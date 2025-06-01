# CNAlign/cli.py

import argparse
import pandas as pd
from CNAlign import CNAlign

def main():
    parser = argparse.ArgumentParser(
        description='CNAlign: use optimal segment alignment to find purity and ploidy values for multi-region sampled bulk tumor genomes'
    )

    # Positional arguments
    parser.add_argument('input', type=str, help='Path to input tab-delimited file (columns: sample, segment, logR, BAF, GC, mb)')
    parser.add_argument('output', type=str, help='Path to output tab-delimited file (columns: Variable, Solution_1, Solution_2, ... Solution_N)')
    parser.add_argument('gurobi_license', type=str, help='Path to Gurobi WLS license file')

    # Optional arguments with defaults
    parser.add_argument('--min_ploidy', type=float, default=1.6, help='Minimum sample ploidy [1.6]')
    parser.add_argument('--max_ploidy', type=float, default=6.0, help='Maximum sample ploidy [6.0]')
    parser.add_argument('--min_purity', type=float, default=0.05, help='Minimum sample purity [0.05]')
    parser.add_argument('--max_purity', type=float, default=0.95, help='Maximum sample purity [0.95]')
    parser.add_argument('--min_aligned_seg_mb', type=float, default=5.0, help='Minimum length (Mb) for a segment in the CNA alignment [5.0]')
    parser.add_argument('--max_homdel_mb', type=float, default=100.0, help='Maximum combined length (Mb) of segments with homozygous-deletions [100.0]')
    parser.add_argument('--delta_tcn_to_int', type=float, default=0.2, help='Maximum distance between a sample\'s TCN value at a given segment and the its nearest integer for the sample to be \"aligned\" [0.2]')
    parser.add_argument('--delta_tcn_to_avg', type=float, default=0.1, help='Maximum distance between a sample\'s TCN value at a given segment and the average of all TCN values at that segment for the sample to be \"aligned\" [0.1]')
    parser.add_argument('--delta_tcnavg_to_int', type=float, default=0.1, help='Maximum distance between the average of all TCN values at a segment and its nearest integer for samples to be \"aligned\" at this segment [0.1]')
    parser.add_argument('--delta_mcn_to_int', type=float, default=0.2, help='Maximum distance between a sample\'s MCN value at a given segment and the its nearest integer for the sample to be \"aligned\" [0.2]')
    parser.add_argument('--delta_mcn_to_avg', type=float, default=0.1, help='Maximum distance between a sample\'s MCN value at a given segment and the average of all MCN values at that segment for the sample to be \"aligned\" [0.1]')
    parser.add_argument('--delta_mcnavg_to_int', type=float, default=0.1, help='Maximum distance between the average of all MCN values at a segment and its nearest integer for samples to be \"aligned\" at this segment [0.1]')
    parser.add_argument('--mcn_weight', type=float, default=0.5, help='Weight of MCN component in Objective 2 such that tcn_weight+mcn_weight=1 [0.5]')
    parser.add_argument('--rho', type=float, default=1.0, help='Minimum fraction of samples with matching TCN and MCN values for a segment to be considered \"aligned\" [1.0]')
    parser.add_argument('--timeout', type=int, default=3*60, help='Time (s) without improvement for optimization to stop [180]')
    parser.add_argument('--min_cna_segments_per_sample', type=int, default=5, help='Minimum number of segments with CNAs per sample for solution to be valid (prevents trivial solutions) [5]')
    parser.add_argument('--obj2_clonalonly', type=bool, default=False, help='Optimize obj2 only among segments with clonal CNAs [False]')
    parser.add_argument('--sol_count', type=int, default=10, help='Top N solutions to return [10]')
    args = parser.parse_args()

    # print out message with input parameters 
    print('CNAlign command-line interface mode:')
    print('- input: '+args.input)
    print('- output: '+args.output)

    # load data
    dat = pd.read_csv(args.input, sep='\t')

    # run CNAlign
    model = CNAlign(
            dat=dat,
            gurobi_license=args.gurobi_license,
            min_ploidy=args.min_ploidy,
            max_ploidy=args.max_ploidy,
            min_purity=args.min_purity,
            max_purity=args.max_purity,
            min_aligned_seg_mb=args.min_aligned_seg_mb,
            max_homdel_mb=args.max_homdel_mb,
            delta_tcn_to_int=args.delta_tcn_to_int,
            delta_tcn_to_avg=args.delta_tcn_to_avg,
            delta_tcnavg_to_int=args.delta_tcnavg_to_int,
            delta_mcn_to_int=args.delta_mcn_to_int,
            delta_mcn_to_avg=args.delta_mcn_to_avg,
            delta_mcnavg_to_int=args.delta_mcnavg_to_int,
            mcn_weight=args.mcn_weight,
            rho=args.rho,
            timeout=args.timeout,
            min_cna_segments_per_sample=args.min_cna_segments_per_sample,
            obj2_clonalonly=args.obj2_clonalonly,
            sol_count=args.sol_count
            )

    n_objectives = 2
    solution_count = model.SolCount
    
    # extract objective values for each solution
    obj1_vals = []
    obj2_vals = []
    for i in range(model.SolCount):
        model.params.SolutionNumber = i
        obj1_vals.append(model._obj1_expr.Xn)
        obj2_vals.append(
            -(1 - args.mcn_weight) * model.getVarByName('tcn_error_clonal').Xn
            - args.mcn_weight * model.getVarByName('mcn_error_clonal').Xn
        )

    df = pd.DataFrame({
        'Obj1': obj1_vals,
        'Obj2': obj2_vals
        }, index=[f"Solution_{i+1}" for i in range(model.SolCount)])
    obj_df = df.T.reset_index()
    obj_df.columns = ['Variable'] + [f"Solution_{i+1}" for i in range(model.SolCount)]

    # extract ploidy for each solution
    pl_vars = [v for v in model.getVars() if v.VarName.startswith("pl[")]
    pl_var_names = [v.VarName for v in pl_vars]
    pl_data = { "Variable": pl_var_names }
    for i in range(solution_count):
        model.params.SolutionNumber = i
        pl_data[f"Solution_{i+1}"] = [v.Xn for v in pl_vars]
    
    pl_df = pd.DataFrame(pl_data) 

    # extract purity for each solution
    pu_vars = [v for v in model.getVars() if v.VarName.startswith("z[")]
    pu_var_names = [v.VarName for v in pu_vars]
    pu_data = { "Variable": pu_var_names }
    for i in range(solution_count):
        model.params.SolutionNumber = i
        pu_data[f"Solution_{i+1}"] = [1/v.Xn for v in pu_vars]
    
    pu_df = pd.DataFrame(pu_data)
    pu_df.Variable = pu_df.Variable.replace('z\\[','pu[',regex=True) 

    # extract allmatch for each solution
    allmatch_vars = [v for v in model.getVars() if v.VarName.startswith("allmatch[")]
    allmatch_var_names = [v.VarName for v in allmatch_vars]
    allmatch_data = { "Variable": allmatch_var_names }
    for i in range(solution_count):
        model.params.SolutionNumber = i
        allmatch_data[f"Solution_{i+1}"] = [v.Xn for v in allmatch_vars]
    
    allmatch_df = pd.DataFrame(allmatch_data)

    # extract tcn for each solution
    tcn_vars = [v for v in model.getVars() if v.VarName.startswith("tcn[")]
    tcn_var_names = [v.VarName for v in tcn_vars]
    tcn_data = { "Variable": tcn_var_names }
    for i in range(solution_count):
        model.params.SolutionNumber = i
        tcn_data[f"Solution_{i+1}"] = [v.Xn for v in tcn_vars]
    
    tcn_df = pd.DataFrame(tcn_data)

    # extract mcn for each solution
    mcn_vars = [v for v in model.getVars() if v.VarName.startswith("mcn[")]
    mcn_var_names = [v.VarName for v in mcn_vars]
    mcn_data = { "Variable": mcn_var_names }
    for i in range(solution_count):
        model.params.SolutionNumber = i
        mcn_data[f"Solution_{i+1}"] = [v.Xn for v in mcn_vars]
    
    mcn_df = pd.DataFrame(mcn_data)

    # merge the output data and write out
    print('Writing output to file: '+args.output)
    out = pd.concat([obj_df, pl_df, pu_df, allmatch_df, tcn_df, mcn_df], axis=0, ignore_index=True)
    out.iloc[:, 1:] = out.iloc[:, 1:].round(6)
    out.to_csv(args.output, sep='\t', index=False)
    print('Done. Have a nice day!')


