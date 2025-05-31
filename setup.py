from setuptools import setup, find_packages

setup(
    name='CNAlign',
    version='0.3.0',
    description='Use optimal segment alignment to fit purity and ploidy across multi-region sampled bulk tumor genomes',
    author='Alex Gorelick',
    author_email='alexander_gorelick@hms.harvard.edu',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'gurobipy',
    ],
    python_requires='>=3.7',
)
