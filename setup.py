from setuptools import setup

requirements = [
    'numpy==1.23*;',
    'pandas==1.5*',
    'matplotlib==3.7*',
    'seaborn==0.12*',
    'scikit-allel==1.3*',
    'scikit-learn==1.2*',
    'psutil==5.9*',
    'adjusttext==0.7*',
    'bcftools==1.17*',
    'plink==1.90*',
    'gemma==0.98.3'
]

setup(
    name='vcf2gwas',
    version='0.8.9',
    description="Python API for comprehensive GWAS analysis using GEMMA",
    license="GNUv3",
    author="Frank Vogt",
    author_email='frvogt@gmail.com',
    url='https://github.com/frankvogt/vcf2gwas',
    packages=['vcf2gwas'],
    package_dir={'vcf2gwas': 'vcf2gwas'},
    package_data={'vcf2gwas': ['*.yml', 'input/example.*', 'GFF_files/GFF*', 'README.*', 'LIC*']},
    entry_points={
        'console_scripts': [
            'vcf2gwas=vcf2gwas.__main__:run_main'
            #'install_vcf2gwas=vcf2gwas.install:main'
        ]
    },
    install_requires=requirements,
    keywords='vcf2gwas',
    classifiers=[
        'Programming Language :: Python :: 3.9'
    ]
)
