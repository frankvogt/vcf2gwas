from setuptools import setup

requirements = [
    'numpy==1.20*',
    'pandas==1.2*',
    'matplotlib==3.4*',
    'seaborn==0.11*',
    'scikit-allel==1.3*',
    'scikit-learn==0.24*',
    'psutil==5.8*',
    'adjusttext==0.7*',
    'bcftools==1.12*',
    'plink==1.90*',
    'gemma==0.98.3'
    # package requirements go here
]

setup(
    name='vcf2gwas',
    version='0.8.4',
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
            'vcf2gwas=vcf2gwas.__main__:main'
            #'install_vcf2gwas=vcf2gwas.install:main'
        ]
    },
    install_requires=requirements,
    keywords='vcf2gwas',
    classifiers=[
        'Programming Language :: Python :: 3.9'
    ]
)
