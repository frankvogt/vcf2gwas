from setuptools import setup

requirements = [
    # package requirements go here
]

setup(
    name='vcf2gwas',
    version='0.5',
    description="Python API for comprehensive GWAS analysis using GEMMA",
    license="GNUv3",
    author="Frank Vogt",
    author_email='frvogt@gmail.com',
    url='https://github.com/frankvogt/vcf2gwas',
    packages=['vcf2gwas'],
    package_dir={'vcf2gwas': 'vcf2gwas'},
    package_data={'vcf2gwas': ['*.yml', 'input/example.*', 'README.*', 'LIC*']},
    entry_points={
        'console_scripts': [
            'vcf2gwas=vcf2gwas.__main__:main',
            'install_vcf2gwas=vcf2gwas.install:main'
            #'test_vcf2gwas=run_test:main'
        ]
    },
    install_requires=requirements,
    keywords='vcf2gwas',
    classifiers=[
        'Programming Language :: Python :: 3.9'
    ]
)
