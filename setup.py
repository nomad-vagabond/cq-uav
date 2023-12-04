from setuptools import setup, find_packages

setup(
    name='cquav',
    version='0.0.2',
    author = "Vadym Pasko",
    author_email = "vadym@splinecloud.com",
    description = "A CaqQuery library for generating UAV components.",
    license = "MIT",
    packages=find_packages(),  
    package_data = {"cquav": ["wing/airfoil/*.json"]},
    include_package_data=True,
    install_requires=[
        'numpy',
        'scipy',
        'pandas',
        'cadquery',
        'requests',
    ],
)
