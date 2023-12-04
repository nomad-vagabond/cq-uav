from setuptools import setup, find_packages

setup(
    name='cquav',
    version='0.0.1',
    author = "Vadym Pasko",
    author_email = "vadym@splinecloud.com",
    description = "A CaqQuery library for eneratingUAV components.",
    license = "MIT",
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'pandas',
        'cadquery',
        'requests',
    ],
)
