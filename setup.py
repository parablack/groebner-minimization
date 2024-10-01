from setuptools import find_packages, setup

setup(
    name='groebner_min',
    version='1.0',
    author='Simon Schwarz, Nicolas Faross',
    author_email='sschwarz@mpi-inf.mpg.de',
    description='Minimization of Boolean functions using Groebner bases.',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'groebner-min=groebner_min.groebner_min:main'
        ]
    },
)