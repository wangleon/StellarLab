#!/usr/bin/env python3
from distutils.core import setup

setup(
    name = 'StellarLab',
    version = '0.1',
    description = '',
    author = 'Liang Wang',
    author_email = 'wang.leon@gmail.com',
    license = 'BSD',
    scripts = [
        'scripts/findstar',
        ],
    packages = [
        'stellarlab',
        'stellarlab/catalog',
        'stellarlab/evolution',
        'stellarlab/extinction',
        'stellarlab/parameter',
        'stellarlab/spectrum',
        'stellarlab/kinetics',
        'stellarlab/utils',
        ],
    )
