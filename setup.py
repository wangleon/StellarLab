#!/usr/bin/env python
from distutils.core import setup

setup(
    name = 'stella',
    version = '0.1',
    description = '',
    author = 'Liang Wang',
    author_email = 'wang.leon@gmail.com',
    license = 'BSD',
    scripts = [
        'scripts/findstar',
        ],
    packages = [
        'stella',
        'stella/catalog',
        'stella/evolution',
        'stella/extinction',
        'stella/parameter',
        'stella/spectroscopy',
        'stella/kinetics',
        'stella/utils',
        ],
    )
