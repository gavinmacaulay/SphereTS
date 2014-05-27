# -*- coding: utf-8 -*-
"""    
    Copyright 2014 Gavin Macaulay

    This file is part of SphereTS.

    SphereTS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SphereTS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SphereTS.  If not, see <http://www.gnu.org/licenses/>.
"""
import os
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = 'sphereTS',
    version = '1.0',
    description = 'Calculates the acoustic target strength of an elastic sphere immersed in a fluid',
    long_description = read('README.txt'),
    author = 'Gavin Macaulay', 
    author_email = 'gavin@macaulay.co.nz',
    scripts = ['sphereTS.py', 'sphereTSGUI.py'],
    packages = find_packages(),
    install_requires = ['numpy', 'traitsui', 'matplotlib', 'traits', 'scipy', 'gsw'],
    keywords = 'acoustic TS target strength sphere calibration',
    url = 'https://bitbucket.org/gjm/spherets',
    classifiers = ('Development Status :: 4 - Beta',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python',
                   'Topic :: Scientific/Engineering'),
    )  
