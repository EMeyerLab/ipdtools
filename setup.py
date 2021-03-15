#!/usr/bin/env python

from setuptools import setup, Extension, find_packages
import os
import sys

setup(
    name='ipdtools',
    version='0.2.1',
    author='DELEVOYE Guillaume',
    license=open('LICENSES.txt').read(),
    packages=find_packages("."),
    python_requires='>=3.6',
    package_data={'ipdtools': ['resources/*.h5']},
    ext_modules=[Extension('ipdtools/tree_predict', ['ipdtools/tree_predict.c'],
                           extra_compile_args=["-O3", "-shared", "-std=c99"],
                           export_symbols=["innerPredict", "innerPredictCtx", "init_native"])],
    install_requires=[
        'numpy',
        'joblib',
        'h5py',
        'pandas',
        'tqdm',
        'h5py',
        'numpy'
    ],
    entry_points={'console_scripts': [
        "ipdtools = ipdtools.launchers.ipdtools_launcher:main",
    ]},
)
