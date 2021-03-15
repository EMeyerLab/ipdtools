#!/usr/bin/env python

from setuptools import setup, Extension, find_packages
import os
import sys

setup(
    name='ipdtools',
    version='0.22',
    description="PacBio codes isolated and repackaged to work on their in-sillico model",
    long_description="PacBio codes isolated and repackaged to work on their in-sillico model. See https://github.com/GDelevoye/ipdtools for more details.",
    author='DELEVOYE Guillaume',
    author_email="delevoye.guillaume@gmail.com",
    packages=find_packages("."),
    url="https://github.com/GDelevoye/ipdtools",
    python_requires='>=3.6',
    package_data={'ipdtools': ['resources/*.h5']},
    ext_modules=[Extension('ipdtools/tree_predict', ['ipdtools/tree_predict.c'],
                           extra_compile_args=["-O3", "-shared", "-std=c99"],
                           export_symbols=["innerPredict", "innerPredictCtx", "init_native"])],
    entry_points={'console_scripts': [
        "ipdtools = ipdtools.launchers.ipdtools_launcher:main",
    ]},
)
