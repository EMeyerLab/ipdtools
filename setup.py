from setuptools import setup, Extension, find_packages
import os
import sys

setup(
    name='ipdtools',
    version='0.1',
    author='DELEVOYE Guillaume',
    license=open('LICENSES.txt').read(),
    packages=find_packages("."),
    package_data={'ipdtools': ['resources/*.h5']},
    ext_modules=[Extension('ipdtools/tree_predict', ['ipdtools/tree_predict.c'],
                           extra_compile_args=["-O3", "-shared", "-std=c99"],
                           export_symbols=["innerPredict", "innerPredictCtx", "init_native"])],
    zip_safe=False,
    install_requires=[
        'numpy'
    ],
    entry_points={'console_scripts': [
        "ipdtools = ipdtools.launchers.ipdtools_launcher:main",
    ]},
)
