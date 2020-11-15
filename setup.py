#!/usr/bin/env python
from setuptools import setup, find_packages
import sys

__version__ = "0.3.3"

long_description = """Interface of TB2J to OpenMX. 
TB2J is a Python package aimed to compute automatically the magnetic interactions (superexchange  and Dzyaloshinskii-Moriya) between atoms of magnetic crystals from DFT Hamiltonian based on Wannier functions or Linear combination of atomic orbitals. It uses the Green's function method and take the local rigid spin rotation as a perturbation. The package can take the output from Wannier90, which is interfaced with many density functional theory codes or from codes based on localised orbitals. A minimal user input is needed, which allows for an easily integration into a high-throughput workflows.
 """

setup(
    name='TB2J_OpenMX',
    version=__version__,
    description=
    'TB2J_OpenMX: TB2J interface to OpenMX',
    long_description=long_description,
    author='Xu He',
    author_email='mailhexu@gmail.com',
    license='GPLv3',
    packages=find_packages(),
    cffi_modules=["TB2J_OpenMX/ffimod.py:ffi"],
    scripts=["scripts/openmx2J.py"],
    install_requires=['numpy', 'scipy', 'matplotlib', 'ase', 'progressbar', 'TB2J>=0.3.3', 'cffi'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3',
        #'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    python_requires='>=3.6',
)
