# TB2J-OpenMX
[![Build Status](https://travis-ci.com/mailhexu/TB2J-OpenMX.svg?branch=master)](https://travis-ci.com/mailhexu/TB2J-OpenMX)

## Description
The TB2J interface to OpenMX.
TB2J is a open source python package for calculating the magnetic interaction parameters in Heisenberg models from DFT. It use the magnetic force theorem and take the local rigid spin rotation as a perturbation in the Green's function method.

The features include:
 - Calculates  parameters in Heisenberg model, including isotropic exchange, anisotropic exchange, Dyzanoshinskii-Moriya interaction.
 - Can use the input from many DFT codes with Wannier90, e.g. Abinit, Quantum Espresso, Siesta, VASP, etc.
 - Can use input from DFT codes with numerical orbitals from Siesta and OpenMX.
 - Calculate magnon band structure from the Heisenberg Hamiltonian.
 - Generate input for spin dynamics/Monte Carlo codes MULTIBINIT.
 - Require only ground state DFT calculation.
 - No need for supercells.
 - Calculate magnetic interaction up to large distance.
 - Minimal user input, which allows for a black-box like experience and automatic workflows.
 - Versatile API on both the input (DFT Hamiltonian) and the output (Heisenberg model) sides.

For more information, see the documentation on
 <https://tb2j.readthedocs.io/en/latest/> .

The section on TB2J-OpenMX can be found at
 <https://tb2j.readthedocs.io/en/latest/src/openmx.html> .

The link to TB2J:
 <https://github.com/mailhexu/TB2J>

