#!/usr/bin/env python3
import argparse
from TB2J.versioninfo import print_license
from TB2J_OpenMX.gen_exchange import gen_exchange

def run_openmx2J():
    print_license()
    print("\n")
    parser = argparse.ArgumentParser(
        description=
        "openmx2J: Using magnetic force theorem to calculate exchange parameter J from openmx Hamiltonian"
    )
    parser.add_argument(
        '--prefix', help="prefix of the openmx files", default='openmx', type=str)
    parser.add_argument(
        '--elements',
        help="elements to be considered in Heisenberg model",
        default=None,
        type=str,
        nargs='*')
    parser.add_argument(
        '--rcut', help='range of R. The default is all the commesurate R to the kmesh', default=None, type=float)
    parser.add_argument(
        '--efermi', help='Fermi energy in eV', default=None, type=float)
    parser.add_argument(
        '--kmesh',
        help='kmesh in the format of kx ky kz. Monkhorst pack. If all the numbers are odd, it is Gamma cenetered. (strongly recommended)',
        type=int,
        nargs='*',
        default=[5, 5, 5])
    parser.add_argument(
        '--emin',
        help='energy minimum below efermi, default -14 eV',
        type=float,
        default=-14.0)
    parser.add_argument(
        '--emax',
        help='energy maximum above efermi, default 0.0 eV',
        type=float,
        default=0.05)
    parser.add_argument(
        '--use_cache',
        help="whether to use disk file for temporary storing wavefunctions and hamiltonian to reduce memory usage. Default: False",
        action='store_true',
        default=False)
    parser.add_argument(
        '--nz', help='number of integration steps, default: 50', default=50, type=int)

    parser.add_argument(
        '--exclude_orbs',
        help=
        "the indices of wannier functions to be excluded from magnetic site. counting start from 0",
        default=[],
        type=int,
        nargs='+')

    parser.add_argument(
        "--description",
        help=
        "add description of the calculatiion to the xml file. Essential information, like the xc functional, U values, magnetic state should be given.",
        type=str,
        default="Calculated with TB2J.\n")

    parser.add_argument(
        "--fname",
        default='exchange.xml',
        type=str,
        help='exchange xml file name. default: exchange.xml')

    args = parser.parse_args()

    if args.elements is None:
        print("Please input the magnetic elements, e.g. --elements Fe Ni")
        exit()
    gen_exchange(
        path='./',
        prefix=args.prefix,
        kmesh=args.kmesh,
        magnetic_elements=args.elements,
        Rcut=args.rcut,
        emin=args.emin,
        emax=args.emax,
        nz=args.nz,
        description=args.description,
        use_cache=args.use_cache,
        exclude_orbs=args.exclude_orbs)

if __name__ == "__main__":
    run_openmx2J()
