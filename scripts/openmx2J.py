#!/usr/bin/env python3
import argparse
import sys

from TB2J.exchange_params import add_exchange_args_to_parser, parser_argument_to_dict
from TB2J.versioninfo import print_license
from TB2J_OpenMX.gen_exchange import gen_exchange


def run_openmx2J():
    print_license()
    print("\n")
    parser = argparse.ArgumentParser(
        description="openmx2J: Using magnetic force theorem to calculate exchange parameter J from openmx Hamiltonian"
    )
    parser.add_argument(
        '--prefix', help="prefix of the openmx files", default='openmx', type=str)

    parser = add_exchange_args_to_parser(parser)

    args = parser.parse_args()

    if args.elements is None and args.index_magnetic_atoms is None:
        print("Please input the magnetic elements, e.g. --elements Fe Ni")
        sys.exit()

    kwargs = parser_argument_to_dict(args)
    gen_exchange(path='./', prefix=args.prefix, **kwargs)


if __name__ == "__main__":
    run_openmx2J()
