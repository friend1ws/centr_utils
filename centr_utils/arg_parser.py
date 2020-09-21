#! /usr/bin/env python

import argparse

from .run import *
from .version import __version__

def create_parser():

    parser = argparse.ArgumentParser(prog = "centr_utils")

    parser.add_argument("--version", action = "version", version = "%(prog)s " + __version__)

    subparsers = parser.add_subparsers()

    ##########


    split_to_monomer = subparsers.add_parser("split_to_monomer",
                                              help = "Split fasta sequences to alpha satellite monomer elements")

    split_to_monomer.add_argument("input_fasta_file", default = None, type = str,
                                    help = "Path to FASTA file")

    split_to_monomer.add_argument("output_file", default = None, type = str,
                                  help = "Patht to output file")

    split_to_monomer.add_argument("hmm_model_fwd", default = None, type = str,
                                   help = "Path to forward HMM file")

    split_to_monomer.add_argument("hmm_model_rev", default = None, type = str,
                                  help = "Patht to reverse HMM file")

    split_to_monomer.add_argument("--min_monomer_len", default = 150, type = int,
                                  help = "Minimum monomer length")

    split_to_monomer.set_defaults(func = split_to_monomer_main)

    ##########

    return parser

