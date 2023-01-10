#!/usr/bin/python

__author__ = 'Abdelrahman Hosny'

import os.path
import datetime
import argparse
import shutil
from pathlib import Path

from cnvsim.fileio import *
from cnvsim.genome_simulator import *

logger = logging.getLogger(__name__)

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.DEBUG)

class CapitalisedHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def add_usage(self, usage, actions, groups, prefix=None):
        if prefix is None:
            prefix = 'Usage: '
            return super(CapitalisedHelpFormatter, self).add_usage(usage, actions, groups, prefix)

def main():
    parser = argparse.ArgumentParser(add_help=True, formatter_class=CapitalisedHelpFormatter, \
                                     description='Generates NGS short reads that encompass copy number variations in whole genome and targeted exome sequencing')
    parser._positionals.title = 'Positional arguments'
    parser._optionals.title = 'Optional arguments'
    parser.add_argument('-v', '--version', action='version', version = 'CNV-Sim v0.9.2', help = "Show program's version number and exit.")

    parser.add_argument("-g", "--genome", type=Path, required=True,
                        help="path to the referece genome file in FASTA format ")

    parser.add_argument("-o", "--output_dir_name",type=str, default="simulation_output", \
                        help="a name to be used to create the output directory (overrides existing directory with the same name).")
    parser.add_argument("-n", "--n_reads", type=int, default=10000, \
                        help="total number of reads without variations")
    parser.add_argument("-l", "--read_length", type=int, default=100, \
                        help="read length (bp)")
    parser.add_argument("--cnv_list", type=Path, required=True, \
                        help="path to a CNV list file in BED format chr | start | end | variation.")
    parser.add_argument("--coverage", type=int, default=1, \
                        help="the integer average depth of coverage of a genome for the reads (only on whole genome simulation)")

    args = parser.parse_args()

    simulation_parameters = {}
    simulation_parameters['genome_file'] = args.genome
    simulation_parameters['output_dir'] = os.path.join(os.getcwd(), args.output_dir_name)
    simulation_parameters['number_of_reads'] = args.n_reads
    simulation_parameters['read_length'] = args.read_length
    if args.cnv_list is not None:
        simulation_parameters['cnv_list_file'] = args.cnv_list
    else:
        simulation_parameters['cnv_list_file'] = None
    simulation_parameters['tmp_dir'] = os.path.join(args.output_dir_name , "tmp")
    simulation_parameters['coverage'] = args.coverage

    simulate_genome_cnv(simulation_parameters)



if __name__ == '__main__':
    main()