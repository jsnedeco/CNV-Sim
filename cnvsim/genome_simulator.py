__author__ = 'Abdelrahman Hosny'

import random
import os
import shlex
import subprocess
import datetime
import shutil
from Bio import SeqIO
from collections import namedtuple

from . import fileio
import logging
logger = logging.getLogger(__name__)

def getScriptPath():
    return os.path.dirname(os.path.realpath(__file__))


def _generateCNVMask(mask_length, p_amplify, p_delete, min_variation, max_variation):
    '''
    This function generates random Copy Number Variations mask list
    :param exons: list of regions.
    :param p_amplify: percentage of amplifications introduced
    :param p_delete: percentage of deletions introduced
    :return: a list of the same length as the exons list. each list item
             indicates the variation to be added to the exon in the same position.
             Positive number: number of amplification
             Zero: no change
             -1: delete this exon
    '''

    number_of_amplifications = int(p_amplify * mask_length)
    number_of_deletions = int(p_delete * mask_length)
    cnv_mask = [0] * mask_length

    # generate CNV mask (a list of amplifications and deletions to be applied to the exons)
    while number_of_amplifications > 0:
        choice = random.randrange(0, mask_length)
        while cnv_mask[choice] != 0:
            choice = random.randrange(0, mask_length)
        cnv_mask[choice] = random.randrange(min_variation, max_variation)     # random amplifications in the range [min_variation, max_variation)
        number_of_amplifications -= 1
    random.shuffle(cnv_mask)
    while number_of_deletions > 0:
        choice = random.randrange(0, mask_length)
        while cnv_mask[choice] != 0:
            choice = random.randrange(0, mask_length)
        cnv_mask[choice] = -1*random.randrange(min_variation, max_variation)  # random deletions in the range [min_variation, max_variation)
        number_of_deletions -= 1
    random.shuffle(cnv_mask)

    return cnv_mask


def _simulateCNV(genome_file, cnv_list, read_length, control_file, cnv_file):
    '''
    Simulate the control genome and the CNV genome
    :param genome: original genome sequence
    :param cnv_list: a list of region variations (chromosome, start, end, variation)
    :return: control_genome, cnv_genome
    '''

    logger.info("loading genome file ..")
    genome = SeqIO.index(genome_file, "fasta")
    logger.info(f"successfully loaded a genome with {len(genome)} chromosomes")
    with open(control_file, 'w') as control_handle:
        with open(cnv_file, 'w') as cnv_handle:
            for cnv in cnv_list:
                header = f">{cnv.chromosome}:{cnv.region_start}:{cnv.region_end}:{cnv.variation}"
                curr_chrom = genome[cnv.chromosome]
                sequence = curr_chrom[cnv.region_start:cnv.region_end].seq
                prefix = curr_chrom[cnv.region_start-read_length: cnv.region_start].seq
                suffix = curr_chrom[cnv.region_end:cnv.region_end+read_length].seq

                if cnv.variation > 0:
                    # amplification
                    amplification = prefix + sequence * cnv.variation + suffix
                    cnv_handle.write(header + "\n")
                    cnv_handle.write(amplification + "\n")
                elif cnv.variation < 0:
                    # deletion
                    deletion = prefix + sequence * abs(cnv.variation) + suffix
                    control_handle.write(header + "\n")
                    control_handle.write(deletion + "\n")

def _callART(genome_file, output_file, read_length, fold_coverage=1):
    '''
    Calls Wessim to generate artificial reads for the targets on the reference genome
    :param genome_file: reference genome file in FASTA format
    :param output_file: output file name
    :param read_length: the read length
    :param fold_coverage: fold coverage for the reads
    :return: None
    '''
    art_path = os.path.join(getScriptPath(), "ART", 'art_illumina')

    cmd = f"""{art_path} \
                -na \
                -i {genome_file} \
                -p \
                -m 200 \
                -s 10 \
                -l {read_length} \
                -f {fold_coverage} \
                -o {output_file}"""
    logger.debug(cmd)

    subprocess.call(shlex.split(cmd), stderr=None)
    os.chdir("..")


def simulate_genome_cnv(simulation_parameters):
    '''
    Simulate copy number variations on the passed reference genome based on the given simulation parameters
    :param simulation_parameters: a dictionary containing parameters for simulation
    :param cnv_list_parameters: a dictionary containing parameters for CNV List creation
    :return: None
    '''
    logger.info('simulation type: whole genome sequencing')

    # create a temporary directory for intermediate files
    if not os.path.exists(simulation_parameters['tmp_dir']):
        os.makedirs(simulation_parameters['tmp_dir'])

    # create output directory for final results
    if not os.path.exists(simulation_parameters['output_dir']):
        os.makedirs(simulation_parameters['output_dir'])

    genome_file = simulation_parameters['genome_file']

    # initialize variables for temporary files
    control_genome_file = os.path.join(simulation_parameters['tmp_dir'], "ControlGenome.fa")
    cnv_genome_file = os.path.join(simulation_parameters['tmp_dir'], "CNVGenome.fa")
    base_reads_file = os.path.join(simulation_parameters['tmp_dir'], "base")
    control_reads_file = os.path.join(simulation_parameters['tmp_dir'], "control")
    cnv_reads_file = os.path.join(simulation_parameters['tmp_dir'], "cnv")

    logger.info("loading CNV list ..")
    CNV = namedtuple("CNV", ["chromosome", "region_start", "region_end", "num_positions", "variation"])


    with open(simulation_parameters['cnv_list_file'], "r") as f:
        cnv_list = []
        lines = f.readlines()
        lines.pop(0)
        for line in lines:
            chromosome, region_start, region_end, num_positions, variation = line.strip().split("\t")
            cnv = CNV(chromosome, int(region_start), int(region_end), int(num_positions), int(variation))
            cnv_list.append(cnv)

        logger.info("successfully loaded CNV list that contains " + str(len(cnv_list)) + " regions ..")

    # call ART to generate reads from the genome file
    logger.info("generating reads for the genome ..")
    logger.info("delegating job to ART ...")
    _callART(genome_file, base_reads_file, simulation_parameters['read_length'], fold_coverage=simulation_parameters['coverage'])

    logger.info("simulating copy number variations (amplifications/deletions)")
    _simulateCNV(genome_file, cnv_list, simulation_parameters['read_length'], control_genome_file, cnv_genome_file)

    logger.info("delegating job to ART ...")
    _callART(control_genome_file, control_reads_file, simulation_parameters['read_length'], fold_coverage=simulation_parameters['coverage'])

    logger.info("delegating job to ART ...")
    _callART(cnv_genome_file, cnv_reads_file, simulation_parameters['read_length'], fold_coverage=simulation_parameters['coverage'])

    logger.info("merging results ..")
    fileio.mergeARTReads(simulation_parameters['tmp_dir'], simulation_parameters['output_dir'])

    logger.info("cleaning temporary files ..")
    fileio.clean(simulation_parameters['tmp_dir'])

    logger.info("simulation completed. find results in " + simulation_parameters['output_dir'])