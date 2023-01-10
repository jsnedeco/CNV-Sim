#!/usr/bin/python

__author__ = 'Abdelrahman Hosny'

import shlex
import subprocess
import os
import logging
import gzip

logger = logging.getLogger(__name__)

def readTargets(filename):
    '''
    Read target file in BED format
    :param filename: target bed file name
    :return: a list of targets (exons). each target is a dictionary of:
             chromosome, start, end, description, score, strand
    '''
    exons = []
    with open(filename, 'r') as tf:
        for line in tf:
            chromosome, start, end = line.strip().split('\t')
            exon = (chromosome, int(start), int(end))
            exons.append(exon)
    return exons

def prepareTargetFile(target_file):
    '''
    sort and merge targets in the target file and writes it to '.sorted.merged'
    :param target_file: target exons in BED format
    :return: sorted and merged file name
    '''
    sorted_file = target_file + ".sorted"
    merged_file = sorted_file + ".merged"

    with open(sorted_file, "w") as f:
        subprocess.call(["sort", "-k1,1", "-k2,2n", target_file], stdout=f)
    with open(merged_file, "w") as f:
        subprocess.call(["bedtools", "merge", "-i", sorted_file], stdout=f)

    return merged_file

def mergeWessimReads(tmp_dir, output_dir):
    '''
    merges the base reads with normal and cnv
    :return: null
    '''
    base_file_1 = os.path.join(tmp_dir, "base_1.fastq")
    base_file_2 = os.path.join(tmp_dir, "base_2.fastq")
    normal_file_1 = os.path.join(tmp_dir, "control_1.fastq")
    normal_file_2 = os.path.join(tmp_dir, "control_2.fastq")
    cnv_file_1 = os.path.join(tmp_dir, "cnv_1.fastq")
    cnv_file_2 = os.path.join(tmp_dir, "cnv_2.fastq")

    cmd1 = f"cat {base_file_1} {normal_file_1} | gzip > {output_dir}/control_1.fastq.gz"
    logger.debug(cmd1)
    subprocess.check_call(shlex.split(cmd1), shell=True)

    cmd2 = f"cat {base_file_2} {normal_file_2} | gzip > {output_dir}/control_2.fastq.gz"
    logger.debug(cmd2)
    subprocess.check_call(shlex.split(cmd2), shell=True)

    cmd3 = f"cat {base_file_1} {cnv_file_1} | gzip > {output_dir}/cnv_1.fastq.fastq.gz"
    logger.debug(cmd3)
    subprocess.check_call(shlex.split(cmd3), shell=True)

    cmd4 = f"cat {base_file_2} {cnv_file_2} | gzip > {output_dir}/cnv_2.fastq.fastq.gz"
    logger.debug(cmd4)
    subprocess.check_call(shlex.split(cmd4), shell=True)

def mergeARTReads(tmp_dir, output_dir):
    '''
    merges the base reads with normal and cnv
    :return: null
    '''
    base_file_1 = os.path.join(tmp_dir, "base1.fq")
    base_file_2 = os.path.join(tmp_dir, "base2.fq")
    normal_file_1 = os.path.join(tmp_dir, "control1.fq")
    normal_file_2 = os.path.join(tmp_dir, "control2.fq")
    cnv_file_1 = os.path.join(tmp_dir, "cnv1.fq")
    cnv_file_2 = os.path.join(tmp_dir, "cnv2.fq")

    cmd1 = f"cat {base_file_1} {normal_file_1} | gzip > {output_dir}/control_1.fastq.gz"
    logger.debug(cmd1)
    subprocess.check_call(shlex.split(cmd1), shell=True)

    cmd2 = f"cat {base_file_2} {normal_file_2} | gzip > {output_dir}/control_2.fastq.gz"
    logger.debug(cmd2)
    subprocess.check_call(shlex.split(cmd2), shell=True)

    cmd3 = f"cat {base_file_1} {cnv_file_1} | gzip > {output_dir}/cnv_1.fastq.fastq.gz"
    logger.debug(cmd3)
    subprocess.check_call(shlex.split(cmd3), shell=True)

    cmd4 = f"cat {base_file_2} {cnv_file_2} | gzip > {output_dir}/cnv_2.fastq.fastq.gz"
    logger.debug(cmd4)
    subprocess.check_call(shlex.split(cmd4), shell=True)

def clean(tmp_dir):
    subprocess.call(["rm", "-rf", tmp_dir + "/*.fq"])
