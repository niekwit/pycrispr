#!/usr/bin/env python3

import os
import click
from ..scripts import utils as utils
import yaml


#load crispr lib info
#script_dir = os.path.abspath(os.path.dirname(__file__))
#yaml_dir = os.path.join(os.path.dirname(os.path.dirname(script_dir)),"yaml")
#with open(os.path.join(yaml_dir,
#                       "crispr.yaml")) as file:
#    crispr_libs = yaml.full_load(file)
#crispr_library_list = list(crispr_libs.keys())


#command line parser
@click.command()
@click.option("--md5sums", is_flag=True, show_default=True, default=False,
              help="Check md5sums of fastq files")
@click.option("--fastqc", is_flag=True, show_default=True, default=False,
              help="Fastqc/Multiqc on fastq files")
@click.option("-r", "--rename", is_flag=True, show_default=True, default=False, 
              help="Rename fastq files according to rename.txt")
@click.option("-t","--threads", default=1, type=int, 
              help="Number of CPU threads used for analysis")
@click.option("-l","--library", 
              help="CRISPR-Cas9 library")
@click.option("-m","--mismatch", default=0, show_default=True, type=int, 
              help="Number of mismatches allowed during alignment")
@click.option("-a","--analysis", default="mageck", show_default=True, 
              type=click.Choice(["mageck","bagel2"]),
              help="Statistical analysis with MAGeCK or BAGEL2")
@click.option("-c","--cnv", default=None, show_default=True, 
              type=int,
              help="Apply CNV correction for MAGeCK/BAGEL2 with given cell line")
@click.option("-f","--fdr", default=0.25, show_default=True, 
              type=float,
              help="FDR cutoff for MAGeCK")
@click.option("--go", is_flag=True, show_default=True, default=False,
              help="Apply gene ontology analysis to MAGeCK or BAGEL2 results")


def cli(md5sum,fastqc,rename,threads,library,mismatch,analysis,cnv,fdr,go):
    """CRISPR-Cas9 screen analysis"""
    if md5sum == True:
        click.echo("Checking md5sums of fastq files in raw-data/")
    else:
        click.echo(f"{library} library selected")
        click.echo(f"Mismatches allowed: {mismatch}")
        if rename == True:
            click.echo("Renaming fastq files according to rename.txt")
            utils.rename()
        



