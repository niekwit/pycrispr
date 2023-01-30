#!/usr/bin/env python3

import subprocess
import os
import click
import sys
import hashlib
import glob
import numpy as np
import pandas as pd
import yaml
from itertools import compress
from tqdm.auto import tqdm

script_dir = os.path.abspath(os.path.dirname(__file__))
work_dir = os.getcwd()

####generic functions
def write2log(command,name=""):
    '''Write bash command to commands.log
    '''
    with open(os.path.join(work_dir,"commands.log"), "a") as file:
        file.write(name)
        print(*command, sep = "", file = file)

def checkInPath(check):
    ''' Check if commands in list are set in $PATH
    '''
    if type(check) != list: #single command (string)
        check = [check]
        
    path = os.environ["PATH"].lower()
    if any([x in path for x in check]) == False:
        bool_list = [x in path for x in check]
        bool_list_rev = [not x for x in bool_list]
        not_in_path= " ".join((list(compress(check, bool_list_rev))))
        click.secho(f"ERROR: {not_in_path} not found in $PATH",fg="red")
        return(False)

####command line parser
@click.group()
def cli():
    """CRISPR-Cas9 screen analysis pipeline
    """
        
@click.command(name='version')
def version():
    """Report the current build and version number
    """
    click.echo()

@click.command(name='show-libs')
def show_libs():
    ''' Displays which sgRNA libraries have been written to crispr.yaml
    '''
    lib = os.path.join(script_dir,"crispr.yaml")
    click.secho("The following sgRNA libraries are available in crispr.yaml:",fg="green")
    subprocess.run(["cat", lib])

@click.command(name='add-lib')
@click.option("-n","--name", required=True,
              help="Library name")
@click.option("-i","--index", required=True,
              help="HISAT2 path")
@click.option("-f","--fasta", required=True,
              help="Fasta file path")
@click.option("-c","--csv", required=True,
              help="CSV file path")
@click.option("--sg-length", required=True, type=int,
              help="sgRNA length")
@click.option("--species", required=True, 
              help="Target genome")
def add_lib(name,index,fasta,csv,sg_length,species):
    ''' Add sgRNA library to crispr.yaml
    '''
    pass

@click.command(name='analysis')
@click.option("--md5sums", is_flag=True, show_default=True, default=False,
              help="Check md5sums of fastq files")
@click.option("--fastqc", is_flag=True, show_default=True, default=False,
              help="Quality control of fastq files")
@click.option("-r", "--rename", is_flag=True, show_default=True, default=False, 
              help="Rename fastq files according to rename.txt")
@click.option("-t","--threads", default=1, type=int, 
              help="Number of CPU threads used during analysis")
@click.option("-l","--library", 
                  help="CRISPR-Cas9 library")
@click.option("-m","--mismatch", default=0, show_default=True, type=int, 
              help="Number of mismatches allowed during alignment")
@click.option("-a","--analysis", default="mageck", show_default=True, 
              type=click.Choice(["mageck","bagel2"]),
              help="Statistical analysis with MAGeCK or BAGEL2")
@click.option("-c","--cnv", default=None, show_default=True, 
              type=str,
              help="Apply CNV correction for MAGeCK/BAGEL2 for given cell line")
@click.option("-f","--fdr", default=0.25, show_default=True, 
              type=float,
              help="FDR cutoff for MAGeCK")
@click.option("--go", is_flag=True, show_default=True, default=False,
              help="Apply gene ontology analysis to MAGeCK and/or BAGEL2 results")
def analysis(md5sums,fastqc,rename,threads,library,mismatch,analysis,cnv,fdr,go):
    ''' Run CRISPR-Cas9 screen analysis pipeline
    '''
    #create output dirs
    os.makedirs(os.path.join(work_dir,"count"),exist_ok=True)
    if analysis == "mageck":
        os.makedirs(os.path.join(work_dir,"mageck"),exist_ok=True)
    elif analysis == "bagel2":
        os.makedirs(os.path.join(work_dir,"mageck"),exist_ok=True)
    
    #run selected functions
    '''
    if md5sums == True:
        click.echo("Checking md5sums of fastq files in raw-data/")
        md5sum_match = utils.md5sums()
        if md5sum_match == False:
            click.secho("ERROR: At least one calculated md5sum did not match the pre-calculated ones\nPlease check md5sums_failed.csv",color="red")
    click.echo(f"{library} library selected")
    click.echo(f"Mismatches allowed: {mismatch}")
    if rename == True:
        click.echo("Renaming fastq files according to rename.txt")
        utils.rename()
    '''
    if fastqc:
        click.echo("Quality control of fastq files using FastQC/MultiQC")
        
        #check if fastqc is in $PATH
        if not checkInPath("fastqc"):
            return()
        
        #create fastqc dir
        fastqc_dir = os.path.join(work_dir,"fastqc")
        os.makedirs(fastqc_dir,exist_ok = True)
        
        #run fastqc
        data_files = os.path.join(work_dir,"raw-data","*.gz")
        fastqc = ["fastqc","--threads",threads,"--quiet","-o",fastqc_dir,data_files]
        write2log(" ".join(fastqc))
        subprocess.run(fastqc)
        
        #run multiqc
        multiqc = ["multiqc","-o",fastqc_dir,fastqc_dir]
        write2log(" ".join(multiqc))
        subprocess.run(multiqc)
    '''
    #align and count sgRNA reads
    utils.count(threads,mismatch,library)
    
    #join count files
    utils.join(library)
    
    #create normalised count file
    utils.normalise()
    
    #apply statistics to count files
    if analysis == "mageck":
        utils.mageck(cnv)
    elif analysis == "bagel2":
        utils.bagel2()

    #run gene ontology analysis
    if go == True:
        utils.go()
    '''
#add subparsers
cli.add_command(show_libs)
cli.add_command(add_lib)
cli.add_command(analysis)
cli.add_command(version)

    
        
    
    