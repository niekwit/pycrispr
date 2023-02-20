#!/usr/bin/env python3

import subprocess
import os
import click
import yaml
import timeit
import time

from ..scripts import utils as utils

#global variables
script_dir = os.path.abspath(os.path.dirname(__file__))
work_dir = os.getcwd()
version = "1.0"


####command line parser
@click.group()
def cli():
    """CRISPR-Cas9 screen analysis pipeline
    """
    utils.logCommandLineArgs()
    
@click.command(name='version')
def version():
    """Report the current build and version number
    """
    click.echo(version)

@click.command(name='show-libs')
def show_libs():
    ''' Displays which sgRNA libraries have been written to crispr.yaml
    '''
    lib = os.path.join(script_dir,"crispr.yaml")
    if os.path.exists(lib):
        click.secho("The following sgRNA libraries are available in crispr.yaml:",fg="green")
        subprocess.run(["cat", lib])
    else:
        click.echo("WARNING: no sgRNA library has been added yet")
        click.echo("Please run `pycrispr add-lib`")
        return()

@click.command(name='md5sums')
def md5sums():
	click.echo("Checking md5sums of fastq files in raw-data/")
	md5sum_match = utils.md5sums()
	if not md5sum_match:
		click.secho("ERROR: At least one calculated md5sum did not match the pre-calculated ones\nPlease check md5sums_failed.csv",color="red")
		return()

@click.command(name='add-lib')
@click.option("-n","--name", required=True,
              help="Library name")
@click.option("-i","--index", required=True,
              help="HISAT2 index path")
@click.option("-f","--fasta", required=True,
              help="Fasta file path")
@click.option("-c","--csv", required=True,
              help="CSV file path")
@click.option("--sg-length", required=True, type=int,
              help="sgRNA length")

def add_lib(name,index,fasta,csv,sg_length,species):
    ''' Add sgRNA library to crispr.yaml
    '''
    #create dictionary keys for yaml dump
    yaml_keys = ["fasta","index","csv","sg_length","species"]
    
    #write/create crispr.yaml file
    yaml_file = os.path.join(script_dir,"crispr.yaml")
    if not os.path.exists(yaml_file):
        #no pre-existing yaml file so start with dict to create one
        doc = {}
        doc[name] = {}
        doc[name]["index"] = index
        doc[name]["fasta"] = fasta
        doc[name]["csv"] = csv
        doc[name]["sg_length"] = sg_length
        
        #write dictionary to crispr.yaml
        with open(yaml_file, "w") as f:
            yaml.dump(doc,f)
    else:
        #open file as dictionary and add library info
        with open(yaml_file) as f:
            doc = yaml.safe_load(f)
            doc[name] = {}
            for i in yaml_keys:
                doc[name][i] = ""
            doc[name]["index"] = index
            doc[name]["fasta"] = fasta
            doc[name]["csv"] = csv
            doc[name]["sg_length"] = sg_length
    
        #write appended dictionary with all sgRNA library info to crispr.yaml
        with open(yaml_file, "w") as f:
            yaml.dump(doc,f)

@click.command(name='analysis')
@click.option("--md5sums", is_flag=True, show_default=True, default=False,
              help="Check md5sums of fastq files")
@click.option("--fastqc", is_flag=True, show_default=True, default=False,
              help="Quality control of fastq files")
@click.option("-r", "--rename", is_flag=True, show_default=True, default=False, 
              help="Rename fastq files according to rename.csv")
@click.option("-t","--threads", default=1, type=int, 
              help="Number of CPU threads used during analysis")
@click.option("-l","--library",
              help="CRISPR-Cas9 library")
@click.option("-m","--mismatch", default=0, show_default=True, type=int, 
              help="Number of mismatches allowed during alignment")
@click.option("-a","--analysis", default="mageck", show_default=True, 
              type=click.Choice(["mageck","bagel2"]),
              help="Statistical analysis with MAGeCK or BAGEL2")

def analysis(md5sums,fastqc,rename,threads,library,mismatch,analysis):
    ''' Run CRISPR-Cas9 screen analysis pipeline
    '''
    #start run timer
    start = timeit.default_timer()
    
    click.secho("CRISPR-Cas9 screen analysis with pycrispr",fg="green")
    
    threads = str(threads)
    mismatch = str(mismatch)
    
    #create output dirs
    os.makedirs(os.path.join(work_dir,"count"),exist_ok=True)
    if analysis == "mageck":
        os.makedirs(os.path.join(work_dir,"mageck"),exist_ok=True)
    elif analysis == "bagel2":
        os.makedirs(os.path.join(work_dir,"mageck"),exist_ok=True)
    
    #run selected functions
        
    click.echo(f"{library} library selected")
    click.echo(f"Mismatches allowed: {mismatch}")
    
    if rename:
        utils.rename()
    
    if fastqc:
       utils.fastqc(threads)
    
    #align and count sgRNA reads
    utils.count(threads,mismatch,library)
    
    #join count files
    utils.join(library)
        
    #apply statistics to count files
    if analysis == "mageck":
        utils.mageck()
    elif analysis == "bagel2":
        utils.bagel2()
    
    #print total run time
    stop = timeit.default_timer()
    total_time = stop - start
    ty_res = time.gmtime(total_time)
    res = time.strftime("%H:%M:%S",ty_res)
    print(f"Total run time: {res}")
    
#add subparsers
cli.add_command(show_libs)
cli.add_command(add_lib)
cli.add_command(analysis)
cli.add_command(version)

    
        
    
    
