#!/usr/bin/env python3

import subprocess
import os
import shutil
import click
import yaml

from ..scripts import utils as utils

#global variables
script_dir = os.path.abspath(os.path.dirname(__file__))
work_dir = os.getcwd()
version = "0.1"


####command line parser####
@click.group()
def cli():
    """Snakemake-based CRISPR-Cas9 screen analysis pipeline
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

@click.command(name='add-lib')
@click.option("-n",
              "--name", 
              required=True,
              help="Library name")
@click.option("-i","--index", 
              required=True,
              help="HISAT2 index path")
@click.option("-f","--fasta", 
              required=True,
              help="Fasta file path")
@click.option("-c","--csv", 
              required=True,
              help="CSV file path")
@click.option("--sg-length", 
              required=True, 
              type=int,
              help="sgRNA length")

def add_lib(name,index,fasta,csv,sg_length):
    ''' Add sgRNA library to crispr.yaml
    '''
    #create dictionary keys for yaml dump
    yaml_keys = ["fasta","index","csv","sg_length","species"]
    
    #write/create crispr.yaml file
    yaml_file = os.path.join(script_dir,"crispr.yaml")
    if not os.path.exists(yaml_file):
        #no pre-existing yaml file so start with dict to create one
        doc = {}
        doc["script_dir"] = script_dir 
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
@click.option("-t","--threads", 
              default=1, 
              show_default=True, 
              help="Total number of CPU threads to use for local analysis")
@click.option("-s","--slurm", 
              default=False, 
              show_default=True, 
              help="Run pipeline on SLURM-based HPC")
@click.option("-d","--dryrun", 
              default=False, 
              show_default=True,
              help="Dry run for running pipeline (helpful for testing if pipeline works)")
@click.option("-v","--verbose", 
              default=False, 
              show_default=True,
              help="Increase verbosity")


def analysis(threads,slurm,dryrun,verbose):
    ''' Run CRISPR-Cas9 screen analysis pipeline
    '''
    click.secho("CRISPR-Cas9 screen analysis with pycrispr",fg="green")
    
    #total threads for local pipeline run
    threads = str(threads)    
      
    #copy snakemake file to work_dir:
    snakemake_file = os.path.join(script_dir,"workflow","snakemake")
    snakemake_copy = os.path.join(work_dir,"snakemake")
    shutil.copyfile(snakemake_file,snakemake_copy)
    
    #plot DAG
    try:
        @click.echo("Plottig snakemake DAG")
        dag = "snakemake --forceall --dag | dot -Tpdf > dag.pdf"
        process=subprocess.check_output(dag,shell=True)
    except subprocess.CalledProcessError:
        print("ERROR: make sure snakemake and graphviz is installed")
    
    #construct snakemake command
    snakemake = "snakemake --use-conda" 
    
    if verbose:
        snakemake = f"{snakemake} -p" #-p prints shell commands
    if dryrun:
        @click.echo("Dry run only")
        snakemake = f"{snakemake} -n"
    if slurm:
        #load slurm default resources
        slurm = utils.loadYaml("slurm")
        account = slurm["account"]
        partition = slurm("partition")
        
        snakemake = f"{snakemake} --slurm --default-resources slurm_account={account} slurm_partition={partition}"
    else:
        snakemake = f"{snakemake} --cores {threads}"
    
    #run snakemake command
    subprocess.run(snakemake, shell = True)
    
       
#add subparsers
cli.add_command(show_libs)
cli.add_command(add_lib)
cli.add_command(analysis)
cli.add_command(version)

    
        
    
    
