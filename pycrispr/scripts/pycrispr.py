#!/usr/bin/env python3

import subprocess
import os
import shutil
import glob
import click

from ..scripts import utils as utils


script_dir, script_file = os.path.split(__file__)
work_dir = os.getcwd()
version = "0.1"


####command line parser####
@click.group()
def cli():
    """Snakemake-based CRISPR-Cas9 screen analysis pipeline
    """
        

@click.command(name='report')

def report():
    ''' Create html report of data
    '''
    #create report for previous pipeline run
    click.secho("Creating report of analysis...",fg="green")
    #copy CU style sheet
    style_ori = os.path.join(script_dir,"src","cam_style.css")
    style = "src/cam_style.css"
    shutil.copyfile(style_ori,style)
    
    #create command
    report = f"snakemake --report report.html --report-stylesheet {style}"
    
    #run command
    subprocess.run(report, shell=True)
    return


@click.command(name='analysis')
@click.option("-t","--threads", 
              default=1, 
              show_default=True, 
              help="Total number of CPU threads to use for local analysis")
@click.option("-s","--slurm", 
              is_flag=True,
              help="Run pipeline on SLURM-based HPC")
@click.option("-d","--dryrun", 
              is_flag=True,
              help="Dry run for running pipeline (helpful for testing if pipeline works)")
@click.option("-v","--verbose", 
              show_default=True,
              is_flag=True,
              help="Increase verbosity")


def analysis(threads,slurm,dryrun,verbose):
    ''' Run CRISPR-Cas9 screen analysis pipeline
    '''
    click.secho("CRISPR-Cas9 screen analysis with pycrispr",fg="green")
    
    #total threads for local pipeline run
    threads = str(threads)    
    
    #check if files need to be renamed
    experiment = utils.loadYaml("experiment")
    if "rename" in experiment:
        utils.rename()
    
    
    #copy snakemake file to work_dir:
    snakemake_file = os.path.join(script_dir,"src","snakefile")
    snakemake_copy = os.path.join(work_dir,"snakefile")
    shutil.copyfile(snakemake_file,snakemake_copy)
    
    #copy utils to work_dir
    utils_file = os.path.join(script_dir,"utils.py")
    utils_copy = os.path.join(work_dir,"utils.py")
    shutil.copyfile(utils_file,utils_copy)
    
    #plot DAG
    click.echo("Plotting snakemake DAG")
    dag = "snakemake --forceall --dag | dot -Tpdf > dag.pdf"
    process=subprocess.check_output(dag,shell=True)
    
    
    #construct snakemake command
    snakemake = "snakemake --use-conda" 
    
    if verbose:
        snakemake = f"{snakemake} -p" #-p prints shell commands
    if dryrun:
        click.echo("Dry run only")
        snakemake = f"{snakemake} -n"
    if slurm:
        #load slurm default resources
        slurm = utils.loadYaml("slurm")
        account = slurm["account"]
        partition = slurm("partition")
        
        snakemake = f"{snakemake} --slurm --default-resources slurm_account={account} slurm_partition={partition}"
    else:
        snakemake = f"{snakemake} --cores {threads}"

    #copy conda envs yamls to work_dir
    conda_envs = glob.glob(os.path.join(script_dir,"src","*.yaml"))
    to_dir = os.path.join(work_dir,"envs")
    os.makedirs(to_dir,exist_ok = True)
    
    [shutil.copyfile(x,os.path.join(to_dir,os.path.basename(x))) for x in conda_envs]
    
    
    #run snakemake command
    subprocess.run(snakemake, shell=True)
    
       
#add subparsers
cli.add_command(report)
cli.add_command(analysis)

    
        
    
    
