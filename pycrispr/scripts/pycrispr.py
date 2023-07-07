#!/usr/bin/env python3

import sys
import subprocess
import os
import shutil
import glob
import click
import yaml

script_dir, script_file = os.path.split(__file__)
work_dir = os.getcwd()


####Python functions####

def loadYaml():
    '''Load experiment.yaml as dictionary
    '''    
    with open("experiment.yaml") as f:
        doc = yaml.safe_load(f)
    return(doc)


def bg_msg(): #background message function
        click.echo("Pipeline will run in the background (nohup)")
        click.echo("Terminal output can be found in logs/terminal.log")


####command line parser####
@click.group()
def cli():
    
    '''Snakemake-based CRISPR-Cas9 screen analysis pipeline
    '''
        

@click.command(name='analysis')
@click.option("-t","--threads", 
              default=1, 
              show_default=True, 
              help="Total number of CPU threads to use for local analysis")
@click.option("-s","--slurm", 
              is_flag=True,
              help="Run pipeline on SLURM-based HPC")
@click.option("-b","--background", 
              is_flag=True,
              help="Run pipeline in the background with nohup (stdin/stderr will be directed to logs/terminal.log)")
@click.option("-c","--noconda", 
              is_flag=True,
              help="Disable jobs from running in a Conda environment (not recommended)")
@click.option("-d","--dryrun", 
              is_flag=True,
              help="Dry run for running pipeline (helpful for testing if pipeline works)")
@click.option("-r","--cleanup", 
              is_flag=True,
              help="Cleanup unused conda environments and packages")
@click.option("-v","--verbose", 
              show_default=True,
              is_flag=True,
              help="Increase verbosity")


def analysis(threads,slurm,background,noconda,dryrun,cleanup,verbose):
    
    ''' Run CRISPR-Cas9 screen analysis pipeline
    '''
    if cleanup:
        click.secho("Cleaning up unused Conda packages and environment...",fg="green")
        subprocess.run("snakemake --cores 1 --conda-cleanup-envs --conda-cleanup-pkgs cache", shell=True)
        click.secho("Done!",fg="green")
        sys.exit(0)
    
    click.secho("CRISPR-Cas9 screen analysis with pycrispr",fg="green")
    
    #total threads for local pipeline run
    threads = str(threads)    
   
    #copy snakemake file to work_dir:
    snakemake_file = os.path.join(script_dir,"src","snakefile")
    snakemake_copy = os.path.join(work_dir,"snakefile")
    shutil.copyfile(snakemake_file,snakemake_copy)
    
    #copy R script to work_dir
    os.makedirs(os.path.join(work_dir,"scripts"),exist_ok=True)
    flute_file = os.path.join(script_dir,"src","flute.R")
    flute_dest = os.path.join(work_dir,"scripts","flute.R")
    shutil.copyfile(flute_file,flute_dest)
    
    #copy Python scripts to work_dir
    python_files = glob.glob(os.path.join(script_dir,"src","python","*.py"))
    [shutil.copyfile(x,os.path.join(work_dir,"scripts",os.path.basename(x))) for x in python_files]
    
    #plot DAG
    if not os.path.exists("dag.pdf"):
        click.echo("Plotting snakemake DAG")
        dag = "snakemake --forceall --dag | grep -v '\-> 0\|0\[label = \"all\"' |dot -Tpdf > dag.pdf"
        process=subprocess.check_output(dag,shell=True)
        
    #construct snakemake command
    snakemake = "snakemake --output-wait 20" 
    if not noconda:
        snakemake = f"{snakemake} --use-conda" 
    if verbose:
        snakemake = f"{snakemake} -p" #prints shell commands
    if dryrun:
        click.echo("Dry run only")
        snakemake = f"{snakemake} -n"
    if slurm:
        click.echo("Submitting pipeline to Slurm workload manager")
        #load slurm default resources
        slurm = loadYaml()
        account = slurm["resources"]["account"]
        partition = slurm["resources"]["partition"]
        max_jobs = slurm["resources"]["max_jobs"]
        
        snakemake = f"{snakemake} --slurm -j {max_jobs} --default-resources slurm_account={account} slurm_partition={partition}" #run snakemake in background with nohup
    else:
        snakemake = f"{snakemake} --cores {threads}"
    
    #copy conda envs yamls to work_dir
    conda_envs = glob.glob(os.path.join(script_dir,"src","*.yaml"))
    os.makedirs("envs",exist_ok = True)
    [shutil.copyfile(x,os.path.join(work_dir,"envs",os.path.basename(x))) for x in conda_envs]
        
    #run snakemake command    
    if not os.path.isdir(".snakemake/"): #this dir does not exist before first run
        if background: #check if it needs to run in the background (nohup)
            os.makedirs("logs", exist_ok=True)
            snakemake = f"nohup {snakemake} >> logs/terminal.log &"
            bg_msg()
        subprocess.run(snakemake, shell=True)
    else: #if it has run before, it probably failed at some step so rerun all failed rules
        snakemake = f"{snakemake} --rerun-incomplete"
        if background:
            os.makedirs("logs", exist_ok=True)
            snakemake = f"nohup {snakemake} >> logs/terminal.log &"
            bg_msg()
        subprocess.run(snakemake, shell=True)


@click.command(name='report')

def report():
    
    ''' Create HTML report of analysis
    '''
    #create report for previous pipeline run
    click.secho("Creating report of analysis...",fg="green")
    
    #copy report files to work_dir
    shutil.copytree(os.path.join(script_dir,"src","report"),
                    os.path.join(work_dir,"report"), 
                    dirs_exist_ok=True)
        
    #create command
    report = "snakemake --report pycrispr-report.html"
    
    #run command
    subprocess.run(report, shell=True)
    return


@click.command(name='version')

def version():
    
    ''' Display version of pycrispr
    '''
    version = "1.0.1"
    
    #create report for previous pipeline run
    click.secho(f"pycrispr v{version}",fg="green")
    
    return
           
#add subparsers
cli.add_command(analysis)
cli.add_command(report)
cli.add_command(version)

    
        
    
    
