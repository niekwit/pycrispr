#!/usr/bin/env python3

import click
from ..scripts import utils as utils


#command line parser
@click.command()
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
@click.option("--load-library", is_flag=True, show_default=True, default=False,
              help="Add sgRNA library to crispr.yaml")


def cli(md5sums,fastqc,rename,threads,library,mismatch,analysis,cnv,fdr,go,load_library):
    """CRISPR-Cas9 screen analysis"""
    click.secho("CRISPR-Cas9 screen analysis with pycrispr",fg="green")
        
    #add sgRNA library info
    if load_library == True:
        utils.loadLibrary()
        click.echo("Please run pycrispr again (sgRNA library info added)")
        return()
    
    #run selected functions
    if md5sums == True:
        click.echo("Checking md5sums of fastq files in raw-data/")
        md5sum_match = utils.md5sums()
        if md5sum_match == False:
            click.echo("ERROR: At least one calculated md5sum did not match the pre-calculated ones\nPlease check md5sums_failed.csv",color="red")
    click.echo(f"{library} library selected")
    click.echo(f"Mismatches allowed: {mismatch}")
    if rename == True:
        click.echo("Renaming fastq files according to rename.txt")
        utils.rename()
    if fastqc == True:
        utils.fastqc(threads)
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
        
        
#if __name__ == '__main__':
#    cli()
        
        