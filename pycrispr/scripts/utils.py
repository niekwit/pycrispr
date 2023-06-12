#!/usr/bin/env python3

''' Functions for pycrispr
'''
import os
import click
import yaml
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('agg')

#set global variables
script_dir = os.path.abspath(os.path.dirname(__file__))
work_dir = os.getcwd()
stdout_log = os.path.join(work_dir,"stdout.log")


#Functions for pycrispr

def loadYaml(name):
    '''Load yaml as dictionary'''    
    with open(f"{name}.yaml") as f:
        doc = yaml.safe_load(f)
    return(doc)


def rename():
    ''' Renames fastq files according experiment.yaml
    '''
    click.echo("Renaming fastq files according to experiment.yaml")
    
    yaml = loadYaml("experiment")
    
    rename = yaml["rename"]

    old_names = list(rename.keys())
    new_names = yaml["rename"]
    
    for old,new in zip(old_names,new_names):
        os.rename(os.path.join(work_dir,"reads",old),
                  os.path.join(work_dir,"reads",new))
    

def join(fasta,file_list):
    ''' Join count files to create MAGeCK/BAGEL2 input
    '''
    counts = {}
       
    for file in file_list:
       #add counts to counts dict
       df = pd.read_csv(file,sep=" ",header = None)
       key = os.path.basename(file).replace(".guidecounts.txt", "")
       df.columns = [key, "sgRNA"]
       df[key].astype(int)
       counts.update({key:df})
  
    df = pd.read_csv(fasta, header = None)
    df.columns = ["sgRNA"]
    
    df = df[df["sgRNA"].str.contains(">")]
    df = df.reset_index(drop = True)
    df["sgRNA"] = df["sgRNA"].str.replace(">", "")
    df["gene"] = df["sgRNA"].str.split(pat = "_",n = 1,expand = True)[0]
    
    #perform left join on count files
    for key, value in counts.items():
        df = pd.merge(df, value, on='sgRNA', how='left')
    df["sgRNA"] = df["sgRNA"].str.split(pat = "_",n = 1,expand = True)[1]
    
    #replace nan with zero
    df = df.fillna(0)
    df = df.sort_values(by = ["sgRNA"])
    
    #convert floats back to int after pandas merge (bug in pandas)
    index_range = range(2, len(df.columns))
    index_list = []
    for i in index_range:
        index_list.append(i)
    df[df.columns[index_list]] = df[df.columns[index_list]].astype(int)
    
    #save data frame to file
    df.to_csv(os.path.join(work_dir,"count",'counts-aggregated.tsv'), 
              sep = '\t',
              index = False)


def plot(df,y_label,save_file):
    '''General plotting function'''
    sns.set_style("white")
    sns.set_style("ticks")
    sns.barplot(x=list(df.keys())[0],
                    y=list(df.keys())[1],
                    data=df,
                    color="royalblue",
                    edgecolor="black",
                    linewidth=1)
    plt.ylabel(y_label)
    plt.xticks(rotation = 'vertical')
    plt.xlabel("")
    plt.tight_layout()
    sns.despine()
    plt.savefig(save_file)
    plt.close()


def plot_alignment_rate(log_files):
    plot_file="count/alignment-rates.pdf"
    
    #create df to store alignment rates
    df = pd.DataFrame(columns=["sample","alignment_rate"],index=np.arange(len(log_files)))
    samples = []
    rates = []

    for i in sorted(log_files):
        #get sample name from file name
        sample = os.path.basename(i).replace(".log","")
        samples.append(sample)
        
        #extract alignment rate from file
        with open(i) as f:
            lines = f.readlines()
        rate = float([x for x in lines if "overall alignment rate" in x][0].replace("% overall alignment rate\n",""))
        rates.append(rate)

    #add values to empty df for plotting
    df["sample"] = samples
    df["alignment_rate"] = rates

    #plot alignment rate
    plot(df,"Overall alignment rate (%)",plot_file)


def plot_coverage(fasta,count_table): #plots coverage per sample after alignment
    plot_file = "count/sequence-coverage.pdf"

    #get number of sgRNAs in CRISPR library
    fasta = pd.read_table(fasta, header = None)
    lib_size = len(fasta) / 2

    #extract number of single mapped aligned reads from counts-aggregated.tsv
    df = pd.read_table("count/counts-aggregated.tsv", sep = "\t")
    column_names = list(df.columns)
    del column_names[0:2] #delete sgRNA and gene columns
    
    counts = {}
    for i in column_names:
        count_sum = []
        count_sum.append(df[i].sum())
        counts[i] = count_sum
    
    #convert counts to sequence coverage
    df = pd.DataFrame(counts)
    df = df / lib_size
    
    #order columns alphabetically
    df = df.reindex(sorted(df.columns), axis=1)
    
    #transpose data frame
    df = df.transpose()
    df["samples"] = df.index
    names = ["coverage", "samples"]
    df.columns = names
    df = df[["samples","coverage"]]
    
    #plot coverage per sample
    plot(df,"Fold sequence coverage per sample",plot_file)









