#!/usr/bin/env python3

''' Functions for pycrispr
'''

import sys
import os
import hashlib
import glob
import subprocess
import pandas as pd
from clint.textui import colored, puts
import yaml
from itertools import compress


def loadLibrary():
    ''' Add sgRNA library to crispr.yaml file
    '''
    print("Adding new sgRNA library to library.yaml")
    script_dir = os.path.abspath(os.path.dirname(__file__))
    
    #ask user for sgRNA library info
    name = input("sgRNA library name? ")
    index = input("sgRNA library HISAT2 index path (if unavailable leave blank and enter fasta file path in next prompt)? ")
    fasta = input("sgRNA library fasta file path (this will be used to create the HISAT2 index. If not available leave blank and enter the path to a csv file with sgRNA name and sequence in two separate column)? ")
    csv = input("CSV file path (sgRNA name and sequence in separate columns)? ")
    #puts(colored.red("ERROR: insufficient information for sgRNA library"))
    sg_length = input("sgRNA length (if sgRNA library contains sgRNAs with varying lengths, enter the shortest length)? ")
    species = input("Species (required for gene ontology analysis)? ")
    
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
        doc[name]["species"] = species
        
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
            doc[name]["species"] = species
    
        #write appended dictionary with all sgRNA library info to crispr.yaml
        with open(yaml_file, "w") as f:
            yaml.dump(doc,f)
    

def md5sums(work_dir):
    """Checks md5sums of fastqc files
    """
    ###to do:add old file names to df so that this function can still be run after renaming files
    
    #get md5sum files
    md5sum_files = glob.glob(os.path.join(work_dir,"raw-data","*.md5sums.txt"))
    if len(md5sum_files) == 0:
        puts(colored.red("ERROR: no *.md5sum.txt files found in raw-data/"))
        return()
    
    #df to store md5sums
    df = pd.DataFrame(columns=["fastq","md5sum","md5sum_calculated","correct"],
                      index = range(0,len(md5sum_files)*2))
    
    #function to cacluate md5sum
    def md5(work_dir,file):
        file = os.path.join(work_dir, "raw-data", file)
        hash_md5 = hashlib.md5()
        with open(file, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return(hash_md5.hexdigest())
    
    #for index,row in df.iterrows():
    for i in range(0,len(md5sum_files)):
        #read file
        file = open(os.path.join(work_dir,"raw-data",md5sum_files[i]), "r")
        lines = file.readlines()
            
        for l in range(0,len(lines)):
            #add all fastq files and their md5sums to df
            md5sum = lines[l].split("  ")[0]
            df.at[i*2+l,"md5sum"] = md5sum
            
            fastq = lines[l].split("  ")[1].replace("\n","")
            df.at[i*2+l,"fastq"] = fastq
    
            #calculate md5sums of actual fastq files and add to df
            md5sum_calc = md5(work_dir,fastq)
            df.at[i*2+l,"md5sum_calculated"] = md5sum_calc
            
    #compare original checksums with calculated ones
    df["correct"] = df["md5sum"] == df["md5sum_calculated"]
    
    #save df to csv
    df.to_csv(os.path.join(work_dir,"raw-data","md5sums_checked.csv"),index=False)      
         
    #check for different check sums and save to separate csv
    df_fail = df[df["correct"] == False]
    
    if len(df_fail.axes[0]) > 0:
        puts(colored.red("ERROR: md5sum of at least one file does not match\nCheck md5sums_failed.csv"))
        df_fail.to_csv(os.path.join(work_dir,"md5sums_failed.csv"),index=False)
        return(False)
    	

def rename(work_dir):
    file = open(os.path.join(work_dir,"rename.txt"), "r")
    lines = file.readlines()
    count = 0
    for line in lines: #removes newline characters
        lines[count] = line.replace("\n","")
        count+=1

    for line in lines:#rename files
        old_name,new_name=line.split(";")
        os.rename(os.path.join(work_dir,"raw-data",old_name),
                  os.path.join(work_dir,"raw-data",new_name))


def fastqc(work_dir,threads):
    '''Quality control of fastq files
    '''
    #check if fastqc is in $PATH
    path = os.environ["PATH"].lower()
    if "fastqc" not in path:
        puts(colored.red("ERROR: fastqc not found in $PATH"))
        return()
    
    #create fastqc dir
    fastqc_dir = os.path.join(work_dir,"fastqc")
    os.makedirs(fastqc_dir,exist_ok = True)
    
    #run fastqc
    data_files = os.path.join(work_dir,"raw-data","*.gz")
    fastqc = f"fastqc --threads {str(threads)} --quiet -o {fastqc_dir} {data_files}"
    write2log(work_dir,fastqc)
    subprocess.run(fastqc, shell = True)


def logCommandLineArgs(work_dir):
    args = sys.argv
    args = " ".join(args)

    if "--help" not in args:
        print(args,file = open(os.path.join(work_dir,"commands.log"), "a"))


def write2log(work_dir,command,name=""):
    '''Write bash command to commands.log
    '''
    with open(os.path.join(work_dir,"commands.log"), "a") as file:
        file.write(name)
        print(*command, sep = "", file = file)


def file_exists(file): #check if file exists/is not size zero
    '''Check if file exists or has 0 file size.
    '''    

    if os.path.exists(file):
            if os.path.getsize(file) > 0:
                print(f"Skipping {file} (already exists/analysed)")
                return(True)
    else:
        return(False)


def loadYaml():
    '''Load crispr.yaml as dictionary
    '''
    script_dir = os.path.abspath(os.path.dirname(__file__))
    yaml_file = os.path.join(script_dir,"crispr.yaml")
    
    with open(yaml_file) as f:
        doc = yaml.safe_load(f)
    
    return(doc)

def count(work_dir,threads,mismatch,library):
    #check if samtools and HISAT2 are in $PATH
    check = ["samtools","hisat2"]
    path = os.environ["PATH"].lower()
    if any([x in path for x in check]) == False:
        bool_list = [x in path for x in check]
        bool_list_rev = [not x for x in bool_list]
        not_in_path= " ".join((list(compress(check, bool_list_rev))))
        puts(colored.red(f"ERROR: {not_in_path} not found in $PATH"))
        return()
    
    #load crispr.yaml
    doc = loadYaml()


def join(work_dir,yaml,library):
    file_list = glob.glob(os.path.join(work_dir,"count","*guidecounts.txt"))
    
    #dictionary to store all guide counts
    counts = {}
    
    for file in file_list:
       #add counts to counts dict
       df = pd.read_csv(file,
                        sep=" ",
                        header = None)
       key = os.path.basename(file).replace(".guidecounts.txt", "")
       df.columns = [key, "sgRNA"]
       df[key].astype(int)
       counts.update({key:df})
       
    #prepare data frame for left join
    fasta = yaml[library]["fasta"]

    df = pd.read_csv(fasta, header = None)
    df.columns = ["sgRNA"]
    
    df = df[df["sgRNA"].str.contains(">")]
    df = df.reset_index(drop = True)
    df["sgRNA"] = df["sgRNA"].str.replace(">", "")
    
    df["gene"] = df["sgRNA"].str.split(pat = "_",
                                       n = 1,
                                       expand = True)[0]
    
    #perform left join
    for key, value in counts.items():
        df = pd.merge(df, value, on='sgRNA', how='left')
        
    df["sgRNA"] = df["sgRNA"].str.split(pat = "_",
                                   n = 1,
                                   expand = True)[1]
    
    #replace nan with zero
    df = df.fillna(0)
    df = df.sort_values(by = ["sgRNA"])
    
    #convert floats back to int after pandas merge (bug in pandas)
    index_range = range(2, len(df.columns))
    index_list = []
    for i in index_range:
        index_list.append(i)
        
    df.iloc[:, index_list] = df.iloc[:, index_list].astype(int)
    
    #save data frame to file
    df.to_csv(os.path.join(work_dir,"count",'counts-aggregated.tsv'), 
              sep = '\t',
              index = False)


def normalise(work_dir):
    '''Normalise counts-aggregated.tsv to total read count per sample
    '''
    
    df = pd.read_table(os.path.join(work_dir,"count","counts-aggregated.tsv"))
    column_range = range(2,len(df.columns))
    for i in column_range:
        column_sum = df.iloc[:,i].sum()
        df.iloc[:,i] = df.iloc[:,i] / column_sum * 1E8
        df.iloc[:,i] = df.iloc[:,i].astype(int)
    df.to_csv(os.path.join(work_dir,"count","counts-aggregated-normalised.csv"),
              index = False,
              header = True)


def mageck(work_dir):
    pass


def bagel2(work_dir):
    pass


def go(work_dir):
    pass







