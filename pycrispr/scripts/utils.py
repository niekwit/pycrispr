#!/usr/bin/env python3

''' Functions for pycrispr
'''
import shutil
import sys
import os
import hashlib
import glob
import subprocess
import numpy as np
import pandas as pd
import click
import yaml
from itertools import compress

#set global variables
script_dir = os.path.abspath(os.path.dirname(__file__))
work_dir = os.getcwd()
stdout_log = os.path.join(work_dir,"stdout.log")


#Functions for pycrispr

def logCommandLineArgs():
    args = sys.argv
    args = " ".join(args)

    if "-h" not in args:
            if "--help" not in args:
                print(args,
                      file = open(os.path.join(work_dir,"commands.log"), "a"))


def md5sums(md5sum_files):
    """Checks md5sums of fastqc files
    """
    ###to do:add old file names to df so that this function can still be run after renaming files
    
    #get md5sum files
    #md5sum_files = glob.glob(os.path.join(work_dir,"raw-data","*.md5sums.txt"))
    if len(md5sum_files) == 0:
        click.secho("ERROR: no *.md5sum.txt files found in raw-data/",fg="red")
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
        click.secho("ERROR: md5sum of at least one file does not match\nCheck md5sums_failed.csv",fg="red")
        df_fail.to_csv(os.path.join(work_dir,"md5sums_failed.csv"),index=False)
        return(False)
    	

def rename():
    ''' Renames fastq files according to rename.csv
    '''
    click.echo("Renaming fastq files according to rename.csv")
    
    file = open(os.path.join(work_dir,"rename.csv"), "r")
    lines = file.readlines()
    count = 0
    for line in lines: #removes newline characters
        lines[count] = line.replace("\n","")
        count+=1

    for line in lines:#rename files
        old_name,new_name=line.split(";")
        os.rename(os.path.join(work_dir,"raw-data",old_name),
                  os.path.join(work_dir,"raw-data",new_name))


def checkInPath(check):
    ''' Check if commands in list are set in $PATH
    '''
    if type(check) != list: #single command (string)
        check = [check]
    
    #make sure everything is lower case
    check = [x.lower() for x in check]
    
    #get $PATH
    path = os.environ["PATH"].lower()
    
    #check if commands are in $PATH and if not show which ones
    if any([x in path for x in check]) == False:
        bool_list = [x in path for x in check]
        bool_list_rev = [not x for x in bool_list]
        not_in_path= list(compress(check, bool_list_rev))
        
        #not in $PATH but but check if in dir such as /home/user/bin that is in $PATH
        not_in_path_final = []
        for n in not_in_path:
            try:
                os.path.exists(shutil.which(n))
            except TypeError:
                not_in_path_final.append(n)
                
        if len(not_in_path_final) != 0:
            not_in_path_final = " & ".join(not_in_path_final)
            click.secho(f"ERROR: {not_in_path_final} not found in $PATH",fg="red")
            return(False)
        else:
            return(True)
    else:
        return(True)


def fastqc(threads):
    '''Quality control of fastq files
    '''
    click.echo("Quality control of fastq files using FastQC/MultiQC")
    
    #check if fastqc is in $PATH
    if not checkInPath("fastqc"):
        click.secho("ERROR: fastqc not found in $PATH",fg="red")
        return()
    
    #create fastqc dir
    fastqc_dir = os.path.join(work_dir,"fastqc")
    os.makedirs(fastqc_dir,exist_ok = True)
    
    #run fastqc
    data_files = glob.glob(os.path.join(work_dir,"raw-data","*.gz"))
    fastqc = ["fastqc","--threads",threads,"--quiet","-o",fastqc_dir]
    fastqc.extend(data_files)
    write2log(" ".join(fastqc))
    subprocess.run(fastqc)
    
    #run multiqc
    multiqc = ["multiqc","-o",fastqc_dir,fastqc_dir]
    write2log(" ".join(multiqc))
    subprocess.run(multiqc)


def write2log(command,name=""):
    '''Write bash command to commands.log
    '''
    with open(os.path.join(work_dir,"commands.log"), "a") as file:
        file.write(name)
        print(*command, sep = "", file = file)


def fileExists(file): #check if file exists/is not size zero
    '''Check if file exists or has 0 file size.
    '''    

    if os.path.exists(file):
            if os.path.getsize(file) > 0:
                print(f"Skipping {file} (already exists/analysed)")
                return(True)
    else:
        return(False)


def loadYaml(name):
    '''Load crispr.yaml as dictionary
    '''
    if name == "crispr"
        yaml_file = os.path.join(script_dir,"crispr.yaml")
    else:
        yaml_file = os.path.join(work_dir,f"{name}.yaml")
    
    with open(yaml_file) as f:
        doc = yaml.safe_load(f)
    
    return(doc)


def csv2fasta(csv):
    ''' Converts CSV file to fasta format
    '''
    df_CSV = pd.read_csv(csv)
    line_number_fasta = len(df_CSV) * 2

    df = pd.DataFrame(columns = ["column"],index = np.arange(line_number_fasta))

    #create fasta df
    df["column"] = df_CSV.stack().reset_index(drop = True)
    df.iloc[0::2, :] = ">"+df.iloc[0::2, :]
    '''
    library_name = os.path.basename(csv)
    library_name = library_name.replace(".csv","")
    fasta_base = library_name + ".fasta"
    fasta_file = os.path.join(script_dir,"index",library_name,fasta_base)
    os.makedirs(os.path.join(script_dir,"index",library_name),
                exist_ok = True)
    df.to_csv(fasta_file,index = False, 
              header = False)

    #add new CRISPR library and fasta file location to library.yaml
    yaml_list = ["clip_seq","fasta","index_path","read_mod","sg_length","species"]
    with open(os.path.join(script_dir, "crispr.yaml")) as f:
        doc = yaml.safe_load(f)
        doc[library_name] = {}
        for i in yaml_list:
            doc[library_name][i] = ""
        doc[library_name]["fasta"] = fasta_file
        with open(os.path.join(script_dir, "yaml" , "crispr-library.yaml"), "w") as f:
            yaml.dump(doc,f)

    #exit message
    sys.exit("Fasta file created and added to library.yaml\nPlease provide more CRISPR library information in this file before first run.")
    '''


def trim(threads,sg_length):
    ''' Remove vector sequence from reads
    '''
    files = glob.glob(os.path.join(work_dir,"raw-data","*.gz"))
    files = [x for x in files if not "trimmed" in x]
    click.echo("Trimming vector sequence from:")
    for file in tqdm(files, position = 0, leave = True):
        trimmed = file.split(".",1)[0] + "_trimmed.fq.gz"
        if not fileExists(trimmed):
            base_file = os.path.basename(file)
            tqdm.write(base_file)
            cutadapt = f"cutadapt -j {threads} --quiet --quality-base 33 -l {sg_length} -o {trimmed} {file} 2>> stdout.log"
            write2log(cutadapt)
            subprocess.run(cutadapt, shell = True)


def buildIndex():
    '''Builds HISAT2 index
    '''
    pass


def align(threads,mismatch,index):
    ''' Align and count trimmed reads with HISAT2
    '''
    os.makedirs(os.path.join(work_dir,"count"),exist_ok=True)
    
    files = glob.glob(os.path.join(work_dir,"raw-data","*_trimmed.fq.gz"))
    click.echo("Aligning samples:")
    for file in tqdm(files, position = 0, leave = True):
        base_file = os.path.basename(file)
        count_file = os.path.join(work_dir,"count",base_file.replace("_trimmed.fq.gz",".guidecounts.txt"))
        
        if not fileExists(count_file):
            '''
            1. Align with HISAT2 to index
            2. Select third field of SAM file that contains sgRNA name using awk
            3. Sort sgRNA names
            4. Count unique sgRNA names
            5. Remove leading white space using sed
            6. Remove line with unmapped reads (first line) (sed)
            '''
            tqdm.write(base_file)
            hisat2 = ["zcat",file,"|","hisat2","-p",threads,"-t","-N",mismatch,"-x",index, 
            "-","2>>",stdout_log,"|","awk","'{print $3}'","|","sort","|","uniq","-c","|",
            "sed","'s/^ *//'","|","sed","'1d",">",count_file]
            write2log(" ".join(hisat2))
            subprocess.run(f'echo "{base_file}" >> {stdout_log}',shell=True)
            subprocess.run(" ".join(hisat2),shell=True)
    
   
def count(threads,mismatch,library):
    
    #check if samtools and HISAT2 are in $PATH
    if not checkInPath(["hisat2","samtools"]):
        return()
    
    #load crispr.yaml
    doc = loadYaml()
    
    #load settings
    index = doc[library]["index"]
    fasta = doc[library]["fasta"]
    csv = doc[library]["csv"]
    sg_length = doc[library]["sg_length"]
    
    #check if index is available
    if index == "":
        #check if fasta is available to build HISAT2 index
        if fasta != "":
            click.secho(f"WARNING: no HISAT2 index found for {library} library, building now...",fg="red")
            index_dir = os.path.dirname(fasta)
            hisat2_build = f""
            
        elif csv != "":
            click.secho(f"WARNING: no HISAT2 index or fasta found for {library} library",fg="red")
            print("Generating fasta file from CSV file")
            csv2fasta(csv)
    
    else:
        #remove vector sequences from reads
        trim(threads,sg_length)
        
        #align and count reads
        align(threads,mismatch,index)


def join(library):
    ''' Join count files to create MAGeCK/BAGEL2 input
    '''
    click.echo("Creating count table for MAGeCK/BAGEL2")
    
    #load crispr.yaml
    doc = loadYaml()
    
    #get all count files
    file_list = glob.glob(os.path.join(work_dir,"count","*guidecounts.txt"))
    
    #dictionary to store all guide counts
    counts = {}
    
    for file in file_list:
       #add counts to counts dict
       df = pd.read_csv(file,sep=" ",header = None)
       key = os.path.basename(file).replace(".guidecounts.txt", "")
       df.columns = [key, "sgRNA"]
       df[key].astype(int)
       counts.update({key:df})
       
    #prepare data frame for left join
    fasta = doc[library]["fasta"]

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


def normalise():
    '''Normalise counts-aggregated.tsv to total read count per sample
    '''
    click.echo("Creating normalised count table")
    
    
    df = pd.read_table(os.path.join(work_dir,"count","counts-aggregated.tsv"))
    column_range = range(2,len(df.columns))
    for i in column_range:
        column_sum = df[df.columns[i]].sum()
        df[df.columns[i]] = df[df.columns[i]] / column_sum * 1E8
        
        df[df.columns[i]] = df[df.columns[i]].astype(int)
    df.to_csv(os.path.join(work_dir,"count","counts-aggregated-normalised.csv"),
              index = False,
              header = True)


def mageck():
    '''Statistical analysis of sgRNA counts with MAGeCK
    '''
    #check if MAGeCK is in $PATH
    if not checkInPath("mageck"):
        return()
    
    #check if stats.txt is available
    stats = os.path.join(work_dir,"stats.csv")
    if not os.path.exists(stats):
        click.secho("ERROR: stats.csv not found (MAGeCK comparisons)",fg="red")
        return()
    
    #load stats.txt
    df = pd.read_csv(os.path.join(work_dir,"stats.csv"))
    
    #run MAGeCK for each line in df
    click.echo("Running MAGeCK for:")
    for index,row in tqdm(df.iterrows(), position = 0, leave = True):
        samples = row.tolist()
        test = samples[0]
        control = samples [1]
        count_table = os.path.join(work_dir,"count","counts-aggregated.tsv")
        exp = f"{test}_vs_{control}"
        tqdm.write(exp.replace("_"," "))
        prefix=os.path.join(work_dir,"mageck",exp,exp)
        
        os.makedirs(os.path.join(work_dir,"mageck",exp), exist_ok=True)
        
        #create MAGeCK command
        mageck = f"mageck test -k {count_table} -t {test} -c {control} -n {prefix} 2>> {stdout_log}"
        
        #run command
        write2log(mageck)
        subprocess.run(mageck,shell=True)
        







