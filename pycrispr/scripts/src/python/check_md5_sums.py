'''Check md5 sums of files in a directory against a file with md5 sums.
'''

import sys
import os
import hashlib
import pandas as pd
import glob
import click


def md5(file):
    
    hash_md5 = hashlib.md5()
    
    with open(file, "rb") as f:
        
        for chunk in iter(lambda: f.read(4096), b""):
            
            hash_md5.update(chunk)
    
    return(hash_md5.hexdigest())


def fetch_md5sums(files):
    
    if any([x for x in files if x.endswith("r_1.fq.gz")]): #Illumina platform data (NovaSeq/HiSeq)
        
        return glob.glob("*.md5sums.txt")
    
    else:
        
        return []
    

def exit():
        
        #create output file so not to crash snakemake if md5sum file does not exist
        with open(snakemake.output[0], "w") as myfile:
            myfile.write("No MD5sums calulated")
        
        #quietly exit 
        sys.exit(0)   
            
#get gz files 
files=glob.glob("*.gz")

#get md5sum files
md5sum_files = fetch_md5sums(files)

#if no md5sum files found, exit
if len(md5sum_files) == 0:
    
    exit()

#check any md5sum file is missing
for file in md5sum_files:
    
    if not os.path.exists(file):
        
        exit()

#df to store md5sums
df = pd.DataFrame(columns=["fastq","md5sum","md5sum_calculated","correct"],
                  index = range(0,len(md5sum_files)))

#add fastq files to df
df["fastq"] = fastq

#calculate md5sums and match to md5sums in md5sum_files
for i, file in enumerate(md5sum_files):
    
    #print(f"checking md5sum for {fastq[i]}...")
    #print(file)
    
    #open md5sum file
    md5sum = open(file, "r")
    md5sum = md5sum.readlines() #read1 md5sum can be on any line
    
    #select md5sum from line with read1
    md5sum = [x for x in md5sum if fastq[i] in x][0].replace("\n","").split(" ")[0]
    
    #add md5sum to df
    df.at[i,"md5sum"] = md5sum
    
    #add calculated md5sum
    df.at[i,"md5sum_calculated"] = md5(fastq[i])

#check if md5sums match
df["correct"] = df["md5sum"] == df["md5sum_calculated"]

#save df to file
df.to_csv(snakemake.output[0], sep="\t", index=False)

#check for different check sums and save to separate csv
df_fail = df[df["correct"] == False]

#print error message if md5sums do not match and raise error
if len(df_fail.axes[0]) > 0:
    
    click.secho(f"ERROR: md5sums do not match for the following files:", fg="red")
    click.echo("\n".join(df_fail["fastq"].tolist()))

    sys.exit(1)




