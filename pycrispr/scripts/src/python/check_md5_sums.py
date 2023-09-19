'''Check md5 sums of files in a directory against a file with md5 sums.
'''

import sys
import os
import hashlib
import pandas as pd
import glob
import click

#function to cacluate md5sum
def md5(file):
    
    hash_md5 = hashlib.md5()
    
    with open(file, "rb") as f:
        
        for chunk in iter(lambda: f.read(4096), b""):
            
            hash_md5.update(chunk)
    
    return(hash_md5.hexdigest())


def fetch_md5sums(files):
    
    match files:
        
        case any([x for x in files if x.endswith(".s_1.r_1.fq.gz")]):
            
            return glob.glob("*.md5sums.txt")
            
    


#get fastq and md5sum files




r1_tag = snakemake.params["r1_tag"]
fastq = glob.glob(f"reads/*{r1_tag}")
md5sum_files = [x.replace(r1_tag,".md5sums.txt") for x in fastq]

#check if md5sum files exist
for file in md5sum_files:
    
    if not os.path.exists(file):
        
        #create output file so not to crash snakemake if md5sum file does not exist
        with open(snakemake.output[0], "w") as myfile:
            myfile.write("No MD5sums calulated")
        
        #quietly exit 
        sys.exit(0)

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




