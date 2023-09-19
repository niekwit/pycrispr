#!/usr/bin/env python

import os
import pandas as pd


fasta = snakemake.params["fa"]
count_table_mageck = snakemake.input[0]
count_table_bagel2 = snakemake.output[0]


###PREPARE COUNT TABLE FOR BAGEL2 FROM counts-aggregated.tsv

#load gRNA sequences and names from fasta
df_fasta = pd.read_csv(fasta, header=None)

#create new data frame: gRNA name and sequence in separate columns
df_name = df_fasta[df_fasta[0].str.contains(">")]
names = df_name.squeeze()#convert to series
names = names.reset_index(drop=True)#reset index
names = names.str.replace(">","")#remove >
names.name = "sgRNA"
df_seq = df_fasta[~df_fasta[0].str.contains(">")]
seq = df_seq.squeeze()#convert to series
seq = seq.reset_index(drop=True)#reset index
seq.name = "SEQUENCE"
df_join = pd.concat([names,seq],axis=1)#create df with names and sequences

#create gene column
df_join["gene"] = df_join["sgRNA"].str.split(pat="_").str[0]

#open MAGeCK formatted count table
df_master = pd.read_csv(count_table_mageck, sep="\t")

#format sgRNA column df_join to same format as df_master
df_join["sgRNA"] = df_join["sgRNA"].str.split(pat="_", n=1).str[1]

#sort data frames
df_join = df_join.sort_values(by="sgRNA")
df_master = df_master.sort_values(by="sgRNA")

#merge data frames and format data frame for BAGEL2
df_join = df_join.drop(["gene"], axis=1)#remove gene column
df_merge = pd.merge(df_master,df_join, on="sgRNA", how="left")
df_merge = df_merge.drop(["sgRNA"], axis=1)#remove gene column
df_merge.rename(columns={"gene":"GENE"}, inplace=True)
cols = list(df_merge)
cols.insert(0,cols.pop(cols.index("SEQUENCE")))
df_merge = df_merge.loc[:,cols]

#remove duplicate lines
df_merge.drop_duplicates(subset = "SEQUENCE", 
                       keep = "first", 
                       inplace = True)

#save df to file
print("BAGEL2: generating count table...")

os.makedirs("bagel2", exist_ok=True)
df_merge.to_csv(count_table_bagel2, sep='\t', index=False)




