#!/usr/bin/env python

import os
import subprocess

b2dir = snakemake.params["b2dir"]
species = snakemake.params["species"]
bf = snakemake.input[0]
pr = snakemake.output[0]
log = snakemake.log[0]


#load gene sets
if species == "hsa":
    eg = f"{b2dir}/CEGv2.txt" #essential genes
    neg = f"{b2dir}/NEGv1.txt" #non-essential genes
elif species == "mmu":
    eg = f"{b2dir}/CEG_mouse.txt" #essential genes
    neg = f"{b2dir}/NEG_mouse.txt" #non-essential genes

#get comparison and control sample
comparison = os.path.basename(bf.replace(".bf",""))

#generate precission-recall command
command = f"python {b2dir}/BAGEL.py pr -i {bf} -o {pr} -e {eg} -n {neg} 2> {log}"

#run command
print(f"BAGEL2: calculating precision-recall for {comparison}...")
subprocess.run(command, shell=True)





