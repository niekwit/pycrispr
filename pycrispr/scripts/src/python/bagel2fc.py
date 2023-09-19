#!/usr/bin/env python

import os
import pandas as pd
import subprocess


b2dir = snakemake.input[0]
count_table_bagel2 = snakemake.input[1]
fc_table = snakemake.output[0]
log = snakemake.log[0]

#create dictionary to store the sample column numbers in count table
#(samples cannot be refered to by column name)
count_table = pd.read_csv(count_table_bagel2, sep="\t")
column_names = list(count_table.columns)
column_dict = {key: i for i, key in enumerate(column_names)}
column_dict = {key: column_dict[key] - 1 for key in column_dict} #first sample column should have value 1

#get comparison and control sample
comparison = os.path.basename(fc_table.replace(".foldchange",""))
control = comparison.split("_vs_")[1]
control_column = column_dict[control]

#generate foldchange command
command = f"python {b2dir}/BAGEL.py fc -i {count_table_bagel2} -o bagel2/{comparison}/{comparison} -c {control_column} 2> {log}"

#run command
print(f"BAGEL2: generating fold change table for {comparison}...")
subprocess.run(command, shell=True)




