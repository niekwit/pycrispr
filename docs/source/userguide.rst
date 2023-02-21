User guide
====================================
Show *pycrispr* help messages
------------------------------------
To display the ``pycrispr`` help message, use this command:

.. code-block:: console

   $ pycrispr --help
   Usage: pycrispr [OPTIONS] COMMAND [ARGS]...
   
   		CRISPR-Cas9 screen analysis pipeline
   
   Options:
     --help  Show this message and exit.
   
   Commands:
     add-lib    Add sgRNA library to crispr.yaml
     analysis   Run CRISPR-Cas9 screen analysis pipeline
     md5sums	Check md5sums of fastq files
     show-libs	Display which sgRNA library info has been stored
     version    Report the current build and version number

This shows all available functionalities, which also have their own help messages, for example:

.. code-block:: console

   $ pycrispr analysis --help
   Usage: pycrispr analysis [OPTIONS]
   
   		Run CRISPR-Cas9 screen analysis pipeline
   
   Options:
     --md5sums                       Check md5sums of fastq files
     --fastqc                        Quality control of fastq files
     -r, --rename                    Rename fastq files according to rename.csv
     -t, --threads INTEGER           Number of CPU threads used during analysis
     -l, --library TEXT              CRISPR-Cas9 library
     -m, --mismatch INTEGER          Number of mismatches allowed during alignment [default: 0]
     -a, --analysis [mageck|bagel2]  Statistical analysis with MAGeCK or BAGEL2 [default: mageck]


   
Add sgRNA library information
------------------------------------
Before using ``pycrispr``, information about the sgRNA library needs to be added:

1. sgRNA library name
2. Path to HISAT2 index (if available, otherwise none)
3. Fasta file that contains the sgRNA names and sequences (if available, otherwise none)
4. CSV file with two columns: sgRNA name and sgRNA sequence

This can be done with the following command:

.. code-block:: console

   $ pycrispr add-lib --name yusa-mouse --index none --fasta /path/to/fasta.fa --csv none --sg-length 20 
   
This will add an sgRNA library called `yusa-mouse <https://www.addgene.org/pooled-library/yusa-crispr-knockout-mouse-v2/>`_, with no HISAT2 index available, a path to the sgRNA library fasta file, and no path to a CSV file that contains sgRNA names and sequences. ``pycrispr`` can now be used to analyse data with this sgRNA library. If no path to a pre-existing HISAT2 index is given, it will be build from the fasta file, or CSV file, which is first converted to a fasta file. The index path will be stored automatically. Adding an sgRNA library will only have to be done once.

Show saved sgRNA library information
------------------------------------
To see which sgRNA libraries have been added to ``pycrispr``, run:

.. code-block:: console

   $ pycrispr show-libs
   
   The following sgRNA libraries are available:
   yusa-mouse:
     csv: none
     fasta: /home/user/Documents/fasta/yusa-mouse.fasta
     index: none
     sg_length: 20

Preparing CRISPR-Cas9 screen data
------------------------------------
Before running ``pycrispr`` an analysis directory has to be created (can be any name or location), and should contain a sub-directory called raw-data. This sub-directory contains all the fastq files from your CRISPR-Cas9 screen experiment::

    analysis_dir
    └── raw-data
    	├── SLXXXXXX1_R1_001.fq
    	├── SLXXXXXX2_R1_001.fq
    	└── SLXXXXXX3_R1_001.fq


.. important::
	Please note that ``pycrispr`` only accepts single-end NGS data, so if your data was sequenced in a paried-end fashion, only include the mate that contains the sgRNA sequence information (most commonly read 1). It also assumes that the first nucleotide sequenced is the first nulceotide of the sgRNA sequence.

Preparing configuration files
------------------------------------
rename.csv (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Depending on the NGS platform, fastq files can have very long file names, and as ``pycrispr`` uses the basename of a file as its sample name, it is advised to rename your fastq files prior to analysis. The existing and new files names can be included in a csv file called ``rename.csv`` as follows::

	existing,new
	SLXXXXXX1_R1_001.fq,S1.fq
	SLXXXXXX2_R1_001.fq,S2.fq
	SLXXXXXX3_R1_001.fq,L1.fq

How to apply this file will be descibed below.

stats.csv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If statistical analysis of sgRNA counts is required, a ``stats.csv`` file is needed with the following content::

	test,control
	S1,L1
	S2,L1
	S1;S2,L1


This will run MAGeCK RRA or BAGEL2 for three pair-wise comparisons:

1. S1 (test) vs L1 (control)
2. S2 (test) vs L1 (control)
3. S1,S2 (combined test samples) vs L1 (control)


As shown in comparison 3, multiple sample can be combined by separating them with a semi-colon. 

The ``rename.csv`` and ``stats.csv`` files should be locatated in the main analysis directory::

	analysis_dir
	├── raw-data
	├── stats.csv
	└── rename.csv

Checking md5sums fastq files
------------------------------------




Analysing CRISPR-Cas9 screen data
------------------------------------
The options for the CRISPR-Cas9 screen analysis are as follows::

	  --md5sums                       Check md5sums of fastq files
	  --fastqc                        Quality control of fastq files
	  -r, --rename                    Rename fastq files according to rename.txt
	  -t, --threads INTEGER           Number of CPU threads used during analysis
	  -l, --library TEXT              CRISPR-Cas9 library
	  -m, --mismatch INTEGER          Number of mismatches allowed during alignment [default: 0]
	  -a, --analysis [mageck|bagel2]  Statistical analysis with MAGeCK or BAGEL2 [default: mageck]

Analysis can be initiated with the following command:

.. code-block:: console

   $ pycrispr pycrispr analysis -r -t 4 -l yusa-mouse 

This will first rename the files according to ``rename.csv``, use four CPU threads, select the yusa-mouse sgRNA library added earlier, and use MAGeCK for pair-wise comparisons between samples according to ``stats.csv``. 

As the HISAT2 index path was not provided, ``pycrispr`` will build this first and include the path with the rest of the information. The index will be stored in a new directory called HISAT2-index in the parent directory of the fasta/csv file. This can be changed by providing an alternative path at the confirmation prompt:




Output files
------------------------------------
























