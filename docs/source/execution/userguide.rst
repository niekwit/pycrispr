User guide
************

Show *pycrispr* help messages
------------------------------------
To display the ``pycrispr`` help message, use this command:

.. code-block:: console
   
   $ pycrispr --help
   Usage: pycrispr [OPTIONS] COMMAND [ARGS]...

   Snakemake-based CRISPR-Cas9 screen analysis pipeline

   Options:
      --help  Show this message and exit.

   Commands:
      analysis  Run CRISPR-Cas9 screen analysis pipeline
      report    Create HTML report of analysis

This shows all available functionalities, which also have their own help messages, for example:

.. code-block:: console

   $ pycrispr analysis --help
   Usage: pycrispr analysis [OPTIONS]

      Run CRISPR-Cas9 screen analysis pipeline

   Options:
      -t, --threads INTEGER  Total number of CPU threads to use for local analysis
                              [default: 1]
      -s, --slurm            Run pipeline on SLURM-based HPC
      -d, --dryrun           Dry run for running pipeline (helpful for testing if
                              pipeline works)
      -v, --verbose          Increase verbosity
      --help                 Show this message and exit.



Getting started with ``pycrispr``
------------------------------------
``pycrispr`` requires a YAML file (experiment.yaml) that contains information of the experiment, and available CRISPR sgRNA libraries:

.. code-block:: yaml

   slurm: False #submit jobs to SLURM-based HPC
   rename: #rename to .fq.gz
      L8_S1_L001_R1_001.fastq.gz: L8.fq.gz
      S8_S3_L001_R1_001.fastq.gz: S8.fq.gz
      S15_S4_L001_R1_001.fastq.gz: S15.fq.gz
   library: dub-only #CRISPR library for current experiment
   lib_info: #all available CRISPR libraries
      dub-only:
         index: /home/user/Documents/references/index/hisat2/dub-only/dub-only-hisat2.index #HISAT2 index path
         fasta: /home/user/Documents/references/fasta/Human/dub-only/DUBonly.fasta
         sg_length: 20
         species: hsa #hsa for human, mmu for mouse
   mismatch: 0 #number of mismatches allowed during alignment with HISAT2
   stats: 
      type: mageck
      comparisons: #test vs control 
         1: S8_vs_L8 #sample names are file names without extension
         2: S15_vs_L8
         3: S8,S15_vs_L8 # samples can be pooled together
   resources:
      account: SLURM_ACCOUNT_NAME #only for HPC use
      partition: PARTITION #only for HPC use
      short:
         cpu: 1
         time: 15 # in minutes; only for HPC use
      trim:
         cpu: 4
         time: 60
      fastqc:
         cpu: 4
         time: 60
      count:
         cpu: 8
         time: 120
      mageck:
         cpu: 1
         time: 60

.. note:: You can delete the rename section if you do not need to rename your files, but please keep in mind that the sample names will be taken from the read files names by removing the file extension. Also, the *comparisons* in the *stats* section should match this.


Preparing CRISPR-Cas9 screen data
------------------------------------
Before running ``pycrispr`` an analysis directory has to be created (can be any name or location), and should contain a sub-directory called *reads*. This sub-directory contains all the fastq files of your CRISPR-Cas9 screen experiment::

   analysis_dir
   └── reads
    	├── L8_S1_L001_R1_001.fastq.gz
    	├── S8_S3_L001_R1_001.fastq.gz
    	└── S15_S4_L001_R1_001.fastq.gz
   └── experiment.yaml 


.. important:: Please note that ``pycrispr`` only accepts single-end NGS data, so if your data was sequenced in a paried-end fashion, only include the mate that contains the sgRNA sequence information (most commonly read 1). It also assumes that the first nucleotide sequenced is the first nulceotide of the sgRNA sequence.


Initiating the pipeline
------------------------------------
To start the analysis, run:

.. code-block:: console

   $ pycrispr analysis -t 24

This will first rename the files according to *experiment.yaml*, use a total of 24 CPU threads, select the *dub-only* sgRNA library, and use MAGeCK for pair-wise comparisons specified in *experiment.yaml*. 


Output files
------------------------------------

Multiple output files will be generated::

   analysis_dir
   └── count
   |   ├── alignment-rates.pdf
   |   ├── counts-aggregated.tsv
   |   ├── L8.guidecounts.txt
   |   ├── S15.guidecounts.txt
   |   ├── S8.guidecounts.txt
   |   └── sequence-coverage.pdf
   └── envs
   |   ├── count.yaml
   |   ├── flute.yaml
   |   ├── join.yaml
   |   ├── mageck.yaml
   |   └── trim.yaml
   └── logs
   |   ├── count
   |   ├── fastqc
   |   ├── mageck
   |   ├── multiqc
   |   └── trim
   └── mageck
   └── mageck_flute
   └── qc
   └── reads
   | 	├── L8.fq.gz
   | 	├── S8.fq.gz
   | 	└── S15.fq.gz
   └── scripts
   |   └── flute.R
   ├── dag.pdf
   ├── experiment.yaml
   ├── snakefile
   └── utils.py



.. figure:: dag.png
   :align: center

   Directed acyclic graph (DAG) for workflow



   