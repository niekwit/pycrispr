Installation
====================================
Installation with pip
------------------------------------

To install **pycrispr** with pip:

.. code-block:: console

   $ pip install pycrispr
   
This command installs **pycrispr** and all the Python dependencies. 

.. important::

    It is recommended to install **pycrispr** in its own `virtual environment <https://docs.python.org/3/library/venv.html>`_ .
    
    **pycrispr** also requires non-Python dependencies that will need to be installed by the user:
    
    * `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_ (optional)
    * `HISAT2 <http://daehwankimlab.github.io/hisat2/>`_
    * `samtools <https://www.htslib.org>`_
    * `MAGeCK (optional) <https://sourceforge.net/p/mageck/wiki/Home/>`_
    * `BAGEL2 (optional) <https://github.com/hart-lab/bagel>`_
    
    These dependencies should be set in $PATH or in a directory that is in $PATH, such as /usr/local/bin. For instructions to add a directory to your $PATH, click `here <https://stackoverflow.com/questions/14637979/how-to-permanently-set-path-on-linux-unix>`_

Statistical analysis can be performed with either MAGeCK or BAGEL2. If these are not available then ``pycrispr`` will only generate sgRNA count tables, and will perform no statistical analysis.

Installation with conda
------------------------------------

To install ``pycrispr`` with conda:

.. code-block:: console

   $ conda install pycrispr
   
The non-Python dependencies will be included in the conda installation.


.. note::

   This is not yet available.











