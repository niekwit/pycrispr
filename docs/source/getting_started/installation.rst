Installation
====================================

Mamba installation
------------------------------------

As ``pycrispr`` is a library that is based on the `snakemake <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_ workflow management software, it is required to install the `mamba <https://mamba.readthedocs.io/en/latest/installation.html>`_ package manager first:

.. code-block:: console
   
   $ curl micro.mamba.pm/install.sh | bash


Install ``snakemake`` into its own virtual environment 
------------------------------------------------------

Next ``snakemake`` can be installed with the following commands:

.. code-block:: console
   
   $ conda activate base
   $ mamba create -c conda-forge -c bioconda -n snakemake snakemake

This will create a new virtual environment named ``snakemake``.

Installation of latest development version of ``pycrispr``
-----------------------------------------------------------
From the ``snakemake`` virtual environment, run the following commands to install the latest development version of ``pycrispr``:

.. code-block:: console

   $  mamba activate snakemake
   $  cd /path/of/choice 
   $  git clone https://github.com/niekwit/pycrispr.git
   $  cd pycrispr
   $  pip install .

To update ``pycrispr`` to the latest version run:

.. code-block:: console

   $  cd /path/of/pycrispr 
   $  git pull
   $  pip install --upgrade .


Installation of ``pycrispr`` via Python Packaging Index (PyPi)
----------------------------------------------------------------

.. attention:: This is not available yet.