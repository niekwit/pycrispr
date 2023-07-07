Installation
====================================


The following steps are required to install ``pycrispr``:

   1. Installation Conda/Mamba (Mamba is highly recommended)
   2. Installation of ``snakemake``
   3. Installation of ``pycrispr``



1. Conda/Mamba installation
------------------------------------

``pycrispr`` requires the `snakemake <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_ workflow management software, it highly recommended to install the `mamba <https://mamba.readthedocs.io/en/latest/installation.html>`_ package manager first:

.. code-block:: console
   
   $ curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh" | bash

.. note:: Although highly recommended, using Conda/Mamba is not an absolute requirement. Using the ``pycrispr`` flags -c/--noconda, the Conda/Mamba requirement can be bypassed. However, it is then up to the user to install all the software and make sure they are set as environment variables.  



2. Installation of snakemake
------------------------------------------------------

Next ``snakemake`` can be installed with the following commands:

.. code-block:: console
   
   $ mamba activate base
   $ mamba create -c conda-forge -c bioconda -n snakemake snakemake

This will create a new virtual environment named ``snakemake``.



3. Installation of pycrispr
-----------------------------------------------------------

Latest development version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

From the ``snakemake`` virtual environment, run the following commands to install the latest development version of ``pycrispr``:

.. code-block:: console

   $  mamba activate snakemake
   $  cd /path/of/choice 
   $  git clone https://github.com/niekwit/pycrispr.git
   $  cd pycrispr
   $  pip install .

To update pycrispr to the latest version run:

.. code-block:: console

   $  cd /path/of/pycrispr 
   $  git pull
   $  pip install --upgrade .


Installation of stable version via Python Packaging Index (PyPi)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: console

   $  pip install pycrispr







