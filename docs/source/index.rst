.. pycrispr documentation master file, created by
   sphinx-quickstart on Sat Jan 28 16:10:33 2023.
   
.. image:: https://readthedocs.org/projects/pycrispr/badge/?version=latest
    :target: https://pycrispr.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

pycrispr
====================================

**pycrispr**: a user-friendly Python library for analysing CRISPR-Cas9 screen data. 

Contents
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. toctree::
   :maxdepth: 2

   /requirements
   /installation
   /userguide
   /about
   Source code <https://github.com/niekwit/pycrispr/>

Quickstart
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Add sgRNA library information:

.. code-block:: console

   $ pycrispr add-lib --name yusa-mouse --index /path/to/hisat2-index --fasta none --csv none --sg-length 20

2. Run analysis:

.. code-block:: console

   $ pycrispr analysis --threads 4 -l yusa-mouse 

Links
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* `Source code <https://github.com/niekwit/pycrispr/>`_
* `Report an issue <https://github.com/niekwit/pycrispr/issues>`_
* `Project page on PyPi <xxx>`_


   
   
   
   
   
   
   
