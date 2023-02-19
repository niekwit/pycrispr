.. pycrispr documentation master file, created by
   sphinx-quickstart on Sat Jan 28 16:10:33 2023.
   
pycrispr
====================================

**pycrispr**: a user-friendly Python library for analysing CRISPR-Cas9 screen data. ``pycrispr`` is being developped by Dr. Niek Wit at the University of Cambridge in the groups of `Prof James Nathan <https://www.jamesnathanlab.com>`_ and `Prof Paul Lehner <https://www.citiid.cam.ac.uk/paul-lehner/>`_. Please cite ``pycrispr`` if you found it useful and use it in a publication.

Index
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. toctree::
   :maxdepth: 2

   /requirements
   /installation
   /userguide
   /cite
   Source code <https://github.com/niekwit/pycrispr/>

Quickstart
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Add sgRNA library information:

.. code-block:: console

   $ pycrispr add-lib --name yusa-mouse --index none --fasta /path/to/fasta.fa --csv none --sg-length 20

2. Run analysis:

.. code-block:: console

   $ pycrispr analysis --threads 4 -l yusa-mouse 

Links
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* `Source code <https://github.com/niekwit/pycrispr/>`_
* `Report an issue <https://github.com/niekwit/pycrispr/issues>`_
* `Project page on PyPi <xxx>`_

.. note::

   This project is under active development.
   
   
   
   
   
   
   
