.. scATS documentation master file, created by
   sphinx-quickstart on Mon Jan 27 18:49:34 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to scATS's documentation! 
=================================
scATS is a tool designed to quantify alternative transcription start site (ATS) events from single-end or paired-end 5'-end single-cell RNA sequencing (scRNA-seq) data. It achieves this by modeling the decreasing trend of short-read counts in exons from the 3' to the 5' ends of transcripts, which is attributable to RNA degradation (RD).
Users are required to specify the genes and cells to be included in the analysis, as well as the minimum number for each gene. With these parameters, scATS can accurately calculate the expected abundance of ATS isoforms and conduct differential analysis for specific groups.

To use scATS, you can follow the instructions below:

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   md/Installation
   md/Usage
   md/Getting Help
   md/Citing scATS


Citing scATS
------------
If you use scATS in your research, please cite **********************************

Contacts
--------
scATS is developed and maintained by `Zijie Xu <http://ramble.3vshej.cn>`_ and `Chao Tang <http://ramble.3vshej.cn>`_ 
from Sichuan University. If you want to contribute or have any questions, 
please leave an issue in our `repository <http://ramble.3vshej.cn>`_.

Thank you !

