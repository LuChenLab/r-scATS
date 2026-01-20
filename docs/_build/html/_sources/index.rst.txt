.. scATS documentation master file, created by
   sphinx-quickstart on Mon Jan 27 18:49:34 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to scATS's documentation! 
=================================
**scATS** is a tool designed to **de novo** quantify alternative transcription start site (ATS) from paired-end 5'-end single-cell RNA sequencing (scRNA-seq) data. It achieves this by modeling the decreasing trend of short-read counts in exons from the 3' to the 5' ends of transcripts, which is attributable to **RNA degradation**.


**scATS** introduces two novel metrics (`α` and `β`) to precisely quantify RNA degradation, and provides both raw (`ψ`) and corrected (`θ`) metrics to correct ATS expression profiles from resulting distortions.

.. image:: _static/method.png
   :alt: scATS Workflow Diagram
   :align: center
   :scale: 40%

Users are required to specify the genes and cells to be included in the analysis, as well as the minimum number for each gene. With these parameters, **scATS** can accurately calculate the expected abundance of ATS isoforms and conduct differential analysis for specific groups.


To use scATS, you can follow the instructions below:

.. toctree::
   :caption: Contents:
   :maxdepth: 2

   md/Installation
   md/Usage
   md/Disease model
   md/Getting Help
   md/Citing scATS


Citing scATS
------------
If you use scATS in your research, please cite **********************************

Contacts
--------
scATS is developed and maintained by `Zijie Xu <https://github.com/kayla-xu>`_ and `Chao Tang <https://github.com/ChaoTang-SCU>`_ 
from Sichuan University. If you want to contribute or have any questions, 
please leave an issue in our `repository <https://github.com/LuChenLab/r-scATS/issues>`_.

Thank you !

