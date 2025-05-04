.. INQuant documentation master file, created by
   sphinx-quickstart on Mon Apr  7 09:01:00 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to INQuant's documentation!
===================================

The INQuant package is a tool for quantitative analysis of proteomics data. Below you'll find the main sections of the documentation.

Contents
========

.. toctree::
   :maxdepth: 2

   user_guide
   modules
   
============
Introduction
============

This is the documentation of the INQuant algorithms, designed for quantitative analysis of InstaNovo proteomics predictions. 

This guide contains examples of usage of both the INQuant class and the Calibration class.

INQuant
---------

The INQuant class contains the algorithms for the quantification of InstaNovo predictions. It loads predictions and experiment (mzML) files and performs the quantification.
There is a number of user options for the quantification in each step of the process. 

.. image:: _static/flow_chart.jpg
    :alt: Flow chart of the INQuant pipeline
    :align: center
    :width: 80em



Calibration
----------------

The Calibration class is designed to be worked into the InstaNovo pipeline to calibrate the predictions in a first pass which will modify the mzML files with the calibrated data before running InstaNovo on the calibrated experiments.
As of now, the Calibration is not implemented in the InstaNovo pipeline, so the current implementation can be used as part of the quantification process.
This is demonstrated in the examples under the User Guide. 


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`