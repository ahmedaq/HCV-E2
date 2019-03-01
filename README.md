## Table of Contents
*  [Overview](#overview)
*  [Details](#details)
*  [Requirements](#requirements)
*  [Usage](#usage)
*  [Troubleshooting](#troubleshooting)

## Overview
This software package comprises scripts that include the main functions for reproducing the results presented in the paper titled "Identifying immunologically-vulnerable regions of the HCV E2 glycoprotein and broadly neutralizing antibodies that target them".


## Details
#### Title of paper
Identifying immunologically-vulnerable regions of the HCV E2 glycoprotein and broadly neutralizing antibodies that target them
#### Authors
Ahmed A. Quadeer, Raymond H. Y. Louie, and Matthew R. McKay

## Requirements
1.  A PC with MATLAB (preferrably v2016b or later) installed on it with the following additional toolboxes:
    * Bioinformatics Toolbox
    * Communications System Toolbox
    * Statistics and Machine Learning Toolbox
    * Curve Fitting Toolbox
    * Parallel Computing Toolbox
    * MATLAB Distributed Computing Server
    * MATLAB supported compiler installed for compiling C files
 
2.  For inferring maximum-entropy model parameters:
    * Minimum probability flow method (MPF-BML), available at https://github.com/raymondlouie/MPF-BML 
    
3.  For computing the mean escape time metric using the available code, extensive computations are required. The reported escape times in the manuscript were computed using ~1 million core-hours on a supercomputer with Intel(R) Xeon(R) CPU E5-2692 v2 @ 2.20GHz.
 
4.  For mapping fitness costs on the available E2 crystal structure
    * Pymol, available at https://pymol.org/     

## Usage
1.  Open MATLAB
2.  Run the script ```main.m```
3.  By default, all data files are provided to reproduce the results in the paper


For visualizing the step-by-step procedure and the corresponding output
1. Download the html folder
2. Open the ```main.html``` file in your browser

## Troubleshooting
For any questions or comments, please email at ahmedaq@gmail.com. 
