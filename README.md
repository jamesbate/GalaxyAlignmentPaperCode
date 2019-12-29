# GalaxyAlignmentPaperCode

This repository contains the python code I used to analyse data from the Horizon AGN "Universe in a Box" simulation for a publication in the Monthly Notices of the Royal Astronomical Society (MNRAS). The paper 
looked at the evolution of galaxy aligment in this simulation. This code was designed to send instructions 
on a computer cluster and to manage up to hundreds of Gigabytes of data. 

## Broad Structure

Each important function has a docstring describing its role

main_iter.py is the central script which controls which functions you want call. main_iter_callfunction.py
is a small script which iteratively sends jobs to the cluster working on different classes of data in parallel. The primary functions called by the script are the following:

* merger_code_4 : extracts galaxy information 
* TidalField2 : loops over extensive tidal field data to gain tidal information at galaxy locations
* Tidal_Plots_2 : plots this galaxy/tidal information
* z_avangle_plot : creates plots as a function of redshift

