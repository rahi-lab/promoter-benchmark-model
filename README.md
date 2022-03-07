# promoter-benchmark-model
This is the code for promoter induction parameters identification using a 2-step gene expression model.
The purpose of the code is to extract the parameters of gene expression dynamics by fitting the experimentally observed data to a model for promoter induction.

Input to the code are the values of single-cell expression data obtained by fluorescence microscopy.
The code fits this data to a simple ODE model of gene expression thus obtaining the values of parameters that describe the leakiness, spead of induction, time-delays upon induction and shut-off as well as degradation rate.

For more details look at the comments in the Main.m file and the accompanying manuscript: https://www.biorxiv.org/content/10.1101/2020.08.16.253310v1

The code is written, run and tested in Matlab 2019a by Vojislav Gligorovski and Sahand Jamal Rahi.
