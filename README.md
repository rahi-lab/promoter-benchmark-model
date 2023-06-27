# promoter-benchmark-model

Matlab code for identification of the parameters of gene expression dynamics by fitting the experimentally observed data to a 2-step model of promoter induction.

Input to the code are the values of single-cell expression data obtained by fluorescence microscopy.
The code fits these data to an ODE model of gene expression thus obtaining the values of parameters that describe the leakiness, speed of induction, time-delays upon induction and shut-off, and the degradation rate.

For more details on how the code functions look at the comments in the Main.m file and the accompanying manuscript: https://www.nature.com/articles/s41467-023-38959-8

The code is written, run and tested in Matlab 2019a by Vojislav Gligorovski and Sahand Jamal Rahi.
