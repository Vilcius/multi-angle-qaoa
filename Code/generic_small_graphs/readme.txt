The main program is "QAOA_BFGS_ManyAngles.f90" 
Adjustable parameters are in "QAOA_parameters_mod.f90"
Supporting subroutines are in "BFGS_mod.f90" and "QAOA_subroutines_mod.f90". 
The BFGS_mod.f90 is adapted from Ref. [52] at https://arxiv.org/pdf/2109.11455.pdf

To compile I use

"gfortran QAOA_BFGS_ManyAngles.f90 -o q.exe -O3" 

and to run "./q.exe"


The program generates output analogous to the n=8 results.  