Code for triangle-free graphs, using the analytical formula of Theorem IV.1 at https://arxiv.org/pdf/2109.11455.pdf

The main program is "QAOA_ma_trianglefree_lowmem.f90"
adjustable parameters are in "QAOA_parameters_lowmem_mod.f90"
The remaining files are supporting subroutines
"QAOA_subroutines_lowmem_mod.f90" - QAOA state calculations
"BFGS_mod.f90" - BFGS optimization from Ref. [52]
"nlopt_mod.f90" - NLopt optimization from Refs. [53,54]

To compile you will need a local install of fortran nlopt. 

I compile with "gfortran QAOA_ma_trianglefree_lowmem.f90 -o q.exe -O3 -I/usr/local/include -L/usr/local/lib -lnlopt -lm"

run with "./q.exe"

The output is analogous the to n=50, n=100 results. 

--
The program "QAOA_ma_trianglefree_lowmem_regqaoa.f90" will do the p=1 optimization of regular QAOA when "regular_qaoa=.true." in the parameters file. 
