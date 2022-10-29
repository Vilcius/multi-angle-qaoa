#!/usr/bin/env bash

runs=5
tdir="tesst"

# make the test directory
mkdir $tdir

# for (( m=1 ; m<=1 ; m++ ));

# BFGS Loops

for (( b=1 ; b<=2 ; b++ ));
do
    bdir="bfgs_${b}_gamma"
    mkdir $tdir/$bdir
    sed -i "12s/.*/integer, parameter :: min_good_loops=$b/g" QAOA_parameters_mod.f90
    do_ma=1

    # if multi-angle or not
    if [[ $do_ma == 1 ]]; then
            max_p=2
            sed -i '34s/.*/logical, parameter :: many_angles=.true./g' QAOA_parameters_mod.f90
            sed -i '35s/.*/logical, parameter :: variable_beta=.true.,variable_gamma=.true./g' QAOA_parameters_mod.f90
        else
            max_p=6
            sed -i '34s/.*/logical, parameter :: many_angles=.false./g' QAOA_parameters_mod.f90
            sed -i '35s/.*/logical, parameter :: variable_beta=.false.,variable_gamma=.false./g' QAOA_parameters_mod.f90
    fi

    # loop through the p's
    for (( i=1 ; i<=$max_p ; i++ ));
    do
        if [[ $do_ma == 1 ]]; then
            qdir="ma_p=$i"
        else
            qdir="normal_p=$i"
        fi

        mkdir $tdir/$bdir/$qdir
        # mkdir $tdir/$qdir

        echo "p=$i"

        # loop through the runs
        for (( ii=1 ; ii<=$runs ; ii++ ));
        do
            sed -i "5s/.*/character*200, parameter :: save_folder='$tdir\/$bdir\/$qdir\/run_$ii\/'/g" QAOA_parameters_mod.f90
            # sed -i "5s/.*/character*200, parameter :: save_folder='$tdir\/$qdir\/run_$ii\/'/g" QAOA_parameters_mod.f90
            sed -i "10s/.*/integer, parameter :: p_max=$i/g" QAOA_parameters_mod.f90

            echo "run $ii"
            gfortran QAOA_BFGS_ManyAngles.f90 -o q.exe && ./q.exe
        done
    done
done