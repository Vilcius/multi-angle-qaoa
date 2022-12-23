#!/usr/bin/env bash
# author: Anthony Wilkie

# The directory where all the results are stored
# if on the server 
rdir="$HOME/papers/ma_qaoa/Results"
# if on local device
# rdir="$HOME/Papers/Maxcut/multi-angle-qaoa/Results"
tdir=$rdir

# File variables
param_file="QAOA_parameters_mod.f90"
maqaoa_file="QAOA_BFGS_ManyAngles.f90"

# choose the parameters for the current test(s)
runs=100
p_max_ma=3
p_max_normal=6
multi_angle=1
more_bfgs=0

# Choose the heuristic for the test
max_deg_beta=0
max_deg_gamma=0
prev_layer=0


# Set parameter file and directory for previous layer heuristic
if [[ $prev_layer == 1 ]]; then
    tdir="${tdir}/prev_layer"
    mkdir $tdir
fi


# Set parameter file and directory for max degree heuristic
if [[ $max_deg_beta == 1 && $max_deg_gamma == 1 ]]; then
    tdir="${tdir}/max_degree_both"
    mkdir $tdir

    sed -i "s/\(set_max_degree_\(beta\|gamma\)=\).false./\1.true./g" $param_file
elif [[ $max_deg_beta == 1 ]];
then
    tdir="${tdir}/max_degree_beta"
    mkdir $tdir

    sed -i "s/\(set_max_degree_beta=\).false./\1.true./g" $param_file
    sed -i "s/\(set_max_degree_gamma=\).true./\1.false./g" $param_file
elif [[ $max_deg_gamma == 1 ]];
then
    tdir="${tdir}/max_degree_gamma"
    mkdir $tdir

    sed -i "s/\(set_max_degree_beta=\).true./\1.false./g" $param_file
    sed -i "s/\(set_max_degree_gamma=\).false./\1.true./g" $param_file
else
    sed -i "s/\(set_max_degree_\(beta\|gamma\)=\).true./\1.false./g" $param_file
fi


# set parameters if doing ma-QAOA
if [[ $multi_angle == 1 ]]; then
    p_max=$p_max_ma
    sed -i "s/\(many_angles=\).false./\1.true./g" $param_file
    sed -i "s/\(variable_\(beta\|gamma\)=\).false./\1.true./g" $param_file
else
    p_max=$p_max_normal
    sed -i "s/\(many_angles=\).true./\1.false./g" $param_file
    sed -i "s/\(variable_\(beta\|gamma\)=\).true./\1.false./g" $param_file
fi


# run the tests for each value of p
for (( i=1 ; i<=$p_max ; i++ ));
do
    if [[ $do_ma == 1 ]]; then
        tdir="${tdir}/ma_p=$i"
        mkdir $tdir
    else
        tdir="${tdir}/normal_p=$i"
        mkdir $tdir
    fi

    # loop through the runs
    for (( ii=1 ; ii<=$runs ; ii++ ));
    do
        tdir="${tdir}/run_${ii}"

        sed -i "s/\(save_folder='\).*/\1$tdir\/'/g" $param_file
        sed -i "s/\(p_max=\)\d\+/\1$i/g" $param_file

        # if we are using angles from smaller p values, then change the corresponding parameters and file name
        if [[ $prev_layer == 1  ]] && [[ $i -ge 2 ]]; then
            sed -i "s/\(angles_from_smaller_p=\).false./\1.true./g" $param_file
            sii=$(($i-1))
            sed -i "s/\( smaller_p=\)\d\+/\1$sii/g" $param_file
            sed -i "s/\(open(11,file=\1.*,/\1$tdir\/QAOA_dat',/g" $maqaoa_file
        else
            sed -i "s/\(angles_from_smaller_p=\).true./\1.false./g" $param_file
        fi

        # print to the terminal what test we are on
        echo "p=$i: run $ii"
        gfortran $maqaoa_file -o q.exe && ./q.exe
    done
done

