include'QAOA_parameters_mod.f90'
include'QAOA_subroutines_mod.f90'
include"BFGS_mod.f90"
Program QAOA
use parameters
use QAOA_subroutines
use BFGS
implicit none

integer :: i,j,k,x,y,z,n,m
integer :: iter,tmp_int
double precision :: fret
double precision :: timef,time0
double precision :: C,p_cz_max,C_test,angles_test(2*p_max),C_tester
double precision :: C_max_c, p_max_c, beta_max_c(p_max), gamma_max_c(p_max)
double precision :: C_max_p, p_max_p, beta_max_p(p_max), gamma_max_p(p_max)
double precision :: prob_cmax
double precision :: angles(2*p_max)
double precision, allocatable :: angles_ma(:)
double precision :: tmp

logical :: more_loops=.true., keepgoing
logical, parameter :: use_random_graph_subset=.true., generate_random_graph_set=.false.
integer, parameter :: n_random_graphs=100
logical :: check_this_graph
integer :: random_graphs(n_random_graphs)
integer, allocatable :: random_graphs_reverse_mapping(:)
double precision :: it_dat(n_random_graphs,min_good_loops,201)
integer :: BFGS_loops=0

double precision :: grad(2*p_max), norm
character*1 :: Mixer='X'
character*3 :: g_string,it_string3
character*2 :: it_string2
character*1 :: it_string1

call cpu_time(time0)

call system('mkdir '//trim(save_folder))
call system('cp QAOA_BFGS_ManyAngles.f90 '//trim(save_folder)//'QAOA_BFGS_ManyAngles.f90')
call system('cp QAOA_parameters_mod.f90 '//trim(save_folder)//'QAOA_parameters_mod.f90')
call system('cp QAOA_subroutines_mod.f90 '//trim(save_folder)//'QAOA_subroutines_mod.f90')
call system('cp BFGS_mod.f90 '//trim(save_folder)//'BFGS_mod.f90')


!set up the basis for calculations
call enumerate_basis
call calc_pair_list_fast

open(1,file=trim(graph_file),status='old')
call Count_num_graphs
close(1)


allocate(C_opt(graph_num_tot),Beta_opt(graph_num_tot,p_max),Gamma_opt(graph_num_tot,p_max), &
		& C0(graph_num_tot),cz_max_save(graph_num_tot),p_C_max(graph_num_tot), &
		& p0_C_max(graph_num_tot),good_loops(graph_num_tot),vertex_degrees(graph_num_tot,0:n_qubits-1), &
		& all_odd_degree(graph_num_tot),all_even_degree(graph_num_tot), &
		& random_graphs_reverse_mapping(graph_num_tot))

if (many_angles) then
    allocate(n_angles_ma(graph_num_tot),beta_opt_ma(graph_num_tot,n_qubits*p_max),&
		& gamma_opt_ma(graph_num_tot,n_edges_max*p_max))
    gamma_opt_ma=0.d0
    beta_opt_ma=0.d0
endif

it_dat=1000.d0

!print *, 'graph_num_tot:', graph_num_tot
good_loops=0
vertex_degrees=0
C_opt=0.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! choose which graphs we are using
!print *, 'random graph set:'
if (use_random_graph_subset) then
	open(2,file='random_graph_set',status='unknown')
	if (generate_random_graph_set) then
		do i = 1,n_random_graphs
			keepgoing=.true.
			do while(keepgoing)  
				call random_number(tmp)
				graph_num=nint(tmp*graph_num_tot+0.5d0)
				keepgoing=.false.
				do j = 1,i
					if (random_graphs(j) .eq. graph_num) then
						keepgoing=.true.
					endif
				enddo
			enddo
			random_graphs(i) = graph_num
			random_graphs_reverse_mapping(random_graphs(i))=i
			write(2,*) i,random_graphs(i)
			print*, i,random_graphs(i)
		enddo
	else
		do i = 1,n_random_graphs
			read(2,*) tmp_int,random_graphs(i)
			random_graphs_reverse_mapping(random_graphs(i))=i
		enddo
	endif
	close(2)
endif


if (single_graph) then
	good_loops=min_good_loops
	good_loops(single_graph_num)=0
endif

if (many_angles) allocate( gamma_vec_ma(0:2**n_qubits-1) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! start doing the algorithm?
do while (more_loops)

	open(1,file=trim(graph_file),status='old')

	graph_num=0
	BFGS_loops=BFGS_loops+1
	
	more_graphs=.true.

    open(11,file='smaller_p/bfgs_3_normal/ma_p=2/run_10/QAOA_dat', status='old')

	do while (more_graphs)
		!loop over the graphs of a given n_qubits

		graph_num = graph_num+1

		if (random_graph) then
			call Random_edges
		else
			call Read_Adjacency_Matrix
		endif

		if (use_random_graph_subset) then
			check_this_graph=.false.
			do i = 1,n_random_graphs
				if (random_graphs(i) .eq. graph_num) then
					check_this_graph=.true.
				endif
			enddo
		else
			check_this_graph=.true.
		endif

		if (more_graphs) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! start to choose angles for algorithm
			if (variable_gamma) then
				n_gamma_angles = p_max*n_edges
			else
				n_gamma_angles = p_max
			endif

			if (variable_beta) then
				n_beta_angles = p_max*n_qubits
			else
				n_beta_angles = p_max
			endif
            if (many_angles) then
                n_angles_ma(graph_num) = n_beta_angles + n_gamma_angles
            endif

			if ( (good_loops(graph_num) .lt. min_good_loops)) then
				if (random_graph) more_graphs=.false.

				call Calc_cz_vec(BFGS_loops)

				if ( check_this_graph ) then
					
					if (many_angles) then
                        allocate(gamma_ma(n_gamma_angles),beta_ma(n_beta_angles),angles_ma(n_angles_ma(graph_num)))
                        call Generate_angles(n_angles_ma(graph_num),angles_ma)
                    else
                        call Generate_angles(2*p_max, angles)
                    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! if we want to optimize <C>
					if (optimize_Expec_C) then
						if (many_angles) then

							call dfpmin(angles_ma,n_angles_ma(graph_num),gtol,iter,fret,QAOA_ExpecC_ma,dfunc)
							C=-fret*cz_max_save(graph_num)
                        else

							call dfpmin(angles,2*p_max,gtol,iter,fret,QAOA_ExpecC,dfunc)
							C=-fret*cz_max_save(graph_num)
						endif
						
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! if <C> > <C>_opt + epsilon, then run ma-QAOA, and set loops to 0
				     	if (C .gt. (C_opt(graph_num)+convergence_tolerance) ) then
					    	C_opt(graph_num) = C
					    	if (many_angles) then
                                beta_opt_ma(graph_num,1:n_beta_angles) = angles_ma(1:n_beta_angles)
                                gamma_opt_ma(graph_num,1:n_angles_ma(graph_num)-n_beta_angles) = angles_ma(&
                                    n_beta_angles+1:n_angles_ma(graph_num))
					    		p_C_max(graph_num) = -QAOA_pmax_ma(n_angles_ma(graph_num),angles_ma)
					    	else
                                beta_opt(graph_num, 1:p_max) = angles(1:p_max)
                                gamma_opt(graph_num, 1:p_max) = angles(p_max+1:2*p_max)
					    		p_C_max(graph_num) = -QAOA_pmax(2*p_max,angles)
					    	endif
					    	good_loops(graph_num) = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! else, don't run ma-QAOA and loop again
					    else
					    	good_loops(graph_num) = good_loops(graph_num)+1
					    endif
					endif

					if ( (single_graph) .and. (graph_num .eq. single_graph_num) .and. (good_loops(graph_num) .eq. min_good_loops) ) then
						call QAOA_algorithm(p_max, beta_opt(graph_num,:), gamma_opt(graph_num,:))
						open(55,file=trim(save_folder)//'pcz',status='unknown')
						tmp = 0.d0
						do z = 0,2**n_qubits-1
							write(55,*) 1.d0/2.d0**n_qubits, dble(psi(z)*dconjg(psi(z))), cz_vec(z)
							tmp = tmp + dble(psi(z)*dconjg(psi(z))) * cz_vec(z)
						enddo
						print *, 'final calculated <C> = ', tmp
					endif

					if (many_angles) deallocate(gamma_ma,beta_ma,angles_ma)

			    endif

			endif


		endif
	enddo

	close(1)
	close(11)

	more_loops=.false.
	do i = 1,graph_num_tot
		if (good_loops(i) .lt. min_good_loops) more_loops=.true.
	enddo

	if (old_angles) more_loops=.false.
	if (fixed_number_iterations .and. (BFGS_loops .eq. min_good_loops) ) more_loops=.false.

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! writing to the QAOA_dat file
open(2,file=trim(save_folder)//'QAOA_dat',status='unknown')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! If we are using a single graph, only print that data
if (single_graph) then
    if (many_angles) then
        write(2,*) single_graph_num, cz_max_save(single_graph_num), C0(single_graph_num), C_opt(single_graph_num), &
            p_C_max(single_graph_num), p_max, (beta_opt_ma(single_graph_num,j)/pi,j=1,n_beta_angles), &
            (gamma_opt_ma(single_graph_num,j)/pi,j=1,n_angles_ma(single_graph_num)-n_beta_angles)
    else
        write(2,*) single_graph_num, cz_max_save(single_graph_num), C0(single_graph_num), C_opt(single_graph_num), &
            p_C_max(single_graph_num), p_max, (beta_opt(single_graph_num,j)/pi,j=1,p_max), &
            (gamma_opt(single_graph_num,j)/pi,j=1,p_max)
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! If we are using a subset of graphs, only print that data
elseif (use_random_graph_subset) then
    do i = 1,n_random_graphs
        graph_num = random_graphs(i)
        if (many_angles) then
            write(2,*) graph_num, cz_max_save(graph_num), C0(graph_num), C_opt(graph_num), p_C_max(graph_num), p_max, &
                (beta_opt_ma(graph_num,j)/pi,j=1,n_beta_angles), &
                (gamma_opt_ma(graph_num,j)/pi,j=1,n_angles_ma(graph_num)-n_beta_angles)
        else
            write(2,*) graph_num, cz_max_save(graph_num), C0(graph_num), C_opt(graph_num), p_C_max(graph_num), p_max, &
                (beta_opt(graph_num,j)/pi,j=1,p_max), &
                (gamma_opt(graph_num,j)/pi,j=1,p_max)
        endif
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Else, we are using all the graphs
else
    do graph_num = 1,graph_num_tot
        if (many_angles) then
            write(2,*) graph_num, cz_max_save(graph_num), C0(graph_num), C_opt(graph_num), p_C_max(graph_num), p_max, &
                (beta_opt_ma(graph_num,j)/pi,j=1,n_beta_angles), &
                (gamma_opt_ma(graph_num,j)/pi,j=1,n_angles_ma(graph_num)-n_beta_angles)
        else
            write(2,*) graph_num, cz_max_save(graph_num), C0(graph_num), C_opt(graph_num), p_C_max(graph_num), p_max, &
                (beta_opt(graph_num,j)/pi,j=1,p_max), &
                (gamma_opt(graph_num,j)/pi,j=1,p_max)
        endif
    enddo
endif

!do graph_num = 1,graph_num_tot
!	if (dabs(c_opt(graph_num) - cz_max_save(graph_num)) .lt. 1.d-4) then
!		if (.not. single_graph) print *, graph_num,'fully optimized!'
!	endif
!enddo
close(2)

! open(2,file=trim(save_folder)//'Conv_dat',status='unknown')
! do i = 1,n_random_graphs
! 	do j = 1,min_good_loops
! 		write(2,*) (it_dat(i,j,k),k=1,201)
! 	enddo
! enddo
! close(2)

print *, 'BFGS_loops:', BFGS_loops
call cpu_time(timef)
print*, 'elapsed time:', timef-time0


end program QAOA

