include'QAOA_parameters_lowmem_mod.f90'
include'QAOA_subroutines_lowmem_mod.f90'
include"BFGS_mod.f90"
include"nlopt_mod.f90"
Program QAOA
use parameters
use QAOA_subroutines
use BFGS
use nlopt
implicit none

integer :: i,j,k,x,y,z,n,m
integer :: iter,tmp_int,its
double precision :: fret
double precision :: timef,time0
double precision :: C,p_cz_max,C_test,angles_test(2*p_max),C_tester
double precision :: angles(2*p_max)
double precision, allocatable :: angles_ma(:)
double precision :: tmp
double precision, allocatable :: C_min(:)

logical, parameter :: use_random_graph_subset=.false., generate_random_graph_set=.false.
integer, parameter :: n_random_graphs=100
integer :: BFGS_loops=0
double precision :: t0,t1

character*1 :: Mixer='X'
character*3 :: g_string,it_string3
character*2 :: it_string2
character*1 :: it_string1

call cpu_time(time0)

call system('mkdir '//trim(save_folder))
call system('cp QAOA_ma_trianglefree_lowmem.f90 '//trim(save_folder)//'QAOA_ma_trianglefree_lowmem.f90')
call system('cp QAOA_parameters_lowmem_mod.f90 '//trim(save_folder)//'QAOA_parameters_lowmem_mod.f90')
call system('cp QAOA_subroutines_lowmem_mod.f90 '//trim(save_folder)//'QAOA_subroutines_lowmem_mod.f90')
call system('cp BFGS_mod.f90 '//trim(save_folder)//'BFGS_mod.f90')

open(1,file=trim(graph_file),status='old')
call Count_num_graphs
close(1)

allocate(C_opt(graph_num_tot),Beta_opt(graph_num_tot,p_max),Gamma_opt(graph_num_tot,p_max),C_min(graph_num_tot), &
		& n_angles_ma(graph_num_tot),beta_opt_ma(graph_num_tot,n_qubits*p_max),&
		& gamma_opt_ma(graph_num_tot,n_edges_max),cz_max_save(graph_num_tot))
C_opt=-10.d0**10
C_min = 10.d0**10
gamma_opt_ma=0.d0
beta_opt_ma=0.d0
call Read_Cz_max

do its = 1,min_good_loops

	open(1,file=trim(graph_file),status='old')
		
	do graph_num = 1,graph_num_tot

		call Read_Adjacency_Matrix	

		if (more_graphs) then

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
			n_angles_ma(graph_num) = n_beta_angles + n_gamma_angles

			if (many_angles) allocate(gamma_ma(n_gamma_angles),beta_ma(n_beta_angles),angles_ma(n_angles_ma(graph_num)))

			call Generate_angles(n_angles_ma(graph_num),angles_ma)

			if (regular_qaoa) then
				do i = 1,2
					call random_number(tmp)
					angles(i) = tmp*2.d0*pi
				enddo
				call cpu_time(t0)
				!call dfpmin(angles,2,gtol,iter,fret,C_trianglefree,dfunc)
				call optimize_nlopt(C_trianglefree_nlopt,2,angles,fret)
				call cpu_time(t1)
			else
				call cpu_time(t0)
				call dfpmin(angles_ma,n_angles_ma(graph_num),gtol,iter,fret,C_ma_trianglefree,dfunc)
				!call dfpmin(angles,2,gtol,iter,fret,C_trianglefree,dfunc)
				call cpu_time(t1)
			endif

			C=-fret*C_rescaling_fac
			
		    !call Check_optimum_degeneracies(n_angles_ma(graph_num),angles_ma,C)
	     	if (C .gt. (C_opt(graph_num)) ) then
		    	C_opt(graph_num) = C
		    	if (regular_qaoa) then
		    		Beta_opt_ma(graph_num,1) = angles(1)
		    		Gamma_opt_ma(graph_num,1) = angles(2)
		    	else
		    		Beta_opt_ma(graph_num,1:n_beta_angles) = angles_ma(1:n_beta_angles)
		    		Gamma_opt_ma(graph_num,1:n_angles_ma(graph_num)-n_beta_angles) = &
		    			& angles_ma(n_beta_angles+1:n_angles_ma(graph_num))
		    	endif
		    endif
		    if (C .lt. C_min(graph_num)) C_min(graph_num) = C

			deallocate(gamma_ma,beta_ma,angles_ma)

		endif

	enddo
	close(1)
	call cpu_time(timef)
enddo

call cpu_time(timef)
print*, 'c_opt:',C_opt(1)
print*, 'elapsed time:', timef-time0
print*, 'time per it:', (timef-time0)/dble(min_good_loops)
open(2,file=trim(save_folder)//'QAOA_dat_ma_tf',status='unknown')

do graph_num = 1,graph_num_tot
	if (regular_qaoa) write(2,*) graph_num, cz_max_save(graph_num), C_min(graph_num), C_opt(graph_num), &
		(beta_opt_ma(graph_num,j)/pi,j=1,1), (gamma_opt_ma(graph_num,j)/pi,j=1,1)
	if (.not. regular_qaoa) write(2,*) graph_num, cz_max_save(graph_num), C_min(graph_num), C_opt(graph_num), &
		(beta_opt_ma(graph_num,j)/pi,j=1,n_beta_angles), (gamma_opt_ma(graph_num,j)/pi,j=1,n_angles_ma(graph_num)-n_beta_angles)
enddo
close(2)
print*, 'elapsed time:', timef-time0


end program QAOA

subroutine Read_Cz_max

	use parameters
	implicit none
	integer :: g
	double precision :: read_dat(2)

	open(531,file=cmax_file,status='old')
	do g = 1,graph_num_tot
		read(531,*) read_dat(1), read_dat(2)
		cz_max_save(g) = read_dat(2)
	enddo
	close(531)

end subroutine Read_Cz_max
