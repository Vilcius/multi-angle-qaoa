
module parameters

implicit none

character*200 :: save_folder='test/'
character*200, parameter :: graph_file='../../Graph_files/Erdos-Renyi/Adj50_tot.txt'
character*200, parameter :: cmax_file='../../Graph_files/Erdos-Renyi/Adj50_cmax.txt'

integer, parameter :: n_qubits=50
integer, parameter :: p_max=1
integer, parameter :: min_good_loops=100

integer, parameter :: n_edges_max=n_qubits*(n_qubits-1)/2
integer :: n_edges, edges(0:n_edges_max-1,0:1)
double precision, parameter :: pi=dacos(-1.d0)
double precision, parameter :: gtol=1.d-14, degen_tol=1.d-8,convergence_tolerance=1.d-12

logical, parameter :: rescale_angles=.true.
double precision, parameter :: angle_rescaling_fac=1.d0
double precision, parameter :: C_rescaling_fac=1.d0

logical, parameter :: regular_qaoa=.false.
logical, parameter :: many_angles=.true.
logical, parameter :: variable_beta=.true.,variable_gamma=.true.
integer :: n_beta_angles=n_qubits,n_gamma_angles
double precision, parameter :: gamma_max = pi, beta_max=pi/2.d0

double precision, allocatable :: gamma_ma(:), beta_ma(:), gamma_vec_ma(:),gamma_opt_ma(:,:),beta_opt_ma(:,:)

logical :: more_graphs=.true.

integer, allocatable :: n_angles_ma(:)
double precision :: beta_a(p_max),gamma_a(p_max),cz_max
double precision, allocatable :: cz_max_save(:), C0(:), C_opt(:), p_C_max(:), beta_opt(:,:), gamma_opt(:,:),p0_C_max(:)

integer :: graph_num_tot,graph_num

end module parameters
