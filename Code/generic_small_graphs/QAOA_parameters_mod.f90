module parameters

implicit none

character*200, parameter :: save_folder='wtf_angles/ma_p=3/run_5/'
character*200, parameter :: graph_file='../../Graph_Files/n=8/graph8c.txt'

integer, parameter :: n_qubits=8
integer, parameter :: smaller_p=1
integer, parameter :: p_max=3

integer, parameter :: min_good_loops=1
integer, parameter :: graph_of_interest=1

integer, parameter :: dim=2**n_qubits, n_edges_max=n_qubits*(n_qubits-1)/2
integer, parameter :: min_it_gamma=0, max_it_gamma=101*40
integer, parameter :: min_it_beta=0, max_it_beta=26*2
double precision, parameter :: pi=dacos(-1.d0)
double precision, parameter :: gtol=1.d-14, degen_tol=1.d-8,convergence_tolerance=1.d-12
logical, parameter :: random_graph=.false.
logical, parameter :: angles_from_smaller_p=.false.
logical, parameter :: old_angles=.false.,directed_angles=.false.
logical, parameter :: force_angle_ordering=.false.
logical, parameter :: even_odd_only=.false., even_odd_angle_dist=.false., not_even_odd_only=.false.
logical, parameter :: no_opt=.false.
logical, parameter :: Optimize_expec_C=.true.,Optimize_p_Cmax=.false.,Optimize_Sq=.false.
logical, parameter :: fix_sign_beta1=.false.,fix_sign_gamma1=.false.,fix_sign_positive=.false.
logical, parameter :: single_graph=.false.
integer, parameter :: single_graph_num = 4115
integer :: graph_num_tot,graph_num
integer :: degeneracy_cz_max
logical, parameter :: fixed_number_iterations=.true.

logical, parameter :: many_angles=.true.
logical, parameter :: variable_beta=.true.,variable_gamma=.true.
integer :: n_beta_angles=n_qubits,n_gamma_angles
!For variable angle QAOA:
double precision, parameter :: gamma_max = pi, beta_max=pi/2.d0

double precision, allocatable :: gamma_ma(:), beta_ma(:), gamma_vec_ma(:),gamma_opt_ma(:,:),beta_opt_ma(:,:)
double precision, allocatable :: Sq_opt(:)

complex*16 :: psi(0:dim-1) 
complex*16 :: Ub(0:1,0:1)

logical :: more_graphs=.true.

double precision :: cz_vec(0:dim-1)

integer :: pairs_list(0:n_qubits-1,0:dim/2-1,0:1),basis(0:dim-1,0:n_qubits-1)
integer :: n_edges, edges(0:n_edges_max-1,0:1)
integer :: n_graphs_even_odd
integer, allocatable :: n_angles_ma(:)
double precision :: beta_a(p_max),gamma_a(p_max),cz_max
double precision, allocatable :: cz_max_save(:), C0(:), C_opt(:), p_C_max(:), beta_opt(:,:), gamma_opt(:,:),p0_C_max(:)
double precision, allocatable :: cz_max_save_weighted(:,:), C0_weighted(:,:), C_opt_weighted(:,:)
double precision, allocatable :: p_C_max_weighted(:,:), beta_opt_weighted(:,:,:), gamma_opt_weighted(:,:,:),p0_C_max_weighted(:,:)
integer :: weight_num

integer, allocatable :: z_max(:)
integer, allocatable :: good_loops(:), vertex_degrees(:,:)
logical, allocatable :: all_odd_degree(:),all_even_degree(:)
logical :: set_max_degree=.false.

end module parameters
