module QAOA_subroutines

implicit none

contains


function Morse(x)

implicit none

double precision :: Morse, x(1)
double precision :: D=1.d0, r0=1.d0, a=1.d0

Morse = D*( exp(-2*a*(x(1)-r0)) - 2*exp(-a*(x(1)-r0)) )

end function Morse


function Morse_test(n,x)

implicit none
integer :: n,m
double precision :: Morse_test, x(2*n)
double precision :: D=1.d0, r0=1.d0, a=1.d0

Morse_test = 0.d0
do m = 1,2*n
	Morse_test = Morse_test + D*( exp(-2*a*(x(m)-r0)) - 2*exp(-a*(x(m)-r0)) )
enddo

end function Morse_test



subroutine QAOA_Algorithm_ma
	use parameters
	implicit none

	integer :: j
	integer :: p

	psi=dcmplx(1.d0/dsqrt(dble(dim)),0.d0)

    do j = 1,p_max
    	call Calc_gamma_vec_ma
    	call Apply_Uc_ma
    	call Apply_Ub_ma
    enddo

end subroutine QAOA_Algorithm_ma


subroutine QAOA_Algorithm(p,B_angles,C_angles)
	use parameters
	implicit none

	integer :: j
	integer :: p 
	double precision :: b_angles(p),c_angles(p)!beta and gamma

	psi=dcmplx(1.d0/dsqrt(dble(dim)),0.d0)

    do j = 1,p
    	call construct_Ub(b_angles(j))
    	call Apply_Uc(c_angles(j))
    	call Apply_Ub
    enddo

end subroutine QAOA_Algorithm

function QAOA_ExpecC_ma(n_angles,angles)

	use parameters 
	integer :: n_angles,p
	double precision :: angles(n_angles)
	double precision :: p_cz_max,QAOA_ExpecC_ma


	beta_ma = angles(1:n_qubits*p_max)
	gamma_ma = angles(n_qubits*p_max+1:n_angles)

	Call QAOA_Algorithm_ma
	call Calc_expec_C(QAOA_ExpecC_ma,p_cz_max)

	QAOA_ExpecC_ma = -QAOA_ExpecC_ma/cz_max_save(graph_num)

end function QAOA_ExpecC_ma


function QAOA_ExpecC(n_angles,angles)

	use parameters 
	integer :: n_angles,p
	double precision :: angles(n_angles),B_angles(n_angles/2),C_angles(n_angles/2)
	double precision :: p_cz_max,QAOA_ExpecC
	p = n_angles/2

	B_angles = angles(1:p)
	C_angles = angles(p+1:2*p)

	Call QAOA_Algorithm(p,B_angles,C_angles)
	call Calc_expec_C(QAOA_ExpecC,p_cz_max)

	QAOA_ExpecC = -QAOA_ExpecC/cz_max_save(graph_num)

	!print *, 'angles:', angles
	!print *, 'QAOA_ExpecC:', QAOA_ExpecC
end function QAOA_ExpecC

function QAOA_Sq(n_angles,angles)

	use parameters 
	integer :: n_angles,p
	double precision :: angles(n_angles),B_angles(n_angles/2),C_angles(n_angles/2)
	double precision :: p_cz_max,tmp,QAOA_Sq
	integer :: z
	p = n_angles/2

	B_angles = angles(1:p)
	C_angles = angles(p+1:2*p)

	Call QAOA_Algorithm(p,B_angles,C_angles)

	QAOA_Sq = 0.d0
	do z = 0,2**n_qubits-1
		tmp = dble(psi(z)*dconjg(psi(z)))
		if (tmp .gt. 1.d-12) QAOA_sq = QAOA_sq - tmp*log2(tmp)
	enddo

end function QAOA_Sq


function C_ma_trianglefree(n_angles,angles)
    use parameters
    implicit none
    integer :: n_angles
    double precision :: angles(n_angles)
    integer :: m,l
    double precision :: C_ma_trianglefree
    double precision :: term0,term1
    integer :: q0,q1,edge_angle_index

    C_ma_trianglefree = 0.d0

    !print*, 'angles:',angles
    !print*, 'edges:',edges
    do m = 0,n_edges-1
    	!check the edge exists in edges array
    	if (edges(m,0) .ne. edges(m,1)) then 
    		q0 = edges(m,0)+1
    		q1 = edges(m,1)+1
    		edge_angle_index=n_beta_angles+m+1
    		term1=dsin(angles(edge_angle_index))*dcos(angles(q0)*2.d0)*dsin(angles(q1)*2.d0)
    		term0=dsin(angles(edge_angle_index))*dcos(angles(q1)*2.d0)*dsin(angles(q0)*2.d0)
    		do l = 0,n_edges-1
    			if (  (l .ne. m) .and. ( edges(l,0) .ne. edges(l,1) )  ) then
    				if ( (edges(l,0) .eq. q0-1) .or. (edges(l,1) .eq. q0-1) ) then
    					edge_angle_index = n_beta_angles+l+1
    					term0 = term0*dcos(angles(edge_angle_index))
    				endif
    				if ( (edges(l,0) .eq. q1-1) .or. (edges(l,1) .eq. q1-1) ) then
    					edge_angle_index = n_beta_angles+l+1
    					term1 = term1*dcos(angles(edge_angle_index))
    				endif
    			endif
    		enddo
    		C_ma_trianglefree = C_ma_trianglefree + 0.5d0*(1.d0 + term0 + term1)
    		!if ( dcos(angles(q0)*2.d0) .gt. 1.d-13 ) then
    		!	derivative(q0) = derivative(q0) -2.d0*dsin(angles(q0)*2.d0)/dcos(angles(q0)*2.d0)*term1
    		!endif
    		!if ( dsin( angles(q0)*2.d0 ) .gt. 1.d-13 ) then
    		!	derivative(q0) = derivative(q0) +2.d0*dcos(angles(q0)*2.d0)/dsin(angles(q0)*2.d0)*term0
    		!endif
    	endif
    enddo

    C_ma_trianglefree=-C_ma_trianglefree/dble(n_qubits)
   ! C_ma_trianglefree = 0.5d0*(1.d0 + term0 + term1)

end function C_ma_trianglefree

function QAOA_pmax(n_angles,angles)

	use parameters 
	implicit none
	integer :: n_angles,p
	double precision :: angles(n_angles),B_angles(n_angles/2),C_angles(n_angles/2)
	double precision :: C,QAOA_pmax

	p = n_angles/2

	B_angles = angles(1:p)
	C_angles = angles(p+1:2*p)

	Call QAOA_Algorithm(p,B_angles,C_angles)
	call Calc_expec_C(C,QAOA_pmax)

	!print *, 'graph_num,C,QAOA_pmax',graph_num,C,QAOA_pmax
	QAOA_pmax = -QAOA_pmax

end function QAOA_pmax

function QAOA_pmax_ma(n_angles,angles)

	use parameters 
	integer :: n_angles,p
	double precision :: angles(n_angles)
	double precision :: C,QAOA_pmax_ma

	beta_ma = angles(1:n_qubits*p_max)
	gamma_ma = angles(n_qubits*p_max+1:n_angles)

	Call QAOA_Algorithm_ma
	call Calc_expec_C(C,QAOA_pmax_ma)

	QAOA_pmax_ma = -QAOA_pmax_ma

end function QAOA_pmax_ma


subroutine Generate_angles(n_angles,angles)

	use parameters
	implicit none

	integer :: i
	integer :: n_angles
	double precision :: junk(6)
	double precision :: angles(n_angles)
	double precision :: tmp

	if (old_angles) then
		read(12,*) (junk(i),i=1,6),(angles(i),i=1,2*p_max)
		angles=angles*pi
	else if (directed_angles) then
		do i = 1,p_max
			tmp=0.d0
			!call random_number(tmp)
			angles(i) = -pi/8 + tmp*0.1d0*pi
			!call random_number(tmp)
			angles(p_max+i) = -3*pi/8+ tmp*0.1d0*pi
		enddo
	else if (angles_from_smaller_p) then
		read(11,*) (junk(i),i=1,6),(angles(i),i=1,smaller_p),(angles(i),i=p_max+1,p_max+smaller_p)
		angles=angles*pi
		!add a small random perturbation t-o the previous angles
		do i = 1,p_max-1
			call random_number(tmp)
			angles(i) = angles(i) + pi/2.d0*(tmp-0.5d0)/2.d0
			call random_number(tmp)
			angles(p_max+i) = angles(p_max+i) + 2.d0*pi*(tmp-0.5d0)/2.d0
		enddo
			call random_number(tmp)
			angles(p_max) = pi/2.d0*(tmp-0.5d0)
			call random_number(tmp)
			angles(2*p_max) = 2*pi*(tmp-0.5d0)
		!print *, 'angles:', (angles(i),i=1,2*p)
	else if (many_angles) then
		do i = 1,n_qubits
    		call random_number(tmp)
    		angles(i) = beta_max*2.d0*(tmp-0.5d0)
    	enddo
    	beta_ma = angles(1:n_qubits*p_max)
    	do i = 1,n_edges
    		call random_number(tmp)
    		angles(n_qubits*p_max+i) = gamma_max*2.d0*(tmp-0.5d0)
    	enddo
    	gamma_ma(1:n_angles-p_max*n_qubits) = angles(n_qubits*p_max+1:n_angles)
    elseif (even_odd_only .and. even_odd_angle_dist) then
    	do i = 1,p_max
    		call random_number(tmp)
    		if (i .eq. 1) then
    			angles(i) = -pi/4.d0*tmp !-0.25pi < \beta_1 < 0
    		else
    			angles(i) = -pi/2.d0*(tmp-0.5d0)
    		endif
    		call random_number(tmp)
    		angles(p_max+i) = 1.d0*pi*(tmp-0.5d0)
    	enddo
    	!print *, 'angles/pi:',angles/pi
	else
    	do i = 1,p_max
    		call random_number(tmp)
    		angles(i) = pi/2.d0*(tmp-0.5d0)
    		call random_number(tmp)
    		angles(p_max+i) = 2.d0*pi*(tmp-0.5d0)
    	enddo
    	!print *, 'angles/pi:',angles/pi
    endif

end subroutine Generate_angles


subroutine calc_vertex_degrees

	use parameters

	implicit none

	integer :: n,q0,q1

	do n = 0,n_edges-1
		q0 = edges(n,0)
		q1 = edges(n,1)
		vertex_degrees(graph_num,q0) = vertex_degrees(graph_num,q0)+1
		vertex_degrees(graph_num,q1) = vertex_degrees(graph_num,q1)+1
	enddo

	all_odd_degree(graph_num)=.true.
	all_even_degree(graph_num)=.true.
	do n = 0,n_qubits-1
		if (mod(vertex_degrees(graph_num,n),2) .eq. 0) all_odd_degree(graph_num)=.false.
		if (mod(vertex_degrees(graph_num,n),2) .eq. 1) all_even_degree(graph_num)=.false.
	enddo
	!if (all_odd_degree(graph_num)) print *, graph_num,'all odd degree'
	!if (all_even_degree(graph_num)) print *, graph_num,'all even degree'

end subroutine calc_vertex_degrees


subroutine dfunc(n,x,grad,func)

	implicit none

	integer :: n,i
	double precision :: x(n), x_dif_plus(n), x_dif_min(n)
	double precision :: nominal_delta = 1.d-6
	double precision :: h,temp
	double precision :: grad(n)
	double precision :: func

	External func

	!set h exact to numerical precision
	!so h can be represented exactly in base 2
	temp =  x(1) + nominal_delta
 	h = temp-x(1)

 	do i = 1,n
  		x_dif_plus = x
 		x_dif_min = x
 		x_dif_plus(i) = x(i) + h
 		x_dif_min(i) = x(i) - h
 		grad(i) = (func(n,x_dif_plus) - func(n,x_dif_min))/(2.d0*h)
	enddo

end subroutine dfunc



double precision function log2(x)
  implicit none
  double precision, intent(in) :: x

  log2 = log(x) / log(2.d0)
end function

subroutine Count_Num_graphs

	use parameters 
	implicit none

	!the below counts an extra time on the last call
	!to read_adjacency_matrix
	graph_num_tot=-1
	do while (more_graphs)
		call Read_Adjacency_Matrix
		graph_num_tot = graph_num_tot+1
	enddo
	more_graphs=.true.

end subroutine Count_num_graphs


subroutine Read_Adjacency_Matrix

	use parameters
	implicit none

	integer :: n,nn,m
	integer :: iostatus,edge_count
	character (len=n_qubits-1) :: line
	logical :: non_graph_line


	non_graph_line=.true.
	edge_count=-1
	edges=0
	!skip lines that aren't the adjacency matrix entries
	do while (non_graph_line)
		read(1,*,IOSTAT = iostatus) line
		if (iostatus .lt. 0) then
			!print*, 'end of file'
			non_graph_line=.false.
			more_graphs=.false.
		else if (iostatus .gt. 0) then
			print*, 'error reading file:', iostatus
			stop
		endif
		if (line(1:1) .eq. '0' .or. line(1:1) .eq. '1') then
			non_graph_line = .false.
		endif
	enddo
	if (more_graphs) then
		do n = 0,n_qubits-2
			if (n .ne. 0) read(1,*) line
			do nn = n+1,n_qubits-1
				if (line(nn-n:nn-n) .eq. '1') then
					edge_count = edge_count+1
					edges(edge_count,0) = n
					edges(edge_count,1) = nn
				endif
			enddo
		enddo
	endif
	n_edges=edge_count+1

end subroutine Read_Adjacency_Matrix


subroutine Random_edges

	!Generate a graph with random edges
	use parameters
	implicit none

	double precision :: tmp
	integer :: i,j,k,l
	logical :: old_edge=.true.

	call random_number(tmp)

	tmp=0.5d0
	n_edges = nint(tmp*n_edges_max)
	if (n_edges .eq. 0) n_edges =1
	print*,'n_edges,n_edges_max:',n_edges,n_edges_max
	!do i = 1,n_edges-1
	!	edges(i,0)=0
	!	edges(i,1)=i
	!enddo
	if (.true.) then
		do i = 0,n_edges-1
			old_edge=.true.
			do while (old_edge)
				call random_number(tmp)
				j = nint(tmp*n_qubits-0.5d0)
				call random_number(tmp)
				k = nint(tmp*n_qubits-0.5d0)
				old_edge=.false.
				do l = 0,i
					if ( j .eq. k .or. (edges(l,0) .eq. j .and. edges(l,1) .eq. k) .or. &
					   	& (edges(l,0) .eq. k .and. edges(l,1) .eq. j) ) then
					   	old_edge=.true.
					endif
				enddo
				if (.not. old_edge) then
					edges(i,0) = j
					edges(i,1) = k
				endif
			enddo
		enddo
	endif

end subroutine Random_edges


SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
END SUBROUTINE init_random_seed


subroutine Calc_gamma_vec_ma

	use parameters
	implicit none

	integer :: z,m

	gamma_vec_ma=0.d0
    do z =0,2**n_qubits-1
        do m = 0,n_edges-1
            if (basis(z,edges(m,0)) .ne. basis(z,edges(m,1))) then
           		gamma_vec_ma(z) = gamma_vec_ma(z) + gamma_ma(m+1)
            endif
        enddo
    enddo


end subroutine Calc_gamma_vec_ma


subroutine Calc_cz_vec(BFGS_loops)
    !calculate a vector c_z
    !with entries c_z[z] that equal the cost function evaluation for
    !the binary expansion of z
    use parameters
    implicit none
    integer :: z,m
    integer :: BFGS_loops
    logical :: calc_zmax=.false.

    cz_vec=0.d0

    do z =0,2**n_qubits-1
        do m = 0,n_edges-1
            if (basis(z,edges(m,0)) .ne. basis(z,edges(m,1))) then
            	cz_vec(z) = cz_vec(z) + 1
            	!if (many_angles) then
            	!	gamma_vec_ma(z) = gamma_vec_ma(z) + gamma_ma(m)
            	!endif
            endif
        enddo
    enddo

	if (BFGS_loops .eq. 1) then
		!calculate C0, the <C> for the initial state
		psi=dcmplx(1.d0/dsqrt(dble(dim)),0.d0)
		call Calc_expec_C(C0(graph_num),p0_C_max(graph_num))
		cz_max_save(graph_num) = maxval(cz_vec)

		call calc_vertex_degrees

		if (all_odd_degree(graph_num)) print *, graph_num,'all odd degree'
		if (all_even_degree(graph_num)) print *, graph_num,'all even degree'

		if (all_even_degree(graph_num) .or. all_odd_degree(graph_num)) n_graphs_even_odd = n_graphs_even_odd+1
		
		if (calc_zmax) then
			P0_C_max(graph_num) = 0.d0
			do z = 0,2**n_qubits-1
				if (cz_vec(z) .eq. maxval(cz_vec)) then
					P0_c_max(graph_num) = P0_C_max(graph_num) +1.d0/dble(dim)
					if (dabs( p0_c_max(graph_num) - 1.d0/dble(dim) ) .lt. 1.d-12) then
						z_max(graph_num) = z
					endif
				endif
			enddo
		endif

	endif

end subroutine Calc_cz_vec


subroutine Calc_Pcz_max(Pczmax)

	use parameters
	implicit none

	integer :: z
	double precision :: Pczmax

	pczmax=0.d0
	do z = 0,2**(n_qubits)-1
		if (dabs(dble(cz_vec(z)) -maxval(cz_vec)) .lt. 1.d-8) then
			pczmax=pczmax+1.d0
		endif
		!print *, 'z,cz_vec(z),cz_max,pczmax',z,cz_vec(z),cz_max,pczmax
	enddo
	pczmax=pczmax/dble(2**n_qubits)

end subroutine calc_Pcz_max

subroutine Construct_Ub(b_angle)

	use parameters
	implicit none
	double precision :: b_angle

	Ub(0,0) = dcmplx(dcos(b_angle),0.d0)
	Ub(0,1) = dcmplx(0.d0,-dsin(b_angle))
	Ub(1,0) = dcmplx(0.d0,-dsin(b_angle))
	Ub(1,1) = dcmplx(dcos(b_angle),0.d0)

end subroutine Construct_Ub



subroutine Apply_Uc_ma

	use parameters
	implicit none
	integer :: z

	if (many_angles) then
		do z = 0,2**n_qubits-1
			psi(z) = psi(z) * dcmplx(dcos(gamma_vec_ma(z)),-dsin(gamma_vec_ma(z)))
		enddo
	else
		print *, 'not many angles!'
		stop
	endif

end subroutine Apply_Uc_ma


subroutine Apply_Uc(c_angle)

	use parameters
	implicit none
	integer :: z
	double precision :: c_angle

    do z = 0,2**n_qubits-1
        psi(z) = psi(z)*dcmplx(dcos(cz_vec(z)*c_angle),-dsin(cz_vec(z)*c_angle))
    enddo

end subroutine Apply_Uc


subroutine Calc_expec_C(C,p_cz_max)

	use parameters
	implicit none
	integer :: z
	double precision :: C,p_cz_max

    C=0.d0
    p_cz_max=0.d0
    do z = 0,2**n_qubits-1
        C = C + realpart(psi(z)*dconjg(psi(z)))*dble(cz_vec(z))  
        if (cz_vec(z) .eq. cz_max_save(graph_num)) then
        	p_cz_max = p_cz_max + realpart(psi(z)*dconjg(psi(z)))
        endif
    enddo
end subroutine Calc_expec_C


subroutine Apply_Ub_ma
	use parameters
    implicit none
    integer :: n,z
    complex*16 :: temporary

    !apply Ub to all qubits
    !for each qubit, loop over coupled pairs of states
    !where the basis state of qubit n changes 0->1
    !and all other qubit basis states are constant.
    !The Ub operators are assumed to commute so the
    !order the n are executed in doesn't matter  

    do n = 0,n_qubits-1

    	call construct_Ub(beta_ma(n+1))
        do z = 0, 2**(n_qubits-1)-1
            temporary = Ub(0,0)*psi(pairs_list(n,z,0)) + Ub(0,1)*psi(pairs_list(n,z,1))
            psi(pairs_list(n,z,1)) = Ub(1,1)*psi(pairs_list(n,z,1)) + Ub(1,0)*psi(pairs_list(n,z,0))
            psi(pairs_list(n,z,0)) = temporary
        enddo

    enddo

end subroutine Apply_Ub_ma


subroutine Apply_Ub
	use parameters
    implicit none
    integer :: n,z
    complex*16 :: temporary

    !apply Ub to all qubits
    !for each qubit, loop over coupled pairs of states
    !where the basis state of qubit n changes 0->1
    !and all other qubit basis states are constant.
    !The Ub operators are assumed to commute so the
    !order the n are executed in doesn't matter  

    do n = 0,n_qubits-1

        do z = 0, 2**(n_qubits-1)-1
            temporary = Ub(0,0)*psi(pairs_list(n,z,0)) + Ub(0,1)*psi(pairs_list(n,z,1))
            psi(pairs_list(n,z,1)) = Ub(1,1)*psi(pairs_list(n,z,1)) + Ub(1,0)*psi(pairs_list(n,z,0))
            psi(pairs_list(n,z,0)) = temporary
        enddo

    enddo

end subroutine Apply_Ub


subroutine Print_binary_expansion_string(z,z_str)
	use parameters
	implicit none

	character*10 :: z_str
	integer :: n,z,m

	!print *, 'basis(z):',(basis(z,n),n=0,n_qubits-1)
	z_str = char(48+basis(z,n_qubits-1))
	!print *, '0, z_str',z_str
	do n = 1,n_qubits-1
		m = n_qubits-1-n
		z_str = trim(z_str)//char(48+basis(z,m))
		!print *, 'm,basis(z,m),z_str',m,basis(z,m),z_str
	enddo

	!print *, 'z,z_str',z,z_str

end subroutine Print_binary_expansion_string


subroutine enumerate_basis
	use parameters
	implicit none
	integer :: n,z
	integer :: tmp

	do z = 0,2**n_qubits-1
		tmp=z
		do n = 0,n_qubits-1
			basis(z,n) = mod(tmp,2)
			tmp=tmp/2
		enddo
		!print *, 'z,basis(z,:):',z,basis(z,:)
	enddo

end subroutine enumerate_basis


subroutine calc_pair_list_fast
	!generate a list of pairs of states that are related by switching single
	!qubit states from 0 -> 1
	!resulting list has 2**(n_qubits-1) entries for each qubit n, with pairs
	!of related states where all other qubits are the same and n =0,1
	use parameters
	implicit none
	integer :: z0,n,j
	integer :: tally(0:n_qubits)

	tally=-1
	do z0=0,2**n_qubits-1
    	do n = 0,n_qubits-1
    		if (basis(z0,n) .eq. 0) then
    			tally(n) = tally(n)+1
    			pairs_list(n,tally(n),0) = z0
    			pairs_list(n,tally(n),1) = z0 + 2**n
    		endif
    	enddo
    enddo

end subroutine calc_pair_list_fast


subroutine Calc_pair_list
    !Make a list of basis state numbers that differ in one qubit
    !These will be used later to calculate single-qubit unitaries
    !acting on the quantum state
    use parameters
    implicit none
    integer :: z0,z1,counts,different,j,n
    integer :: tally(0:n_qubits-1)


    tally=-1
    do z0 =0,2**n_qubits-1
        do z1 =z0+1,2**n_qubits-1
            counts=0
            do n = 0,n_qubits-1
                if (basis(z0,n) .ne. basis(z1,n)) then
                    counts=counts+1
                    different=n
                endif
            enddo
            if (counts .eq. 1) then
            	tally(different) = tally(different)+1
                pairs_list(different,tally(different),0) = z0
                pairs_list(different,tally(different),1) = z1
                !print*, 'different,z0,z1,basis(z0,:),basis(z1,:)'
                !print*, different,z0,z1
                !print*, (basis(z0,j),j=1,n_qubits)
               	!print*, (basis(z1,j),j=1,n_qubits)
            endif
        enddo
    enddo

end subroutine Calc_pair_list


end module QAOA_subroutines