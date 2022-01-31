module QAOA_subroutines

implicit none

contains

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

    do m = 0,n_edges-1
    	!check the edge exists in edges array
    	if (edges(m,0) .ne. edges(m,1)) then 
    		!if (m .eq. 1) print*, '*',edges(m,0),edges(m,1)
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
    		!print*, 'edge term:',0.5d0*(1.d0 + term0 + term1)
    		C_ma_trianglefree = C_ma_trianglefree + 0.5d0*(1.d0 + term0 + term1)
    	endif
    enddo

    C_ma_trianglefree=-C_ma_trianglefree/dble(n_qubits)

end function C_ma_trianglefree

function C_ma_trianglefree_rescaleangles(n_angles,angles)
    use parameters
    implicit none
    integer :: n_angles
    double precision :: angles(n_angles)
    integer :: m,l
    double precision :: c_ma_trianglefree_rescaleangles
    double precision :: term0,term1,t_0,t_1
    integer :: q0,q1,edge_angle_index

   ! call cpu_time(t_0)
    
    angles = angles/angle_rescaling_fac
    c_ma_trianglefree_rescaleangles = 0.d0
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
    		c_ma_trianglefree_rescaleangles = c_ma_trianglefree_rescaleangles + 0.5d0*(1.d0 + term0 + term1)
    	endif
    enddo

    c_ma_trianglefree_rescaleangles=-c_ma_trianglefree_rescaleangles/C_rescaling_fac
 	angles = angles*angle_rescaling_fac
    !print*, 'c_ma_trianglefree_rescaleangles',c_ma_trianglefree_rescaleangles
    !stop
    !call cpu_time(t_1)
    !print*, t_1-t_0
end function C_ma_trianglefree_rescaleangles



subroutine analytic_grad_C_ma_trianglefree_rescaleangles(n_angles,angles,grad,func)
    use parameters
    implicit  none
    integer :: n_angles
    double precision :: angles(n_angles)
    integer :: j,m,l,counts(n_angles,0:1)
    double precision :: c_ma_trianglefree_rescaleangles
    double precision :: term0,term1,t_0,t_1
    double precision :: term0_array(n_angles),term1_array(n_angles)
    integer :: q0,q1,edge_angle_index

    double precision :: grad(n_angles)

    External func

    call cpu_time(t_0)
    grad=0.d0
    angles = angles/angle_rescaling_fac
    do m = 0,n_edges-1
    	!check the edge exists in edges array
    	if (edges(m,0) .ne. edges(m,1)) then 
    		q0 = edges(m,0)+1
    		q1 = edges(m,1)+1
    		edge_angle_index=n_beta_angles+m+1
    		counts=0
    		term1_array=0.d0
    		term0_array=0.d0
    		do j = 1,n_angles
    			if (j .eq. q0) then
    				term1_array(j)=&
    				&dsin(angles(edge_angle_index))*(-2.d0)*dsin(angles(q0)*2.d0)*dsin(angles(q1)*2.d0)
    				term0_array(j)=&
    				&dsin(angles(edge_angle_index))*dcos(angles(q1)*2.d0)*(2.d0)*dcos(angles(q0)*2.d0)
    				counts(j,:) = counts(j,:) +1
    				!if (j .eq. 1) print*, 'q0,term1_array(j),term0_array(j)',q0,term1_array(j),term0_array(j)
    			else if (j .eq. q1) then
    				term1_array(j)=&
    				&dsin(angles(edge_angle_index))*dcos(angles(q0)*2.d0)*(2.d0)*dcos(angles(q1)*2.d0)
   					term0_array(j)=&
   					&dsin(angles(edge_angle_index))*(-2.d0)*dsin(angles(q1)*2.d0)*dsin(angles(q0)*2.d0)
   					counts(j,:) = counts(j,:) +1
   					!print*, 'q1,term1_array(j),term0_array(j)',q1,term1_array(j),term0_array(j)
   					!print*, '*****'
   				else if (j .eq. edge_angle_index) then
   					term1_array(j)=dcos(angles(edge_angle_index))*dcos(angles(q0)*2.d0)*dsin(angles(q1)*2.d0)
   					term0_array(j)=dcos(angles(edge_angle_index))*dcos(angles(q1)*2.d0)*dsin(angles(q0)*2.d0)
   					counts(j,:) = counts(j,:) +1
   				else
    				term1_array(j)=dsin(angles(edge_angle_index))*dcos(angles(q0)*2.d0)*dsin(angles(q1)*2.d0)
    				term0_array(j)=dsin(angles(edge_angle_index))*dcos(angles(q1)*2.d0)*dsin(angles(q0)*2.d0)
    			endif
    		enddo
    		do l = 0,n_edges-1
    			if (  (l .ne. m) .and. ( edges(l,0) .ne. edges(l,1) )  ) then
    				if ( (edges(l,0) .eq. q0-1) .or. (edges(l,1) .eq. q0-1) ) then
    					edge_angle_index = n_beta_angles+l+1
    					!term0 = term0*dcos(angles(edge_angle_index))
    					do j =1,n_angles
    						if (j .eq. edge_angle_index) then
    							term0_array(j) = term0_array(j)*(-1.d0)*dsin(angles(edge_angle_index))
    							counts(j,0)=counts(j,0)+1
    						else
    							term0_array(j)=term0_array(j)*dcos(angles(edge_angle_index))
    						endif
    					enddo
    				endif
    				if ( (edges(l,0) .eq. q1-1) .or. (edges(l,1) .eq. q1-1) ) then
    					edge_angle_index = n_beta_angles+l+1
    					!term1 = term1*dcos(angles(edge_angle_index))
    					do j=1,n_angles
    						if (j .eq. edge_angle_index) then
    							term1_array(j) = term1_array(j)*(-1.d0)*dsin(angles(edge_angle_index))
    							counts(j,1)=counts(j,1)+1
    						else
    							term1_array(j) = term1_array(j)*dcos(angles(edge_angle_index))
    						endif
    					enddo
    				endif
    			endif
    		enddo

    		do j = 1,n_angles
    			!print*, 'j,counts(j,0),counts(j,1):',j,counts(j,0),counts(j,1)
    			if ( counts(j,0) .eq. 1 ) then
    				grad(j) = grad(j) + 0.5d0*term0_array(j)
    				!print*, 'j,grad(j),term0_array(j)',j,grad(j),term0_array(j)
    			endif
    			if ( counts(j,1) .eq. 1 ) then
    				grad(j) = grad(j) + 0.5d0*term1_array(j)
    				!print*, 'j,grad(j),term1_array(j)',j,grad(j),term1_array(j)
    			endif
    		enddo

    	endif
    enddo

    grad=-grad/C_rescaling_fac
    angles = angles*angle_rescaling_fac
    !print*, 'grad:',grad
   ! stop
    call cpu_time(t_1)
    print*, '*',t_1-t_0
end subroutine analytic_grad_C_ma_trianglefree_rescaleangles

subroutine analytic_grad_C_ma_trianglefree_rescaleangles2(n_angles,angles,grad)
    use parameters
    implicit  none
    integer :: n_angles
    double precision :: angles(n_angles)
    integer :: j,m,l,counts(n_angles,0:1)
    double precision :: c_ma_trianglefree_rescaleangles
    double precision :: term0,term1
    double precision :: term0_array(n_angles),term1_array(n_angles)
    integer :: q0,q1,edge_angle_index

    double precision :: grad(n_angles)

    External func

    grad=0.d0
    angles = angles/angle_rescaling_fac
    !stop
    do m = 0,n_edges-1
    	!check the edge exists in edges array
    	if (edges(m,0) .ne. edges(m,1)) then 
    		q0 = edges(m,0)+1
    		q1 = edges(m,1)+1
    		edge_angle_index=n_beta_angles+m+1
    		counts=0
    		term1_array=0.d0
    		term0_array=0.d0
    		do j = 1,n_angles
    			if (j .eq. q0) then
    				term1_array(j)=&
    				&dsin(angles(edge_angle_index))*(-2.d0)*dsin(angles(q0)*2.d0)*dsin(angles(q1)*2.d0)
    				term0_array(j)=&
    				&dsin(angles(edge_angle_index))*dcos(angles(q1)*2.d0)*(2.d0)*dcos(angles(q0)*2.d0)
    				counts(j,:) = counts(j,:) +1
    				!if (j .eq. 1) print*, 'q0,term1_array(j),term0_array(j)',q0,term1_array(j),term0_array(j)
    			else if (j .eq. q1) then
    				term1_array(j)=&
    				&dsin(angles(edge_angle_index))*dcos(angles(q0)*2.d0)*(2.d0)*dcos(angles(q1)*2.d0)
   					term0_array(j)=&
   					&dsin(angles(edge_angle_index))*(-2.d0)*dsin(angles(q1)*2.d0)*dsin(angles(q0)*2.d0)
   					counts(j,:) = counts(j,:) +1
   					!print*, 'q1,term1_array(j),term0_array(j)',q1,term1_array(j),term0_array(j)
   					!print*, '*****'
   				else if (j .eq. edge_angle_index) then
   					term1_array(j)=dcos(angles(edge_angle_index))*dcos(angles(q0)*2.d0)*dsin(angles(q1)*2.d0)
   					term0_array(j)=dcos(angles(edge_angle_index))*dcos(angles(q1)*2.d0)*dsin(angles(q0)*2.d0)
   					counts(j,:) = counts(j,:) +1
   				else
    				term1_array(j)=dsin(angles(edge_angle_index))*dcos(angles(q0)*2.d0)*dsin(angles(q1)*2.d0)
    				term0_array(j)=dsin(angles(edge_angle_index))*dcos(angles(q1)*2.d0)*dsin(angles(q0)*2.d0)
    			endif
    		enddo
    		do l = 0,n_edges-1
    			if (  (l .ne. m) .and. ( edges(l,0) .ne. edges(l,1) )  ) then
    				if ( (edges(l,0) .eq. q0-1) .or. (edges(l,1) .eq. q0-1) ) then
    					edge_angle_index = n_beta_angles+l+1
    					!term0 = term0*dcos(angles(edge_angle_index))
    					do j =1,n_angles
    						if (j .eq. edge_angle_index) then
    							term0_array(j) = term0_array(j)*(-1.d0)*dsin(angles(edge_angle_index))
    							counts(j,0)=counts(j,0)+1
    						else
    							term0_array(j)=term0_array(j)*dcos(angles(edge_angle_index))
    						endif
    					enddo
    				endif
    				if ( (edges(l,0) .eq. q1-1) .or. (edges(l,1) .eq. q1-1) ) then
    					edge_angle_index = n_beta_angles+l+1
    					!term1 = term1*dcos(angles(edge_angle_index))
    					do j=1,n_angles
    						if (j .eq. edge_angle_index) then
    							term1_array(j) = term1_array(j)*(-1.d0)*dsin(angles(edge_angle_index))
    							counts(j,1)=counts(j,1)+1
    						else
    							term1_array(j) = term1_array(j)*dcos(angles(edge_angle_index))
    						endif
    					enddo
    				endif
    			endif
    		enddo

    		do j = 1,n_angles
    			!print*, 'j,counts(j,0),counts(j,1):',j,counts(j,0),counts(j,1)
    			if ( counts(j,0) .eq. 1 ) then
    				grad(j) = grad(j) + 0.5d0*term0_array(j)
    				!print*, 'j,grad(j),term0_array(j)',j,grad(j),term0_array(j)
    			endif
    			if ( counts(j,1) .eq. 1 ) then
    				grad(j) = grad(j) + 0.5d0*term1_array(j)
    				!print*, 'j,grad(j),term1_array(j)',j,grad(j),term1_array(j)
    			endif
    		enddo

    	endif
    enddo

    grad=-grad/C_rescaling_fac
    angles = angles*angle_rescaling_fac
    !print*, 'grad:',grad
   ! stop
end subroutine analytic_grad_C_ma_trianglefree_rescaleangles2

subroutine C_ma_trianglefree_rescaleangles_nlopt(val, n_angles, angles, grad, need_gradient, f_data)

    use parameters 
    implicit none

    integer :: n_angles,need_gradient,f_data
    double precision :: angles(n_angles),grad(n_angles),val

    !print*, 'angles nlopt:',angles

    val=C_ma_trianglefree_rescaleangles(n_angles,angles)

    if (need_gradient .ne. 0) then
    	call analytic_grad_C_ma_trianglefree_rescaleangles2(n_angles,angles,grad)
    	!call dfunc(n_angles,angles,grad,C_ma_trianglefree_rescaleangles)
    	!print*, 'grad:',grad
    endif

end subroutine C_ma_trianglefree_rescaleangles_nlopt


subroutine C_trianglefree_nlopt(val, n_angles, angles, grad, need_gradient, f_data)

    use parameters 
    implicit none

    integer :: n_angles,need_gradient,f_data
    double precision :: angles(n_angles),grad(n_angles),val

    !print*, 'angles nlopt:',angles

    val=C_trianglefree(n_angles,angles)

    if (need_gradient .ne. 0) then
    	!call analytic_grad_C_ma_trianglefree_rescaleangles2(n_angles,angles,grad)
    	call dfunc(n_angles,angles,grad,C_trianglefree)
    	!print*, 'grad:',grad
    endif

end subroutine C_trianglefree_nlopt


function C_trianglefree(n_angles,angles)
    use parameters
    implicit none
    integer :: n_angles
    double precision :: angles(n_angles)
    integer :: m,l
    double precision :: C_trianglefree
    double precision :: term0,term1
    integer :: q0,q1,edge_angle_index

    C_trianglefree = 0.d0

    do m = 0,n_edges-1
    	!check the edge exists in edges array
    	if (edges(m,0) .ne. edges(m,1)) then 
    		q0 = edges(m,0)+1
    		q1 = edges(m,1)+1
    		!edge_angle_index=n_beta_angles+m+1
    		term1=dsin(angles(2))*dcos(angles(1)*2.d0)*dsin(angles(1)*2.d0)
    		term0=dsin(angles(2))*dcos(angles(1)*2.d0)*dsin(angles(1)*2.d0)
    		do l = 0,n_edges-1
    			if (  (l .ne. m) .and. ( edges(l,0) .ne. edges(l,1) )  ) then
    				if ( (edges(l,0) .eq. q0-1) .or. (edges(l,1) .eq. q0-1) ) then
    					edge_angle_index = n_beta_angles+l+1
    					term0 = term0*dcos(angles(2))
    				endif
    				if ( (edges(l,0) .eq. q1-1) .or. (edges(l,1) .eq. q1-1) ) then
    					edge_angle_index = n_beta_angles+l+1
    					term1 = term1*dcos(angles(2))
    				endif
    			endif
    		enddo
    		C_trianglefree = C_trianglefree + 0.5d0*(1.d0 + term0 + term1)
    	endif
    enddo

    C_trianglefree=-C_trianglefree/C_rescaling_fac

end function C_trianglefree


subroutine dfunc(n,x,grad,func)

	implicit none

	integer :: n,i
	double precision :: x(n), x_dif_plus(n), x_dif_min(n)
	double precision :: nominal_delta = 1.d-5
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
 		!print*, 'grad(i):',grad(i)
	enddo

	!stop

end subroutine dfunc


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
			!print*, 'line:',line
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


subroutine Generate_angles(n_angles,angles)

	use parameters
	implicit none

	integer :: i
	integer :: n_angles
	double precision :: junk(6)
	double precision :: angles(n_angles)
	double precision :: tmp

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

end subroutine Generate_angles


end module QAOA_subroutines