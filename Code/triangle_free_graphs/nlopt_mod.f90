module nlopt

contains

subroutine optimize_nlopt(func,dim,parameters,minf)
      external func
      integer :: dim,graph_number,weight_number,i
      double precision :: parameters(dim),lb(dim),ub(dim)
      integer*8 :: opt
      double precision :: minf
      integer :: ires, f_data
      include 'nlopt.f'
      f_data=0
      opt=0

      !optimization algorithms:  https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/
      !call nlo_create(opt, NLOPT_LD_LBFGS, dim)
      call nlo_create(opt, NLOPT_LD_MMA, dim)

      !minimize func
      call nlo_set_min_objective(ires, opt, func, f_data)
      if (ires .lt. 0) print*, 'fail! nlo_set_max_objective'

      !bounds +/- infinity
      call nlo_get_lower_bounds(ires, opt, lb)
      call nlo_get_upper_bounds(ires, opt, ub)

      !tolerance
      call nlo_set_xtol_rel(ires, opt, 1.d-7)
      if (ires .lt. 0) print*, 'fail! nlo_set_xtol_rel'
      call nlo_set_ftol_rel(ires, opt, 1.d-7)
      if (ires .lt. 0) print*, 'fail! nlo_set_ftol_rel'

      !optimize
      call nlo_optimize(ires, opt, parameters, minf)
      if (ires.lt.0) then
        write(*,*) 'nlopt failed!', ires, parameters, minf
      endif

      call nlo_destroy(opt)

end subroutine

end module nlopt
