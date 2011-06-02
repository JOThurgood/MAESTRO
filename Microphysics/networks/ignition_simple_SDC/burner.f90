module burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module
  use network

  private
  public :: burner



contains

  subroutine burner(dens, temp, Xin, dt, Xout, rho_omegadot, rho_Hnuc, sdc_rho, sdc_X, p0)

    ! outputs:
    !   Xout are the mass fractions after burning through timestep dt
    !   rho_omegadot = rho dX/dt
    !   rho_Hnuc = - sum_k q_k rho_omegadot_k  [erg / cm^3 / s]

    use burner_aux_module, only : sdc_rho_pass, sdc_X_pass, p0_pass

    implicit none

    real(kind=dp_t), intent(inout) :: dens
    real(kind=dp_t), intent(in   ) :: temp, Xin(nspec), dt
    real(kind=dp_t), intent(  out) :: Xout(nspec), rho_omegadot(nspec), rho_Hnuc
    real(kind=dp_t), intent(in   ) :: sdc_rho, sdc_X(nspec)
    real(kind=dp_t), intent(in   ) :: p0

    integer :: n
    real(kind=dp_t) :: enuc, dX

    logical, parameter :: verbose = .false.

    ! set the number of independent variables -- this should be density
    ! + the number of species
    integer, parameter :: NEQ = 1 + nspec
  

    ! allocate storage for the input state
    real(kind=dp_t), dimension(NEQ) :: y


    ! we will always refer to the species by integer indices that come from
    ! the network module -- this makes things robust to a shuffling of the 
    ! species ordering
    integer, save :: ic12, io16, img24

    ! our problem is stiff, tell ODEPACK that. 21 means stiff, jacobian 
    ! function is supplied, 22 means stiff, figure out my jacobian through 
    ! differencing
    integer, parameter :: MF_ANALYTIC_JAC = 21, MF_NUMERICAL_JAC = 22


    ! tolerance parameters:
    !
    !  itol specifies whether to use an single absolute tolerance for
    !  all variables (1), or to pass an array of absolute tolerances, one
    !  for each variable with a scalar relative tol (2), a scalar absolute
    !  and array of relative tolerances (3), or arrays for both (4)
    !  
    !  The error is determined as e(i) = rtol*abs(y(i)) + atol, and must
    !  be > 0.  Since we have some compositions that may be 0 initially,
    !  we will specify both an absolute and a relative tolerance.
    !
    ! We will use arrays for both the absolute and relative tolerances, 
    ! since we want to be easier on the temperature than the species
    integer, parameter :: ITOL = 4
    real(kind=dp_t), dimension(NEQ) :: atol, rtol


    real(kind=dp_t) :: time
    

    ! we want to do a normal computation, and get the output values of y(t)
    ! after stepping though dt
    integer, PARAMETER :: ITASK = 1
  

    ! istate determines the state of the calculation.  A value of 1 meeans
    ! this is the first call to the problem -- this is what we will want.
    ! Note, istate is changed over the course of the calculation, so it
    ! cannot be a parameter
    integer :: istate


    ! we will override the maximum number of steps, so turn on the 
    ! optional arguments flag
    integer, parameter :: IOPT = 1
    
    ! declare a real work array of size 22 + 9*NEQ + 2*NEQ**2 and an
    ! integer work array of since 30 + NEQ
    integer, parameter :: LRW = 22 + 9*NEQ + 2*NEQ**2
    real(kind=dp_t), dimension(LRW) :: rwork
    
    integer, parameter :: LIW = 30 + NEQ
    integer, dimension(LIW) :: iwork
    

    real(kind=dp_t) :: rpar
    integer :: ipar

    EXTERNAL jac, f_rhs
    
    logical, save :: firstCall = .true.

    if (firstCall) then

       if (.NOT. network_initialized) then
          call bl_error("ERROR in burner: must initialize network first")
       endif
     
       ic12 = network_species_index("carbon-12")
       io16 = network_species_index("oxygen-16")
       img24 = network_species_index("magnesium-24")
       
       if (ic12 < 0 .OR. io16 < 0 .OR. img24 < 0) then
          call bl_error("ERROR in burner: species undefined")
       endif
       
       firstCall = .false.
    endif

    ! set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.  
    atol(1:nspec) = 1.d-12    ! mass fractions
    atol(nspec+1) = 1.d-8     ! rho
       
    rtol(1:nspec) = 1.d-12    ! mass fractions
    rtol(nspec+1) = 1.d-5     ! rho
    

    ! we want VODE to re-initialize each time we call it
    istate = 1
    
    rwork(:) = ZERO
    iwork(:) = 0
    
    
    ! set the maximum number of steps allowed (the VODE default is 500)
    iwork(6) = 15000
    
    
    ! initialize the integration time
    time = ZERO
    
    
    ! abundances are the first nspec values and density is the last
    y(ic12) = Xin(ic12)
    y(io16) = Xin(io16)
    y(img24) = Xin(img24)
    y(nspec+1) = dens

    ! sdc source terms (sdc_rho and sdc_X) are needed
    ! in the righthand side routine, so we will pass these in through the
    ! burner_aux module.
    !
    ! Since we are only integrating C12, we will need the O16 mass fraction
    ! in the RHS routine to compute the screening (and we know that the
    ! Mg24 abundance is constraint so things add to 1).

    sdc_rho_pass = sdc_rho
    sdc_X_pass(:) = sdc_X(:)
    p0_pass = p0

    ! call the integration routine
    call dvode(f_rhs, NEQ, y, time, dt, ITOL, rtol, atol, ITASK, &
               istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_NUMERICAL_JAC, &
               rpar, ipar)

    if (istate < 0) then
       print *, 'ERROR: integration failed in net'
       print *, 'istate = ', istate
       print *, 'time = ', time
       call bl_error("ERROR in burner: integration failed")
    endif


    ! store the new mass fractions -- note, we discard the temperature
    ! here and instead compute the energy release from the binding
    ! energy -- make sure that they are positive
    Xout(ic12)  = max(y(ic12), ZERO)
    Xout(io16)  = Xin(io16)
    Xout(img24) = ONE - Xout(ic12) - Xout(io16)
        
    dens = y(nspec+1)

    ! compute the energy release.  Our convention is that the binding 
    ! energies are negative, so the energy release is
    ! - sum_k { (Xout(k) - Xin(k)) ebin(k) }
    !
    ! also compute the density-weighted creation rates, rho_omegadot
    enuc = 0.0_dp_t
    do n = 1, nspec
       dX = Xout(n) - Xin(n) - dt * sdc_X(n)
       enuc = enuc - ebin(n) * dX
       rho_omegadot(n) = dens * dX / dt
    enddo

    rho_Hnuc = dens*enuc/dt

    if (verbose) then
       
       ! print out some integration statistics, if desired
       print *, 'integration summary: '
       print *, 'dens: ', dens, ' temp: ', temp
       print *, 'number of steps taken: ', iwork(11)
       print *, 'number of f evaluations: ', iwork(12)
    endif

  end subroutine burner

end module burner_module