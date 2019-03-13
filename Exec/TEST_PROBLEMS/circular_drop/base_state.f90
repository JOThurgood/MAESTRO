module base_state_module

  use bl_types
  use network, only: nspec
  
  implicit none

  private

  public :: init_base_state

contains

  subroutine init_base_state(n,model_file,s0_init,p0_init,dx)

    use bl_prof_module
    use parallel
    use bl_error_module
    use bl_constants_module
    use eos_module, only: eos_input_rp, eos
    use eos_type_module
    use probin_module, only : rho_1, p0_base, do_smallscale, &
                              prob_lo, grav_const
 
    use variables, only: rho_comp, rhoh_comp, temp_comp, spec_comp, trac_comp, ntrac
    use geometry, only: dr, spherical, nr
    use inlet_bc_module, only: set_inlet_bcs
    
    integer           , intent(in   ) :: n
    character(len=256), intent(in   ) :: model_file
    real(kind=dp_t)   , intent(inout) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(inout) :: p0_init(0:)
    real(kind=dp_t)   , intent(in   ) :: dx(:)

    ! local
    integer         :: r
    real(kind=dp_t) :: rloc
    real(kind=dp_t) :: d_ambient,t_ambient,p_ambient,xn_ambient(nspec)
    real(kind=dp_t) :: t_guess
    real(kind=dp_t), parameter :: TINY = 1.0d-10
    real(kind=dp_t), parameter :: SMALL = 1.d-12

    type (eos_t) :: eos_state

    if (spherical .eq. 1) then
       call bl_error("ERROR: circular_drop base_state is not valid for spherical")
    endif

    ! set a guess for the temperature for the EOS calls
    t_guess = 1.e-8

    do r = 0, nr(n)-1
      ! height above the bottom of the domain
      rloc = (dble(r) + HALF)*dr(n)
 
      d_ambient = rho_1
      if (do_smallscale) then
         p_ambient = p0_base
      else
         p_ambient = p0_base + rho_1*grav_const*(rloc - prob_lo(size(dx)))
      end if
      t_ambient = t_guess
      xn_ambient(:) = ONE-SMALL
      ! not sure if -SMALL is necessary - monkey see monkey do

      ! use the EOS to make the state consistent
      eos_state%T     = t_ambient
      eos_state%rho   = d_ambient
      eos_state%p     = p_ambient
      eos_state%xn(:) = xn_ambient(:)

      ! (rho,p) --> T, h
      call eos(eos_input_rp, eos_state)

      s0_init(r, rho_comp) = d_ambient
      s0_init(r,rhoh_comp) = d_ambient * eos_state%h
      s0_init(r,spec_comp:spec_comp+nspec-1) = d_ambient * xn_ambient(1:nspec)
      if (do_smallscale) then
         p0_init(r) = p0_base
      else
         p0_init(r) = eos_state%p
      end if
      
      s0_init(r,temp_comp) = eos_state%T

      if (ntrac .gt. 0) then
         s0_init(r,trac_comp:trac_comp+ntrac-1) = ZERO
      end if

    end do ! r = 0, rn(n)-1

    ! initialize any inlet BC parameters
    call set_inlet_bcs()

  end subroutine init_base_state

end module base_state_module
