module init_scalar_module

  use multifab_module
  use ml_layout_module
  use bl_constants_module
  use ml_restrict_fill_module
  use variables
  use network
  use fill_3d_module, only: put_1d_array_on_cart_3d_sphr

  implicit none

  private
  public :: initscalardata, initscalardata_on_level

contains

  subroutine initscalardata(s,s0_init,p0_init,dx,the_bc_level,mla)

    use geometry, only: spherical

    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: sop(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim),ng
    integer :: i,n,dm,nlevs
    
    dm = mla%dim
    nlevs = mla%nlevel

    ng = nghost(s(1))

    do n=1,nlevs
      do i = 1, nfabs(s(n))
        sop => dataptr(s(n),i)
        lo =  lwb(get_box(s(n),i))
        hi =  upb(get_box(s(n),i))

        select case (dm)
        case (1)
          call bl_error("ERROR: initscalardata not support in 1d")
        case (2)
          call initscalardata_2d(sop(:,:,1,:), lo, hi, ng, dx(n,:), s0_init(n,:,:), &
                                 p0_init(n,:))
        case (3)
          if (spherical .eq. 1) then
            print *, 'spherical not allowed for circular_drop'
            print *, 'STOP from initscaldata.f90'
            STOP
          else
            print *, '3d not implemented circular_drop'
            print *, 'STOP from initscaldata.f90'
            STOP
           end if
        end select
      end do
    enddo

    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,s,mla%mba%rr,the_bc_level, &
                              icomp=rho_comp, &
                              bcomp=dm+rho_comp, &
                              nc=nscal, &
                              ng=s(1)%ng)

  end subroutine initscalardata

  subroutine initscalardata_on_level(n,s,s0_init,p0_init,dx,the_bc_level)

    use geometry, only: spherical

    integer        , intent(in   ) :: n
    type(multifab) , intent(inout) :: s
    real(kind=dp_t), intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_level) , intent(in   ) :: the_bc_level

    ! local
    integer                  :: ng,i,dm
    integer                  :: lo(get_dim(s)),hi(get_dim(s))
    real(kind=dp_t), pointer :: sop(:,:,:,:)

    dm = get_dim(s)

    ng = nghost(s)

    do i = 1, nfabs(s)
       sop => dataptr(s,i)
       lo =  lwb(get_box(s,i))
       hi =  upb(get_box(s,i))
       select case (dm)
       case (1)
          call bl_error("ERROR: initscalardata not supported in 1d")
       case (2)
          call initscalardata_2d(sop(:,:,1,:),lo,hi,ng,dx,s0_init,p0_init)
       case (3)
          if (spherical .eq. 1) then
            print *, 'spherical not allowed for circular_drop'
            print *, 'STOP from initscaldata.f90'
            STOP
          else
            call initscalardata_3d(sop(:,:,:,:),lo,hi,ng,dx,s0_init,p0_init)
            print *, '3d not implemented circular_drop'
            print *, 'STOP from initscaldata.f90'
            STOP
          end if
       end select
    end do

  end subroutine initscalardata_on_level

  subroutine initscalardata_2d(s,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only: prob_lo, radius, rho_2
    use eos_module, only: eos_input_rp, eos
    use eos_type_module

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    ! Local variables
    integer         :: i,j
    real(kind=dp_t) :: x,y, rloc
    type (eos_t) :: eos_state

    ! initial the domain with the base state
    s = ZERO

    ! initialize the whole 2d field to the 1d background
    do j = lo(2), hi(2)
      do i = lo(1), hi(1)
        s(i,j,rho_comp)  = s0_init(j,rho_comp)
        s(i,j,rhoh_comp) = s0_init(j,rhoh_comp)
        s(i,j,temp_comp) = s0_init(j,temp_comp)
        s(i,j,spec_comp:spec_comp+nspec-1) = &
             s0_init(j,spec_comp:spec_comp+nspec-1)
        s(i,j,trac_comp:trac_comp+ntrac-1) = &
             s0_init(j,trac_comp:trac_comp+ntrac-1)
      enddo
    enddo

    ! where necessary, change state to be different
    ! density and set rest of variables for 
    ! thermodynamic consistency

    do j = lo(2), hi(2)
      do i = lo(1), hi(1)

        x = prob_lo(1) + (dble(i)+HALF) * dx(1)
        y = prob_lo(2) + (dble(j)+HALF) * dx(2)
        x = x - ( prob_lo(1) + 0.5d0  )
        y = y - ( prob_lo(2) + 0.75d0 )
        rloc = sqrt(x**2 + y**2)

        if (rloc < radius ) then

          eos_state%rho   = rho_2 ! should replace with rt param 
          eos_state%p     = p0_init(j)
          eos_state%T     = s0_init(j,temp_comp) ! not necessarily consistent yet?
          eos_state%xn(:) = ONE ! single fluid
 
          ! (rho,p) --> T, h
          call eos(eos_input_rp, eos_state)

          s(i,j,rho_comp)  = eos_state%rho
          s(i,j,rhoh_comp) = eos_state%rho * eos_state%h
          s(i,j,temp_comp) = eos_state%T
          s(i,j,spec_comp:spec_comp+nspec-1) = eos_state%rho * eos_state%xn ! single fluid
        
          ! (do nothing to trac_comp)

        endif
      enddo
    enddo

  end subroutine initscalardata_2d

  subroutine initscalardata_3d(s,lo,hi,ng,dx,s0_init,p0_init)

    use init_perturb_module
    
    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    s = ZERO
    
  end subroutine initscalardata_3d

end module init_scalar_module