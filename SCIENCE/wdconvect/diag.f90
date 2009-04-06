! This is a general interface for doing runtime diagnostics on the state.
! It is called at the end of advance

module diag_module

  use bl_types
  use bl_IO_module
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: diag

contains

  subroutine diag(time,dt,dx,s,rho_Hnuc,rho_Hext, &
                  rho0,rhoh0,p0,tempbar,gamma1bar,div_coeff, &
                  u,w0,normal, &
                  mla,the_bc_tower)

    use bl_prof_module
    use geometry, only: dm, nlevs, spherical
    use bl_constants_module
    use variables, only: foextrap_comp
    use fill_3d_module
    use probin_module, only: prob_lo_x, prob_lo_y, prob_lo_z, &
                             prob_hi_x, prob_hi_y, prob_hi_z, &
                             job_name,edge_nodal_flag

    real(kind=dp_t), intent(in   ) :: dt,dx(:,:),time
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(in   ) :: rho_Hnuc(:)
    type(multifab) , intent(in   ) :: rho_Hext(:)    
    type(multifab) , intent(in   ) :: u(:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) ::      rho0(:,0:)
    real(kind=dp_t), intent(in   ) ::     rhoh0(:,0:)
    real(kind=dp_t), intent(in   ) ::        p0(:,0:)
    real(kind=dp_t), intent(in   ) ::   tempbar(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(in   ) :: div_coeff(:,0:)
    real(kind=dp_t), intent(in   ) ::        w0(:,0:)
    type(ml_layout), intent(in   ) :: mla
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! Local
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: rhnp(:,:,:,:)
    real(kind=dp_t), pointer :: rhep(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)
    real(kind=dp_t), pointer :: w0rp(:,:,:,:)
    real(kind=dp_t), pointer :: w0xp(:,:,:,:)
    real(kind=dp_t), pointer :: w0yp(:,:,:,:)
    real(kind=dp_t), pointer :: w0zp(:,:,:,:)
    logical,         pointer :: mp(:,:,:,:)

    type(multifab) :: w0r_cart(mla%nlevel)
    type(multifab) ::    w0mac(mla%nlevel,dm)

    real(kind=dp_t) :: vr(dm),    vr_level(dm),    vr_local(dm)
    real(kind=dp_t) :: vr_max,    vr_max_level,    vr_max_local
    real(kind=dp_t) :: rhovr(dm), rhovr_level(dm), rhovr_local(dm)
    real(kind=dp_t) :: mass,      mass_level,      mass_local
    real(kind=dp_t) :: nzones,    nzones_level,    nzones_local
    real(kind=dp_t) :: T_max,     T_max_level,     T_max_local
    real(kind=dp_t) :: enuc_max,  enuc_max_level,  enuc_max_local
    real(kind=dp_t) :: kin_ener,  kin_ener_level,  kin_ener_local
    real(kind=dp_t) :: U_max,     U_max_level,     U_max_local

    real(kind=dp_t) :: vr_favre(dm)

    integer :: lo(dm),hi(dm)
    integer :: ng_s,ng_u,ng_n,ng_w,ng_wm,ng_rhn,ng_rhe
    integer :: i,n, comp
    integer :: un,un2,un3,un4
    logical :: lexist

    logical, save :: firstCall = .true.

    type(bl_prof_timer), save :: bpt

    call build(bpt, "diagnostics")

    if (spherical .eq. 1) then

       do n=1,nlevs

          do comp=1,dm
             ! w0mac will contain an edge-centered w0 on a Cartesian grid,                             
             ! for use in computing divergences.                                                       
             call multifab_build(w0mac(n,comp), mla%la(n),1,1,nodal=edge_nodal_flag(comp,:))
             call setval(w0mac(n,comp), ZERO, all=.true.)
          enddo

          ! w0r_cart is w0 but onto a Cartesian grid in cell-centered as
          ! a scalar.  Since w0 is the radial expansion velocity, w0r_cart
          ! is the radial w0 in a zone
          call multifab_build(w0r_cart(n), mla%la(n),1,0)
          call setval(w0r_cart(n), ZERO, all=.true.)
       end do

       ! put w0 on Cartesian edges as a vector                                                         
       call put_w0_on_edges(mla,w0,w0mac,dx,div_coeff,the_bc_tower)


       ! put w0 in Cartesian cell-centers as a scalar (the radial 
       ! expansion velocity)
       call put_1d_array_on_cart(w0,w0r_cart,foextrap_comp,.true.,.false.,dx, &
                                 the_bc_tower%bc_tower_array,mla,normal=normal)

    else

       call bl_error("ERROR: wdconvect/diag.f90: geometry not spherical")

    endif

    ng_s = s(1)%ng
    ng_u = u(1)%ng
    ng_n = normal(1)%ng
    ng_w = w0r_cart(1)%ng
    ng_wm = w0mac(1,1)%ng
    ng_rhn = rho_Hnuc(1)%ng
    ng_rhe = rho_Hext(1)%ng

    !=========================================================================
    ! initialize
    !=========================================================================
    vr(:)    = ZERO
    rhovr(:) = ZERO
    mass     = ZERO
    nzones   = ZERO
    vr_max   = ZERO
    T_max    = ZERO
    enuc_max = ZERO
    kin_ener = ZERO
    U_max    = ZERO


    !=========================================================================
    ! loop over the levels and compute the global quantities
    !=========================================================================
    do n = 1, nlevs

       vr_level(:) = ZERO
       vr_local(:) = ZERO
       
       rhovr_level(:) = ZERO
       rhovr_local(:) = ZERO
       
       mass_level = ZERO
       mass_local = ZERO
       
       nzones_level = ZERO
       nzones_local = ZERO
       
       vr_max_level = ZERO
       vr_max_local = ZERO
       
       T_max_level = ZERO
       T_max_local = ZERO
       
       enuc_max_level = ZERO
       enuc_max_local = ZERO

       kin_ener_level = ZERO
       kin_ener_local = ZERO

       U_max_level = ZERO
       U_max_local = ZERO


       !----------------------------------------------------------------------
       ! loop over boxes in a given level
       !----------------------------------------------------------------------
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle
          sp => dataptr(s(n) , i)
          rhnp => dataptr(rho_Hnuc(n), i)
          rhep => dataptr(rho_Hext(n), i)
          up => dataptr(u(n) , i)
          np => dataptr(normal(n) , i)
          w0rp => dataptr(w0r_cart(n), i)
          w0xp => dataptr(w0mac(n,1), i)
          w0yp => dataptr(w0mac(n,2), i)
          w0zp => dataptr(w0mac(n,3), i)

          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))

          select case (dm)
          case (2)
             call bl_error("ERROR: 2-d diag not implmented")
          case (3)
             if (n .eq. nlevs) then
                call diag_3d(n,time,dt,dx(n,:), &
                             sp(:,:,:,:),ng_s, &
                             rhnp(:,:,:,1), ng_rhn, &
                             rhep(:,:,:,1), ng_rhe, &
                             up(:,:,:,:),ng_u, &
                             w0rp(:,:,:,1), ng_w, &
                             w0xp(:,:,:,1),w0yp(:,:,:,1),w0zp(:,:,:,1),ng_wm, &                             
                             np(:,:,:,:),ng_n, &
                             lo,hi, &
                             nzones_local, &
                             vr_local(1),vr_local(2),vr_local(3),vr_max_local, &
                             rhovr_local(1), rhovr_local(2), rhovr_local(3), mass_local, &
                             T_max_local, enuc_max_local, kin_ener_local, U_max_local)
             else
                mp => dataptr(mla%mask(n), i)
                call diag_3d(n,time,dt,dx(n,:), &
                             sp(:,:,:,:),ng_s, &
                             rhnp(:,:,:,1), ng_rhn, &
                             rhep(:,:,:,1), ng_rhe, &
                             up(:,:,:,:),ng_u, &
                             w0rp(:,:,:,1), ng_w, &
                             w0xp(:,:,:,1),w0yp(:,:,:,1),w0zp(:,:,:,1),ng_wm, &                             
                             np(:,:,:,:),ng_n, &
                             lo,hi, &
                             nzones_local, &
                             vr_local(1),vr_local(2),vr_local(3),vr_max_local, &
                             rhovr_local(1), rhovr_local(2), rhovr_local(3), mass_local, &
                             T_max_local, enuc_max_local, kin_ener_local, U_max_local, &
                             mp(:,:,:,1))
             end if
          end select
       end do

       !----------------------------------------------------------------------
       ! do the appropriate parallel reduction for the current level
       !----------------------------------------------------------------------
       call parallel_reduce(vr_level, vr_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(rhovr_level, rhovr_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(mass_level, mass_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(nzones_level, nzones_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(vr_max_level, vr_max_local, MPI_MAX, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(T_max_level, T_max_local, MPI_MAX, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(enuc_max_level, enuc_max_local, MPI_MAX, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(kin_ener_level, kin_ener_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(U_max_level, U_max_local, MPI_MAX, &
                            proc = parallel_IOProcessorNode())

       vr       = vr     + vr_level
       rhovr    = rhovr  + rhovr_level
       mass     = mass   + mass_level
       nzones   = nzones + nzones_level
       vr_max   = max(vr_max,   vr_max_level)
       T_max    = max(T_max,    T_max_level)
       enuc_max = max(enuc_max, enuc_max_level)     
       kin_ener = kin_ener + kin_ener_level
       U_max    = max(U_max,    U_max_level)

    end do


    !=========================================================================
    ! normalize
    !=========================================================================
    vr(:) = vr(:)/nzones
    vr_favre(:) = rhovr(:)/mass    ! note, the common dV normalization cancels
    mass = mass*dx(1,1)*dx(1,2)*dx(1,3)
    kin_ener = kin_ener*dx(1,1)*dx(1,2)*dx(1,3)


    !=========================================================================
    ! output
    !=========================================================================
 999 format("# job name: ",a)
1000 format(1x,10(g20.10,1x))
1001 format("#",10(a20,1x))

    if (parallel_IOProcessor()) then

       ! open the diagnostic files for output, taking care not to overwrite
       ! an existing file
       un = unit_new()
       inquire(file="wdconvect_radvel_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un, file="wdconvect_radvel_diag.out", &
               status="old", position="append")
       else
          open(unit=un, file="wdconvect_radvel_diag.out", status="new")
       endif

       un2 = unit_new()
       inquire(file="wdconvect_temp_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un2, file="wdconvect_temp_diag.out", &
               status="old", position="append")
       else
          open(unit=un2, file="wdconvect_temp_diag.out", status="new")
       endif

       un3 = unit_new()
       inquire(file="wdconvect_enuc_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un3, file="wdconvect_enuc_diag.out", &
               status="old", position="append")
       else
          open(unit=un3, file="wdconvect_enuc_diag.out", status="new")
       endif

       un4 = unit_new()
       inquire(file="wdconvect_vel_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un4, file="wdconvect_vel_diag.out", &
               status="old", position="append")
       else
          open(unit=un4, file="wdconvect_vel_diag.out", status="new")
       endif


       ! write out the headers
       if (firstCall) then
          
          write (un, *) " "
          write (un, 999) trim(job_name)
          write (un, 1001) "time", "<vr_x>", "<vr_y>", "<vr_z>", "<vr>", &
                           "max{|vr|}", &
                           "int{rhovr_x}/mass", "int{rhovr_y}/mass", "int{rhovr_z}/mass", &
                           "mass"

          write (un2, *) " "
          write (un2, 999) trim(job_name)
          write (un2,1001) "time", "max{T}"

          write (un3, *) " "
          write (un3, 999) trim(job_name)
          write (un3,1001) "time", "max{enuc}"

          write (un4, *) " "
          write (un4, 999) trim(job_name)
          write (un4,1001) "time", "max{|U|}", "tot. kin. energy"

          firstCall = .false.
       endif

       ! write out the data
       write (un,1000) time, vr(1), vr(2), vr(3), &
            sqrt(vr(1)**2 + vr(2)**2 + vr(3)**2), vr_max, &
            vr_favre(1), vr_favre(2), vr_favre(3), mass
       
       write (un2,1000) time, T_max

       write (un3,1000) time, enuc_max

       write (un4,1000) time, U_max, kin_ener

       close(un)
       close(un2)
       close(un3)
       close(un4)
    endif

    if (spherical .eq. 1) then
       do n=1,nlevs
          call destroy(w0r_cart(n))
       end do
    end if

    call destroy(bpt)

  end subroutine diag


  !===========================================================================
  ! the actual diagnostic routine
  !===========================================================================
  subroutine diag_3d(n,time,dt,dx, &
                     s,ng_s, &
                     rho_Hnuc,ng_rhn, &
                     rho_Hext,ng_rhe, &
                     u,ng_u, &
                     w0r,ng_w, &
                     w0macx,w0macy,w0macz,ng_wm, &
                     normal,ng_n, &
                     lo,hi, &
                     nzones, &
                     vr_x,vr_y,vr_z,vr_max, &
                     rhovr_x,rhovr_y,rhovr_z,mass, &
                     T_max,enuc_max,kin_ener,U_max, &
                     mask)

    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp
    use bl_constants_module
    use network, only: nspec
    use geometry, only: spherical
    use probin_module, only: base_cutoff_density

    integer,          intent(in   ) :: n,lo(:),hi(:),ng_s,ng_u,ng_n,ng_w,ng_wm,ng_rhn,ng_rhe
    real (kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:  ,lo(2)-ng_s:  ,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho_Hnuc(lo(1)-ng_rhn:,lo(2)-ng_rhn:,lo(3)-ng_rhn:)
    real (kind=dp_t), intent(in   ) :: rho_Hext(lo(1)-ng_rhe:,lo(2)-ng_rhe:,lo(3)-ng_rhe:)
    real (kind=dp_t), intent(in   ) ::        u(lo(1)-ng_u:  ,lo(2)-ng_u:  ,lo(3)-ng_u:,:)
    real (kind=dp_t), intent(in   ) ::      w0r(lo(1)-ng_w:  ,lo(2)-ng_w:  ,lo(3)-ng_w:)
    real (kind=dp_t), intent(in   ) ::   w0macx(lo(1)-ng_wm: ,lo(2)-ng_wm: ,lo(3)-ng_wm:)
    real (kind=dp_t), intent(in   ) ::   w0macy(lo(1)-ng_wm: ,lo(2)-ng_wm: ,lo(3)-ng_wm:)
    real (kind=dp_t), intent(in   ) ::   w0macz(lo(1)-ng_wm: ,lo(2)-ng_wm: ,lo(3)-ng_wm:)
    real (kind=dp_t), intent(in   ) ::   normal(lo(1)-ng_n:  ,lo(2)-ng_n:  ,lo(3)-ng_n:,:)
    real (kind=dp_t), intent(in   ) :: time, dt, dx(:)
    real (kind=dp_t), intent(inout) :: vr_x, vr_y, vr_z, vr_max
    real (kind=dp_t), intent(inout) :: T_max, enuc_max, kin_ener, U_max
    real (kind=dp_t), intent(inout) :: rhovr_x, rhovr_y, rhovr_z, mass, nzones
    logical,          intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)

    !     Local variables
    integer            :: i, j, k
    real (kind=dp_t)   :: velr, vel, weight
    logical            :: cell_valid

    weight = 1.d0 / 8.d0**(n-1)

    if (.not. spherical == 1) then
       call bl_error("ERROR: geometry not spherical in diag")
    endif

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1) 

             cell_valid = .true.
             if (present(mask)) then
                if ( (.not. mask(i,j,k)) ) cell_valid = .false.
             end if

             if (cell_valid .and. s(i,j,k,rho_comp) > base_cutoff_density) then
                   
                velr = u(i,j,k,1)*normal(i,j,k,1) + &
                       u(i,j,k,2)*normal(i,j,k,2) + &
                       u(i,j,k,3)*normal(i,j,k,3) + w0r(i,j,k)

                vel = sqrt( (u(i,j,k,1)+HALF*(w0macx(i,j,k)+w0macx(i+1,j,k)))**2 + &
                            (u(i,j,k,2)+HALF*(w0macy(i,j,k)+w0macy(i,j+1,k)))**2 + &
                            (u(i,j,k,3)+HALF*(w0macz(i,j,k)+w0macz(i,j,k+1)))**2)
                
                vr_max = max(vr_max,abs(velr))
                
                vr_x = vr_x + weight*velr*normal(i,j,k,1)
                vr_y = vr_y + weight*velr*normal(i,j,k,2)
                vr_z = vr_z + weight*velr*normal(i,j,k,3)
                
                rhovr_x = rhovr_x + weight*s(i,j,k,rho_comp)*velr*normal(i,j,k,1)
                rhovr_y = rhovr_y + weight*s(i,j,k,rho_comp)*velr*normal(i,j,k,2)
                rhovr_z = rhovr_z + weight*s(i,j,k,rho_comp)*velr*normal(i,j,k,3)
                
                mass = mass + weight*s(i,j,k,rho_comp)
                nzones = nzones + weight
                
                T_max = max(T_max,s(i,j,k,temp_comp))
                enuc_max = max(enuc_max,rho_Hnuc(i,j,k)/s(i,j,k,rho_comp))

                kin_ener = kin_ener + weight*s(i,j,k,rho_comp)*vel**2

                U_max = max(U_max,vel)

             endif

          enddo
       enddo
    enddo

  end subroutine diag_3d

end module diag_module
