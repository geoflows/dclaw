!! D-Claw specific core file
!! This file is a modified version of
!! clawpack/geoclaw/src/2d/shallow/qinit.f90 
!!
subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

    use geoclaw_module, only: sea_level,grav,coordinate_system,dry_tolerance
    use amr_module, only: t0
    use qinit_module, only: qinit_type,add_perturbation
    use qinit_module, only: variable_eta_init
    use qinit_module, only: force_dry,use_force_dry,mx_fdry, my_fdry
    use qinit_module, only: xlow_fdry, ylow_fdry
    use qinit_module, only: dxqinit,dyqinit,i0qinit,iqinit,mqinit
    use qinit_module, only: mxqinit,myqinit,xhiqinit,xlowqinit
    use qinit_module, only: yhiqinit,ylowqinit
    use qinit_module, only: dx_fdry, dy_fdry
    use qinit_module, only: tend_force_dry
    use qinit_module, only: qinitwork,mqinitfiles


    use digclaw_module, only: qfix,calc_pmin,calc_taudir
    use digclaw_module, only: bed_normal,chi0,init_ptype,segregation
    use digclaw_module, only: i_theta,m0,rho_f,rho_s,init_pmin_ratio
    use digclaw_module, only: i_h,i_hu,i_hv,i_hm,i_pb,i_hchi,i_hf,i_hs,i_ent

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Locals
    integer :: i,j,m,ii,jj,i1,i2,iend,istart,jstart,jend,mf
    integer :: initchi,initm,initpv,initu,initv,inithf
    real(kind=8) :: x,y,xc,xhigher,xim,ximc,xinthi,xintlow,xip,xipc
    real(kind=8) :: yc,yhigher,yinthi,yintlow,yjm,yjmc,yjp,yjpc
    real(kind=8) :: veta(1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8) :: ddxy
    real(kind=8) :: u,v,sv,dq,p_ratioij,rho,gz,chi
    ! Topography integral function
    real(kind=8) :: topointegral

    if (variable_eta_init) then
        ! Set initial surface eta based on eta_init
        call set_eta_init(mbc,mx,my,xlower,ylower,dx,dy,t0,veta)
        ! DIG: this is redundant to current method of passing a file.
    else
        veta = sea_level  ! same value everywhere
    endif

    q(:,:,:) = 0.d0   ! set all to zero (set further below)

    ! 1/12/24: To make consistent with dclaw 4 we needed
    ! to change the typical geoclaw behavior and set sea
    ! level in ghost cells (rather than just in interior cells).
    ! This is necessary because eta in ghost cells is used to
    ! calculate gradients in eta used by (at least) taudir_x and
    ! taudir_y.
    ! 1/30/24: Confirmed with DLG that this is the behavior we want.

    forall(i=1-mbc:mx+mbc, j=1-mbc:my+mbc)
        q(i_h,i,j) = max(0.d0, veta(i,j) - aux(1,i,j))
    end forall

    if (use_force_dry .and. (t0 <= tend_force_dry)) then
     ! only use the force_dry if it specified on a grid that matches the
     ! resolution of this patch, since we only check the cell center:
     ddxy = max(abs(dx-dx_fdry), abs(dy-dy_fdry))
     if (ddxy < 0.01d0*min(dx_fdry,dy_fdry)) then
       do i=1,mx
          x = xlower + (i-0.5d0)*dx
          ii = int((x - xlow_fdry + 1d-7) / dx_fdry)
          do j=1,my
              y = ylower + (j-0.5d0)*dy
              jj = int((y - ylow_fdry + 1d-7) / dy_fdry)
              jj = my_fdry - jj  ! since index 1 corresponds to north edge
              if ((ii>=1) .and. (ii<=mx_fdry) .and. &
                  (jj>=1) .and. (jj<=my_fdry)) then
                  ! grid cell lies in region covered by force_dry,
                  ! check if this cell is forced to be dry
                  ! Otherwise don't change value set above:
                  if (force_dry(ii,jj) == 1) then
                      q(i_h,i,j) = 0.d0
                      endif
                  endif
          enddo ! loop on j
       enddo ! loop on i
       endif ! dx and dy agree with dx_fdry, dy_fdry
    endif ! use_force_dry


    ! Add perturbation to initial conditions
    if (qinit_type > 0) then
        call add_perturbation(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    endif

    if (.false.) then
        open(23, file='fort.aux',status='unknown',form='formatted')
        print *,'Writing out aux arrays'
        print *,' '
        do j=1,my
            do i=1,mx
                write(23,*) i,j,(q(m,i,j),m=1,meqn)
            enddo
        enddo
        close(23)
    endif

    !===========================================
    ! D-Claw specific initialization of q
    !===========================================

    ! set i_hs to i_ent
    ! and set i_hf to 0.0d0
    do i=1,mx
      do j=1,my
         q(i_hs,i,j) = aux(i_ent,i,j)
         q(i_hf,i,j) = 0.d0
      enddo
    enddo

      xhigher = xlower + (mx-0.5d0)*dx ! DIG: higher is cell center but lower is cell edge. Correct?
      yhigher = ylower + (my-0.5d0)*dy ! below, we are comparing against xlow/hi, ylow/hi which are cell edges.

      ! loop through qinitfiles and fill in q. topo has already been set.
      do mf =1,mqinitfiles

        if ((xlower.le.xhiqinit(mf).and.xhigher.ge.xlowqinit(mf)).and. &
           (ylower.le.yhiqinit(mf).and.yhigher.ge.ylowqinit(mf))) then

            xintlow = dmax1(xlower,xlowqinit(mf))
            xinthi  = dmin1(xhigher,xhiqinit(mf))
            istart  = max(1-mbc,int(0.5d0 + (xintlow-xlower)/dx)-mbc)
            iend    = min(mx+mbc,int(1.0d0 + (xinthi-xlower)/dx)+mbc)

            yintlow = dmax1(ylower,ylowqinit(mf))
            yinthi  = dmin1(yhigher,yhiqinit(mf))
            jstart  = max(1-mbc,int(0.5d0 + (yintlow-ylower)/dy)-mbc)
            jend    = min(my+mbc,int(1.0d0 + (yinthi-ylower)/dy)+mbc)

            do i=istart,iend
               x = xlower + (i-0.5d0)*dx
               xim = x - 0.5d0*dx
               xip = x + 0.5d0*dx
               do j=jstart,jend
                  y = ylower + (j-0.5d0)*dy
                  yjm = y - 0.5d0*dy
                  yjp = y + 0.5d0*dy

                  if (xip.gt.xlowqinit(mf).and.xim.lt.xhiqinit(mf) &
                     .and.yjp.gt.ylowqinit(mf) &
                     .and.yjm.lt.yhiqinit(mf)) then

                     xipc=dmin1(xip,xhiqinit(mf))
                     ximc=dmax1(xim,xlowqinit(mf))
                     xc=0.5d0*(xipc+ximc)

                     yjpc=dmin1(yjp,yhiqinit(mf))
                     yjmc=dmax1(yjm,ylowqinit(mf))
                     yc=0.5d0*(yjmc+yjpc)

                     i1 = i0qinit(mf)
                     i2 = i0qinit(mf)+mqinit(mf)-1
                     dq = topointegral(ximc,xipc,yjmc,yjpc, &
                        xlowqinit(mf),ylowqinit(mf), &
                        dxqinit(mf),dyqinit(mf), &
                        mxqinit(mf),myqinit(mf), &
                        qinitwork(i1:i2), 1)
                         ! qinitwork(i0qinit(mf):i0qinit(mf)+mqinit(mf)-1) ,1)

                     if (coordinate_system == 2) then
                        dq=dq/((xipc-ximc)*(yjpc-yjmc)*aux(2,i,j))
                     else
                        dq=dq/((xipc-ximc)*(yjpc-yjmc))
                     endif


                     if (iqinit(mf).le.meqn) then
                        q(iqinit(mf),i,j) = q(iqinit(mf),i,j) + dq
                     else
                        ! if a file is provide for eta, set h.
                        ! if files are provided for both h and eta (
                        ! not recommended, the second will take precedence.
                        ! if eta is provided as less than topo, then
                        ! h will be zero.
                        if (dq-aux(1,i,j).gt.0.d0) then
                          q(i_h,i,j) = dmax1(q(i_h,i,j),dq-aux(1,i,j))
                        endif
                     endif

                  endif

               enddo
            enddo
         endif
      enddo

      ! adjust q values.
      initu = 0
      initv = 0
      initm = 0
      initchi = 0
      initpv = 0
      inithf = 0
      ! hs cannot be user-specified at initialization
      do mf =1,mqinitfiles
         if (iqinit(mf).eq.i_hu) initu=1
         if (iqinit(mf).eq.i_hv) initv = 1
         if (iqinit(mf).eq.i_hm) initm = 1
         if (iqinit(mf).eq.i_pb) initpv = 1
         if (iqinit(mf).eq.i_hchi) initchi = 1
         if (iqinit(mf).eq.i_hf) inithf = 1
      enddo


      do j=1-mbc,my+mbc
         do i=1-mbc,mx+mbc
               if (initm.eq.0) then
                  if (dabs((q(i_h,i,j) + aux(1,i,j)) - veta(i,j)).lt.1d-6) then
                    q(i_hm,i,j) = 0.d0 ! If eta is sea level (stored in veta), assume m is zero
                    ! DIG: This may need to change if treatment of sea level does not use veta.
                    ! potentially have a m0fill and m0perm
                  else
                    q(i_hm,i,j) = m0*q(i_h,i,j)
                  endif
               else
                  q(i_hm,i,j) = q(i_h,i,j)*q(i_hm,i,j)
               endif
               if (segregation.eq.1) then
                  ! only initialize chi if segregation is enabled
                  if (initchi.eq.0) then
                    q(i_hchi,i,j) = chi0*q(i_h,i,j)
                  else
                    q(i_hchi,i,j) = q(i_h,i,j)*q(i_hchi,i,j)
                  endif
               endif
               if (initu.eq.1) then
                  q(i_hu,i,j) = q(i_h,i,j)*q(i_hu,i,j)
               endif
               if (initv.eq.1) then
                  q(i_hv,i,j) = q(i_h,i,j)*q(i_hv,i,j)
               endif
               if (initpv.eq.1) then
                  q(i_pb,i,j) = q(i_h,i,j)*q(i_pb,i,j)
               endif


               if (q(i_h,i,j).le.dry_tolerance) then
                  do m = 1,meqn
                     if ((m.eq.i_hf).or.(m.eq.i_hs)) continue
                     q(m,i,j) = 0.d0
                  enddo
               endif
         enddo
      enddo

      ! ========= Pressure initialization for Mobilization Modeling======
      call calc_pmin(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

      select case (init_ptype)
         case (-1)
            !p should already be 0 or set by qinit file
         case (0)
            !set to hydrostatic
            do i=1-mbc,mx+mbc
               do j=1-mbc,my+mbc
                 if (bed_normal.eq.1) then
                     gz = grav*dcos(aux(i_theta,i,j))
                 else
                     gz=grav
                 endif
                 q(i_pb,i,j) = rho_f*gz*q(i_h,i,j)
               enddo
            enddo

         case(1:2) ! DIG: Not yet tested
            !set to failure
            do i=1-mbc,mx+mbc
               do j=1-mbc,my+mbc
                  p_ratioij = init_pmin_ratio

                  if (bed_normal.eq.1) then
                     gz = grav*dcos(aux(i_theta,i,j))
                  else
                     gz = grav
                  endif

                  if (q(i_h,i,j).le.dry_tolerance) then
                     q(i_pb,i,j) = init_pmin_ratio*rho_f*gz*q(i_h,i,j)
                     cycle
                  endif


                  call qfix(q(i_h,i,j),q(i_hu,i,j),q(i_hv,i,j), &
                       q(i_hm,i,j),q(i_pb,i,j),q(i_hchi,i,j), &
                          q(i_hf,i,j),u,v,sv,chi,rho,gz)

                  if (bed_normal.eq.1) then
                     p_ratioij = init_pmin_ratio &
                         + (init_pmin_ratio - 1.0)*aux(1,i,j)/q(i_h,i,j)
                  endif

                  q(i_pb,i,j) = p_ratioij*rho*gz*q(i_h,i,j)

               enddo
            enddo

      end select

      !===============set factor of safety=============================
      call calc_taudir(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

end subroutine qinit
