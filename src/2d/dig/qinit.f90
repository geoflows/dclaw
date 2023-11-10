
subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use geoclaw_module, only: sea_level
    use amr_module, only: t0
    use qinit_module, only: qinit_type,add_perturbation
    use qinit_module, only: variable_eta_init
    use qinit_module, only: force_dry,use_force_dry,mx_fdry, my_fdry
    use qinit_module, only: xlow_fdry, ylow_fdry, xhi_fdry, yhi_fdry
    use qinit_module, only: dx_fdry, dy_fdry
    use qinit_module, only: tend_force_dry
    use qinit_module, only: qinitwork
    use digclaw_module  !!DIG specify which variables

    
    !implicit none  !!DIG need to declare variables properly
    implicit double precision (a-h,o-z)  !!DIG temporary
    
    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Locals
    integer :: i,j,m, ii,jj
    real(kind=8) :: x,y
    real(kind=8) :: veta(1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8) :: ddxy
    
    
    if (variable_eta_init) then
        ! Set initial surface eta based on eta_init
        call set_eta_init(mbc,mx,my,xlower,ylower,dx,dy,t0,veta)
      else
        veta = sea_level  ! same value everywhere
      endif

    q(2:3,:,:) = 0.d0   ! set momenta to zero

    forall(i=1:mx, j=1:my)
        q(1,i,j) = max(0.d0, veta(i,j) - aux(1,i,j))
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
                      q(1,i,j) = 0.d0
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
    
    !=============
    ! DIG modifications:
    

      xhigher = xlower + (mx-0.5d0)*dx
      yhigher = ylower + (my-0.5d0)*dy

      do mf =1,mqinitfiles

        if ((xlower.le.xhiqinit(mf).and.xhigher.ge.xlowqinit(mf)).and. &
           (ylower.le.yhiqinit(mf).and.yhigher.ge.ylowqinit(mf))) then

            xintlow = dmax1(xlower,xlowqinit(mf))
            xinthi  = dmin1(xhigher,xhiqinit(mf))
            istart  = max(1-mbc,int(0.5 + (xintlow-xlower)/dx)-mbc)
            iend    = min(mx+mbc,int(1.0 + (xinthi-xlower)/dx)+mbc)

            yintlow = dmax1(ylower,ylowqinit(mf))
            yinthi  = dmin1(yhigher,yhiqinit(mf))
            jstart  = max(1-mbc,int(0.5 + (yintlow-ylower)/dy)-mbc)
            jend    = min(my+mbc,int(1.0 + (yinthi-ylower)/dy)+mbc)

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
                     dq = topointegral(ximc,xc,xipc,yjmc,yc,yjpc, &
                        xlowqinit(mf),ylowqinit(mf), &
                        dxqinit(mf),dyqinit(mf), &
                        mxqinit(mf),myqinit(mf), &
                        qinitwork(i1:i2), 1)
                         ! qinitwork(i0qinit(mf):i0qinit(mf)+mqinit(mf)-1) ,1)
                     dq=dq/((xipc-ximc)*(yjpc-yjmc)*aux(2,i,j))

                     if (iqinit(mf).le.meqn) then
                        q(iqinit(mf),i,j) = q(iqinit(mf),i,j) + dq
                     else
                        if (dq-aux(1,i,j).gt.0.d0) then
                          q(1,i,j) = dmax1(q(1,i,j),dq-aux(1,i,j))
                        endif
                     endif

                  endif

               enddo
            enddo
         endif
      enddo

      initm = 0
      initchi = 0
      do mf =1,mqinitfiles
         if (iqinit(mf).eq.2) initu=1
         if (iqinit(mf).eq.3) initv = 1
         if (iqinit(mf).eq.4) initm = 1
         if (iqinit(mf).eq.5) initpv = 1
         if (iqinit(mf).eq.6) initchi = 1
      enddo

      do i=1-mbc,mx+mbc
         do j=1-mbc,my+mbc
               if (initm.eq.0) then
                  q(4,i,j) = m0*q(1,i,j)
               else
                  q(4,i,j) = q(1,i,j)*q(4,i,j)
               endif
               if (initchi.eq.0) then
                  q(6,i,j) = 0.5*q(1,i,j)
               else
                  q(6,i,j) = q(1,i,j)*q(6,i,j)
               endif
               if (initu.eq.1) then
                  q(2,i,j) = q1(1,i,j)*q(2,i,j)
               endif
               if (initv.eq.1) then
                  q(3,i,j) = q(1,i,j)*q(3,i,j)
               endif
               if (initpv.eq.1) then
                  q(5,i,j) = q(1,i,j)*q(5,i,j)
               endif
               if (q(1,i,j).le.drytolerance) then
                  do m = 1,meqn
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
            p_initialized = 1
         case (0)
            !set to hydrostatic
            do i=1-mbc,mx+mbc
               do j=1-mbc,my+mbc
                 if (bed_normal.eq.1) gmod = grav*cos(aux(i_theta,i,j))
                 q(5,i,j) = rho_f*gmod*q(1,i,j)
               enddo
            enddo
            p_initialized = 1
         case(1:2)
            !set to failure
            do i=1-mbc,mx+mbc
               do j=1-mbc,my+mbc
                  p_ratioij = init_pmin_ratio
                  if (q(1,i,j).le.drytolerance) then
                     q(5,i,j) = init_pmin_ratio*rho_f*gmod*q(1,i,j)
                     cycle
                  endif
                  call admissibleq(q(1,i,j),q(2,i,j),q(3,i,j), &
                             q(4,i,j),q(5,i,j),u,v,sv,aux(i_theta,i,j))
                  if (bed_normal.eq.1) then
                     gmod = grav*cos(aux(i_theta,i,j))
                     p_ratioij = init_pmin_ratio &
                         + (init_pmin_ratio - 1.0)*aux(1,i,j)/q(1,i,j)
                  endif

                  rho_dig = sv*rho_s + (1.d0-sv)*rho_f !!DIG rho array in geoclaw_module conflicts 
                  q(5,i,j) = p_ratioij*rho_dig*gmod*q(1,i,j)
               enddo
            enddo
            p_initialized = 1
         case(3:4)
            !p will be updated in b4step2
            do i=1-mbc,mx+mbc
               do j=1-mbc,my+mbc
                 p_ratioij = init_pmin_ratio
                 if (q(1,i,j).le.drytolerance) then
                     q(5,i,j) = init_pmin_ratio*rho_f*gmod*q(1,i,j)
                     cycle
                 endif
                 call admissibleq(q(1,i,j),q(2,i,j),q(3,i,j), &
                             q(4,i,j),q(5,i,j),u,v,sv,aux(i_theta,i,j))
                 if (bed_normal.eq.1) then
                     gmod = grav*cos(aux(i,j,i_theta))
                     p_ratioij = init_pmin_ratio &
                         + (init_pmin_ratio - 1.0)*aux(1,i,j)/q(1,i,j)
                 endif
                 rho = sv*rho_s + (1.0-sv)*rho_f
                     pfail = p_ratioij*rho_dig*gmod*q(1,i,j)
                     q(5,i,j) = pfail - abs(pfail)
               enddo
            enddo
            p_initialized = 0
      end select
      !write(*,*) 'qinit:init,dx,dy:',init_pmin_ratio,dx,dy

      !===============set factor of safety=============================
      call calc_taudir(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

end subroutine qinit
