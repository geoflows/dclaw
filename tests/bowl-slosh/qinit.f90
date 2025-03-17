
subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use geoclaw_module, only: sea_level
    use amr_module, only: t0
    !use qinit_module, only: qinit_type,add_perturbation
    !use qinit_module, only: variable_eta_init
    !use qinit_module, only: force_dry,use_force_dry,mx_fdry, my_fdry
    !use qinit_module, only: xlow_fdry, ylow_fdry, xhi_fdry, yhi_fdry
    !use qinit_module, only: dx_fdry, dy_fdry
    !use qinit_module, only: tend_force_dry
    !use qinit_module, only: qinitwork
    use qinit_module
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
    

    ! Parameters for problem
    real(kind=8), parameter :: a = 1.d0
    real(kind=8), parameter :: sigma = 0.5d0
    real(kind=8), parameter :: h0 = 0.1d0

    ! Other storage
    real(kind=8) :: omega,eta
    
    omega = sqrt(2.d0 * grav * h0) / a
    
    q(:,:,:) = 0.d0  

    do i=1-mbc,mx+mbc
        x = xlower + (i - 0.5d0)*dx
        do j=1-mbc,my+mbc
            y = ylower + (j - 0.5d0) * dy
            eta = sigma * h0 / a**2 * (2.d0 * x - sigma)
            
            q(1,i,j) = max(0.d0,eta - aux(1,i,j))
            q(2,i,j) = 0.d0
            q(3,i,j) = sigma * omega * q(1,i,j)
        enddo
    enddo
    
    
    !=============
    ! DIG modifications:
    

      xhigher = xlower + (mx-0.5d0)*dx
      yhigher = ylower + (my-0.5d0)*dy

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
      do mf =1,mqinitfiles
         if (iqinit(mf).eq.i_hu) initu=1
         if (iqinit(mf).eq.i_hv) initv = 1
         if (iqinit(mf).eq.i_hm) initm = 1
         if (iqinit(mf).eq.i_pb) initpv = 1
         if (iqinit(mf).eq.i_hchi) initchi = 1
      enddo

      do i=1-mbc,mx+mbc
         do j=1-mbc,my+mbc
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
                 if (bed_normal.eq.1) then
                     gz = grav*dcos(aux(i_theta,i,j))
                 else
                     gz=grav
                 endif
                 q(i_pb,i,j) = rho_f*gz*q(i_h,i,j)
               enddo
            enddo
            p_initialized = 1
         case(1:2)
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
                          u,v,sv,chi,rho,gz)

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

      end select
      !write(*,*) 'qinit:init,dx,dy:',init_pmin_ratio,dx,dy

      !===============set factor of safety=============================
      call calc_taudir(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

end subroutine qinit
