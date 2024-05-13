! This code is the main program of ``stationary droplet" accompanying the following book.
!
!   An Introduction to the Lattice Boltzmann Method: A Numerical Method for Complex Boundary and Moving Boundary Flows
!   T. Inamuro, M. Yoshino, K. Suzuki
!
!  Copyright @ 2021 Kosuke Suzuki. All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program laplace
 implicit none

 !====================parameters for LBM=============================
 real(8), parameter:: ds = 1.d0		!lattice spacing in lattice unit

 integer, parameter:: imax = 100 	!number of lattice points in x direction
 integer, parameter:: jmax = 100 	!number of lattice points in y direction
 integer, parameter:: kmax = 100 	!number of lattice points in z direction

 real(8), parameter:: D = 40.d0*ds	!initial droplet diameter

 real(8), parameter:: h0 = 5.d0*ds		!characteristic length

 !initial droplet center
 real(8), parameter:: xc = 0.5d0*dble(imax)	
 real(8), parameter:: yc = 0.5d0*dble(jmax)
 real(8), parameter:: zc = D/2.d0+9.d0
				
 !particle velocities (D3Q15) in integer
 integer, parameter:: ci(15) = (/ 0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
 integer, parameter:: cj(15) = (/ 0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/) 
 integer, parameter:: ck(15) = (/ 0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/) 

 !particle velocities (D3Q15) in double precision
 real(8):: cr(1:3, 1:15)	

 !coefficients
 real(8), parameter:: E(15) = (/ 2.d0/9.d0, 1.d0/9.d0, 1.d0/9.d0, 1.d0/9.d0, 1.d0/9.d0, 1.d0/9.d0, 1.d0/9.d0, &
                                1.d0/72.d0, 1.d0/72.d0, 1.d0/72.d0, 1.d0/72.d0, &
                                1.d0/72.d0, 1.d0/72.d0, 1.d0/72.d0, 1.d0/72.d0 /)
 real(8), parameter:: H(15) = (/ 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
                                 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
 real(8), parameter:: F(15) = (/ -7.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, &
                                1.d0/24.d0, 1.d0/24.d0, 1.d0/24.d0, 1.d0/24.d0, &
                                1.d0/24.d0, 1.d0/24.d0, 1.d0/24.d0, 1.d0/24.d0 /)
 !====================================================================
 !====================governing parameter=============================
 real(8), parameter:: ratio_rho	= 800.d0	!density ratio
 real(8), parameter:: ratio_mu = 100.d0		!viscosity ratio

 real(8), parameter:: rhog = 1.d0		!gas density
 real(8), parameter:: rhol = ratio_rho*rhog	!liquid density

 real(8), parameter:: mug = 1.2d-2*ds		!gas viscosity
 real(8), parameter:: mul = ratio_mu*mug	!liquid viscosity

 real(8), parameter:: sigma  = 4.2d-1*ds	!interfacial tension
 real(8), parameter:: U_initial  = 0.0075	!initial speed
 real(8), parameter:: g_acc  = 1.84d-6*ds	!grivity acc
 !====================================================================
 !=====================parameters for index function phi==============
 real(8), parameter:: a = 1.d0
 real(8), parameter:: b = 1.d0
 real(8), parameter:: T = 2.93d-1
 real(8), parameter:: phi1 = 2.638d-1
 real(8), parameter:: phi2 = 4.031d-1

 real(8), parameter:: kappaf = 0.06d0*ds**2
 real(8), parameter:: C = 0.d0
 real(8), parameter:: M = (0.5d0-C/3.d0)*ds		!mobility
 !====================================================================
 !======================variables fo index function phi===============
 real(8), dimension(1:15,0:imax,0:jmax,0:kmax):: feq		!equilibrium distribution functions

 real(8), dimension(0:imax,0:jmax,0:kmax):: phi		!index function

 !intermediate variables for index function
 real(8), dimension(0:imax,0:jmax,0:kmax):: lap_phi
 real(8), dimension(1:3, 0:imax,0:jmax,0:kmax):: grad_phi
 real(8), dimension(0:imax,0:jmax,0:kmax):: p0
 real(8), dimension(1:3, 1:3, 0:imax,0:jmax,0:kmax):: gphi
 real(8), dimension(1:3, 1:3, 0:imax,0:jmax,0:kmax):: pcap
 real(8), dimension(1:3, 0:imax,0:jmax,0:kmax):: div_pcap

 real(8):: phi_min, phi_max			!maximum and minimum values of phi
 !====================================================================
 !=======parameters for pressure and velocity calculations============
 integer, parameter:: itr_max = 2			!number of iteration
 real(8), parameter:: omega_max = rhol/dble(itr_max)	!maximum value of acceleration parameter
 real(8), parameter:: lambda = 1.d0			!coefficient of stabilization term
 !====================================================================
 !========variables for pressure and velocity calculations=============
 real(8), dimension(1:15,0:imax,0:jmax,0:kmax):: geq		!equilibrium distribution functions

 real(8), dimension(0:imax,0:jmax,0:kmax):: rho			!density
 real(8), dimension(0:imax,0:jmax,0:kmax):: p			!pressure
 real(8), dimension(0:imax,0:jmax,0:kmax):: u, v, w		!velocity

 real(8), dimension(0:imax,0:jmax,0:kmax):: omega		!acceleration parameter
 real(8), dimension(1:15,0:imax,0:jmax,0:kmax):: delta_p	!intermediate variable for pressure calculation

 real(8), dimension(0:imax,0:jmax,0:kmax):: rhonext		!density at next time step
 real(8), dimension(0:imax,0:jmax,0:kmax):: unext, vnext, wnext	!velocity at next time step

 real(8), dimension(0:imax,0:jmax,0:kmax):: mu			!viscosity
 real(8), dimension(0:imax,0:jmax,0:kmax)::Au			!parameter determining viscosity

 real(8), dimension(1:3, 0:imax,0:jmax,0:kmax)::grad_mu		!gradient of viscosity
 real(8), dimension(1:3,1:3, 0:imax,0:jmax,0:kmax)::grad_u	!gradient of velocity
 real(8), dimension(1:3,0:imax,0:jmax,0:kmax):: vcap		!V_alpha

 real(8), dimension(0:imax,0:jmax,0:kmax)::lap_u		!laplacian of u
 real(8), dimension(0:imax,0:jmax,0:kmax)::lap_v		!laplacian of v
 real(8), dimension(0:imax,0:jmax,0:kmax)::lap_w		!laplacian of w
 real(8), dimension(0:imax,0:jmax,0:kmax)::laplap_u		!laplacian of lap_u
 real(8), dimension(0:imax,0:jmax,0:kmax)::laplap_v		!laplacian of lap_v
 real(8), dimension(0:imax,0:jmax,0:kmax)::laplap_w		!laplacian of lap_w
 
 real(8), dimension(1:3, 0:imax,0:jmax,0:kmax)::grad_rho	!gradient of density
 real(8), dimension(1:3, 0:imax,0:jmax,0:kmax)::normal		!normal direction
 real(8), dimension(0:imax,0:jmax,0:kmax)::chi			!curvature 
 real(8), dimension(1:3, 0:imax,0:jmax,0:kmax)::fsv		!interfacial tension (CSF)
 !====================================================================
 !=====================other parameters and variables=================
 integer, parameter:: step = 20000		!maximum time step

 integer, parameter:: start_sigma = 5000	!time step starting to apply interfacial tension
 integer, parameter:: end_sigma = 15000		!time step when interfacial tension reaches its desired value
 real(8):: sigma_temp				!value of interfacial tension at each time step

 real(8):: epsilon = 1.d-12			!threshold value for avoiding division by zero

 real(8):: krone(1:3,1:3)			!Kronecker's delta
 
 real(8):: gtemp				!temporary variable

 real(8), parameter:: pi = dacos(-1.d0)		!circle ratio

 real(8), parameter:: rinput = 0.5d0*D			!input radius
 real(8), parameter:: vinput = 4.d0/3.d0*pi*rinput**3	!input volume

 real(8):: veff, reff				!effective droplet volume, radius
 real(8):: dp_th 				!theoretical value of Laplace pressure
 real(8):: pg, pl, dp_calc			!calculated values of gas and liquid pressures and their difference
 real(8):: err					!error from Laplace pressure

 integer:: i, j, k, l, n, alpha, beta, itr	!index
 character*6:: num				!characters for file number
 !====================================================================
 !-----------substituting--------------------------
 do l=1,15
  cr(1,l) = dble(ci(l))
  cr(2,l) = dble(cj(l))
  cr(3,l) = dble(ck(l))
 end do

 krone(1,1) = 1.d0; krone(1,2) = 0.d0; krone(1,3) = 0.d0
 krone(2,1) = 0.d0; krone(2,2) = 1.d0; krone(2,3) = 0.d0
 krone(3,1) = 0.d0; krone(3,2) = 0.d0; krone(3,3) = 1.d0
 !------------------------------------------------
 !--output of computational parameters------------
  open(8,file='parameter.txt')
  write(8,'(11a20)') 'h0','ratio_rho','ratio_mu','mug','mul','sigma','kappaf','C','M','itr_max','omega_max'
  write(8,'(11f20.10)') h0, ratio_rho, ratio_mu, mug, mul, sigma, kappaf, C, M, dble(itr_max), omega_max
  close(8)
 !------------------------------------------------
 !--------------file open-------------------------
  open(8,file='dp.txt')
  open(9,file='reff.txt')
 !------------------------------------------------


 !======================initial condition============================
 do k=0, kmax
  do j=0,jmax
   do i=0,imax

    u(i,j,k) = 0.d0
    v(i,j,k) = 0.d0
    w(i,j,k) = 0.d0
    p(i,j,k) = 1.d0/3.d0
    phi(i,j,k) = phi1

!    if(k <= h0) then
!      phi(i,j,k) = phi2
    if ((i*ds-xc)**2+(j*ds-yc)**2+(k*ds-zc)**2 <= (0.5d0*D)**2) then
      phi(i,j,k) = phi2 
      ! w(i,j,k)=-U_initial    
    end if

    !please select linear or trigonometric interpolation for function form of density (c.f. line 480)
    !rho(i,j,k) = (phi(i,j,k)-phi1)/(phi2-phi1)*(rhol-rhog) + rhog		!linear
    rho(i,j,k) = 0.5d0*(rhol-rhog)*(dsin(pi*(phi(i,j,k)-0.5d0*(phi2+phi1))/(phi2-phi1)) + 1.d0) + rhog	!trigonometric

   end do
  end do
 end do
 !====================================================================


 !======================time evolution================================
 DO n=1,step
 !----------------------------------------------------------phi calculation
 !-----------lap_phi (all sides are periodic)-------
 do k=0,kmax
  do j=0,jmax
   do i=0,imax
    lap_phi(i,j,k) = - 14.d0*phi(i,j,k)

    !bulk
    if((i/=0).and.(i/=imax).and.(j/=0).and.(j/=jmax).and.(k/=0).and.(k/=kmax)) then
     do l=2,15
      lap_phi(i,j,k) = lap_phi(i,j,k) + phi(i+ci(l),j+cj(l),k+ck(l))
     end do

    !boundary (periodic)
    else
     do l=2,15
      lap_phi(i,j,k) = lap_phi(i,j,k) + phi(perx(i+ci(l)),pery(j+cj(l)),perz(k+ck(l)))
     end do
    end if

    lap_phi(i,j,k) = lap_phi(i,j,k) / (5.d0*ds**2)

   end do
  end do
 end do
 !-----------------------------------------------
 !-----------grad_phi (all sides are periodic)---
 do k=0,kmax
  do j=0,jmax
   do i=0,imax

    do alpha=1,3
     grad_phi(alpha,i,j,k) = 0.d0

     !bulk
     if((i/=0).and.(i/=imax).and.(j/=0).and.(j/=jmax).and.(k/=0).and.(k/=kmax)) then
      do l=2,15
       grad_phi(alpha,i,j,k) = grad_phi(alpha,i,j,k) + cr(alpha,l)*phi(i+ci(l),j+cj(l),k+ck(l))
      end do

     !boundary (periodic)
     else
      do l=2,15
       grad_phi(alpha,i,j,k) = grad_phi(alpha,i,j,k) + cr(alpha,l)*phi(perx(i+ci(l)),pery(j+cj(l)),perz(k+ck(l)))
      end do
     end if

     grad_phi(alpha,i,j,k) = grad_phi(alpha,i,j,k) / (10.d0*ds)

    end do

   end do
  end do
 end do
 !-----------------------------------------------
 !-----------p0----------------------------------
 do k=0,kmax
  do j=0,jmax
   do i=0,imax
    p0(i,j,k) = phi(i,j,k)*T/(1.d0-b*phi(i,j,k)) - a*phi(i,j,k)**2
   end do
  end do
 end do
 !-----------------------------------------------
 !-----------gphi--------------------------------
 do k=0,kmax
  do j=0,jmax
   do i=0,imax
    
    do beta=1,3
     do alpha=1,3
      gphi(alpha,beta,i,j,k) = 4.5d0*grad_phi(alpha,i,j,k)*grad_phi(beta,i,j,k) &
                      - 1.5d0*(grad_phi(1,i,j,k)**2 + grad_phi(2,i,j,k)**2 + grad_phi(3,i,j,k)**2)&
                             *krone(alpha,beta)
     end do
    end do

   end do
  end do
 end do
 !-----------------------------------------------
 !-----------pcap-------------------------------
 do k=0,kmax
  do j=0,jmax
   do i=0,imax
    
    do beta=1,3
     do alpha=1,3
      pcap(alpha,beta,i,j,k) = (p0(i,j,k)-kappaf*phi(i,j,k)*lap_phi(i,j,k) &
                                -0.5d0*kappaf*(grad_phi(1,i,j,k)**2 &
                                             + grad_phi(2,i,j,k)**2 &
                                             + grad_phi(3,i,j,k)**2))*krone(alpha,beta) &
                                + kappaf*grad_phi(alpha,i,j,k)*grad_phi(beta,i,j,k)
     end do
    end do

   end do
  end do
 end do
 !-----------------------------------------------
 !------div_pcap (all sides are periodic)--------
 do k=0,kmax
  do j=0,jmax
   do i=0,imax
    
    do alpha=1,3
     div_pcap(alpha,i,j,k) = 0.d0

     !bulk
     if((i/=0).and.(i/=imax).and.(j/=0).and.(j/=jmax).and.(k/=0).and.(k/=kmax)) then
      do l=2,15
       do beta = 1,3
        div_pcap(alpha,i,j,k) = div_pcap(alpha,i,j,k) + cr(beta,l)*pcap(alpha,beta,i+ci(l),j+cj(l),k+ck(l))
       end do
      end do

     !boundary (periodic)
     else
      do l=2,15
       do beta = 1,3
        div_pcap(alpha,i,j,k) = div_pcap(alpha,i,j,k) &
                              + cr(beta,l)*pcap(alpha,beta,perx(i+ci(l)),pery(j+cj(l)),perz(k+ck(l)))
       end do
      end do
     end if
 
     div_pcap(alpha,i,j,k) = div_pcap(alpha,i,j,k) / (10.d0*ds)

    end do

   end do
  end do
 end do
 !-----------------------------------------------
 !-----------feq---------------------------------
 do k=0,kmax
  do j=0,jmax
   do i=0,imax

    do l=1,15

     gtemp = 0.d0

     do beta=1,3
      do alpha=1,3
       gtemp = gtemp + gphi(alpha,beta,i,j,k)*cr(alpha,l)*cr(beta,l)
      end do
     end do

     feq(l,i,j,k) = H(l)*phi(i,j,k) &
                  + F(l)*(p0(i,j,k)-kappaf*phi(i,j,k)*lap_phi(i,j,k)&
                          -kappaf/6.d0*(grad_phi(1,i,j,k)**2 + grad_phi(2,i,j,k)**2 + grad_phi(3,i,j,k)**2)) &
                  + 3.d0*E(l)*phi(i,j,k)*(cr(1,l)*u(i,j,k)+cr(2,l)*v(i,j,k)+cr(3,l)*w(i,j,k)) &
                  + E(l)*kappaf*gtemp &
                  + E(l)*C*(div_pcap(1,i,j,k)*cr(1,l)+div_pcap(2,i,j,k)*cr(2,l)+div_pcap(3,i,j,k)*cr(3,l))*ds
    end do

   end do
  end do
 end do
 !-----------------------------------------------
 !***time evolution of phi (all sides are periodic)***
 do k=0,kmax
  do j=0,jmax
   do i=0,imax

    phi(i,j,k) = 0.d0

    !bulk
    if((i/=0).and.(i/=imax).and.(j/=0).and.(j/=jmax).and.(k/=0).and.(k/=kmax)) then
     do l=1,15
      phi(i,j,k) = phi(i,j,k) + feq(l,i-ci(l), j-cj(l), k-ck(l))
     end do

    !boundary (periodic)
    else
     do l=1,15
      phi(i,j,k) = phi(i,j,k) + feq(l,perx(i-ci(l)), pery(j-cj(l)), perz(k-ck(l)))
     end do
    end if

   end do
  end do
 end do
 !************************************************
 !----------------------------------------------------------
 !----------------------------------------------------------pressure calculation
 !--------------omega-----------------------------
 do k=0,kmax
  do j=0,jmax
   do i=0,imax
    omega(i,j,k) = 1.d0 + (rho(i,j,k)-rhog)/(rhol-rhog)*(omega_max-1.d0)
   end do
  end do
 end do
 !------------------------------------------------
 !--------geq-------------------------------------
 do k=0,kmax
  do j=0,jmax
   do i=0,imax
    do l=1,15
     geq(l,i,j,k) = E(l)*(3.d0*(cr(1,l)*u(i,j,k)+cr(2,l)*v(i,j,k)+cr(3,l)*w(i,j,k)) &
                        - 1.5d0*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2) &
                        + 4.5d0*(cr(1,l)*u(i,j,k)+cr(2,l)*v(i,j,k)+cr(3,l)*w(i,j,k))**2)!-3*E(l)*ck(l)*(1-(rhog/rho(i,j,k)))*g_acc*ds
    end do
   end do
  end do
 end do
 !------------------------------------------------
 !********iteration********************************
 do itr = 1, itr_max

  !--------delta_p (all sides are periodic)----
  do k=0,kmax
   do j=0,jmax
    do i=0,imax

     !bulk
     if((i/=0).and.(i/=imax).and.(j/=0).and.(j/=jmax).and.(k/=0).and.(k/=kmax)) then
      do l=1,15
       delta_p(l,i,j,k) = 1.5d0*E(l)*(1.d0/rho(i-ci(l),j-cj(l),k-ck(l)) + 1.d0/rho(i,j,k)) &
                                    * (p(i-ci(l),j-cj(l),k-ck(l)) - p(i,j,k))
      end do

     !boundary (periodic)
     else
      do l=1,15
       delta_p(l,i,j,k) = 1.5d0*E(l)*(1.d0/rho(perx(i-ci(l)),pery(j-cj(l)),perz(k-ck(l))) + 1.d0/rho(i,j,k)) &
                                    *(p(perx(i-ci(l)),pery(j-cj(l)),perz(k-ck(l))) - p(i,j,k))
      end do
     end if

    end do
   end do
  end do
  !-------------------------------------
  !--------p (all sides are periodic)----
  do k=0,kmax
   do j=0,jmax
    do i=0,imax

     !bulk
     if((i/=0).and.(i/=imax).and.(j/=0).and.(j/=jmax).and.(k/=0).and.(k/=kmax)) then
      do l=1,15
       p(i,j,k) = p(i,j,k) + omega(i,j,k)/3.d0*(delta_p(l,i,j,k) + geq(l,i-ci(l), j-cj(l), k-ck(l)))
      end do

     !boundary (periodic)
     else
      do l=1,15
       p(i,j,k) = p(i,j,k) &
                 + omega(i,j,k)/3.d0*(delta_p(l,i,j,k) + geq(l,perx(i-ci(l)), pery(j-cj(l)), perz(k-ck(l))))
      end do
     end if

    end do
   end do
  end do
  !-------------------------------------

 end do
 !************************************************
 !----------------------------------------------------------
 !----------------------------------------------------------velocity calculation
 !---------rho at next time step------------------
 phi_min = phi2
 phi_max = phi1

 do k=0,kmax
  do j=0,jmax
   do i=0,imax
    if(phi_min > phi(i,j,k)) phi_min = phi(i,j,k)
    if(phi_max < phi(i,j,k)) phi_max = phi(i,j,k)
   end do
  end do
 end do


 do k=0,kmax
  do j=0,jmax
   do i=0,imax

    !please select linear or trigonometric interpolation for function form of density (c.f. line 186)
    !rhonext(i,j,k) = (phi(i,j,k)-phi_min)/(phi_max-phi_min) + rhog		!linear
    rhonext(i,j,k) = 0.5d0*(rhol-rhog)*(dsin(pi*(phi(i,j,k)-0.5d0*(phi_max+phi_min))/(phi_max-phi_min)) + 1.d0)&
                   + rhog							!trigonometric

   end do
  end do
 end do
 !-----------------------------------------------
 !--------delta_p (all sides are periodic)-------
  do k=0,kmax
   do j=0,jmax
    do i=0,imax

     !bulk
     if((i/=0).and.(i/=imax).and.(j/=0).and.(j/=jmax).and.(k/=0).and.(k/=kmax)) then
      do l=1,15
       delta_p(l,i,j,k) = 1.5d0*E(l)*(1.d0/rhonext(i-ci(l),j-cj(l),k-ck(l)) + 1.d0/rhonext(i,j,k)) &
                                    *(p(i-ci(l),j-cj(l),k-ck(l)) - p(i,j,k))
      end do

     !boundary (periodic)
     else
      do l=1,15
       delta_p(l,i,j,k) = 1.5d0*E(l)*(1.d0/rhonext(perx(i-ci(l)),pery(j-cj(l)),perz(k-ck(l))) &
                                    + 1.d0/rhonext(i,j,k)) &
                                    *(p(perx(i-ci(l)),pery(j-cj(l)),perz(k-ck(l))) - p(i,j,k))
      end do
     end if

    end do
   end do
  end do
 !------------------------------------------------
 !---------mu (linear interpolation)--------------
 do k=0,kmax
  do j=0,jmax
   do i=0,imax
    mu(i,j,k) = mug + (rho(i,j,k)-rhog)/(rhol-rhog)*(mul-mug)
   end do
  end do
 end do
 !------------------------------------------------
 !---------Au-------------------------------------
 do k=0,kmax
  do j=0,jmax
   do i=0,imax
    Au(i,j,k) = 1.d0 - 6.d0*mu(i,j,k)/rho(i,j,k)/ds 
   end do
  end do
 end do
 !------------------------------------------------
 !---------grad_mu (all sides are periodic)--------
 do k=0,kmax
  do j=0,jmax
   do i=0,imax

    do alpha=1,3
     grad_mu(alpha,i,j,k) = 0.d0

     !bulk
     if((i/=0).and.(i/=imax).and.(j/=0).and.(j/=jmax).and.(k/=0).and.(k/=kmax)) then
      do l=2,15
       grad_mu(alpha,i,j,k) = grad_mu(alpha,i,j,k) + cr(alpha,l)*mu(i+ci(l),j+cj(l),k+ck(l))
      end do

     !boundary (periodic)
     else
      do l=2,15
       grad_mu(alpha,i,j,k) = grad_mu(alpha,i,j,k) + cr(alpha,l)*mu(perx(i+ci(l)),pery(j+cj(l)),perz(k+ck(l)))
      end do
     end if

     grad_mu(alpha,i,j,k) = grad_mu(alpha,i,j,k) / (10.d0*ds)

    end do

   end do
  end do
 end do
 !-----------------------------------------------
 !--------grad_u (all sides are periodic)--------
 do k=0,kmax
  do j=0,jmax
   do i=0,imax

    do beta=1,3
     grad_u(1,beta,i,j,k) = 0.d0
     grad_u(2,beta,i,j,k) = 0.d0
     grad_u(3,beta,i,j,k) = 0.d0

     !bulk
     if((i/=0).and.(i/=imax).and.(j/=0).and.(j/=jmax).and.(k/=0).and.(k/=kmax)) then
      do l=2,15
       grad_u(1,beta,i,j,k) = grad_u(1,beta,i,j,k) + cr(beta,l)*u(i+ci(l),j+cj(l),k+ck(l))
       grad_u(2,beta,i,j,k) = grad_u(2,beta,i,j,k) + cr(beta,l)*v(i+ci(l),j+cj(l),k+ck(l))
       grad_u(3,beta,i,j,k) = grad_u(3,beta,i,j,k) + cr(beta,l)*w(i+ci(l),j+cj(l),k+ck(l))
      end do

     !boundary (periodic)
     else
      do l=2,15
       grad_u(1,beta,i,j,k) = grad_u(1,beta,i,j,k) + cr(beta,l)*u(perx(i+ci(l)),pery(j+cj(l)),perz(k+ck(l)))
       grad_u(2,beta,i,j,k) = grad_u(2,beta,i,j,k) + cr(beta,l)*v(perx(i+ci(l)),pery(j+cj(l)),perz(k+ck(l)))
       grad_u(3,beta,i,j,k) = grad_u(3,beta,i,j,k) + cr(beta,l)*w(perx(i+ci(l)),pery(j+cj(l)),perz(k+ck(l)))
      end do
     end if

     grad_u(1,beta,i,j,k) = grad_u(1,beta,i,j,k) / (10.d0*ds)
     grad_u(2,beta,i,j,k) = grad_u(2,beta,i,j,k) / (10.d0*ds)
     grad_u(3,beta,i,j,k) = grad_u(3,beta,i,j,k) / (10.d0*ds)

    end do

   end do
  end do
 end do
 !------------------------------------------------
 !--------vcap------------------------------------
 do k=0,kmax
  do j=0,jmax
   do i=0,imax

    do alpha=1,3
     vcap(alpha,i,j,k) = 0.d0

     do beta = 1,3
      vcap(alpha,i,j,k) = vcap(alpha,i,j,k) &
                         + 1.d0/rho(i,j,k)*grad_mu(beta,i,j,k)&
                                          *(grad_u(alpha,beta,i,j,k) + grad_u(beta,alpha,i,j,k))*ds
     end do

    end do

   end do
  end do
 end do
 !------------------------------------------------
 !-------lap_u etc. (all sides are periodic)-------
 do k=0,kmax
  do j=0,jmax
   do i=0,imax
    lap_u(i,j,k) = - 14.d0*u(i,j,k)
    lap_v(i,j,k) = - 14.d0*v(i,j,k)
    lap_w(i,j,k) = - 14.d0*w(i,j,k)

     !bulk
     if((i/=0).and.(i/=imax).and.(j/=0).and.(j/=jmax).and.(k/=0).and.(k/=kmax)) then
      do l=2,15
       lap_u(i,j,k) = lap_u(i,j,k) + u(i+ci(l),j+cj(l),k+ck(l))
       lap_v(i,j,k) = lap_v(i,j,k) + v(i+ci(l),j+cj(l),k+ck(l))
       lap_w(i,j,k) = lap_w(i,j,k) + w(i+ci(l),j+cj(l),k+ck(l))
      end do

     !boundary (periodic)
     else
      do l=2,15
       lap_u(i,j,k) = lap_u(i,j,k) + u(perx(i+ci(l)),pery(j+cj(l)),perz(k+ck(l)))
       lap_v(i,j,k) = lap_v(i,j,k) + v(perx(i+ci(l)),pery(j+cj(l)),perz(k+ck(l)))
       lap_w(i,j,k) = lap_w(i,j,k) + w(perx(i+ci(l)),pery(j+cj(l)),perz(k+ck(l)))
      end do
     end if

    lap_u(i,j,k) = lap_u(i,j,k) / (5.d0*ds**2)
    lap_v(i,j,k) = lap_v(i,j,k) / (5.d0*ds**2)
    lap_w(i,j,k) = lap_w(i,j,k) / (5.d0*ds**2)

   end do
  end do
 end do
 !------------------------------------------------
 !----laplap_u etc. (all sides are periodic)------
 do k=0,kmax
  do j=0,jmax
   do i=0,imax
    laplap_u(i,j,k) = - 14.d0*lap_u(i,j,k)
    laplap_v(i,j,k) = - 14.d0*lap_v(i,j,k)
    laplap_w(i,j,k) = - 14.d0*lap_w(i,j,k)

     !bulk
     if((i/=0).and.(i/=imax).and.(j/=0).and.(j/=jmax).and.(k/=0).and.(k/=kmax)) then
      do l=2,15
       laplap_u(i,j,k) = laplap_u(i,j,k) + lap_u(i+ci(l),j+cj(l),k+ck(l))
       laplap_v(i,j,k) = laplap_v(i,j,k) + lap_v(i+ci(l),j+cj(l),k+ck(l))
       laplap_w(i,j,k) = laplap_w(i,j,k) + lap_w(i+ci(l),j+cj(l),k+ck(l))
      end do

     !boundary (periodic)
     else
      do l=2,15
       laplap_u(i,j,k) = laplap_u(i,j,k) + lap_u(perx(i+ci(l)),pery(j+cj(l)),perz(k+ck(l)))
       laplap_v(i,j,k) = laplap_v(i,j,k) + lap_v(perx(i+ci(l)),pery(j+cj(l)),perz(k+ck(l)))
       laplap_w(i,j,k) = laplap_w(i,j,k) + lap_w(perx(i+ci(l)),pery(j+cj(l)),perz(k+ck(l)))
      end do
     end if

    laplap_u(i,j,k) = laplap_u(i,j,k) / (5.d0*ds**2)
    laplap_v(i,j,k) = laplap_v(i,j,k) / (5.d0*ds**2)
    laplap_w(i,j,k) = laplap_w(i,j,k) / (5.d0*ds**2)

   end do
  end do
 end do
 !------------------------------------------------
 !---------grad_rho  (all sides are periodic)-----
 do k=0,kmax
  do j=0,jmax
   do i=0,imax

    do alpha=1,3
     grad_rho(alpha,i,j,k) = 0.d0

     !bulk
     if((i/=0).and.(i/=imax).and.(j/=0).and.(j/=jmax).and.(k/=0).and.(k/=kmax)) then
      do l=2,15
       grad_rho(alpha,i,j,k) = grad_rho(alpha,i,j,k) + cr(alpha,l)*rho(i+ci(l),j+cj(l),k+ck(l))
      end do

     !boundary (periodic)
     else
      do l=2,15
       grad_rho(alpha,i,j,k) = grad_rho(alpha,i,j,k) &
                             + cr(alpha,l)*rho(perx(i+ci(l)),pery(j+cj(l)),perz(k+ck(l)))
      end do
     end if

     grad_rho(alpha,i,j,k) = grad_rho(alpha,i,j,k) / (10.d0*ds)

    end do

   end do
  end do
 end do
 !------------------------------------------------
 !---------normal---------------------------------
 do k=0,kmax
  do j=0,jmax
   do i=0,imax

    gtemp = dsqrt(grad_rho(1,i,j,k)**2+grad_rho(2,i,j,k)**2+grad_rho(3,i,j,k)**2)

     do alpha = 1, 3
      normal(alpha,i,j,k) = grad_rho(alpha,i,j,k)/(gtemp + epsilon)
     end do

   end do
  end do
 end do
 !------------------------------------------------
 !-----------chi  (all sides are periodic)--------
 do k=0,kmax
  do j=0,jmax
   do i=0,imax
    
     chi(i,j,k) = 0.d0

     !bulk
     if((i/=0).and.(i/=imax).and.(j/=0).and.(j/=jmax).and.(k/=0).and.(k/=kmax)) then
      do l=2,15
       do beta = 1,3
        chi(i,j,k) = chi(i,j,k) - cr(beta,l)*normal(beta,i+ci(l),j+cj(l),k+ck(l))
       end do
      end do

     !boundary (periodic)
     else
      do l=2,15
       do beta = 1,3
        chi(i,j,k) = chi(i,j,k) - cr(beta,l)*normal(beta,perx(i+ci(l)),pery(j+cj(l)),perz(k+ck(l)))
       end do
      end do
     end if
 
     chi(i,j,k) = chi(i,j,k) / (10.d0*ds)

   end do
  end do
 end do
 !------------------------------------------------
 !--gradually applying interfacial tension---------
 if(n<start_sigma) then
  sigma_temp = 0.d0
 else if(n>=start_sigma.and.n<=end_sigma) then
  sigma_temp = dble(n-start_sigma)/dble(end_sigma-start_sigma)*sigma
 else
  sigma_temp = sigma
 end if
 !------------------------------------------------
 !-----------fsv----------------------------------
 do k=0,kmax
  do j=0,jmax
   do i=0,imax
    
    do alpha = 1,3
     fsv(alpha,i,j,k) = sigma_temp*chi(i,j,k)*grad_rho(alpha,i,j,k)*rho(i,j,k)/((rhol-rhog)*(0.5d0*(rhog+rhol)))
    end do
 
   end do
  end do
 end do
 !------------------------------------------------
 !********velocity (all sides are periodic)********
  do k=0,kmax
   do j=0,jmax
    do i=0,imax

     unext(i,j,k) = -(1.d0-Au(i,j,k))/6.d0*lambda*ds**4*laplap_u(i,j,k) + fsv(1,i,j,k)/rho(i,j,k)*ds
     vnext(i,j,k) = -(1.d0-Au(i,j,k))/6.d0*lambda*ds**4*laplap_v(i,j,k) + fsv(2,i,j,k)/rho(i,j,k)*ds
     wnext(i,j,k) = -(1.d0-Au(i,j,k))/6.d0*lambda*ds**4*laplap_w(i,j,k) + fsv(3,i,j,k)/rho(i,j,k)*ds

     !bulk
     if((i/=0).and.(i/=imax).and.(j/=0).and.(j/=jmax).and.(k/=0).and.(k/=kmax)) then
      do l=1,15
       unext(i,j,k) = unext(i,j,k) + cr(1,l)*(delta_p(l,i,j,k) + geq(l,i-ci(l), j-cj(l), k-ck(l)) &
 				     +3.d0*Au(i,j,k)*E(l)*(cr(1,l)*(u(i,j,k)-u(i-ci(l), j-cj(l), k-ck(l))) &
                                                          +cr(2,l)*(v(i,j,k)-v(i-ci(l), j-cj(l), k-ck(l))) &
                                                          +cr(3,l)*(w(i,j,k)-w(i-ci(l), j-cj(l), k-ck(l)))) &
                               +3.d0*E(l)*(cr(1,l)*vcap(1,i,j,k)+cr(2,l)*vcap(2,i,j,k)+cr(3,l)*vcap(3,i,j,k)))


       vnext(i,j,k) = vnext(i,j,k) + cr(2,l)*(delta_p(l,i,j,k) + geq(l,i-ci(l), j-cj(l), k-ck(l)) &
  				     +3.d0*Au(i,j,k)*E(l)*(cr(1,l)*(u(i,j,k)-u(i-ci(l), j-cj(l), k-ck(l))) &
                                                          +cr(2,l)*(v(i,j,k)-v(i-ci(l), j-cj(l), k-ck(l))) &
                                                          +cr(3,l)*(w(i,j,k)-w(i-ci(l), j-cj(l), k-ck(l)))) &
                               +3.d0*E(l)*(cr(1,l)*vcap(1,i,j,k)+cr(2,l)*vcap(2,i,j,k)+cr(3,l)*vcap(3,i,j,k)))

       wnext(i,j,k) = wnext(i,j,k) + cr(3,l)*(delta_p(l,i,j,k) + geq(l,i-ci(l), j-cj(l), k-ck(l)) &
				     +3.d0*Au(i,j,k)*E(l)*(cr(1,l)*(u(i,j,k)-u(i-ci(l), j-cj(l), k-ck(l))) &
                                                          +cr(2,l)*(v(i,j,k)-v(i-ci(l), j-cj(l), k-ck(l))) &
                                                          +cr(3,l)*(w(i,j,k)-w(i-ci(l), j-cj(l), k-ck(l)))) &
                               +3.d0*E(l)*(cr(1,l)*vcap(1,i,j,k)+cr(2,l)*vcap(2,i,j,k)+cr(3,l)*vcap(3,i,j,k))) 
      end do

     !boundary (periodic)
     else
      do l=1,15
       unext(i,j,k) = unext(i,j,k) &
                     + cr(1,l)*(delta_p(l,i,j,k) + geq(l,perx(i-ci(l)), pery(j-cj(l)), perz(k-ck(l))) &
 				     +3.d0*Au(i,j,k)*E(l)*(cr(1,l)*(u(i,j,k)-u(perx(i-ci(l)), pery(j-cj(l)), perz(k-ck(l)))) &
                                          +cr(2,l)*(v(i,j,k)-v(perx(i-ci(l)), pery(j-cj(l)), perz(k-ck(l)))) &
                                          +cr(3,l)*(w(i,j,k)-w(perx(i-ci(l)), pery(j-cj(l)), perz(k-ck(l))))) &
                               +3.d0*E(l)*(cr(1,l)*vcap(1,i,j,k)+cr(2,l)*vcap(2,i,j,k)+cr(3,l)*vcap(3,i,j,k))) 

       vnext(i,j,k) = vnext(i,j,k) &
                    + cr(2,l)*(delta_p(l,i,j,k) + geq(l,perx(i-ci(l)), pery(j-cj(l)), perz(k-ck(l))) &
  				     +3.d0*Au(i,j,k)*E(l)*(cr(1,l)*(u(i,j,k)-u(perx(i-ci(l)), pery(j-cj(l)), perz(k-ck(l)))) &
                                          +cr(2,l)*(v(i,j,k)-v(perx(i-ci(l)), pery(j-cj(l)), perz(k-ck(l)))) &
                                          +cr(3,l)*(w(i,j,k)-w(perx(i-ci(l)), pery(j-cj(l)), perz(k-ck(l))))) &
                               +3.d0*E(l)*(cr(1,l)*vcap(1,i,j,k)+cr(2,l)*vcap(2,i,j,k)+cr(3,l)*vcap(3,i,j,k))) 

       wnext(i,j,k) = wnext(i,j,k) &
                    + cr(3,l)*(delta_p(l,i,j,k) + geq(l,perx(i-ci(l)), pery(j-cj(l)), perz(k-ck(l))) &
				     +3.d0*Au(i,j,k)*E(l)*(cr(1,l)*(u(i,j,k)-u(perx(i-ci(l)), pery(j-cj(l)), perz(k-ck(l)))) &
                                          +cr(2,l)*(v(i,j,k)-v(perx(i-ci(l)), pery(j-cj(l)), perz(k-ck(l)))) &
                                          +cr(3,l)*(w(i,j,k)-w(perx(i-ci(l)), pery(j-cj(l)), perz(k-ck(l))))) &
                               +3.d0*E(l)*(cr(1,l)*vcap(1,i,j,k)+cr(2,l)*vcap(2,i,j,k)+cr(3,l)*vcap(3,i,j,k))) 
      end do
     end if


    end do
   end do
  end do
 !------------counter back u----------------------
!  do i=0,imax
!    do j=0,jmax
!      if (wnext(i,j,0)<0) then
!        wnext(i,j,0)=-wnext(i,j,0)
!      end if
!      if (wnext(i,j,kmax)>0) then
!        wnext(i,j,kmax)=-wnext(i,j,kmax)
!      end if
!    end do
!  end do
!
!  do j=0, jmax
!    do k=0, kmax
!      if (unext(0,j,k)<0) then
!        unext(0,j,k)=-unext(0,j,k)
!      end if
!      if (unext(imax,j,k)>0) then
!        unext(imax,j,k)=-unext(imax,j,k)
!      end if
!    end do
!  end do
!
!  do i=0,imax
!    do k=0,kmax
!      if (vnext(i,0,k)<0) then
!        vnext(i,0,k)=-vnext(i,0,k)
!      end if
!      if (vnext(i,jmax,k)>0) then
!        vnext(i,jmax,k)=-vnext(i,jmax,k)
!      end if
!    end do
!  end do
 !************************************************
 !------------updating rho and u------------------
  do k=0,kmax
   do j=0,jmax
    do i=0,imax
     rho(i,j,k) = rhonext(i,j,k)
     u(i,j,k) = unext(i,j,k)
     v(i,j,k) = vnext(i,j,k)
     w(i,j,k) = wnext(i,j,k)
    end do
   end do
  end do
 !------------------------------------------------
 !---------------------------------------------------------
 !---------------output------------------------------------
 !---------pressure difference--------------------
 phi_min = phi2
 phi_max = phi1
 veff = 0.d0

 do k=0,kmax
  do j=0,jmax
   do i=0,imax

    if(phi_min > phi(i,j,k)) then
     phi_min = phi(i,j,k)
     pg = p(i,j,k)
    end if

    if(phi_max < phi(i,j,k)) then
     phi_max = phi(i,j,k)
     pl = p(i,j,k)
    end if

    if(phi(i,j,k) >= 0.5d0*(phi1+phi2)) then
     veff = veff + ds**3
    end if

   end do
  end do
 end do

 reff = (3.d0/4.d0*veff / pi)**(1.d0/3.d0)
 dp_th = 2.d0*sigma/reff
 dp_calc = pl-pg
 err = dabs(dp_calc-dp_th)/dp_th

 write(8,'(1i10, 1x, 3f20.15)') n, dp_calc, dp_th, err
 !-------------------------------------------
 !-----effective droplet volume, radius------
 write(9,'(1i10, 1x, 2f20.15)') n, veff/vinput, reff/rinput
 !-------------------------------------------
 !distribution of phi, density, pressure, and velocity for gnuplot
  if(mod(n, step/20000) == 0) then

    num='      '
    write(num,'(1i6)') n
    do j=1,6
     if(num(j:j)==' ') num(j:j)='0'
    end do

    open(11,file='flow'//num//'.txt')


    do j=0, jmax
     do i=0, imax
      do k=0, kmax
      write(11,'(3i20, 1x, 6f20.8)') i, j, k, phi(i,j,k), rho(i,j,k), p(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k)
      end do
      write(11,*)
     end do
    write(11,*) 
    end do

    close(11)
  end if
 !------------------------------------------
 !----------------------------------------------------------


 END DO
 !====================================================================

 close(8)
 close(9)


 contains

 !====================================================================function for periodic condition
 !-------------------------------------
 function perx(i)
  integer:: perx
  integer, intent(in):: i

  if(i<0) then
   perx = i+imax
  else if(i>imax) then
   perx = i-imax
  else
   perx = i
  end if

 end function perx
 !-------------------------------------
 !-------------------------------------
 function pery(j)
  integer:: pery
  integer, intent(in):: j

  if(j<0) then
   pery = j+jmax
  else if(j>jmax) then
   pery = j-jmax
  else
   pery = j
  end if

 end function pery
 !-------------------------------------
 !-------------------------------------
 function perz(k)
  integer:: perz
  integer, intent(in):: k

  if(k<0) then
   perz = k+kmax
  else if(k>kmax) then
   perz = k-kmax
  else
   perz = k
  end if

 end function perz
 !-------------------------------------

 !====================================================================

end program laplace
