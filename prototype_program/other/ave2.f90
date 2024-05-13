program main
implicit none
    integer,parameter:: xmax = 256 !ｘ方向格子数
    integer,parameter:: ymax = 256 !ｙ方向格子数
    integer,parameter:: zmax = 256 !ｚ方向格子数
    ! integer,parameter:: xmax = 64 !ｘ方向格子数
    ! integer,parameter:: ymax = 64 !ｙ方向格子数
    ! integer,parameter:: zmax = 64 !ｚ方向格子数
    ! integer,parameter:: xmax = 32 !ｘ方向格子数
    ! integer,parameter:: ymax = 32 !ｙ方向格子数
    ! integer,parameter:: zmax = 32 !ｚ方向格子数
    ! integer,parameter:: xmax = 16 !ｘ方向格子数
    ! integer,parameter:: ymax = 16 !ｙ方向格子数
    ! integer,parameter:: zmax = 16 !ｚ方向格子数
    real(8),allocatable ::  f(:)
    real(8) dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7
    integer i
    real(8) dkx, dky, dkz, dk, kmax_abs
    real(8),parameter:: pi = acos(-1.0d0) !円周率
    real(8),parameter:: Dv = 127.5d0
    real(8),parameter:: umax = 0.1d0
    ! real(8),parameter:: nu = 5.3125d-4

    real(8),parameter:: Dv2 = 7.5d0
    real(8),parameter:: umax2 = 0.05d0
    ! real(8),parameter:: nu2 = 6.458d-5

    dkx = 2.0d0*pi/dble(xmax)
    dky = 2.0d0*pi/dble(ymax)
    dkz = 2.0d0*pi/dble(zmax)
    dk = sqrt( (dkx)**2 + (dky)**2 + (dkz)**2 )
    kmax_abs = sqrt( (dkx*dble(xmax)/2.0d0)**2 + (dky*dble(ymax)/2.0d0)**2 + (dkz*dble(zmax)/2.0d0)**2 )
    allocate(f( floor(kmax_abs / dk) + 1))

    open(10, file="enesupe_re12000.d")
    ! open(10, file="DNS.d")
    do i = 1, size(f)
        read(10,*) dummy1, f(i)
    enddo
    close(10)

    open(11,file="enesupe_re12000_2.d")
    ! open(11,file="DNS_3.d")
    do i = 1, size(f)
        write(11, "(2es24.16)") dble(i)-1.0d0,  f(i)
        ! write(11, "(2es24.16)") (dble(i-1) + 0.5d0)*dk, f(i) / (umax*umax*Dv)
    enddo
end program