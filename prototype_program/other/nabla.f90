program main
    implicit none
    real(8),parameter :: ds = 1.0d0
    integer,parameter :: xmax = 128
    integer,parameter :: ymax = 128
    integer,parameter :: zmax = 128
    !粒子速度（整数）
    integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
    integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
    integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)
    real(8),parameter:: pi = acos(-1.0d0) !円周率
    real(8),parameter :: dx = 2.0d0*pi/dble(xmax)
    real(8):: cr(1:3, 1:15)  !粒子速度（実数）
    real(8) u(-1:xmax+1, -1:ymax+1, -1:zmax+1)
    real(8) grad_u_lbm(1:3, 0:xmax, 0:ymax, 0:zmax), grad_u_dif(0:xmax, 0:ymax, 0:zmax)

    integer xi, yi, zi, i, beta
    real(8) err, err2, err_sum, err_sum2

    do i=1,15
        cr(1,i) = dble(cx(i))
        cr(2,i) = dble(cy(i))
        cr(3,i) = dble(cz(i))
    enddo

    !初期条件
    do zi = -1, zmax+1
        do yi = -1, ymax+1
            do xi = -1, xmax+1
                u(xi,yi,zi) = sin(2.0d0*pi*dble(xi)/dble(xmax))*cos(2.0d0*pi*dble(yi)/dble(ymax))
            enddo
        enddo
    enddo


    !中心差分
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                ! grad_u_dif(xi,yi,zi) = (u(xi+1,yi,zi) - u(xi-1,yi,zi)) / (2.0d0*ds)
                grad_u_dif(xi,yi,zi) = (u(xi,yi+1,zi) - u(xi,yi-1,zi)) / (2.0d0*ds)
            enddo
        enddo
    enddo
    !LBM
    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                do beta=1,3
                    grad_u_lbm(beta,xi,yi,zi) = 0.0d0
                    do i=2,15
                        grad_u_lbm(beta,xi,yi,zi) = grad_u_lbm(beta,xi,yi,zi)+cr(beta,i)*u(xi+cx(i),yi+cy(i),zi+cz(i))
                    enddo
                    grad_u_lbm(beta,xi,yi,zi) = grad_u_lbm(beta,xi,yi,zi)/(10.0d0*ds)
                enddo
            enddo
        enddo
    enddo

    open(10,file="./err_128_case2.d")
    zi = (zmax) / 2
    xi = (xmax) / 2
    err_sum = 0.0d0
    err_sum2 = 0.0d0
    do yi = 0, ymax
        ! err = (abs( 2.0d0*pi/dble(xmax)*cos(2.0d0*pi*dble(xi)/dble(xmax))  -  grad_u_lbm(1,xi,yi,zi) ) )/ abs( 2.0d0*pi/dble(xmax)*cos(2.0d0*pi*dble(xi)/dble(xmax)) )
        ! err2 = (abs( 2.0d0*pi/dble(xmax)*cos(2.0d0*pi*dble(xi)/dble(xmax))  -  grad_u_dif(xi,yi,zi) ) )/ abs( 2.0d0*pi/dble(xmax)*cos(2.0d0*pi*dble(xi)/dble(xmax)) )

        err = (abs( -2.0d0*pi/dble(xmax)*sin(2.0d0*pi*dble(xi)/dble(xmax))*sin(2.0d0*pi*dble(yi)/dble(ymax))  -  grad_u_lbm(2,xi,yi,zi) ) )/ abs( -2.0d0*pi/dble(xmax)*sin(2.0d0*pi*dble(xi)/dble(xmax))*sin(2.0d0*pi*dble(yi)/dble(ymax)) )
        err2 = (abs( -2.0d0*pi/dble(xmax)*sin(2.0d0*pi*dble(xi)/dble(xmax))*sin(2.0d0*pi*dble(yi)/dble(ymax))  -  grad_u_dif(xi,yi,zi) ) )/ abs( -2.0d0*pi/dble(xmax)*sin(2.0d0*pi*dble(xi)/dble(xmax))*sin(2.0d0*pi*dble(yi)/dble(ymax)) )
        err_sum = err_sum + err
        err_sum2 = err_sum2 + err2

        write(10,"(6es24.16)") 2.0d0*pi*dble(yi)/dble(ymax), grad_u_lbm(2,xi,yi,zi), grad_u_dif(xi,yi,zi), -2.0d0*pi/dble(xmax)*sin(2.0d0*pi*dble(xi)/dble(xmax))*sin(2.0d0*pi*dble(yi)/dble(ymax)), err, err2
    enddo
    close(10)

    open(11,file="./err_sum.d")
    write(11,"(1es24.16)") err_sum/(dble(xmax+1)), err_sum2/(dble(xmax+1))
    close(11)

end program