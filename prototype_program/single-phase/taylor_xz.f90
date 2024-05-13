module globals
        real(8),parameter:: ds = 1.0d0 !格子間隔（lattice unit）
        integer,parameter:: xmax = 39 !ｘ方向格子数
        integer,parameter:: ymax = 1 !ｙ方向格子数
        integer,parameter:: zmax = 39 !ｚ方向格子数
        integer,parameter:: step = 1
        integer i, n, xi, yi, zi
        real(8) err_u1,err_up_u1,err_down_u1,err_u3,err_up_u3,err_down_u3, &
                err_p,err_up_p,err_down_p
        real(8) ut1,ut3,pt
    
        !粒子速度（整数）
        integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
        integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
        integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)
        real(8):: cr(1:3, 1:15)  !粒子速度（実数）
        
        !パラメータ
        real(8),parameter:: rho = 1.0d0 !液体の密度
        real(8),parameter:: tau = 0.8d0
        real(8),parameter:: nu = (tau - 0.5d0)/3.0d0 !液体の動粘度
        real(8),parameter:: u0 = 0.01d0 !流速の初期値
        real(8),parameter:: p0 = rho / 3.0d0 !圧力の初期値
        real(8),parameter:: ep = 1.0d0 / tau
        real(8),parameter:: pi = acos(-1.0d0)
        real(8),parameter:: k1 = 2.0d0*pi/dble(xmax+1)
        real(8),parameter:: k2 = 2.0d0*pi/dble(ymax+1)
        real(8),parameter:: k3 = 2.0d0*pi/dble(zmax+1)
        real(8),parameter:: E(15) = (/ 2.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, &
                                    1.0d0/9.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, &
                                    1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0 /)
    contains
        subroutine parv(cx,cy,cz,cr)
            integer, intent(in):: cx(15),cy(15),cz(15)
            real(8), intent(out):: cr(1:3,1:15)
                do i=1,15
                    cr(1,i) = dble(cx(i))
                    cr(2,i) = dble(cy(i))
                    cr(3,i) = dble(cz(i))
                enddo
        end subroutine parv
    
        subroutine initial(p,u1,u2,u3,f)
            real(8), intent(inout):: p(-1:xmax+1, -1:ymax+1, -1:zmax+1),f(1:15, -1:xmax+1, -1:ymax+1, -1:zmax+1)
            real(8), intent(inout):: u1(-1:xmax+1, -1:ymax+1, -1:zmax+1),u2(-1:xmax+1, -1:ymax+1, -1:zmax+1), &
                                        u3(-1:xmax+1, -1:ymax+1, -1:zmax+1)
            open(88,file="ini_u_ori.d")
            do zi=0,zmax
                do yi=0,ymax
                    do xi=0,xmax
                        p(xi,yi,zi) = p0 - rho * u0**2 / 4.0d0 * (cos(2.0d0*k1*dble(xi)) + &
                         k1**2 / k3**2 * cos(2.0d0*k3*dble(zi)))
                        u1(xi,yi,zi) = -u0*cos(k1*dble(xi))*sin(k3*dble(zi))
                        u2(xi,yi,zi) = 0.0d0
                        u3(xi,yi,zi) = k1/k3*u0*sin(k1*dble(xi))*cos(k3*dble(zi))
                        write(88,"(7es23.16)") dble(xi),dble(yi),dble(zi),u1(xi,yi,zi),u2(xi,yi,zi),&
                                    u3(xi,yi,zi),p(xi,yi,zi)
                    enddo
                enddo
            enddo
            close(88)
            open(89,file="ini_ori.d")
            do zi=0,zmax
                do yi=0,ymax
                    do xi=0,xmax
                        do i =1,15
                            f(i,xi,yi,zi) = E(i)*(3.0d0*p(xi,yi,zi)+3.0d0*(cr(1,i)*u1(xi,yi,zi)+cr(2,i)*u2(xi,yi,zi) &
                            +cr(3,i)*u3(xi,yi,zi))+4.5d0*(cr(1,i)*u1(xi,yi,zi)+cr(2,i)*u2(xi,yi,zi)+cr(3,i)*u3(xi,yi,zi))**2 &
                            -1.5d0*(u1(xi,yi,zi)**2+u2(xi,yi,zi)**2+u3(xi,yi,zi)**2))
                            write(89,"(5es23.16)") dble(xi),dble(yi),dble(zi),dble(i),f(i,xi,yi,zi)
                        enddo
                    enddo
                enddo
            enddo
            close(89)
        end subroutine initial

        subroutine periodic(u1,u2,u3,p,f)
            real(8), intent(inout):: p(-1:xmax+1, -1:ymax+1, -1:zmax+1),f(1:15, -1:xmax+1, -1:ymax+1, -1:zmax+1)
            real(8), intent(inout):: u1(-1:xmax+1, -1:ymax+1, -1:zmax+1),u2(-1:xmax+1, -1:ymax+1, -1:zmax+1), &
                                        u3(-1:xmax+1, -1:ymax+1, -1:zmax+1)
            do zi=0,zmax
                do yi=-1,ymax+1
                    do xi=0,xmax
                        u1(xi,-1,zi) = u1(xi,ymax,zi)
                        u2(xi,-1,zi) = u2(xi,ymax,zi)
                        u3(xi,-1,zi) = u3(xi,ymax,zi)
                        p(xi,-1,zi) = p(xi,ymax,zi)
                        u1(xi,ymax+1,zi) = u1(xi,0,zi)
                        u2(xi,ymax+1,zi) = u2(xi,0,zi)
                        u3(xi,ymax+1,zi) = u3(xi,0,zi)
                        p(xi,ymax+1,zi) = p(xi,0,zi)
                        u1(-1,yi,zi) = u1(xmax,yi,zi)
                        u2(-1,yi,zi) = u2(xmax,yi,zi)
                        u3(-1,yi,zi) = u3(xmax,yi,zi)
                        p(-1,yi,zi) = p(xmax,yi,zi)
                        u1(xmax+1,yi,zi) = u1(0,yi,zi)
                        u2(xmax+1,yi,zi) = u2(0,yi,zi)
                        u3(xmax+1,yi,zi) = u3(0,yi,zi)
                        p(xmax+1,yi,zi) = p(0,yi,zi)

                        do i=1,15
                            f(i,xi,-1,zi) = f(i,xi,ymax,zi)
                            f(i,xi,ymax+1,zi) = f(i,xi,0,zi)
                            f(i,-1,yi,zi) = f(i,xmax,yi,zi)
                            f(i,xmax+1,yi,zi) = f(i,0,yi,zi)
                        enddo
                    enddo
                enddo
            enddo
            
            do yi=-1,ymax+1
                do xi=-1,xmax+1
                    u1(xi,yi,-1) = u1(xi,yi,zmax)
                    u2(xi,yi,-1) = u2(xi,yi,zmax)
                    u3(xi,yi,-1) = u3(xi,yi,zmax)

                    u1(xi,yi,zmax+1) = u1(xi,yi,0)
                    u2(xi,yi,zmax+1) = u2(xi,yi,0)
                    u3(xi,yi,zmax+1) = u3(xi,yi,0)
                    
                    p(xi,yi,-1) = p(xi,yi,zmax)
                    p(xi,yi,zmax+1) = p(xi,yi,0)

                    do i=1,15
                        f(i,xi,yi,-1) = f(i,xi,yi,zmax)
                        f(i,xi,yi,zmax+1) = f(i,xi,yi,0)
                    enddo
                enddo
            enddo
        end subroutine periodic

        subroutine renew(f,fnext)
            real(8),intent(in) :: fnext(1:15,-1:xmax+1,-1:ymax+1,-1:zmax+1)
            real(8),intent(out) :: f(1:15,-1:xmax+1,-1:ymax+1,-1:zmax+1)
            do zi=0,zmax
                do yi=0,ymax
                    do xi=0,xmax
                        do i=1,15
                            f(i,xi,yi,zi) = fnext(i,xi,yi,zi)
                        enddo
                    enddo
                enddo
            enddo
        end subroutine renew
    
        subroutine reset(p,u1,u2,u3)
            real(8),intent(inout) :: p(-1:xmax+1, -1:ymax+1, -1:zmax+1)
            real(8),intent(inout) :: u1(-1:xmax+1, -1:ymax+1, -1:zmax+1),u2(-1:xmax+1, -1:ymax+1, -1:zmax+1),&
                                        u3(-1:xmax+1, -1:ymax+1, -1:zmax+1)
            do zi=0,zmax
                do yi=0,ymax
                    do xi=0,xmax
                        p(xi,yi,zi) = 0.0d0
                        u1(xi,yi,zi) = 0.0d0
                        u2(xi,yi,zi) = 0.0d0
                        u3(xi,yi,zi) = 0.0d0
                    enddo
                enddo
            enddo
        end subroutine reset
    
        subroutine pressure_cal(f,p)
            real(8),intent(in) :: f(1:15,-1:xmax+1,-1:ymax+1,-1:zmax+1)
            real(8),intent(out) :: p(-1:xmax+1, -1:ymax+1, -1:zmax+1)
            do zi=0,zmax
                do yi=0,ymax
                    do xi=0,xmax
                        do i =1,15
                            p(xi,yi,zi) = p(xi,yi,zi) + f(i,xi,yi,zi)
                        enddo
                        p(xi,yi,zi) = p(xi,yi,zi) / 3.0d0
                    enddo
                enddo
            enddo
        end subroutine pressure_cal
    
        subroutine velocity_cal(f,u1,u2,u3)
            real(8),intent(in) :: f(1:15,-1:xmax+1,-1:ymax+1,-1:zmax+1)
            real(8),intent(out) :: u1(-1:xmax+1, -1:ymax+1, -1:zmax+1),u2(-1:xmax+1, -1:ymax+1, -1:zmax+1),&
                                    u3(-1:xmax+1, -1:ymax+1, -1:zmax+1)
            open(91,file="sum.d")
            do zi=0,zmax
                do yi=0,ymax
                    do xi=0,xmax
                        do i =1,15
                            u1(xi,yi,zi) = u1(xi,yi,zi) + f(i,xi,yi,zi) * cr(1,i)
                            u2(xi,yi,zi) = u2(xi,yi,zi) + f(i,xi,yi,zi) * cr(2,i)
                            u3(xi,yi,zi) = u3(xi,yi,zi) + f(i,xi,yi,zi) * cr(3,i)
                            write(91,*) xi,yi,zi,i,u1(xi,yi,zi)
                        enddo
                    enddo
                enddo
            enddo
            close(91)
        end subroutine velocity_cal

        subroutine output(u1,u2,u3,p,f)
            real(8),intent(in) :: u1(-1:xmax+1, -1:ymax+1, -1:zmax+1), u2(-1:xmax+1, -1:ymax+1, -1:zmax+1),&
                                    u3(-1:xmax+1, -1:ymax+1, -1:zmax+1)
            real(8),intent(in) :: p(-1:xmax+1, -1:ymax+1, -1:zmax+1),f(1:15, -1:xmax+1, -1:ymax+1, -1:zmax+1)
        open(94,file="para_original.d")
        write(94,"(5es23.16)") k1,k2,k3,pi,p0
        close(94)
            open(10,file = "taylor_u.d")
            do zi=0,zmax
                do yi=0,ymax
                    do xi=0,xmax
                        write(10, "(7es23.16)") dble(xi),dble(yi),dble(zi),u1(xi,yi,zi),u2(xi,yi,zi),u3(xi,yi,zi),p(xi,yi,zi)
                    enddo
                enddo
            enddo
            close(10)

            open(11,file = "taylor_theory.d")
            n = step
            do zi=0,zmax
                do yi=0,ymax
                    do xi=0,xmax
                        ut1 = -u0*exp(-nu*dble(n)*(k1**2+k3**2))*cos(k1*dble(xi))*sin(k3*dble(zi))
                        ut3 = k1/k3*u0*exp(-nu*dble(n)*(k1**2+k3**2))*sin(k1*dble(xi))*cos(k3*dble(zi))
                        pt = p0 -rho/4.0d0*exp(-2.0d0*nu*dble(n)*(k1**2+k3**2))*(u0**2*cos(2.0d0*k1*dble(xi))+ &
                        k1**2/k3**2*u0**2*cos(2.0d0*k3*dble(zi)))
                        write(11,*) xi,yi,zi,ut1,ut2,pt
                    enddo
                enddo
            enddo
            close(11)

            open(12,file = "taylor_err.d")
            write(12,*) "err_u1   err_u3   err_p"
            err_u1 = 0.0d0
            err_up_u1 = 0.0d0
            err_down_u1 = 0.0d0
            err_u3 = 0.0d0
            err_up_u3 = 0.0d0
            err_down_u3 = 0.0d0
            err_p = 0.0d0
            err_up_p = 0.0d0
            err_down_p = 0.0d0
            n = step
            yi = 1
            do zi=0,zmax
                do xi=0,xmax
                    ut1 = -u0*exp(-nu*dble(n)*(k1**2+k3**2))*cos(k1*dble(xi))*sin(k3*dble(zi))
                    ut3 = k1/k3*u0*exp(-nu*dble(n)*(k1**2+k3**2))*sin(k1*dble(xi))*cos(k3*dble(zi))
                    pt = p0 -rho/4.0d0*exp(-2.0d0*nu*dble(n)*(k1**2+k3**2))*(u0**2*cos(2.0d0*k1*dble(xi))+ &
                    k1**2/k3**2*u0**2*cos(2.0d0*k3*dble(zi)))
                    err_up_u1 = err_up_u1 + abs(ut1 - u1(xi,yi,zi))
                    err_down_u1 = err_down_u1 + abs(ut1)
                    err_up_u3 = err_up_u3 + abs(ut3 - u3(xi,yi,zi))
                    err_down_u3 = err_down_u3 + abs(ut3)
                    err_up_p = err_up_p + abs(pt - p(xi,yi,zi))
                    err_down_p = err_down_p + abs(pt)
                enddo
            enddo
            err_u1 = (err_up_u1 / err_down_u1)
            err_u3 = (err_up_u3 / err_down_u3)
            err_p = (err_up_p / err_down_p)
            write(12,*) err_u1,err_u3,err_p
            close(12)

            open(64, file = 'fp_normal.d')   
            do zi=0,zmax
                do yi=0,ymax
                    do xi=0,xmax
                        do i=1,15
                            write(64, "(5es23.16)") dble(i),dble(xi),dble(yi),dble(zi),f(i,xi,yi,zi)
                        enddo
                    enddo
                enddo
            enddo
            close(64)

            open(65, file = 'data1.d')    
            do zi=0,zmax
                ! do yi=0,ymax
                    do xi=0,xmax
                        ! do i=1,15
                            write(65, "(1es24.17)") u1(xi,0,zi)
                        ! enddo
                    enddo
                ! enddo
            enddo
            close(65)

        end subroutine output
    
        function feq(i,p, u1, u2, u3) result(feqr)
            integer,intent(in) :: i
            real(8),intent(in) :: p,u1,u2,u3
            real(8) feqr
            feqr = E(i)*(3.0d0*p + 3.0d0*(u1*dble(cx(i)) + u2*dble(cy(i)) + u3*dble(cz(i))) &
                    + 4.5d0*(u1*dble(cx(i)) + u2*dble(cy(i)) + u3*dble(cz(i)))**2 -1.5d0*(u1**2 + u2**2 + u3**2))
        end function feq
    
end module globals
    
program main
use globals
    implicit none
    integer in,jn,kn
    real(8) f(1:15, -1:xmax+1, -1:ymax+1, -1:zmax+1), fnext(1:15, -1:xmax+1, -1:ymax+1, -1:zmax+1) !速度分布関数
    real(8) p(-1:xmax+1, -1:ymax+1, -1:zmax+1)  !圧力
    real(8) u1(-1:xmax+1, -1:ymax+1, -1:zmax+1), u2(-1:xmax+1, -1:ymax+1, -1:zmax+1), u3(-1:xmax+1, -1:ymax+1, -1:zmax+1) !流速
    call parv(cx,cy,cz,cr)
    !初期値の設定
    call initial(p,u1,u2,u3,f)
!===========時間発展=====================
    DO n=1,step
        call periodic(u1,u2,u3,p,f)
        !次の時刻の速度分布関数fnextを計算
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    do i=1,15
                        fnext(i,xi,yi,zi) = f(i,xi-cx(i),yi-cy(i),zi-cz(i)) - ep * (f(i,xi-cx(i),yi-cy(i),zi-cz(i)) &
                            - feq(i,p(xi-cx(i),yi-cy(i),zi-cz(i)),u1(xi-cx(i),yi-cy(i),zi-cz(i)),&
                                u2(xi-cx(i),yi-cy(i),zi-cz(i)),u3(xi-cx(i),yi-cy(i),zi-cz(i))))
                    enddo
                enddo
            enddo
        enddo
        call renew(f,fnext)
        !初期化
        call reset(p,u1,u2,u3)
        !圧力の計算
        call pressure_cal(f,p)
        !流速の計算
        call velocity_cal(f,u1,u2,u3)
        write(*,*) "step = ", n
    ENDDO
    !出力
    call output(u1,u2,u3,p,f)
end program main
