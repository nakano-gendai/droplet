module globals
        real(8),parameter:: ds = 1.0d0 !格子間隔（lattice unit）
        integer,parameter:: xmax = 3 !ｘ方向格子数
        integer,parameter:: ymax = 3 !ｙ方向格子数
        integer,parameter:: zmax = 99 !ｚ方向格子数
        integer,parameter:: step = 100000000
        integer i, n, xi, yi, zi
        real(8) err,ut,err_up,err_down
    
        !粒子速度（整数）
        integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
        integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
        integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)
        real(8):: cr(1:3, 1:15)  !粒子速度（実数）
        
        !パラメータ
        real(8),parameter:: rho = 1.0d0 !液体の密度
        real(8),parameter:: nu = 0.05d0 !液体の動粘度
        real(8),parameter:: U(3) = (/0.01d0, 0.0d0, 0.0d0/) !壁の速度
        real(8),parameter:: dp = 0.0d0
        real(8) c
        real(8),parameter:: tau = 0.5d0 + 3.0d0 * nu / ds
        real(8),parameter:: ep = 1.0d0 / tau
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
            do zi=0,zmax
                do yi=0,ymax
                    do xi=0,xmax
                        p(xi,yi,zi) = rho / 3.0d0
                        u1(xi,yi,zi) = 0.0d0
                        u2(xi,yi,zi) = 0.0d0
                        u3(xi,yi,zi) = 0.0d0
                    enddo
                enddo
            enddo
            do zi=0,zmax
                do yi=0,ymax
                    do xi=0,xmax
                        do i =1,15
                            f(i,xi,yi,zi) = E(i)*(3.0d0*p(xi,yi,zi)+3.0d0*(cr(1,i)*u1(xi,yi,zi)+cr(2,i)*u2(xi,yi,zi) &
                            +cr(3,i)*u3(xi,yi,zi))+4.5d0*(cr(1,i)*u1(xi,yi,zi)+cr(2,i)*u2(xi,yi,zi)+cr(3,i)*u3(xi,yi,zi))**2 &
                            -1.5d0*(u1(xi,yi,zi)**2+u2(xi,yi,zi)**2+u3(xi,yi,zi)**2))
                        enddo
                    enddo
                enddo
            enddo
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

        subroutine periodic_pressure(fnext)
            real(8),intent(inout) :: fnext(1:15, -1:xmax+1, -1:ymax+1, -1:zmax+1)
            do zi=0,zmax
                do yi=0,ymax
                    c = dp - (fnext(1,0,yi,zi)-fnext(1,xmax,yi,zi) + fnext(3,0,yi,zi)-fnext(3,xmax,yi,zi) &
                            + fnext(4,0,yi,zi)-fnext(4,xmax,yi,zi) + fnext(6,0,yi,zi)-fnext(6,xmax,yi,zi) &
                            + fnext(7,0,yi,zi)-fnext(7,xmax,yi,zi)) / 3.0d0
                    !左境界
                    fnext(2,0,yi,zi) = fnext(2,xmax,yi,zi) + c
                    fnext(8,0,yi,zi) = fnext(8,xmax,yi,zi) + 0.125d0 * c
                    fnext(10,0,yi,zi) = fnext(10,xmax,yi,zi) + 0.125d0 * c
                    fnext(11,0,yi,zi) = fnext(11,xmax,yi,zi) + 0.125d0 * c
                    fnext(13,0,yi,zi) = fnext(13,xmax,yi,zi) + 0.125d0 * c
                    !右境界
                    fnext(5,xmax,yi,zi) = fnext(5,0,yi,zi) - c
                    fnext(9,xmax,yi,zi) = fnext(9,0,yi,zi) - 0.125d0 * c
                    fnext(12,xmax,yi,zi) = fnext(12,0,yi,zi) - 0.125d0 * c
                    fnext(14,xmax,yi,zi) = fnext(14,0,yi,zi) - 0.125d0 * c
                    fnext(15,xmax,yi,zi) = fnext(15,0,yi,zi) - 0.125d0 * c
                enddo
            enddo
        end subroutine periodic_pressure

        subroutine bounce_back(cr,fnext)
            real(8),intent(in) :: cr(1:3,1:15)
            real(8),intent(inout) :: fnext(1:15, -1:xmax+1, -1:ymax+1, -1:zmax+1)
            !下の壁（zi=0）
            zi = 0
            do yi=0,ymax
                do xi=0,xmax
                    fnext(4,xi,yi,zi) = fnext(7,xi,yi,zi) - 6.0d0 * E(7) * ((-U(1))*cr(1,7)+(-U(2))*cr(2,7)+(-U(3))*cr(3,7))
                    fnext(8,xi,yi,zi) = fnext(12,xi,yi,zi) - 6.0d0 * E(12) * ((-U(1))*cr(1,12)+(-U(2))*cr(2,12)+(-U(3))*cr(3,12))
                    fnext(9,xi,yi,zi) = fnext(13,xi,yi,zi) - 6.0d0 * E(13) * ((-U(1))*cr(1,13)+(-U(2))*cr(2,13)+(-U(3))*cr(3,13))
                    fnext(10,xi,yi,zi) = fnext(14,xi,yi,zi) - 6.0d0 * E(14) * ((-U(1))*cr(1,14)+(-U(2))*cr(2,14)+(-U(3))*cr(3,14))
                    fnext(15,xi,yi,zi) = fnext(11,xi,yi,zi) - 6.0d0 * E(11) * ((-U(1))*cr(1,11)+(-U(2))*cr(2,11)+(-U(3))*cr(3,11))
                enddo
            enddo
            !上の壁（zi=zmax）
            zi = zmax
            do yi=0,ymax
                do xi=0,xmax
                    fnext(7,xi,yi,zi) = fnext(4,xi,yi,zi) - 6.0d0 * E(4) * ((U(1))*cr(1,4)+(U(2))*cr(2,4)+(U(3))*cr(3,4))
                    fnext(11,xi,yi,zi) = fnext(15,xi,yi,zi) - 6.0d0 * E(15) * ((U(1))*cr(1,15)+(U(2))*cr(2,15)+(U(3))*cr(3,15))
                    fnext(12,xi,yi,zi) = fnext(8,xi,yi,zi) - 6.0d0 * E(8) * ((U(1))*cr(1,8)+(U(2))*cr(2,8)+(U(3))*cr(3,8))
                    fnext(13,xi,yi,zi) = fnext(9,xi,yi,zi) - 6.0d0 * E(9) * ((U(1))*cr(1,9)+(U(2))*cr(2,9)+(U(3))*cr(3,9))
                    fnext(14,xi,yi,zi) = fnext(10,xi,yi,zi) - 6.0d0 * E(10) * ((U(1))*cr(1,10)+(U(2))*cr(2,10)+(U(3))*cr(3,10))
                enddo
            enddo
        end subroutine bounce_back

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
            do zi=0,zmax
                do yi=0,ymax
                    do xi=0,xmax
                        do i =1,15
                            u1(xi,yi,zi) = u1(xi,yi,zi) + f(i,xi,yi,zi) * cr(1,i)
                            u2(xi,yi,zi) = u2(xi,yi,zi) + f(i,xi,yi,zi) * cr(2,i)
                            u3(xi,yi,zi) = u3(xi,yi,zi) + f(i,xi,yi,zi) * cr(3,i)
                        enddo
                    enddo
                enddo
            enddo
        end subroutine velocity_cal

        subroutine output(u1,u2,u3,p,f)
            real(8),intent(in) :: u1(-1:xmax+1, -1:ymax+1, -1:zmax+1), u2(-1:xmax+1, -1:ymax+1, -1:zmax+1),&
                                    u3(-1:xmax+1, -1:ymax+1, -1:zmax+1)
            real(8),intent(in) :: p(-1:xmax+1, -1:ymax+1, -1:zmax+1), f(1:15,-1:xmax+1,-1:ymax+1,-1:zmax+1)
            open(10,file = "u_z.d")
            do zi=0,zmax
                do yi=0,ymax
                    do xi=0,xmax
                        ! write(10,"(7es16.8)") dble(xi)/dble(xmax),dble(yi)/dble(ymax),dble(zi)/dble(zmax), &
                        ! u1(xi,yi,zi)/sqrt(U(1)**2+U(2)**2+U(3)**2),u2(xi,yi,zi)/sqrt(U(1)**2+U(2)**2+U(3)**2) &
                        ! ,u3(xi,yi,zi)/sqrt(U(1)**2+U(2)**2+U(3)**2),p(xi,yi,zi)
                        write(10,*) xi,yi,zi,u1(xi,yi,zi),u2(xi,yi,zi),u3(xi,yi,zi),p(xi,yi,zi)
                    enddo
                enddo
            enddo
            close(10)

            open(11,file="err_z.d")
            err_up = 0.0d0
            err_down = 0.0d0
            err = 0.0d0
            do zi=0,zmax
                ut = 2.0d0 * sqrt(U(1)**2+U(2)**2+U(3)**2) * dble(zi) / dble(zmax) - sqrt(U(1)**2+U(2)**2+U(3)**2)
                err_up = err_up + abs((ut) - u1(1,1,zi))
                err_down = err_down + abs(ut)
                err = abs((ut) - u1(1,1,zi)) / abs(ut)
                write(11,*) zi,u1(1,1,zi),ut,err
            enddo
            err = err_up / err_down
            write(11,*) err
            close(11)

            open(12,file = "udis_z.d")
            do zi=0,zmax
                write(12,*) zi, u1(2,2,zi)
            enddo
            close(12)

            open(14,file="fz.d")
            do zi=0,zmax
                do yi=0,ymax
                    do xi=0,xmax
                        do i=1,15
                            write(14,*) xi,yi,zi,i,f(i,xi,yi,zi)
                        enddo
                    enddo
                enddo
            enddo
            close(14)
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
        !左右（Y-Z平面）周期境界条件（dp=0）
        call periodic_pressure(fnext)
        !壁（zi=0,zmax）
        call bounce_back(cr,fnext)
        !速度分布関数の更新
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
    do zi=0,zmax
        write(*,"(7es16.8)") dble(2)/dble(xmax),dble(2)/dble(ymax),dble(zi)/dble(zmax), &
        u1(1,1,zi)/sqrt(U(1)**2+U(2)**2+U(3)**2),u2(1,1,zi)/sqrt(U(1)**2+U(2)**2+U(3)**2), &
        u3(1,1,zi)/sqrt(U(1)**2+U(2)**2+U(3)**2),p(1,1,zi)
        ! do yi=0,ymax
        !     do xi=0,xmax
        !         do i=1,15
        !             write(*,*) xi,yi,zi,i,f(i,xi,yi,zi)
        !         enddo
        !     enddo
        ! enddo
    enddo
end program main
