module globals
    real(8),parameter:: ds = 1.0d0 !格子間隔（lattice unit）
    integer,parameter:: xmax = 3 !ｘ方向格子数
    integer,parameter:: ymax = 29 !ｙ方向格子数
    integer,parameter:: step = 30000
    integer i, n, xi, yi
    real(8) err,ut,errs
    
    !粒子速度（整数）
    integer,parameter:: cx(9) = (/0, 1, 0, -1, 0, 1, -1, -1, 1/)
    integer,parameter:: cy(9) = (/0, 0, 1, 0, -1, 1, 1, -1, -1/)
    !粒子速度（実数）
    real(8):: cr(1:2, 1:9)
        
    !パラメータ
    real(8),parameter:: rho = 1.0d0 !液体の密度
    real(8),parameter:: nu = 1.0d0 / 6.0d0 !液体の動粘度
    real(8),parameter:: U(1:2) = (/0.01d0, 0.0d0/) !壁の速度
    real(8),parameter:: dp = 0.0d0
    real(8) c
    real(8),parameter:: tau = 0.5d0 + 3.0d0 * nu / ds
    real(8),parameter:: ep = 1.0d0 / tau
    real(8),parameter:: E(9) = (/ 4.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/36.0d0, &
                        1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0 /)

contains
    subroutine parv(cx,cy,cr)
        integer, intent(in):: cx(9),cy(9)
        real(8), intent(out):: cr(1:2,1:9)
            do i=1,9
                cr(1,i) = dble(cx(i))
                cr(2,i) = dble(cy(i))
            enddo
    end subroutine parv

    subroutine initial(p,u1,u2,f)
        real(8), intent(inout):: p(-1:xmax+1, -1:ymax+1),f(1:9, -1:xmax+1, -1:ymax+1)
        real(8), intent(inout):: u1(-1:xmax+1, -1:ymax+1),u2(-1:xmax+1, -1:ymax+1)
        do yi=0,ymax
            do xi=0,xmax
                p(xi,yi) = rho / 3.0d0
                u1(xi,yi) = 0.0d0
                u2(xi,yi) = 0.0d0
            enddo
        enddo
        do yi=0,ymax
            do xi=0,xmax
                do i =1,9
                    f(i,xi,yi) = E(i)*(3.0d0*p(xi,yi)+3.0d0*(cr(1,i)*u1(xi,yi)+cr(2,i)*u2(xi,yi)) &
                    +4.5d0*(cr(1,i)*u1(xi,yi)+cr(2,i)*u2(xi,yi))**2 &
                    -1.5d0*(u1(xi,yi)**2+u2(xi,yi)**2))
                enddo
            enddo
        enddo
    end subroutine initial

    subroutine periodic(u1,u2,p,f)
            real(8), intent(inout):: p(-1:xmax+1, -1:ymax+1),u1(-1:xmax+1, -1:ymax+1),&
                                u2(-1:xmax+1, -1:ymax+1),f(1:9, -1:xmax+1, -1:ymax+1)
            do yi=-1,ymax+1
                do xi=0,xmax
                    u1(xi,-1) = u1(xi,ymax)
                    u2(xi,-1) = u2(xi,ymax)
                    p(xi,-1) = p(xi,ymax)
                    u1(xi,ymax+1) = u1(xi,0)
                    u2(xi,ymax+1) = u2(xi,0)
                    p(xi,ymax+1) = p(xi,0)
                    u1(-1,yi) = u1(xmax,yi)
                    u2(-1,yi) = u2(xmax,yi)
                    p(-1,yi) = p(xmax,yi)
                    u1(xmax+1,yi) = u1(0,yi)
                    u2(xmax+1,yi) = u2(0,yi)
                    p(xmax+1,yi) = p(0,yi)
                    do i=1,9
                        f(i,xi,-1) = f(i,xi,ymax)
                        f(i,xi,ymax+1) = f(i,xi,0)
                        f(i,-1,yi) = f(i,xmax,yi)
                        f(i,xmax+1,yi) = f(i,0,yi)
                    enddo
                enddo
            enddo
        end subroutine periodic

    subroutine periodic_pressure(fnext)
        real(8),intent(inout) :: fnext(1:9, -1:xmax+1, -1:ymax+1)
        do yi=0,ymax
            c = dp - (fnext(1,0,yi)-fnext(1,xmax,yi) + fnext(3,0,yi)-fnext(3,xmax,yi) + fnext(5,0,yi)-fnext(5,xmax,yi))/3.0d0
            !左境界
            fnext(2,0,yi) = fnext(2,xmax,yi) + c
            fnext(6,0,yi) = fnext(6,xmax,yi) + 0.25d0 * c
            fnext(9,0,yi) = fnext(9,xmax,yi) + 0.25d0 * c
            !右境界
            fnext(4,xmax,yi) = fnext(4,0,yi) - c
            fnext(7,xmax,yi) = fnext(7,0,yi) - 0.25d0 * c
            fnext(8,xmax,yi) = fnext(8,0,yi) - 0.25d0 * c
        enddo
    end subroutine periodic_pressure

    subroutine bounce_back(cr,fnext)
        real(8),intent(in) :: cr(1:2,1:9)
        real(8),intent(inout) :: fnext(1:9, -1:xmax+1, -1:ymax+1)
        !下の壁
        yi = 0
        do xi=0,xmax
            fnext(3,xi,yi) = fnext(5,xi,yi) - 6.0d0 * E(5) * ((-U(1))*cr(1,5)+(-U(2))*cr(2,5))
            fnext(6,xi,yi) = fnext(8,xi,yi) - 6.0d0 * E(8) * ((-U(1))*cr(1,8)+(-U(2))*cr(2,8))
            fnext(7,xi,yi) = fnext(9,xi,yi) - 6.0d0 * E(9) * ((-U(1))*cr(1,9)+(-U(2))*cr(2,9))
        enddo
        !上の壁
        yi = ymax
        do xi=0,xmax
            fnext(5,xi,yi) = fnext(3,xi,yi) - 6.0d0 * E(3) * ((U(1))*cr(1,3)+(U(2))*cr(2,3))
            fnext(8,xi,yi) = fnext(6,xi,yi) - 6.0d0 * E(6) * ((U(1))*cr(1,6)+(U(2))*cr(2,6))
            fnext(9,xi,yi) = fnext(7,xi,yi) - 6.0d0 * E(7) * ((U(1))*cr(1,7)+(U(2))*cr(2,7))
        enddo
    end subroutine bounce_back

    subroutine renew(f,fnext)
        real(8),intent(in) :: fnext(1:9,-1:xmax+1,-1:ymax+1)
        real(8),intent(out) :: f(1:9,-1:xmax+1,-1:ymax+1)
        do yi=0,ymax
            do xi=0,xmax
                do i=1,9
                    f(i,xi,yi) = fnext(i,xi,yi)
                enddo
            enddo
        enddo
    end subroutine renew

    subroutine reset(p,u1,u2)
        real(8),intent(inout) :: p(-1:xmax+1, -1:ymax+1)
        real(8),intent(inout) :: u1(-1:xmax+1, -1:ymax+1),u2(-1:xmax+1, -1:ymax+1)
        do yi=0,ymax
            do xi=0,xmax
                p(xi,yi) = 0.0d0
                u1(xi,yi) = 0.0d0
                u2(xi,yi) = 0.0d0
            enddo
        enddo
    end subroutine reset

    subroutine pressure_cal(f,p)
        real(8),intent(in) :: f(1:9,-1:xmax+1,-1:ymax+1)
        real(8),intent(out) :: p(-1:xmax+1, -1:ymax+1)
        do yi=0,ymax
            do xi=0,xmax
                do i =1,9
                    p(xi,yi) = p(xi,yi) + f(i,xi,yi)
                enddo
                p(xi,yi) = p(xi,yi) / 3.0d0
            enddo
        enddo
    end subroutine pressure_cal

    subroutine velocity_cal(f,u1,u2)
        real(8),intent(in) :: f(1:9,-1:xmax+1,-1:ymax+1)
        real(8),intent(out) :: u1(-1:xmax+1, -1:ymax+1),u2(-1:xmax+1, -1:ymax+1)
        do yi=0,ymax
            do xi=0,xmax
                do i =1,9
                    u1(xi,yi) = u1(xi,yi) + f(i,xi,yi) * cr(1,i)
                    u2(xi,yi) = u2(xi,yi) + f(i,xi,yi) * cr(2,i)
                enddo
            enddo
        enddo
    end subroutine velocity_cal

    subroutine output(u1,u2,p,f)
        real(8),intent(in) :: u1(-1:xmax+1, -1:ymax+1), u2(-1:xmax+1, -1:ymax+1)
        real(8),intent(in) :: p(-1:xmax+1, -1:ymax+1), f(1:9,-1:xmax+1,-1:ymax+1)
        open(10,file = "u_y.d")
        do yi=0,ymax
            do xi=0,xmax
                ! write(10,"(7es16.8)") dble(xi)/dble(xmax),dble(yi)/dble(ymax),dble(zi)/dble(zmax), &
                ! u1(xi,yi,zi)/sqrt(U(1)**2+U(2)**2+U(3)**2),u2(xi,yi,zi)/sqrt(U(1)**2+U(2)**2+U(3)**2) &
                ! ,u3(xi,yi,zi)/sqrt(U(1)**2+U(2)**2+U(3)**2),p(xi,yi,zi)
                write(10,*) xi,yi,u1(xi,yi),u2(xi,yi),p(xi,yi)
            enddo
        enddo
        close(10)

        open(11,file="err_y.d")
        do yi=0,ymax
            ut = 2.0d0 * sqrt(U(1)**2+U(2)**2) * dble(yi) / dble(ymax) - sqrt(U(1)**2+U(2)**2)
            err = (abs(ut) - sqrt(u1(2,yi)*u1(2,yi)+u2(2,yi)*u2(2,yi))) / abs(ut)
            write(11,*) dble(2)/dble(xmax),dble(yi)/dble(ymax),err
        enddo
        close(11)

        open(12,file = "udis_y.d")
        do yi=0,ymax
            write(12,"(6es16.8)") dble(2)/dble(xmax),dble(yi)/dble(ymax), &
            u1(2,yi),u2(2,yi)
        enddo
        close(12)

        open(14,file="fy.d")
        do yi=0,ymax
            do xi=0,xmax
                do i=1,9
                    write(14,*) xi,yi,i,f(i,xi,yi)
                enddo
            enddo
        enddo
        close(14)
    end subroutine output

    function feq(i,p, u1, u2) result(feqr)
        integer,intent(in) :: i
        real(8),intent(in) :: p,u1,u2
        real(8) feqr
        feqr = E(i)*(3.0d0*p + 3.0d0*(u1*dble(cx(i)) + u2*dble(cy(i))) &
                + 4.5d0*(u1*dble(cx(i)) + u2*dble(cy(i)))**2 -1.5d0*(u1**2 + u2**2))
    end function feq
    
end module globals
    
program main
use globals
    implicit none
    real(8) f(1:9, -1:xmax+1, -1:ymax+1), fnext(1:9, -1:xmax+1, -1:ymax+1) !速度分布関数
    real(8) p(-1:xmax+1, -1:ymax+1)  !圧力
    real(8) u1(-1:xmax+1, -1:ymax+1), u2(-1:xmax+1, -1:ymax+1) !流速
    real(8) time1,time2
    call cpu_time(time1)

    call parv(cx,cy,cr)
    !初期値の設定
    call initial(p,u1,u2,f)
!===========時間発展=====================
    DO n=1,step
        call periodic(u1,u2,p,f)
        !次の時刻の速度分布関数fnextを計算
        do yi=0,ymax
            do xi=0,xmax
                do i=1,9
                    fnext(i,xi,yi) = f(i,xi-cx(i),yi-cy(i)) - ep * (f(i,xi-cx(i),yi-cy(i)) &
                        - feq(i,p(xi-cx(i),yi-cy(i)),u1(xi-cx(i),yi-cy(i)),u2(xi-cx(i),yi-cy(i))))
                enddo
            enddo
        enddo
        
        !左右 周期境界条件（dp=0）
        ! call periodic_pressure(fnext)
        call bounce_back(cr,fnext)
        ! do yi=-1,ymax+1
        !     do xi=0,xmax
        !         do i=1,9
        !             fnext(i,xi,-1) = fnext(i,xi,ymax)
        !             fnext(i,xi,ymax+1) = fnext(i,xi,0)
        !             fnext(i,-1,yi) = fnext(i,xmax,yi)
        !             fnext(i,xmax+1,yi) = fnext(i,0,yi)
        !         enddo
        !     enddo
        ! enddo
        !壁（yi=0,ymax）
        ! call bounce_back(cr,fnext)
        ! 速度分布関数の更新
        call renew(f,fnext)
        !初期化
        call reset(p,u1,u2)
        !圧力の計算
        call pressure_cal(f,p)
        !流速の計算
        call velocity_cal(f,u1,u2)
        write(*,*) "step = ", n
    ENDDO
    !出力
    call output(u1,u2,p,f)
    do yi=0,ymax
        write(*,*) yi,u1(1,yi),u2(1,yi),p(2,yi)
    enddo
    call cpu_time(time2)
    write(*,*) time2-time1
end program main
