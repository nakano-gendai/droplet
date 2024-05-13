module globals
    real(8),parameter:: ds = 1.0d0 !格子間隔（lattice unit）
    integer,parameter:: xmax = 3 !ｘ方向格子数
    integer,parameter:: ymax = 49 !ｙ方向格子数
    integer,parameter:: step = 30000
    integer i, n, xi, yi, alpha,beta
    real(8) err_u1,err_up_u1,err_down_u1,err_u2,err_up_u2,err_down_u2, &
                err_p,err_up_p,err_down_p
    real(8) ut1,ut2,pt
    
    !粒子速度（整数）
    integer,parameter:: cx(9) = (/0, 1, 0, -1, 0, 1, -1, -1, 1/)
    integer,parameter:: cy(9) = (/0, 0, 1, 0, -1, 1, 1, -1, -1/)
    !粒子速度（実数）
    real(8):: cr(1:2, 1:9)
        
    !パラメータ
    real(8),parameter:: rho = 1.0d0 !液体の密度
    real(8),parameter:: nu = 0.01d0 !液体の動粘度
    real(8),parameter:: U(1:2) = (/0.01d0, 0.0d0/) !壁の速度
    real(8),parameter:: dp = 0.0d0
    real(8) c
    ! real(8),parameter:: A = 4.5d0*(1.0d0/6.0d0-nu/ds)
    real(8),parameter:: A = 0.0d0
    real(8),parameter:: E(9) = (/ 4.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/36.0d0, &
                        1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0 /)

    ! taylor
    real(8),parameter:: u0 = 0.01d0 !流速の初期値
    real(8),parameter:: p0 = rho / 3.0d0 !圧力の初期値
    real(8),parameter:: pi = acos(-1.0d0)
    real(8),parameter:: k1 = 2.0d0*pi/dble(xmax+1)
    real(8),parameter:: k2 = 2.0d0*pi/dble(ymax+1)

contains
    subroutine parv(cx,cy,cr)
        integer, intent(in):: cx(9),cy(9)
        real(8), intent(out):: cr(1:2,1:9)
            do i=1,9
                cr(1,i) = dble(cx(i))
                cr(2,i) = dble(cy(i))
            enddo
    end subroutine parv

    subroutine initial(p,u1,u2)
        real(8), intent(inout):: p(-1:xmax+1, -1:ymax+1)
        real(8), intent(inout):: u1(-1:xmax+1, -1:ymax+1),u2(-1:xmax+1, -1:ymax+1)
        do yi=0,ymax
            do xi=0,xmax
                ! p(xi,yi) = p0 -rho/4.0d0*(u0**2*sin(2.0d0*k1*dble(xi))+ &
                ! k1**2/k2**2*u0**2*cos(2.0d0*k2*dble(yi)))
                ! u1(xi,yi) = -u0*cos(k1*dble(xi))*sin(k2*dble(yi))
                ! u2(xi,yi) = (k1/k2)*u0*sin(k1*dble(xi))*cos(k2*dble(yi))
                p(xi,yi) = rho / 3.0d0
                u1(xi,yi) = 0.0d0
                u2(xi,yi) = 0.0d0
            enddo
        enddo
    end subroutine initial

    subroutine periodic(u1,u2,p)
        real(8), intent(inout):: p(-1:xmax+1, -1:ymax+1),u1(-1:xmax+1, -1:ymax+1),&
                            u2(-1:xmax+1, -1:ymax+1)
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
            enddo
        enddo
    end subroutine periodic

    subroutine periodic_f(feq)
        real(8),intent(inout):: feq(1:9,-1:xmax+1,-1:ymax+1)
        do yi=-1,ymax+1
            do xi=0,xmax
                do i=1,9
                    feq(i,xi,-1) = feq(i,xi,ymax)
                    feq(i,xi,ymax+1) = feq(i,xi,0)
                    feq(i,-1,yi) = feq(i,xmax,yi)
                    feq(i,xmax+1,yi) = feq(i,0,yi)
                enddo
            enddo
        enddo
    endsubroutine periodic_f 

    subroutine periodic_pressure(fnext)
        real(8),intent(inout) :: fnext(1:9, 0:xmax, 0:ymax)
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

    ! subroutine bounce_back(cr,feq)
    !     real(8),intent(in) :: cr(1:2,1:9)
    !     real(8),intent(inout) :: feq(1:9, -1:xmax+1, -1:ymax+1)
    !     ! !下の壁
    !     ! yi = 0
    !     ! do xi=0,xmax
    !     !     feq(3,xi-cx(3),yi-cy(3)) = feq(5,xi-cx(5),yi-cy(5)) - 6.0d0 * E(5) * ((-U(1))*cr(1,5)+(-U(2))*cr(2,5))
    !     !     feq(6,xi-cx(6),yi-cy(6)) = feq(8,xi-cx(8),yi-cy(8)) - 6.0d0 * E(8) * ((-U(1))*cr(1,8)+(-U(2))*cr(2,8))
    !     !     feq(7,xi-cx(7),yi-cy(7)) = feq(9,xi-cx(9),yi-cy(9)) - 6.0d0 * E(9) * ((-U(1))*cr(1,9)+(-U(2))*cr(2,9))
    !     ! enddo
    !     ! !上の壁
    !     ! yi = ymax
    !     ! do xi=0,xmax
    !     !     feq(5,xi-cx(5),yi-cy(5)) = feq(3,xi-cx(3),yi-cy(3)) - 6.0d0 * E(3) * ((U(1))*cr(1,3)+(U(2))*cr(2,3))
    !     !     feq(8,xi-cx(8),yi-cy(8)) = feq(6,xi-cx(6),yi-cy(6)) - 6.0d0 * E(6) * ((U(1))*cr(1,6)+(U(2))*cr(2,6))
    !     !     feq(9,xi-cx(9),yi-cy(9)) = feq(7,xi-cx(7),yi-cy(7)) - 6.0d0 * E(7) * ((U(1))*cr(1,7)+(U(2))*cr(2,7))
    !     ! enddo

    !     !下の壁
    !     yi = 0
    !     do xi=0,xmax
    !         feq(3,xi-cx(3),yi-cy(3)) = feq(5,xi-cx(3),yi-cy(3)) - 6.0d0 * E(5) * ((-U(1))*cr(1,5)+(-U(2))*cr(2,5))
    !         feq(6,xi-cx(6),yi-cy(6)) = feq(8,xi-cx(6),yi-cy(6)) - 6.0d0 * E(8) * ((-U(1))*cr(1,8)+(-U(2))*cr(2,8))
    !         feq(7,xi-cx(7),yi-cy(7)) = feq(9,xi-cx(7),yi-cy(7)) - 6.0d0 * E(9) * ((-U(1))*cr(1,9)+(-U(2))*cr(2,9))
    !     enddo
    !     !上の壁
    !     yi = ymax
    !     do xi=0,xmax
    !         feq(5,xi-cx(5),yi-cy(5)) = feq(3,xi-cx(5),yi-cy(5)) - 6.0d0 * E(3) * ((U(1))*cr(1,3)+(U(2))*cr(2,3))
    !         feq(8,xi-cx(8),yi-cy(8)) = feq(6,xi-cx(8),yi-cy(8)) - 6.0d0 * E(6) * ((U(1))*cr(1,6)+(U(2))*cr(2,6))
    !         feq(9,xi-cx(9),yi-cy(9)) = feq(7,xi-cx(9),yi-cy(9)) - 6.0d0 * E(7) * ((U(1))*cr(1,7)+(U(2))*cr(2,7))
    !     enddo
    ! end subroutine bounce_back


    subroutine bounce_back(cr,f)
        real(8),intent(in) :: cr(1:2,1:9)
        real(8),intent(inout) :: f(1:9, 0:xmax, 0:ymax)
        !下の壁
        yi = 0
        do xi=0,xmax
            f(3,xi,yi) = f(5,xi,yi) - 6.0d0 * E(5) * ((-U(1))*cr(1,5)+(-U(2))*cr(2,5))
            f(6,xi,yi) = f(8,xi,yi) - 6.0d0 * E(8) * ((-U(1))*cr(1,8)+(-U(2))*cr(2,8))
            f(7,xi,yi) = f(9,xi,yi) - 6.0d0 * E(9) * ((-U(1))*cr(1,9)+(-U(2))*cr(2,9))
        enddo
        !上の壁
        yi = ymax
        do xi=0,xmax
            f(5,xi,yi) = f(3,xi,yi) - 6.0d0 * E(3) * ((U(1))*cr(1,3)+(U(2))*cr(2,3))
            f(8,xi,yi) = f(6,xi,yi) - 6.0d0 * E(6) * ((U(1))*cr(1,6)+(U(2))*cr(2,6))
            f(9,xi,yi) = f(7,xi,yi) - 6.0d0 * E(7) * ((U(1))*cr(1,7)+(U(2))*cr(2,7))
        enddo
    end subroutine bounce_back


    subroutine physics(p,u1,u2,f)
        real(8),intent(inout):: p(-1:xmax+1,-1:ymax+1)
        real(8),intent(inout):: u1(-1:xmax+1,-1:ymax+1)
        real(8),intent(inout):: u2(-1:xmax+1,-1:ymax+1)
        ! real(8),intent(in):: f(1:9, -1:xmax+1, -1:ymax+1)
        real(8),intent(in):: f(1:9, 0:xmax, 0:ymax)
        
        do yi=0,ymax
            do xi=0,xmax
                p(xi,yi) = 0.0d0
                u1(xi,yi) = 0.0d0
                u2(xi,yi) = 0.0d0
                do i=1,9
                    ! p(xi,yi) = p(xi,yi) + f(i,xi-cx(i),yi-cy(i)) / 3.0d0
                    ! u1(xi,yi) = u1(xi,yi) + dble(cx(i))*f(i,xi-cx(i),yi-cy(i))
                    ! u2(xi,yi) = u2(xi,yi) + dble(cy(i))*f(i,xi-cx(i),yi-cy(i))
                    p(xi,yi) = p(xi,yi) + f(i,xi,yi) / 3.0d0
                    u1(xi,yi) = u1(xi,yi) + dble(cx(i))*f(i,xi,yi)
                    u2(xi,yi) = u2(xi,yi) + dble(cy(i))*f(i,xi,yi)
                enddo
            enddo
        enddo

    end subroutine physics

    subroutine grad_u_cal(grad_u,u1,u2)
        real(8),intent(inout):: grad_u(1:2,1:2,0:xmax,0:ymax)
        real(8),intent(in):: u1(-1:xmax+1,-1:ymax+1)
        real(8),intent(in):: u2(-1:xmax+1,-1:ymax+1)
        
        do yi=0,ymax
            do xi=0,xmax
                do beta=1,2
                    grad_u(1,beta,xi,yi)=0.0d0
                    grad_u(2,beta,xi,yi)=0.0d0
                    do i=1,9
                        grad_u(1,beta,xi,yi)=grad_u(1,beta,xi,yi)+cr(beta,i)*u1(xi+cx(i),yi+cy(i))
                        grad_u(2,beta,xi,yi)=grad_u(2,beta,xi,yi)+cr(beta,i)*u2(xi+cx(i),yi+cy(i))
                    enddo
                    grad_u(1,beta,xi,yi)=grad_u(1,beta,xi,yi)/(6.0d0*ds)
                    grad_u(2,beta,xi,yi)=grad_u(2,beta,xi,yi)/(6.0d0*ds)
                enddo
            enddo
        enddo
    end subroutine grad_u_cal

    subroutine output(u1,u2,p)
        real(8),intent(in) :: u1(-1:xmax+1, -1:ymax+1), u2(-1:xmax+1, -1:ymax+1)
        real(8),intent(in) :: p(-1:xmax+1, -1:ymax+1)
        open(10,file = "taylor_u_ori.d")
        do xi=0,xmax
            do yi=0,ymax
                write(10,*) xi,yi,u1(xi,yi),u2(xi,yi),p(xi,yi)
            enddo
        enddo
        close(10)

        open(11,file = "taylor_theory_ori.d")
        n = step
        do xi=0,xmax
            do yi=0,ymax
                ut1 = -u0*exp(-nu*dble(n)*(k1**2+k2**2))*cos(k1*dble(xi))*sin(k2*dble(yi))
                ut2 = k1/k2*u0*exp(-nu*n*(k1**2+k2**2))*sin(k1*dble(xi))*cos(k2*dble(yi))
                pt = p0 -rho/4.0d0*exp(-2.0d0*nu*n*(k1**2+k2**2))*(u0**2*sin(2.0d0*k1*dble(xi))+ &
                k1**2/k2**2*u0**2*cos(2.0d0*k2*dble(yi)))
                write(11,*) xi,yi,ut1,ut2,pt
            enddo
        enddo
        close(11)

        open(12,file = "taylor_err_ori.d")
        err_u1 = 0.0d0
        err_up_u1 = 0.0d0
        err_down_u1 = 0.0d0
        err_u2 = 0.0d0
        err_up_u2 = 0.0d0
        err_down_u2 = 0.0d0
        err_p = 0.0d0
        err_up_p = 0.0d0
        err_down_p = 0.0d0
        n = step
        do xi=0,xmax
            do yi=0,ymax
                ut1 = -u0*exp(-nu*dble(n)*(k1**2+k2**2))*cos(k1*dble(xi))*sin(k2*dble(yi))
                ut2 = k1/k2*u0*exp(-nu*dble(n)*(k1**2+k2**2))*sin(k1*dble(xi))*cos(k2*dble(yi))
                pt = p0 -rho/4.0d0*exp(-2.0d0*nu*dble(n)*(k1**2+k2**2))*(u0**2*sin(2.0d0*k1*dble(xi))+ &
                k1**2/k2**2*u0**2*cos(2.0d0*k2*dble(yi)))
                err_up_u1 = err_up_u1 + abs(ut1 - u1(xi,yi))
                err_down_u1 = err_down_u1 + abs(ut1)
                err_up_u2 = err_up_u2 + abs(ut2 - u2(xi,yi))
                err_down_u2 = err_down_u2 + abs(ut2)
                err_up_p = err_up_p + abs(pt - p(xi,yi))
                err_down_p = err_down_p + abs(pt)
            enddo
        enddo
        err_u1 = (err_up_u1 / err_down_u1)
        err_u2 = (err_up_u2 / err_down_u2)
        err_p = (err_up_p / err_down_p)
        write(12,*) err_u1,err_u2,err_p
        close(12)

    end subroutine output
    
end module globals
    
program main
use globals
    implicit none
    real(8) feq(1:9, -1:xmax+1, -1:ymax+1), f(1:9, 0:xmax, 0:ymax) !速度分布関数
    real(8) p(-1:xmax+1, -1:ymax+1)  !圧力
    real(8) u1(-1:xmax+1, -1:ymax+1), u2(-1:xmax+1, -1:ymax+1) !流速
    real(8) grad_u(1:2,1:2,0:xmax,0:ymax)
    real(8) time1,time2
    real(8) ftemp
    call cpu_time(time1)

    call parv(cx,cy,cr)
    !初期値の設定
    call initial(p,u1,u2)
!===========時間発展=====================
    DO n=1,step
        call periodic(u1,u2,p)
        call grad_u_cal(grad_u,u1,u2)
        !速度分布関数feqを計算
        do yi=0,ymax
            do xi=0,xmax
                do i=1,9
                    ftemp = 0.0d0
                    do beta=1,2
                        do alpha=1,2
                            ftemp = ftemp + (grad_u(beta,alpha,xi,yi)+grad_u(alpha,beta,xi,yi))*cr(alpha,i)*cr(beta,i)
                        enddo
                    enddo
                    feq(i,xi,yi)=E(i)*(3.0d0*p(xi,yi) &
                                +3.0d0*(dble(cx(i))*u1(xi,yi)+dble(cy(i))*u2(xi,yi)) &
                                +4.5d0*(dble(cx(i))*u1(xi,yi)+dble(cy(i))*u2(xi,yi))**2 &
                                -1.5d0*(u1(xi,yi)**2+u2(xi,yi)**2)) + A*E(i)*ds*ftemp
                enddo
            enddo
        enddo
        
        !周期境界
        call periodic_f(feq)
        do yi=0,ymax
            do xi=0,xmax
                do i=1,9
                    f(i,xi,yi) = feq(i,xi-cx(i),yi-cy(i))
                enddo
            enddo
        enddo
        !壁（yi=0,ymax）
        call bounce_back(cr,f)
        !phisics
        call physics(p,u1,u2,f)
        write(*,*) "step = ", n
    ENDDO

    ! call output(u1,u2,p)

    open(10,file="couette.d")
    do yi=0,ymax
        write(10,*) yi,u1(1,yi),u2(1,yi),p(2,yi)
    enddo
    close(10)
    call cpu_time(time2)
    write(*,*) time2-time1
end program main
