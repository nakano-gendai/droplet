module globals
    real(8),parameter:: ds = 1.0d0
    integer,parameter:: xmax = 50
    integer,parameter:: ymax = 50
    integer,parameter:: step = 10500
    integer i, n, xi, yi
    real(8) tmp, norm, veff, reff
    character :: filename*20
    !粒子速度（整数）
    integer,parameter:: cx(9) = (/0, 1, 0, -1, 0, 1, -1, -1, 1/)
    integer,parameter:: cy(9) = (/0, 0, 1, 0, -1, 1, 1, -1, -1/)
    real(8):: cr(1:2, 1:9)  !粒子速度（実数）
    !パラメータ
    real(8),parameter:: E(9) = (/ 4.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/36.0d0, &
                            1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0 /)
    real(8),parameter:: pi = acos(-1.0d0)
    real(8),parameter:: phi1 = 1.0d0
    real(8),parameter:: phi2 = -1.0d0
    real(8),parameter:: phi0 = (phi1-phi2) / 2.0d0
    real(8),parameter:: tauf = 0.7d0
    real(8),parameter:: taug = 0.7d0
    real(8),parameter:: sigma = 0.0001d0
    real(8),parameter:: width = 5.0d0
    real(8),parameter:: beta = 3.0d0/4.0d0*sigma/width*phi0**4
    real(8),parameter:: kappa = 3.0d0/8.0d0*sigma*width/phi0**2
    real(8),parameter:: gamma = 10.0d0
    real(8),parameter:: H = 2.0d0*taug - 1.0d0
    real(8),parameter:: xc = dble(xmax)/2.0d0
    real(8),parameter:: yc = dble(ymax)/2.0d0
    real(8),parameter:: R = 12.5d0
    real(8) phi_max, phi_min
    
contains
    subroutine par(cx,cy,cr)
        integer,intent(in):: cx(9), cy(9)
        real(8),intent(out):: cr(1:2,1:9)
        do i=1,9
            cr(1,i) = dble(cx(i))
            cr(2,i) = dble(cy(i))
        enddo
    end subroutine par

    subroutine periodic(fun)
        real(8),intent(inout):: fun(-1:xmax+1,-1:ymax+1)
        do yi=-1,ymax+1
            do xi=0,xmax
                fun(xi,-1) = fun(xi,ymax)
                fun(xi,ymax+1) = fun(xi,0)
                fun(-1,yi) = fun(xmax,yi)
                fun(xmax+1,yi) = fun(0,yi)
            enddo
        enddo
    end subroutine periodic

    subroutine periodic_fg(fun_fg)
        real(8),intent(inout):: fun_fg(1:9,-1:xmax+1,-1:ymax+1)
        do yi=-1,ymax+1
            do xi=0,xmax
                do i=1,9
                    fun_fg(i,xi,-1) = fun_fg(i,xi,ymax)
                    fun_fg(i,xi,ymax+1) = fun_fg(i,xi,0)
                    fun_fg(i,-1,yi) = fun_fg(i,xmax,yi)
                    fun_fg(i,xmax+1,yi) = fun_fg(i,0,yi)
                enddo
            enddo
        enddo
    end subroutine periodic_fg

    function feq_1_cal(in,rho,u1,u2,phi,che) result(feq_1_result)
        integer,intent(in):: in
        real(8),intent(in):: rho,u1,u2,phi,che
        real(8) feq_1_result
        feq_1_result = E(in)*rho*(1.0d0-1.5d0*(u1**2+u2**2))-15.0d0/4.0d0*E(in)*phi*che
    end function feq_1_cal

    function feq_cal(in,rho,u1,u2,phi,che) result(feq_result)
        integer,intent(in):: in
        real(8),intent(in):: rho,u1,u2,phi,che
        real(8) feq_result
        feq_result = E(in)*rho*(1.0d0+3.0d0*(dble(cx(in))*u1+dble(cy(in))*u2) &
                -1.5d0*(u1**2+u2**2)+4.5d0*(dble(cx(in))*u1+dble(cy(in))*u2)**2)+E(in)*3.0d0*phi*che
    end function feq_cal

    function geq_1_cal(in,che,phi) result(geq_1_result)
        integer,intent(in):: in
        real(8),intent(in):: che,phi
        real(8) geq_1_result
        geq_1_result = phi - (1.0d0-E(in))*gamma*che
    end function geq_1_cal

    function geq_cal(in,che,phi,u1,u2) result(geq_result)
        integer,intent(in):: in
        real(8),intent(in):: che,phi,u1,u2
        real(8) geq_result
        geq_result = E(in)*(gamma*che+3.0d0*phi*(dble(cx(in))*u1+dble(cy(in))*u2))
    end function geq_cal
end module globals

program main
use globals
    implicit none
    real(8) phi(-1:xmax+1,-1:ymax+1)
    real(8) u1(-1:xmax+1,-1:ymax+1), u2(-1:xmax+1,-1:ymax+1)
    real(8) u1n(-1:xmax+1,-1:ymax+1), u2n(-1:xmax+1,-1:ymax+1)
    real(8) che(-1:xmax+1,-1:ymax+1)
    real(8) rho(-1:xmax+1,-1:ymax+1)
    real(8) g(1:9,-1:xmax+1,-1:ymax+1), f(1:9,-1:xmax+1,-1:ymax+1)
    real(8) gnext(1:9,-1:xmax+1,-1:ymax+1), fnext(1:9,-1:xmax+1,-1:ymax+1)
    real(8) fi(1:9,-1:xmax+1,-1:ymax+1)
    real(8) fx(0:xmax,0:ymax), fy(0:xmax,0:ymax)

!===========初期条件=======================================================
    do yi=0,ymax
        do xi=0,xmax
            u1(xi,yi) = 0.0d0
            u2(xi,yi) = 0.0d0
            u1n(xi,yi) = 0.0d0
            u2n(xi,yi) = 0.0d0
            rho(xi,yi) = 1.0d0
            phi(xi,yi) = phi0*tanh(2.0d0*(R-sqrt((dble(xi)-xc)**2+(dble(yi)-yc)**2))/width)
        enddo
    enddo
    call periodic(u1)
    call periodic(u2)
    call periodic(rho)
    call periodic(phi)
    do yi=0,ymax
        do xi=0,xmax
            tmp = (phi(xi+1,yi)+phi(xi-1,yi)+phi(xi,yi+1)+phi(xi,yi-1)-4.0d0*phi(xi,yi))/ds**2
            che(xi,yi) = 4.0d0*beta*(phi(xi,yi)**3-phi(xi,yi)*phi0**2)-kappa*tmp
        enddo
    enddo
    call periodic(che)
    do yi=0,ymax
        do xi=0,xmax
            i = 1
            f(i,xi,yi) = feq_1_cal(i,rho(xi,yi),u1(xi,yi),u2(xi,yi),phi(xi,yi),che(xi,yi))
            g(i,xi,yi) = geq_1_cal(i,che(xi,yi),phi(xi,yi))
            do i=2,9
                f(i,xi,yi) = feq_cal(i,rho(xi,yi),u1(xi,yi),u2(xi,yi),phi(xi,yi),che(xi,yi))
                g(i,xi,yi) = geq_cal(i,che(xi,yi),phi(xi,yi),u1(xi,yi),u2(xi,yi))
            enddo
        enddo
    enddo
    call periodic_fg(f)
    call periodic_fg(g)
!===========================================================================
!========時間発展スタート====================================================
DO n=1,step
    call periodic(u1)
    call periodic(u2)
    call periodic(rho)
    call periodic(phi)
    call periodic(che)
    call periodic_fg(f)
    call periodic_fg(g)

    do yi=0,ymax
        do xi=0,xmax
            u1n(xi,yi) = u1(xi,yi)
            u2n(xi,yi) = u2(xi,yi)
        enddo
    enddo
    !=======外力項の計算================================================
    do yi=0,ymax
        do xi=0,xmax
            fx(xi,yi) = 0.0d0
            fy(xi,yi) = 0.0d0
            fx(xi,yi) = che(xi,yi)*(phi(xi+1,yi)-phi(xi-1,yi))*0.5d0
            fy(xi,yi) = che(xi,yi)*(phi(xi,yi+1)-phi(xi,yi-1))*0.5d0

            do i=1,9
                fi(i,xi,yi) = E(i)*(1.0d0-0.5d0/tauf) &
                *(3.0d0*(dble(cx(i))*fx(xi,yi)+dble(cy(i))*fy(xi,yi)) &
                -3.0d0*(u1(xi,yi)*fx(xi,yi)+u2(xi,yi)*fy(xi,yi)) &
                +9.0d0*(u1(xi,yi)*fx(xi,yi)*dble(cx(i))*dble(cx(i)) &
                        +u1(xi,yi)*fy(xi,yi)*dble(cx(i))*dble(cy(i)) &
                        +u2(xi,yi)*fx(xi,yi)*dble(cx(i))*dble(cy(i)) &
                        +u2(xi,yi)*fy(xi,yi)*dble(cy(i))*dble(cy(i))))
            enddo
        enddo
    enddo
    call periodic_fg(fi)
    !========SRT===========================================================
    do yi=0,ymax
        do xi=0,xmax
            i = 1
            fnext(i,xi,yi) = f(i,xi-cx(i),yi-cy(i))-(f(i,xi-cx(i),yi-cy(i)) &
                -feq_1_cal(i,rho(xi-cx(i),yi-cy(i)),u1(xi-cx(i),yi-cy(i)) &
                ,u2(xi-cx(i),yi-cy(i)),phi(xi-cx(i),yi-cy(i)),che(xi-cx(i),yi-cy(i))))/tauf &
                +fi(i,xi-cx(i),yi-cy(i))
            
            gnext(i,xi,yi) = g(i,xi-cx(i),yi-cy(i))-(g(i,xi-cx(i),yi-cy(i)) &
                -geq_1_cal(i,che(xi-cx(i),yi-cy(i)),phi(xi-cx(i),yi-cy(i))))/taug
            do i=2,9
                fnext(i,xi,yi) = f(i,xi-cx(i),yi-cy(i))-(f(i,xi-cx(i),yi-cy(i)) &
                -feq_cal(i,rho(xi-cx(i),yi-cy(i)),u1(xi-cx(i),yi-cy(i)) &
                ,u2(xi-cx(i),yi-cy(i)),phi(xi-cx(i),yi-cy(i)),che(xi-cx(i),yi-cy(i))))/tauf &
                +fi(i,xi-cx(i),yi-cy(i))

                gnext(i,xi,yi) = g(i,xi-cx(i),yi-cy(i))-(g(i,xi-cx(i),yi-cy(i)) &
                -geq_cal(i,che(xi-cx(i),yi-cy(i)),phi(xi-cx(i),yi-cy(i)),u1(xi-cx(i),yi-cy(i)),u2(xi-cx(i),yi-cy(i))))/taug
            enddo
        enddo
    enddo
    do yi=0,ymax
        do xi=0,xmax
            do i=1,9
                f(i,xi,yi) = fnext(i,xi,yi)
                g(i,xi,yi) = gnext(i,xi,yi)
            enddo
        enddo
    enddo
    !========速度・密度==========================================================
    do yi=0,ymax
        do xi=0,xmax
            rho(xi,yi) = 0.0d0
            u1(xi,yi) = 0.0d0
            u2(xi,yi) = 0.0d0
            phi(xi,yi) = 0.0d0
            do i=1,9
                rho(xi,yi) = rho(xi,yi)+f(i,xi,yi)
                u1(xi,yi) = u1(xi,yi)+f(i,xi,yi)*dble(cx(i))
                u2(xi,yi) = u2(xi,yi)+f(i,xi,yi)*dble(cy(i))
                phi(xi,yi) = phi(xi,yi)+g(i,xi,yi)
            enddo
            u1(xi,yi) = (u1(xi,yi)+0.5d0*fx(xi,yi))/rho(xi,yi)
            u2(xi,yi) = (u2(xi,yi)+0.5d0*fy(xi,yi))/rho(xi,yi)
        enddo
    enddo
    !=========化学ポテンシャル====================================================
    call periodic(phi)
    do yi=0,ymax
        do xi=0,xmax
            tmp = (phi(xi+1,yi)+phi(xi-1,yi)+phi(xi,yi+1)+phi(xi,yi-1)-4.0d0*phi(xi,yi))/ds**2
            che(xi,yi) = 4.0d0*beta*(phi(xi,yi)**3-phi(xi,yi)*phi0**2)-kappa*tmp
        enddo
    enddo
    !========前の時間との比較====================================================
    norm = 0.0d0
    do yi=0,ymax
        do xi=0,xmax
            tmp = sqrt((u1(xi,yi)-u1n(xi,yi))**2+(u2(xi,yi)-u2n(xi,yi))**2)
            if(tmp > norm) then
                norm = tmp
            endif
        enddo
    enddo
    !========出力==============================================================
    ! 圧力差を求める phi1 = 2.211 phi2 = 4.895
    veff = 0.0d0
    
    do yi=0,ymax
        do xi=0,xmax
            if(phi(xi,yi) >= 0.5d0*(phi1+phi2)) then
                veff = veff + ds**3
            end if
        enddo
    enddo
    
    reff = sqrt(veff/pi)
    write(*,*) "step, norm, lap, dp =", n,norm,sigma/R,(rho(xmax/2,ymax/2)-rho(0,0))/3.0d0,veff,reff
ENDDO
    open(20,file="lbm_u.txt")
    do xi=0,xmax
        do yi=0,ymax
            write(20,*) xi,yi,1000000.0d0*u1(xi,yi),1000000.0d0*u2(xi,yi) &
            ,sqrt(u1(xi,yi)**2+u2(xi,yi)**2)
        enddo
        write(20,*)
    enddo
    close(20)
end program main