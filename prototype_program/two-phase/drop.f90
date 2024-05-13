module globals
    real(8),parameter:: ds = 1.0d0 !格子間隔（lattice unit）
    integer,parameter:: xmax = 39 !ｘ方向格子数
    integer,parameter:: ymax = 39 !ｙ方向格子数
    integer,parameter:: zmax = 39 !ｚ方向格子数
    integer,parameter:: step = 1 !計算時間step
    integer,parameter:: start = 10000 !壁を動かす時間step
    integer i, n, xi, yi, zi, alpha, beta
    character :: filename*20
    character(*),parameter :: datadir = "./data_p1/"

    !粒子速度（整数）
    integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
    integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
    integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)
    real(8):: cr(1:3, 1:15)  !粒子速度（実数）
    real(8) krone(1:3,1:3) !クロネッカーのデルタ

    
    !パラメータ
    real(8),parameter:: E(15) = (/ 2.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, &
                                1.0d0/9.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, &
                                1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0 /)
    real(8),parameter:: H(15) = (/ 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
                                0.0d0, 0.0d0, 0.0d0 /)
    real(8),parameter:: F(15) = (/ -7.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, &
                                1.d0/24.d0, 1.d0/24.d0, 1.d0/24.d0, 1.d0/24.d0, &
                                1.d0/24.d0, 1.d0/24.d0, 1.d0/24.d0, 1.d0/24.d0 /)
    real(8),parameter:: pi = acos(-1.0d0)

    !支配パラメータ
    real(8),parameter:: D = 20.0d0 * ds !設置する液滴のサイズ
    real(8),parameter:: nu1 = 0.1d0 !液体の動粘度
    real(8),parameter:: nu2 = 0.1d0 !液滴の動粘度
    real(8),parameter:: uw = 9.375d-3
    real(8),parameter:: kappaf = 0.01d0*ds**2
    real(8),parameter:: kappag = 1.0d-6
    real(8),parameter:: H1 = dble(zmax)
    real(8),parameter:: Re = 2.0d0*uw*D**2/(4.0d0*H1*nu1)
    real(8),parameter:: Ca = 2.0d0*uw*nu1*D/(2.0d0*H1*2.5d-3)
    real(8),parameter:: phi1 = 2.211
    real(8),parameter:: phi2 = 4.895
    real(8),parameter:: a = 9.0d0/49.0d0
    real(8),parameter:: b = 2.0d0/21.0d0
    real(8),parameter:: T = 0.55d0
    real(8),parameter:: C = 0.0d0
    real(8),parameter:: M = (0.5d0-C/3.0d0)*ds

    !初期液滴中心
    real(8), parameter:: xc = 0.5d0*dble(xmax)	
    real(8), parameter:: yc = 0.5d0*dble(ymax)
    real(8), parameter:: zc = 0.5d0*dble(zmax)
    
    !その他パラメータ
    real(8) phi_min, phi_max, p_in, p_out, veff, reff, sigma, seff, distance
    real(8) dp_cal, dp_th
    real(8) L1, B1, Dd
    real(8) gtemp
    real(8) ftemp

contains
    subroutine par(cx,cy,cz,cr,krone,U)
        integer, intent(in):: cx(15),cy(15),cz(15)
        real(8), intent(out):: cr(1:3,1:15), krone(1:3,1:3), U(1:3)
        do i=1,15
            cr(1,i) = dble(cx(i))
            cr(2,i) = dble(cy(i))
            cr(3,i) = dble(cz(i))
        enddo
        krone(1,1) = 1.d0; krone(1,2) = 0.d0; krone(1,3) = 0.d0
        krone(2,1) = 0.d0; krone(2,2) = 1.d0; krone(2,3) = 0.d0
        krone(3,1) = 0.d0; krone(3,2) = 0.d0; krone(3,3) = 1.d0
        U(1) = 0.0d0
        U(2) = 0.0d0
        U(3) = 0.0d0
    end subroutine par

    subroutine ini(phi,u1,u2,u3,p,nu,Anu)
        real(8),intent(inout):: phi(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: u1(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: u2(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: u3(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: p(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: nu(0:xmax,0:ymax,0:zmax)
        real(8),intent(inout):: Anu(0:xmax,0:ymax,0:zmax)
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    phi(xi,yi,zi) = phi1
                    u1(xi,yi,zi) = 0.0d0
                    u2(xi,yi,zi) = 0.0d0
                    u3(xi,yi,zi) = 0.0d0
                    p(xi,yi,zi) = 0.0d0

                    if((dble(xi)*ds-xc)**2+(dble(yi)*ds-yc)**2+(dble(zi)*ds-zc)**2 <= (0.5d0*D)**2) then
                        phi(xi,yi,zi) = phi2
                    endif

                    nu(xi,yi,zi) = 0.5d0*(nu2-nu1)*(dsin(pi*(phi(xi,yi,zi)-0.5d0*(phi2+phi1))/(phi2-phi1)) + 1.d0) + nu1
                    Anu(xi,yi,zi) = 4.5d0*(1.0d0/6.0d0-nu(xi,yi,zi)/ds)
                enddo
            enddo
        enddo
    end subroutine ini

    subroutine periodic(fun)
        real(8),intent(inout):: fun(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        do zi=0,zmax
            do yi=-1,ymax+1
                do xi=0,xmax
                    fun(xi,-1,zi) = fun(xi,ymax,zi)
                    fun(xi,ymax+1,zi) = fun(xi,0,zi)
                    fun(-1,yi,zi) = fun(xmax,yi,zi)
                    fun(xmax+1,yi,zi) = fun(0,yi,zi)
                enddo
            enddo
        enddo
        do yi=-1,ymax+1
            do xi=-1,xmax+1
                fun(xi,yi,-1) = fun(xi,yi,zmax)
                fun(xi,yi,zmax+1) = fun(xi,yi,0)
            enddo
        enddo
    end subroutine periodic

    subroutine periodic_fg(fun_fg)
        real(8),intent(inout):: fun_fg(1:15,-1:xmax+1,-1:ymax+1,-1:zmax+1)
        do zi=0,zmax
            do yi=-1,ymax+1
                do xi=0,xmax
                    do i=1,15
                        fun_fg(i,xi,-1,zi) = fun_fg(i,xi,ymax,zi)
                        fun_fg(i,xi,ymax+1,zi) = fun_fg(i,xi,0,zi)
                        fun_fg(i,-1,yi,zi) = fun_fg(i,xmax,yi,zi)
                        fun_fg(i,xmax+1,yi,zi) = fun_fg(i,0,yi,zi)
                    enddo
                enddo
            enddo
        enddo
        do yi=-1,ymax+1
            do xi=-1,xmax+1
                do i=1,15
                    fun_fg(i,xi,yi,-1) = fun_fg(i,xi,yi,zmax)
                    fun_fg(i,xi,yi,zmax+1) = fun_fg(i,xi,yi,0)
                enddo
            enddo
        enddo
    end subroutine periodic_fg

    subroutine feq_cal(gphi,phi,p0,lap_phi,grad_phi,u1,u2,u3,div_pcap,feq)
        real(8),intent(in):: gphi(1:3,1:3,0:xmax,0:ymax,0:zmax)
        real(8),intent(in):: phi(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(in):: p0(0:xmax,0:ymax,0:zmax)
        real(8),intent(in):: lap_phi(0:xmax,0:ymax,0:zmax)
        real(8),intent(in):: grad_phi(1:3,0:xmax,0:ymax,0:zmax)
        real(8),intent(in):: u1(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(in):: u2(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(in):: u3(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(in):: div_pcap(1:3,0:xmax,0:ymax,0:zmax)
        real(8),intent(inout):: feq(1:15,-1:xmax+1,-1:ymax+1,-1:zmax+1)
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    do i=1,15
                        gtemp = 0.0d0
                        do beta=1,3
                            do alpha=1,3
                                gtemp = gtemp + gphi(alpha,beta,xi,yi,zi)*cr(alpha,i)*cr(beta,i)
                            enddo
                        enddo
                        feq(i,xi,yi,zi) = H(i)*phi(xi,yi,zi) &
                                        + F(i)*(p0(xi,yi,zi)-kappaf*phi(xi,yi,zi)*lap_phi(xi,yi,zi) &
                                        -kappaf/6.0d0*(grad_phi(1,xi,yi,zi)**2+grad_phi(2,xi,yi,zi)**2+ &
                                                grad_phi(3,xi,yi,zi)**2)) &
                                        + 3.0d0*E(i)*phi(xi,yi,zi)*(cr(1,i)*u1(xi,yi,zi)+ &
                                                cr(2,i)*u2(xi,yi,zi)+cr(3,i)*u3(xi,yi,zi)) &
                                        + E(i)*kappaf*gtemp &
                                        + E(i)*C*(div_pcap(1,xi,yi,zi)*cr(1,i)+ &
                                                div_pcap(2,xi,yi,zi)*cr(2,i)+div_pcap(3,xi,yi,zi)*cr(3,i))*ds
                    enddo
                enddo
            enddo
        enddo
    end subroutine feq_cal

    subroutine geq_cal(gphi,grad_u,p,u1,u2,u3,Anu,geq)
        real(8),intent(in):: gphi(1:3,1:3,0:xmax,0:ymax,0:zmax)
        real(8),intent(in):: grad_u(1:3,1:3,0:xmax,0:ymax,0:zmax)
        real(8),intent(in):: p(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(in):: u1(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(in):: u2(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(in):: u3(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(in):: Anu(0:xmax,0:ymax,0:zmax)
        real(8),intent(inout):: geq(1:15,-1:xmax+1,-1:ymax+1,-1:zmax+1)
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    do i=1,15
                        gtemp = 0.0d0
                        ftemp = 0.0d0
                        do beta=1,3
                            do alpha=1,3
                                gtemp = gtemp + gphi(alpha,beta,xi,yi,zi)*cr(alpha,i)*cr(beta,i)
    
                                ftemp = ftemp + (grad_u(beta,alpha,xi,yi,zi)+grad_u(alpha,beta,xi,yi,zi))*cr(alpha,i)*cr(beta,i)
                            enddo
                        enddo
                        geq(i,xi,yi,zi) = E(i)*(3.0d0*p(xi,yi,zi) + 3.0d0*(u1(xi,yi,zi)*dble(cx(i)) &
                        + u2(xi,yi,zi)*dble(cy(i)) + u3(xi,yi,zi)*dble(cz(i))) &
                        + 4.5d0*(u1(xi,yi,zi)*dble(cx(i)) + u2(xi,yi,zi)*dble(cy(i)) + u3(xi,yi,zi)*dble(cz(i)))**2 &
                        -1.5d0*(u1(xi,yi,zi)**2 + u2(xi,yi,zi)**2 + u3(xi,yi,zi)**2) &
                        +Anu(xi,yi,zi)*ds*ftemp) &
                        + E(i)*kappag*gtemp
                    enddo
                enddo
            enddo
        enddo
    end subroutine geq_cal

    subroutine bounce_back(fun,U)
        real(8),intent(inout):: fun(1:15,-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: U(1:3)
        !安定したら壁を動かす
        if((n >= start)) then
            U(1) = uw
            U(2) = 0.0d0
            U(3) = 0.0d0
        endif
        zi = 0
        do yi=0,ymax
            do xi=0,xmax
                fun(4,xi,yi,zi) = fun(7,xi,yi,zi) - 6.0d0 * E(7) * ((-U(1))*cr(1,7)+(-U(2))*cr(2,7)+(-U(3))*cr(3,7))
                fun(8,xi,yi,zi) = fun(12,xi,yi,zi) - 6.0d0 * E(12) * ((-U(1))*cr(1,12)+(-U(2))*cr(2,12)+(-U(3))*cr(3,12))
                fun(9,xi,yi,zi) = fun(13,xi,yi,zi) - 6.0d0 * E(13) * ((-U(1))*cr(1,13)+(-U(2))*cr(2,13)+(-U(3))*cr(3,13))
                fun(10,xi,yi,zi) = fun(14,xi,yi,zi) - 6.0d0 * E(14) * ((-U(1))*cr(1,14)+(-U(2))*cr(2,14)+(-U(3))*cr(3,14))
                fun(15,xi,yi,zi) = fun(11,xi,yi,zi) - 6.0d0 * E(11) * ((-U(1))*cr(1,11)+(-U(2))*cr(2,11)+(-U(3))*cr(3,11))
            enddo
        enddo
        zi = zmax
        do yi=0,ymax
            do xi=0,xmax
                fun(7,xi,yi,zi) = fun(4,xi,yi,zi) - 6.0d0 * E(4) * ((U(1))*cr(1,4)+(U(2))*cr(2,4)+(U(3))*cr(3,4))
                fun(11,xi,yi,zi) = fun(15,xi,yi,zi) - 6.0d0 * E(15) * ((U(1))*cr(1,15)+(U(2))*cr(2,15)+(U(3))*cr(3,15))
                fun(12,xi,yi,zi) = fun(8,xi,yi,zi) - 6.0d0 * E(8) * ((U(1))*cr(1,8)+(U(2))*cr(2,8)+(U(3))*cr(3,8))
                fun(13,xi,yi,zi) = fun(9,xi,yi,zi) - 6.0d0 * E(9) * ((U(1))*cr(1,9)+(U(2))*cr(2,9)+(U(3))*cr(3,9))
                fun(14,xi,yi,zi) = fun(10,xi,yi,zi) - 6.0d0 * E(10) * ((U(1))*cr(1,10)+(U(2))*cr(2,10)+(U(3))*cr(3,10))
            enddo
        enddo
    end subroutine bounce_back

    subroutine physics(phi,p,u1,u2,u3,feq,geq,nu,Anu)
        real(8),intent(inout):: phi(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: p(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: u1(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: u2(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: u3(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(in):: feq(1:15,-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(in):: geq(1:15,-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: nu(0:xmax,0:ymax,0:zmax)
        real(8),intent(inout):: Anu(0:xmax,0:ymax,0:zmax)
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    phi(xi,yi,zi) = 0.0d0
                    p(xi,yi,zi) = 0.0d0
                    u1(xi,yi,zi) = 0.0d0
                    u2(xi,yi,zi) = 0.0d0
                    u3(xi,yi,zi) = 0.0d0
                    do i=1,15
                        phi(xi,yi,zi) = phi(xi,yi,zi) + feq(i,xi-cx(i),yi-cy(i),zi-cz(i))
                        p(xi,yi,zi) = p(xi,yi,zi) + geq(i,xi-cx(i),yi-cy(i),zi-cz(i)) / 3.0d0
                        u1(xi,yi,zi) = u1(xi,yi,zi) + dble(cx(i))*geq(i,xi-cx(i),yi-cy(i),zi-cz(i))
                        u2(xi,yi,zi) = u2(xi,yi,zi) + dble(cy(i))*geq(i,xi-cx(i),yi-cy(i),zi-cz(i))
                        u3(xi,yi,zi) = u3(xi,yi,zi) + dble(cz(i))*geq(i,xi-cx(i),yi-cy(i),zi-cz(i))
                    enddo
                    nu(xi,yi,zi) = 0.5d0*(nu2-nu1)*(dsin(pi*(phi(xi,yi,zi)-0.5d0*(phi2+phi1))/(phi2-phi1)) + 1.d0) + nu1
                    Anu(xi,yi,zi) = 4.5d0*(1.0d0/6.0d0-nu(xi,yi,zi)/ds)
                enddo
            enddo
        enddo
    end subroutine physics

    subroutine lap_cal(lap_f,fun)
        real(8),intent(inout):: lap_f(0:xmax,0:ymax,0:zmax)
        real(8),intent(in):: fun(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    lap_f(xi,yi,zi) = -14.0d0*fun(xi,yi,zi)
                    do i=2,15
                        lap_f(xi,yi,zi) = lap_f(xi,yi,zi) + fun(xi+cx(i),yi+cy(i),zi+cz(i))
                    enddo
                    lap_f(xi,yi,zi) = lap_f(xi,yi,zi) / (5.0d0*ds**2)
                enddo
            enddo
        enddo
    end subroutine lap_cal

    subroutine grad_cal(grad_f,fun)
        real(8),intent(inout):: grad_f(1:3,0:xmax,0:ymax,0:zmax)
        real(8),intent(in):: fun(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    do alpha=1,3
                        grad_f(alpha,xi,yi,zi) = 0.0d0
                        do i = 2,15
                            grad_f(alpha,xi,yi,zi) = grad_f(alpha,xi,yi,zi) + cr(alpha,i)*fun(xi+cx(i),yi+cy(i),zi+cz(i))
                        enddo
                        grad_f(alpha,xi,yi,zi) = grad_f(alpha,xi,yi,zi) / (10.0d0*ds)
                    enddo
                enddo
            enddo
        enddo
    end subroutine grad_cal

    subroutine grad_u_cal(grad_u,u1,u2,u3)
        real(8),intent(inout):: grad_u(1:3,1:3,0:xmax,0:ymax,0:zmax)
        real(8),intent(in):: u1(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(in):: u2(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(in):: u3(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    do beta=1,3
                        grad_u(1,beta,xi,yi,zi)=0.0d0
                        grad_u(2,beta,xi,yi,zi)=0.0d0
                        grad_u(3,beta,xi,yi,zi)=0.0d0
                        do i=2,15
                            grad_u(1,beta,xi,yi,zi)=grad_u(1,beta,xi,yi,zi)+cr(beta,i)*u1(xi+cx(i),yi+cy(i),zi+cz(i))
                            grad_u(2,beta,xi,yi,zi)=grad_u(2,beta,xi,yi,zi)+cr(beta,i)*u2(xi+cx(i),yi+cy(i),zi+cz(i))
                            grad_u(3,beta,xi,yi,zi)=grad_u(3,beta,xi,yi,zi)+cr(beta,i)*u3(xi+cx(i),yi+cy(i),zi+cz(i))
                        enddo
                        grad_u(1,beta,xi,yi,zi)=grad_u(1,beta,xi,yi,zi)/(10.0d0*ds)
                        grad_u(2,beta,xi,yi,zi)=grad_u(2,beta,xi,yi,zi)/(10.0d0*ds)
                        grad_u(3,beta,xi,yi,zi)=grad_u(3,beta,xi,yi,zi)/(10.0d0*ds)
                    enddo
                enddo
            enddo
        enddo
    end subroutine grad_u_cal

    subroutine p0_cal(p0,phi)
        real(8),intent(inout):: p0(0:xmax,0:ymax,0:zmax)
        real(8),intent(in):: phi(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    p0(xi,yi,zi) = phi(xi,yi,zi)*T/(1.0d0-b*phi(xi,yi,zi)) - a*phi(xi,yi,zi)**2
                enddo
            enddo
        enddo
    end subroutine p0_cal

    subroutine gphi_cal(gphi,grad_phi)
        real(8),intent(inout):: gphi(1:3,1:3,0:xmax,0:ymax,0:zmax)
        real(8),intent(in):: grad_phi(1:3,0:xmax,0:ymax,0:zmax)
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    do beta=1,3
                        do alpha=1,3
                            gphi(alpha,beta,xi,yi,zi) = 4.5d0*grad_phi(alpha,xi,yi,zi)*grad_phi(beta,xi,yi,zi)&
                            -1.5d0*(grad_phi(1,xi,yi,zi)**2+grad_phi(2,xi,yi,zi)**2+grad_phi(3,xi,yi,zi)**2)&
                            *krone(alpha,beta)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine gphi_cal

    subroutine pcap_cal(pcap,p0,phi,lap_phi,grad_phi)
        real(8),intent(inout):: pcap(1:3,1:3,0:xmax,0:ymax,0:zmax)
        real(8),intent(in):: p0(0:xmax,0:ymax,0:zmax)
        real(8),intent(in):: phi(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(in):: lap_phi(0:xmax,0:ymax,0:zmax)
        real(8),intent(in):: grad_phi(1:3,0:xmax,0:ymax,0:zmax)
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    do beta=1,3
                        do alpha=1,3
                            pcap(alpha,beta,xi,yi,zi) = (p0(xi,yi,zi)-kappaf*phi(xi,yi,zi)*lap_phi(xi,yi,zi)&
                                    -0.5d0*kappaf*(grad_phi(1,xi,yi,zi)**2+grad_phi(2,xi,yi,zi)**2&
                                    +grad_phi(3,xi,yi,zi)**2))*krone(alpha,beta)&
                                    +kappaf*grad_phi(alpha,xi,yi,zi)*grad_phi(beta,xi,yi,zi)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine pcap_cal

    subroutine div_pcap_cal(div_pcap,pcap)
        real(8),intent(inout):: div_pcap(1:3,0:xmax,0:ymax,0:zmax)
        real(8),intent(in):: pcap(1:3,1:3,0:xmax,0:ymax,0:zmax)
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    do alpha=1,3
                        div_pcap(alpha,xi,yi,zi) = 0.0d0
                        do i=2,15
                            do beta=1,3
                                div_pcap(alpha,xi,yi,zi) = div_pcap(alpha,xi,yi,zi) + cr(beta,i)* &
                                                            pcap(alpha,beta,xi+cx(i),yi+cy(i),zi+cz(i))
                            enddo
                        enddo
                        div_pcap(alpha,xi,yi,zi) = div_pcap(alpha,xi,yi,zi) / (10.0d0*ds)
                    enddo
                enddo
            enddo
        enddo
    end subroutine div_pcap_cal

    subroutine laplace(phi,p)
        real(8),intent(in):: phi(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(in):: p(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        !圧力差を求める phi1 = 2.211 phi2 = 4.895
        veff = 0.0d0
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    if(phi(xi,yi,zi) >= 0.5d0*(phi1+phi2)) then
                        veff = veff + ds**3
                    end if
                enddo
            enddo
        enddo
        reff = (3.d0/4.d0*veff / pi)**(1.d0/3.d0)
        p_in = p(xmax/2,ymax/2,zmax/2)
        p_out = p(0,ymax/2,zmax/2)

        dp_cal = p_in - p_out
        sigma = dp_cal*reff/2.0d0
        write(*,"(7es16.8)") dble(n),dp_cal,veff,reff,sigma,kappag
    end subroutine laplace

    subroutine break_up(phi)
        real(8),intent(in):: phi(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        !楕円の長辺・短辺をもとめる
        seff = 0.0d0
        distance = 0.0d0
        yi = ymax/2
        do zi=0,zmax
            do xi=0,xmax
                if(phi(xi,yi,zi) >= 0.5d0*(phi1+phi2)) then
                    seff = seff + ds**2
                    if( distance <= sqrt((dble(xi)-xc)**2 + (dble(zi)-zc)**2) ) then
                        distance = sqrt((dble(xi)-xc)**2 + (dble(zi)-zc)**2)
                    endif
                endif
            enddo
        enddo

        L1 = 2.0d0 * distance
        B1 = (4.0d0*seff) / (pi*L1)
        Dd = (L1-B1) / (L1+B1)
        write(*,"(6es16.8)") "step=",dble(n),veff,seff,L1,B1,Dd
        open(21,file="Dd-time.d")
        write(21,"(2es16.8)") (dble(n)-dble(start))*2.0d0*uw/H1, Dd
        close(21)
    end subroutine break_up

    subroutine output(phi,p,u1,u2,u3)
        real(8),intent(inout):: phi(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: p(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: u1(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: u2(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: u3(-1:xmax+1,-1:ymax+1,-1:zmax+1)

        write(filename,*) n !i->filename 変換
        filename=datadir//'step='//trim(adjustl(filename))//'.d' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
        print *, filename !表示してみる
        open(20,file=filename, status='replace') !使う
        zi = zmax/2
        do xi=0,xmax
            do yi=0,ymax
                write(20,"(8es16.8)") dble(xi),dble(yi),dble(zi),u1(xi,yi,zi),u2(xi,yi,zi),u3(xi,yi,zi),p(xi,yi,zi),phi(xi,yi,zi)
            enddo
            write(20,*)
        enddo
        close(20)
    end subroutine output

    !ディレクトリ作成
    subroutine mk_dirs(outdir)
        implicit none
        character(len=*), intent(in) :: outdir
        character(len=256) command
        write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
        call system(command)
    end subroutine mk_dirs

end module globals

program main
use globals
    implicit none
    real(8) feq(1:15, -1:xmax+1, -1:ymax+1, -1:zmax+1) !局所平衡分布関数（phi計算用）
    real(8) geq(1:15, -1:xmax+1, -1:ymax+1, -1:zmax+1) !局所平衡分布関数（速度・圧力計算用）
    real(8) phi(-1:xmax+1, -1:ymax+1, -1:zmax+1) !秩序パラメータ
    real(8) u1(-1:xmax+1, -1:ymax+1, -1:zmax+1), u2(-1:xmax+1, -1:ymax+1, -1:zmax+1),&
            u3(-1:xmax+1, -1:ymax+1, -1:zmax+1) !流速
    real(8) p(-1:xmax+1, -1:ymax+1, -1:zmax+1) !圧力
    real(8) Anu(0:xmax,0:ymax,0:zmax)
    real(8) nu(0:xmax,0:ymax,0:zmax)
    real(8) U(1:3)
    real(8),dimension(0:xmax,0:ymax,0:zmax):: p0
    real(8),dimension(0:xmax,0:ymax,0:zmax):: lap_phi
    real(8),dimension(1:3, 0:xmax,0:ymax,0:zmax):: grad_phi
    real(8),dimension(1:3, 1:3, 0:xmax,0:ymax,0:zmax):: gphi
    real(8),dimension(1:3, 1:3, 0:xmax,0:ymax,0:zmax):: pcap
    real(8),dimension(1:3, 0:xmax,0:ymax,0:zmax):: div_pcap
    real(8),dimension(1:3,1:3,0:xmax,0:ymax,0:zmax):: grad_u

    open(7,file="para.d")
    write(7,*) "Re, Ca, nu1, nu2, kappaf, kappag, C, H1, D, uw"
    write(7,"(10es13.5)") Re,Ca,nu1,nu2,kappaf,kappag,C,H1,D,uw
    write(7,*) "phi1, phi2, a, b, T"
    write(7,"(5es13.5)") phi1,phi2,a,b,T
    close(7)
    call mk_dirs(datadir)
    call par(cx,cy,cz,cr,krone,U)
    call ini(phi,u1,u2,u3,p,nu,Anu)

DO n=1,step
!======================のりしろ境界==========================================================
    call periodic(phi)
    call periodic(u1)
    call periodic(u2)
    call periodic(u3)
    call periodic(p)
!=======================div/gradの計算============================================================
    call lap_cal(lap_phi,phi)
    call grad_cal(grad_phi,phi)
    call p0_cal(p0,phi)
    call gphi_cal(gphi,grad_phi)
    call pcap_cal(pcap,p0,phi,lap_phi,grad_phi)
    call div_pcap_cal(div_pcap,pcap)
    call grad_u_cal(grad_u,u1,u2,u3)
!=======================平衡分布関数の計算========================================================
    call feq_cal(gphi,phi,p0,lap_phi,grad_phi,u1,u2,u3,div_pcap,feq)
    call geq_cal(gphi,grad_u,p,u1,u2,u3,Anu,geq)
    ! call bounce_back(feq,U)
    ! call bounce_back(geq,U)
    call periodic_fg(feq)
    call periodic_fg(geq)
!==========================物理量の計算==========================================================
    call physics(phi,p,u1,u2,u3,feq,geq,nu,Anu)
!============================静止液滴======================================================
    call laplace(phi,p)
!===========================液滴分裂========================================================
    ! call break_up(phi)
    if ( mod(n+1,100)==0 ) then
        ! call output(phi,p,u1,u2,u3)
    endif
ENDDO

    open(41,file="data1.d")
    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                do i=1,15
                    write(41,"(1es24.17)") feq(i,xi,yi,zi)
                enddo
            enddo
        enddo
    enddo
    close(41)
end program main