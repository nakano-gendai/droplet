module globals
    include "mpif.h"
    real(8),parameter:: ds = 1.0d0 !格子間隔（lattice unit）
    integer,parameter:: xmax = 68 !ｘ方向格子数
    integer,parameter:: ymax = 68 !ｙ方向格子数
    integer,parameter:: zmax = 68 !ｚ方向格子数
    integer,parameter:: step = 10000000 !計算時間step
    integer,parameter:: start = 8000 !壁を動かす時間step
    integer i, n, xi, yi, zi, alpha, beta, k
    character :: filename*20
    character(*),parameter :: datadir = "./data_b/"

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
    real(8),parameter:: D = 32.0d0 * ds !設置する液滴径
    real(8),parameter:: nu1 = 0.1d0 !液体の動粘度
    real(8),parameter:: nu2 = 0.1d0 !液滴の動粘度
    real(8),parameter:: uw = 2.66d-3 !壁の移動速度
    real(8),parameter:: kappaf = 0.01d0*ds**2 !界面厚さを決めるパラメータ
    real(8),parameter:: kappag = 1.62d-4 !界面張力を決めるパラメータ
    real(8),parameter:: H1 = dble(zmax) !代表長さ
    real(8),parameter:: Re = 2.0d0*uw*D**2/(4.0d0*H1*nu1) !レイノルズ数
    real(8),parameter:: Ca = 2.0d0*uw*nu1*D/(2.0d0*H1*1.25d-3) !キャピラリー数
    real(8),parameter:: phi1 = 2.211d0
    real(8),parameter:: phi2 = 4.895d0
    real(8),parameter:: a = 9.0d0/49.0d0
    real(8),parameter:: b = 2.0d0/21.0d0
    real(8),parameter:: T = 0.55d0
    real(8),parameter:: C = 0.0d0
    real(8),parameter:: M = (0.5d0-C/3.0d0)*ds !モビリティー

    !初期液滴中心
    real(8), parameter:: xc = 0.5d0*dble(xmax)	
    real(8), parameter:: yc = 0.5d0*dble(ymax)
    real(8), parameter:: zc = 0.5d0*dble(zmax)
    
    !その他変数
    real(8) phi_min, phi_max, sigma
    real(8) veff, reff, seff, distance
    real(8) p_in, p_out, dp_cal
    real(8) L1, B1, Dd
    real(8) ftemp
    real(8) gtemp
    real(8) time1,time2
    real(8) kappag_tmp

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

end module globals

program main
use globals
    implicit none
    real(8) feq(1:15, -1:xmax+1, -1:ymax+1, -1:zmax+1) !局所平衡分布関数（phi計算用）
    real(8) geq(1:15, -1:xmax+1, -1:ymax+1, -1:zmax+1) !局所平衡分布関数（速度・圧力計算用）
    real(8) g(1:15, -1:xmax+1, -1:ymax+1, -1:zmax+1) !平衡分布関数（速度・圧力計算用）
    real(8) gnext(1:15, 0:xmax, 0:ymax, 0:zmax)
    real(8) phi(-1:xmax+1,-1:ymax+1,-1:zmax+1) !秩序パラメータ
    real(8) u1(0:xmax,0:ymax,0:zmax), u2(0:xmax,0:ymax,0:zmax),&
            u3(0:xmax,0:ymax,0:zmax) !流速
    real(8) p(0:xmax,0:ymax,0:zmax) !圧力
    real(8) U(1:3)
    real(8) grad_phi(1:3,0:xmax,0:ymax,0:zmax)
    real(8) gphi(1:3,1:3,0:xmax,0:ymax,0:zmax)
    real(8),dimension(1:3, 1:3, 0:xmax,0:ymax,0:zmax):: pcap
    real(8),dimension(1:3, 0:xmax,0:ymax,0:zmax):: div_pcap
    real(8),dimension(0:xmax,0:ymax,0:zmax):: p0
    real(8),dimension(0:xmax,0:ymax,0:zmax):: lap_phi
    real(8) nu(0:xmax,0:ymax,0:zmax)
    real(8) taug(-1:xmax+1,-1:ymax+1,-1:zmax+1)

    call par(cx,cy,cz,cr,krone,U)
!============初期条件======================================
    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                phi(xi,yi,zi) = phi1
                u1(xi,yi,zi) = 0.0d0
                u2(xi,yi,zi) = 0.0d0
                u3(xi,yi,zi) = 0.0d0
                p(xi,yi,zi) = 0.0d0
                !中央に液滴生成
                if((dble(xi)*ds-xc)**2+(dble(yi)*ds-yc)**2+(dble(zi)*ds-zc)**2 <= (0.5d0*D)**2) then
                    phi(xi,yi,zi) = phi2
                endif
            enddo
        enddo
    enddo
    !phiの周期境界
    do zi=0,zmax
        do yi=-1,ymax+1
            do xi=0,xmax
                phi(xi,-1,zi) = phi(xi,ymax,zi)
                phi(xi,ymax+1,zi) = phi(xi,0,zi)
                phi(-1,yi,zi) = phi(xmax,yi,zi)
                phi(xmax+1,yi,zi) = phi(0,yi,zi)
            enddo
        enddo
    enddo
    do yi=-1,ymax+1
        do xi=-1,xmax+1
            phi(xi,yi,-1) = phi(xi,yi,zmax)
            phi(xi,yi,zmax+1) = phi(xi,yi,0)
        enddo
    enddo
    !grad_phiの計算
    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                do alpha=1,3
                    grad_phi(alpha,xi,yi,zi) = 0.0d0
                    do i = 2,15
                        grad_phi(alpha,xi,yi,zi) = grad_phi(alpha,xi,yi,zi) + cr(alpha,i)*phi(xi+cx(i),yi+cy(i),zi+cz(i))
                    enddo
                    grad_phi(alpha,xi,yi,zi) = grad_phi(alpha,xi,yi,zi) / (10.0d0*ds)
                enddo
            enddo
        enddo
    enddo
    !gphiの計算
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
    !gの初期条件
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
                    g(i,xi,yi,zi) = E(i)*(3.0d0*p(xi,yi,zi) + 3.0d0*(u1(xi,yi,zi)*dble(cx(i)) &
                    + u2(xi,yi,zi)*dble(cy(i)) + u3(xi,yi,zi)*dble(cz(i))) &
                    + 4.5d0*(u1(xi,yi,zi)*dble(cx(i)) + u2(xi,yi,zi)*dble(cy(i)) + u3(xi,yi,zi)*dble(cz(i)))**2 &
                    -1.5d0*(u1(xi,yi,zi)**2 + u2(xi,yi,zi)**2 + u3(xi,yi,zi)**2)) 
                    ! + E(i)*kappag*gtemp
                enddo
            enddo
        enddo
    enddo

end program