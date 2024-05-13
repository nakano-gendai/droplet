module globals
    real(8),parameter:: ds = 1.0d0 !格子間隔（lattice unit）
    integer,parameter:: xmax = 90 !ｘ方向格子数
    integer,parameter:: ymax = 90 !ｙ方向格子数
    integer,parameter:: zmax = 90 !ｚ方向格子数
    integer,parameter:: step = 10 !計算時間step
    integer i, j, k, n, xi, yi, zi, alpha, beta
    character :: filename*200
    character(*),parameter :: datadir = "/home/nakano/droplet/program"
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
    real(8),parameter:: D = 40.0d0 * ds !設置する液滴径
    real(8),parameter:: nu1 = 0.4d0 !液体の動粘度
    real(8),parameter:: nu2 = 0.4d0 !液滴の動粘度
    real(8),parameter:: uw = 0.02d0 !壁の移動速度
    real(8),parameter:: kappaf = 0.01d0*ds**2 !界面厚さを決めるパラメータ
    real(8),parameter:: kappag = 2.82d-3 !界面張力を決めるパラメータ
    real(8),parameter:: H1 = dble(zmax) !代表長さ
    real(8),parameter:: phi1 = -0.5d0
    real(8),parameter:: phi2 = 0.5d0
    real(8),parameter:: Anu = 0.0d0
    real(8),parameter:: A = 0.9d0
    real(8),parameter:: M = (1.0d0 - A) * ds / 6.0d0
    real(8),parameter:: W = 4.0d0 * ds
    real(8),parameter:: eps = 1.0d-12

    !初期液滴中心
    real(8), parameter:: xc = 0.5d0*dble(xmax)	
    real(8), parameter:: yc = 0.5d0*dble(ymax)
    real(8), parameter:: zc = 0.5d0*dble(zmax)

    !その他変数
    real(8) phi_min, phi_max, sigma, sigma_th, wa
    real(8) veff, seff, reff, distance
    real(8) p_in, p_out, dp_cal
    real(8) L1, B1, Dd
    real(8) gtemp, ftemp
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
    real(8) feq(1:15,-1:xmax+1,-1:ymax+1,-1:zmax+1) !平衡分布関数（オーダーパラメータ）
    real(8) geq(1:15,0:xmax,0:ymax,0:zmax) !平衡分布関数（速度・圧力計算用）
    real(8) phi(-1:xmax+1,-1:ymax+1,-1:zmax+1) !秩序パラメータ
    real(8) u1(-1:xmax+1,-1:ymax+1,-1:zmax+1), u2(-1:xmax+1,-1:ymax+1,-1:zmax+1), u3(-1:xmax+1,-1:ymax+1,-1:zmax+1) !流速
    real(8) p(0:xmax,0:ymax,0:zmax) !圧力
    real(8) U(1:3) !壁の移動速度
    real(8) grad_phi(1:3,0:xmax,0:ymax,0:zmax)
    real(8) gphi(1:3,1:3,0:xmax,0:ymax,0:zmax)
    real(8) grad_u(1:3,1:3,0:xmax,0:ymax,0:zmax)
    real(8) n_vector(1:3,0:xmax,0:ymax,0:zmax)
    real(8) phinext(0:xmax,0:ymax,0:zmax)

    call par(cx,cy,cz,cr,krone,U)
    !========================初期条件==========================================================================
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
    !周期境界
    do zi=0,zmax
        do yi=-1,ymax+1
            do xi=0,xmax
                phi(xi,-1,zi) = phi(xi,ymax,zi)
                phi(xi,ymax+1,zi) = phi(xi,0,zi)
                phi(-1,yi,zi) = phi(xmax,yi,zi)
                phi(xmax+1,yi,zi) = phi(0,yi,zi)

                u1(xi,-1,zi) = u1(xi,ymax,zi)
                u1(xi,ymax+1,zi) = u1(xi,0,zi)
                u1(-1,yi,zi) = u1(xmax,yi,zi)
                u1(xmax+1,yi,zi) = u1(0,yi,zi)

                u2(xi,-1,zi) = u2(xi,ymax,zi)
                u2(xi,ymax+1,zi) = u2(xi,0,zi)
                u2(-1,yi,zi) = u2(xmax,yi,zi)
                u2(xmax+1,yi,zi) = u2(0,yi,zi)

                u3(xi,-1,zi) = u3(xi,ymax,zi)
                u3(xi,ymax+1,zi) = u3(xi,0,zi)
                u3(-1,yi,zi) = u3(xmax,yi,zi)
                u3(xmax+1,yi,zi) = u3(0,yi,zi)
            enddo
        enddo
    enddo
    do yi=-1,ymax+1
        do xi=-1,xmax+1
            phi(xi,yi,-1) = phi(xi,yi,zmax)
            phi(xi,yi,zmax+1) = phi(xi,yi,0)

            u1(xi,yi,-1) = u1(xi,yi,zmax)
            u1(xi,yi,zmax+1) = u1(xi,yi,0)

            u2(xi,yi,-1) = u2(xi,yi,zmax)
            u2(xi,yi,zmax+1) = u2(xi,yi,0)

            u3(xi,yi,-1) = u3(xi,yi,zmax)
            u3(xi,yi,zmax+1) = u3(xi,yi,0)
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
    !n_vectorの計算
    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                do alpha=1,3
                    n_vector(alpha,xi,yi,zi) = grad_phi(alpha,xi,yi,zi) / &
                    (sqrt(grad_phi(1,xi,yi,zi)**2 + grad_phi(2,xi,yi,zi)**2 + grad_phi(3,xi,yi,zi)**2) + eps)
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
    !grad_uの計算
    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                do beta=1,3
                    grad_u(1,beta,xi,yi,zi) = 0.0d0
                    grad_u(2,beta,xi,yi,zi) = 0.0d0
                    grad_u(3,beta,xi,yi,zi) = 0.0d0
                    do i=2,15
                        grad_u(1,beta,xi,yi,zi) = grad_u(1,beta,xi,yi,zi)+cr(beta,i)*u1(xi+cx(i),yi+cy(i),zi+cz(i))
                        grad_u(2,beta,xi,yi,zi) = grad_u(2,beta,xi,yi,zi)+cr(beta,i)*u2(xi+cx(i),yi+cy(i),zi+cz(i))
                        grad_u(3,beta,xi,yi,zi) = grad_u(3,beta,xi,yi,zi)+cr(beta,i)*u3(xi+cx(i),yi+cy(i),zi+cz(i))
                    enddo
                    grad_u(1,beta,xi,yi,zi) = grad_u(1,beta,xi,yi,zi)/(10.0d0*ds)
                    grad_u(2,beta,xi,yi,zi) = grad_u(2,beta,xi,yi,zi)/(10.0d0*ds)
                    grad_u(3,beta,xi,yi,zi) = grad_u(3,beta,xi,yi,zi)/(10.0d0*ds)
                enddo
            enddo
        enddo
    enddo
    !fの初期条件
    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                do i=1,15
                    feq(i,xi,yi,zi) = E(i)*phi(xi,yi,zi)*(1.0d0+3.0d0*(cr(1,i) &
                    *u1(xi,yi,zi)+cr(2,i)*u2(xi,yi,zi)+cr(3,i)*u3(xi,yi,zi))) &
                                    +3.0d0*E(i)*M*(1.0d0-4.0d0*phi(xi,yi,zi)**2) &
                                    *(cr(1,i)*n_vector(1,xi,yi,zi)+cr(2,i)*n_vector &
                                    (2,xi,yi,zi)+cr(3,i)*n_vector(3,xi,yi,zi))/W
                enddo
            enddo
        enddo
    enddo
    !====================計算=========================================
DO n=1,step
    !周期境界
    do zi=0,zmax
        do yi=-1,ymax+1
            do xi=0,xmax
                phi(xi,-1,zi) = phi(xi,ymax,zi)
                phi(xi,ymax+1,zi) = phi(xi,0,zi)
                phi(-1,yi,zi) = phi(xmax,yi,zi)
                phi(xmax+1,yi,zi) = phi(0,yi,zi)

                u1(xi,-1,zi) = u1(xi,ymax,zi)
                u1(xi,ymax+1,zi) = u1(xi,0,zi)
                u1(-1,yi,zi) = u1(xmax,yi,zi)
                u1(xmax+1,yi,zi) = u1(0,yi,zi)

                u2(xi,-1,zi) = u2(xi,ymax,zi)
                u2(xi,ymax+1,zi) = u2(xi,0,zi)
                u2(-1,yi,zi) = u2(xmax,yi,zi)
                u2(xmax+1,yi,zi) = u2(0,yi,zi)

                u3(xi,-1,zi) = u3(xi,ymax,zi)
                u3(xi,ymax+1,zi) = u3(xi,0,zi)
                u3(-1,yi,zi) = u3(xmax,yi,zi)
                u3(xmax+1,yi,zi) = u3(0,yi,zi)
            enddo
        enddo
    enddo
    do yi=-1,ymax+1
        do xi=-1,xmax+1
            phi(xi,yi,-1) = phi(xi,yi,zmax)
            phi(xi,yi,zmax+1) = phi(xi,yi,0)

            u1(xi,yi,-1) = u1(xi,yi,zmax)
            u1(xi,yi,zmax+1) = u1(xi,yi,0)

            u2(xi,yi,-1) = u2(xi,yi,zmax)
            u2(xi,yi,zmax+1) = u2(xi,yi,0)

            u3(xi,yi,-1) = u3(xi,yi,zmax)
            u3(xi,yi,zmax+1) = u3(xi,yi,0)
        enddo
    enddo


    !grad_phiの計算
    open(65,file="grad.d")
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
    yi=ymax/2
    zi=zmax/2
    if(n==10)then
    do xi=0,xmax
        write(65,*) alpha,xi,yi,zi,grad_phi(alpha,xi,yi,zi)
    enddo
    endif
    close(65)
    !n_vectorの計算
    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                do alpha=1,3
                    n_vector(alpha,xi,yi,zi) = grad_phi(alpha,xi,yi,zi) &
                    / (sqrt(grad_phi(1,xi,yi,zi)**2 + grad_phi(2,xi,yi,zi)**2 + grad_phi(3,xi,yi,zi)**2) + eps)
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
    !grad_uの計算
    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                do beta=1,3
                    grad_u(1,beta,xi,yi,zi) = 0.0d0
                    grad_u(2,beta,xi,yi,zi) = 0.0d0
                    grad_u(3,beta,xi,yi,zi) = 0.0d0
                    do i=2,15
                        grad_u(1,beta,xi,yi,zi) = grad_u(1,beta,xi,yi,zi)+cr(beta,i)*u1(xi+cx(i),yi+cy(i),zi+cz(i))
                        grad_u(2,beta,xi,yi,zi) = grad_u(2,beta,xi,yi,zi)+cr(beta,i)*u2(xi+cx(i),yi+cy(i),zi+cz(i))
                        grad_u(3,beta,xi,yi,zi) = grad_u(3,beta,xi,yi,zi)+cr(beta,i)*u3(xi+cx(i),yi+cy(i),zi+cz(i))
                    enddo
                    grad_u(1,beta,xi,yi,zi) = grad_u(1,beta,xi,yi,zi)/(10.0d0*ds)
                    grad_u(2,beta,xi,yi,zi) = grad_u(2,beta,xi,yi,zi)/(10.0d0*ds)
                    grad_u(3,beta,xi,yi,zi) = grad_u(3,beta,xi,yi,zi)/(10.0d0*ds)
                enddo
            enddo
        enddo
    enddo


    !feq
    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                do i=1,15
                    feq(i,xi,yi,zi) = E(i)*phi(xi,yi,zi)*(1.0d0+3.0d0*(cr(1,i)* &
                    u1(xi,yi,zi)+cr(2,i)*u2(xi,yi,zi)+cr(3,i)*u3(xi,yi,zi))) &
                                    +3.0d0*E(i)*M*(1.0d0-4.0d0*phi(xi,yi,zi)**2) &
                                    *(cr(1,i)*n_vector(1,xi,yi,zi)+cr(2,i)* &
                                    n_vector(2,xi,yi,zi)+cr(3,i)*n_vector(3,xi,yi,zi))/W
                enddo
            enddo
        enddo
    enddo
    !周期境界
    do zi=0,zmax
        do yi=-1,ymax+1
            do xi=0,xmax
                do i=1,15
                    feq(i,xi,-1,zi) = feq(i,xi,ymax,zi)
                    feq(i,xi,ymax+1,zi) = feq(i,xi,0,zi)
                    feq(i,-1,yi,zi) = feq(i,xmax,yi,zi)
                    feq(i,xmax+1,yi,zi) = feq(i,0,yi,zi)
                enddo
            enddo
        enddo
    enddo
    do yi=-1,ymax+1
        do xi=-1,xmax+1
            do i=1,15
                feq(i,xi,yi,-1) = feq(i,xi,yi,zmax)
                feq(i,xi,yi,zmax+1) = feq(i,xi,yi,0)
            enddo
        enddo
    enddo
    !phinext
    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                phinext(xi,yi,zi) = 0.0d0
                do i=1,15
                    phinext(xi,yi,zi) = phinext(xi,yi,zi) + feq(i,xi-cx(i),yi-cy(i),zi-cz(i)) &
                     + A*E(i)*(phi(xi,yi,zi)-phi(xi-cx(i),yi-cy(i),zi-cz(i)))
                enddo
            enddo
        enddo
    enddo
    !renew
    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                phi(xi,yi,zi) = phinext(xi,yi,zi)
                u1(xi,yi,zi) = 0.0d0
                u2(xi,yi,zi) = 0.0d0
                u3(xi,yi,zi) = 0.0d0
                p(xi,yi,zi) = 0.0d0
            enddo
        enddo
    enddo


    ! write(filename,*) n !i->filename 変換
    ! filename=datadir//'step='//trim(adjustl(filename))//'.d' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
    ! print *, filename !表示してみる
    ! open(20,file=filename, status='replace') !使う

    ! yi=(ymax)/2
    ! zi=(zmax)/2
    ! do xi=0,xmax
    !     write(20,"(2es16.8)") dble(xi),phi(xi,yi,zi)
    ! enddo
    ! close(20)
ENDDO
end program main