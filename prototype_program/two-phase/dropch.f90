! 【Cahn-Hilliard 方程式】
! このプログラムは、界面を「Cahn-Hilliard 方程式」で表現している。
! 液滴を長時間保持できないが、系全体の保存性には優れている。
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module globals
    include "mpif.h"
    real(8),parameter:: ds = 1.0d0 !格子間隔（lattice unit）
    integer,parameter:: xmax = 200 !ｘ方向格子数
    integer,parameter:: ymax = 200 !ｙ方向格子数
    integer,parameter:: zmax = 200 !ｚ方向格子数
    integer,parameter:: step = 5000 !計算時間step
    integer,parameter:: start = 5001 !壁を動かし始める時間step
    integer i, j, k, n, xi, yi, zi, alpha, beta
    character :: filename*200
    character(*),parameter :: datadir = "/data/n/n517/sigmatest/"
    ! character(*),parameter :: datadir = "/home/nakano/aomsin/data_aomsin/"
    ! character(*),parameter :: datadir = "/data/2022/nakano/multitest/"
    ! character(*),parameter :: datadir = "/data/2022/nakano/rechange_caconst/re150ca0.2/"
    

    !MPI用変数
    integer ierr,comm_procs,comm_rank
    integer Nx,Ny !xi、yi方向の分割数
    integer x_procs,y_procs
    integer key_x, group_x, key_y, group_y
    integer newy_comm_world,newy_procs,newy_rank
    integer newx_comm_world,newx_procs,newx_rank
    integer rectangle_type,rectangle_type_f,xtype,xtype_f,output_type
    integer next_rank_x,former_rank_x,next_rank_y,former_rank_y
    integer req1s,req1r,req2s,req2r,req3s,req3r,req4s,req4r,req5s,req5r,req6s,req6r,req7s,req7r,req8s,req8r,req9s,req9r,req10s,req10r
    integer, dimension(MPI_STATUS_SIZE) :: sta1s,sta1r,sta2s,sta2r,sta3s,sta3r,sta4s,sta4r,sta5s,sta5r,sta6s,sta6r,sta7s,sta7r,sta8s,sta8r,sta9s,sta9r,sta10s,sta10r
    integer Nyy,Nxx
    real(8), allocatable:: tmp1(:,:,:,:),tmp2(:,:,:,:),tmp3(:,:,:,:),tmp4(:,:,:,:),tmp5(:,:,:,:),tmpf(:,:,:,:,:),tmpff(:,:,:,:,:),tmp(:,:,:),tmp_f(:,:,:,:),tmp_nu(:,:,:,:)
    integer tags, tagr, recv_rank
    character(2) chmyrank

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
    real(8),parameter:: D = 50.0d0 * ds !設置する液滴径
    real(8),parameter:: nu1 = 0.01d0 !液体の動粘度
    real(8),parameter:: nu2 = 0.01d0 !液滴の動粘度
    real(8),parameter:: uw = 0.02d0 !壁の移動速度
    real(8),parameter:: kappaf = 0.01d0*ds**2 !界面厚さを決めるパラメータ
    real(8),parameter:: kappag = 1.8d-2 !界面張力を決めるパラメータ
    real(8),parameter:: H1 = dble(zmax) !代表長さ
    real(8),parameter:: Re = 2.0d0*uw*D**2/(4.0d0*H1*nu1) !レイノルズ数
    real(8),parameter:: Ca = 2.0d0*uw*nu1*D/(2.0d0*H1*0.02d0) !キャピラリー数
    real(8),parameter:: phi1 = 2.211d0
    real(8),parameter:: phi2 = 4.895d0
    ! real(8),parameter:: phi2 = 2.211d0
    real(8),parameter:: a = 9.0d0/49.0d0
    real(8),parameter:: b = 2.0d0/21.0d0
    real(8),parameter:: T = 0.55d0
    real(8),parameter:: tauf = 0.7d0
    real(8),parameter:: Anu = 0.0d0

    !初期液滴中心
    real(8), parameter:: xc = 0.5d0*dble(xmax)	
    real(8), parameter:: yc = 0.5d0*dble(ymax)
    real(8), parameter:: zc = 0.5d0*dble(zmax)
    
    !その他変数
    real(8) phi_min, phi_max, sigma, sigma_th, wa
    real(8) min,max
    real(8) veff, seff, reff, distance
    real(8) p_in, p_out, dp_cal
    real(8) L1, B1, Dd
    real(8) gtemp, ftemp
    real(8) time1,time2
    real(8) kappag_tmp
    !重心計算用変数
    real(8) xg, yg, zg
    real(8) phi_tmp
    !慣性行列の固有値・固有ベクトル計算用変数
    real(8) xx_wa, yy_wa, zz_wa, xy_wa, yz_wa, zx_wa
    real(8) inertia_matrix(1:3,1:3), inertia_matrix_tmp(1:3,1:3)
    real(8) principal_moment_1, principal_moment_2, principal_moment_3
    real(8) inertia_principal_axis_1(1:3), inertia_principal_axis_2(1:3), inertia_principal_axis_3(1:3)
    real(8) w(1:3,1:3), w_t(1:3,1:3), t0(1:3,1:3), t_tmp(1:3,1:3)
    real(8) theta, mx
    integer pp, qq, check 
    !線形補間用変数
    real(8), parameter :: weight = 0.01d0
    real(8) xj, zj, phi_hokann
    integer x1, x2, z1, z2

    integer x_tmp, z_tmp
    real(8),parameter :: zz = 0.1d0
    real(8),parameter :: xx = 0.1d0
    

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

    subroutine ini(phi,u1,u2,u3,p,geq,g_procs,gnext_procs,feq,f_procs,fnext_procs,feq_procs,geq_procs,phi_procs,p_procs,u1_procs,u2_procs,u3_procs,grad_phi,gphi,lap_phi,p0,grad_u)
        real(8),intent(inout):: phi(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: u1(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: u2(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: u3(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: p(0:xmax,0:ymax,0:zmax)
        real(8),intent(inout):: feq(1:15,0:xmax,0:ymax,0:zmax)
        real(8),intent(inout):: geq(1:15,0:xmax,0:ymax,0:zmax)
        real(8),intent(inout):: grad_phi(1:3,0:xmax,0:ymax,0:zmax)
        real(8),intent(inout):: lap_phi(0:xmax,0:ymax,0:zmax)
        real(8),intent(inout):: gphi(1:3,1:3,0:xmax,0:ymax,0:zmax)
        real(8),intent(inout):: p0(0:xmax,0:ymax,0:zmax)
        real(8),intent(inout):: grad_u(1:3,1:3,0:xmax,0:ymax,0:zmax)
        
        real(8),allocatable :: f_procs(:,:,:,:), fnext_procs(:,:,:,:), feq_procs(:,:,:,:)
        real(8),allocatable :: g_procs(:,:,:,:), gnext_procs(:,:,:,:), geq_procs(:,:,:,:)
        real(8),allocatable :: phi_procs(:,:,:), p_procs(:,:,:)
        real(8),allocatable :: u1_procs(:,:,:), u2_procs(:,:,:), u3_procs(:,:,:)

        !配列のallocate(xi:0とx_procs+1がのりしろ)(yi:0とymax+2がのりしろ)(zi:0とy_procs+1がのりしろ)
        Nx = 3 !ｘ方向の並列数（ただし，Nx/=comm_procs）
        Ny = comm_procs / Nx !ｙ方向の並列数
        x_procs = (xmax+1) / Nx
        y_procs = (zmax+1) / Ny
        !以下はのりしろ有りの変数
        allocate(phi_procs(0:x_procs+1,0:ymax+2,0:y_procs+1))
        allocate(p_procs(0:x_procs+1,0:ymax+2,0:y_procs+1))
        allocate(u1_procs(0:x_procs+1,0:ymax+2,0:y_procs+1))
        allocate(u2_procs(0:x_procs+1,0:ymax+2,0:y_procs+1))
        allocate(u3_procs(0:x_procs+1,0:ymax+2,0:y_procs+1))
        allocate(f_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1))
        allocate(fnext_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1))
        allocate(feq_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1))
        allocate(g_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1))
        allocate(gnext_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1))
        allocate(geq_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1))
        allocate(tmp1(0:comm_procs-1,1:x_procs,1:ymax+1,1:y_procs))
        allocate(tmp2(0:comm_procs-1,1:x_procs,1:ymax+1,1:y_procs))
        allocate(tmp3(0:comm_procs-1,1:x_procs,1:ymax+1,1:y_procs))
        allocate(tmp4(0:comm_procs-1,1:x_procs,1:ymax+1,1:y_procs))
        allocate(tmp5(0:comm_procs-1,1:x_procs,1:ymax+1,1:y_procs))
        allocate(tmpf(0:comm_procs-1,1:15,1:x_procs,1:ymax+1,1:y_procs))
        allocate(tmpff(0:comm_procs-1,1:15,1:x_procs,1:ymax+1,1:y_procs))
        allocate(tmp(1:x_procs,1:ymax+1,1:y_procs))
        allocate(tmp_nu(0:comm_procs-1,1:x_procs,1:ymax+1,1:y_procs))
        
        !初期化
        phi_procs(:,:,:) = 0.0d0
        p_procs(:,:,:) = 0.0d0
        u1_procs(:,:,:) = 0.0d0
        u2_procs(:,:,:) = 0.0d0
        u3_procs(:,:,:) = 0.0d0
        feq_procs(:,:,:,:) = 0.0d0
        f_procs(:,:,:,:) = 0.0d0
        fnext_procs(:,:,:,:) = 0.0d0
        geq_procs(:,:,:,:) = 0.0d0
        g_procs(:,:,:,:) = 0.0d0
        gnext_procs(:,:,:,:) = 0.0d0
    !===================================初期条件の設定================================================
        !phi, u1, u2, u3, p の初期条件
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
        !lap_phiの計算
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    lap_phi(xi,yi,zi) = -14.0d0*phi(xi,yi,zi)
                    do i=2,15
                        lap_phi(xi,yi,zi) = lap_phi(xi,yi,zi) + phi(xi+cx(i),yi+cy(i),zi+cz(i))
                    enddo
                    lap_phi(xi,yi,zi) = lap_phi(xi,yi,zi) / (5.0d0*ds**2)
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
        !p0の計算
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    p0(xi,yi,zi) = phi(xi,yi,zi)*T/(1.0d0-b*phi(xi,yi,zi)) - a*phi(xi,yi,zi)**2
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
                        gtemp = 0.0d0
                        do beta=1,3
                            do alpha=1,3
                                gtemp = gtemp + gphi(alpha,beta,xi,yi,zi)*cr(alpha,i)*cr(beta,i)
                            enddo
                        enddo
                        feq(i,xi,yi,zi) = H(i)*phi(xi,yi,zi) &
                                        + F(i)*(p0(xi,yi,zi)-kappaf*phi(xi,yi,zi)*lap_phi(xi,yi,zi) &
                                        -kappaf/6.0d0*(grad_phi(1,xi,yi,zi)**2+grad_phi(2,xi,yi,zi)**2+grad_phi(3,xi,yi,zi)**2)) &
                                        + E(i)*phi(xi,yi,zi)*3.0d0*(cr(1,i)*u1(xi,yi,zi)+cr(2,i)*u2(xi,yi,zi)+cr(3,i)*u3(xi,yi,zi)) &
                                        + E(i)*kappaf*gtemp
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
                        - 1.5d0*(u1(xi,yi,zi)**2 + u2(xi,yi,zi)**2 + u3(xi,yi,zi)**2) &
                        + Anu*ds*ftemp) &
                        + E(i)*kappag*gtemp
                    enddo
                enddo
            enddo
        enddo
    !=======================================================================================================================
    !============================================初期条件の分割==============================================================
        n = 0
        do Nyy=0,Ny-1
            do Nxx=0,Nx-1
                do zi=1,y_procs
                    do yi=1,ymax+1
                        do xi=1,x_procs
                            tmp1(n,xi,yi,zi) = u1((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs))
                            tmp2(n,xi,yi,zi) = u2((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs))
                            tmp3(n,xi,yi,zi) = u3((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs))
                            tmp4(n,xi,yi,zi) = p((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs))
                            tmp5(n,xi,yi,zi) = phi((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs))
                            do i=1,15
                                tmpf(n,i,xi,yi,zi) = geq(i,(xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs))
                                tmpff(n,i,xi,yi,zi) = feq(i,(xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs))
                            enddo
                        enddo
                    enddo
                enddo
                n = n + 1
            enddo
        enddo
        do n=0,comm_procs-1
            if(n == comm_rank) then
                do zi=1,y_procs
                    do yi=1,ymax+1
                        do xi=1,x_procs
                            u1_procs(xi,yi,zi) = tmp1(n,xi,yi,zi)
                            u2_procs(xi,yi,zi) = tmp2(n,xi,yi,zi)
                            u3_procs(xi,yi,zi) = tmp3(n,xi,yi,zi)
                            p_procs(xi,yi,zi) = tmp4(n,xi,yi,zi)
                            phi_procs(xi,yi,zi) = tmp5(n,xi,yi,zi)
                            do i=1,15
                                g_procs(i,xi,yi,zi) = tmpf(n,i,xi,yi,zi)
                                f_procs(i,xi,yi,zi) = tmpff(n,i,xi,yi,zi)
                            enddo
                        enddo
                    enddo
                enddo
            endif
        enddo
    !====================================================================================================================
    !===========================================コミュニケータ・通信先設定==================================================
        !z方向にコミュニケータを分割する
        key_y = comm_rank
        group_y = comm_rank / Nx
        call MPI_Comm_Split(MPI_COMM_WORLD,group_y,key_y,newy_comm_world,ierr)
        call MPI_Comm_Size(newy_comm_world,newy_procs,ierr)
        call MPI_Comm_Rank(newy_comm_world,newy_rank,ierr)

        !x方向にコミュニケータを分割する
        key_x = comm_rank
        group_x = mod(comm_rank,Nx)
        call MPI_Comm_Split(MPI_COMM_WORLD,group_x,key_x,newx_comm_world,ierr)
        call MPI_Comm_Size(newx_comm_world,newx_procs,ierr)
        call MPI_Comm_Rank(newx_comm_world,newx_rank,ierr)

        !のりしろ境界の通信先設定
        next_rank_x = newx_rank + 1
        former_rank_x = newx_rank - 1
        if(newx_rank == 0) then
            former_rank_x = Ny - 1
        else if(newx_rank == Ny - 1) then
            next_rank_x = 0
        endif
        next_rank_y = newy_rank + 1
        former_rank_y = newy_rank - 1
        if(newy_rank == 0) then
            former_rank_y = Nx - 1
        else if(newy_rank == Nx - 1) then
            next_rank_y = 0
        endif

        !X世界での受け渡しをする際の型作成
        call MPI_Type_Vector(ymax+3,x_procs,x_procs+2,MPI_REAL8,xtype,ierr)
        call MPI_Type_Commit(xtype,ierr)
        call MPI_Type_Vector(ymax+3,15*x_procs,15*(x_procs+2),MPI_REAL8,xtype_f,ierr)
        call MPI_Type_Commit(xtype_f,ierr)

        !Y世界での受け渡しをする際の型作成
        call MPI_Type_Vector((y_procs+2)*(ymax+3),1,x_procs+2,MPI_REAL8,rectangle_type,ierr)
        call MPI_Type_Commit(rectangle_type,ierr)
        call MPI_Type_Vector((y_procs+2)*(ymax+3),15,15*(x_procs+2),MPI_REAL8,rectangle_type_f,ierr)
        call MPI_Type_Commit(rectangle_type_f,ierr)
    !===============================================================================================================
    ! open(77,file="origin.d")
    ! do zi=0,zmax
    !     do yi=0,ymax
    !         do xi=0,xmax
    !             write(77,"(8es24.16)") dble(xi),dble(yi),dble(zi),phi(xi,yi,zi)
    !         enddo
    !     enddo
    ! enddo
    ! close(77)
    end subroutine ini

    subroutine ini_op(taug_procs,nu_procs,p0_procs,lap_phi_procs,grad_phi_procs,gphi_procs,grad_u_procs,phi_procs)
        real(8),intent(in):: phi_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),allocatable :: taug_procs(:,:,:), nu_procs(:,:,:)
        real(8),allocatable :: p0_procs(:,:,:), lap_phi_procs(:,:,:), grad_phi_procs(:,:,:,:), gphi_procs(:,:,:,:,:), grad_u_procs(:,:,:,:,:)
        !以下はのりしろ無しの変数
        allocate(taug_procs(0:x_procs+1,0:ymax+2,0:y_procs+1))
        allocate(nu_procs(1:x_procs,1:ymax+1,1:y_procs))
        allocate(p0_procs(1:x_procs,1:ymax+1,1:y_procs))
        allocate(lap_phi_procs(1:x_procs,1:ymax+1,1:y_procs))
        allocate(grad_phi_procs(1:3,1:x_procs,1:ymax+1,1:y_procs))
        allocate(gphi_procs(1:3,1:3,1:x_procs,1:ymax+1,1:y_procs))
        allocate(grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:y_procs))

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    nu_procs(xi,yi,zi) = (phi_procs(xi,yi,zi)-phi1)/(phi2-phi1)*(nu2-nu1) + nu1
                    taug_procs(xi,yi,zi) = 0.5d0 + 3.0d0*nu_procs(xi,yi,zi)/ds + 2.0d0/3.0d0*Anu
                enddo
            enddo
        enddo
    end subroutine ini_op

    subroutine glue(var)
        real(8),intent(inout) :: var(0:x_procs+1,0:ymax+2,0:y_procs+1)
        do zi=1,y_procs
            do xi=1,x_procs
                var(xi,0,zi) = var(xi,ymax+1,zi)
                var(xi,ymax+2,zi) = var(xi,1,zi)
            enddo
        enddo
        !X世界でののりしろ通信
        call MPI_Isend(var(1,0,y_procs),1,xtype,next_rank_x,1,newx_comm_world,req1s,ierr)
        call MPI_Irecv(var(1,0,0),1,xtype,former_rank_x,1,newx_comm_world,req1r,ierr)

        call MPI_Isend(var(1,0,1),1,xtype,former_rank_x,2,newx_comm_world,req2s,ierr)
        call MPI_Irecv(var(1,0,y_procs+1),1,xtype,next_rank_x,2,newx_comm_world,req2r,ierr)

        call MPI_Wait(req1s,sta1s,ierr)
        call MPI_Wait(req1r,sta1r,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
        call MPI_Wait(req2r,sta2r,ierr)
        !Y世界でののりしろ通信
        call MPI_Isend(var(x_procs,0,0),1,rectangle_type,next_rank_y,1,newy_comm_world,req1s,ierr)
        call MPI_Irecv(var(0,0,0),1,rectangle_type,former_rank_y,1,newy_comm_world,req1r,ierr)

        call MPI_Isend(var(1,0,0),1,rectangle_type,former_rank_y,2,newy_comm_world,req2s,ierr)
        call MPI_Irecv(var(x_procs+1,0,0),1,rectangle_type,next_rank_y,2,newy_comm_world,req2r,ierr)

        call MPI_Wait(req1s,sta1s,ierr)
        call MPI_Wait(req1r,sta1r,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
        call MPI_Wait(req2r,sta2r,ierr)
    end subroutine glue

    subroutine glue_f(var_f)
        real(8),intent(inout) :: var_f(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        do zi=1,y_procs
            do xi=1,x_procs
                do i=1,15
                    var_f(i,xi,0,zi) = var_f(i,xi,ymax+1,zi)
                    var_f(i,xi,ymax+2,zi) = var_f(i,xi,1,zi)
                enddo
            enddo
        enddo
        !X世界でののりしろ通信
        call MPI_Isend(var_f(1,1,0,y_procs),1,xtype_f,next_rank_x,1,newx_comm_world,req1s,ierr)
        call MPI_Irecv(var_f(1,1,0,0),1,xtype_f,former_rank_x,1,newx_comm_world,req1r,ierr)

        call MPI_Isend(var_f(1,1,0,1),1,xtype_f,former_rank_x,2,newx_comm_world,req2s,ierr)
        call MPI_Irecv(var_f(1,1,0,y_procs+1),1,xtype_f,next_rank_x,2,newx_comm_world,req2r,ierr)

        call MPI_Wait(req1s,sta1s,ierr)
        call MPI_Wait(req1r,sta1r,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
        call MPI_Wait(req2r,sta2r,ierr)

        !Y世界でののりしろ通信
        call MPI_Isend(var_f(1,x_procs,0,0),1,rectangle_type_f,next_rank_y,1,newy_comm_world,req1s,ierr)
        call MPI_Irecv(var_f(1,0,0,0),1,rectangle_type_f,former_rank_y,1,newy_comm_world,req1r,ierr)

        call MPI_Isend(var_f(1,1,0,0),1,rectangle_type_f,former_rank_y,2,newy_comm_world,req2s,ierr)
        call MPI_Irecv(var_f(1,x_procs+1,0,0),1,rectangle_type_f,next_rank_y,2,newy_comm_world,req2r,ierr)

        call MPI_Wait(req1s,sta1s,ierr)
        call MPI_Wait(req1r,sta1r,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
        call MPI_Wait(req2r,sta2r,ierr)
    end subroutine glue_f

    subroutine MPI_boundary(phi_procs,u1_procs,u2_procs,u3_procs,f_procs,g_procs,taug_procs)
        real(8),intent(inout) :: phi_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(inout) :: u1_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(inout) :: u2_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(inout) :: u3_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(inout) :: f_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(inout) :: g_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(inout) :: taug_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        call glue(phi_procs)
        call glue(u1_procs)
        call glue(u2_procs)
        call glue(u3_procs)
        call glue(taug_procs)
        call glue_f(f_procs)
        call glue_f(g_procs)
    endsubroutine MPI_boundary

    subroutine MPI_boundary_fg(feq_procs,geq_procs)
        real(8),intent(inout) :: feq_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(inout) :: geq_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        call glue_f(feq_procs)
        call glue_f(geq_procs)
    endsubroutine MPI_boundary_fg

    subroutine feq_cal(gphi_procs,phi_procs,p0_procs,lap_phi_procs,grad_phi_procs,u1_procs,u2_procs,u3_procs,feq_procs)
        real(8),intent(in):: gphi_procs(1:3,1:3,1:x_procs,1:ymax+1,1:y_procs)
        real(8),intent(in):: phi_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(in):: p0_procs(1:x_procs,1:ymax+1,1:y_procs)
        real(8),intent(in):: lap_phi_procs(1:x_procs,1:ymax+1,1:y_procs)
        real(8),intent(in):: grad_phi_procs(1:3,1:x_procs,1:ymax+1,1:y_procs)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(out):: feq_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        gtemp = 0.0d0
                        do beta=1,3
                            do alpha=1,3
                                gtemp = gtemp + gphi_procs(alpha,beta,xi,yi,zi)*cr(alpha,i)*cr(beta,i)
                            enddo
                        enddo
                        feq_procs(i,xi,yi,zi) = H(i)*phi_procs(xi,yi,zi) &
                                        + F(i)*(p0_procs(xi,yi,zi)-kappaf*phi_procs(xi,yi,zi)*lap_phi_procs(xi,yi,zi) &
                                        -kappaf/6.0d0*(grad_phi_procs(1,xi,yi,zi)**2+grad_phi_procs(2,xi,yi,zi)**2+grad_phi_procs(3,xi,yi,zi)**2)) &
                                        + E(i)*phi_procs(xi,yi,zi)*3.0d0*(cr(1,i)*u1_procs(xi,yi,zi)+cr(2,i)*u2_procs(xi,yi,zi)+cr(3,i)*u3_procs(xi,yi,zi)) &
                                        + E(i)*kappaf*gtemp
                    enddo
                enddo
            enddo
        enddo
    end subroutine feq_cal

    subroutine geq_cal(gphi_procs,p_procs,u1_procs,u2_procs,u3_procs,geq_procs,cr,grad_u_procs)
        real(8),intent(in):: gphi_procs(1:3,1:3,1:x_procs,1:ymax+1,1:y_procs)
        real(8),intent(in):: p_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(out):: geq_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(in):: grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:y_procs)
        real(8),intent(in):: cr(1:3,1:15)
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        gtemp = 0.0d0
                        ftemp = 0.0d0
                        do beta=1,3
                            do alpha=1,3
                                gtemp = gtemp + gphi_procs(alpha,beta,xi,yi,zi)*cr(alpha,i)*cr(beta,i)

                                ftemp = ftemp + (grad_u_procs(beta,alpha,xi,yi,zi)+grad_u_procs(alpha,beta,xi,yi,zi))*cr(alpha,i)*cr(beta,i)
                            enddo
                        enddo
                        geq_procs(i,xi,yi,zi) = E(i)*(3.0d0*p_procs(xi,yi,zi) + 3.0d0*(u1_procs(xi,yi,zi)*dble(cx(i)) &
                        + u2_procs(xi,yi,zi)*dble(cy(i)) + u3_procs(xi,yi,zi)*dble(cz(i))) &
                        + 4.5d0*(u1_procs(xi,yi,zi)*dble(cx(i)) + u2_procs(xi,yi,zi)*dble(cy(i)) + u3_procs(xi,yi,zi)*dble(cz(i)))**2 &
                        - 1.5d0*(u1_procs(xi,yi,zi)**2 + u2_procs(xi,yi,zi)**2 + u3_procs(xi,yi,zi)**2) &
                        + Anu*ds*ftemp) &
                        + E(i)*kappag_tmp*gtemp
                    enddo
                enddo
            enddo
        enddo
    end subroutine geq_cal

    subroutine f_cal(fnext_procs,f_procs,feq_procs)
        real(8),intent(inout):: fnext_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(in):: f_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(in):: feq_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        fnext_procs(i,xi,yi,zi) = f_procs(i,xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                - (f_procs(i,xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                - feq_procs(i,xi-cx(i),yi-cy(i),zi-cz(i))) / tauf
                    enddo
                enddo
            enddo
        enddo
    end subroutine f_cal

    subroutine g_cal(gnext_procs,g_procs,geq_procs,taug_procs)
        real(8),intent(inout):: gnext_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(in):: g_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(in):: geq_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(in):: taug_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        gnext_procs(i,xi,yi,zi) = g_procs(i,xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                - (g_procs(i,xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                - geq_procs(i,xi-cx(i),yi-cy(i),zi-cz(i))) / taug_procs(xi-cx(i),yi-cy(i),zi-cz(i))
                    enddo
                enddo
            enddo
        enddo
    end subroutine g_cal

    subroutine bounce_back_LBM(fun,U)
        real(8),intent(inout):: fun(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(inout):: U(1:3)
        !安定したら壁を動かす
        if((n >= start)) then
            U(1) = uw
            U(2) = 0.0d0
            U(3) = 0.0d0
        endif
        if(group_y == 0) then
            zi = 1
            do yi=1,ymax+1
                do xi=1,x_procs
                    fun(4,xi,yi,zi) = fun(7,xi,yi,zi) - 6.0d0 * E(7) * ((-U(1))*cr(1,7)+(-U(2))*cr(2,7)+(-U(3))*cr(3,7))
                    fun(8,xi,yi,zi) = fun(12,xi,yi,zi) - 6.0d0 * E(12) * ((-U(1))*cr(1,12)+(-U(2))*cr(2,12)+(-U(3))*cr(3,12))
                    fun(9,xi,yi,zi) = fun(13,xi,yi,zi) - 6.0d0 * E(13) * ((-U(1))*cr(1,13)+(-U(2))*cr(2,13)+(-U(3))*cr(3,13))
                    fun(10,xi,yi,zi) = fun(14,xi,yi,zi) - 6.0d0 * E(14) * ((-U(1))*cr(1,14)+(-U(2))*cr(2,14)+(-U(3))*cr(3,14))
                    fun(15,xi,yi,zi) = fun(11,xi,yi,zi) - 6.0d0 * E(11) * ((-U(1))*cr(1,11)+(-U(2))*cr(2,11)+(-U(3))*cr(3,11))
                enddo
            enddo
        else if(group_y == Ny-1) then
            zi = y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    fun(7,xi,yi,zi) = fun(4,xi,yi,zi) - 6.0d0 * E(4) * ((U(1))*cr(1,4)+(U(2))*cr(2,4)+(U(3))*cr(3,4))
                    fun(11,xi,yi,zi) = fun(15,xi,yi,zi) - 6.0d0 * E(15) * ((U(1))*cr(1,15)+(U(2))*cr(2,15)+(U(3))*cr(3,15))
                    fun(12,xi,yi,zi) = fun(8,xi,yi,zi) - 6.0d0 * E(8) * ((U(1))*cr(1,8)+(U(2))*cr(2,8)+(U(3))*cr(3,8))
                    fun(13,xi,yi,zi) = fun(9,xi,yi,zi) - 6.0d0 * E(9) * ((U(1))*cr(1,9)+(U(2))*cr(2,9)+(U(3))*cr(3,9))
                    fun(14,xi,yi,zi) = fun(10,xi,yi,zi) - 6.0d0 * E(10) * ((U(1))*cr(1,10)+(U(2))*cr(2,10)+(U(3))*cr(3,10))
                enddo
            enddo
        endif
    end subroutine bounce_back_LBM

    subroutine physics(phi_procs,p_procs,u1_procs,u2_procs,u3_procs,f_procs,g_procs,nu_procs,taug_procs)
        real(8),intent(inout):: phi_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(inout):: p_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(inout):: u1_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(inout):: u2_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(inout):: u3_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(in):: f_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(in):: g_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(inout):: nu_procs(1:x_procs,1:ymax+1,1:y_procs)
        real(8),intent(inout):: taug_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    phi_procs(xi,yi,zi) = 0.0d0
                    p_procs(xi,yi,zi) = 0.0d0
                    u1_procs(xi,yi,zi) = 0.0d0
                    u2_procs(xi,yi,zi) = 0.0d0
                    u3_procs(xi,yi,zi) = 0.0d0
                    do i=1,15
                        phi_procs(xi,yi,zi) = phi_procs(xi,yi,zi) + f_procs(i,xi,yi,zi)
                        p_procs(xi,yi,zi) = p_procs(xi,yi,zi) + g_procs(i,xi,yi,zi) / 3.0d0
                        u1_procs(xi,yi,zi) = u1_procs(xi,yi,zi) + dble(cx(i))*g_procs(i,xi,yi,zi)
                        u2_procs(xi,yi,zi) = u2_procs(xi,yi,zi) + dble(cy(i))*g_procs(i,xi,yi,zi)
                        u3_procs(xi,yi,zi) = u3_procs(xi,yi,zi) + dble(cz(i))*g_procs(i,xi,yi,zi)
                    enddo
                enddo
            enddo
        enddo
        min = phi2
        max = phi1
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    if(phi_procs(xi,yi,zi) > max) then
                        max = phi_procs(xi,yi,zi)
                    elseif(phi_procs(xi,yi,zi) < min) then
                        min = phi_procs(xi,yi,zi)
                    endif
                enddo
            enddo
        enddo
        call MPI_Allreduce(max,phi_max,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(min,phi_min,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ierr)
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    nu_procs(xi,yi,zi) = (phi_procs(xi,yi,zi)-phi_min)/(phi_max-phi_min)*(nu2-nu1) + nu1
                    taug_procs(xi,yi,zi) = 0.5d0 + 3.0d0*nu_procs(xi,yi,zi)/ds + 2.0d0/3.0d0*Anu
                enddo
            enddo
        enddo
    end subroutine physics

    subroutine renew(fun_procs,funnext_procs)
        real(8),intent(in) :: funnext_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(out) :: fun_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        fun_procs(i,xi,yi,zi) = funnext_procs(i,xi,yi,zi)
                    enddo
                enddo
            enddo
        enddo
    end subroutine renew

    subroutine lap_cal(lap_f,fun)
        real(8),intent(out):: lap_f(1:x_procs,1:ymax+1,1:y_procs)
        real(8),intent(in):: fun(0:x_procs+1,0:ymax+2,0:y_procs+1)
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
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
        real(8),intent(out):: grad_f(1:3,1:x_procs,1:ymax+1,1:y_procs)
        real(8),intent(in):: fun(0:x_procs+1,0:ymax+2,0:y_procs+1)
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
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

    subroutine p0_cal(p0_procs,phi_procs)
        real(8),intent(out):: p0_procs(1:x_procs,1:ymax+1,1:y_procs)
        real(8),intent(in):: phi_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    p0_procs(xi,yi,zi) = phi_procs(xi,yi,zi)*T/(1.0d0-b*phi_procs(xi,yi,zi)) - a*phi_procs(xi,yi,zi)**2
                enddo
            enddo
        enddo
    end subroutine p0_cal

    subroutine gphi_cal(gphi_procs,grad_phi_procs)
        real(8),intent(out):: gphi_procs(1:3,1:3,1:x_procs,1:ymax+1,1:y_procs)
        real(8),intent(in):: grad_phi_procs(1:3,1:x_procs,1:ymax+1,1:y_procs)
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do beta=1,3
                        do alpha=1,3
                            gphi_procs(alpha,beta,xi,yi,zi) = 4.5d0*grad_phi_procs(alpha,xi,yi,zi)*grad_phi_procs(beta,xi,yi,zi) &
                            -1.5d0*(grad_phi_procs(1,xi,yi,zi)**2 + grad_phi_procs(2,xi,yi,zi)**2 + grad_phi_procs(3,xi,yi,zi)**2) &
                            *krone(alpha,beta)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine gphi_cal

    subroutine grad_u_cal(grad_u_procs,u1_procs,u2_procs,u3_procs)
        real(8),intent(out):: grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:y_procs)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do beta=1,3
                        grad_u_procs(1,beta,xi,yi,zi) = 0.0d0
                        grad_u_procs(2,beta,xi,yi,zi) = 0.0d0
                        grad_u_procs(3,beta,xi,yi,zi) = 0.0d0
                        do i=2,15
                            grad_u_procs(1,beta,xi,yi,zi) = grad_u_procs(1,beta,xi,yi,zi)+cr(beta,i)*u1_procs(xi+cx(i),yi+cy(i),zi+cz(i))
                            grad_u_procs(2,beta,xi,yi,zi) = grad_u_procs(2,beta,xi,yi,zi)+cr(beta,i)*u2_procs(xi+cx(i),yi+cy(i),zi+cz(i))
                            grad_u_procs(3,beta,xi,yi,zi) = grad_u_procs(3,beta,xi,yi,zi)+cr(beta,i)*u3_procs(xi+cx(i),yi+cy(i),zi+cz(i))
                        enddo
                        grad_u_procs(1,beta,xi,yi,zi) = grad_u_procs(1,beta,xi,yi,zi)/(10.0d0*ds)
                        grad_u_procs(2,beta,xi,yi,zi) = grad_u_procs(2,beta,xi,yi,zi)/(10.0d0*ds)
                        grad_u_procs(3,beta,xi,yi,zi) = grad_u_procs(3,beta,xi,yi,zi)/(10.0d0*ds)
                    enddo
                enddo
            enddo
        enddo
    end subroutine grad_u_cal

    subroutine break_up3(phi,p,grad_phi)
        real(8),intent(in):: p(0:xmax,0:ymax,0:zmax)
        real(8),intent(inout):: phi(-1:xmax+1,-1:ymax+1,-1:zmax+1), grad_phi(1:3,0:xmax,0:ymax,0:zmax)
        yi = ymax/2  !yi = ymax / 2　断面を見たい
        !==========================重心座標を求める===============================================
        xg = 0.0d0
        yg = 0.0d0
        zg = 0.0d0
        phi_tmp = 0.0d0
        do zi=0,zmax
            ! do yi=0,ymax
                do xi=0,xmax
                    xg = xg + dble(xi) * (phi(xi,yi,zi) - phi1)
                    yg = yg + dble(yi) * (phi(xi,yi,zi) - phi1)
                    zg = zg + dble(zi) * (phi(xi,yi,zi) - phi1)
                    phi_tmp = phi_tmp + (phi(xi,yi,zi) - phi1)
                enddo
            ! enddo
        enddo
        xg = xg / phi_tmp
        yg = yg / phi_tmp
        zg = zg / phi_tmp

        !回転角を計算
        distance = 0.0d0
        yi = (ymax+1)/2
        do zi=0,zmax
            do xi=0,xmax
                if(phi(xi,yi,zi) >= 0.5d0*(phi1+phi2)) then
                    if( distance <= sqrt((dble(xi)-xg)**2 + (dble(zi)-zg)**2) ) then
                        distance = sqrt((dble(xi)-xg)**2 + (dble(zi)-zg)**2)
                        x_tmp = xi
                        z_tmp = zi
                    endif
                endif
            enddo
        enddo
        theta = acos(abs((dble(x_tmp) - xg) / distance))
        if(((dble(x_tmp)-xg) > 0.0d0).and.((dble(z_tmp)-zg) > 0.0d0)) then !第一象限
            theta = theta
        else if(((dble(x_tmp)-xg) < 0.0d0).and.((dble(z_tmp)-zg) > 0.0d0)) then !第二象限
            theta = pi - theta
        else if(((dble(x_tmp)-xg) < 0.0d0).and.((dble(z_tmp)-zg) < 0.0d0)) then !第三象限
            theta = pi + theta
        else if(((dble(x_tmp)-xg) > 0.0d0).and.((dble(z_tmp)-zg) < 0.0d0)) then !第四象限
            theta = 2.0d0 * pi - theta 
        endif
        !短辺を計算
        k = 0
        do
            k = k + 1
            zj = (zz*dble(k)) / (tan(theta)*sin(theta) + cos(theta))
            xj = -zj*tan(theta)
            xj = xj + xg
            zj = zj + zg
            x1 = int(xj)
            x2 = int(xj) + 1
            z1 = int(zj)
            z2 = int(zj) + 1
            !線形補間
            phi_hokann = ((dble(x2)-xj)*(dble(z2)-zj)*phi(x1,yi,z1) + (xj-dble(x1))*(dble(z2)-zj)*phi(x2,yi,z1) &
                            + (dble(x2)-xj)*(zj-dble(z1))*phi(x1,yi,z2) + (xj-dble(x1))*(zj-dble(z1))*phi(x2,yi,z2)) &
                            / ((dble(x2)-dble(x1))*(dble(z2)-dble(z1)))
            if(phi_hokann >= 0.5d0*(phi1+phi2)) then
                B1 = sqrt((xj-xg)**2 + (zj-zg)**2) * 2.0d0
            endif
            if(phi_hokann < 0.5d0*(phi1+phi2)) exit
        enddo
        !長辺を計算
        k = 0
        do
            k = k + 1
            xj = (xx*dble(k)) / (tan(theta)*sin(theta) + cos(theta))
            zj = xj*tan(theta)
            xj = xj + xg
            zj = zj + zg
            x1 = int(xj)
            x2 = int(xj) + 1
            z1 = int(zj)
            z2 = int(zj) + 1
            !線形補間
            phi_hokann = ((dble(x2)-xj)*(dble(z2)-zj)*phi(x1,yi,z1) + (xj-dble(x1))*(dble(z2)-zj)*phi(x2,yi,z1) &
                            + (dble(x2)-xj)*(zj-dble(z1))*phi(x1,yi,z2) + (xj-dble(x1))*(zj-dble(z1))*phi(x2,yi,z2)) &
                            / ((dble(x2)-dble(x1))*(dble(z2)-dble(z1)))
            if(phi_hokann >= 0.5d0*(phi1+phi2)) then
                L1 = sqrt((xj-xg)**2 + (zj-zg)**2) * 2.0d0
            endif
            if(phi_hokann < 0.5d0*(phi1+phi2)) exit
        enddo

        !変形度
        Dd = (L1-B1) / (L1+B1)

        !=============================界面張力===============================================
        if(n <= start) then
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
            !表面張力の計算
            yi = (ymax+1) / 2
            zi = (zmax+1) / 2
            sigma = p(0,0,0) - p(1,1,1)
            sigma = 0.0d0
            do xi=0,(xmax)/2
                sigma = sigma + grad_phi(1,xi,yi,zi)**2
            enddo
            sigma = sigma * kappag_tmp
            if(n == 10) then
                open(10,file="./stable5.d")
                write(10,"(4es16.8)") dble(n),L1,sigma,kappag_tmp
                close(10)
            else
                open(10,file="./stable5.d",action="write",position="append")
                write(10,"(4es16.8)") dble(n),L1,sigma,kappag_tmp
                close(10)
            endif
            
        endif

        ! write(9,"(5es16.8)") dble(n),(dble(n)-dble(start))*2.0d0*uw/H1,L1,B1,Dd
    end subroutine break_up3

    subroutine output(phi,u1,u2,u3,p)
        real(8),intent(inout):: phi(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: p(0:xmax,0:ymax,0:zmax)
        real(8),intent(inout):: u1(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: u2(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: u3(-1:xmax+1,-1:ymax+1,-1:zmax+1)

        write(filename,*) n !i->filename 変換
        filename=datadir//trim(adjustl(filename))//'.bin' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
        print *, filename !表示してみる
        open(20,file=filename, form = "unformatted",status='replace') 

        ! yi = (ymax+1) /2 
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    write(20) dble(xi),dble(yi),dble(zi),phi(xi,yi,zi),u1(xi,yi,zi),u2(xi,yi,zi),u3(xi,yi,zi),p(xi,yi,zi)
                enddo
            enddo
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
    real(8) feq(1:15,0:xmax,0:ymax,0:zmax) !平衡分布関数（オーダーパラメータ）
    real(8) geq(1:15,0:xmax,0:ymax,0:zmax) !平衡分布関数（速度・圧力計算用）
    real(8) phi(-1:xmax+1,-1:ymax+1,-1:zmax+1) !秩序パラメータ
    real(8) u1(-1:xmax+1,-1:ymax+1,-1:zmax+1), u2(-1:xmax+1,-1:ymax+1,-1:zmax+1), u3(-1:xmax+1,-1:ymax+1,-1:zmax+1) !流速
    real(8) p(0:xmax,0:ymax,0:zmax) !圧力
    real(8) U(1:3) !壁の移動速度
    real(8) grad_phi(1:3,0:xmax,0:ymax,0:zmax)
    real(8) lap_phi(0:xmax,0:ymax,0:zmax)
    real(8) gphi(1:3,1:3,0:xmax,0:ymax,0:zmax)
    real(8) p0(0:xmax,0:ymax,0:zmax)
    real(8) grad_u(1:3,1:3,0:xmax,0:ymax,0:zmax)
    ! real(8) nu(-1:xmax+1,-1:ymax+1,-1:zmax+1)

    real(8),allocatable :: f_procs(:,:,:,:), fnext_procs(:,:,:,:), feq_procs(:,:,:,:)
    real(8),allocatable :: g_procs(:,:,:,:), gnext_procs(:,:,:,:), geq_procs(:,:,:,:)
    real(8),allocatable :: phi_procs(:,:,:), p_procs(:,:,:)
    real(8),allocatable :: u1_procs(:,:,:), u2_procs(:,:,:), u3_procs(:,:,:)

    real(8),allocatable :: taug_procs(:,:,:), nu_procs(:,:,:)
    real(8),allocatable :: p0_procs(:,:,:), lap_phi_procs(:,:,:), grad_phi_procs(:,:,:,:), gphi_procs(:,:,:,:,:), grad_u_procs(:,:,:,:,:)
    
    call cpu_time(time1)

    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, comm_procs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, comm_rank, ierr)

    call mk_dirs(datadir)
    call par(cx,cy,cz,cr,krone,U)
    call ini(phi,u1,u2,u3,p,geq,g_procs,gnext_procs,feq,f_procs,fnext_procs,feq_procs,geq_procs,phi_procs,p_procs,u1_procs,u2_procs,u3_procs,grad_phi,gphi,lap_phi,p0,grad_u)
    call ini_op(taug_procs,nu_procs,p0_procs,lap_phi_procs,grad_phi_procs,gphi_procs,grad_u_procs,phi_procs)

    
    ! open(9,file="break2.d")

DO n=1,step
    kappag_tmp = kappag
!======================のりしろ境界==========================================================
    call MPI_boundary(phi_procs,u1_procs,u2_procs,u3_procs,f_procs,g_procs,taug_procs)
!=======================div/gradの計算============================================================
    call lap_cal(lap_phi_procs,phi_procs)
    call grad_cal(grad_phi_procs,phi_procs)
    call p0_cal(p0_procs,phi_procs)
    call gphi_cal(gphi_procs,grad_phi_procs)
    call grad_u_cal(grad_u_procs,u1_procs,u2_procs,u3_procs)
!=======================局所平衡分布関数の計算========================================================
    call feq_cal(gphi_procs,phi_procs,p0_procs,lap_phi_procs,grad_phi_procs,u1_procs,u2_procs,u3_procs,feq_procs)
    call geq_cal(gphi_procs,p_procs,u1_procs,u2_procs,u3_procs,geq_procs,cr,grad_u_procs)
    call MPI_boundary_fg(feq_procs,geq_procs)
!=======================平衡分布関数gの計算=============================================================
    call f_cal(fnext_procs,f_procs,feq_procs)
    call g_cal(gnext_procs,g_procs,geq_procs,taug_procs)
!=======================境界条件=============================================================
    call bounce_back_LBM(gnext_procs,U)
    call renew(f_procs,fnext_procs)
    call renew(g_procs,gnext_procs)
!==========================物理量の計算==========================================================
    call physics(phi_procs,p_procs,u1_procs,u2_procs,u3_procs,f_procs,g_procs,nu_procs,taug_procs)
!============================物理量をまとめる======================================================
if ( mod(n,10)==0 ) then
    !p
    tmp(:,:,:) = 0.0d0
    do zi=1,y_procs
        do yi=1,ymax+1
            do xi=1,x_procs
                tmp(xi,yi,zi) = p_procs(xi,yi,zi)
            enddo
        enddo
    enddo
    if(comm_rank == 0) then
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp4(0,xi,yi,zi) = p_procs(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 1) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,1,MPI_COMM_WORLD,req1s,ierr)
        call MPI_Wait(req1s,sta1s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,1,1,MPI_COMM_WORLD,req1r,ierr)
        call MPI_Wait(req1r,sta1r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp4(1,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 2) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,2,MPI_COMM_WORLD,req2s,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,2,2,MPI_COMM_WORLD,req2r,ierr)
        call MPI_Wait(req2r,sta2r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp4(2,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 3) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,3,MPI_COMM_WORLD,req3s,ierr)
        call MPI_Wait(req3s,sta3s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,3,3,MPI_COMM_WORLD,req3r,ierr)
        call MPI_Wait(req3r,sta3r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp4(3,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 4) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,4,MPI_COMM_WORLD,req4s,ierr)
        call MPI_Wait(req4s,sta4s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,4,4,MPI_COMM_WORLD,req4r,ierr)
        call MPI_Wait(req4r,sta4r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp4(4,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 5) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,5,MPI_COMM_WORLD,req5s,ierr)
        call MPI_Wait(req5s,sta5s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,5,5,MPI_COMM_WORLD,req5r,ierr)
        call MPI_Wait(req5r,sta5r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp4(5,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 6) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,6,MPI_COMM_WORLD,req6s,ierr)
        call MPI_Wait(req6s,sta6s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,6,6,MPI_COMM_WORLD,req6r,ierr)
        call MPI_Wait(req6r,sta6r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp4(6,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 7) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,7,MPI_COMM_WORLD,req7s,ierr)
        call MPI_Wait(req7s,sta7s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,7,7,MPI_COMM_WORLD,req7r,ierr)
        call MPI_Wait(req7r,sta7r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp4(7,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 8) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,8,MPI_COMM_WORLD,req8s,ierr)
        call MPI_Wait(req8s,sta8s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,8,8,MPI_COMM_WORLD,req8r,ierr)
        call MPI_Wait(req8r,sta8r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp4(8,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif


    !phi
    tmp(:,:,:) = 0.0d0
    do zi=1,y_procs
        do yi=1,ymax+1
            do xi=1,x_procs
                tmp(xi,yi,zi) = phi_procs(xi,yi,zi)
            enddo
        enddo
    enddo
    if(comm_rank == 0) then
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp5(0,xi,yi,zi) = phi_procs(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 1) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,1,MPI_COMM_WORLD,req1s,ierr)
        call MPI_Wait(req1s,sta1s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,1,1,MPI_COMM_WORLD,req1r,ierr)
        call MPI_Wait(req1r,sta1r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp5(1,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 2) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,2,MPI_COMM_WORLD,req2s,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,2,2,MPI_COMM_WORLD,req2r,ierr)
        call MPI_Wait(req2r,sta2r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp5(2,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 3) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,3,MPI_COMM_WORLD,req3s,ierr)
        call MPI_Wait(req3s,sta3s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,3,3,MPI_COMM_WORLD,req3r,ierr)
        call MPI_Wait(req3r,sta3r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp5(3,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 4) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,4,MPI_COMM_WORLD,req4s,ierr)
        call MPI_Wait(req4s,sta4s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,4,4,MPI_COMM_WORLD,req4r,ierr)
        call MPI_Wait(req4r,sta4r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp5(4,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 5) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,5,MPI_COMM_WORLD,req5s,ierr)
        call MPI_Wait(req5s,sta5s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,5,5,MPI_COMM_WORLD,req5r,ierr)
        call MPI_Wait(req5r,sta5r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp5(5,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 6) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,6,MPI_COMM_WORLD,req6s,ierr)
        call MPI_Wait(req6s,sta6s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,6,6,MPI_COMM_WORLD,req6r,ierr)
        call MPI_Wait(req6r,sta6r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp5(6,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 7) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,7,MPI_COMM_WORLD,req7s,ierr)
        call MPI_Wait(req7s,sta7s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,7,7,MPI_COMM_WORLD,req7r,ierr)
        call MPI_Wait(req7r,sta7r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp5(7,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 8) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,8,MPI_COMM_WORLD,req8s,ierr)
        call MPI_Wait(req8s,sta8s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,8,8,MPI_COMM_WORLD,req8r,ierr)
        call MPI_Wait(req8r,sta8r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp5(8,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    ! !u1
    ! tmp(:,:,:) = 0.0d0
    ! do zi=1,y_procs
    !     do yi=1,ymax+1
    !         do xi=1,x_procs
    !             tmp(xi,yi,zi) = u1_procs(xi,yi,zi)
    !         enddo
    !     enddo
    ! enddo
    ! if(comm_rank == 0) then
    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp1(0,xi,yi,zi) = u1_procs(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 1) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,1,MPI_COMM_WORLD,req1s,ierr)
    !     call MPI_Wait(req1s,sta1s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,1,1,MPI_COMM_WORLD,req1r,ierr)
    !     call MPI_Wait(req1r,sta1r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp1(1,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 2) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,2,MPI_COMM_WORLD,req2s,ierr)
    !     call MPI_Wait(req2s,sta2s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,2,2,MPI_COMM_WORLD,req2r,ierr)
    !     call MPI_Wait(req2r,sta2r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp1(2,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 3) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,3,MPI_COMM_WORLD,req3s,ierr)
    !     call MPI_Wait(req3s,sta3s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,3,3,MPI_COMM_WORLD,req3r,ierr)
    !     call MPI_Wait(req3r,sta3r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp1(3,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 4) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,4,MPI_COMM_WORLD,req4s,ierr)
    !     call MPI_Wait(req4s,sta4s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,4,4,MPI_COMM_WORLD,req4r,ierr)
    !     call MPI_Wait(req4r,sta4r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp1(4,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 5) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,5,MPI_COMM_WORLD,req5s,ierr)
    !     call MPI_Wait(req5s,sta5s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,5,5,MPI_COMM_WORLD,req5r,ierr)
    !     call MPI_Wait(req5r,sta5r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp1(5,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 6) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,6,MPI_COMM_WORLD,req6s,ierr)
    !     call MPI_Wait(req6s,sta6s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,6,6,MPI_COMM_WORLD,req6r,ierr)
    !     call MPI_Wait(req6r,sta6r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp1(6,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 7) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,7,MPI_COMM_WORLD,req7s,ierr)
    !     call MPI_Wait(req7s,sta7s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,7,7,MPI_COMM_WORLD,req7r,ierr)
    !     call MPI_Wait(req7r,sta7r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp1(7,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 8) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,8,MPI_COMM_WORLD,req8s,ierr)
    !     call MPI_Wait(req8s,sta8s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,8,8,MPI_COMM_WORLD,req8r,ierr)
    !     call MPI_Wait(req8r,sta8r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp1(8,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    
    ! !u2
    ! tmp(:,:,:) = 0.0d0
    ! do zi=1,y_procs
    !     do yi=1,ymax+1
    !         do xi=1,x_procs
    !             tmp(xi,yi,zi) = u2_procs(xi,yi,zi)
    !         enddo
    !     enddo
    ! enddo
    ! if(comm_rank == 0) then
    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp2(0,xi,yi,zi) = u2_procs(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 1) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,1,MPI_COMM_WORLD,req1s,ierr)
    !     call MPI_Wait(req1s,sta1s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,1,1,MPI_COMM_WORLD,req1r,ierr)
    !     call MPI_Wait(req1r,sta1r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp2(1,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 2) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,2,MPI_COMM_WORLD,req2s,ierr)
    !     call MPI_Wait(req2s,sta2s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,2,2,MPI_COMM_WORLD,req2r,ierr)
    !     call MPI_Wait(req2r,sta2r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp2(2,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 3) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,3,MPI_COMM_WORLD,req3s,ierr)
    !     call MPI_Wait(req3s,sta3s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,3,3,MPI_COMM_WORLD,req3r,ierr)
    !     call MPI_Wait(req3r,sta3r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp2(3,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 4) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,4,MPI_COMM_WORLD,req4s,ierr)
    !     call MPI_Wait(req4s,sta4s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,4,4,MPI_COMM_WORLD,req4r,ierr)
    !     call MPI_Wait(req4r,sta4r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp2(4,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 5) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,5,MPI_COMM_WORLD,req5s,ierr)
    !     call MPI_Wait(req5s,sta5s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,5,5,MPI_COMM_WORLD,req5r,ierr)
    !     call MPI_Wait(req5r,sta5r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp2(5,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 6) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,6,MPI_COMM_WORLD,req6s,ierr)
    !     call MPI_Wait(req6s,sta6s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,6,6,MPI_COMM_WORLD,req6r,ierr)
    !     call MPI_Wait(req6r,sta6r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp2(6,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 7) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,7,MPI_COMM_WORLD,req7s,ierr)
    !     call MPI_Wait(req7s,sta7s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,7,7,MPI_COMM_WORLD,req7r,ierr)
    !     call MPI_Wait(req7r,sta7r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp2(7,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 8) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,8,MPI_COMM_WORLD,req8s,ierr)
    !     call MPI_Wait(req8s,sta8s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,8,8,MPI_COMM_WORLD,req8r,ierr)
    !     call MPI_Wait(req8r,sta8r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp2(8,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    
    ! !u3
    ! tmp(:,:,:) = 0.0d0
    ! do zi=1,y_procs
    !     do yi=1,ymax+1
    !         do xi=1,x_procs
    !             tmp(xi,yi,zi) = u3_procs(xi,yi,zi)
    !         enddo
    !     enddo
    ! enddo
    ! if(comm_rank == 0) then
    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp3(0,xi,yi,zi) = u3_procs(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 1) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,1,MPI_COMM_WORLD,req1s,ierr)
    !     call MPI_Wait(req1s,sta1s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,1,1,MPI_COMM_WORLD,req1r,ierr)
    !     call MPI_Wait(req1r,sta1r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp3(1,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 2) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,2,MPI_COMM_WORLD,req2s,ierr)
    !     call MPI_Wait(req2s,sta2s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,2,2,MPI_COMM_WORLD,req2r,ierr)
    !     call MPI_Wait(req2r,sta2r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp3(2,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 3) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,3,MPI_COMM_WORLD,req3s,ierr)
    !     call MPI_Wait(req3s,sta3s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,3,3,MPI_COMM_WORLD,req3r,ierr)
    !     call MPI_Wait(req3r,sta3r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp3(3,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 4) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,4,MPI_COMM_WORLD,req4s,ierr)
    !     call MPI_Wait(req4s,sta4s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,4,4,MPI_COMM_WORLD,req4r,ierr)
    !     call MPI_Wait(req4r,sta4r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp3(4,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 5) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,5,MPI_COMM_WORLD,req5s,ierr)
    !     call MPI_Wait(req5s,sta5s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,5,5,MPI_COMM_WORLD,req5r,ierr)
    !     call MPI_Wait(req5r,sta5r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp3(5,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 6) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,6,MPI_COMM_WORLD,req6s,ierr)
    !     call MPI_Wait(req6s,sta6s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,6,6,MPI_COMM_WORLD,req6r,ierr)
    !     call MPI_Wait(req6r,sta6r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp3(6,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 7) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,7,MPI_COMM_WORLD,req7s,ierr)
    !     call MPI_Wait(req7s,sta7s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,7,7,MPI_COMM_WORLD,req7r,ierr)
    !     call MPI_Wait(req7r,sta7r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp3(7,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 8) then
    !     call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,8,MPI_COMM_WORLD,req8s,ierr)
    !     call MPI_Wait(req8s,sta8s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,8,8,MPI_COMM_WORLD,req8r,ierr)
    !     call MPI_Wait(req8r,sta8r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 tmp3(8,xi,yi,zi) = tmp(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    ! endif


    !u1,u2,u3,p,phiをまとめるぞ！
    if(comm_rank == 0) then
        k = 0
        do Nyy=0,Ny-1
            do Nxx=0,Nx-1
                do zi=1,y_procs
                    do yi=1,ymax+1
                        do xi=1,x_procs
                            ! u1((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs)) = tmp1(k,xi,yi,zi)
                            ! u2((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs)) = tmp2(k,xi,yi,zi)
                            ! u3((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs)) = tmp3(k,xi,yi,zi)
                            p((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs)) = tmp4(k,xi,yi,zi)
                            phi((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs)) = tmp5(k,xi,yi,zi)
                        enddo
                    enddo
                enddo
                k = k + 1
            enddo
        enddo
    endif
    if(comm_rank == 0) then
        call break_up3(phi,p,grad_phi)
        if ( mod(n,30000)==0 ) then
            call output(phi,u1,u2,u3,p)
        endif
    endif
endif
!============================================================================================
write(*,*) "time = ", n
ENDDO
    ! close(10)
    ! close(9)

    call MPI_Finalize(ierr)
    call cpu_time(time2)
end program main