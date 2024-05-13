! 【Multi-times simulations】
! このプログラムは、多数回の液滴分裂の数値シミュレーションを実行する。
! このプログラムは、界面を「Cahn-Hilliard 方程式」で表現している。
! 液滴を長時間保持できないが、系全体の保存性には優れている。
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module globals
    include "mpif.h"
    real(8),parameter:: ds = 1.0d0 !格子間隔（lattice unit）
    integer,parameter:: xmax = 299 !ｘ方向格子数（０から数える）
    integer,parameter:: ymax = 128 !ｙ方向格子数（０から数える）
    integer,parameter:: zmax = 128 !ｚ方向格子数（０から数える）
    integer,parameter:: Nxall = 1 !x方向の分割数（Nxall*Nzall=全体の計算数）
    integer,parameter:: Nzall = 4 !z方向の分割数（Nxall*Nzall=全体の計算数）
    integer,parameter:: xall = (xmax + 1) * Nxall !全体のx方向格子数
    integer,parameter:: zall = (zmax + 1) * Nzall !全体のz方向格子数
    integer,parameter:: step = 200000 !計算時間step
    integer,parameter:: start = 5000 !壁を動かし始める時間step
    integer,parameter:: dirstep = 100 !ディレクトリーを作成するstep
    integer i, j, k, n, nall, nstar, xi, yi, zi, alpha, beta
    character :: filename*200
    character :: filename2*200
    ! character(*),parameter :: datadir = "/data/group1/z40137r/rechange_caconst_eta1/"
    character(*),parameter :: datadir = "/data/n/n517/crit/re150_2/"

    !支配パラメータ
    real(8) D !設置する液滴径
    real(8) uw !壁の移動速度
    real(8),parameter:: We = 2.0d0 !ウエーバー数
    real(8) nu1  !液体の動粘度
    real(8) nu2  !液滴の動粘度
    real(8),parameter:: kappaf = 0.01d0*ds**2 !界面厚さを決めるパラメータ
    real(8) kappag  !界面張力を決めるパラメータ
    ! real(8),parameter:: nu1 = 0.8d0 !液体の動粘度
    ! real(8),parameter:: nu2 = 0.08d0 !液滴の動粘度
    ! real(8),parameter:: kappaf = 0.01d0*ds**2 !界面厚さを決めるパラメータ
    ! real(8),parameter:: kappag = 9.2d-3 !界面張力を決めるパラメータ
    real(8),parameter:: H1 = dble(zmax) !代表長さ
    real(8),parameter:: phi1 = 2.211d0
    real(8),parameter:: phi2 = 4.895d0
    real(8),parameter:: a = 9.0d0/49.0d0
    real(8),parameter:: b = 2.0d0/21.0d0
    real(8),parameter:: T = 0.55d0
    real(8),parameter:: tauf = 0.7d0
    real(8),parameter:: Anu = 0.0d0

    !初期液滴中心
    real(8), parameter:: xc = 0.5d0*dble(xmax)	
    real(8), parameter:: yc = 0.5d0*dble(ymax)
    real(8), parameter:: zc = 0.5d0*dble(zmax)

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

    !その他変数
    real(8) phi_min, phi_max, min, max
    real(8) gtemp, ftemp

    !MPI用変数
    integer ierr, comm_procs, comm_rank
    integer Nx, Nz !xi、yi方向の分割数
    integer x_procs, z_procs
    integer key_new, group_new, key_x, group_x, key_z, group_z
    integer new_comm_world, new_procs, new_rank
    integer newx_comm_world, newx_procs, newx_rank
    integer newz_comm_world, newz_procs, newz_rank
    integer rectangle_type,rectangle_type_f,xtype,xtype_f,output_type
    integer next_rank_x,former_rank_x,next_rank_z,former_rank_z
    integer req1s,req1r,req2s,req2r,req3s,req3r,req4s,req4r,req5s,req5r,req6s,req6r,req7s,req7r,req8s,req8r,req9s,req9r,req10s,req10r
    integer, dimension(MPI_STATUS_SIZE) :: sta1s,sta1r,sta2s,sta2r,sta3s,sta3r,sta4s,sta4r,sta5s,sta5r,sta6s,sta6r,sta7s,sta7r,sta8s,sta8r,sta9s,sta9r,sta10s,sta10r
    integer Nzz,Nxx
    real(8), allocatable:: tmp1(:,:,:,:),tmpf(:,:,:,:,:),tmp(:,:,:)
    integer tags, tagr, recv_rank
    character(2) chmyrank


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
        real(8),intent(inout):: p(-1:xmax+1,-1:ymax+1,-1:zmax+1)
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

    !========================並列数・コミュニケータ分割・通信先設定========================
        !配列のallocate(xi:0とx_procs+1がのりしろ)(yi:0とymax+2がのりしろ)(zi:0とz_procs+1がのりしろ)
        Nx = 6 !x方向の並列数（ただし，Nx/=comm_procs）
        Nz = comm_procs / (Nx * Nxall * Nzall) !z方向の並列数
        x_procs = (xmax+1) / Nx
        z_procs = (zmax+1) / Nz

        !各条件でコミュニケータを分割する
        key_new = comm_rank
        group_new = comm_rank / (Nx * Nz)
        call MPI_Comm_Split(MPI_COMM_WORLD,group_new,key_new,new_comm_world,ierr)
        call MPI_Comm_Size(new_comm_world,new_procs,ierr)
        call MPI_Comm_Rank(new_comm_world,new_rank,ierr)
        !z方向にコミュニケータを分割する
        key_z = new_rank
        group_z = new_rank / Nx
        call MPI_Comm_Split(new_comm_world,group_z,key_z,newz_comm_world,ierr)
        call MPI_Comm_Size(newz_comm_world,newz_procs,ierr)
        call MPI_Comm_Rank(newz_comm_world,newz_rank,ierr)

        !x方向にコミュニケータを分割する
        key_x = new_rank
        group_x = mod(new_rank,Nx)
        call MPI_Comm_Split(new_comm_world,group_x,key_x,newx_comm_world,ierr)
        call MPI_Comm_Size(newx_comm_world,newx_procs,ierr)
        call MPI_Comm_Rank(newx_comm_world,newx_rank,ierr)

        !のりしろ境界の通信先設定
        next_rank_x = newx_rank + 1
        former_rank_x = newx_rank - 1
        if(newx_rank == 0) then
            former_rank_x = Nz - 1
        else if(newx_rank == Nz - 1) then
            next_rank_x = 0
        endif
        next_rank_z = newz_rank + 1
        former_rank_z = newz_rank - 1
        if(newz_rank == 0) then
            former_rank_z = Nx - 1
        else if(newz_rank == Nx - 1) then
            next_rank_z = 0
        endif

        !x世界での受け渡しをする際の型作成
        call MPI_Type_Vector(ymax+3,x_procs,x_procs+2,MPI_REAL8,xtype,ierr)
        call MPI_Type_Commit(xtype,ierr)
        call MPI_Type_Vector(ymax+3,15*x_procs,15*(x_procs+2),MPI_REAL8,xtype_f,ierr)
        call MPI_Type_Commit(xtype_f,ierr)

        !z世界での受け渡しをする際の型作成
        call MPI_Type_Vector((z_procs+2)*(ymax+3),1,x_procs+2,MPI_REAL8,rectangle_type,ierr)
        call MPI_Type_Commit(rectangle_type,ierr)
        call MPI_Type_Vector((z_procs+2)*(ymax+3),15,15*(x_procs+2),MPI_REAL8,rectangle_type_f,ierr)
        call MPI_Type_Commit(rectangle_type_f,ierr)
    !===============================================================================================================
        !以下はのりしろ有りの変数
        allocate(phi_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(p_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(f_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(fnext_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(feq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(g_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(gnext_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(geq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(tmp1(0:new_procs-1,1:x_procs,1:ymax+1,1:z_procs))
        ! allocate(tmp2(0:new_procs-1,1:x_procs,1:ymax+1,1:z_procs))
        ! allocate(tmp3(0:new_procs-1,1:x_procs,1:ymax+1,1:z_procs))
        ! allocate(tmp4(0:new_procs-1,1:x_procs,1:ymax+1,1:z_procs))
        ! allocate(tmp5(0:new_procs-1,1:x_procs,1:ymax+1,1:z_procs))
        allocate(tmpf(0:new_procs-1,1:15,1:x_procs,1:ymax+1,1:z_procs))
        ! allocate(tmpff(0:new_procs-1,1:15,1:x_procs,1:ymax+1,1:z_procs))
        allocate(tmp(1:x_procs,1:ymax+1,1:z_procs))
        ! allocate(tmp_nu(0:new_procs-1,1:x_procs,1:ymax+1,1:z_procs))
        
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

        if(group_new == 0) then
            uw = 0.07d0
            nu1 = 0.00746d0
            nu2 = 0.00746d0
            kappag = 3.02d-3
            !液滴径をランダムに発生させる（後でやる）
            D = 64.0d0
        elseif(group_new == 1) then
            uw = 0.07d0
            nu1 = 0.00746d0
            nu2 = 0.00746d0
            kappag = 2.52d-3
            !液滴径をランダムに発生させる（後でやる）
            D = 64.0d0
        elseif(group_new == 2) then
            uw = 0.07d0
            nu1 = 0.00746d0
            nu2 = 0.00746d0
            kappag = 2.16d-3
            !液滴径をランダムに発生させる（後でやる）
            D = 64.0d0
        elseif(group_new == 3) then
            uw = 0.07d0
            nu1 = 0.00746d0
            nu2 = 0.00746d0
            kappag = 1.89d-3
            !液滴径をランダムに発生させる（後でやる）
            D = 64.0d0
        ! elseif(group_new == 4) then
        !     uw = 0.07d0
        !     nu1 = 0.00746d0
        !     nu2 = 0.00746d0
        !     kappag = 1.89d-3
        !     !液滴径をランダムに発生させる（後でやる）
        !     D = 64.0d0
        ! elseif(group_new == 5) then
        !     uw = 0.07d0
        !     nu1 = 0.00746d0
        !     nu2 = 0.00746d0
        !     kappag = 1.51d-3
        !     !液滴径をランダムに発生させる（後でやる）
        !     D = 64.0d0
        endif
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

        call tmpfunction(u1, u1_procs)
        call tmpfunction(u2, u2_procs)
        call tmpfunction(u3, u3_procs)
        call tmpfunction(p, p_procs)
        call tmpfunction(phi, phi_procs)
        call tmpfunction2(geq, g_procs)
        call tmpfunction2(feq, f_procs)
        ! nstar = 0
        ! do Nzz=0,Nz-1
        !     do Nxx=0,Nx-1
        !         do zi=1,z_procs
        !             do yi=1,ymax+1
        !                 do xi=1,x_procs
        !                     tmp1(nstar,xi,yi,zi) = u1((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*(z_procs))
        !                     ! tmp2(nstar,xi,yi,zi) = u2((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*(z_procs))
        !                     ! tmp3(nstar,xi,yi,zi) = u3((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*(z_procs))
        !                     ! tmp4(nstar,xi,yi,zi) = p((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*(z_procs))
        !                     ! tmp5(nstar,xi,yi,zi) = phi((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*(z_procs))
        !                     ! do i=1,15
        !                     !     tmpf(nstar,i,xi,yi,zi) = geq(i,(xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*(z_procs))
        !                     !     tmpff(nstar,i,xi,yi,zi) = feq(i,(xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*(z_procs))
        !                     ! enddo
        !                 enddo
        !             enddo
        !         enddo
        !         nstar = nstar + 1
        !     enddo
        ! enddo

        ! do n = 0, new_procs-1
        !     if(n == new_rank) then
        !         do zi=1,z_procs
        !             do yi=1,ymax+1
        !                 do xi=1,x_procs
        !                     u1_procs(xi,yi,zi) = tmp1(n,xi,yi,zi)
        !                     ! u2_procs(xi,yi,zi) = tmp2(n,xi,yi,zi)
        !                     ! u3_procs(xi,yi,zi) = tmp3(n,xi,yi,zi)
        !                     ! p_procs(xi,yi,zi) = tmp4(n,xi,yi,zi)
        !                     ! phi_procs(xi,yi,zi) = tmp5(n,xi,yi,zi)
        !                     ! do i=1,15
        !                     !     g_procs(i,xi,yi,zi) = tmpf(n,i,xi,yi,zi)
        !                     !     f_procs(i,xi,yi,zi) = tmpff(n,i,xi,yi,zi)
        !                     ! enddo
        !                 enddo
        !             enddo
        !         enddo
        !     endif
        ! enddo
    end subroutine ini

    subroutine tmpfunction(fun,fun_procs)
        real(8),intent(inout):: fun(-1:xmax+1,-1:ymax+1,-1:zmax+1)
        real(8),intent(inout):: fun_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)

        nstar = 0
        do Nzz=0,Nz-1
            do Nxx=0,Nx-1
                do zi=1,z_procs
                    do yi=1,ymax+1
                        do xi=1,x_procs
                            tmp1(nstar,xi,yi,zi) = fun((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*(z_procs))
                        enddo
                    enddo
                enddo
                nstar = nstar + 1
            enddo
        enddo

        do n = 0, new_procs-1
            if(n == new_rank) then
                do zi=1,z_procs
                    do yi=1,ymax+1
                        do xi=1,x_procs
                            fun_procs(xi,yi,zi) = tmp1(n,xi,yi,zi)
                        enddo
                    enddo
                enddo
            endif
        enddo
    endsubroutine tmpfunction

    subroutine tmpfunction2(fun,fun_procs)
        real(8),intent(inout):: fun(1:15,0:xmax,0:ymax,0:zmax)
        real(8),intent(inout):: fun_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)

        nstar = 0
        do Nzz=0,Nz-1
            do Nxx=0,Nx-1
                do zi=1,z_procs
                    do yi=1,ymax+1
                        do xi=1,x_procs
                            do i=1,15
                                tmpf(nstar,i,xi,yi,zi) = fun(i,(xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*(z_procs))
                            enddo
                        enddo
                    enddo
                enddo
                nstar = nstar + 1
            enddo
        enddo

        do n = 0, new_procs-1
            if(n == new_rank) then
                do zi=1,z_procs
                    do yi=1,ymax+1
                        do xi=1,x_procs
                            do i=1,15
                                fun_procs(i,xi,yi,zi) = tmpf(n,i,xi,yi,zi)
                            enddo
                        enddo
                    enddo
                enddo
            endif
        enddo
    endsubroutine tmpfunction2


    subroutine ini_op(taug_procs,nu_procs,p0_procs,lap_phi_procs,grad_phi_procs,gphi_procs,grad_u_procs,phi_procs)
        real(8),intent(in):: phi_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),allocatable :: taug_procs(:,:,:), nu_procs(:,:,:)
        real(8),allocatable :: p0_procs(:,:,:), lap_phi_procs(:,:,:), grad_phi_procs(:,:,:,:), gphi_procs(:,:,:,:,:), grad_u_procs(:,:,:,:,:)
        !以下はのりしろ無しの変数
        allocate(taug_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(nu_procs(1:x_procs,1:ymax+1,1:z_procs))
        allocate(p0_procs(1:x_procs,1:ymax+1,1:z_procs))
        allocate(lap_phi_procs(1:x_procs,1:ymax+1,1:z_procs))
        allocate(grad_phi_procs(1:3,1:x_procs,1:ymax+1,1:z_procs))
        allocate(gphi_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs))
        allocate(grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs))

        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    nu_procs(xi,yi,zi) = (phi_procs(xi,yi,zi)-phi1)/(phi2-phi1)*(nu2-nu1) + nu1
                    taug_procs(xi,yi,zi) = 0.5d0 + 3.0d0*nu_procs(xi,yi,zi)/ds + 2.0d0/3.0d0*Anu
                enddo
            enddo
        enddo
    end subroutine ini_op

    subroutine glue(var)
        real(8),intent(inout) :: var(0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
            do xi=1,x_procs
                var(xi,0,zi) = var(xi,ymax+1,zi)
                var(xi,ymax+2,zi) = var(xi,1,zi)
            enddo
        enddo
        !X世界でののりしろ通信
        call MPI_Isend(var(1,0,z_procs),1,xtype,next_rank_x,1,newx_comm_world,req1s,ierr)
        call MPI_Irecv(var(1,0,0),1,xtype,former_rank_x,1,newx_comm_world,req1r,ierr)

        call MPI_Isend(var(1,0,1),1,xtype,former_rank_x,2,newx_comm_world,req2s,ierr)
        call MPI_Irecv(var(1,0,z_procs+1),1,xtype,next_rank_x,2,newx_comm_world,req2r,ierr)

        call MPI_Wait(req1s,sta1s,ierr)
        call MPI_Wait(req1r,sta1r,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
        call MPI_Wait(req2r,sta2r,ierr)
        !Y世界でののりしろ通信
        call MPI_Isend(var(x_procs,0,0),1,rectangle_type,next_rank_z,1,newz_comm_world,req1s,ierr)
        call MPI_Irecv(var(0,0,0),1,rectangle_type,former_rank_z,1,newz_comm_world,req1r,ierr)

        call MPI_Isend(var(1,0,0),1,rectangle_type,former_rank_z,2,newz_comm_world,req2s,ierr)
        call MPI_Irecv(var(x_procs+1,0,0),1,rectangle_type,next_rank_z,2,newz_comm_world,req2r,ierr)

        call MPI_Wait(req1s,sta1s,ierr)
        call MPI_Wait(req1r,sta1r,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
        call MPI_Wait(req2r,sta2r,ierr)
    end subroutine glue

    subroutine glue_f(var_f)
        real(8),intent(inout) :: var_f(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
            do xi=1,x_procs
                do i=1,15
                    var_f(i,xi,0,zi) = var_f(i,xi,ymax+1,zi)
                    var_f(i,xi,ymax+2,zi) = var_f(i,xi,1,zi)
                enddo
            enddo
        enddo
        !X世界でののりしろ通信
        call MPI_Isend(var_f(1,1,0,z_procs),1,xtype_f,next_rank_x,1,newx_comm_world,req1s,ierr)
        call MPI_Irecv(var_f(1,1,0,0),1,xtype_f,former_rank_x,1,newx_comm_world,req1r,ierr)

        call MPI_Isend(var_f(1,1,0,1),1,xtype_f,former_rank_x,2,newx_comm_world,req2s,ierr)
        call MPI_Irecv(var_f(1,1,0,z_procs+1),1,xtype_f,next_rank_x,2,newx_comm_world,req2r,ierr)

        call MPI_Wait(req1s,sta1s,ierr)
        call MPI_Wait(req1r,sta1r,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
        call MPI_Wait(req2r,sta2r,ierr)

        !Y世界でののりしろ通信
        call MPI_Isend(var_f(1,x_procs,0,0),1,rectangle_type_f,next_rank_z,1,newz_comm_world,req1s,ierr)
        call MPI_Irecv(var_f(1,0,0,0),1,rectangle_type_f,former_rank_z,1,newz_comm_world,req1r,ierr)

        call MPI_Isend(var_f(1,1,0,0),1,rectangle_type_f,former_rank_z,2,newz_comm_world,req2s,ierr)
        call MPI_Irecv(var_f(1,x_procs+1,0,0),1,rectangle_type_f,next_rank_z,2,newz_comm_world,req2r,ierr)

        call MPI_Wait(req1s,sta1s,ierr)
        call MPI_Wait(req1r,sta1r,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
        call MPI_Wait(req2r,sta2r,ierr)
    end subroutine glue_f

    subroutine MPI_boundary(phi_procs,u1_procs,u2_procs,u3_procs,f_procs,g_procs,taug_procs)
        real(8),intent(inout) :: phi_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout) :: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout) :: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout) :: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout) :: f_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout) :: g_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout) :: taug_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        call glue(phi_procs)
        call glue(u1_procs)
        call glue(u2_procs)
        call glue(u3_procs)
        call glue(taug_procs)
        call glue_f(f_procs)
        call glue_f(g_procs)
    endsubroutine MPI_boundary

    subroutine MPI_boundary_fg(feq_procs,geq_procs)
        real(8),intent(inout) :: feq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout) :: geq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        call glue_f(feq_procs)
        call glue_f(geq_procs)
    endsubroutine MPI_boundary_fg

    subroutine feq_cal(gphi_procs,phi_procs,p0_procs,lap_phi_procs,grad_phi_procs,u1_procs,u2_procs,u3_procs,feq_procs)
        real(8),intent(in):: gphi_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: phi_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: p0_procs(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: lap_phi_procs(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: grad_phi_procs(1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(out):: feq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
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
        real(8),intent(in):: gphi_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: p_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(out):: geq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: cr(1:3,1:15)
        do zi=1,z_procs
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
                        + E(i)*kappag*gtemp
                    enddo
                enddo
            enddo
        enddo
    end subroutine geq_cal

    subroutine f_cal(fnext_procs,f_procs,feq_procs)
        real(8),intent(inout):: fnext_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: f_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: feq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
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
        real(8),intent(inout):: gnext_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: g_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: geq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: taug_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
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
        real(8),intent(inout):: fun(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: U(1:3)
        !安定したら壁を動かす
        if((n >= start)) then
            U(1) = uw
            U(2) = 0.0d0
            U(3) = 0.0d0
        endif
        if(group_z == 0) then
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
        else if(group_z == Nz-1) then
            zi = z_procs
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
        real(8),intent(inout):: phi_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: p_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: f_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: g_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: nu_procs(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(inout):: taug_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
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
        do zi=1,z_procs
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
        call MPI_Allreduce(max,phi_max,1,MPI_REAL8,MPI_MAX,new_comm_world,ierr)
        call MPI_Allreduce(min,phi_min,1,MPI_REAL8,MPI_MIN,new_comm_world,ierr)
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    nu_procs(xi,yi,zi) = (phi_procs(xi,yi,zi)-phi_min)/(phi_max-phi_min)*(nu2-nu1) + nu1
                    taug_procs(xi,yi,zi) = 0.5d0 + 3.0d0*nu_procs(xi,yi,zi)/ds + 2.0d0/3.0d0*Anu
                enddo
            enddo
        enddo
    end subroutine physics

    subroutine renew(fun_procs,funnext_procs)
        real(8),intent(in) :: funnext_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(out) :: fun_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
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
        real(8),intent(out):: lap_f(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: fun(0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
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
        real(8),intent(out):: grad_f(1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: fun(0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
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
        real(8),intent(out):: p0_procs(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: phi_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    p0_procs(xi,yi,zi) = phi_procs(xi,yi,zi)*T/(1.0d0-b*phi_procs(xi,yi,zi)) - a*phi_procs(xi,yi,zi)**2
                enddo
            enddo
        enddo
    end subroutine p0_cal

    subroutine gphi_cal(gphi_procs,grad_phi_procs)
        real(8),intent(out):: gphi_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: grad_phi_procs(1:3,1:x_procs,1:ymax+1,1:z_procs)
        do zi=1,z_procs
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
        real(8),intent(out):: grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
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

    subroutine output(phi)
        real(8),intent(inout):: phi(-1:xmax+1,-1:ymax+1,-1:zmax+1)

        write(filename,*) group_new !i->filename 変換
        write(filename2,*) n
        filename=datadir//trim(adjustl(filename))//'_'//trim(adjustl(filename2))//'.d' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
        print *, filename !表示してみる
        open(20,file=filename, status='replace') 

        yi = ymax / 2 
        do zi=0,zmax
            ! do yi=0,ymax
                do xi=0,xmax
                    write(20,"(3es16.8)") dble(xi),dble(zi),phi(xi,yi,zi)
                enddo
                write(20,*)
            ! enddo
        enddo
        close(20)
    end subroutine output

    subroutine laplace(phi,grad_phi,sigma)
        real(8),intent(inout):: phi(-1:xmax+1,-1:ymax+1,-1:zmax+1), grad_phi(1:3,0:xmax,0:ymax,0:zmax), sigma

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
        yi = (ymax) / 2
        zi = (zmax) / 2
        sigma = 0.0d0
        do xi=0,(xmax)/2
            sigma = sigma + grad_phi(1,xi,yi,zi)**2
        enddo
        sigma = sigma * kappag
    end subroutine laplace
    
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
    real(8) p(-1:xmax+1,-1:ymax+1,-1:zmax+1) !圧力
    real(8) U(1:3) !壁の移動速度
    real(8) grad_phi(1:3,0:xmax,0:ymax,0:zmax)
    real(8) lap_phi(0:xmax,0:ymax,0:zmax)
    real(8) gphi(1:3,1:3,0:xmax,0:ymax,0:zmax)
    real(8) p0(0:xmax,0:ymax,0:zmax)
    real(8) grad_u(1:3,1:3,0:xmax,0:ymax,0:zmax)

    real(8) sigma !界面張力

    real(8),allocatable :: f_procs(:,:,:,:), fnext_procs(:,:,:,:), feq_procs(:,:,:,:)
    real(8),allocatable :: g_procs(:,:,:,:), gnext_procs(:,:,:,:), geq_procs(:,:,:,:)
    real(8),allocatable :: phi_procs(:,:,:), p_procs(:,:,:)
    real(8),allocatable :: u1_procs(:,:,:), u2_procs(:,:,:), u3_procs(:,:,:)

    real(8),allocatable :: taug_procs(:,:,:), nu_procs(:,:,:)
    real(8),allocatable :: p0_procs(:,:,:), lap_phi_procs(:,:,:), grad_phi_procs(:,:,:,:), gphi_procs(:,:,:,:,:), grad_u_procs(:,:,:,:,:)
    
    real(8),allocatable :: phiout(:,:,:,:)

    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, comm_procs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, comm_rank, ierr)

    call mk_dirs(datadir)
    call par(cx,cy,cz,cr,krone,U)
    call ini(phi,u1,u2,u3,p,geq,g_procs,gnext_procs,feq,f_procs,fnext_procs,feq_procs,geq_procs,phi_procs,p_procs,u1_procs,u2_procs,u3_procs,grad_phi,gphi,lap_phi,p0,grad_u)
    call ini_op(taug_procs,nu_procs,p0_procs,lap_phi_procs,grad_phi_procs,gphi_procs,grad_u_procs,phi_procs)
    
    allocate(phiout(1:x_procs,1:ymax+1,1:z_procs,0:new_procs-1))
    if(new_rank == 0) then
        write(chmyrank, "(i2.2)") group_new
        open(10,file="./laplace_"//chmyrank//".d")
    endif
DO n=1,step
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

    if( mod(n,100)==0 ) then
        !phiを出力する
        tmp(:,:,:) = 0.0d0
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp(xi,yi,zi) = phi_procs(xi,yi,zi)
                enddo
            enddo
        enddo
        do i=1,new_procs-1
            if(new_rank == i) then
                tags = i
                call MPI_Isend(tmp(1,1,1),x_procs*z_procs*(ymax+1),MPI_REAL8,0,tags,new_comm_world,req1s,ierr)
                call MPI_Wait(req1s,sta1s,ierr)
            endif
            if(new_rank == 0) then
                tagr = i
                recv_rank = i
                call MPI_Irecv(phiout(1,1,1,i),x_procs*z_procs*(ymax+1),MPI_REAL8,recv_rank,tagr,new_comm_world,req1r,ierr)
                call MPI_Wait(req1r,sta1r,ierr)
            endif
        enddo
        if(new_rank == 0) then
            do zi=1,z_procs
                do yi=1,ymax+1
                    do xi=1,x_procs
                        phiout(xi,yi,zi,0) = phi_procs(xi,yi,zi)
                    enddo
                enddo
            enddo
        endif
        !物理量をまとめる
        if(new_rank == 0) then
            k = 0
            do Nzz=0,Nz-1
                do Nxx=0,Nx-1
                    do zi=1,z_procs
                        do yi=1,ymax+1
                            do xi=1,x_procs
                                phi((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*z_procs) = phiout(xi,yi,zi,k)
                            enddo
                        enddo
                    enddo
                    k = k + 1
                enddo
            enddo
            call laplace(phi,grad_phi,sigma)
            call output(phi)
            write(*,*) "step=", n, sigma
            write(10,*) "Myrank :", group_new, n, sigma
        endif
    endif
ENDDO
    close(10)
    call MPI_Finalize(ierr)
end program main