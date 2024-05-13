!【LES】
! Smagorinskyモデルを使用し、渦粘性の効果を導入している。
! 【Multi-times simulations】
! このプログラムは、多数回の液滴分裂の数値シミュレーションを実行する。
! このプログラムは、界面を「Cahn-Hilliard 方程式」で表現している。
! 液滴を長時間保持できないが、系全体の保存性には優れている。
! 【高速計算】
! このプログラムは、MPI + OpenMP のハイブリット並列により高速計算を可能としている。
! 【系について】
! 単相のシミュレーションを基本的には行う。
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module globals
    !$ use omp_lib
    include "mpif.h"
    !計算領域
    real(8),parameter:: ds = 1.0d0 !格子間隔（lattice unit）
    integer,parameter:: xmax = 255 !ｘ方向格子数（０から数える）
    integer,parameter:: ymax = 255 !ｙ方向格子数（０から数える）
    integer,parameter:: zmax = 255 !ｚ方向格子数（０から数える）
    !MPI並列数
    integer,parameter:: MPI_X = 8 !x方向並列数
    !パラメータ並列数
    integer,parameter:: Nxall = 1 !x方向の分割数（Nxall*Nzall=全体の計算数）
    integer,parameter:: Nzall = 1 !z方向の分割数（Nxall*Nzall=全体の計算数）
    integer,parameter:: xall = (xmax + 1) * Nxall !全体のx方向格子数
    integer,parameter:: zall = (zmax + 1) * Nzall !全体のz方向格子数
    !時間に関するパラメータ
    integer,parameter:: step_start = 0 !計算開始ステップ数
    integer,parameter:: step_end = 10000000 !計算終了ステップ数
    !乱流場の生成(Frag = 0)，液滴分裂の計算(Frag = 1)，ファイルを入力して途中から計算(Frag = 2)
    integer,parameter:: Frag = 0
    !出力ディレクトリ
    character(*),parameter :: datadir = "/data/sht/nakanog/taylor_re12000_LES_2/"
    character(*),parameter :: datadir2 = "/data/sht/nakanog/taylor_re12000_LES_2/fg/"

    !LESパラメータ
    real(8),parameter:: delta = 8.0d0*ds
    real(8),parameter:: Cs = 0.173d0

    !無次元数
    real(8),parameter:: We = 3.0d0 !粒子ウエーバー数
    real(8),parameter:: Re = 12000.0d0 !粒子レイノルズ数
    real(8),parameter:: eta = 1.0d0 !粘度比（nu2/nu1）

    !撹乱（乱数）のオーダー
    real(8), parameter :: ran = 1.0d-3

    !支配パラメータ
    real(8),parameter:: pi = acos(-1.0d0) !円周率
    real(8),parameter:: D = 0.0d0 !設置する液滴直径
    real(8),parameter:: D_vortex = 127.5d0 !渦の大きさ
    real(8),parameter:: kw = pi/D_vortex !波数
    real(8),parameter:: umax = 0.05d0 !最大流速
    real(8),parameter:: nu1 = umax*D_vortex/Re !連続相の粘性係数
    real(8),parameter:: nu2 = eta*nu1 !分散相の粘性係数
    real(8),parameter:: sigma = D*umax*umax/We !界面張力
    ! real(8),parameter:: kappag = (sigma/1.7039d0)**(1.0d0/0.9991d0)  !界面張力を決めるパラメータ
    real(8),parameter:: kappag = 1.0d-5  !界面張力を決めるパラメータ
    real(8),parameter:: kappaf = 0.01d0*ds**2 !界面厚さを決めるパラメータ
    real(8),parameter:: phi1 = 2.211d0 !連続相のオーダーパラメータ
    real(8),parameter:: phi2 = 4.895d0 !分散相のオーダーパラメータ
    real(8),parameter:: a = 9.0d0/49.0d0
    real(8),parameter:: b = 2.0d0/21.0d0
    real(8),parameter:: T = 0.55d0
    real(8),parameter:: tauf = 0.7d0 !緩和時間
    real(8),parameter:: Anu = 0.0d0 !粘性係数が小さいときに値を入れると良い（O）

    !初期液滴中心
    real(8), parameter:: xc = 0.5d0*dble(xmax)	
    real(8), parameter:: yc = 0.5d0*dble(ymax)
    real(8), parameter:: zc = 0.5d0*dble(zmax)

    !粒子速度（整数）
    integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
    integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
    integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)
    real(8) cr(1:3, 1:15)  !粒子速度（実数）
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

    !その他変数
    real(8) phi_min, phi_max, min, max
    real(8) gtemp, ftemp
    real(8) dif
    integer grobalx, grobaly, grobalz
    real(8) dummy
    real(8) time1, time2 !計算時間測定用
    
    character :: filename*200
    character :: filename2*200
    character :: filename3*200

    !MPI用変数
    integer ierr, comm_procs, comm_rank
    integer Nx, Nz 
    integer x_procs, z_procs
    integer key_new, group_new, key_x, group_x, key_z, group_z
    integer new_comm_world, new_procs, new_rank
    integer newx_comm_world, newx_procs, newx_rank
    integer newz_comm_world, newz_procs, newz_rank
    integer rectangle_type, rectangle_type_f, xtype, xtype_f
    integer next_rank_x, former_rank_x, next_rank_z, former_rank_z
    integer req1s, req1r, req2s, req2r
    integer, dimension(MPI_STATUS_SIZE) :: sta1s, sta1r, sta2s, sta2r
    integer Nxx, Nzz
    integer tags, tagr, recv_rank
    character(2) chmyrank

contains
    subroutine par(cx,cy,cz,cr,krone)
        integer, intent(in):: cx(15),cy(15),cz(15)
        real(8), intent(out):: cr(1:3,1:15), krone(1:3,1:3)
        integer i
        do i=1,15
            cr(1,i) = dble(cx(i))
            cr(2,i) = dble(cy(i))
            cr(3,i) = dble(cz(i))
        enddo
        krone(1,1) = 1.d0; krone(1,2) = 0.d0; krone(1,3) = 0.d0
        krone(2,1) = 0.d0; krone(2,2) = 1.d0; krone(2,3) = 0.d0
        krone(3,1) = 0.d0; krone(3,2) = 0.d0; krone(3,3) = 1.d0
    end subroutine par

    subroutine mk_random(x_m,y_m,eps,eps2)
        real(8),allocatable :: x_m(:,:,:), y_m(:,:,:), eps(:,:,:), eps2(:,:,:)
        integer seedsize
        integer, allocatable :: seeds(:), seeds2(:)
        integer seed, seed2

        integer i, xi, yi, zi

        allocate(x_m(1:x_procs,1:ymax+1,1:z_procs))
        allocate(y_m(1:x_procs,1:ymax+1,1:z_procs))
        allocate(eps(1:x_procs,1:ymax+1,1:z_procs))
        allocate(eps2(1:x_procs,1:ymax+1,1:z_procs))

        !================================乱数の発生（標準正規分布）=====================================
        call random_seed(size=seedsize)
        allocate(seeds(1:seedsize))
        allocate(seeds2(1:seedsize))

        call system_clock(seed)
        seed = seed + comm_rank
        do i=1,seedsize
            seeds(i) = seed + i
        enddo
        call random_seed(put=seeds)
        call random_number(eps)

        call system_clock(seed2)
        seed2 = seed + seed2 + (2 * comm_rank)
        do i=1,seedsize
            seeds2(i) = seed2 + i
        enddo
        call random_seed(put=seeds2)
        call random_number(eps2)

        !$omp parallel
        !$omp do
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    x_m(xi,yi,zi) = sqrt(-2.0d0*log(eps(xi,yi,zi)))*cos(2.0d0*pi*eps2(xi,yi,zi))
                    y_m(xi,yi,zi) = sqrt(-2.0d0*log(eps(xi,yi,zi)))*sin(2.0d0*pi*eps2(xi,yi,zi))
                    if(x_m(xi,yi,zi) > 3.0d0) then
                        x_m(xi,yi,zi) = 3.0d0
                    elseif(x_m(xi,yi,zi) < -3.0d0) then
                        x_m(xi,yi,zi) = -3.0d0
                    endif

                    if(y_m(xi,yi,zi) > 3.0d0) then
                        y_m(xi,yi,zi) = 3.0d0
                    elseif(y_m(xi,yi,zi) < -3.0d0) then
                        y_m(xi,yi,zi) = -3.0d0
                    endif
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine mk_random

    subroutine ini(g_procs,gnext_procs,f_procs,fnext_procs,feq_procs,geq_procs,phi_procs,p_procs,u1_procs,u2_procs,u3_procs,grad_phi_procs,gphi_procs,lap_phi_procs,p0_procs,grad_u_procs,x_m,y_m,taug_procs,nu_procs,strain_procs,nue)
        real(8),allocatable :: f_procs(:,:,:,:), fnext_procs(:,:,:,:), feq_procs(:,:,:,:)
        real(8),allocatable :: g_procs(:,:,:,:), gnext_procs(:,:,:,:), geq_procs(:,:,:,:)
        real(8),allocatable :: phi_procs(:,:,:), p_procs(:,:,:)
        real(8),allocatable :: u1_procs(:,:,:), u2_procs(:,:,:), u3_procs(:,:,:)
        real(8),allocatable :: p0_procs(:,:,:), lap_phi_procs(:,:,:), grad_phi_procs(:,:,:,:)
        real(8),allocatable :: gphi_procs(:,:,:,:,:), grad_u_procs(:,:,:,:,:)
        real(8),allocatable :: taug_procs(:,:,:), nu_procs(:,:,:)
        real(8),allocatable :: strain_procs(:,:,:,:,:)
        real(8),allocatable :: nue(:,:,:)
        real(8),intent(in) :: x_m(1:x_procs,1:ymax+1,1:z_procs), y_m(1:x_procs,1:ymax+1,1:z_procs)

        integer i, xi, yi, zi, alpha, beta

    !========================並列数・コミュニケータ分割・通信先設定========================
        !配列のallocate(xi:0とx_procs+1がのりしろ)(yi:0とymax+2がのりしろ)(zi:0とz_procs+1がのりしろ)
        Nx = MPI_X !x方向の並列数（ただし，Nx/=comm_procs）
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
        allocate(p0_procs(1:x_procs,1:ymax+1,1:z_procs))
        allocate(lap_phi_procs(1:x_procs,1:ymax+1,1:z_procs))
        allocate(grad_phi_procs(1:3,1:x_procs,1:ymax+1,1:z_procs))
        allocate(gphi_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs))
        allocate(grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs))
        allocate(taug_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(nu_procs(1:x_procs,1:ymax+1,1:z_procs))
        allocate(strain_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs))
        allocate(nue(1:x_procs,1:ymax+1,1:z_procs))

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
        p0_procs(:,:,:) = 0.0d0
        lap_phi_procs(:,:,:) = 0.0d0
        grad_phi_procs(:,:,:,:) = 0.0d0
        gphi_procs(:,:,:,:,:) = 0.0d0
        grad_u_procs(:,:,:,:,:) = 0.0d0
        taug_procs(:,:,:) = 0.0d0
        nu_procs(:,:,:) = 0.0d0
        strain_procs(:,:,:,:,:) = 0.0d0
        nue(:,:,:) = 0.0d0
        !===================================初期条件の設定================================================
        if(Frag == 0) then
            do zi=1,z_procs
                do yi=1,ymax+1
                    do xi=1,x_procs
                        phi_procs(xi,yi,zi) = phi1
                        p_procs(xi,yi,zi) = 0.0d0
                        u1_procs(xi,yi,zi) = x_m(xi,yi,zi) * ran
                        u2_procs(xi,yi,zi) = y_m(xi,yi,zi) * ran
                        u3_procs(xi,yi,zi) = 0.5d0 * ( x_m(xi,yi,zi) + y_m(xi,yi,zi) ) * ran
                    enddo
                enddo
            enddo
        else if(Frag == 1) then
            do zi=1,z_procs
                do yi=1,ymax+1
                    do xi=1,x_procs
                        grobalx = (xi-1) + group_x * x_procs
                        grobaly = yi-1
                        grobalz = (zi-1) + group_z * z_procs
                        dif = (dble(grobalx)*ds-xc)**2 + (dble(grobaly)*ds-yc)**2 + (dble(grobalz)*ds-zc)**2
                        if(dif <= (0.5d0*D)**2) then
                            phi_procs(xi,yi,zi) = phi2
                        else
                            phi_procs(xi,yi,zi) = phi1
                        endif
                        p_procs(xi,yi,zi) = 0.0d0
                        u1_procs(xi,yi,zi) = 0.0d0
                        u2_procs(xi,yi,zi) = 0.0d0
                        u3_procs(xi,yi,zi) = 0.0d0
                    enddo
                enddo
            enddo
        endif

        call glue(phi_procs)
        call glue(u1_procs)
        call glue(u2_procs)
        call glue(u3_procs)
        call various_cal(lap_phi_procs,grad_phi_procs,grad_u_procs,p0_procs,phi_procs,u1_procs,u2_procs,u3_procs)
        call gphi_cal(gphi_procs,grad_phi_procs)
        call equilibrium_cal(gphi_procs,phi_procs,p0_procs,lap_phi_procs,grad_phi_procs,grad_u_procs,p_procs,u1_procs,u2_procs,u3_procs,f_procs,g_procs)
        call strain_cal(grad_u_procs,strain_procs)
    
        !$omp parallel
        !$omp do
        do zi = 1, z_procs
            do yi = 1, ymax+1
                do xi = 1, x_procs
                    do beta = 1, 3
                        do alpha = 1, 3
                            nue(xi,yi,zi) = nue(xi,yi,zi) + strain_procs(alpha,beta,xi,yi,zi)**2.0d0
                        enddo
                    enddo
                    nue(xi,yi,zi) = ( 2.0d0 * nue(xi,yi,zi) )**0.5d0
                    nue(xi,yi,zi) = Cs * delta * delta * nue(xi,yi,zi)

                    nu_procs(xi,yi,zi) = (phi_procs(xi,yi,zi)-phi1)/(phi2-phi1)*(nu2-nu1) + nu1
                    taug_procs(xi,yi,zi) = 0.5d0 + 3.0d0 * (nu_procs(xi,yi,zi) + nue(xi,yi,zi)) / ds + 2.0d0/3.0d0 * Anu
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine ini

    subroutine ini_restart(g_procs,gnext_procs,f_procs,fnext_procs,feq_procs,geq_procs,phi_procs,p_procs,u1_procs,u2_procs,u3_procs,grad_phi_procs,gphi_procs,lap_phi_procs,p0_procs,grad_u_procs,taug_procs,nu_procs,strain_procs,nue)
        real(8),allocatable :: f_procs(:,:,:,:), fnext_procs(:,:,:,:), feq_procs(:,:,:,:)
        real(8),allocatable :: g_procs(:,:,:,:), gnext_procs(:,:,:,:), geq_procs(:,:,:,:)
        real(8),allocatable :: phi_procs(:,:,:), p_procs(:,:,:)
        real(8),allocatable :: u1_procs(:,:,:), u2_procs(:,:,:), u3_procs(:,:,:)
        real(8),allocatable :: p0_procs(:,:,:), lap_phi_procs(:,:,:), grad_phi_procs(:,:,:,:)
        real(8),allocatable :: gphi_procs(:,:,:,:,:), grad_u_procs(:,:,:,:,:)
        real(8),allocatable :: taug_procs(:,:,:), nu_procs(:,:,:)
        real(8),allocatable :: strain_procs(:,:,:,:,:)
        real(8),allocatable :: nue(:,:,:)

        integer i, xi, yi, zi

    !========================並列数・コミュニケータ分割・通信先設定========================
        !配列のallocate(xi:0とx_procs+1がのりしろ)(yi:0とymax+2がのりしろ)(zi:0とz_procs+1がのりしろ)
        Nx = MPI_X !x方向の並列数（ただし，Nx/=comm_procs）
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
        allocate(p0_procs(1:x_procs,1:ymax+1,1:z_procs))
        allocate(lap_phi_procs(1:x_procs,1:ymax+1,1:z_procs))
        allocate(grad_phi_procs(1:3,1:x_procs,1:ymax+1,1:z_procs))
        allocate(gphi_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs))
        allocate(grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs))
        allocate(taug_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(nu_procs(1:x_procs,1:ymax+1,1:z_procs))
        allocate(strain_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs))
        allocate(nue(1:x_procs,1:ymax+1,1:z_procs))

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
        p0_procs(:,:,:) = 0.0d0
        lap_phi_procs(:,:,:) = 0.0d0
        grad_phi_procs(:,:,:,:) = 0.0d0
        gphi_procs(:,:,:,:,:) = 0.0d0
        grad_u_procs(:,:,:,:,:) = 0.0d0
        taug_procs(:,:,:) = 0.0d0
        nu_procs(:,:,:) = 0.0d0

        !===============================途中から計算再開=================================================
        write(filename,*) group_new !i->filename 変換
        write(filename2,*) step_start
        write(filename3,*) new_rank
        filename=datadir2//trim(adjustl(filename))//'_'//trim(adjustl(filename3))//'_'//trim(adjustl(filename2))//'_fg.bin' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
        print *, filename !表示してみる
        open(99, file=filename, form="unformatted")
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        read(99) f_procs(i,xi,yi,zi),g_procs(i,xi,yi,zi)
                    enddo
                enddo
            enddo
        enddo
        close(99)
        
        !$omp parallel
        !$omp do
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
        !$omp end do
        !$omp end parallel

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

        call grad_u_cal(grad_u_procs,u1_procs,u2_procs,u3_procs)
        call strain_cal(grad_u_procs,strain_procs)
        nue(:,:,:) = 0.0d0
        !$omp parallel
        !$omp do
        do zi = 1, z_procs
            do yi = 1, ymax+1
                do xi = 1, x_procs
                    do beta = 1, 3
                        do alpha = 1, 3
                            nue(xi,yi,zi) = nue(xi,yi,zi) + strain_procs(alpha,beta,xi,yi,zi)**2.0d0
                        enddo
                    enddo
                    nue(xi,yi,zi) = ( 2.0d0 * nue(xi,yi,zi) )**0.5d0
                    nue(xi,yi,zi) = Cs * delta * delta * nue(xi,yi,zi)

                    nu_procs(xi,yi,zi) = (phi_procs(xi,yi,zi)-phi_min)/(phi_max-phi_min)*(nu2-nu1) + nu1
                    taug_procs(xi,yi,zi) = 0.5d0 + 3.0d0 * (nu_procs(xi,yi,zi) + nue(xi,yi,zi)) / ds + 2.0d0/3.0d0 * Anu
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine ini_restart

    subroutine glue(var)
        real(8),intent(inout) :: var(0:x_procs+1,0:ymax+2,0:z_procs+1)
        integer xi, zi

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
        !Z世界でののりしろ通信
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
        integer i, xi, zi

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

        !Z世界でののりしろ通信
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

    subroutine grad_u_cal(grad_u_procs,u1_procs,u2_procs,u3_procs)
        real(8),intent(out):: grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        integer xi, yi, zi, i, beta

        !$omp parallel
        !$omp do
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
        !$omp end do
        !$omp end parallel
    end subroutine grad_u_cal

    subroutine strain_cal(grad_u_procs,strain_procs)
        real(8),intent(in):: grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(out):: strain_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        integer xi, yi, zi, alpha, beta 

        !$omp parallel
        !$omp do
        do zi = 1, z_procs
            do yi = 1, ymax+1
                do xi = 1, x_procs
                    do beta = 1, 3
                        do alpha = 1, 3
                            strain_procs(alpha,beta,xi,yi,zi) = 0.5d0*( grad_u_procs(alpha,beta,xi,yi,zi) + grad_u_procs(beta,alpha,xi,yi,zi) )
                        enddo
                    enddo
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine strain_cal

    subroutine equilibrium_cal(gphi_procs,phi_procs,p0_procs,lap_phi_procs,grad_phi_procs,grad_u_procs,p_procs,u1_procs,u2_procs,u3_procs,feq_procs,geq_procs)
        real(8),intent(in):: gphi_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: phi_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: p0_procs(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: lap_phi_procs(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: grad_phi_procs(1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: p_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(out):: feq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(out):: geq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        integer i, alpha, beta, xi, yi, zi

        !$omp parallel
        !$omp do
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        ftemp = 0.0d0
                        gtemp = 0.0d0
                        do beta=1,3
                            do alpha=1,3
                                ftemp = ftemp + (grad_u_procs(beta,alpha,xi,yi,zi) &
                                                    + grad_u_procs(alpha,beta,xi,yi,zi))*cr(alpha,i)*cr(beta,i)

                                gtemp = gtemp + gphi_procs(alpha,beta,xi,yi,zi)*cr(alpha,i)*cr(beta,i)
                            enddo
                        enddo
                        feq_procs(i,xi,yi,zi) = H(i)*phi_procs(xi,yi,zi) &
                                        + F(i) * (p0_procs(xi,yi,zi) &
                                        - kappaf*phi_procs(xi,yi,zi)*lap_phi_procs(xi,yi,zi) &
                                        - kappaf/6.0d0*(grad_phi_procs(1,xi,yi,zi)**2 &
                                        + grad_phi_procs(2,xi,yi,zi)**2+grad_phi_procs(3,xi,yi,zi)**2)) &
                                        + E(i)*phi_procs(xi,yi,zi)*3.0d0*(cr(1,i)*u1_procs(xi,yi,zi) &
                                        + cr(2,i)*u2_procs(xi,yi,zi)+cr(3,i)*u3_procs(xi,yi,zi)) &
                                        + E(i)*kappaf*gtemp

                        geq_procs(i,xi,yi,zi) = E(i)*(3.0d0*p_procs(xi,yi,zi) &
                                        + 3.0d0*(u1_procs(xi,yi,zi)*dble(cx(i)) &
                                        + u2_procs(xi,yi,zi)*dble(cy(i)) + u3_procs(xi,yi,zi)*dble(cz(i))) &
                                        + 4.5d0*(u1_procs(xi,yi,zi)*dble(cx(i)) &
                                        + u2_procs(xi,yi,zi)*dble(cy(i)) + u3_procs(xi,yi,zi)*dble(cz(i)))**2 &
                                        - 1.5d0*(u1_procs(xi,yi,zi)**2 &
                                        + u2_procs(xi,yi,zi)**2 + u3_procs(xi,yi,zi)**2) &
                                        + Anu*ds*ftemp) &
                                        + E(i)*kappag*gtemp
                    enddo
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine equilibrium_cal

    subroutine fg_cal(fnext_procs,f_procs,feq_procs,gnext_procs,g_procs,geq_procs,taug_procs,forcex,forcey,forcez)
        real(8),intent(inout):: fnext_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: f_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: feq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)

        real(8),intent(inout):: gnext_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: g_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: geq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: taug_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: forcex(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: forcey(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: forcez(0:x_procs+1,0:ymax+2,0:z_procs+1)
        integer i, xi, yi, zi

        !$omp parallel
        !$omp do
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        fnext_procs(i,xi,yi,zi) = f_procs(i,xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                - (f_procs(i,xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                - feq_procs(i,xi-cx(i),yi-cy(i),zi-cz(i))) / tauf

                        gnext_procs(i,xi,yi,zi) = g_procs(i,xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                - (g_procs(i,xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                - geq_procs(i,xi-cx(i),yi-cy(i),zi-cz(i))) &
                                                / taug_procs(xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                + 3.0d0 * ds * E(i) &
                                                * (dble(cx(i))*forcex(xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                + dble(cy(i))*forcey(xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                + dble(cz(i))*forcez(xi-cx(i),yi-cy(i),zi-cz(i)))
                    enddo
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine fg_cal

    subroutine physics(phi_procs,p_procs,u1_procs,u2_procs,u3_procs,f_procs,g_procs,nu_procs,taug_procs,fnext_procs,gnext_procs,strain_procs,grad_u_procs,nue)
        real(8),intent(inout):: phi_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: p_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: f_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: g_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: nu_procs(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(inout):: taug_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in) :: fnext_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in) :: gnext_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout) :: strain_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(inout) :: grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(inout) :: nue(1:x_procs,1:ymax+1,1:z_procs)
        integer i, xi, yi, zi, beta, alpha

        !$omp parallel
        !$omp do
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        f_procs(i,xi,yi,zi) = fnext_procs(i,xi,yi,zi)
                        g_procs(i,xi,yi,zi) = gnext_procs(i,xi,yi,zi)
                    enddo

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
        !$omp end do
        !$omp end parallel

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

        call grad_u_cal(grad_u_procs,u1_procs,u2_procs,u3_procs)
        call strain_cal(grad_u_procs,strain_procs)
        nue(:,:,:) = 0.0d0
        !$omp parallel
        !$omp do
        do zi = 1, z_procs
            do yi = 1, ymax+1
                do xi = 1, x_procs
                    do beta = 1, 3
                        do alpha = 1, 3
                            nue(xi,yi,zi) = nue(xi,yi,zi) + strain_procs(alpha,beta,xi,yi,zi)**2.0d0
                        enddo
                    enddo
                    nue(xi,yi,zi) = ( 2.0d0 * nue(xi,yi,zi) )**0.5d0
                    nue(xi,yi,zi) = Cs * delta * delta * nue(xi,yi,zi)
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel

        !$omp parallel
        !$omp do
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    nu_procs(xi,yi,zi) = (phi_procs(xi,yi,zi)-phi_min)/(phi_max-phi_min)*(nu2-nu1) + nu1
                    taug_procs(xi,yi,zi) = 0.5d0 + 3.0d0 * (nu_procs(xi,yi,zi) + nue(xi,yi,zi)) / ds + 2.0d0/3.0d0 * Anu
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine physics

    subroutine various_cal(lap_f,grad_f,grad_u_procs,p0_procs,fun,u1_procs,u2_procs,u3_procs)
        real(8),intent(out):: lap_f(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(out):: grad_f(1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(out):: grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(out):: p0_procs(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: fun(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        integer alpha, i, xi, yi, zi

        !$omp parallel
        !$omp do
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    p0_procs(xi,yi,zi) = fun(xi,yi,zi)*T/(1.0d0-b*fun(xi,yi,zi)) - a*fun(xi,yi,zi)**2

                    lap_f(xi,yi,zi) = -14.0d0*fun(xi,yi,zi)

                    do i=2,15
                        lap_f(xi,yi,zi) = lap_f(xi,yi,zi) + fun(xi+cx(i),yi+cy(i),zi+cz(i))
                    enddo
                    lap_f(xi,yi,zi) = lap_f(xi,yi,zi) / (5.0d0*ds**2)

                    do alpha=1,3
                        grad_f(alpha,xi,yi,zi) = 0.0d0
                        grad_u_procs(1,alpha,xi,yi,zi) = 0.0d0
                        grad_u_procs(2,alpha,xi,yi,zi) = 0.0d0
                        grad_u_procs(3,alpha,xi,yi,zi) = 0.0d0
                        do i = 2,15
                            grad_f(alpha,xi,yi,zi) = grad_f(alpha,xi,yi,zi) &
                                                    + cr(alpha,i)*fun(xi+cx(i),yi+cy(i),zi+cz(i))
                            grad_u_procs(1,alpha,xi,yi,zi) = grad_u_procs(1,alpha,xi,yi,zi) &
                                                            + cr(alpha,i)*u1_procs(xi+cx(i),yi+cy(i),zi+cz(i))
                            grad_u_procs(2,alpha,xi,yi,zi) = grad_u_procs(2,alpha,xi,yi,zi) &
                                                            + cr(alpha,i)*u2_procs(xi+cx(i),yi+cy(i),zi+cz(i))
                            grad_u_procs(3,alpha,xi,yi,zi) = grad_u_procs(3,alpha,xi,yi,zi) &
                                                            + cr(alpha,i)*u3_procs(xi+cx(i),yi+cy(i),zi+cz(i))
                        enddo
                        grad_f(alpha,xi,yi,zi) = grad_f(alpha,xi,yi,zi) / (10.0d0*ds)
                        grad_u_procs(1,alpha,xi,yi,zi) = grad_u_procs(1,alpha,xi,yi,zi)/(10.0d0*ds)
                        grad_u_procs(2,alpha,xi,yi,zi) = grad_u_procs(2,alpha,xi,yi,zi)/(10.0d0*ds)
                        grad_u_procs(3,alpha,xi,yi,zi) = grad_u_procs(3,alpha,xi,yi,zi)/(10.0d0*ds)
                    enddo
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine various_cal

    subroutine gphi_cal(gphi_procs,grad_phi_procs)
        real(8),intent(out):: gphi_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: grad_phi_procs(1:3,1:x_procs,1:ymax+1,1:z_procs)
        integer alpha, beta, xi, yi, zi

        !$omp parallel
        !$omp do
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do beta=1,3
                        do alpha=1,3
                            gphi_procs(alpha,beta,xi,yi,zi) = 4.5d0*grad_phi_procs(alpha,xi,yi,zi)*grad_phi_procs(beta,xi,yi,zi) &
                                                            -1.5d0*(grad_phi_procs(1,xi,yi,zi)**2 &
                                                            + grad_phi_procs(2,xi,yi,zi)**2 &
                                                            + grad_phi_procs(3,xi,yi,zi)**2) &
                                                            *krone(alpha,beta)
                        enddo
                    enddo
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine gphi_cal

    subroutine externalforce(forcex,forcey,forcez)
        real(8),allocatable :: forcex(:,:,:), forcey(:,:,:), forcez(:,:,:)
        integer xi, yi, zi

        allocate(forcex(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(forcey(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(forcez(0:x_procs+1,0:ymax+2,0:z_procs+1))
        forcex(:,:,:) = 0.0d0
        forcey(:,:,:) = 0.0d0
        forcez(:,:,:) = 0.0d0

        !$omp parallel
        !$omp do
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    grobalx = (xi-1) + group_x * x_procs
                    grobalz = (zi-1) + group_z * z_procs

                    forcex(xi,yi,zi) = 2.0d0 * pi * pi * umax * umax &
                                        / (D_vortex * Re) * sin(kw*dble(grobalx))  * cos(kw*dble(grobalz))
                    forcey(xi,yi,zi) = 0.0d0
                    forcez(xi,yi,zi) = -2.0d0 * pi * pi * umax * umax &
                                        / (D_vortex * Re) * cos(kw*dble(grobalx)) * sin(kw*dble(grobalz))
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    endsubroutine externalforce

    subroutine output(phi_procs,u1_procs,u2_procs,u3_procs,p_procs,n)
        real(8),intent(inout):: phi_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: p_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        integer,intent(in):: n
        integer xi, yi, zi

        write(filename,*) group_new !i->filename 変換
        write(filename2,*) n
        write(filename3,*) new_rank
        filename=datadir//trim(adjustl(filename))//'_'//trim(adjustl(filename3))//'_'//trim(adjustl(filename2))//'.bin' 
        print *, filename !表示してみる
        open(100,file=filename, form='unformatted',status='replace')
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    write(100) u1_procs(xi,yi,zi),u2_procs(xi,yi,zi),u3_procs(xi,yi,zi), p_procs(xi,yi,zi), phi_procs(xi,yi,zi)
                enddo
            enddo
        enddo
        close(100)
    end subroutine output

    subroutine outputphi(phi_procs,n)
        real(8),intent(inout):: phi_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        integer,intent(in):: n
        integer xi, yi, zi

        write(filename,*) group_new !i->filename 変換
        write(filename2,*) n
        write(filename3,*) new_rank
        filename=datadir//trim(adjustl(filename))//'_'//trim(adjustl(filename3))//'_'//trim(adjustl(filename2))//'.bin' 
        print *, filename !表示してみる
        open(101,file=filename, form='unformatted',status='replace') 
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    write(101) phi_procs(xi,yi,zi)
                enddo
            enddo
        enddo
        close(101)
    end subroutine outputphi

    subroutine outputfg(f_procs, g_procs, n)
        real(8),intent(inout):: f_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: g_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        integer,intent(in):: n
        integer i, xi, yi, zi

        write(filename,*) group_new !i->filename 変換
        write(filename2,*) n
        write(filename3,*) new_rank
        filename=datadir2//trim(adjustl(filename))//'_'//trim(adjustl(filename3))//'_'//trim(adjustl(filename2))//'_fg.bin' 
        print *, filename !表示してみる
        open(102,file=filename, form='unformatted',status='replace') 
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        write(102) f_procs(i,xi,yi,zi), g_procs(i,xi,yi,zi)
                    enddo
                enddo
            enddo
        enddo
        close(102)
    end subroutine outputfg
    
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
    !LBMの計算に使用する変数
    real(8),allocatable :: f_procs(:,:,:,:), fnext_procs(:,:,:,:), feq_procs(:,:,:,:)
    real(8),allocatable :: g_procs(:,:,:,:), gnext_procs(:,:,:,:), geq_procs(:,:,:,:)
    real(8),allocatable :: phi_procs(:,:,:), p_procs(:,:,:)
    real(8),allocatable :: u1_procs(:,:,:), u2_procs(:,:,:), u3_procs(:,:,:)
    real(8),allocatable :: taug_procs(:,:,:), nu_procs(:,:,:)
    real(8),allocatable :: p0_procs(:,:,:), lap_phi_procs(:,:,:), grad_phi_procs(:,:,:,:), gphi_procs(:,:,:,:,:)
    real(8),allocatable :: grad_u_procs(:,:,:,:,:)
    real(8),allocatable :: forcex(:,:,:), forcey(:,:,:), forcez(:,:,:)
    real(8),allocatable :: strain_procs(:,:,:,:,:)
    real(8),allocatable :: nue(:,:,:)
    integer n, i, xi, yi, zi

    !乱数変数
    real(8),allocatable :: eps(:,:,:), eps2(:,:,:), x_m(:,:,:), y_m(:,:,:)

    !解析
    real(8),allocatable ::  kinetic_procs(:,:,:)
    real(8) kinetic_sum, kinetic_ave

!==================================MPI並列開始===============================================
    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, comm_procs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, comm_rank, ierr)
!================================ディレクトリの作成============================================
    call mk_dirs(datadir)
    call mk_dirs(datadir2)
!================================値の代入======================================
    do i = 1, 15
        cr(1,i) = dble(cx(i))
        cr(2,i) = dble(cy(i))
        cr(3,i) = dble(cz(i))
    enddo
    krone(1,1) = 1.d0; krone(1,2) = 0.d0; krone(1,3) = 0.d0
    krone(2,1) = 0.d0; krone(2,2) = 1.d0; krone(2,3) = 0.d0
    krone(3,1) = 0.d0; krone(3,2) = 0.d0; krone(3,3) = 1.d0
!=================================初期条件=======================================
    if(Frag == 0) then
        call ini(g_procs,gnext_procs,f_procs,fnext_procs,feq_procs,geq_procs,phi_procs,p_procs,u1_procs,u2_procs,u3_procs,grad_phi_procs,gphi_procs,lap_phi_procs,p0_procs,grad_u_procs,x_m,y_m,taug_procs,nu_procs,strain_procs,nue)
        call externalforce(forcex,forcey,forcez)
        call glue(forcex)
        call glue(forcey)
        call glue(forcez)
    else if(Frad == 1) then
        call ini(g_procs,gnext_procs,f_procs,fnext_procs,feq_procs,geq_procs,phi_procs,p_procs,u1_procs,u2_procs,u3_procs,grad_phi_procs,gphi_procs,lap_phi_procs,p0_procs,grad_u_procs,x_m,y_m,taug_procs,nu_procs,strain_procs,nue)
    else if(Frag == 2) then
        call ini_restart(g_procs,gnext_procs,f_procs,fnext_procs,feq_procs,geq_procs,phi_procs,p_procs,u1_procs,u2_procs,u3_procs,grad_phi_procs,gphi_procs,lap_phi_procs,p0_procs,grad_u_procs,taug_procs,nu_procs,strain_procs,nue)
        call externalforce(forcex,forcey,forcez)
        call glue(forcex)
        call glue(forcey)
        call glue(forcez)
    endif
!=========================時間発展開始========================================================
    call MPI_Barrier(new_comm_world, ierr)
    time1 = MPI_Wtime()

    allocate(kinetic_procs(1:x_procs,1:ymax+1,1:z_procs))
DO n = step_start+1, step_end
!========================のりしろ境界=======================================================
    call MPI_boundary(phi_procs,u1_procs,u2_procs,u3_procs,f_procs,g_procs,taug_procs)
!=======================諸変数の計算(OMP)======================================================
    call various_cal(lap_phi_procs,grad_phi_procs,grad_u_procs,p0_procs,phi_procs,u1_procs,u2_procs,u3_procs)
    call gphi_cal(gphi_procs,grad_phi_procs)
!=======================局所平衡分布関数の計算(OMP)==============================================
    call equilibrium_cal(gphi_procs,phi_procs,p0_procs,lap_phi_procs,grad_phi_procs,grad_u_procs,p_procs,u1_procs,u2_procs,u3_procs,feq_procs,geq_procs)
    call MPI_boundary_fg(feq_procs,geq_procs)
!=======================平衡分布関数fgの計算(OMP)=================================================
    call fg_cal(fnext_procs,f_procs,feq_procs,gnext_procs,g_procs,geq_procs,taug_procs,forcex,forcey,forcez)
!==========================物理量の計算(OMP)=====================================================
    call physics(phi_procs,p_procs,u1_procs,u2_procs,u3_procs,f_procs,g_procs,nu_procs,taug_procs,fnext_procs,gnext_procs,strain_procs,grad_u_procs,nue)
!================================出力==================================
    if((mod(n,4000)==0)) then
        call output(phi_procs,u1_procs,u2_procs,u3_procs,p_procs,n)
    endif

    if(mod(n,100000)==0) then
        call outputfg(f_procs,g_procs,n)
    endif

    call MPI_Barrier(new_comm_world, ierr)
    time2 = MPI_Wtime()
    if(comm_rank == 0) then
        if(mod(n,100)==0) then
            if(n == 1) then
                open(10,file="./time_omp.d")
                write(10,*) n, time2-time1, u1_procs(1,1,1)
                close(10)
            else
                open(10,file="./time_omp.d",action="write",position="append")
                write(10,*) n, time2-time1, u1_procs(1,1,1)
                close(10)
            endif
        endif
    endif

    if(mod(n,4000)==0) then
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    kinetic_procs(xi,yi,zi) = 0.5d0 * (u1_procs(xi,yi,zi)*u1_procs(xi,yi,zi) &
                                            + u2_procs(xi,yi,zi)*u2_procs(xi,yi,zi) &
                                            + u3_procs(xi,yi,zi)*u3_procs(xi,yi,zi))
                enddo
            enddo
        enddo
        kinetic_sum = 0.0d0
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    kinetic_sum = kinetic_sum + kinetic_procs(xi,yi,zi)
                enddo
            enddo
        enddo
        call MPI_Reduce(kinetic_sum, kinetic_ave, 1, MPI_REAL8, MPI_SUM, 0, new_comm_world, ierr)
        if(new_rank == 0) then
            kinetic_ave = kinetic_ave / (dble(xmax+1)*dble(ymax+1)*dble(zmax+1))
            if(n == 4000) then
                open(11,file="./kinetic.d")
                write(11,*) dble(n), dble(n)*umax/D_vortex, kinetic_ave
                close(11)
            else
                open(11,file="./kinetic.d",action="write",position="append")
                write(11,*) dble(n), dble(n)*umax/D_vortex, kinetic_ave
                close(11)
            endif
        endif
    endif
ENDDO
!================MPI並列終わり=======================================
    call MPI_Finalize(ierr)
end program main