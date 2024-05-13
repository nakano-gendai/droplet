! 【高速計算】
! このプログラムは、MPI + OpenMP のハイブリット並列により高速計算を可能としている。
! 【系について】
! 単相のシミュレーションを行う。
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module globals
    !$ use omp_lib
    include "mpif.h"
    !計算領域
    real,parameter:: ds = 1.0d0 !格子間隔（lattice unit）
    integer,parameter:: xmax = 255 !ｘ方向格子数（０から数える）
    integer,parameter:: ymax = 255 !ｙ方向格子数（０から数える）
    integer,parameter:: zmax = 255 !ｚ方向格子数（０から数える）
    !並行して計算する数
    integer,parameter:: Nxall = 1 !x方向の分割数（Nxall*Nzall=全体の計算数）
    integer,parameter:: Nzall = 1 !z方向の分割数（Nxall*Nzall=全体の計算数）
    integer,parameter:: xall = (xmax + 1) * Nxall !全体のx方向格子数
    integer,parameter:: zall = (zmax + 1) * Nzall !全体のz方向格子数
    !時間に関するパラメータ
    integer,parameter:: step = 100000000 !計算時間step
    !出力ディレクトリ
    character(*),parameter :: datadir = "/data/sht/nakanog/test_real/"
    character(*),parameter :: datadir2 = "/data/sht/nakanog/test_real/fg/"

    !無次元数
    real,parameter:: Re = 80000.0e0 !粒子レイノルズ数

    !LESパラメータ
    real,parameter:: delta = 1.0e0*ds
    real,parameter:: Cs = 0.0e0

    !撹乱（乱数）のオーダー
    real, parameter :: ran = 1.0e-3

    !支配パラメータ
    real,parameter:: pi = acos(-1.0e0) !円周率
    real,parameter:: D_vortex = 127.5e0 !渦の大きさ
    real,parameter:: kw = pi/D_vortex !波数
    real,parameter:: umax = 0.1e0 !最大流速
    real,parameter:: nu = umax*D_vortex/Re !動粘性係数
    real,parameter:: taug = 3.0e0*nu/ds + 0.5e0

    !粒子速度（整数）
    integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
    integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
    integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)
    real :: cr(1:3, 1:15)  !粒子速度（実数）
    real krone(1:3,1:3) !クロネッカーのデルタ

    !パラメータ
    real,parameter:: E(15) = (/ 2.0e0/9.0e0, 1.0e0/9.0e0, 1.0e0/9.0e0, 1.0e0/9.0e0, 1.0e0/9.0e0, 1.0e0/9.0e0, &
                                1.0e0/9.0e0, 1.0e0/72.0e0, 1.0e0/72.0e0, 1.0e0/72.0e0, 1.0e0/72.0e0, 1.0e0/72.0e0, &
                                1.0e0/72.0e0, 1.0e0/72.0e0, 1.0e0/72.0e0 /)

    !その他変数
    real dif
    integer grobalx, grobaly, grobalz
    real dummy
    real time1, time2 !計算時間測定用
    
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
        real, intent(out):: cr(1:3,1:15), krone(1:3,1:3)
        integer i
        do i=1,15
            cr(1,i) = dble(cx(i))
            cr(2,i) = dble(cy(i))
            cr(3,i) = dble(cz(i))
        enddo
        krone(1,1) = 1.e0; krone(1,2) = 0.e0; krone(1,3) = 0.e0
        krone(2,1) = 0.e0; krone(2,2) = 1.e0; krone(2,3) = 0.e0
        krone(3,1) = 0.e0; krone(3,2) = 0.e0; krone(3,3) = 1.e0
    end subroutine par

    subroutine ini(g_procs,gnext_procs,geq_procs,p_procs,u1_procs,u2_procs,u3_procs,grad_u_procs,x_m,y_m,eps,eps2)
        real,allocatable :: g_procs(:,:,:,:), gnext_procs(:,:,:,:), geq_procs(:,:,:,:)
        real,allocatable :: p_procs(:,:,:)
        real,allocatable :: u1_procs(:,:,:), u2_procs(:,:,:), u3_procs(:,:,:)
        real,allocatable :: grad_u_procs(:,:,:,:,:)

        real,allocatable :: x_m(:,:,:), y_m(:,:,:), eps(:,:,:), eps2(:,:,:)
        integer seedsize
        integer, allocatable :: seeds(:), seeds2(:)
        integer seed, seed2

        integer i, xi, yi, zi

    !========================並列数・コミュニケータ分割・通信先設定========================
        !配列のallocate(xi:0とx_procs+1がのりしろ)(yi:0とymax+2がのりしろ)(zi:0とz_procs+1がのりしろ)
        Nx = 8 !x方向の並列数（ただし，Nx/=comm_procs）
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
        call MPI_Type_Vector(ymax+3,x_procs,x_procs+2,MPI_REAL,xtype,ierr)
        call MPI_Type_Commit(xtype,ierr)
        call MPI_Type_Vector(ymax+3,15*x_procs,15*(x_procs+2),MPI_REAL,xtype_f,ierr)
        call MPI_Type_Commit(xtype_f,ierr)

        !z世界での受け渡しをする際の型作成
        call MPI_Type_Vector((z_procs+2)*(ymax+3),1,x_procs+2,MPI_REAL,rectangle_type,ierr)
        call MPI_Type_Commit(rectangle_type,ierr)
        call MPI_Type_Vector((z_procs+2)*(ymax+3),15,15*(x_procs+2),MPI_REAL,rectangle_type_f,ierr)
        call MPI_Type_Commit(rectangle_type_f,ierr)
    !===============================================================================================================
        !以下はのりしろ有りの変数
        allocate(p_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(g_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(gnext_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(geq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs))

        allocate(x_m(1:x_procs,1:ymax+1,1:z_procs))
        allocate(y_m(1:x_procs,1:ymax+1,1:z_procs))
        allocate(eps(1:x_procs,1:ymax+1,1:z_procs))
        allocate(eps2(1:x_procs,1:ymax+1,1:z_procs))

        !初期化
        p_procs(:,:,:) = 0.0e0
        u1_procs(:,:,:) = 0.0e0
        u2_procs(:,:,:) = 0.0e0
        u3_procs(:,:,:) = 0.0e0
        geq_procs(:,:,:,:) = 0.0e0
        g_procs(:,:,:,:) = 0.0e0
        gnext_procs(:,:,:,:) = 0.0e0
        grad_u_procs(:,:,:,:,:) = 0.0e0

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
                    x_m(xi,yi,zi) = sqrt(-2.0e0*log(eps(xi,yi,zi)))*cos(2.0e0*pi*eps2(xi,yi,zi))
                    y_m(xi,yi,zi) = sqrt(-2.0e0*log(eps(xi,yi,zi)))*sin(2.0e0*pi*eps2(xi,yi,zi))
                    if(x_m(xi,yi,zi) > 3.0e0) then
                        x_m(xi,yi,zi) = 3.0e0
                    elseif(x_m(xi,yi,zi) < -3.0e0) then
                        x_m(xi,yi,zi) = -3.0e0
                    endif

                    if(y_m(xi,yi,zi) > 3.0e0) then
                        y_m(xi,yi,zi) = 3.0e0
                    elseif(y_m(xi,yi,zi) < -3.0e0) then
                        y_m(xi,yi,zi) = -3.0e0
                    endif
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
        !===================================初期条件の設定================================================
        !$omp parallel
        !$omp do
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    u1_procs(xi,yi,zi) = x_m(xi,yi,zi) * ran
                    u2_procs(xi,yi,zi) = y_m(xi,yi,zi) * ran
                    u3_procs(xi,yi,zi) = 0.5e0 * ( x_m(xi,yi,zi) + y_m(xi,yi,zi) ) * ran
                    p_procs(xi,yi,zi) = 1.0e0 / 3.0e0
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel

        call glue(u1_procs)
        call glue(u2_procs)
        call glue(u3_procs)
        
        call equilibrium_cal(p_procs,u1_procs,u2_procs,u3_procs,g_procs)
    end subroutine ini

    subroutine ini_op(kinetic_procs,forcex,forcey,forcez,strain_procs,grad_u_procs,nue,taug_procs,u1_procs,u2_procs,u3_procs)
        real,intent(inout):: grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real,intent(in):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(in):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(in):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,allocatable :: kinetic_procs(:,:,:)
        real,allocatable :: forcex(:,:,:), forcey(:,:,:), forcez(:,:,:)
        real,allocatable :: taug_procs(:,:,:)
        real,allocatable :: strain_procs(:,:,:,:,:)
        real,allocatable :: nue(:,:,:)
        integer xi, yi, zi, alpha, beta

        !以下はのりしろ無しの変数
        allocate(kinetic_procs(1:x_procs,1:ymax+1,1:z_procs))
        allocate(forcex(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(forcey(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(forcez(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(taug_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(strain_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs))
        allocate(nue(1:x_procs,1:ymax+1,1:z_procs))

        kinetic_procs(:,:,:) = 0.0e0
        forcex(:,:,:) = 0.0e0
        forcey(:,:,:) = 0.0e0
        forcez(:,:,:) = 0.0e0
        taug_procs(:,:,:) = 0.0e0
        strain_procs(:,:,:,:,:) = 0.0e0
        nue(:,:,:) = 0.0e0

        call grad_cal(grad_u_procs,u1_procs,u2_procs,u3_procs)
        call strain_cal(grad_u_procs,strain_procs)

        !$omp parallel
        !$omp do
        do zi = 1, z_procs
            do yi = 1, ymax+1
                do xi = 1, x_procs
                    do beta = 1, 3
                        do alpha = 1, 3
                            nue(xi,yi,zi) = nue(xi,yi,zi) + strain_procs(alpha,beta,xi,yi,zi)**2.0e0
                        enddo
                    enddo
                    nue(xi,yi,zi) = ( 2.0e0 * nue(xi,yi,zi) )**0.5e0
                    nue(xi,yi,zi) = Cs * Cs * delta * delta * nue(xi,yi,zi)
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
                    taug_procs(xi,yi,zi) = 0.5e0 + 3.0e0 * (nu + nue(xi,yi,zi)) / ds
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine ini_op

    subroutine glue(var)
        real,intent(inout) :: var(0:x_procs+1,0:ymax+2,0:z_procs+1)
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
        real,intent(inout) :: var_f(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
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

    subroutine MPI_boundary(u1_procs,u2_procs,u3_procs,g_procs,taug_procs)
        real,intent(inout) :: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(inout) :: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(inout) :: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(inout) :: g_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(inout) :: taug_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)

        call glue(u1_procs)
        call glue(u2_procs)
        call glue(u3_procs)
        call glue(taug_procs)
        call glue_f(g_procs)
    endsubroutine MPI_boundary

    subroutine MPI_boundary_g(geq_procs)
        real,intent(inout) :: geq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)

        call glue_f(geq_procs)
    endsubroutine MPI_boundary_g

    subroutine grad_cal(grad_u_procs,u1_procs,u2_procs,u3_procs)
        real,intent(out):: grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real,intent(in):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(in):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(in):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        integer alpha, i, xi, yi, zi

        !$omp parallel
        !$omp do
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do alpha=1,3
                        grad_u_procs(1,alpha,xi,yi,zi) = 0.0e0
                        grad_u_procs(2,alpha,xi,yi,zi) = 0.0e0
                        grad_u_procs(3,alpha,xi,yi,zi) = 0.0e0
                        do i = 2,15
                            grad_u_procs(1,alpha,xi,yi,zi) = grad_u_procs(1,alpha,xi,yi,zi) &
                                                            + cr(alpha,i)*u1_procs(xi+cx(i),yi+cy(i),zi+cz(i))
                            grad_u_procs(2,alpha,xi,yi,zi) = grad_u_procs(2,alpha,xi,yi,zi) &
                                                            + cr(alpha,i)*u2_procs(xi+cx(i),yi+cy(i),zi+cz(i))
                            grad_u_procs(3,alpha,xi,yi,zi) = grad_u_procs(3,alpha,xi,yi,zi) &
                                                            + cr(alpha,i)*u3_procs(xi+cx(i),yi+cy(i),zi+cz(i))
                        enddo
                        grad_u_procs(1,alpha,xi,yi,zi) = grad_u_procs(1,alpha,xi,yi,zi)/(10.0e0*ds)
                        grad_u_procs(2,alpha,xi,yi,zi) = grad_u_procs(2,alpha,xi,yi,zi)/(10.0e0*ds)
                        grad_u_procs(3,alpha,xi,yi,zi) = grad_u_procs(3,alpha,xi,yi,zi)/(10.0e0*ds)
                    enddo
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine grad_cal

    subroutine strain_cal(grad_u_procs,strain_procs)
        real,intent(in):: grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real,intent(out):: strain_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        integer xi, yi, zi, alpha, beta 

        !$omp parallel
        !$omp do
        do zi = 1, z_procs
            do yi = 1, ymax+1
                do xi = 1, x_procs
                    do beta = 1, 3
                        do alpha = 1, 3
                            strain_procs(alpha,beta,xi,yi,zi) = 0.5e0*( grad_u_procs(alpha,beta,xi,yi,zi) + grad_u_procs(beta,alpha,xi,yi,zi) )
                        enddo
                    enddo
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine strain_cal

    subroutine equilibrium_cal(p_procs,u1_procs,u2_procs,u3_procs,geq_procs)
        real,intent(in):: p_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(in):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(in):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(in):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(out):: geq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        integer i, xi, yi, zi

        !$omp parallel
        !$omp do
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        geq_procs(i,xi,yi,zi) = E(i)*(3.0e0*p_procs(xi,yi,zi) &
                                        + 3.0e0*(u1_procs(xi,yi,zi)*real(cx(i)) &
                                        + u2_procs(xi,yi,zi)*real(cy(i)) + u3_procs(xi,yi,zi)*real(cz(i))) &
                                        + 4.5e0*(u1_procs(xi,yi,zi)*real(cx(i)) &
                                        + u2_procs(xi,yi,zi)*real(cy(i)) + u3_procs(xi,yi,zi)*real(cz(i)))**2 &
                                        - 1.5e0*(u1_procs(xi,yi,zi)**2 &
                                        + u2_procs(xi,yi,zi)**2 + u3_procs(xi,yi,zi)**2))
                    enddo
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine equilibrium_cal

    subroutine g_cal(gnext_procs,g_procs,geq_procs,forcex,forcey,forcez,taug_procs)
        real,intent(inout):: gnext_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(in):: g_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(in):: geq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(in):: forcex(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(in):: forcey(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(in):: forcez(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(in):: taug_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        integer i, xi, yi, zi

        !$omp parallel
        !$omp do
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        gnext_procs(i,xi,yi,zi) = g_procs(i,xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                - (g_procs(i,xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                - geq_procs(i,xi-cx(i),yi-cy(i),zi-cz(i))) / taug_procs(xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                + 3.0e0 * ds * E(i) &
                                                * (real(cx(i))*forcex(xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                + real(cy(i))*forcey(xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                + real(cz(i))*forcez(xi-cx(i),yi-cy(i),zi-cz(i)))
                    enddo
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine g_cal

    subroutine physics(p_procs,u1_procs,u2_procs,u3_procs,g_procs,gnext_procs,strain_procs,grad_u_procs,nue,taug_procs)
        real,intent(inout):: p_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(inout):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(inout):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(inout):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(inout):: g_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(in) :: gnext_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(inout) :: strain_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real,intent(inout) :: grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real,intent(inout) :: nue(1:x_procs,1:ymax+1,1:z_procs)
        real,intent(inout) :: taug_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        integer i, xi, yi, zi, alpha, beta

        !$omp parallel
        !$omp do
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        g_procs(i,xi,yi,zi) = gnext_procs(i,xi,yi,zi)
                    enddo

                    p_procs(xi,yi,zi) = 0.0e0
                    u1_procs(xi,yi,zi) = 0.0e0
                    u2_procs(xi,yi,zi) = 0.0e0
                    u3_procs(xi,yi,zi) = 0.0e0
                    do i=1,15
                        p_procs(xi,yi,zi) = p_procs(xi,yi,zi) + g_procs(i,xi,yi,zi) / 3.0e0
                        u1_procs(xi,yi,zi) = u1_procs(xi,yi,zi) + real(cx(i))*g_procs(i,xi,yi,zi)
                        u2_procs(xi,yi,zi) = u2_procs(xi,yi,zi) + real(cy(i))*g_procs(i,xi,yi,zi)
                        u3_procs(xi,yi,zi) = u3_procs(xi,yi,zi) + real(cz(i))*g_procs(i,xi,yi,zi)
                    enddo
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel

        nue(:,:,:) = 0.0e0

        call grad_cal(grad_u_procs,u1_procs,u2_procs,u3_procs)
        call strain_cal(grad_u_procs,strain_procs)

        !$omp parallel
        !$omp do
        do zi = 1, z_procs
            do yi = 1, ymax+1
                do xi = 1, x_procs
                    do beta = 1, 3
                        do alpha = 1, 3
                            nue(xi,yi,zi) = nue(xi,yi,zi) + strain_procs(alpha,beta,xi,yi,zi)**2.0e0
                        enddo
                    enddo
                    nue(xi,yi,zi) = ( 2.0e0 * nue(xi,yi,zi) )**0.5e0
                    nue(xi,yi,zi) = Cs * Cs * delta * delta * nue(xi,yi,zi)
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
                    taug_procs(xi,yi,zi) = 0.5e0 + 3.0e0 * (nu + nue(xi,yi,zi)) / ds
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine physics

    subroutine externalforce(forcex,forcey,forcez)
        real,intent(inout):: forcex(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(inout):: forcey(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(inout):: forcez(0:x_procs+1,0:ymax+2,0:z_procs+1)
        integer xi, yi, zi

        !$omp parallel
        !$omp do
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    grobalx = (xi-1) + group_x * x_procs
                    grobalz = (zi-1) + group_z * z_procs

                    forcex(xi,yi,zi) = 2.0e0 * pi * pi * umax * umax &
                                        / (D_vortex * Re) * sin(kw*real(grobalx))  * cos(kw*real(grobalz))
                    forcey(xi,yi,zi) = 0.0e0
                    forcez(xi,yi,zi) = -2.0e0 * pi * pi * umax * umax &
                                        / (D_vortex * Re) * cos(kw*real(grobalx)) * sin(kw*real(grobalz))
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    endsubroutine externalforce

    subroutine output(u1_procs,u2_procs,u3_procs,p_procs,n)
        real,intent(inout):: p_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(inout):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(inout):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real,intent(inout):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
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
                    write(100) u1_procs(xi,yi,zi),u2_procs(xi,yi,zi),u3_procs(xi,yi,zi), p_procs(xi,yi,zi)
                enddo
            enddo
        enddo
        close(100)
    end subroutine output

    subroutine outputg(g_procs, n)
        real,intent(inout):: g_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
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
                        write(102) g_procs(i,xi,yi,zi)
                    enddo
                enddo
            enddo
        enddo
        close(102)
    end subroutine outputg
    
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
    real,allocatable :: g_procs(:,:,:,:), gnext_procs(:,:,:,:), geq_procs(:,:,:,:)
    real,allocatable :: p_procs(:,:,:)
    real,allocatable :: u1_procs(:,:,:), u2_procs(:,:,:), u3_procs(:,:,:)
    real,allocatable :: forcex(:,:,:), forcey(:,:,:), forcez(:,:,:)

    !LESの計算に使用する変数
    real,allocatable :: grad_u_procs(:,:,:,:,:), strain_procs(:,:,:,:,:)
    real,allocatable :: nue(:,:,:)
    real,allocatable :: taug_procs(:,:,:)

    !乱数変数
    real,allocatable :: eps(:,:,:), eps2(:,:,:), x_m(:,:,:), y_m(:,:,:)

    !解析変数
    real,allocatable ::  kinetic_procs(:,:,:)
    real kinetic_sum, kinetic_ave

    !その他変数
    integer n, xi, yi, zi

!==================================MPI並列開始===============================================
    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, comm_procs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, comm_rank, ierr)
!================================ディレクトリの作成============================================
    call mk_dirs(datadir)
    call mk_dirs(datadir2)
!================================初期条件・配列のallocate======================================
    call par(cx,cy,cz,cr,krone)
    call ini(g_procs,gnext_procs,geq_procs,p_procs,u1_procs,u2_procs,u3_procs,grad_u_procs,x_m,y_m,eps,eps2)
    call ini_op(kinetic_procs,forcex,forcey,forcez,strain_procs,grad_u_procs,nue,taug_procs,u1_procs,u2_procs,u3_procs)
    call externalforce(forcex,forcey,forcez)
    call glue(forcex)
    call glue(forcey)
    call glue(forcez)
!=========================時間発展開始========================================================
    call MPI_Barrier(new_comm_world, ierr)
    time1 = MPI_Wtime()
DO n=1,step
!========================のりしろ境界=======================================================
    call MPI_boundary(u1_procs,u2_procs,u3_procs,g_procs,taug_procs)
!=======================局所平衡分布関数の計算(OMP)==============================================
    call equilibrium_cal(p_procs,u1_procs,u2_procs,u3_procs,geq_procs)
    call MPI_boundary_g(geq_procs)
!=======================平衡分布関数fgの計算(OMP)=================================================
    call g_cal(gnext_procs,g_procs,geq_procs,forcex,forcey,forcez,taug_procs)
!==========================物理量の計算(OMP)=====================================================
    call physics(p_procs,u1_procs,u2_procs,u3_procs,g_procs,gnext_procs,strain_procs,grad_u_procs,nue,taug_procs)
!================================出力==================================
    if((mod(n,4000)==0)) then
        call output(u1_procs,u2_procs,u3_procs,p_procs,n)
    endif

    if(mod(n,100000)==0) then
        call outputg(g_procs,n)
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
        call MPI_Reduce(kinetic_sum, kinetic_ave, 1, MPI_REAL, MPI_SUM, 0, new_comm_world, ierr)
        if(new_rank == 0) then
            kinetic_ave = kinetic_ave / (real(xmax+1)*real(ymax+1)*real(zmax+1))
            if(n == 4000) then
                open(11,file="./kinetic.d")
                write(11,*) real(n), real(n)*umax/D_vortex, kinetic_ave
                close(11)
            else
                open(11,file="./kinetic.d",action="write",position="append")
                write(11,*) real(n), real(n)*umax/D_vortex, kinetic_ave
                close(11)
            endif
        endif
    endif
ENDDO
    !================MPI並列終わり=======================================
    call MPI_Finalize(ierr)
end program main