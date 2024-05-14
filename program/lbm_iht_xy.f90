! 【高速計算】
! このプログラムは、MPI + OpenMP のハイブリット並列により高速計算を可能としている。
! 【系について】
! 単相のシミュレーションを行う。
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module globals
    !$ use omp_lib
    use decomp_2d
    use decomp_2d_fft
    use decomp_2d_io
    use glassman
    include "mpif.h"
    include 'fftw3.f'
    !計算領域
    real(8),parameter:: ds = 1.0d0 !格子間隔（lattice unit）
    integer,parameter:: xmax = 255 !ｘ方向格子数（０から数える）
    integer,parameter:: ymax = 255 !ｙ方向格子数（０から数える）
    integer,parameter:: zmax = 255 !ｚ方向格子数（０から数える）
    !並行して計算する数
    integer,parameter:: Nxall = 1 !x方向の分割数（Nxall*Nzall=全体の計算数）
    integer,parameter:: Nyall = 1 !z方向の分割数（Nxall*Nzall=全体の計算数）
    integer,parameter:: xall = (xmax + 1) * Nxall !全体のx方向格子数
    integer,parameter:: yall = (ymax + 1) * Nyall !全体のz方向格子数
    !時間に関するパラメータ
    integer,parameter:: step = 10000000 !計算時間step
    !出力ディレクトリ
    character(*),parameter :: datadir = "/data/sht/nakanog/DNS_turbulence_256_IHT_7/"
    character(*),parameter :: datadir2 = "/data/sht/nakanog/DNS_turbulence_256_IHT_7/fg/"
    integer,parameter:: step_output = 4000
    integer,parameter:: step_putput_fg = 100000

    !LESパラメータ
    real(8),parameter:: delta = 0.0d0*ds
    real(8),parameter:: Cs = 0.0d0

    !ローパスフィルター
    integer,parameter:: k_low = 8

    !撹乱（乱数）のオーダー
    real(8), parameter :: ran = 1.0d-4

    !支配パラメータ
    real(8),parameter:: pi = acos(-1.0d0) !円周率

    !カットオフ波数（IHT）
    integer,parameter:: kc = 3
    !動粘性係数（IHT）
    real(8),parameter:: nu = 0.001d0
    !緩和時間（IHT）
    real(8),parameter:: taug = 3.0d0*nu/ds + 0.5d0
    !Kolmogorovスケール (IHT)
    real(8),parameter:: Kolmogorov_scale = 1.5d0 * ds
    !エネルギー散逸率（IHT）
    real(8),parameter:: epsilon = (nu**3.0d0) / (Kolmogorov_scale**4.0d0)

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

    !その他変数
    real(8) dif
    integer grobalx, grobaly, grobalz
    real(8) dummy
    real(8) time1, time2 !計算時間測定用
    
    character :: filename*200
    character :: filename2*200
    character :: filename3*200

    !MPI用変数
    integer ierr, comm_procs, comm_rank
    integer comm_rank_x, comm_rank_y
    integer Nx, Ny 
    integer x_procs, y_procs
    integer key_new, group_new, key_x, group_x, key_y, group_y
    integer new_comm_world, new_procs, new_rank
    integer newx_comm_world, newx_procs, newx_rank
    integer newy_comm_world, newy_procs, newy_rank
    integer rectangle_type, rectangle_type_f, xtype, xtype_f
    integer next_rank_x, former_rank_x, next_rank_y, former_rank_y
    integer req1s, req1r, req2s, req2r
    integer, dimension(MPI_STATUS_SIZE) :: sta1s, sta1r, sta2s, sta2r
    integer tags, tagr, recv_rank
    character(2) chmyrank
    !2decomp用変数
    integer, save:: sta(1:3), last(1:3), sized(1:3) 

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

    subroutine ini(g_procs,gnext_procs,geq_procs,p_procs,u1_procs,u2_procs,u3_procs,grad_u_procs,x_m,y_m,eps,eps2)
        real(8),allocatable :: g_procs(:,:,:,:), gnext_procs(:,:,:,:), geq_procs(:,:,:,:)
        real(8),allocatable :: p_procs(:,:,:)
        real(8),allocatable :: u1_procs(:,:,:), u2_procs(:,:,:), u3_procs(:,:,:)
        real(8),allocatable :: grad_u_procs(:,:,:,:,:)

        real(8),allocatable :: x_m(:,:,:), y_m(:,:,:), eps(:,:,:), eps2(:,:,:)
        integer seedsize
        integer, allocatable :: seeds(:), seeds2(:)
        integer seed, seed2

        integer i, xi, yi, zi

    !========================並列数・コミュニケータ分割・通信先設定========================
        !配列のallocate(xi:0とx_procs+1がのりしろ)(yi:0とymax+2がのりしろ)(zi:0とz_procs+1がのりしろ)
        Nx = 32 !x方向の並列数（ただし，Nx/=comm_procs）
        Ny = comm_procs / (Nx * Nxall * Nyall) !y方向の並列数
        x_procs = (xmax+1) / Nx
        y_procs = (ymax+1) / Ny

        comm_rank_x = comm_rank / Ny
        comm_rank_y = comm_rank - comm_rank_x * Ny

        !X方向プロセスの通信番号の取得
        if((comm_rank >= 0) .and. (comm_rank <= Ny - 1)) then
            former_rank_x = comm_rank + Ny*(Nx-1)
        else
            former_rank_x = comm_rank - Ny
        endif
        if((comm_rank >= comm_procs - Ny) .and. (comm_rank <= comm_procs -1)) then
            next_rank_x = comm_rank - Ny*(Nx-1)
        else
            next_rank_x = comm_rank + Ny
        end if

        !!!Y方向プロセスの通信番号の取得
        if(mod(comm_rank,Ny) == 0) then
            former_rank_y = comm_rank + (Ny - 1)
        else
            former_rank_y = comm_rank - 1
        endif
        if(mod(comm_rank,Ny) == Ny - 1) then
            next_rank_y = comm_rank - (Ny - 1)
        else
            next_rank_y = comm_rank + 1
        end if

        !x世界での受け渡しをする際の型作成
        call MPI_Type_Vector((y_procs+2)*(zmax+3),1,x_procs+2,MPI_REAL8,xtype,ierr)
        call MPI_Type_Commit(xtype,ierr)
        call MPI_Type_Vector((y_procs+2)*(zmax+3),15,15*(x_procs+2),MPI_REAL8,xtype_f,ierr)
        call MPI_Type_Commit(xtype_f,ierr)

        !y世界での受け渡しをする際の型作成
        call MPI_Type_Vector(zmax+3,x_procs+2,(x_procs+2)*(y_procs+2),MPI_REAL8,rectangle_type,ierr)
        call MPI_Type_Commit(rectangle_type,ierr)
        call MPI_Type_Vector(zmax+3,15*(x_procs+2),15*(x_procs+2)*(y_procs+2),MPI_REAL8,rectangle_type_f,ierr)
        call MPI_Type_Commit(rectangle_type_f,ierr)
    !===============================================================================================================
        !以下はのりしろ有りの変数
        allocate(p_procs(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(g_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(gnext_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(geq_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(grad_u_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1))

        allocate(x_m(1:x_procs,1:y_procs,1:zmax+1))
        allocate(y_m(1:x_procs,1:y_procs,1:zmax+1))
        allocate(eps(1:x_procs,1:y_procs,1:zmax+1))
        allocate(eps2(1:x_procs,1:y_procs,1:zmax+1))

        !初期化
        p_procs(:,:,:) = 0.0d0
        u1_procs(:,:,:) = 0.0d0
        u2_procs(:,:,:) = 0.0d0
        u3_procs(:,:,:) = 0.0d0
        geq_procs(:,:,:,:) = 0.0d0
        g_procs(:,:,:,:) = 0.0d0
        gnext_procs(:,:,:,:) = 0.0d0
        grad_u_procs(:,:,:,:,:) = 0.0d0

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
        do zi=1,zmax+1
            do yi=1,y_procs
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
        !===================================初期条件の設定================================================
        !$omp parallel
        !$omp do
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    grobalx = (xi-1) + comm_rank_x * x_procs
                    grobaly = (yi-1) + comm_rank_y * y_procs
                    grobalz = (zi-1)

                    u1_procs(xi,yi,zi) = -0.0001d0*cos(2.0d0*pi/dble(xmax+1)*dble(grobalx))*sin(2.0d0*pi/dble(zmax+1)*dble(grobalz))
                    u2_procs(xi,yi,zi) = 0.0d0
                    u3_procs(xi,yi,zi) = (2.0d0*pi/dble(xmax+1))/(2.0d0*pi/dble(zmax+1))*0.0001d0*sin(2.0d0*pi/dble(xmax+1)*dble(grobalx))*cos(2.0d0*pi/dble(zmax+1)*dble(grobalz))
                    p_procs(xi,yi,zi) = 1.0d0/3.0d0 - 0.0001d0**2 / 4.0d0 * (cos(2.0d0*2.0d0*pi/dble(xmax+1)*dble(grobalx)) + &
                                        (2.0d0*pi/dble(xmax+1))**2 / (2.0d0*pi/dble(zmax+1))**2 * cos(2.0d0*2.0d0*pi/dble(zmax+1)*dble(grobalz)))

                    ! u1_procs(xi,yi,zi) = x_m(xi,yi,zi) * ran
                    ! u2_procs(xi,yi,zi) = y_m(xi,yi,zi) * ran
                    ! u3_procs(xi,yi,zi) = 0.5d0 * ( x_m(xi,yi,zi) + y_m(xi,yi,zi) ) * ran
                    ! p_procs(xi,yi,zi) = 1.0d0 / 3.0d0
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
        real(8),intent(inout):: grad_u_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),allocatable :: kinetic_procs(:,:,:)
        real(8),allocatable :: forcex(:,:,:), forcey(:,:,:), forcez(:,:,:)
        real(8),allocatable :: taug_procs(:,:,:)
        real(8),allocatable :: strain_procs(:,:,:,:,:)
        real(8),allocatable :: nue(:,:,:)
        integer xi, yi, zi, alpha, beta

        !以下はのりしろ無しの変数
        allocate(kinetic_procs(1:x_procs,1:y_procs,1:zmax+1))
        allocate(forcex(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(forcey(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(forcez(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(taug_procs(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(strain_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1))
        allocate(nue(1:x_procs,1:y_procs,1:zmax+1))

        kinetic_procs(:,:,:) = 0.0d0
        forcex(:,:,:) = 0.0d0
        forcey(:,:,:) = 0.0d0
        forcez(:,:,:) = 0.0d0
        taug_procs(:,:,:) = 0.0d0
        strain_procs(:,:,:,:,:) = 0.0d0
        nue(:,:,:) = 0.0d0

        call grad_cal(grad_u_procs,u1_procs,u2_procs,u3_procs)
        call strain_cal(grad_u_procs,strain_procs)

        !$omp parallel
        !$omp do
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    do beta=1,3
                        do alpha=1,3
                            nue(xi,yi,zi) = nue(xi,yi,zi) + strain_procs(alpha,beta,xi,yi,zi)**2.0d0
                        enddo
                    enddo
                    nue(xi,yi,zi) = ( 2.0d0 * nue(xi,yi,zi) )**0.5d0
                    nue(xi,yi,zi) = Cs * Cs * delta * delta * nue(xi,yi,zi)
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel


        !$omp parallel
        !$omp do
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    taug_procs(xi,yi,zi) = 0.5d0 + 3.0d0 * (nu + nue(xi,yi,zi)) / ds
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine ini_op

    subroutine ini_fft(u1_hat,u2_hat,u3_hat,u1_procs_re,u2_procs_re,u3_procs_re,forcex_hat,forcey_hat,forcez_hat,forcex_re,forcey_re,forcez_re)
        complex(kind(0d0)), allocatable :: u1_hat(:,:,:), u2_hat(:,:,:), u3_hat(:,:,:)
        real(8),allocatable :: u1_procs_re(:,:,:), u2_procs_re(:,:,:), u3_procs_re(:,:,:)
        complex(kind(0d0)), allocatable :: forcex_hat(:,:,:), forcey_hat(:,:,:), forcez_hat(:,:,:)
        real(8),allocatable :: forcex_re(:,:,:), forcey_re(:,:,:), forcez_re(:,:,:)

        call decomp_2d_init(xmax+1, ymax+1, zmax+1, Nx, Ny)
        call decomp_2d_fft_init(PHYSICAL_IN_Z)
        call decomp_2d_fft_get_size(sta, last, sized)

        allocate(u1_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax))
        allocate(u2_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax))
        allocate(u3_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax))
        allocate(u1_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1))
        allocate(u2_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1))
        allocate(u3_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1))

        allocate(forcex_re(0:x_procs-1, 0:y_procs-1, 0:zmax))
        allocate(forcey_re(0:x_procs-1, 0:y_procs-1, 0:zmax))
        allocate(forcez_re(0:x_procs-1, 0:y_procs-1, 0:zmax))
        allocate(forcex_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1))
        allocate(forcey_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1))
        allocate(forcez_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1))

        u1_procs_re(:,:,:) = 0.0d0
        u2_procs_re(:,:,:) = 0.0d0
        u3_procs_re(:,:,:) = 0.0d0
        u1_hat(:,:,:) = (0.0d0, 0.0d0)
        u2_hat(:,:,:) = (0.0d0, 0.0d0)
        u3_hat(:,:,:) = (0.0d0, 0.0d0)

        forcex_re(:,:,:) = 0.0d0
        forcey_re(:,:,:) = 0.0d0
        forcez_re(:,:,:) = 0.0d0
        forcex_hat(:,:,:) = (0.0d0, 0.0d0)
        forcey_hat(:,:,:) = (0.0d0, 0.0d0)
        forcez_hat(:,:,:) = (0.0d0, 0.0d0)
    end subroutine ini_fft

    subroutine glue(var)
        real(8),intent(inout) :: var(0:x_procs+1,0:y_procs+1,0:zmax+2)
        integer xi, yi

        do yi=1,y_procs
            do xi=1,x_procs
                var(xi,yi,0) = var(xi,yi,zmax+1)
                var(xi,yi,zmax+2) = var(xi,yi,1)
            enddo
        enddo

        !X世界でののりしろ通信
        call MPI_Isend(var(x_procs,0,0),1,xtype,next_rank_x,1,MPI_COMM_WORLD,req1s,ierr)
        call MPI_Irecv(var(0,0,0),1,xtype,former_rank_x,1,MPI_COMM_WORLD,req1r,ierr)

        call MPI_Isend(var(1,0,0),1,xtype,former_rank_x,2,MPI_COMM_WORLD,req2s,ierr)
        call MPI_Irecv(var(x_procs+1,0,0),1,xtype,next_rank_x,2,MPI_COMM_WORLD,req2r,ierr)

        call MPI_Wait(req1s,sta1s,ierr)
        call MPI_Wait(req1r,sta1r,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
        call MPI_Wait(req2r,sta2r,ierr)
        !Y世界でののりしろ通信
        call MPI_Isend(var(0,y_procs,0),1,rectangle_type,next_rank_y,1,MPI_COMM_WORLD,req1s,ierr)
        call MPI_Irecv(var(0,0,0),1,rectangle_type,former_rank_y,1,MPI_COMM_WORLD,req1r,ierr)

        call MPI_Isend(var(0,1,0),1,rectangle_type,former_rank_y,2,MPI_COMM_WORLD,req2s,ierr)
        call MPI_Irecv(var(0,y_procs+1,0),1,rectangle_type,next_rank_y,2,MPI_COMM_WORLD,req2r,ierr)

        call MPI_Wait(req1s,sta1s,ierr)
        call MPI_Wait(req1r,sta1r,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
        call MPI_Wait(req2r,sta2r,ierr)
    end subroutine glue

    subroutine glue_f(var_f)
        real(8),intent(inout) :: var_f(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        integer i, xi, yi

        do yi=1,y_procs
            do xi=1,x_procs
                do i=1,15
                    var_f(i,xi,yi,0) = var_f(i,xi,yi,zmax+1)
                    var_f(i,xi,yi,zmax+2) = var_f(i,xi,yi,1)
                enddo
            enddo
        enddo

        !X世界でののりしろ通信
        call MPI_Isend(var_f(1,x_procs,0,0),1,xtype_f,next_rank_x,1,MPI_COMM_WORLD,req1s,ierr)
        call MPI_Irecv(var_f(1,0,0,0),1,xtype_f,former_rank_x,1,MPI_COMM_WORLD,req1r,ierr)

        call MPI_Isend(var_f(1,1,0,0),1,xtype_f,former_rank_x,2,MPI_COMM_WORLD,req2s,ierr)
        call MPI_Irecv(var_f(1,x_procs+1,0,0),1,xtype_f,next_rank_x,2,MPI_COMM_WORLD,req2r,ierr)

        call MPI_Wait(req1s,sta1s,ierr)
        call MPI_Wait(req1r,sta1r,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
        call MPI_Wait(req2r,sta2r,ierr)

        !Y世界でののりしろ通信
        call MPI_Isend(var_f(1,0,y_procs,0),1,rectangle_type_f,next_rank_y,1,MPI_COMM_WORLD,req1s,ierr)
        call MPI_Irecv(var_f(1,0,0,0),1,rectangle_type_f,former_rank_y,1,MPI_COMM_WORLD,req1r,ierr)

        call MPI_Isend(var_f(1,0,1,0),1,rectangle_type_f,former_rank_y,2,MPI_COMM_WORLD,req2s,ierr)
        call MPI_Irecv(var_f(1,0,y_procs+1,0),1,rectangle_type_f,next_rank_y,2,MPI_COMM_WORLD,req2r,ierr)

        call MPI_Wait(req1s,sta1s,ierr)
        call MPI_Wait(req1r,sta1r,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
        call MPI_Wait(req2r,sta2r,ierr)
    end subroutine glue_f

    subroutine MPI_boundary(u1_procs,u2_procs,u3_procs,g_procs,taug_procs,forcex,forcey,forcez)
        real(8),intent(inout) :: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout) :: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout) :: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout) :: g_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout) :: taug_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout) :: forcex(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout) :: forcey(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout) :: forcez(0:x_procs+1,0:y_procs+1,0:zmax+2)

        call glue(u1_procs)
        call glue(u2_procs)
        call glue(u3_procs)
        call glue(taug_procs)
        call glue(forcex)
        call glue(forcey)
        call glue(forcez)
        call glue_f(g_procs)
    endsubroutine MPI_boundary

    subroutine MPI_boundary_g(geq_procs)
        real(8),intent(inout) :: geq_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)

        call glue_f(geq_procs)
    endsubroutine MPI_boundary_g

    subroutine grad_cal(grad_u_procs,u1_procs,u2_procs,u3_procs)
        real(8),intent(out):: grad_u_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        integer alpha, i, xi, yi, zi

        !$omp parallel
        !$omp do
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    do alpha=1,3
                        grad_u_procs(1,alpha,xi,yi,zi) = 0.0d0
                        grad_u_procs(2,alpha,xi,yi,zi) = 0.0d0
                        grad_u_procs(3,alpha,xi,yi,zi) = 0.0d0
                        do i = 2,15
                            grad_u_procs(1,alpha,xi,yi,zi) = grad_u_procs(1,alpha,xi,yi,zi) &
                                                            + cr(alpha,i)*u1_procs(xi+cx(i),yi+cy(i),zi+cz(i))
                            grad_u_procs(2,alpha,xi,yi,zi) = grad_u_procs(2,alpha,xi,yi,zi) &
                                                            + cr(alpha,i)*u2_procs(xi+cx(i),yi+cy(i),zi+cz(i))
                            grad_u_procs(3,alpha,xi,yi,zi) = grad_u_procs(3,alpha,xi,yi,zi) &
                                                            + cr(alpha,i)*u3_procs(xi+cx(i),yi+cy(i),zi+cz(i))
                        enddo
                        grad_u_procs(1,alpha,xi,yi,zi) = grad_u_procs(1,alpha,xi,yi,zi)/(10.0d0*ds)
                        grad_u_procs(2,alpha,xi,yi,zi) = grad_u_procs(2,alpha,xi,yi,zi)/(10.0d0*ds)
                        grad_u_procs(3,alpha,xi,yi,zi) = grad_u_procs(3,alpha,xi,yi,zi)/(10.0d0*ds)
                    enddo
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine grad_cal

    subroutine strain_cal(grad_u_procs,strain_procs)
        real(8),intent(in):: grad_u_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(out):: strain_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)
        integer xi, yi, zi, alpha, beta 

        !$omp parallel
        !$omp do
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    do beta=1,3
                        do alpha=1,3
                            strain_procs(alpha,beta,xi,yi,zi) = 0.5d0*( grad_u_procs(alpha,beta,xi,yi,zi) + grad_u_procs(beta,alpha,xi,yi,zi) )
                        enddo
                    enddo
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine strain_cal

    subroutine equilibrium_cal(p_procs,u1_procs,u2_procs,u3_procs,geq_procs)
        real(8),intent(in):: p_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(out):: geq_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        integer i, xi, yi, zi

        !$omp parallel
        !$omp do
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    do i=1,15
                        geq_procs(i,xi,yi,zi) = E(i)*(3.0d0*p_procs(xi,yi,zi) &
                                        + 3.0d0*(u1_procs(xi,yi,zi)*dble(cx(i)) &
                                        + u2_procs(xi,yi,zi)*dble(cy(i)) + u3_procs(xi,yi,zi)*dble(cz(i))) &
                                        + 4.5d0*(u1_procs(xi,yi,zi)*dble(cx(i)) &
                                        + u2_procs(xi,yi,zi)*dble(cy(i)) + u3_procs(xi,yi,zi)*dble(cz(i)))**2 &
                                        - 1.5d0*(u1_procs(xi,yi,zi)**2 &
                                        + u2_procs(xi,yi,zi)**2 + u3_procs(xi,yi,zi)**2))
                    enddo
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine equilibrium_cal

    subroutine g_cal(gnext_procs,g_procs,geq_procs,forcex,forcey,forcez,taug_procs)
        real(8),intent(inout):: gnext_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: g_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: geq_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: forcex(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: forcey(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: forcez(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: taug_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        integer i, xi, yi, zi

        !$omp parallel
        !$omp do
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    do i=1,15
                        gnext_procs(i,xi,yi,zi) = g_procs(i,xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                - (g_procs(i,xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                - geq_procs(i,xi-cx(i),yi-cy(i),zi-cz(i))) / taug_procs(xi-cx(i),yi-cy(i),zi-cz(i)) &
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
    end subroutine g_cal

subroutine physics(p_procs,u1_procs,u2_procs,u3_procs,g_procs,gnext_procs,strain_procs,grad_u_procs,nue,taug_procs)
        real(8),intent(inout):: p_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: g_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in) :: gnext_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout) :: strain_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(inout) :: grad_u_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(inout) :: nue(1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(inout) :: taug_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        integer i, xi, yi, zi, alpha, beta

        !$omp parallel
        !$omp do
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    ! do i=1,15
                    !     g_procs(i,xi,yi,zi) = gnext_procs(i,xi,yi,zi)
                    ! enddo

                    p_procs(xi,yi,zi) = 0.0d0
                    u1_procs(xi,yi,zi) = 0.0d0
                    u2_procs(xi,yi,zi) = 0.0d0
                    u3_procs(xi,yi,zi) = 0.0d0
                    do i=1,15
                        g_procs(i,xi,yi,zi) = gnext_procs(i,xi,yi,zi)
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

        nue(:,:,:) = 0.0d0

        ! call grad_cal(grad_u_procs,u1_procs,u2_procs,u3_procs)
        ! call strain_cal(grad_u_procs,strain_procs)

        !$omp parallel
        !$omp do
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    do beta=1,3
                        do alpha=1,3
                            nue(xi,yi,zi) = nue(xi,yi,zi) + strain_procs(alpha,beta,xi,yi,zi)**2.0d0
                        enddo
                    enddo
                    nue(xi,yi,zi) = ( 2.0d0 * nue(xi,yi,zi) )**0.5d0
                    nue(xi,yi,zi) = Cs * Cs * delta * delta * nue(xi,yi,zi)
                    taug_procs(xi,yi,zi) = 0.5d0 + 3.0d0 * (nu + nue(xi,yi,zi)) / ds
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine physics

    subroutine externalforce(forcex,forcey,forcez,forcex_re,forcey_re,forcez_re,forcex_hat,forcey_hat,forcez_hat,u1_procs,u2_procs,u3_procs,u1_procs_re,u2_procs_re,u3_procs_re,u1_hat,u2_hat,u3_hat)
        real(8),intent(inout):: forcex(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: forcey(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: forcez(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: forcex_re(0:x_procs-1, 0:y_procs-1, 0:zmax)
        real(8),intent(inout):: forcey_re(0:x_procs-1, 0:y_procs-1, 0:zmax)
        real(8),intent(inout):: forcez_re(0:x_procs-1, 0:y_procs-1, 0:zmax)
        real(8),intent(inout):: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u1_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax)
        real(8),intent(inout):: u2_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax)
        real(8),intent(inout):: u3_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax)
        complex(kind(0d0)), intent(inout) :: forcex_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)
        complex(kind(0d0)), intent(inout) :: forcey_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)
        complex(kind(0d0)), intent(inout) :: forcez_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)

        complex(kind(0d0)), intent(inout) :: u1_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)
        complex(kind(0d0)), intent(inout) :: u2_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)
        complex(kind(0d0)), intent(inout) :: u3_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)
        integer xi, yi, zi, i
        integer k1, k2, k3, k(1:3), k_index
        real(8) k_abs
        real(8) dkx, dky, dkz, dk
        real(8) enesupe((xmax+1)/2 + 1), enesupe_sum((xmax+1)/2 + 1)
        real(8) energy, energy_procs

        dkx = 2.0d0*pi/dble(xmax+1)
        dky = 2.0d0*pi/dble(ymax+1)
        dkz = 2.0d0*pi/dble(zmax+1)
        dk = ( (dkx)**2.0d0 + (dky)**2.0d0 + (dkz)**2.0d0 )**0.5d0

        !u_hatの配列数と合わせる
        !$omp parallel
        !$omp do
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    u1_procs_re(xi-1,yi-1,zi-1) = u1_procs(xi,yi,zi)
                    u2_procs_re(xi-1,yi-1,zi-1) = u2_procs(xi,yi,zi)
                    u3_procs_re(xi-1,yi-1,zi-1) = u3_procs(xi,yi,zi)
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel

        !フーリエ変換（実数→複素数）
        call fft_r2c(u1_procs_re, u1_hat)
        call fft_r2c(u2_procs_re, u2_hat)
        call fft_r2c(u3_procs_re, u3_hat)

        enesupe(:) = 0.0d0
        enesupe_sum(:) = 0.0d0
        energy_procs = 0.0d0
        do k3=sta(3)-1,last(3)-1
            k(3) =( k3 - judge(k3, zmax+1) )
            do k2=sta(2)-1,last(2)-1
                k(2) = ( k2 - judge(k2, ymax+1) )
                do k1=sta(1)-1,last(1)-1
                    k(1) = ( k1 - judge(k1, xmax+1) )

                    k_abs = ( dble(k(1)**2 + k(2)**2 + k(3)**2) )**0.5d0
                    k_index = int( k_abs ) + 1

                    if((k_index > 1) .and. (k_index < kc + 1)) then
                        energy_procs = energy_procs + 0.5d0 * (abs(u1_hat(k1,k2,k3))**2 + abs(u2_hat(k1,k2,k3))**2 + abs(u3_hat(k1,k2,k3))**2)
                    endif
                enddo
            enddo
        enddo

        energy = 0.0d0
        call MPI_Allreduce(energy_procs, energy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        do k3=sta(3)-1,last(3)-1
            k(3) = k3 - judge(k3, zmax+1)
            do k2=sta(2)-1,last(2)-1
                k(2) = k2 - judge(k2, ymax+1)
                do k1=sta(1)-1,last(1)-1
                    k(1) = k1 - judge(k1, xmax+1)

                    k_abs = ( dble(k(1))**2.0d0 + dble(k(2))**2.0d0 + dble(k(3))**2.0d0 )**0.5d0
                    k_index = int( k_abs ) + 1

                    if((k_index > 1) .and. (k_index < kc + 1)) then
                        forcex_hat(k1,k2,k3) = 0.5d0 * epsilon / energy * u1_hat(k1,k2,k3)
                        forcey_hat(k1,k2,k3) = 0.5d0 * epsilon / energy * u2_hat(k1,k2,k3)
                        forcez_hat(k1,k2,k3) = 0.5d0 * epsilon / energy * u3_hat(k1,k2,k3)
                    else
                        forcex_hat(k1,k2,k3) = (0.0d0, 0.0d0)
                        forcey_hat(k1,k2,k3) = (0.0d0, 0.0d0)
                        forcez_hat(k1,k2,k3) = (0.0d0, 0.0d0)
                    endif
                enddo
            enddo
        enddo

        forcex_re(:,:,:) = 0.0d0
        forcey_re(:,:,:) = 0.0d0
        forcez_re(:,:,:) = 0.0d0
        !逆フーリエ変換（複素数→実数）
        call fft_c2r(forcex_hat, forcex_re)
        call fft_c2r(forcey_hat, forcey_re)
        call fft_c2r(forcez_hat, forcez_re)
        !元の速度の配列に戻す
        !$omp parallel
        !$omp do
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    forcex(xi,yi,zi) = forcex_re(xi-1,yi-1,zi-1)
                    forcey(xi,yi,zi) = forcey_re(xi-1,yi-1,zi-1)
                    forcez(xi,yi,zi) = forcez_re(xi-1,yi-1,zi-1)
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    endsubroutine externalforce

    subroutine kinetic_cal(kinetic_procs,u1_procs,u2_procs,u3_procs,n)
        real(8),intent(inout):: kinetic_procs(1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        integer,intent(in):: n
        real(8) kinetic_sum, kinetic_ave
        integer xi, yi, zi

        !$omp parallel
        !$omp do
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    kinetic_procs(xi,yi,zi) = 0.5d0 * (u1_procs(xi,yi,zi)*u1_procs(xi,yi,zi) &
                                            + u2_procs(xi,yi,zi)*u2_procs(xi,yi,zi) &
                                            + u3_procs(xi,yi,zi)*u3_procs(xi,yi,zi))
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
        kinetic_sum = 0.0d0
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    kinetic_sum = kinetic_sum + kinetic_procs(xi,yi,zi)
                enddo
            enddo
        enddo
        call MPI_Reduce(kinetic_sum, kinetic_ave, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        !出力
        if(comm_rank == 0) then
            kinetic_ave = kinetic_ave / (dble(xmax+1)*dble(ymax+1)*dble(zmax+1))
            if(n == step_output) then
                open(11,file="./kinetic.d")
                write(11,*) dble(n), kinetic_ave
                close(11)
            else
                open(11,file="./kinetic.d",action="write",position="append")
                write(11,*) dble(n), kinetic_ave
                close(11)
            endif
        endif
    end subroutine kinetic_cal

    subroutine output(u1_procs,u2_procs,u3_procs,p_procs,n)
        real(8),intent(inout):: p_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        integer,intent(in):: n
        integer xi, yi, zi

        ! write(filename,*) group_new !i->filename 変換
        write(filename2,*) n
        write(filename3,*) comm_rank
        ! filename=datadir//trim(adjustl(filename))//'_'//trim(adjustl(filename3))//'_'//trim(adjustl(filename2))//'.bin' 
        filename=datadir//'0_'//trim(adjustl(filename3))//'_'//trim(adjustl(filename2))//'.bin' 
        print *, filename !表示してみる
        open(100,file=filename, form='unformatted',status='replace')
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    write(100) u1_procs(xi,yi,zi), u2_procs(xi,yi,zi), u3_procs(xi,yi,zi), p_procs(xi,yi,zi)
                enddo
            enddo
        enddo
        close(100)
    end subroutine output

    subroutine outputg(g_procs, n)
        real(8),intent(inout):: g_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        integer,intent(in):: n
        integer i, xi, yi, zi

        ! write(filename,*) group_new !i->filename 変換
        write(filename2,*) n
        write(filename3,*) comm_rank
        filename=datadir2//'0_'//trim(adjustl(filename3))//'_'//trim(adjustl(filename2))//'_fg.bin' 
        print *, filename !表示してみる
        open(102,file=filename, form='unformatted',status='replace') 
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    do i=1,15
                        write(102) g_procs(i,xi,yi,zi)
                    enddo
                enddo
            enddo
        enddo
        close(102)
    end subroutine outputg

    !フーリエ変換
    subroutine fft_r2c(q, q_hat)
        real(8), intent(inout)            :: q(0:x_procs-1, 0:y_procs-1, 0:zmax)
        complex(kind(0d0)), intent(inout) :: q_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)
        call decomp_2d_fft_3d(q, q_hat)
        q_hat(:,:,:) = q_hat(:,:,:) / (dble(xmax+1)*dble(ymax+1)*dble(zmax+1))
    end subroutine fft_r2c

    !逆フーリエ変換
    subroutine fft_c2r(q_hat, q)
        real(8), intent(inout)            :: q(0:x_procs-1, 0:y_procs-1, 0:zmax)
        complex(kind(0d0)), intent(inout) :: q_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)
        call decomp_2d_fft_3d(q_hat, q)
    end subroutine fft_c2r

    !ディレクトリ作成
    subroutine mk_dirs(outdir)
        implicit none
        character(len=*), intent(in) :: outdir
        character(len=256) command
        write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
        call system(command)
    end subroutine mk_dirs

    function judge(i, N)
        integer :: i, N, judge
        if(i >= N / 2 + 1)then
            judge = N
        else
            judge = 0
        endif
    end function

end module globals

program main
use globals
use decomp_2d
use decomp_2d_fft
use decomp_2d_io
use glassman
    implicit none
    !LBMの計算に使用する変数
    real(8),allocatable :: g_procs(:,:,:,:), gnext_procs(:,:,:,:), geq_procs(:,:,:,:)
    real(8),allocatable :: p_procs(:,:,:)
    real(8),allocatable :: u1_procs(:,:,:), u2_procs(:,:,:), u3_procs(:,:,:)
    real(8),allocatable :: forcex(:,:,:), forcey(:,:,:), forcez(:,:,:)

    !LESの計算に使用する変数
    real(8),allocatable :: grad_u_procs(:,:,:,:,:), strain_procs(:,:,:,:,:)
    real(8),allocatable :: nue(:,:,:)
    real(8),allocatable :: taug_procs(:,:,:)

    !乱数変数
    real(8),allocatable :: eps(:,:,:), eps2(:,:,:), x_m(:,:,:), y_m(:,:,:)

    !解析変数
    real(8),allocatable ::  kinetic_procs(:,:,:)

    !その他変数
    integer n

    !波数空間の変数
    complex(kind(0d0)), allocatable :: u1_hat(:,:,:), u2_hat(:,:,:), u3_hat(:,:,:)
    real(8),allocatable :: u1_procs_re(:,:,:), u2_procs_re(:,:,:), u3_procs_re(:,:,:)
    complex(kind(0d0)), allocatable :: forcex_hat(:,:,:), forcey_hat(:,:,:), forcez_hat(:,:,:)
    real(8),allocatable :: forcex_re(:,:,:), forcey_re(:,:,:), forcez_re(:,:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !MPI並列開始
    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, comm_procs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, comm_rank, ierr)
    !ディレクトリの作成
    call mk_dirs(datadir)
    call mk_dirs(datadir2)
    !初期条件・配列のallocate
    call par(cx,cy,cz,cr,krone)
    call ini(g_procs,gnext_procs,geq_procs,p_procs,u1_procs,u2_procs,u3_procs,grad_u_procs,x_m,y_m,eps,eps2)
    call ini_op(kinetic_procs,forcex,forcey,forcez,strain_procs,grad_u_procs,nue,taug_procs,u1_procs,u2_procs,u3_procs)
    call ini_fft(u1_hat,u2_hat,u3_hat,u1_procs_re,u2_procs_re,u3_procs_re,forcex_hat,forcey_hat,forcez_hat,forcex_re,forcey_re,forcez_re)
    call externalforce(forcex,forcey,forcez,forcex_re,forcey_re,forcez_re,forcex_hat,forcey_hat,forcez_hat, &
    u1_procs,u2_procs,u3_procs,u1_procs_re,u2_procs_re,u3_procs_re,u1_hat,u2_hat,u3_hat)
    call glue(forcex)
    call glue(forcey)
    call glue(forcez)
!!!!!!!!!!!!!!!!!!!!!!!時間発展開始!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    time1 = MPI_Wtime()
DO n=1,step
    !外力
    call externalforce(forcex,forcey,forcez,forcex_re,forcey_re,forcez_re,forcex_hat,forcey_hat,forcez_hat, &
    u1_procs,u2_procs,u3_procs,u1_procs_re,u2_procs_re,u3_procs_re,u1_hat,u2_hat,u3_hat)
    !のりしろ境界
    call MPI_boundary(u1_procs,u2_procs,u3_procs,g_procs,taug_procs,forcex,forcey,forcez)
    !局所平衡分布関数の計算(OMP)
    call equilibrium_cal(p_procs,u1_procs,u2_procs,u3_procs,geq_procs)
    call MPI_boundary_g(geq_procs)
    !平衡分布関数fgの計算(OMP)
    call g_cal(gnext_procs,g_procs,geq_procs,forcex,forcey,forcez,taug_procs)
    !物理量の計算(OMP)
    call physics(p_procs,u1_procs,u2_procs,u3_procs,g_procs,gnext_procs,strain_procs,grad_u_procs,nue,taug_procs)
    !出力
    if((mod(n,step_output)==0)) then
        call output(u1_procs,u2_procs,u3_procs,p_procs,n)
        call kinetic_cal(kinetic_procs,u1_procs,u2_procs,u3_procs,n)
    endif

    if(mod(n,step_putput_fg)==0) then
        call outputg(g_procs,n)
    endif
    !時間計測出力
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    time2 = MPI_Wtime()
    if(comm_rank == 0) then
        if(mod(n,100)==0) then
            if(n == 100) then
                open(10,file="./time_xy.d")
                write(10,*) n, time2-time1, u1_procs(1,1,1)
                close(10)
            else
                open(10,file="./time_xy.d",action="write",position="append")
                write(10,*) n, time2-time1, u1_procs(1,1,1)
                close(10)
            endif
        endif
    endif
ENDDO
!!!!!!!!!!!!!!!!!!MPI並列とFFT終わり!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call decomp_2d_fft_finalize
    call decomp_2d_finalize
    call MPI_Finalize(ierr)
end program main