! 【高速計算】
! このプログラムは、MPI + OpenMP のハイブリット並列により高速計算を可能としている。
! 【系について】
! 二相のシミュレーションを行う。
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
    integer,parameter:: step = 1000000 !計算時間step
    integer,parameter:: step_input = 5000 !速度場入力時間step
    integer,parameter:: step_input_file_num = 400000 !入力する乱流場のstep
    !入力ディレクトリ
    character(*),parameter :: datadir_input = "/data/sht/nakanog/DNS_turbulence_256_IHT/fg/"
    !出力ディレクトリ
    character(*),parameter :: datadir_output = "/data/sht/nakanog/DNS_turbulence_256_IHT/case1/"
    character(*),parameter :: datadir_output_fg = "/data/sht/nakanog/DNS_turbulence_256_IHT/case1/fg/"
    integer,parameter:: step_output = 1000
    integer,parameter:: step_putput_fg = 100000

    !無次元数
    real(8),parameter:: We = 1.0d0 !ウェーバー数
    real(8),parameter:: eta = 1.0d0 !粘度比（nu2/nu1）

    !撹乱（乱数）のオーダー
    real(8), parameter :: ran = 1.0d-3

    !支配パラメータ
    real(8),parameter:: pi = acos(-1.0d0) !円周率
    real(8),parameter:: D = 40.0d0 !設置する液滴直径
    real(8),parameter:: nu1 = 0.001d0 !連続相の粘性係数
    real(8),parameter:: nu2 = eta*nu1 !分散相の粘性係数
    real(8),parameter:: sigma = 5.63d-4 !界面張力
    real(8),parameter:: kappaf = 0.01d0*ds**2 !界面厚さを決めるパラメータ
    ! real(8),parameter:: phi1 = 2.211d0 !連続相のオーダーパラメータ
    ! real(8),parameter:: phi2 = 4.895d0 !分散相のオーダーパラメータ
    ! real(8),parameter:: a = 9.0d0/49.0d0
    ! real(8),parameter:: b = 2.0d0/21.0d0
    ! real(8),parameter:: T = 0.55d0
    ! real(8),parameter:: kappag = (sigma/1.7039d0)**(1.0d0/0.9991d0)  !界面張力を決めるパラメータ
    real(8),parameter:: phi1 = 2.638d-1 !連続相のオーダーパラメータ
    real(8),parameter:: phi2 = 4.031d-1 !分散相のオーダーパラメータ
    real(8),parameter:: a = 1.0d0
    real(8),parameter:: b = 1.0d0
    real(8),parameter:: T = 2.93d-1
    real(8),parameter:: kappag = sigma / (6.22d-3)
    real(8),parameter:: tauf = 0.7d0 !緩和時間
    real(8),parameter:: Anu = 0.0d0 !粘性係数が小さいときに値を入れると良い

    !エネルギー注入率一定の外力
    integer,parameter:: kc = 3 !カットオフ波数
    real(8),parameter:: kolmogorov_scale = 1.0d0*ds !Kolmogorovスケール（想定値）
    real(8),parameter:: epsilon = (nu1**3.0d0) / (kolmogorov_scale**4.0d0) !エネルギー注入率

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

    subroutine ini(g_procs,gnext_procs,f_procs,fnext_procs,feq_procs,geq_procs,phi_procs,p_procs,u1_procs,u2_procs,u3_procs,grad_phi_procs,gphi_procs,lap_phi_procs,p0_procs,grad_u_procs)
        real(8),allocatable :: f_procs(:,:,:,:), fnext_procs(:,:,:,:), feq_procs(:,:,:,:)
        real(8),allocatable :: g_procs(:,:,:,:), gnext_procs(:,:,:,:), geq_procs(:,:,:,:)
        real(8),allocatable :: phi_procs(:,:,:), p_procs(:,:,:)
        real(8),allocatable :: u1_procs(:,:,:), u2_procs(:,:,:), u3_procs(:,:,:)
        real(8),allocatable :: p0_procs(:,:,:), lap_phi_procs(:,:,:), grad_phi_procs(:,:,:,:)
        real(8),allocatable :: gphi_procs(:,:,:,:,:), grad_u_procs(:,:,:,:,:)

        integer xi, yi, zi

    !========================並列数・コミュニケータ分割・通信先設定========================
        !配列のallocate(xi:0とx_procs+1がのりしろ)(yi:0とymax+2がのりしろ)(zi:0とz_procs+1がのりしろ)
        Nx = 32 !x方向の並列数（ただし，Nx/=comm_procs）
        Ny = comm_procs / (Nx * Nxall * Nyall) !z方向の並列数
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
        allocate(phi_procs(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(p_procs(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(f_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(fnext_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(feq_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(g_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(gnext_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(geq_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(grad_u_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1))
        allocate(p0_procs(1:x_procs,1:y_procs,1:zmax+1))
        allocate(lap_phi_procs(1:x_procs,1:y_procs,1:zmax+1))
        allocate(grad_phi_procs(1:3,1:x_procs,1:y_procs,1:zmax+1))
        allocate(gphi_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1))

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

!===================================初期条件の設定================================================
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    grobalx = (xi-1) + comm_rank_x * x_procs
                    grobaly = (yi-1) + comm_rank_y * y_procs
                    grobalz = (zi-1)
                    dif = (dble(grobalx)*ds-xc)**2 + (dble(grobaly)*ds-yc)**2 + (dble(grobalz)*ds-zc)**2
                    if(dif <= (0.5d0*D)**2) then
                        phi_procs(xi,yi,zi) = phi2
                    else
                        phi_procs(xi,yi,zi) = phi1
                    endif
                    u1_procs(xi,yi,zi) = 0.0d0
                    u2_procs(xi,yi,zi) = 0.0d0
                    u3_procs(xi,yi,zi) = 0.0d0
                    p_procs(xi,yi,zi) = 1.0d0 / 3.0d0
                enddo
            enddo
        enddo

        call glue(phi_procs)
        call glue(u1_procs)
        call glue(u2_procs)
        call glue(u3_procs)
        call various_cal(lap_phi_procs,grad_phi_procs,grad_u_procs,p0_procs,phi_procs,u1_procs,u2_procs,u3_procs)
        call gphi_cal(gphi_procs,grad_phi_procs)
        
        call equilibrium_cal(gphi_procs,phi_procs,p0_procs,lap_phi_procs,grad_phi_procs,grad_u_procs,p_procs,u1_procs,u2_procs,u3_procs,f_procs,g_procs)
    end subroutine ini

    subroutine ini_op(taug_procs,nu_procs,phi_procs,forcex,forcey,forcez)
        real(8),intent(in):: phi_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),allocatable :: taug_procs(:,:,:), nu_procs(:,:,:)
        real(8),allocatable :: forcex(:,:,:), forcey(:,:,:), forcez(:,:,:)
        integer xi, yi, zi

        !以下はのりしろ無しの変数
        allocate(forcex(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(forcey(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(forcez(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(taug_procs(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(nu_procs(1:x_procs,1:y_procs,1:zmax+1))

        taug_procs(:,:,:) = 0.0d0
        nu_procs(:,:,:) = 0.0d0
        forcex(:,:,:) = 0.0d0
        forcey(:,:,:) = 0.0d0
        forcez(:,:,:) = 0.0d0

        !$omp parallel
        !$omp do
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    nu_procs(xi,yi,zi) = (phi_procs(xi,yi,zi)-phi1)/(phi2-phi1)*(nu2-nu1) + nu1
                    taug_procs(xi,yi,zi) = 0.5d0 + 3.0d0 * nu_procs(xi,yi,zi)/ds + 2.0d0/3.0d0 * Anu
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

    subroutine MPI_boundary(phi_procs,u1_procs,u2_procs,u3_procs,f_procs,g_procs,taug_procs)
        real(8),intent(inout) :: phi_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout) :: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout) :: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout) :: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout) :: f_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout) :: g_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout) :: taug_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)

        call glue(phi_procs)
        call glue(u1_procs)
        call glue(u2_procs)
        call glue(u3_procs)
        call glue(taug_procs)
        call glue_f(f_procs)
        call glue_f(g_procs)
    endsubroutine MPI_boundary

    subroutine MPI_boundary_fg(feq_procs,geq_procs)
        real(8),intent(inout) :: feq_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout) :: geq_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)

        call glue_f(feq_procs)
        call glue_f(geq_procs)
    endsubroutine MPI_boundary_fg

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

    subroutine equilibrium_cal(gphi_procs,phi_procs,p0_procs,lap_phi_procs,grad_phi_procs,grad_u_procs,p_procs,u1_procs,u2_procs,u3_procs,feq_procs,geq_procs)
        real(8),intent(in):: gphi_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: phi_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: p0_procs(1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: lap_phi_procs(1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: grad_phi_procs(1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: grad_u_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: p_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(out):: feq_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(out):: geq_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        integer i, alpha, beta, xi, yi, zi

        !$omp parallel
        !$omp do
        do zi=1,zmax+1
            do yi=1,y_procs
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
        real(8),intent(inout):: fnext_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: f_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: feq_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: gnext_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: g_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: geq_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: taug_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: forcex(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: forcey(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: forcez(0:x_procs+1,0:y_procs+1,0:zmax+2)
        integer i, xi, yi, zi

        !$omp parallel
        !$omp do
        do zi=1,zmax+1
            do yi=1,y_procs
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

    subroutine physics(phi_procs,p_procs,u1_procs,u2_procs,u3_procs,f_procs,g_procs,nu_procs,taug_procs,fnext_procs,gnext_procs)
        real(8),intent(inout):: phi_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: p_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: f_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: g_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: nu_procs(1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(inout):: taug_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in) :: fnext_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in) :: gnext_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        integer i, xi, yi, zi

        !$omp parallel
        !$omp do
        do zi=1,zmax+1
            do yi=1,y_procs
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
        do zi=1,zmax+1
            do yi=1,y_procs
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

        !$omp parallel
        !$omp do
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    nu_procs(xi,yi,zi) = (phi_procs(xi,yi,zi)-phi_min)/(phi_max-phi_min)*(nu2-nu1) + nu1
                    taug_procs(xi,yi,zi) = 0.5d0 + 3.0d0*nu_procs(xi,yi,zi)/ds + 2.0d0/3.0d0*Anu
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine physics

    subroutine various_cal(lap_f,grad_f,grad_u_procs,p0_procs,fun,u1_procs,u2_procs,u3_procs)
        real(8),intent(out):: lap_f(1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(out):: grad_f(1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(out):: grad_u_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(out):: p0_procs(1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: fun(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        integer alpha, i, xi, yi, zi

        !$omp parallel
        !$omp do
        do zi=1,zmax+1
            do yi=1,y_procs
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
        real(8),intent(out):: gphi_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: grad_phi_procs(1:3,1:x_procs,1:y_procs,1:zmax+1)
        integer alpha, beta, xi, yi, zi

        !$omp parallel
        !$omp do
        do zi=1,zmax+1
            do yi=1,y_procs
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
        integer xi, yi, zi
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

        ! enesupe(:) = 0.0d0
        ! enesupe_sum(:) = 0.0d0
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

    subroutine output(phi_procs,u1_procs,u2_procs,u3_procs,p_procs,n)
        real(8),intent(inout):: phi_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
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
        filename=datadir_output//'0_'//trim(adjustl(filename3))//'_'//trim(adjustl(filename2))//'.bin' 
        print *, filename !表示してみる
        open(100,file=filename, form='unformatted',status='replace')
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    write(100) phi_procs(xi,yi,zi), u1_procs(xi,yi,zi), u2_procs(xi,yi,zi), u3_procs(xi,yi,zi), p_procs(xi,yi,zi)
                enddo
            enddo
        enddo
        close(100)
    end subroutine output

    subroutine outputfg(f_procs, g_procs, n)
        real(8),intent(inout):: f_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: g_procs(1:15,0:x_procs+1,0:y_procs+1,0:zmax+2)
        integer,intent(in):: n
        integer i, xi, yi, zi

        ! write(filename,*) group_new !i->filename 変換
        write(filename2,*) n
        write(filename3,*) comm_rank
        filename=datadir_output_fg//'0_'//trim(adjustl(filename3))//'_'//trim(adjustl(filename2))//'_fg.bin' 
        print *, filename !表示してみる
        open(102,file=filename, form='unformatted',status='replace') 
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    do i=1,15
                        write(102) f_procs(i,xi,yi,zi), g_procs(i,xi,yi,zi)
                    enddo
                enddo
            enddo
        enddo
        close(102)
    end subroutine outputfg

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
    real(8),allocatable :: f_procs(:,:,:,:), fnext_procs(:,:,:,:), feq_procs(:,:,:,:)
    real(8),allocatable :: g_procs(:,:,:,:), gnext_procs(:,:,:,:), geq_procs(:,:,:,:)
    real(8),allocatable :: phi_procs(:,:,:), p_procs(:,:,:)
    real(8),allocatable :: u1_procs(:,:,:), u2_procs(:,:,:), u3_procs(:,:,:)
    real(8),allocatable :: taug_procs(:,:,:), nu_procs(:,:,:)
    real(8),allocatable :: p0_procs(:,:,:), lap_phi_procs(:,:,:), grad_phi_procs(:,:,:,:), gphi_procs(:,:,:,:,:)
    real(8),allocatable :: grad_u_procs(:,:,:,:,:)
    real(8),allocatable :: forcex(:,:,:), forcey(:,:,:), forcez(:,:,:)

    !その他変数
    integer n, xi, yi, zi, i

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
    call mk_dirs(datadir_output)
    call mk_dirs(datadir_output_fg)
    !初期条件・配列のallocate
    call par(cx,cy,cz,cr,krone)
    call ini(g_procs,gnext_procs,f_procs,fnext_procs,feq_procs,geq_procs,phi_procs,p_procs,u1_procs,u2_procs,u3_procs,grad_phi_procs,gphi_procs,lap_phi_procs,p0_procs,grad_u_procs)
    call ini_op(taug_procs,nu_procs,phi_procs,forcex,forcey,forcez)
    call ini_fft(u1_hat,u2_hat,u3_hat,u1_procs_re,u2_procs_re,u3_procs_re,forcex_hat,forcey_hat,forcez_hat,forcex_re,forcey_re,forcez_re)
!!!!!!!!!!!!!!!!!!!!!!!時間発展開始!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    time1 = MPI_Wtime()
DO n=1,step
!============================流れの入力（撹乱を添加）=========================================
    if(n == step_input) then
        write(filename2,*) step_input_file_num
        write(filename3,*) comm_rank
        filename=datadir_input//'0_'//trim(adjustl(filename3))//'_'//trim(adjustl(filename2))//'_fg.bin' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
        print *, filename !表示してみる
        open(103, file=filename, form="unformatted")
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    do i=1,15
                        read(103) g_procs(i,xi,yi,zi)
                        ! read(103) dummy, g_procs(i,xi,yi,zi)
                    enddo
                enddo
            enddo
        enddo
        close(103)

        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    p_procs(xi,yi,zi) = 0.0d0
                    u1_procs(xi,yi,zi) = 0.0d0
                    u2_procs(xi,yi,zi) = 0.0d0
                    u3_procs(xi,yi,zi) = 0.0d0
                    do i=1,15
                        p_procs(xi,yi,zi) = p_procs(xi,yi,zi) + g_procs(i,xi,yi,zi) / 3.0d0
                        u1_procs(xi,yi,zi) = u1_procs(xi,yi,zi) + dble(cx(i))*g_procs(i,xi,yi,zi)
                        u2_procs(xi,yi,zi) = u2_procs(xi,yi,zi) + dble(cy(i))*g_procs(i,xi,yi,zi)
                        u3_procs(xi,yi,zi) = u3_procs(xi,yi,zi) + dble(cz(i))*g_procs(i,xi,yi,zi)
                    enddo
                enddo
            enddo
        enddo
    endif

    if(n >= step_input) then
        !外力
        call externalforce(forcex,forcey,forcez,forcex_re,forcey_re,forcez_re,forcex_hat,forcey_hat,forcez_hat, &
        u1_procs,u2_procs,u3_procs,u1_procs_re,u2_procs_re,u3_procs_re,u1_hat,u2_hat,u3_hat)
        call glue(forcex)
        call glue(forcey)
        call glue(forcez)
    endif
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
    call physics(phi_procs,p_procs,u1_procs,u2_procs,u3_procs,f_procs,g_procs,nu_procs,taug_procs,fnext_procs,gnext_procs)
!================================出力==================================
    if((mod(n,step_output)==0)) then
        call output(phi_procs,u1_procs,u2_procs,u3_procs,p_procs,n)
    endif

    if(mod(n,step_putput_fg)==0) then
        call outputfg(f_procs,g_procs,n)
    endif

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    time2 = MPI_Wtime()
    if(comm_rank == 0) then
        if(mod(n,100)==0) then
            if(n == 100) then
                open(10,file="./time.d")
                write(10,*) n, time2-time1, u1_procs(1,1,1)
                close(10)
            else
                open(10,file="./time.d",action="write",position="append")
                write(10,*) n, time2-time1, u1_procs(1,1,1)
                close(10)
            endif
        endif
    endif
ENDDO
!!!!!!!!!!!!!!!!!!MPI並列とFFT終わり!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call MPI_Finalize(ierr)
end program main
