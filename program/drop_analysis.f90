! 【Multi-times simulations】
! このプログラムは、多数回の液滴分裂の数値シミュレーションを実行する。
! このプログラムは、界面を「Cahn-Hilliard 方程式」で表現している。
! 液滴を長時間保持できないが、系全体の保存性には優れている。
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module globals
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
    integer,parameter:: step_start = 5000
    integer,parameter:: step_end = 22000
    integer,parameter:: step_bin = 1000
    integer,parameter:: step_num2 = (step_end - step_start) / step_bin + 1 

    !読み込みディレクトリ
    character(*),parameter :: datadir_input = "/data/sht/nakanog/DNS_turbulence_256_IHT/case4/"
    !出力ディレクトリ
    character(*),parameter :: datadir_output = "/data/sht/nakanog/DNS_turbulence_256_IHT/case4/large/"

    !粒子速度（整数）
    integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
    integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
    integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)

    !無次元数
    real(8),parameter:: We = 1.4d0 !粒子ウエーバー数
    real(8),parameter:: eta = 1.0d0 !粘度比（nu2/nu1）

    real(8),parameter:: epsilon = 1.32d-9 !エネルギー散逸率
    real(8),parameter:: kolmogorov = 0.93d0 !コルモゴロフ長
    real(8),parameter:: nu1 = 0.001d0 !連続相の粘性係数
    real(8),parameter:: nu2 = eta*nu1 !分散相の粘性係数

    real(8),parameter:: D = 40.0d0 !設置する液滴直径
    real(8),parameter:: phi1 = 2.638d-1 !連続相のオーダーパラメータ
    real(8),parameter:: phi2 = 4.031d-1 !分散相のオーダーパラメータ
    real(8),parameter:: a = 1.0d0
    real(8),parameter:: b = 1.0d0
    real(8),parameter:: T = 2.93d-1
    real(8),parameter:: sigma = epsilon**(2.0d0/3.0d0)*D**(5.0d0/3.0d0)/We !界面張力
    real(8),parameter:: kappag = sigma / (2.95d-3) !kappaf=0.06のときの近似式

    real(8),parameter:: pi = acos(-1.0d0) !円周率

    integer,parameter:: k_high = 3
    integer,parameter:: k_low = 1

    !その他変数
    real(8) dummy
    real(8) time1, time2 !計算時間測定用
    integer i, j, k, n, xi, yi, zi, alpha, beta
    character :: filename*200
    character :: filename2*200
    character :: filename3*200
    integer step_num

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
    subroutine par(cr)
        real(8), intent(out):: cr(1:3,1:15)
        do i=1,15
            cr(1,i) = dble(cx(i))
            cr(2,i) = dble(cy(i))
            cr(3,i) = dble(cz(i))
        enddo
    end subroutine par

    subroutine ini(p_procs,u1_procs,u2_procs,u3_procs,phi_procs)
        real(8),allocatable :: p_procs(:,:,:), u1_procs(:,:,:), u2_procs(:,:,:), u3_procs(:,:,:), phi_procs(:,:,:)

    !========================並列数・コミュニケータ分割・通信先設定========================
        !配列のallocate(xi:0とx_procs+1がのりしろ)(yi:0とymax+2がのりしろ)(zi:0とz_procs+1がのりしろ)
        Nx = 32 !x方向の並列数（ただし，Nx/=comm_procs）
        Ny = comm_procs / (Nx) !z方向の並列数
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
        allocate(phi_procs(0:x_procs+1,0:y_procs+1,0:zmax+2))

        !初期化
        p_procs(:,:,:) = 0.0d0
        u1_procs(:,:,:) = 0.0d0
        u2_procs(:,:,:) = 0.0d0
        u3_procs(:,:,:) = 0.0d0
        phi_procs(:,:,:) = 0.0d0
    end subroutine ini

    subroutine ini_fft(chemical_potencial_hat,advection_phi_hat)
        complex(kind(0d0)), allocatable :: chemical_potencial_hat(:,:,:), advection_phi_hat(:,:,:)
        
        call decomp_2d_init(xmax+1, ymax+1, zmax+1, Nx, Ny)
        call decomp_2d_fft_init(PHYSICAL_IN_Z)
        call decomp_2d_fft_get_size(sta, last, sized)

        allocate(chemical_potencial_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1))
        allocate(advection_phi_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1))

        chemical_potencial_hat(:,:,:) = (0.0d0, 0.0d0)
        advection_phi_hat(:,:,:) = (0.0d0, 0.0d0)
    end subroutine ini_fft

    subroutine variable(grad_phi_procs,chemical_potencial,free_energy,grad_free_energy,lap_phi_procs,advection_phi)
        real(8),allocatable :: grad_phi_procs(:,:,:,:)
        real(8),allocatable :: chemical_potencial(:,:,:)
        real(8),allocatable :: free_energy(:,:,:)
        real(8),allocatable :: grad_free_energy(:,:,:,:)
        real(8),allocatable :: lap_phi_procs(:,:,:)
        real(8),allocatable :: advection_phi(:,:,:)

        allocate(grad_phi_procs(1:3,1:x_procs,1:y_procs,1:zmax+1))
        allocate(chemical_potencial(1:x_procs,1:y_procs,1:zmax+1))
        allocate(free_energy(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(grad_free_energy(1:3,1:x_procs,1:y_procs,1:zmax+1))
        allocate(lap_phi_procs(1:x_procs,1:y_procs,1:zmax+1))
        allocate(advection_phi(1:x_procs,1:y_procs,1:zmax+1))

        grad_phi_procs(:,:,:,:) = 0.0d0
        chemical_potencial(:,:,:) = 0.0d0
        free_energy(:,:,:) = 0.0d0
        grad_free_energy(:,:,:,:) = 0.0d0
        lap_phi_procs(:,:,:) = 0.0d0
        advection_phi(:,:,:) = 0.0d0
    end subroutine variable

    subroutine input(u1_procs,u2_procs,u3_procs,p_procs,phi_procs)
        real(8),intent(inout):: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: p_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: phi_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)

        !==============================データ読み込み=================================================
        write(filename2,*) n
        write(filename3,*) comm_rank
        filename=datadir_input//'0_'//trim(adjustl(filename3))//'_'//trim(adjustl(filename2))//'.bin' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
        print *, filename !表示してみる
        open(10, file=filename, form="unformatted")
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    read(10) phi_procs(xi,yi,zi), u1_procs(xi,yi,zi), u2_procs(xi,yi,zi), u3_procs(xi,yi,zi), p_procs(xi,yi,zi)
                enddo
            enddo
        enddo
        close(10)
    end subroutine input

    subroutine glue(var)
        real(8),intent(inout) :: var(0:x_procs+1,0:y_procs+1,0:zmax+2)

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

    subroutine grad_cal(varout,var,cr)
        real(8),intent(out):: varout(1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: var(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: cr(1:3, 1:15)

        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    do alpha=1,3
                        varout(alpha,xi,yi,zi) = 0.0d0
                        do i = 2,15
                            varout(alpha,xi,yi,zi) = varout(alpha,xi,yi,zi) &
                                                    + cr(alpha,i)*var(xi+cx(i),yi+cy(i),zi+cz(i))
                        enddo
                        varout(alpha,xi,yi,zi) = varout(alpha,xi,yi,zi) / (10.0d0*ds)
                    enddo
                enddo
            enddo
        enddo
    end subroutine grad_cal

    subroutine lap_cal(lap_f,fun)
        real(8),intent(out):: lap_f(1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: fun(0:x_procs+1,0:y_procs+1,0:zmax+2)

        do zi=1,zmax+1
            do yi=1,y_procs
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

    subroutine free_energy_cal(free_energy,phi_procs)
        real(8),intent(out):: free_energy(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: phi_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)

        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    free_energy(xi,yi,zi) = phi_procs(xi,yi,zi)*T*log(phi_procs(xi,yi,zi)/(1.0d0-b*phi_procs(xi,yi,zi))) - a*phi_procs(xi,yi,zi)*phi_procs(xi,yi,zi)
                enddo
            enddo
        enddo
    end subroutine free_energy_cal

    subroutine chemical_potencial_cal(chemical_potencial,phi_procs,lap_phi_procs)
        real(8),intent(inout):: chemical_potencial(1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: phi_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: lap_phi_procs(1:x_procs,1:y_procs,1:zmax+1)

        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    chemical_potencial(xi,yi,zi) = ( T*log( phi_procs(xi,yi,zi)/(1.0d0-b*phi_procs(xi,yi,zi)) ) + T/(1.0d0-b*phi_procs(xi,yi,zi)) - 2.0d0*a*phi_procs(xi,yi,zi) ) - kappag*lap_phi_procs(xi,yi,zi)
                enddo
            enddo
        enddo
    end subroutine chemical_potencial_cal

    subroutine advection_cal(advection_phi,grad_phi_procs,u1_procs,u2_procs,u3_procs)
        real(8),intent(inout):: advection_phi(1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: grad_phi_procs(1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)

        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    advection_phi(xi,yi,zi) = u1_procs(xi,yi,zi)*grad_phi_procs(1,xi,yi,zi) &
                                            + u2_procs(xi,yi,zi)*grad_phi_procs(2,xi,yi,zi) &
                                            + u3_procs(xi,yi,zi)*grad_phi_procs(3,xi,yi,zi)
                enddo
            enddo
        enddo
    end subroutine advection_cal

    subroutine At_cal(grad_phi_procs,At_ini)
        real(8),intent(in):: grad_phi_procs(1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(inout):: At_ini
        real(8) At, At_sum

        At = 0.0d0
        At_sum = 0.0d0
        do zi = 1, zmax+1
            do yi = 1, y_procs
                do xi = 1, x_procs
                    At = At + sqrt( grad_phi_procs(1,xi,yi,zi)**2 + grad_phi_procs(2,xi,yi,zi)**2 + grad_phi_procs(3,xi,yi,zi)**2 )
                enddo
            enddo
        enddo

        call MPI_Reduce(At, At_sum, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        
        if(comm_rank == 0) then
            if(n == step_start) then
                At_ini = At_sum
                open(103,file="./At_case29.d")
                write(103,*) dble(n), At_sum / At_ini
                close(103)
            else
                open(103,file="./At_case29.d",action="write",position="append")
                write(103,*) dble(n), At_sum / At_ini
                close(103)
            endif
        endif
    end subroutine At_cal

    subroutine contribution(chemical_potencial,advection_phi,chemical_potencial_hat,advection_phi_hat,enesupe_phi,enesupe_phi_sum,enesupe_phi_result)
        real(8),intent(inout):: chemical_potencial(1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(inout):: advection_phi(1:x_procs,1:y_procs,1:zmax+1)

        complex(kind(0d0)), intent(inout) :: chemical_potencial_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)
        complex(kind(0d0)), intent(inout) :: advection_phi_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)

        real(8),intent(inout):: enesupe_phi((xmax+1)/2 + 1), enesupe_phi_sum((xmax+1)/2 + 1), enesupe_phi_result((xmax+1)/2 + 1)

        complex(kind(0d0)) comp_tmp
        integer k1, k2, k3, k(1:3), k_index
        real(8) k_abs
        integer kk((xmax+1)/2 + 1), kk_sum((xmax+1)/2 + 1)

        kk (:) = 0
        kk_sum(:) = 0

        !フーリエ変換（実数→複素数）
        call fft_r2c(chemical_potencial, chemical_potencial_hat)
        call fft_r2c(advection_phi, advection_phi_hat)

        enesupe_phi(:) = 0.0d0
        enesupe_phi_sum(:) = 0.0d0
        do k3=sta(3)-1,last(3)-1
            k(3) =( k3 - judge(k3, zmax+1) )
            do k2=sta(2)-1,last(2)-1
                k(2) = ( k2 - judge(k2, ymax+1) )
                do k1=sta(1)-1,last(1)-1
                    k(1) = ( k1 - judge(k1, xmax+1) )

                    k_abs = ( dble(k(1)**2 + k(2)**2 + k(3)**2) )**0.5d0
                    k_index = int( k_abs ) + 1

                    if(k3 == 0) then
                        comp_tmp = advection_phi_hat(k1,k2,k3) * conjg(chemical_potencial_hat(k1,k2,k3))
                        enesupe_phi(k_index) = enesupe_phi(k_index) + REAL(comp_tmp)

                        kk(k_index) = kk(k_index) + 1
                    else 
                        comp_tmp = advection_phi_hat(k1,k2,k3) * conjg(chemical_potencial_hat(k1,k2,k3))
                        enesupe_phi(k_index) = enesupe_phi(k_index) + 2.0d0 * REAL(comp_tmp) 

                        kk(k_index) = kk(k_index) + 2
                    endif
                enddo
            enddo
        enddo

        call MPI_Reduce(enesupe_phi(1), enesupe_phi_sum(1), (xmax+1)/2 + 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

        call MPI_Reduce(kk(1), kk_sum(1), (xmax+1)/2 + 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

        if(comm_rank == 0) then
            do i = 1, (xmax+1)/2 + 1
                enesupe_phi_result(i) = enesupe_phi_result(i) + enesupe_phi_sum(i) / dble(kk_sum(i))
                ! enesupe_phi_result(i) = enesupe_phi_result(i) + enesupe_phi_sum(i)
            enddo
        endif

        if(n == step_end) then
            if(comm_rank==0) then
                open(37,file ="./enesupe_phi_40we1.4_2.d")
                do i=1,(xmax+1)/2+1
                    enesupe_phi_result(i) = enesupe_phi_result(i) / dble(step_num2)
                    enesupe_phi_result(i) = -4.0d0*pi*(dble(i)-1.0d0)**2 * enesupe_phi_result(i)
                    write(37,"(4es16.8)") (dble(i)-1.0d0), enesupe_phi_result(i)
                enddo
                close(37)
            endif
        endif
    endsubroutine contribution

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
    real(8),allocatable :: p_procs(:,:,:), u1_procs(:,:,:), u2_procs(:,:,:), u3_procs(:,:,:), phi_procs(:,:,:)
    real(8),allocatable :: grad_phi_procs(:,:,:,:)
    real(8),allocatable :: chemical_potencial(:,:,:)
    real(8),allocatable :: free_energy(:,:,:), grad_free_energy(:,:,:,:)
    real(8),allocatable :: lap_phi_procs(:,:,:)
    real(8),allocatable :: advection_phi(:,:,:)
    real(8):: cr(1:3, 1:15)  !粒子速度（実数）
    real(8) At_ini

    !波数空間の変数
    complex(kind(0d0)), allocatable :: chemical_potencial_hat(:,:,:), advection_phi_hat(:,:,:)
    real(8),allocatable :: enesupe_phi(:), enesupe_phi_sum(:), enesupe_phi_result(:)

    At_ini = 0.0d0
    cr(:,:) = 0.0d0
    allocate(enesupe_phi((xmax+1)/2 + 1))
    allocate(enesupe_phi_sum((xmax+1)/2 + 1))
    allocate(enesupe_phi_result((xmax+1)/2 + 1))
    enesupe_phi(:) = 0.0d0
    enesupe_phi_sum(:) = 0.0d0
    enesupe_phi_result(:) = 0.0d0

!==================================MPI並列開始===============================================
    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, comm_procs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, comm_rank, ierr)
!================================ディレクトリの作成============================================
    ! call mk_dirs(datadir_output)
    call par(cr)
    call ini(p_procs,u1_procs,u2_procs,u3_procs,phi_procs)
    call ini_fft(chemical_potencial_hat,advection_phi_hat)
    call variable(grad_phi_procs,chemical_potencial,free_energy,grad_free_energy,lap_phi_procs,advection_phi)

DO n = step_start, step_end, step_bin
    call input(u1_procs,u2_procs,u3_procs,p_procs,phi_procs)
    call glue(phi_procs)
    call glue(u1_procs)
    call glue(u2_procs)
    call glue(u3_procs)
    !化学ポテンシャル
    call lap_cal(lap_phi_procs,phi_procs)
    call chemical_potencial_cal(chemical_potencial,phi_procs,lap_phi_procs)
    !移流項
    call grad_cal(grad_phi_procs,phi_procs,cr)
    call advection_cal(advection_phi,grad_phi_procs,u1_procs,u2_procs,u3_procs)
    !各スケールの界面への寄与
    call contribution(chemical_potencial,advection_phi,chemical_potencial_hat,advection_phi_hat,enesupe_phi,enesupe_phi_sum,enesupe_phi_result)
ENDDO

    ! if(comm_rank==0) then
    !     open(37,file ="./enesupe_phi.d")
    !     do i=1,(xmax+1)/2+1
    !         enesupe_phi_result(i) = enesupe_phi_result(i) / dble(step_num2)
    !         write(37,"(2es16.8)") dble(i)-1.0d0, enesupe_phi_result(i)
    !     enddo
    !     close(37)
    ! endif

!液滴の変形率の計算
! DO n = step_start, step_end, step_bin
!     call input(u1_procs,u2_procs,u3_procs,p_procs,phi_procs)
!     call glue(phi_procs)
!     call grad_cal(grad_phi_procs,phi_procs,cr)
!     call At_cal(grad_phi_procs,At_ini)
! ENDDO

!================MPI並列終わり=======================================
    call decomp_2d_fft_finalize
    call decomp_2d_finalize
    call MPI_Finalize(ierr)
end program main
