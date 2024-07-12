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
    integer,parameter:: step_end = 50000
    integer,parameter:: step_bin = 1000
    integer,parameter:: step_num2 = (step_end - step_start) / step_bin + 1 

    !読み込みディレクトリ
    character(*),parameter :: datadir_input = "/data/sht/nakanog/IHT_drop_d70_we5/"
    character(*),parameter :: datadir_input2 = "/data/sht/nakanog/IHT_drop_d70_we5/u/"
    !出力ディレクトリ
    character(*),parameter :: datadir_output = "/data/sht/nakanog/IHT_drop_d70_we5/contribution/"

    !粒子速度（整数）
    integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
    integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
    integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)

    !無次元数
    real(8),parameter:: We = 5.0d0 !粒子ウエーバー数
    real(8),parameter:: eta = 1.0d0 !粘度比（nu2/nu1）

    real(8),parameter:: epsilon = 9.15d-10 !エネルギー散逸率
    real(8),parameter:: kolmogorov = 1.0d0 !コルモゴロフ長
    real(8),parameter:: nu1 = 0.001d0 !連続相の粘性係数
    real(8),parameter:: nu2 = eta*nu1 !分散相の粘性係数

    real(8),parameter:: D = 70.0d0 !設置する液滴直径
    real(8),parameter:: phi1 = 2.638d-1 !連続相のオーダーパラメータ
    real(8),parameter:: phi2 = 4.031d-1 !分散相のオーダーパラメータ
    real(8),parameter:: a = 1.0d0
    real(8),parameter:: b = 1.0d0
    real(8),parameter:: T = 2.93d-1
    real(8),parameter:: sigma = epsilon**(2.0d0/3.0d0)*D**(5.0d0/3.0d0)/We !界面張力
    real(8),parameter:: kappag = sigma / (2.95d-3) !kappaf=0.06のときの近似式
    real(8),parameter:: kappaf = 0.06d0

    real(8),parameter:: pi = acos(-1.0d0) !円周率

    integer,parameter:: k_high = 127
    integer,parameter:: k_low = 9

    !その他変数
    real(8) dummy
    real(8) time1, time2 !計算時間測定用
    integer i, j, k, n, xi, yi, zi, alpha, beta
    integer i_input, Nxx, Nyy
    character :: filename*200
    character :: filename2*200
    character :: filename3*200
    character :: filenameu1*200
    character :: filenameu2*200
    character :: filenameu3*200
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

    subroutine ini_fft(chemical_potencial_hat,advection_phi_hat,grad_phi1_hat,grad_phi2_hat,grad_phi3_hat,grad_phi1_re,grad_phi2_re,grad_phi3_re,u1_hat,u2_hat,u3_hat,u1_procs_re,u2_procs_re,u3_procs_re,u1_procs_scale,u2_procs_scale,u3_procs_scale)
        complex(kind(0d0)), allocatable :: chemical_potencial_hat(:,:,:), advection_phi_hat(:,:,:)
        complex(kind(0d0)), allocatable :: grad_phi1_hat(:,:,:), grad_phi2_hat(:,:,:), grad_phi3_hat(:,:,:)
        real(8),allocatable :: grad_phi1_re(:,:,:), grad_phi2_re(:,:,:), grad_phi3_re(:,:,:)
        complex(kind(0d0)), allocatable :: u1_hat(:,:,:), u2_hat(:,:,:), u3_hat(:,:,:)
        real(8),allocatable :: u1_procs_re(:,:,:), u2_procs_re(:,:,:), u3_procs_re(:,:,:)
        real(8),allocatable :: u1_procs_scale(:,:,:), u2_procs_scale(:,:,:), u3_procs_scale(:,:,:)

        call decomp_2d_init(xmax+1, ymax+1, zmax+1, Nx, Ny)
        call decomp_2d_fft_init(PHYSICAL_IN_Z)
        call decomp_2d_fft_get_size(sta, last, sized)

        allocate(chemical_potencial_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1))
        allocate(advection_phi_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1))
        allocate(grad_phi1_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1))
        allocate(grad_phi2_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1))
        allocate(grad_phi3_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1))
        allocate(grad_phi1_re(0:x_procs-1, 0:y_procs-1, 0:zmax))
        allocate(grad_phi2_re(0:x_procs-1, 0:y_procs-1, 0:zmax))
        allocate(grad_phi3_re(0:x_procs-1, 0:y_procs-1, 0:zmax))

        chemical_potencial_hat(:,:,:) = (0.0d0, 0.0d0)
        advection_phi_hat(:,:,:) = (0.0d0, 0.0d0)
        grad_phi1_hat(:,:,:) = (0.0d0, 0.0d0)
        grad_phi2_hat(:,:,:) = (0.0d0, 0.0d0)
        grad_phi3_hat(:,:,:) = (0.0d0, 0.0d0)

        grad_phi1_re(:,:,:) = 0.0d0
        grad_phi2_re(:,:,:) = 0.0d0
        grad_phi3_re(:,:,:) = 0.0d0

        allocate(u1_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax))
        allocate(u2_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax))
        allocate(u3_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax))
        allocate(u1_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1))
        allocate(u2_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1))
        allocate(u3_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1))
        allocate(u1_procs_scale(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(u2_procs_scale(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(u3_procs_scale(0:x_procs+1,0:y_procs+1,0:zmax+2))

        u1_procs_re(:,:,:) = 0.0d0
        u2_procs_re(:,:,:) = 0.0d0
        u3_procs_re(:,:,:) = 0.0d0
        u1_hat(:,:,:) = (0.0d0, 0.0d0)
        u2_hat(:,:,:) = (0.0d0, 0.0d0)
        u3_hat(:,:,:) = (0.0d0, 0.0d0)
        u1_procs_scale(:,:,:) = 0.0d0
        u2_procs_scale(:,:,:) = 0.0d0
        u3_procs_scale(:,:,:) = 0.0d0
    end subroutine ini_fft

    subroutine variable(grad_phi_procs,chemical_potencial,free_energy,grad_free_energy,lap_phi_procs,advection_phi,grad_u_procs,strain_procs,korteweg_procs,interaction_procs)
        real(8),allocatable :: grad_phi_procs(:,:,:,:)
        real(8),allocatable :: chemical_potencial(:,:,:)
        real(8),allocatable :: free_energy(:,:,:)
        real(8),allocatable :: grad_free_energy(:,:,:,:)
        real(8),allocatable :: lap_phi_procs(:,:,:)
        real(8),allocatable :: advection_phi(:,:,:)
        real(8),allocatable :: grad_u_procs(:,:,:,:,:)
        real(8),allocatable :: strain_procs(:,:,:,:,:)
        real(8),allocatable :: korteweg_procs(:,:,:,:,:)
        real(8),allocatable :: interaction_procs(:,:,:)

        allocate(grad_u_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1))
        allocate(strain_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1))
        allocate(korteweg_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1))
        allocate(interaction_procs(1:x_procs,1:y_procs,1:zmax+1))

        allocate(grad_phi_procs(1:3,1:x_procs,1:y_procs,1:zmax+1))
        allocate(chemical_potencial(1:x_procs,1:y_procs,1:zmax+1))
        allocate(free_energy(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(grad_free_energy(1:3,1:x_procs,1:y_procs,1:zmax+1))
        allocate(lap_phi_procs(1:x_procs,1:y_procs,1:zmax+1))
        allocate(advection_phi(1:x_procs,1:y_procs,1:zmax+1))

        grad_u_procs(:,:,:,:,:) = 0.0d0
        strain_procs(:,:,:,:,:) = 0.0d0
        grad_phi_procs(:,:,:,:) = 0.0d0
        chemical_potencial(:,:,:) = 0.0d0
        free_energy(:,:,:) = 0.0d0
        grad_free_energy(:,:,:,:) = 0.0d0
        lap_phi_procs(:,:,:) = 0.0d0
        advection_phi(:,:,:) = 0.0d0
        korteweg_procs(:,:,:,:,:) = 0.0d0
    end subroutine variable

    subroutine ini_energy3(tensor11,tensor12,tensor13,tensor21,tensor22,tensor23,tensor31,tensor32,tensor33)
        real(8),allocatable :: tensor11(:,:,:), tensor12(:,:,:), tensor13(:,:,:)
        real(8),allocatable :: tensor21(:,:,:), tensor22(:,:,:), tensor23(:,:,:)
        real(8),allocatable :: tensor31(:,:,:), tensor32(:,:,:), tensor33(:,:,:)

        !以下はのりしろ無しの変数
        allocate(tensor11(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(tensor12(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(tensor13(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(tensor21(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(tensor22(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(tensor23(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(tensor31(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(tensor32(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(tensor33(0:x_procs+1,0:y_procs+1,0:zmax+2))

        tensor11(:,:,:) = 0.0d0
        tensor12(:,:,:) = 0.0d0
        tensor13(:,:,:) = 0.0d0
        tensor21(:,:,:) = 0.0d0
        tensor22(:,:,:) = 0.0d0
        tensor23(:,:,:) = 0.0d0
        tensor31(:,:,:) = 0.0d0
        tensor32(:,:,:) = 0.0d0
        tensor33(:,:,:) = 0.0d0
    end subroutine ini_energy3

    subroutine ini_energy5(grad_tensor11,grad_tensor12,grad_tensor13,grad_tensor21,grad_tensor22,grad_tensor23,grad_tensor31,grad_tensor32,grad_tensor33)
        real(8),allocatable :: grad_tensor11(:,:,:), grad_tensor12(:,:,:), grad_tensor13(:,:,:)
        real(8),allocatable :: grad_tensor21(:,:,:), grad_tensor22(:,:,:), grad_tensor23(:,:,:)
        real(8),allocatable :: grad_tensor31(:,:,:), grad_tensor32(:,:,:), grad_tensor33(:,:,:)

        !以下はのりしろ無しの変数
        allocate(grad_tensor11(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(grad_tensor12(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(grad_tensor13(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(grad_tensor21(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(grad_tensor22(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(grad_tensor23(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(grad_tensor31(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(grad_tensor32(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(grad_tensor33(0:x_procs+1,0:y_procs+1,0:zmax+2))

        grad_tensor11(:,:,:) = 0.0d0
        grad_tensor12(:,:,:) = 0.0d0
        grad_tensor13(:,:,:) = 0.0d0
        grad_tensor21(:,:,:) = 0.0d0
        grad_tensor22(:,:,:) = 0.0d0
        grad_tensor23(:,:,:) = 0.0d0
        grad_tensor31(:,:,:) = 0.0d0
        grad_tensor32(:,:,:) = 0.0d0
        grad_tensor33(:,:,:) = 0.0d0
    end subroutine ini_energy5

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

        real(8) energy_large, energy_medium, energy_small, energy_all

        kk (:) = 0
        kk_sum(:) = 0

        !フーリエ変換（実数→複素数）
        call fft_r2c(chemical_potencial, chemical_potencial_hat)
        call fft_r2c(advection_phi, advection_phi_hat)

        enesupe_phi(:) = 0.0d0
        enesupe_phi_sum(:) = 0.0d0
        enesupe_phi_result(:) = 0.0d0
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

        energy_all = 0.0d0
        energy_large = 0.0d0
        energy_medium = 0.0d0
        energy_small = 0.0d0
        if(comm_rank == 0) then
            do i = 2, 4
                energy_large = energy_large + enesupe_phi_result(i)
                energy_all = energy_all + enesupe_phi_result(i)
            enddo

            do i = 5, 10
                energy_medium = energy_medium + enesupe_phi_result(i)
                energy_all = energy_all + enesupe_phi_result(i)
            enddo

            do i = 11, (xmax+1)/2 + 1
                energy_small = energy_small + enesupe_phi_result(i)
                energy_all = energy_all + enesupe_phi_result(i)
            enddo

            energy_all = energy_all * (-1.0d0)
            energy_large = energy_large * (-1.0d0)
            energy_medium = energy_medium * (-1.0d0)
            energy_small = energy_small * (-1.0d0)

            if(n == step_start) then
                open(108,file="./energy_time_d70we1.4.d")
                write(108,*) dble(n), energy_all, energy_large, energy_medium, energy_small
                close(108)
            else
                open(108,file="./energy_time_d70we1.4.d",action="write",position="append")
                write(108,*) dble(n), energy_all, energy_large, energy_medium, energy_small
                close(108)
            endif
        endif

        ! if(n == step_end) then
        !     if(comm_rank==0) then
        !         open(37,file ="./enesupe_phi_40we1.4_2.d")
        !         do i=1,(xmax+1)/2+1
        !             enesupe_phi_result(i) = enesupe_phi_result(i) / dble(step_num2)
        !             enesupe_phi_result(i) = -4.0d0*pi*(dble(i)-1.0d0)**2 * enesupe_phi_result(i)
        !             write(37,"(4es16.8)") (dble(i)-1.0d0), enesupe_phi_result(i)
        !         enddo
        !         close(37)
        !     endif
        ! endif
    endsubroutine contribution

    subroutine phi_supectrum(grad_phi1_hat,grad_phi2_hat,grad_phi3_hat,grad_phi_procs,grad_phi1_re,grad_phi2_re,grad_phi3_re,phisupe,phisupe_sum,free_energy,phi_procs)
        complex(kind(0d0)), intent(inout) :: grad_phi1_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)
        complex(kind(0d0)), intent(inout) :: grad_phi2_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)
        complex(kind(0d0)), intent(inout) :: grad_phi3_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)
        real(8),intent(inout):: phi_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: grad_phi_procs(1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(inout):: grad_phi1_re(0:x_procs-1, 0:y_procs-1, 0:zmax) 
        real(8),intent(inout):: grad_phi2_re(0:x_procs-1, 0:y_procs-1, 0:zmax) 
        real(8),intent(inout):: grad_phi3_re(0:x_procs-1, 0:y_procs-1, 0:zmax) 

        real(8),intent(inout):: free_energy(0:x_procs+1,0:y_procs+1,0:zmax+2)

        real(8),intent(inout):: phisupe((xmax+1)/2 + 1), phisupe_sum((xmax+1)/2 + 1)
        integer k1, k2, k3, k(1:3), k_index
        real(8) k_abs

        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    grad_phi1_re(xi-1,yi-1,zi-1) = phi_procs(xi,yi,zi)
                    ! grad_phi2_re(xi-1,yi-1,zi-1) = grad_phi_procs(2,xi,yi,zi)*grad_phi_procs(2,xi,yi,zi)
                    ! grad_phi3_re(xi-1,yi-1,zi-1) = grad_phi_procs(3,xi,yi,zi)*grad_phi_procs(3,xi,yi,zi)
                enddo
            enddo
        enddo

        !フーリエ変換（実数→複素数）
        call fft_r2c(grad_phi1_re, grad_phi1_hat)
        ! call fft_r2c(grad_phi2_re, grad_phi2_hat)
        ! call fft_r2c(grad_phi3_re, grad_phi3_hat)

        phisupe(:) = 0.0d0
        phisupe_sum(:) = 0.0d0
        do k3=sta(3)-1,last(3)-1
            k(3) =( k3 - judge(k3, zmax+1) )
            do k2=sta(2)-1,last(2)-1
                k(2) = ( k2 - judge(k2, ymax+1) )
                do k1=sta(1)-1,last(1)-1
                    k(1) = ( k1 - judge(k1, xmax+1) )

                    k_abs = ( dble(k(1)**2 + k(2)**2 + k(3)**2) )**0.5d0
                    k_index = int( k_abs ) + 1

                    if(k3 == 0) then
                        ! phisupe(k_index) = phisupe(k_index) + ( abs(grad_phi1_hat(k1,k2,k3))+abs(grad_phi2_hat(k1,k2,k3))+abs(grad_phi3_hat(k1,k2,k3)) )
                        phisupe(k_index) = phisupe(k_index) + abs(grad_phi1_hat(k1,k2,k3))
                    else 
                        ! phisupe(k_index) = phisupe(k_index) + ( abs(grad_phi1_hat(k1,k2,k3))+abs(grad_phi2_hat(k1,k2,k3))+abs(grad_phi3_hat(k1,k2,k3)) ) * 2.0d0
                        phisupe(k_index) = phisupe(k_index) + abs(grad_phi1_hat(k1,k2,k3)) * 2.0d0
                    endif
                enddo
            enddo
        enddo

        call MPI_Reduce(phisupe(1), phisupe_sum(1), (xmax+1)/2 + 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    end subroutine phi_supectrum

    subroutine scale_cal(u1_procs,u2_procs,u3_procs,u1_procs_re,u2_procs_re,u3_procs_re,u1_hat,u2_hat,u3_hat)
        real(8),intent(inout):: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u1_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax)
        real(8),intent(inout):: u2_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax)
        real(8),intent(inout):: u3_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax)
        complex(kind(0d0)), intent(inout) :: u1_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)
        complex(kind(0d0)), intent(inout) :: u2_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)
        complex(kind(0d0)), intent(inout) :: u3_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)
        integer xi, yi, zi
        integer k1, k2, k3, k(1:3), k_index
        real(8) k_abs

        !u_hatの配列数と合わせる
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    u1_procs_re(xi-1,yi-1,zi-1) = u1_procs(xi,yi,zi)
                    u2_procs_re(xi-1,yi-1,zi-1) = u2_procs(xi,yi,zi)
                    u3_procs_re(xi-1,yi-1,zi-1) = u3_procs(xi,yi,zi)
                enddo
            enddo
        enddo

        !フーリエ変換（実数→複素数）
        call fft_r2c(u1_procs_re, u1_hat)
        call fft_r2c(u2_procs_re, u2_hat)
        call fft_r2c(u3_procs_re, u3_hat)

        !バンドパスフィルターをかける
        do k3=sta(3)-1,last(3)-1
            k(3) = k3 - judge(k3, zmax+1)
            do k2=sta(2)-1,last(2)-1
                k(2) = k2 - judge(k2, ymax+1)
                do k1=sta(1)-1,last(1)-1
                    k(1) = k1 - judge(k1, xmax+1)

                    k_abs = ( dble(k(1))**2.0d0 + dble(k(2))**2.0d0 + dble(k(3))**2.0d0 )**0.5d0
                    k_index = int( k_abs ) + 1

                    if(k_index < k_low+1) then
                        u1_hat(k1,k2,k3) = (0.0d0, 0.0d0)
                        u2_hat(k1,k2,k3) = (0.0d0, 0.0d0)
                        u3_hat(k1,k2,k3) = (0.0d0, 0.0d0)
                    endif

                    if(k_index >= k_high+1) then
                        u1_hat(k1,k2,k3) = (0.0d0, 0.0d0)
                        u2_hat(k1,k2,k3) = (0.0d0, 0.0d0)
                        u3_hat(k1,k2,k3) = (0.0d0, 0.0d0)
                    endif
                enddo
            enddo
        enddo

        u1_procs_re(:,:,:) = 0.0d0
        u2_procs_re(:,:,:) = 0.0d0
        u3_procs_re(:,:,:) = 0.0d0
        u1_procs(:,:,:) = 0.0d0
        u2_procs(:,:,:) = 0.0d0
        u3_procs(:,:,:) = 0.0d0
        !逆フーリエ変換（複素数→実数）
        call fft_c2r(u1_hat, u1_procs_re)
        call fft_c2r(u2_hat, u2_procs_re)
        call fft_c2r(u3_hat, u3_procs_re)

        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    u1_procs(xi,yi,zi) = u1_procs_re(xi-1,yi-1,zi-1)
                    u2_procs(xi,yi,zi) = u2_procs_re(xi-1,yi-1,zi-1)
                    u3_procs(xi,yi,zi) = u3_procs_re(xi-1,yi-1,zi-1)
                enddo
            enddo
        enddo
    end subroutine scale_cal

    subroutine scale_cal2(u1_procs,u2_procs,u3_procs,u1_procs_re,u2_procs_re,u3_procs_re,u1_hat,u2_hat,u3_hat,k_analysis,u1_procs_scale,u2_procs_scale,u3_procs_scale)
        real(8),intent(inout):: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u1_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax)
        real(8),intent(inout):: u2_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax)
        real(8),intent(inout):: u3_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax)
        complex(kind(0d0)), intent(inout) :: u1_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)
        complex(kind(0d0)), intent(inout) :: u2_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)
        complex(kind(0d0)), intent(inout) :: u3_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)
        integer,intent(in) :: k_analysis
        real(8),intent(inout):: u1_procs_scale(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u2_procs_scale(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u3_procs_scale(0:x_procs+1,0:y_procs+1,0:zmax+2)
        integer xi, yi, zi
        integer k1, k2, k3, k(1:3), k_index
        real(8) k_abs

        !u_hatの配列数と合わせる
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    u1_procs_re(xi-1,yi-1,zi-1) = u1_procs(xi,yi,zi)
                    u2_procs_re(xi-1,yi-1,zi-1) = u2_procs(xi,yi,zi)
                    u3_procs_re(xi-1,yi-1,zi-1) = u3_procs(xi,yi,zi)
                enddo
            enddo
        enddo

        !フーリエ変換（実数→複素数）
        call fft_r2c(u1_procs_re, u1_hat)
        call fft_r2c(u2_procs_re, u2_hat)
        call fft_r2c(u3_procs_re, u3_hat)

        !バンドパスフィルターをかける
        do k3=sta(3)-1,last(3)-1
            k(3) = k3 - judge(k3, zmax+1)
            do k2=sta(2)-1,last(2)-1
                k(2) = k2 - judge(k2, ymax+1)
                do k1=sta(1)-1,last(1)-1
                    k(1) = k1 - judge(k1, xmax+1)

                    k_abs = ( dble(k(1))**2.0d0 + dble(k(2))**2.0d0 + dble(k(3))**2.0d0 )**0.5d0
                    k_index = int( k_abs ) + 1

                    if(k_index - 1 == k_analysis) then
                        u1_hat(k1,k2,k3) = u1_hat(k1,k2,k3)
                        u2_hat(k1,k2,k3) = u2_hat(k1,k2,k3)
                        u3_hat(k1,k2,k3) = u3_hat(k1,k2,k3)
                    else
                        u1_hat(k1,k2,k3) = (0.0d0, 0.0d0)
                        u2_hat(k1,k2,k3) = (0.0d0, 0.0d0)
                        u3_hat(k1,k2,k3) = (0.0d0, 0.0d0)
                    endif
                enddo
            enddo
        enddo

        u1_procs_re(:,:,:) = 0.0d0
        u2_procs_re(:,:,:) = 0.0d0
        u3_procs_re(:,:,:) = 0.0d0
        !逆フーリエ変換（複素数→実数）
        call fft_c2r(u1_hat, u1_procs_re)
        call fft_c2r(u2_hat, u2_procs_re)
        call fft_c2r(u3_hat, u3_procs_re)

        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    u1_procs_scale(xi,yi,zi) = u1_procs_re(xi-1,yi-1,zi-1)
                    u2_procs_scale(xi,yi,zi) = u2_procs_re(xi-1,yi-1,zi-1)
                    u3_procs_scale(xi,yi,zi) = u3_procs_re(xi-1,yi-1,zi-1)
                enddo
            enddo
        enddo
    end subroutine scale_cal2

    subroutine grad_u_cal(grad_u_procs,u1_procs,u2_procs,u3_procs,cr)
        real(8),intent(inout):: grad_u_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: cr(1:3, 1:15)

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
    end subroutine grad_u_cal

    subroutine strain_cal(grad_u_procs,strain_procs)
        real(8),intent(in):: grad_u_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(out):: strain_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)
        integer xi, yi, zi, alpha, beta 

        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    do beta=1,3
                        do alpha=1,3
                            strain_procs(alpha,beta,xi,yi,zi) = 0.5d0 * (grad_u_procs(alpha,beta,xi,yi,zi) + grad_u_procs(beta,alpha,xi,yi,zi))
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine strain_cal

    subroutine korteweg_cal(grad_phi_procs,korteweg_procs)
        real(8),intent(in):: grad_phi_procs(1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(inout):: korteweg_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)

        do zi = 1, zmax + 1
            do yi = 1, y_procs
                do xi = 1, x_procs
                    do beta = 1, 3
                        do alpha = 1, 3
                            korteweg_procs(alpha,beta,xi,yi,zi) = grad_phi_procs(alpha,xi,yi,zi)*grad_phi_procs(beta,xi,yi,zi)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine korteweg_cal

    subroutine interaction_cal(interaction_procs,strain_procs,korteweg_procs)
        real(8),intent(inout):: interaction_procs(1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: strain_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: korteweg_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)

        do zi = 1, zmax + 1
            do yi = 1, y_procs
                do xi = 1, x_procs
                    interaction_procs(xi,yi,zi) = 0.0d0
                    do beta = 1, 3
                        do alpha = 1, 3
                            interaction_procs(xi,yi,zi) = interaction_procs(xi,yi,zi) + korteweg_procs(alpha,beta,xi,yi,zi)*strain_procs(alpha,beta,xi,yi,zi)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine interaction_cal

    subroutine output(interaction_procs)
        real(8),intent(inout):: interaction_procs(1:x_procs,1:y_procs,1:zmax+1)

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
                    write(100) interaction_procs(xi,yi,zi)
                enddo
            enddo
        enddo
        close(100)
    end subroutine output

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
    real(8),allocatable :: u1_procs_scale(:,:,:), u2_procs_scale(:,:,:), u3_procs_scale(:,:,:)
    real(8),allocatable :: grad_phi_procs(:,:,:,:)
    real(8),allocatable :: chemical_potencial(:,:,:)
    real(8),allocatable :: free_energy(:,:,:), grad_free_energy(:,:,:,:)
    real(8),allocatable :: lap_phi_procs(:,:,:)
    real(8),allocatable :: advection_phi(:,:,:)
    real(8):: cr(1:3, 1:15)  !粒子速度（実数）
    real(8) At_ini
    real(8),allocatable :: grad_u_procs(:,:,:,:,:)
    real(8),allocatable :: strain_procs(:,:,:,:,:)
    real(8),allocatable :: korteweg_procs(:,:,:,:,:)
    real(8),allocatable :: interaction_procs(:,:,:)

    real(8),allocatable :: tensor11(:,:,:), tensor12(:,:,:), tensor13(:,:,:)
    real(8),allocatable :: tensor21(:,:,:), tensor22(:,:,:), tensor23(:,:,:)
    real(8),allocatable :: tensor31(:,:,:), tensor32(:,:,:), tensor33(:,:,:)

    real(8),allocatable :: grad_tensor11(:,:,:), grad_tensor12(:,:,:), grad_tensor13(:,:,:)
    real(8),allocatable :: grad_tensor21(:,:,:), grad_tensor22(:,:,:), grad_tensor23(:,:,:)
    real(8),allocatable :: grad_tensor31(:,:,:), grad_tensor32(:,:,:), grad_tensor33(:,:,:)

    real(8) interaction_sum, interaction_all
    real(8) free_energy_sum, free_energy_all
    real(8) kinetic_energy_sum, kinetic_energy_all

    integer k_analysis
    real(8),allocatable :: contribution_each_scale(:)

    !波数空間の変数
    complex(kind(0d0)), allocatable :: chemical_potencial_hat(:,:,:), advection_phi_hat(:,:,:)
    real(8),allocatable :: enesupe_phi(:), enesupe_phi_sum(:), enesupe_phi_result(:)

    complex(kind(0d0)), allocatable :: grad_phi1_hat(:,:,:), grad_phi2_hat(:,:,:), grad_phi3_hat(:,:,:)
    real(8),allocatable :: grad_phi1_re(:,:,:), grad_phi2_re(:,:,:), grad_phi3_re(:,:,:)
    real(8),allocatable :: phisupe(:), phisupe_sum(:), phisupe_result(:)

    complex(kind(0d0)), allocatable :: u1_hat(:,:,:), u2_hat(:,:,:), u3_hat(:,:,:)
    real(8),allocatable :: u1_procs_re(:,:,:), u2_procs_re(:,:,:), u3_procs_re(:,:,:)

    real(8) phi(0:xmax,0:ymax,0:zmax), u1(0:xmax,0:ymax,0:zmax), u2(0:xmax,0:ymax,0:zmax), u3(0:xmax,0:ymax,0:zmax)
    real(8),allocatable :: tmp1(:,:,:,:), tmp2(:,:,:,:), tmp3(:,:,:,:), tmp4(:,:,:,:)

    At_ini = 0.0d0
    cr(:,:) = 0.0d0
    allocate(enesupe_phi((xmax+1)/2 + 1))
    allocate(enesupe_phi_sum((xmax+1)/2 + 1))
    allocate(enesupe_phi_result((xmax+1)/2 + 1))
    enesupe_phi(:) = 0.0d0
    enesupe_phi_sum(:) = 0.0d0
    enesupe_phi_result(:) = 0.0d0

    allocate(phisupe((xmax+1)/2 + 1))
    allocate(phisupe_sum((xmax+1)/2 + 1))
    allocate(phisupe_result((xmax+1)/2 + 1))
    phisupe(:) = 0.0d0
    phisupe_sum(:) = 0.0d0
    phisupe_result(:) = 0.0d0

    allocate(contribution_each_scale((xmax+1)/2 + 1))
    contribution_each_scale(:) = 0.0d0

!==================================MPI並列開始===============================================
    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, comm_procs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, comm_rank, ierr)
!================================ディレクトリの作成============================================
    call mk_dirs(datadir_output)
    call par(cr)
    call ini(p_procs,u1_procs,u2_procs,u3_procs,phi_procs)
    call ini_fft(chemical_potencial_hat,advection_phi_hat,grad_phi1_hat,grad_phi2_hat,grad_phi3_hat,grad_phi1_re,grad_phi2_re,grad_phi3_re,u1_hat,u2_hat,u3_hat,u1_procs_re,u2_procs_re,u3_procs_re,u1_procs_scale,u2_procs_scale,u3_procs_scale)
    call variable(grad_phi_procs,chemical_potencial,free_energy,grad_free_energy,lap_phi_procs,advection_phi,grad_u_procs,strain_procs,korteweg_procs,interaction_procs)
    call ini_energy3(tensor11,tensor12,tensor13,tensor21,tensor22,tensor23,tensor31,tensor32,tensor33)
    call ini_energy5(grad_tensor11,grad_tensor12,grad_tensor13,grad_tensor21,grad_tensor22,grad_tensor23,grad_tensor31,grad_tensor32,grad_tensor33)

    allocate(tmp1(0:comm_procs-1,0:x_procs+1,0:y_procs+1,0:zmax+2))
    allocate(tmp2(0:comm_procs-1,0:x_procs+1,0:y_procs+1,0:zmax+2))
    allocate(tmp3(0:comm_procs-1,0:x_procs+1,0:y_procs+1,0:zmax+2))
    allocate(tmp4(0:comm_procs-1,0:x_procs+1,0:y_procs+1,0:zmax+2))

DO n = step_start, step_end, step_bin
    ! call input(u1_procs,u2_procs,u3_procs,p_procs,phi_procs)
    write(filename2,*) n
    filename=datadir_input//'1_'//trim(adjustl(filename2))//'.bin'
    filenameu1=datadir_input2//'1_'//trim(adjustl(filename2))//'_u1.bin'
    filenameu2=datadir_input2//'1_'//trim(adjustl(filename2))//'_u2.bin'
    filenameu3=datadir_input2//'1_'//trim(adjustl(filename2))//'_u3.bin'
    open(10, file=filename, form="unformatted")
    open(11, file=filenameu1, form="unformatted")
    open(12, file=filenameu2, form="unformatted")
    open(13, file=filenameu3, form="unformatted")
    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                read(10) phi(xi,yi,zi)
                read(11) u1(xi,yi,zi)
                read(12) u2(xi,yi,zi)
                read(13) u3(xi,yi,zi)
            enddo
        enddo
    enddo
    close(10)
    close(11)
    close(12)
    close(13)

    i_input = 0
    do Nxx=0,Nx-1
        do Nyy=0,Ny-1
            do zi=1,zmax+1
                do yi=1,y_procs
                    do xi=1,x_procs
                        tmp1(i_input,xi,yi,zi) = u1((xi-1)+Nxx*x_procs,(yi-1)+Nyy*y_procs,zi-1)
                        tmp2(i_input,xi,yi,zi) = u2((xi-1)+Nxx*x_procs,(yi-1)+Nyy*y_procs,zi-1)
                        tmp3(i_input,xi,yi,zi) = u3((xi-1)+Nxx*x_procs,(yi-1)+Nyy*y_procs,zi-1)
                        tmp4(i_input,xi,yi,zi) = phi((xi-1)+Nxx*x_procs,(yi-1)+Nyy*y_procs,zi-1)
                    enddo
                enddo
            enddo
            i_input = i_input + 1
        enddo
    enddo
    do i_input=0,comm_procs-1
        if(i_input == comm_rank)then
            do zi=1,zmax+1
                do yi=1,y_procs
                    do xi=1,x_procs
                        u1_procs(xi,yi,zi) = tmp1(i_input,xi,yi,zi)
                        u2_procs(xi,yi,zi) = tmp2(i_input,xi,yi,zi)
                        u3_procs(xi,yi,zi) = tmp3(i_input,xi,yi,zi)
                        phi_procs(xi,yi,zi) = tmp4(i_input,xi,yi,zi)
                    enddo
                enddo
            enddo
        endif
    enddo
    call glue(phi_procs)
    call glue(u1_procs)
    call glue(u2_procs)
    call glue(u3_procs)
    !化学ポテンシャル
    ! call lap_cal(lap_phi_procs,phi_procs)
    ! call chemical_potencial_cal(chemical_potencial,phi_procs,lap_phi_procs)
    !移流項
    ! call grad_cal(grad_phi_procs,phi_procs,cr)
    ! call advection_cal(advection_phi,grad_phi_procs,u1_procs,u2_procs,u3_procs)
    !各スケールの界面への寄与
    ! call contribution(chemical_potencial,advection_phi,chemical_potencial_hat,advection_phi_hat,enesupe_phi,enesupe_phi_sum,enesupe_phi_result)
    !自由エネルギースペクトル
    ! if((n == step_start)) then
        ! call free_energy_cal(free_energy,phi_procs)
    !     call phi_supectrum(grad_phi1_hat,grad_phi2_hat,grad_phi3_hat,grad_phi_procs,grad_phi1_re,grad_phi2_re,grad_phi3_re,phisupe,phisupe_sum,free_energy,phi_procs)
    !     if(comm_rank == 0) then
    !         do i = 1, (xmax+1)/2 + 1
    !             phisupe_result(i) = phisupe_result(i) + phisupe_sum(i)
    !         enddo

    !         open(37,file ="./phisupe_d70we1.4_phionly.d")
    !         do i=1,(xmax+1)/2+1
    !             write(37,"(2es16.8)") dble(i)-1.0d0, phisupe_result(i)
    !         enddo
    !         close(37)
    !     endif
    ! endif

    !運動エネルギー
    ! kinetic_energy_sum = 0.0d0
    ! kinetic_energy_all = 0.0d0
    ! do zi=1,zmax+1
    !     do yi=1,y_procs
    !         do xi=1,x_procs
    !             kinetic_energy_sum = kinetic_energy_sum + 0.5d0*(u1_procs(xi,yi,zi)**2 + u2_procs(xi,yi,zi)**2 + u3_procs(xi,yi,zi)**2)
    !         enddo
    !     enddo
    ! enddo
    ! call MPI_Reduce(kinetic_energy_sum, kinetic_energy_all, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    !Korteweg応力テンソル
    call grad_cal(grad_phi_procs,phi_procs,cr)
    call korteweg_cal(grad_phi_procs,korteweg_procs)
    do zi = 1, zmax+1
        do yi = 1, y_procs
            do xi = 1, x_procs
                tensor11(xi,yi,zi) = korteweg_procs(1,1,xi,yi,zi)
                tensor21(xi,yi,zi) = korteweg_procs(2,1,xi,yi,zi)
                tensor31(xi,yi,zi) = korteweg_procs(3,1,xi,yi,zi)

                tensor12(xi,yi,zi) = korteweg_procs(1,2,xi,yi,zi)
                tensor22(xi,yi,zi) = korteweg_procs(2,2,xi,yi,zi)
                tensor32(xi,yi,zi) = korteweg_procs(3,2,xi,yi,zi)

                tensor13(xi,yi,zi) = korteweg_procs(1,3,xi,yi,zi)
                tensor23(xi,yi,zi) = korteweg_procs(2,3,xi,yi,zi)
                tensor33(xi,yi,zi) = korteweg_procs(3,3,xi,yi,zi)
            enddo
        enddo
    enddo
    call glue(tensor11)
    call glue(tensor21)
    call glue(tensor31)

    call glue(tensor12)
    call glue(tensor22)
    call glue(tensor32)

    call glue(tensor13)
    call glue(tensor23)
    call glue(tensor33)

    ! alpha = 1
    do zi=1,zmax+1
        do yi=1,y_procs
            do xi=1,x_procs
                grad_tensor11(xi,yi,zi) = 0.0d0
                grad_tensor21(xi,yi,zi) = 0.0d0
                grad_tensor31(xi,yi,zi) = 0.0d0
                do i = 2,15
                    grad_tensor11(xi,yi,zi) = grad_tensor11(xi,yi,zi) &
                                            + cr(1,i)*tensor11(xi+cx(i),yi+cy(i),zi+cz(i))
                    grad_tensor21(xi,yi,zi) = grad_tensor21(xi,yi,zi) &
                                            + cr(1,i)*tensor21(xi+cx(i),yi+cy(i),zi+cz(i))
                    grad_tensor31(xi,yi,zi) = grad_tensor31(xi,yi,zi) &
                                            + cr(1,i)*tensor31(xi+cx(i),yi+cy(i),zi+cz(i))
                enddo
                grad_tensor11(xi,yi,zi) = grad_tensor11(xi,yi,zi)/(10.0d0*ds)
                grad_tensor21(xi,yi,zi) = grad_tensor21(xi,yi,zi)/(10.0d0*ds)
                grad_tensor31(xi,yi,zi) = grad_tensor31(xi,yi,zi)/(10.0d0*ds)
            enddo
        enddo
    enddo

    ! alpha = 2
    do zi=1,zmax+1
        do yi=1,y_procs
            do xi=1,x_procs
                grad_tensor12(xi,yi,zi) = 0.0d0
                grad_tensor22(xi,yi,zi) = 0.0d0
                grad_tensor32(xi,yi,zi) = 0.0d0
                do i = 2,15
                    grad_tensor12(xi,yi,zi) = grad_tensor12(xi,yi,zi) &
                                            + cr(2,i)*tensor12(xi+cx(i),yi+cy(i),zi+cz(i))
                    grad_tensor22(xi,yi,zi) = grad_tensor22(xi,yi,zi) &
                                            + cr(2,i)*tensor22(xi+cx(i),yi+cy(i),zi+cz(i))
                    grad_tensor32(xi,yi,zi) = grad_tensor32(xi,yi,zi) &
                                            + cr(2,i)*tensor32(xi+cx(i),yi+cy(i),zi+cz(i))
                enddo
                grad_tensor12(xi,yi,zi) = grad_tensor12(xi,yi,zi)/(10.0d0*ds)
                grad_tensor22(xi,yi,zi) = grad_tensor22(xi,yi,zi)/(10.0d0*ds)
                grad_tensor32(xi,yi,zi) = grad_tensor32(xi,yi,zi)/(10.0d0*ds)
            enddo
        enddo
    enddo

    ! alpha = 3
    do zi=1,zmax+1
        do yi=1,y_procs
            do xi=1,x_procs
                grad_tensor13(xi,yi,zi) = 0.0d0
                grad_tensor23(xi,yi,zi) = 0.0d0
                grad_tensor33(xi,yi,zi) = 0.0d0
                do i = 2,15
                    grad_tensor13(xi,yi,zi) = grad_tensor13(xi,yi,zi) &
                                            + cr(3,i)*tensor13(xi+cx(i),yi+cy(i),zi+cz(i))
                    grad_tensor23(xi,yi,zi) = grad_tensor23(xi,yi,zi) &
                                            + cr(3,i)*tensor23(xi+cx(i),yi+cy(i),zi+cz(i))
                    grad_tensor33(xi,yi,zi) = grad_tensor33(xi,yi,zi) &
                                            + cr(3,i)*tensor33(xi+cx(i),yi+cy(i),zi+cz(i))
                enddo
                grad_tensor13(xi,yi,zi) = grad_tensor13(xi,yi,zi)/(10.0d0*ds)
                grad_tensor23(xi,yi,zi) = grad_tensor23(xi,yi,zi)/(10.0d0*ds)
                grad_tensor33(xi,yi,zi) = grad_tensor33(xi,yi,zi)/(10.0d0*ds)
            enddo
        enddo
    enddo

    ! do zi=1,zmax+1
    !     do yi=1,y_procs
    !         do xi=1,x_procs
    !             interaction_procs(xi,yi,zi) = u1_procs(xi,yi,zi)*grad_tensor11(xi,yi,zi) + u2_procs(xi,yi,zi)*grad_tensor21(xi,yi,zi) + u3_procs(xi,yi,zi)*grad_tensor31(xi,yi,zi) &
    !                                     + u1_procs(xi,yi,zi)*grad_tensor12(xi,yi,zi) + u2_procs(xi,yi,zi)*grad_tensor22(xi,yi,zi) + u3_procs(xi,yi,zi)*grad_tensor32(xi,yi,zi) &
    !                                     + u1_procs(xi,yi,zi)*grad_tensor13(xi,yi,zi) + u2_procs(xi,yi,zi)*grad_tensor23(xi,yi,zi) + u3_procs(xi,yi,zi)*grad_tensor33(xi,yi,zi)
                
    !             interaction_procs(xi,yi,zi) = interaction_procs(xi,yi,zi) * kappag
    !         enddo
    !     enddo
    ! enddo

    !スケール分解して渦の液滴への寄与を見る
    do k_analysis = 1, (xmax+1)/2 + 1
        ! k_analysis = 2
        call scale_cal2(u1_procs,u2_procs,u3_procs,u1_procs_re,u2_procs_re,u3_procs_re,u1_hat,u2_hat,u3_hat,k_analysis,u1_procs_scale,u2_procs_scale,u3_procs_scale)
        call glue(u1_procs)
        call glue(u2_procs)
        call glue(u3_procs)

        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    interaction_procs(xi,yi,zi) = u1_procs_scale(xi,yi,zi)*grad_tensor11(xi,yi,zi) + u2_procs_scale(xi,yi,zi)*grad_tensor21(xi,yi,zi) + u3_procs_scale(xi,yi,zi)*grad_tensor31(xi,yi,zi) &
                                            + u1_procs_scale(xi,yi,zi)*grad_tensor12(xi,yi,zi) + u2_procs_scale(xi,yi,zi)*grad_tensor22(xi,yi,zi) + u3_procs_scale(xi,yi,zi)*grad_tensor32(xi,yi,zi) &
                                            + u1_procs_scale(xi,yi,zi)*grad_tensor13(xi,yi,zi) + u2_procs_scale(xi,yi,zi)*grad_tensor23(xi,yi,zi) + u3_procs_scale(xi,yi,zi)*grad_tensor33(xi,yi,zi)
                    
                    interaction_procs(xi,yi,zi) = interaction_procs(xi,yi,zi) * kappag
                enddo
            enddo
        enddo

        interaction_sum = 0.0d0
        interaction_all = 0.0d0
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    interaction_sum = interaction_sum + interaction_procs(xi,yi,zi)
                enddo
            enddo
        enddo
        call MPI_Reduce(interaction_sum, interaction_all, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

        if(comm_rank == 0) then
            contribution_each_scale(k_analysis) = interaction_all / (dble(xmax+1)*dble(ymax+1)*dble(zmax+1))
            ! contribution_each_scale(k_analysis) = interaction_all
        endif

    enddo


    ! call grad_u_cal(grad_u_procs,u1_procs,u2_procs,u3_procs,cr)
    ! call strain_cal(grad_u_procs,strain_procs)

    ! call interaction_cal(interaction_procs,strain_procs,korteweg_procs)

    ! interaction_sum = 0.0d0
    ! interaction_all = 0.0d0
    ! do zi=1,zmax+1
    !     do yi=1,y_procs
    !         do xi=1,x_procs
    !             interaction_sum = interaction_sum + interaction_procs(xi,yi,zi)
    !         enddo
    !     enddo
    ! enddo
    ! call MPI_Reduce(interaction_sum, interaction_all, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    !自由エネルギー
    ! call free_energy_cal(free_energy,phi_procs)
    ! free_energy_sum = 0.0d0
    ! free_energy_all = 0.0d0
    ! do zi=1,zmax+1
    !     do yi=1,y_procs
    !         do xi=1,x_procs
    !             free_energy_sum = free_energy_sum + (free_energy(xi,yi,zi) + 0.5d0*kappag*(grad_phi_procs(1,xi,yi,zi)**2 + grad_phi_procs(2,xi,yi,zi)**2 + grad_phi_procs(3,xi,yi,zi)**2))
    !             ! free_energy_sum = free_energy_sum + (free_energy(xi,yi,zi))
    !             ! free_energy_sum = free_energy_sum + (0.5d0*kappag*(grad_phi_procs(1,xi,yi,zi)**2 + grad_phi_procs(2,xi,yi,zi)**2 + grad_phi_procs(3,xi,yi,zi)**2))
    !         enddo
    !     enddo
    ! enddo
    ! call MPI_Reduce(free_energy_sum, free_energy_all, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    ! if(comm_rank == 0) then
    !     if(mod(n,1000)==0) then
    !         interaction_all = interaction_all * (-1.0d0) * kappag
    !         if(n == 1000) then
    !             ! open(18,file="./contribution_s.d",status='replace')
    !             ! write(18,"(2es16.8)") dble(n), interaction_all / (dble(xmax+1)*dble(ymax+1)*dble(zmax+1))
    !             ! close(18)

    !             open(19,file="./energy_all.d",status='replace')
    !             write(19,"(3es16.8)") dble(n), free_energy_all, kinetic_energy_all
    !             close(19)
    !         else
    !             ! open(18,file="./contribution_s.d",action="write",position="append")
    !             ! write(18,"(2es16.8)") dble(n), interaction_all / (dble(xmax+1)*dble(ymax+1)*dble(zmax+1))
    !             ! close(18)

    !             open(19,file="./energy_all.d",action="write",position="append")
    !             write(19,"(3es16.8)") dble(n), free_energy_all, kinetic_energy_all
    !             close(19)
    !         endif
    !     endif
    ! endif

    ! call output(interaction_procs)

    if(comm_rank == 0) then
        if(mod(n,1000)==0) then
            write(filename2,*) n
            filename=datadir_output//trim(adjustl(filename2))//'.d' 
            print *, filename !表示してみる
            open(100,file=filename, form='formatted',status='replace')
            do k_analysis=1, (xmax+1)/2+1
                write(100,*) k_analysis, contribution_each_scale(k_analysis)
            enddo
            close(100)

            if(n == 5000) then
                k_analysis = 1
                open(20,file="./energy_1_10.d",status='replace')
                write(20,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
                close(20)
            else
                k_analysis = 1
                open(20,file="./energy_1_10.d",action="write",position="append")
                write(20,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
                close(20)
            endif

            ! if(n == 5000) then
            !     k_analysis = 11
            !     open(21,file="./energy_11_20.d",status='replace')
            !     write(21,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(21)
            ! else
            !     k_analysis = 11
            !     open(21,file="./energy_11_20.d",action="write",position="append")
            !     write(21,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(21)
            ! endif

            ! if(n == 5000) then
            !     k_analysis = 21
            !     open(22,file="./energy_21_30.d",status='replace')
            !     write(22,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(22)
            ! else
            !     k_analysis = 21
            !     open(22,file="./energy_21_30.d",action="write",position="append")
            !     write(22,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(22)
            ! endif

            ! if(n == 5000) then
            !     k_analysis = 31
            !     open(23,file="./energy_31_40.d",status='replace')
            !     write(23,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(23)
            ! else
            !     k_analysis = 31
            !     open(23,file="./energy_31_40.d",action="write",position="append")
            !     write(23,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(23)
            ! endif

            ! if(n == 5000) then
            !     k_analysis = 41
            !     open(24,file="./energy_41_50.d",status='replace')
            !     write(24,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(24)
            ! else
            !     k_analysis = 41
            !     open(24,file="./energy_41_50.d",action="write",position="append")
            !     write(24,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(24)
            ! endif

            ! if(n == 5000) then
            !     k_analysis = 51
            !     open(25,file="./energy_51_60.d",status='replace')
            !     write(25,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(25)
            ! else
            !     k_analysis = 51
            !     open(25,file="./energy_51_60.d",action="write",position="append")
            !     write(25,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(25)
            ! endif

            ! if(n == 5000) then
            !     k_analysis = 61
            !     open(26,file="./energy_61_70.d",status='replace')
            !     write(26,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(26)
            ! else
            !     k_analysis = 61
            !     open(26,file="./energy_61_70.d",action="write",position="append")
            !     write(26,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(26)
            ! endif

            ! if(n == 5000) then
            !     k_analysis = 71
            !     open(27,file="./energy_71_80.d",status='replace')
            !     write(27,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(27)
            ! else
            !     k_analysis = 71
            !     open(27,file="./energy_71_80.d",action="write",position="append")
            !     write(27,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(27)
            ! endif

            ! if(n == 5000) then
            !     k_analysis = 81
            !     open(28,file="./energy_81_90.d",status='replace')
            !     write(28,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(28)
            ! else
            !     k_analysis = 81
            !     open(28,file="./energy_81_90.d",action="write",position="append")
            !     write(28,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(28)
            ! endif

            ! if(n == 5000) then
            !     k_analysis = 91
            !     open(29,file="./energy_91_100.d",status='replace')
            !     write(29,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(29)
            ! else
            !     k_analysis = 91
            !     open(29,file="./energy_91_100.d",action="write",position="append")
            !     write(29,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(29)
            ! endif

            ! if(n == 5000) then
            !     k_analysis = 101
            !     open(30,file="./energy_101_110.d",status='replace')
            !     write(30,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(30)
            ! else
            !     k_analysis = 101
            !     open(30,file="./energy_101_110.d",action="write",position="append")
            !     write(30,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(30)
            ! endif

            ! if(n == 5000) then
            !     k_analysis = 111
            !     open(31,file="./energy_111_120.d",status='replace')
            !     write(31,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(31)
            ! else
            !     k_analysis = 111
            !     open(31,file="./energy_111_120.d",action="write",position="append")
            !     write(31,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(31)
            ! endif

            ! if(n == 5000) then
            !     k_analysis = 121
            !     open(32,file="./energy_121_128.d",status='replace')
            !     write(32,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(32)
            ! else
            !     k_analysis = 121
            !     open(32,file="./energy_121_128.d",action="write",position="append")
            !     write(32,"(11es16.8)") dble(n), contribution_each_scale(k_analysis), contribution_each_scale(k_analysis+1), contribution_each_scale(k_analysis+2), contribution_each_scale(k_analysis+3), contribution_each_scale(k_analysis+4), contribution_each_scale(k_analysis+5), contribution_each_scale(k_analysis+6), contribution_each_scale(k_analysis+7), contribution_each_scale(k_analysis+8), contribution_each_scale(k_analysis+9)
            !     close(32)
            ! endif

        endif
    endif

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
