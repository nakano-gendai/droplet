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
    integer,parameter:: step_start = 452000
    integer,parameter:: step_end = 652000
    integer,parameter:: step_bin = 4000
    integer,parameter:: step_num2 = (step_end - step_start) / step_bin + 1 

    !読み込みディレクトリ
    character(*),parameter :: datadir_input = "/data/sht/nakanog/DNS_turbulence_256_IHT_new/"
    !出力ディレクトリ
    character(*),parameter :: datadir_output = "/data/sht/nakanog/DNS_turbulence_256_IHT/eddy/small/"

    !粒子速度（整数）
    integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
    integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
    integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)
    real(8):: cr(1:3, 1:15)  !粒子速度（実数）

    !無次元数
    real(8),parameter:: eta = 1.0d0 !粘度比（nu2/nu1）

    real(8),parameter:: D = 40.0d0 !設置する液滴直径
    real(8),parameter:: nu1 = 0.001d0 !連続相の粘性係数
    real(8),parameter:: nu2 = eta*nu1 !分散相の粘性係数
    real(8),parameter:: pi = acos(-1.0d0) !円周率

    integer,parameter:: k_high = 27
    integer,parameter:: k_low = 9

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
    subroutine par(cx,cy,cz,cr)
        integer, intent(in):: cx(15),cy(15),cz(15)
        real(8), intent(out):: cr(1:3,1:15)
        do i=1,15
            cr(1,i) = dble(cx(i))
            cr(2,i) = dble(cy(i))
            cr(3,i) = dble(cz(i))
        enddo
    end subroutine par

    subroutine ini(p_procs,u1_procs,u2_procs,u3_procs)
        real(8),allocatable :: p_procs(:,:,:), u1_procs(:,:,:), u2_procs(:,:,:), u3_procs(:,:,:)

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

        !初期化
        p_procs(:,:,:) = 0.0d0
        u1_procs(:,:,:) = 0.0d0
        u2_procs(:,:,:) = 0.0d0
        u3_procs(:,:,:) = 0.0d0
    end subroutine ini

    subroutine ini_fft(u1_hat,u2_hat,u3_hat,u1_procs_re,u2_procs_re,u3_procs_re)
        complex(kind(0d0)), allocatable :: u1_hat(:,:,:), u2_hat(:,:,:), u3_hat(:,:,:)
        real(8),allocatable :: u1_procs_re(:,:,:), u2_procs_re(:,:,:), u3_procs_re(:,:,:)

        call decomp_2d_init(xmax+1, ymax+1, zmax+1, Nx, Ny)
        call decomp_2d_fft_init(PHYSICAL_IN_Z)
        call decomp_2d_fft_get_size(sta, last, sized)

        allocate(u1_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax))
        allocate(u2_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax))
        allocate(u3_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax))
        allocate(u1_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1))
        allocate(u2_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1))
        allocate(u3_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1))

        u1_procs_re(:,:,:) = 0.0d0
        u2_procs_re(:,:,:) = 0.0d0
        u3_procs_re(:,:,:) = 0.0d0
        u1_hat(:,:,:) = (0.0d0, 0.0d0)
        u2_hat(:,:,:) = (0.0d0, 0.0d0)
        u3_hat(:,:,:) = (0.0d0, 0.0d0)
    end subroutine ini_fft

    subroutine variable(grad_u_procs,vorticity_procs,enstrophy_procs,u1_procs_ave,u2_procs_ave,u3_procs_ave,u1_procs_fluctuation,u2_procs_fluctuation,u3_procs_fluctuation,strain_procs,u_rms,u_grad,vtensor_procs,Q_procs)
        real(8),allocatable :: grad_u_procs(:,:,:,:,:)
        real(8),allocatable :: vorticity_procs(:,:,:,:), enstrophy_procs(:,:,:)
        real(8),allocatable :: u1_procs_ave(:,:,:), u2_procs_ave(:,:,:), u3_procs_ave(:,:,:)
        real(8),allocatable :: u1_procs_fluctuation(:,:,:), u2_procs_fluctuation(:,:,:), u3_procs_fluctuation(:,:,:)
        real(8),allocatable :: strain_procs(:,:,:,:,:), vtensor_procs(:,:,:,:,:)
        real(8),allocatable :: u_rms(:,:,:), u_grad(:,:,:)
        real(8),allocatable :: Q_procs(:,:,:)

        allocate(grad_u_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1))
        allocate(vorticity_procs(1:3,1:x_procs,1:y_procs,1:zmax+1))
        allocate(enstrophy_procs(1:x_procs,1:y_procs,1:zmax+1))

        allocate(u1_procs_ave(1:x_procs,1:y_procs,1:zmax+1))
        allocate(u2_procs_ave(1:x_procs,1:y_procs,1:zmax+1))
        allocate(u3_procs_ave(1:x_procs,1:y_procs,1:zmax+1))

        allocate(u1_procs_fluctuation(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(u2_procs_fluctuation(0:x_procs+1,0:y_procs+1,0:zmax+2))
        allocate(u3_procs_fluctuation(0:x_procs+1,0:y_procs+1,0:zmax+2))

        allocate(strain_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1))
        allocate(vtensor_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1))

        allocate(u_rms(1:x_procs,1:y_procs,1:zmax+1))
        allocate(u_grad(1:x_procs,1:y_procs,1:zmax+1))

        allocate(Q_procs(1:x_procs,1:y_procs,1:zmax+1))

        grad_u_procs(:,:,:,:,:) = 0.0d0
        vorticity_procs(:,:,:,:) = 0.0d0
        enstrophy_procs(:,:,:) = 0.0d0
        u1_procs_ave(:,:,:) = 0.0d0
        u2_procs_ave(:,:,:) = 0.0d0
        u3_procs_ave(:,:,:) = 0.0d0
        u1_procs_fluctuation(:,:,:) = 0.0d0
        u2_procs_fluctuation(:,:,:) = 0.0d0
        u3_procs_fluctuation(:,:,:) = 0.0d0
        strain_procs(:,:,:,:,:) = 0.0d0
        vtensor_procs(:,:,:,:,:) = 0.0d0
        u_rms(:,:,:) = 0.0d0
        u_grad(:,:,:) = 0.0d0
        Q_procs(:,:,:) = 0.0d0
    end subroutine variable

    subroutine input(u1_procs,u2_procs,u3_procs,p_procs)
        real(8),intent(inout):: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: p_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8) s

        !==============================データ読み込み=================================================
        write(filename2,*) n
        write(filename3,*) comm_rank
        filename=datadir_input//'0_'//trim(adjustl(filename3))//'_'//trim(adjustl(filename2))//'.bin' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
        print *, filename !表示してみる
        open(10, file=filename, form="unformatted")
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    read(10) u1_procs(xi,yi,zi), u2_procs(xi,yi,zi), u3_procs(xi,yi,zi), p_procs(xi,yi,zi)
                    ! read(10) s, u1_procs(xi,yi,zi), u2_procs(xi,yi,zi), u3_procs(xi,yi,zi), p_procs(xi,yi,zi)
                enddo
            enddo
        enddo
        close(10)
    end subroutine input

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

    subroutine grad_u_cal(grad_u_procs,u1_procs,u2_procs,u3_procs)
        real(8),intent(out):: grad_u_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        integer alpha, i, xi, yi, zi

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
                            strain_procs(alpha,beta,xi,yi,zi) = grad_u_procs(alpha,beta,xi,yi,zi) + grad_u_procs(beta,alpha,xi,yi,zi)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine strain_cal

    subroutine vtensor_cal(grad_u_procs,vtensor_procs)
        real(8),intent(in):: grad_u_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(out):: vtensor_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)

        do zi = 1, zmax + 1
            do yi = 1, y_procs
                do xi = 1, x_procs
                    do beta = 1, 3
                        do alpha = 1, 3
                            vtensor_procs(alpha,beta,xi,yi,zi) = grad_u_procs(alpha,beta,xi,yi,zi) - grad_u_procs(beta,alpha,xi,yi,zi)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine vtensor_cal

    subroutine u_rms_cal(u_rms,u_rms_ave,u1_procs_fluctuation,u2_procs_fluctuation,u3_procs_fluctuation)
        real(8),intent(out):: u_rms(1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(inout):: u_rms_ave
        real(8),intent(in):: u1_procs_fluctuation(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u2_procs_fluctuation(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u3_procs_fluctuation(0:x_procs+1,0:y_procs+1,0:zmax+2)

        do zi = 1, zmax + 1
            do yi = 1, y_procs
                do xi = 1, x_procs
                    u_rms(xi,yi,zi) = (u1_procs_fluctuation(xi,yi,zi)**2.0d0 &
                                    + u2_procs_fluctuation(xi,yi,zi)**2.0d0 &
                                    + u3_procs_fluctuation(xi,yi,zi)**2.0d0)
                enddo
            enddo
        enddo
        u_rms_ave = 0.0d0
        do zi = 1, zmax + 1
            do yi = 1, y_procs
                do xi = 1, x_procs
                    u_rms_ave = u_rms_ave + u_rms(xi,yi,zi)
                enddo
            enddo
        enddo
    end subroutine u_rms_cal

    subroutine u_grad_cal(u_grad,u_grad_ave,grad_u_procs)
        real(8),intent(out):: u_grad(1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(inout):: u_grad_ave
        real(8),intent(in):: grad_u_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)

        do zi = 1, zmax + 1
            do yi = 1, y_procs
                do xi = 1, x_procs
                    u_grad(xi,yi,zi) = (grad_u_procs(1,1,xi,yi,zi)**2.0d0 &
                                    + grad_u_procs(2,2,xi,yi,zi)**2.0d0 &
                                    + grad_u_procs(3,3,xi,yi,zi)**2.0d0)
                enddo
            enddo
        enddo
        u_grad_ave = 0.0d0
        do zi = 1, zmax + 1
            do yi = 1, y_procs
                do xi = 1, x_procs
                    u_grad_ave = u_grad_ave + u_grad(xi,yi,zi)
                enddo
            enddo
        enddo
    end subroutine u_grad_cal

    subroutine taylor_re_cal(u_rms_ave,u_rms_ave_sum,u_grad_ave,u_grad_ave_sum,taylor_re,taylor_length,u_rms,u_grad)
        real(8),intent(inout):: u_rms_ave
        real(8),intent(inout):: u_grad_ave

        real(8),intent(in):: u_rms(1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: u_grad(1:x_procs,1:y_procs,1:zmax+1)

        real(8),intent(inout):: u_rms_ave_sum
        real(8),intent(inout):: u_grad_ave_sum
        real(8),intent(out):: taylor_re
        real(8),intent(out):: taylor_length

        u_rms_ave_sum = 0.0d0
        u_grad_ave_sum = 0.0d0

        call MPI_Reduce(u_rms_ave, u_rms_ave_sum, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_Reduce(u_grad_ave, u_grad_ave_sum, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        
        if(comm_rank == 0) then
            u_rms_ave_sum = ( u_rms_ave_sum / (dble(xmax+1)*dble(ymax+1)*dble(zmax+1)) )**0.5d0
            u_grad_ave_sum = ( u_grad_ave_sum / (dble(xmax+1)*dble(ymax+1)*dble(zmax+1)) )**0.5d0

            taylor_length = u_rms_ave_sum / u_grad_ave_sum
            taylor_re = u_rms_ave_sum * taylor_length / nu1

            if(n == step_start) then
                open(103,file="./taylor_para.d")
                write(103,*) dble(n), dble(n)*umax/D_vortex, taylor_length, taylor_re, u_rms_ave_sum
                close(103)
            else
                open(103,file="./taylor_para.d",action="write",position="append")
                write(103,*) dble(n), dble(n)*umax/D_vortex, taylor_length, taylor_re, u_rms_ave_sum
                close(103)
            endif
        endif
    end subroutine taylor_re_cal

    subroutine energy_dissipation_cal(strain_procs,energy_dissipation,energy_dissipation_sum,Kolmogorov_scale)
        real(8),intent(in):: strain_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(out):: energy_dissipation, energy_dissipation_sum, Kolmogorov_scale

        energy_dissipation = 0.0d0
        energy_dissipation_sum = 0.0d0
        Kolmogorov_scale = 0.0d0
        do zi = 1, zmax + 1
            do yi = 1, y_procs
                do xi = 1, x_procs
                    do beta = 1, 3
                        do alpha = 1, 3
                            energy_dissipation_sum = energy_dissipation_sum + strain_procs(alpha,beta,xi,yi,zi)**2
                        enddo
                    enddo
                enddo
            enddo
        enddo

        call MPI_Reduce(energy_dissipation_sum, energy_dissipation, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        
        if(comm_rank == 0) then
            energy_dissipation = energy_dissipation / ( dble(xmax+1)*dble(ymax+1)*dble(zmax+1) )
            energy_dissipation = energy_dissipation * 0.5d0 * nu1

            Kolmogorov_scale = nu1**(3.0d0/4.0d0) * energy_dissipation**(-1.0d0/4.0d0)

            if(n == step_start) then
                open(102,file="./kolmogorov.d")
                write(102,*) dble(n), dble(n)*umax/D_vortex, Kolmogorov_scale, energy_dissipation, nu1
                close(102)
            else
                open(102,file="./kolmogorov.d",action="write",position="append")
                write(102,*) dble(n), dble(n)*umax/D_vortex, Kolmogorov_scale, energy_dissipation, nu1
                close(102)
            endif
        endif
    end subroutine energy_dissipation_cal

    subroutine u_sum_step_cal(u1_procs,u2_procs,u3_procs,u1_procs_ave,u2_procs_ave,u3_procs_ave)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(out):: u1_procs_ave(1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(out):: u2_procs_ave(1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(out):: u3_procs_ave(1:x_procs,1:y_procs,1:zmax+1)

        do zi = 1, zmax + 1
            do yi = 1, y_procs
                do xi = 1, x_procs
                    u1_procs_ave(xi,yi,zi) = u1_procs_ave(xi,yi,zi) + u1_procs(xi,yi,zi)
                    u2_procs_ave(xi,yi,zi) = u2_procs_ave(xi,yi,zi) + u2_procs(xi,yi,zi)
                    u3_procs_ave(xi,yi,zi) = u3_procs_ave(xi,yi,zi) + u3_procs(xi,yi,zi)
                enddo
            enddo
        enddo
    end subroutine u_sum_step_cal

    subroutine vorticity_cal(grad_u_procs,vorticity_procs)
        real(8),intent(in):: grad_u_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(out):: vorticity_procs(1:3,1:x_procs,1:y_procs,1:zmax+1)
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    vorticity_procs(1,xi,yi,zi) = grad_u_procs(3,2,xi,yi,zi) - grad_u_procs(2,3,xi,yi,zi)
                    vorticity_procs(2,xi,yi,zi) = grad_u_procs(1,3,xi,yi,zi) - grad_u_procs(3,1,xi,yi,zi)
                    vorticity_procs(3,xi,yi,zi) = grad_u_procs(2,1,xi,yi,zi) - grad_u_procs(1,2,xi,yi,zi)
                enddo
            enddo
        enddo
    end subroutine vorticity_cal

    subroutine enstrophy_cal(vorticity_procs,enstrophy_procs)
        real(8),intent(in):: vorticity_procs(1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(out):: enstrophy_procs(1:x_procs,1:y_procs,1:zmax+1)
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    enstrophy_procs(xi,yi,zi) = 0.5d0 * (vorticity_procs(1,xi,yi,zi)**2 + &
                                                        vorticity_procs(2,xi,yi,zi)**2 + &
                                                        vorticity_procs(3,xi,yi,zi)**2)
                enddo
            enddo
        enddo
    end subroutine enstrophy_cal

    subroutine Q_cal(Q_procs,strain_procs,vtensor_procs)
        real(8),intent(inout):: Q_procs(1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: strain_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8),intent(in):: vtensor_procs(1:3,1:3,1:x_procs,1:y_procs,1:zmax+1)
        real(8) ftmp, gtmp

        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    ftmp = 0.0d0
                    gtmp = 0.0d0
                    
                    ftmp = strain_procs(1,1,xi,yi,zi)*strain_procs(2,2,xi,yi,zi)*strain_procs(3,3,xi,yi,zi) &
                            + strain_procs(1,2,xi,yi,zi)*strain_procs(3,3,xi,yi,zi)*strain_procs(3,1,xi,yi,zi) &
                            + strain_procs(1,3,xi,yi,zi)*strain_procs(2,1,xi,yi,zi)*strain_procs(3,2,xi,yi,zi) &
                            - strain_procs(1,3,xi,yi,zi)*strain_procs(2,2,xi,yi,zi)*strain_procs(3,1,xi,yi,zi) &
                            - strain_procs(1,1,xi,yi,zi)*strain_procs(2,3,xi,yi,zi)*strain_procs(3,2,xi,yi,zi) &
                            - strain_procs(1,2,xi,yi,zi)*strain_procs(2,1,xi,yi,zi)*strain_procs(3,3,xi,yi,zi)
                    
                    gtmp = vtensor_procs(1,1,xi,yi,zi)*vtensor_procs(2,2,xi,yi,zi)*vtensor_procs(3,3,xi,yi,zi) &
                            + vtensor_procs(1,2,xi,yi,zi)*vtensor_procs(3,3,xi,yi,zi)*vtensor_procs(3,1,xi,yi,zi) &
                            + vtensor_procs(1,3,xi,yi,zi)*vtensor_procs(2,1,xi,yi,zi)*vtensor_procs(3,2,xi,yi,zi) &
                            - vtensor_procs(1,3,xi,yi,zi)*vtensor_procs(2,2,xi,yi,zi)*vtensor_procs(3,1,xi,yi,zi) &
                            - vtensor_procs(1,1,xi,yi,zi)*vtensor_procs(2,3,xi,yi,zi)*vtensor_procs(3,2,xi,yi,zi) &
                            - vtensor_procs(1,2,xi,yi,zi)*vtensor_procs(2,1,xi,yi,zi)*vtensor_procs(3,3,xi,yi,zi)
                    Q_procs(xi,yi,zi) = gtmp*gtmp - ftmp*ftmp
                enddo
            enddo
        enddo
    end subroutine Q_cal

    subroutine output(Q_procs)
        real(8),intent(inout):: Q_procs(1:x_procs,1:y_procs,1:zmax+1)

        write(filename3,*) comm_rank
        write(filename2,*) n
        filename=datadir_output//'0_'//trim(adjustl(filename3))//'_'//trim(adjustl(filename2))//'.bin' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
        print *, filename !表示してみる
        open(101,file=filename, form='unformatted',status='replace') 
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    write(101) Q_procs(xi,yi,zi)
                enddo
            enddo
        enddo
        close(101)
    end subroutine output

    subroutine output_ens(enstrophy_procs)
        real(8),intent(inout):: enstrophy_procs(1:x_procs,1:y_procs,1:zmax+1)

        write(filename3,*) comm_rank
        write(filename2,*) n
        filename=datadir_output//'0_'//trim(adjustl(filename3))//'_'//trim(adjustl(filename2))//'.bin' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
        print *, filename !表示してみる
        open(101,file=filename, form='unformatted',status='replace') 
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    write(101) enstrophy_procs(xi,yi,zi)
                enddo
            enddo
        enddo
        close(101)
    end subroutine output_ens
    
    !ディレクトリ作成
    subroutine mk_dirs(outdir)
        implicit none
        character(len=*), intent(in) :: outdir
        character(len=256) command
        write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
        call system(command)
    end subroutine mk_dirs

    subroutine enesupe_cal(u1_procs,u2_procs,u3_procs,u1_procs_re,u2_procs_re,u3_procs_re,u1_hat,u2_hat,u3_hat,enesupe,enesupe_sum)
        real(8),intent(inout):: u1_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u2_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u3_procs(0:x_procs+1,0:y_procs+1,0:zmax+2)
        real(8),intent(inout):: u1_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax)
        real(8),intent(inout):: u2_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax)
        real(8),intent(inout):: u3_procs_re(0:x_procs-1, 0:y_procs-1, 0:zmax)
        complex(kind(0d0)), intent(inout) :: u1_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)
        complex(kind(0d0)), intent(inout) :: u2_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)
        complex(kind(0d0)), intent(inout) :: u3_hat(sta(1)-1:last(1)-1, sta(2)-1:last(2)-1, sta(3)-1:last(3)-1)
        real(8),intent(inout):: enesupe((xmax+1)/2 + 1), enesupe_sum((xmax+1)/2 + 1)
        integer k1, k2, k3, k(1:3), k_index
        real(8) k_abs
        real(8) dkx, dky, dkz, dk

        dkx = 2.0d0*pi/dble(xmax+1)
        dky = 2.0d0*pi/dble(ymax+1)
        dkz = 2.0d0*pi/dble(zmax+1)
        dk = ( (dkx)**2.0d0 + (dky)**2.0d0 + (dkz)**2.0d0 )**0.5d0
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

        enesupe(:) = 0.0d0
        enesupe_sum(:) = 0.0d0
        do k3=sta(3)-1,last(3)-1
            k(3) =( k3 - judge(k3, zmax+1) )
            do k2=sta(2)-1,last(2)-1
                k(2) = ( k2 - judge(k2, ymax+1) )
                do k1=sta(1)-1,last(1)-1
                    k(1) = ( k1 - judge(k1, xmax+1) )

                    k_abs = ( dble(k(1)**2 + k(2)**2 + k(3)**2) )**0.5d0
                    k_index = int( k_abs ) + 1

                    if(k3 == 0) then
                        enesupe(k_index) = enesupe(k_index) + 0.5d0 * (abs(u1_hat(k1,k2,k3))**2 + abs(u2_hat(k1,k2,k3))**2 + abs(u3_hat(k1,k2,k3))**2)
                    else 
                        enesupe(k_index) = enesupe(k_index) + 0.5d0 * (abs(u1_hat(k1,k2,k3))**2 + abs(u2_hat(k1,k2,k3))**2 + abs(u3_hat(k1,k2,k3))**2) * 2.0d0
                    endif
                enddo
            enddo
        enddo

        call MPI_Reduce(enesupe(1), enesupe_sum(1), (xmax+1)/2 + 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    end subroutine enesupe_cal

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
    real(8),allocatable :: p_procs(:,:,:), u1_procs(:,:,:), u2_procs(:,:,:), u3_procs(:,:,:)
    real(8),allocatable :: grad_u_procs(:,:,:,:,:)
    real(8),allocatable :: vorticity_procs(:,:,:,:), enstrophy_procs(:,:,:)
    real(8) enstrophy_sum, enstrophy_ave
    real(8),allocatable :: u1_procs_ave(:,:,:), u2_procs_ave(:,:,:), u3_procs_ave(:,:,:)
    real(8),allocatable :: u1_procs_fluctuation(:,:,:), u2_procs_fluctuation(:,:,:), u3_procs_fluctuation(:,:,:)
    real(8),allocatable :: strain_procs(:,:,:,:,:), vtensor_procs(:,:,:,:,:)
    real(8) Kolmogorov_scale, energy_dissipation, energy_dissipation_sum
    real(8) u_abs, uave_abs, uflu_abs
    real(8),allocatable :: u_rms(:,:,:), u_grad(:,:,:)
    real(8) u_rms_ave, u_grad_ave, u_rms_ave_sum, u_grad_ave_sum,taylor_re, taylor_length
    real(8),allocatable ::  kinetic_procs(:,:,:)
    real(8) kinetic_sum, kinetic_ave
    real(8),allocatable :: Q_procs(:,:,:)
    real(8) integral_scale, integral_top, integral_bottom

    !波数空間の変数
    complex(kind(0d0)), allocatable :: u1_hat(:,:,:), u2_hat(:,:,:), u3_hat(:,:,:)
    real(8),allocatable :: u1_procs_re(:,:,:), u2_procs_re(:,:,:), u3_procs_re(:,:,:)
    real(8),allocatable :: enesupe(:), enesupe_sum(:), enesupe_result(:)
    real(8) dkx, dky, dkz, dk, kmax_abs
!==================================MPI並列開始===============================================
    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, comm_procs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, comm_rank, ierr)
    call cpu_time(time1)
!================================ディレクトリの作成============================================
    call mk_dirs(datadir_output)
    call par(cx,cy,cz,cr)
    call ini(p_procs,u1_procs,u2_procs,u3_procs)
    call ini_fft(u1_hat,u2_hat,u3_hat,u1_procs_re,u2_procs_re,u3_procs_re)
    call variable(grad_u_procs,vorticity_procs,enstrophy_procs,u1_procs_ave,u2_procs_ave,u3_procs_ave,u1_procs_fluctuation,u2_procs_fluctuation,u3_procs_fluctuation,strain_procs,u_rms,u_grad,vtensor_procs,Q_procs)
    allocate(kinetic_procs(1:x_procs,1:y_procs,1:zmax+1))
    allocate(enesupe((xmax+1)/2 + 1))
    allocate(enesupe_sum((xmax+1)/2 + 1))
    allocate(enesupe_result((xmax+1)/2 + 1))
    enesupe(:) = 0.0d0
    enesupe_sum(:) = 0.0d0
    enesupe_result(:) = 0.0d0
!==============時間平均の計算==================
    ! step_num = 0
    ! DO n = step_start, step_end, step_bin
    !     step_num = step_num + 1
    !     call input(u1_procs,u2_procs,u3_procs,p_procs)
    !     call u_sum_step_cal(u1_procs,u2_procs,u3_procs,u1_procs_ave,u2_procs_ave,u3_procs_ave)
    ! ENDDO
    ! do zi = 1, zmax + 1 
    !     do yi = 1, y_procs
    !         do xi = 1, x_procs
    !             u1_procs_ave(xi,yi,zi) = u1_procs_ave(xi,yi,zi) / (dble(step_num))
    !             u2_procs_ave(xi,yi,zi) = u2_procs_ave(xi,yi,zi) / (dble(step_num))
    !             u3_procs_ave(xi,yi,zi) = u3_procs_ave(xi,yi,zi) / (dble(step_num))
    !         enddo
    !     enddo
    ! enddo
!================================物理量計算================================
    DO n = step_start, step_end, step_bin
        !入力ファイル読み込み
        call input(u1_procs,u2_procs,u3_procs,p_procs)
        !変動速度の計算
        u1_procs_fluctuation(:,:,:) = 0.0d0
        u2_procs_fluctuation(:,:,:) = 0.0d0
        u3_procs_fluctuation(:,:,:) = 0.0d0
        do zi = 1, zmax + 1 
            do yi = 1, y_procs
                do xi = 1, x_procs
                    ! u1_procs_fluctuation(xi,yi,zi) = u1_procs(xi,yi,zi) - u1_procs_ave(xi,yi,zi)
                    ! u2_procs_fluctuation(xi,yi,zi) = u2_procs(xi,yi,zi) - u2_procs_ave(xi,yi,zi)
                    ! u3_procs_fluctuation(xi,yi,zi) = u3_procs(xi,yi,zi) - u3_procs_ave(xi,yi,zi)
                    u1_procs_fluctuation(xi,yi,zi) = u1_procs(xi,yi,zi)
                    u2_procs_fluctuation(xi,yi,zi) = u2_procs(xi,yi,zi)
                    u3_procs_fluctuation(xi,yi,zi) = u3_procs(xi,yi,zi)
                enddo
            enddo
        enddo
        ! !コルモゴロフスケールの計算
        ! call glue(u1_procs_fluctuation)
        ! call glue(u2_procs_fluctuation)
        ! call glue(u3_procs_fluctuation)
        ! call grad_u_cal(grad_u_procs,u1_procs_fluctuation,u2_procs_fluctuation,u3_procs_fluctuation)
        ! call strain_cal(grad_u_procs,strain_procs)
        ! call energy_dissipation_cal(strain_procs,energy_dissipation,energy_dissipation_sum,Kolmogorov_scale)

        ! !テイラー長レイノルズ数の計算
        ! call u_rms_cal(u_rms,u_rms_ave,u1_procs_fluctuation,u2_procs_fluctuation,u3_procs_fluctuation)
        ! call u_grad_cal(u_grad,u_grad_ave,grad_u_procs)
        ! call taylor_re_cal(u_rms_ave,u_rms_ave_sum,u_grad_ave,u_grad_ave_sum,taylor_re,taylor_length,u_rms,u_grad)

        !エネルギースペクトル
        call enesupe_cal(u1_procs,u2_procs,u3_procs,u1_procs_re,u2_procs_re,u3_procs_re,u1_hat,u2_hat,u3_hat,enesupe,enesupe_sum)
        if(comm_rank == 0) then
            do i = 1, (xmax+1)/2 + 1
                enesupe_result(i) = enesupe_result(i) + enesupe_sum(i)
            enddo
        endif
    ENDDO

    if(comm_rank==0) then
        open(37,file ="./enesupe.d")
        do i=1,(xmax+1)/2+1
            enesupe_result(i) = enesupe_result(i) / dble(step_num2)
            write(37,"(2es16.8)") dble(i)-1.0d0, enesupe_result(i)
        enddo
        close(37)
    endif

    ! !積分長
    ! if(comm_rank == 0) then
    !     integral_top = 0.0d0
    !     integral_bottom = 0.0d0
    !     do i = 2, (xmax+1)/2 + 1
    !         integral_top = integral_top + enesupe_result(i) / (dble(i) - 1.0d0)
    !         integral_bottom = integral_bottom + enesupe_result(i)
    !     enddo
    !     integral_scale = 3.0d0 / 4.0d0 * pi * integral_top / integral_bottom
    !     open(38,file="./integral_scal.d")
    !     write(38,*) integral_scale
    !     close(38)
    ! endif

    !=======エンストロフィー＝＝＝＝＝＝＝＝＝＝＝＝＝
    ! DO n = step_start, step_end, step_bin
        ! n = 400000
        ! call input(u1_procs,u2_procs,u3_procs,p_procs)
        ! call scale_cal(u1_procs,u2_procs,u3_procs,u1_procs_re,u2_procs_re,u3_procs_re,u1_hat,u2_hat,u3_hat)
        ! call glue(u1_procs)
        ! call glue(u2_procs)
        ! call glue(u3_procs)
        ! call grad_u_cal(grad_u_procs, u1_procs, u2_procs, u3_procs)
        ! call vorticity_cal(grad_u_procs,vorticity_procs)
        ! call enstrophy_cal(vorticity_procs,enstrophy_procs)
        ! enstrophy_sum = 0.0d0
        ! do zi=1,zmax+1
        !     do yi=1,y_procs
        !         do xi=1,x_procs
        !             enstrophy_sum = enstrophy_sum + enstrophy_procs(xi,yi,zi)
        !         enddo
        !     enddo
        ! enddo
        ! call MPI_Reduce(enstrophy_sum, enstrophy_ave, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        ! if(comm_rank == 0) then
        !     enstrophy_ave = enstrophy_ave / (dble(xmax+1)*dble(ymax+1)*dble(zmax+1))
        !     open(100,file = "./enstrophy_ave_small.d",action="write",position="append")
        !     write(100,*) enstrophy_ave
        !     close(100)
        ! endif
        ! call output(enstrophy_procs)
    ! ENDDO

    !Q
    ! DO n = 5000, 100000, 1000
        ! n = 3000000
        ! call input(u1_procs,u2_procs,u3_procs,p_procs)
        ! call glue(u1_procs)
        ! call glue(u2_procs)
        ! call glue(u3_procs)
        ! call grad_u_cal(grad_u_procs,u1_procs,u2_procs,u3_procs)
        ! call strain_cal(grad_u_procs,strain_procs)
        ! call vtensor_cal(grad_u_procs,vtensor_procs)
        ! call Q_cal(Q_procs,strain_procs,vtensor_procs)
        ! call output(Q_procs)
    ! ENDDO
!================MPI並列終わり=======================================
    call decomp_2d_fft_finalize
    call decomp_2d_finalize
    call MPI_Finalize(ierr)
end program main
