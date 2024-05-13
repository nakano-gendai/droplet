! 【Multi-times simulations】
! このプログラムは、多数回の液滴分裂の数値シミュレーションを実行する。
! このプログラムは、界面を「Cahn-Hilliard 方程式」で表現している。
! 液滴を長時間保持できないが、系全体の保存性には優れている。
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module globals
    include "mpif.h"
    !計算領域
    real(8),parameter:: ds = 1.0d0 !格子間隔（lattice unit）
    integer,parameter:: xmax = 511 !ｘ方向格子数（０から数える）
    integer,parameter:: ymax = 511 !ｙ方向格子数（０から数える）
    integer,parameter:: zmax = 511 !ｚ方向格子数（０から数える）
    integer,parameter:: step_start = 5000
    integer,parameter:: step_end = 3800000
    !並行して計算する数
    integer,parameter:: Nxall = 1 !x方向の分割数（Nxall*Nzall=全体の計算数）
    integer,parameter:: Nzall = 1 !z方向の分割数（Nxall*Nzall=全体の計算数）
    integer,parameter:: xall = (xmax + 1) * Nxall !全体のx方向格子数
    integer,parameter:: zall = (zmax + 1) * Nzall !全体のz方向格子数
    !読み込みディレクトリ
    ! character(*),parameter :: datadir_input = "/data/n/n517/taylor_re1000_255_ran/"
    ! !出力ディレクトリ
    ! character(*),parameter :: datadir_output = "/data/n/n517/taylor_re1000_255_ran/output/"

    !読み込みディレクトリ
    character(*),parameter :: datadir_input = "/data/sht/nakanog/taylor_512_dsd/"
    !出力ディレクトリ
    character(*),parameter :: datadir_output = "/data/sht/nakanog/DNS_large_256_re80000/eddy/"

    !粒子速度（整数）
    integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
    integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
    integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)
    real(8):: cr(1:3, 1:15)  !粒子速度（実数）

    !無次元数
    real(8),parameter:: We = 3.0d0 !粒子ウエーバー数
    real(8),parameter:: Re = 160000.0d0 !粒子レイノルズ数
    real(8),parameter:: eta = 1.0d0 !粘度比（nu2/nu1）

    real(8),parameter:: D = 255.5d0 !設置する液滴直径
    real(8),parameter:: D_vortex = 255.5d0 !渦の大きさ
    real(8),parameter:: umax = 0.1d0 !最大流速
    real(8),parameter:: nu1 = umax*D_vortex/Re !連続相の粘性係数
    real(8),parameter:: nu2 = eta*nu1 !分散相の粘性係数

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
    integer Nx, Nz 
    integer x_procs, z_procs
    integer key_new, group_new, key_x, group_x, key_z, group_z
    integer new_comm_world, new_procs, new_rank
    integer newx_comm_world, newx_procs, newx_rank
    integer newz_comm_world, newz_procs, newz_rank
    integer rectangle_type, xtype
    integer next_rank_x, former_rank_x, next_rank_z, former_rank_z
    integer req1s, req1r, req2s, req2r
    integer, dimension(MPI_STATUS_SIZE) :: sta1s, sta1r, sta2s, sta2r
    integer Nxx, Nzz
    integer tags, tagr, recv_rank
    character(2) chmyrank

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

        !z世界での受け渡しをする際の型作成
        call MPI_Type_Vector((z_procs+2)*(ymax+3),1,x_procs+2,MPI_REAL8,rectangle_type,ierr)
        call MPI_Type_Commit(rectangle_type,ierr)
    !===============================================================================================================
        !以下はのりしろ有りの変数
        allocate(p_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))

        !初期化
        p_procs(:,:,:) = 0.0d0
        u1_procs(:,:,:) = 0.0d0
        u2_procs(:,:,:) = 0.0d0
        u3_procs(:,:,:) = 0.0d0
    end subroutine ini

    subroutine variable(grad_u_procs,vorticity_procs,enstrophy_procs,u1_procs_ave,u2_procs_ave,u3_procs_ave,u1_procs_fluctuation,u2_procs_fluctuation,u3_procs_fluctuation,strain_procs,u_rms,u_grad,vtensor_procs,Q_procs)
        real(8),allocatable :: grad_u_procs(:,:,:,:,:)
        real(8),allocatable :: vorticity_procs(:,:,:,:), enstrophy_procs(:,:,:)
        real(8),allocatable :: u1_procs_ave(:,:,:), u2_procs_ave(:,:,:), u3_procs_ave(:,:,:)
        real(8),allocatable :: u1_procs_fluctuation(:,:,:), u2_procs_fluctuation(:,:,:), u3_procs_fluctuation(:,:,:)
        real(8),allocatable :: strain_procs(:,:,:,:,:), vtensor_procs(:,:,:,:,:)
        real(8),allocatable :: u_rms(:,:,:), u_grad(:,:,:)
        real(8),allocatable :: Q_procs(:,:,:)

        allocate(grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs))
        allocate(vorticity_procs(1:3,1:x_procs,1:ymax+1,1:z_procs))
        allocate(enstrophy_procs(1:x_procs,1:ymax+1,1:z_procs))

        allocate(u1_procs_ave(1:x_procs,1:ymax+1,1:z_procs))
        allocate(u2_procs_ave(1:x_procs,1:ymax+1,1:z_procs))
        allocate(u3_procs_ave(1:x_procs,1:ymax+1,1:z_procs))

        allocate(u1_procs_fluctuation(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(u2_procs_fluctuation(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(u3_procs_fluctuation(0:x_procs+1,0:ymax+2,0:z_procs+1))

        allocate(strain_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs))
        allocate(vtensor_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs))

        allocate(u_rms(1:x_procs,1:ymax+1,1:z_procs))
        allocate(u_grad(1:x_procs,1:ymax+1,1:z_procs))

        allocate(Q_procs(1:x_procs,1:ymax+1,1:z_procs))

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
        real(8),intent(inout):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: p_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8) s

        !==============================データ読み込み=================================================
        write(filename,*) group_new !i->filename 変換
        write(filename2,*) n
        write(filename3,*) new_rank
        filename=datadir_input//trim(adjustl(filename))//'_'//trim(adjustl(filename3))//'_'//trim(adjustl(filename2))//'.bin' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
        print *, filename !表示してみる
        open(10, file=filename, form="unformatted")
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    read(10) s, u1_procs(xi,yi,zi), u2_procs(xi,yi,zi), u3_procs(xi,yi,zi), p_procs(xi,yi,zi)
                    ! read(10) u1_procs(xi,yi,zi), u2_procs(xi,yi,zi), u3_procs(xi,yi,zi), p_procs(xi,yi,zi)
                enddo
            enddo
        enddo
        close(10)
    end subroutine input

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

    subroutine strain_cal(grad_u_procs,strain_procs)
        real(8),intent(in):: grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(out):: strain_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)

        do zi = 1, z_procs
            do yi = 1, ymax+1
                do xi = 1, x_procs
                    do beta = 1, 3
                        do alpha = 1, 3
                            strain_procs(alpha,beta,xi,yi,zi) = grad_u_procs(alpha,beta,xi,yi,zi) + grad_u_procs(beta,alpha,xi,yi,zi)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine strain_cal

    subroutine vtensor_cal(grad_u_procs,vtensor_procs)
        real(8),intent(in):: grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(out):: vtensor_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)

        do zi = 1, z_procs
            do yi = 1, ymax+1
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
        real(8),intent(out):: u_rms(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(inout):: u_rms_ave
        real(8),intent(in):: u1_procs_fluctuation(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u2_procs_fluctuation(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u3_procs_fluctuation(0:x_procs+1,0:ymax+2,0:z_procs+1)

        do zi = 1, z_procs
            do yi = 1, ymax+1
                do xi = 1, x_procs
                    u_rms(xi,yi,zi) = (u1_procs_fluctuation(xi,yi,zi)**2.0d0 &
                                    + u2_procs_fluctuation(xi,yi,zi)**2.0d0 &
                                    + u3_procs_fluctuation(xi,yi,zi)**2.0d0)
                enddo
            enddo
        enddo
        u_rms_ave = 0.0d0
        do zi = 1, z_procs
            do yi = 1, ymax+1
                do xi = 1, x_procs
                    u_rms_ave = u_rms_ave + u_rms(xi,yi,zi)
                enddo
            enddo
        enddo
    end subroutine u_rms_cal

    subroutine u_grad_cal(u_grad,u_grad_ave,grad_u_procs)
        real(8),intent(out):: u_grad(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(inout):: u_grad_ave
        real(8),intent(in):: grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)

        do zi = 1, z_procs
            do yi = 1, ymax+1
                do xi = 1, x_procs
                    u_grad(xi,yi,zi) = (grad_u_procs(1,1,xi,yi,zi)**2.0d0 &
                                    + grad_u_procs(2,2,xi,yi,zi)**2.0d0 &
                                    + grad_u_procs(3,3,xi,yi,zi)**2.0d0)
                enddo
            enddo
        enddo
        u_grad_ave = 0.0d0
        do zi = 1, z_procs
            do yi = 1, ymax+1
                do xi = 1, x_procs
                    u_grad_ave = u_grad_ave + u_grad(xi,yi,zi)
                enddo
            enddo
        enddo
    end subroutine u_grad_cal

    subroutine taylor_re_cal(u_rms_ave,u_rms_ave_sum,u_grad_ave,u_grad_ave_sum,taylor_re,taylor_length,u_rms,u_grad)
        real(8),intent(inout):: u_rms_ave
        real(8),intent(inout):: u_grad_ave

        real(8),intent(in):: u_rms(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: u_grad(1:x_procs,1:ymax+1,1:z_procs)

        real(8),intent(inout):: u_rms_ave_sum
        real(8),intent(inout):: u_grad_ave_sum
        real(8),intent(out):: taylor_re
        real(8),intent(out):: taylor_length

        u_rms_ave_sum = 0.0d0
        u_grad_ave_sum = 0.0d0

        call MPI_Reduce(u_rms_ave, u_rms_ave_sum, 1, MPI_REAL8, MPI_SUM, 0, new_comm_world, ierr)
        call MPI_Reduce(u_grad_ave, u_grad_ave_sum, 1, MPI_REAL8, MPI_SUM, 0, new_comm_world, ierr)
        
        if(new_rank == 0) then
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
        real(8),intent(in):: strain_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(out):: energy_dissipation, energy_dissipation_sum, Kolmogorov_scale

        energy_dissipation = 0.0d0
        energy_dissipation_sum = 0.0d0
        Kolmogorov_scale = 0.0d0
        do zi = 1, z_procs
            do yi = 1, ymax+1
                do xi = 1, x_procs
                    do beta = 1, 3
                        do alpha = 1, 3
                            energy_dissipation_sum = energy_dissipation_sum + strain_procs(alpha,beta,xi,yi,zi)**2
                        enddo
                    enddo
                enddo
            enddo
        enddo

        call MPI_Reduce(energy_dissipation_sum, energy_dissipation, 1, MPI_REAL8, MPI_SUM, 0, new_comm_world, ierr)
        
        if(new_rank == 0) then
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
        real(8),intent(in):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(out):: u1_procs_ave(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(out):: u2_procs_ave(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(out):: u3_procs_ave(1:x_procs,1:ymax+1,1:z_procs)

        do zi = 1, z_procs
            do yi = 1, ymax+1
                do xi = 1, x_procs
                    u1_procs_ave(xi,yi,zi) = u1_procs_ave(xi,yi,zi) + u1_procs(xi,yi,zi)
                    u2_procs_ave(xi,yi,zi) = u2_procs_ave(xi,yi,zi) + u2_procs(xi,yi,zi)
                    u3_procs_ave(xi,yi,zi) = u3_procs_ave(xi,yi,zi) + u3_procs(xi,yi,zi)
                enddo
            enddo
        enddo
    end subroutine u_sum_step_cal

    subroutine vorticity_cal(grad_u_procs,vorticity_procs)
        real(8),intent(in):: grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(out):: vorticity_procs(1:3,1:x_procs,1:ymax+1,1:z_procs)
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    vorticity_procs(1,xi,yi,zi) = grad_u_procs(3,2,xi,yi,zi) - grad_u_procs(2,3,xi,yi,zi)
                    vorticity_procs(2,xi,yi,zi) = grad_u_procs(1,3,xi,yi,zi) - grad_u_procs(3,1,xi,yi,zi)
                    vorticity_procs(3,xi,yi,zi) = grad_u_procs(2,1,xi,yi,zi) - grad_u_procs(1,2,xi,yi,zi)
                enddo
            enddo
        enddo
    end subroutine vorticity_cal

    subroutine enstrophy_cal(vorticity_procs,enstrophy_procs)
        real(8),intent(in):: vorticity_procs(1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(out):: enstrophy_procs(1:x_procs,1:ymax+1,1:z_procs)
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    enstrophy_procs(xi,yi,zi) = 0.5d0 * (vorticity_procs(1,xi,yi,zi)**2 + &
                                                        vorticity_procs(2,xi,yi,zi)**2 + &
                                                        vorticity_procs(3,xi,yi,zi)**2)
                enddo
            enddo
        enddo
    end subroutine enstrophy_cal

    subroutine Q_cal(Q_procs,strain_procs,vtensor_procs)
        real(8),intent(inout):: Q_procs(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: strain_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: vtensor_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8) ftmp, gtmp

        do zi=1,z_procs
            do yi=1,ymax+1
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
        real(8),intent(inout):: Q_procs(1:x_procs,1:ymax+1,1:z_procs)

        write(filename,*) group_new !i->filename 変換
        write(filename3,*) new_rank
        write(filename2,*) n
        filename=datadir_output//trim(adjustl(filename))//'_'//trim(adjustl(filename3))//'_'//trim(adjustl(filename2))//'.bin' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
        print *, filename !表示してみる
        open(101,file=filename, form='unformatted',status='replace') 
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    write(101) Q_procs(xi,yi,zi)
                enddo
            enddo
        enddo
        close(101)
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
!==================================MPI並列開始===============================================
    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, comm_procs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, comm_rank, ierr)
    call cpu_time(time1)
!================================ディレクトリの作成============================================
    call mk_dirs(datadir_input)
    call mk_dirs(datadir_output)
    call par(cx,cy,cz,cr)
    call ini(p_procs,u1_procs,u2_procs,u3_procs)
    call variable(grad_u_procs,vorticity_procs,enstrophy_procs,u1_procs_ave,u2_procs_ave,u3_procs_ave,u1_procs_fluctuation,u2_procs_fluctuation,u3_procs_fluctuation,strain_procs,u_rms,u_grad,vtensor_procs,Q_procs)
    allocate(kinetic_procs(1:x_procs,1:ymax+1,1:z_procs))
!==============時間平均の計算==================
    step_num = 0
    DO n = step_start, step_end, 5000
        step_num = step_num + 1
        call input(u1_procs,u2_procs,u3_procs,p_procs)
        call u_sum_step_cal(u1_procs,u2_procs,u3_procs,u1_procs_ave,u2_procs_ave,u3_procs_ave)
    ENDDO
    do zi = 1, z_procs
        do yi = 1, ymax+1
            do xi = 1, x_procs
                u1_procs_ave(xi,yi,zi) = u1_procs_ave(xi,yi,zi) / (dble(step_num))
                u2_procs_ave(xi,yi,zi) = u2_procs_ave(xi,yi,zi) / (dble(step_num))
                u3_procs_ave(xi,yi,zi) = u3_procs_ave(xi,yi,zi) / (dble(step_num))
            enddo
        enddo
    enddo
!================================物理量計算================================
    DO n = step_start, step_end, 5000
        !入力ファイル読み込み
        call input(u1_procs,u2_procs,u3_procs,p_procs)
        !変動速度の計算
        u1_procs_fluctuation(:,:,:) = 0.0d0
        u2_procs_fluctuation(:,:,:) = 0.0d0
        u3_procs_fluctuation(:,:,:) = 0.0d0
        do zi = 1, z_procs
            do yi = 1, ymax+1
                do xi = 1, x_procs
                    u1_procs_fluctuation(xi,yi,zi) = u1_procs(xi,yi,zi) - u1_procs_ave(xi,yi,zi)
                    u2_procs_fluctuation(xi,yi,zi) = u2_procs(xi,yi,zi) - u2_procs_ave(xi,yi,zi)
                    u3_procs_fluctuation(xi,yi,zi) = u3_procs(xi,yi,zi) - u3_procs_ave(xi,yi,zi)
                enddo
            enddo
        enddo
        !コルモゴロフスケールの計算
        call glue(u1_procs_fluctuation)
        call glue(u2_procs_fluctuation)
        call glue(u3_procs_fluctuation)
        call grad_u_cal(grad_u_procs,u1_procs_fluctuation,u2_procs_fluctuation,u3_procs_fluctuation)
        call strain_cal(grad_u_procs,strain_procs)
        call energy_dissipation_cal(strain_procs,energy_dissipation,energy_dissipation_sum,Kolmogorov_scale)

        !テイラー長レイノルズ数の計算
        call u_rms_cal(u_rms,u_rms_ave,u1_procs_fluctuation,u2_procs_fluctuation,u3_procs_fluctuation)
        call u_grad_cal(u_grad,u_grad_ave,grad_u_procs)
        call taylor_re_cal(u_rms_ave,u_rms_ave_sum,u_grad_ave,u_grad_ave_sum,taylor_re,taylor_length,u_rms,u_grad)


        ! do zi=1,z_procs
        !     do yi=1,ymax+1
        !         do xi=1,x_procs
        !             kinetic_procs(xi,yi,zi) = 0.5d0 * (u1_procs(xi,yi,zi)*u1_procs(xi,yi,zi) &
        !                                     + u2_procs(xi,yi,zi)*u2_procs(xi,yi,zi) &
        !                                     + u3_procs(xi,yi,zi)*u3_procs(xi,yi,zi))
        !         enddo
        !     enddo
        ! enddo
        ! kinetic_sum = 0.0d0
        ! do zi=1,z_procs
        !     do yi=1,ymax+1
        !         do xi=1,x_procs
        !             kinetic_sum = kinetic_sum + kinetic_procs(xi,yi,zi)
        !         enddo
        !     enddo
        ! enddo
        ! call MPI_Reduce(kinetic_sum, kinetic_ave, 1, MPI_REAL8, MPI_SUM, 0, new_comm_world, ierr)
        ! if(new_rank == 0) then
        !     kinetic_ave = kinetic_ave / (dble(xmax+1)*dble(ymax+1)*dble(zmax+1))
        !     if(n == 4000) then
        !         open(11,file="./kinetic.d")
        !         write(11,*) dble(n), dble(n)*umax/D_vortex, kinetic_ave
        !         close(11)
        !     else
        !         open(11,file="./kinetic.d",action="write",position="append")
        !         write(11,*) dble(n), dble(n)*umax/D_vortex, kinetic_ave
        !         close(11)
        !     endif
        ! endif
    ENDDO

    !=======エンストロフィー＝＝＝＝＝＝＝＝＝＝＝＝＝
    ! DO n = 600000, 1300000, 4000
    !     ! n = 1300000
    !     call input(u1_procs,u2_procs,u3_procs,p_procs)
    !     call glue(u1_procs)
    !     call glue(u2_procs)
    !     call glue(u3_procs)
    !     call grad_u_cal(grad_u_procs, u1_procs, u2_procs, u3_procs)
    !     call vorticity_cal(grad_u_procs,vorticity_procs)
    !     call enstrophy_cal(vorticity_procs,enstrophy_procs)
    !     enstrophy_sum = 0.0d0
    !     do zi=1,z_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 enstrophy_sum = enstrophy_sum + enstrophy_procs(xi,yi,zi)
    !             enddo
    !         enddo
    !     enddo
    !     call MPI_Reduce(enstrophy_sum, enstrophy_ave, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    !     if(comm_rank == 0) then
    !         enstrophy_ave = enstrophy_ave / (dble(xmax+1)*dble(ymax+1)*dble(zmax+1))
    !         open(100,file = "./enstrophy_ave_3.d",action="write",position="append")
    !         write(100,*) enstrophy_ave
    !         close(100)
    !     endif
    !     call output(enstrophy_procs)
    ! ENDDO

    !Q
    ! DO n = 5000, 100000, 1000
    !     call input(u1_procs,u2_procs,u3_procs,p_procs)
    !     call glue(u1_procs)
    !     call glue(u2_procs)
    !     call glue(u3_procs)
    !     call grad_u_cal(grad_u_procs,u1_procs,u2_procs,u3_procs)
    !     call strain_cal(grad_u_procs,strain_procs)
    !     call vtensor_cal(grad_u_procs,vtensor_procs)
    !     call Q_cal(Q_procs,strain_procs,vtensor_procs)
    !     call output(Q_procs)
    ! ENDDO
!================MPI並列終わり=======================================
    call MPI_Finalize(ierr)
end program main

!===============sub program=========================================

! DO n = step_start, step_end, 4000
    !=======エンストロフィー＝＝＝＝＝＝＝＝＝＝＝＝＝
    ! call input(u1_procs,u2_procs,u3_procs,p_procs)
    ! call glue(u1_procs)
    ! call glue(u2_procs)
    ! call glue(u3_procs)
    ! call grad_u_cal(grad_u_procs, u1_procs, u2_procs, u3_procs)
    ! call vorticity_cal(grad_u_procs,vorticity_procs)
    ! call enstrophy_cal(vorticity_procs,enstrophy_procs)
    ! enstrophy_sum = 0.0d0
    ! do zi=1,z_procs
    !     do yi=1,ymax+1
    !         do xi=1,x_procs
    !             enstrophy_sum = enstrophy_sum + enstrophy_procs(xi,yi,zi)
    !         enddo
    !     enddo
    ! enddo
    ! call MPI_Reduce(enstrophy_sum, enstrophy_ave, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    ! if(comm_rank == 0) then
    !     enstrophy_ave = enstrophy_ave / (dble(xmax+1)*dble(ymax+1)*dble(zmax+1))
    !     open(100,file = "./enstrophy_ave.d")
    !     write(100,*) enstrophy_ave
    !     close(100)
    ! endif
    ! call output(enstrophy_procs)
! ENDDO