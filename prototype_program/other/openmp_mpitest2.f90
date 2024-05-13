module globals
    include "mpif.h"
    !変数の定義
    integer,parameter :: xmax = 99
    integer,parameter :: ymax = 99
    integer,parameter :: zmax = 99
    integer xi, yi, zi
    real(8) time1, time2

    !MPI用変数
    integer ierr, comm_procs, comm_rank
    integer x_procs, z_procs
    integer Nx, Nz
    integer key_new, group_new, key_x, group_x, key_z, group_z
    integer newx_comm_world, newx_procs, newx_rank
    integer newz_comm_world, newz_procs, newz_rank
    integer rectangle_type, xtype
    integer next_rank_x, former_rank_x, next_rank_z, former_rank_z
    integer req1s, req1r, req2s, req2r
    integer, dimension(MPI_STATUS_SIZE) :: sta1s, sta1r, sta2s, sta2r

contains
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

end module globals

program main
use globals
    implicit none
    real(8), allocatable :: f(:, :, :)

    !MPI計算の開始
    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, comm_procs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, comm_rank, ierr)
    !並列数の定義
    Nx = 5
    Nz = comm_procs / Nx
    x_procs = (xmax+1) / Nx
    z_procs = (zmax+1) / Nz

    !z方向にコミュニケータを分割する
    key_z = comm_rank
    group_z = comm_rank / Nx
    call MPI_Comm_Split(MPI_COMM_WORLD,group_z,key_z,newz_comm_world,ierr)
    call MPI_Comm_Size(newz_comm_world,newz_procs,ierr)
    call MPI_Comm_Rank(newz_comm_world,newz_rank,ierr)

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

    !配列の定義
    allocate(f(0:x_procs+1,0:ymax+2,0:z_procs+1))
    !初期化
    f(:,:,:) = 0.0d0

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    time1 = MPI_Wtime()

    !$omp parallel
    do zi = 1, z_procs
        do yi = 1, ymax
            do xi = 1, x_procs
                f(xi, yi, zi) = 100.0d0 * dble(xi) + 10.0d0 * dble(yi) + dble(zi)
            enddo
        enddo
    enddo
    !$omp end parallel

    !MPI通信（のりしろ境界条件）
    call glue(f)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    time2 = MPI_Wtime()
    if(comm_rank == 0) then
        open(10,file="./time2.d")
        write(10,*) comm_rank, time2 - time1
        close(10)
    endif

    call MPI_Finalize(ierr)
end program main