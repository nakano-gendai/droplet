program main
    implicit none
    include "mpif.h"

    integer,parameter :: xmax = 100
    integer,parameter :: ymax = 100
    integer,parameter :: zmax = 100
    real(8), allocatable :: f(:, :, :)
    integer xi, yi, zi
    real(8) time1, time2

    !MPI用変数
    integer ierr, comm_procs, comm_rank
    integer x_procs, z_procs
    integer Nx, Nz

    !MPI計算の開始
    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, comm_procs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, comm_rank, ierr)
    !並列数の定義
    Nx = 4
    Nz = comm_procs / Nx
    x_procs = xmax / Nx
    z_procs = zmax / Nz
    !配列の定義
    allocate(f(x_procs, ymax, z_procs))
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
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    time2 = MPI_Wtime()
    if(comm_rank == 0) then
        open(10,file="./time2.d")
        write(10,*) comm_rank, time2 - time1
        close(10)
    endif

    call MPI_Finalize(ierr)
end program main