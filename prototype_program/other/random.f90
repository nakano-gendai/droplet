program main
    implicit none
    include 'mpif.h'
    integer, parameter :: NX = 4, NY = NX, NZ = NX
    real(8), allocatable :: U_procs(:, :, :), U2_procs(:,:,:)
    real(8), allocatable :: x_m(:, :, :), y_m(:,:,:)
    character(2) chmyrank
    integer N_procs
    integer ierr, procs, myrank
    integer seed, seed2
    integer, allocatable :: seeds(:), seeds2(:)
    integer xi, yi, zi, k
    real(8), parameter :: pi = 3.14d0
    integer seedsize
    character(8) date
    character(10) time
    character(5) zone
    integer values(8)

    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, procs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, myrank, ierr)

    N_procs = NZ / procs
    allocate(U_procs(0:NX+1, 0:NY+1, 0:N_procs+1))
    allocate(U2_procs(0:NX+1, 0:NY+1, 0:N_procs+1))
    allocate(x_m(0:NX+1, 0:NY+1, 0:N_procs+1))
    allocate(y_m(0:NX+1, 0:NY+1, 0:N_procs+1))

    call random_seed(size=seedsize)
    allocate(seeds(seedsize))
    call random_seed(put=seeds)
    call date_and_time(date, time, zone, values)
    seeds(seedsize) = values(8) + myrank
    ! seeds(1) = seed
    ! call random_seed(put=seeds)  ! 入力するシード値は配列
    
    call random_number(U_procs)

    seedsize = 0
    call random_seed(size=seedsize)
    allocate(seeds2(seedsize))
    call random_seed(put=seeds2)
    call date_and_time(date, time, zone, values)
    seeds2(seedsize) = values(8) + myrank
    ! seeds(1) = seed
    ! call random_seed(put=seeds2)  ! 入力するシード値は配列
    
    call random_number(U2_procs)

    ! call system_clock(seed2)  ! 現在時刻をシード値として使用する
    ! seed2 = seed2 + myrank  ! プロセスごとにシード値を変える
    ! seeds2(1) = seed2
    ! call random_seed(put=seeds2)  ! 入力するシード値は配列
    
    ! call random_number(U2_procs)

    do zi=0,N_procs+1
        do yi=0,NY+1
            do xi=0,NX+1
                x_m(xi,yi,zi) = sqrt(-2.0d0*log(U_procs(xi,yi,zi)))*cos(2.0d0*pi*U2_procs(xi,yi,zi))
                y_m(xi,yi,zi) = sqrt(-2.0d0*log(U_procs(xi,yi,zi)))*sin(2.0d0*pi*U2_procs(xi,yi,zi))
            enddo
        enddo
    enddo


    write(chmyrank, '(I2.2)') myrank
    open(10, file = './'//chmyrank//'.d')
    write(10,*) seeds(:)
    write(10,*) seeds2(:)
    do zi=0,N_procs+1
        do yi=0,NY+1
            do xi=0,NX+1
                write(10,*) U_procs(xi,yi,zi), U2_procs(xi,yi,zi)
            enddo
        enddo
    enddo
    close(10)


    call MPI_Finalize(ierr)
    
end program main