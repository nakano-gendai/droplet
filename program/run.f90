program main
implicit none
    integer,parameter:: xmax = 255 !ｘ方向格子数（０から数える）
    integer,parameter:: ymax = 255 !ｙ方向格子数（０から数える）
    integer,parameter:: zmax = 255 !ｚ方向格子数（０から数える）
    real(8) eps, eps2
    integer seedsize
    integer, allocatable :: seeds(:), seeds2(:)
    integer seed, seed2
    integer i, case_num
    real(8) xc, yc, zc
    integer,parameter :: case_initial_num = 1
    integer,parameter :: case_end_num = 2

    call random_seed(size=seedsize)
    allocate(seeds(1:seedsize))
    allocate(seeds2(1:seedsize))

    DO  case_num = case_initial_num, case_end_num
        !乱数生成
        seed = 0
        call system_clock(seed)
        do i=1,seedsize
            seeds(i) = seed + i*case_num*1
        enddo
        call random_seed(put=seeds)
        call random_number(eps)
        xc = dble(xmax)*eps

        seed = 0
        call system_clock(seed)
        do i=1,seedsize
            seeds(i) = seed + i*case_num*1000
        enddo
        call random_seed(put=seeds)
        call random_number(eps)
        yc = dble(ymax)*eps

        seed = 0
        call system_clock(seed)
        do i=1,seedsize
            seeds(i) = seed + i*case_num*10000
        enddo
        call random_seed(put=seeds)
        call random_number(eps)
        zc = dble(zmax)*eps

        write(*,*) xc/dble(xmax), yc/dble(ymax), zc/dble(zmax), xc, yc, zc
    ENDDO

end program main