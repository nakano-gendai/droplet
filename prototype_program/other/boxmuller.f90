program main
implicit none
    integer seed1, seed2
    integer seeds1(1), seeds2(1)
    real(8) eps1(10000)
    real(8) eps2(10000)
    real(8) x_m(10000), y_m(10000)
    integer i
    real(8), parameter :: pi = 3.14d0

    integer kazu(6)
    integer other
    real(8),parameter :: bin = 1.0d0
    integer kk
    real(8) saidai, saisyo

    call system_clock(seed1)
    seed2 = seed1 + 1
    seeds1(1) = seed1
    seeds2(1) = seed2

    call random_seed(put=seeds1)
    call random_number(eps1)
    call random_seed(put=seeds2)
    call random_number(eps2)

    do i = 1, 10000
        x_m(i) = sqrt(-2.0d0*log(eps1(i)))*cos(2.0d0*pi*eps2(i))
        y_m(i) = sqrt(-2.0d0*log(eps1(i)))*sin(2.0d0*pi*eps2(i))
    enddo

    kazu(:) = 0
    other = 0
    do i = 1, 10000
        if((x_m(i) >= -3.0d0) .and. (x_m(i) < -2.0d0)) then
            kazu(1) = kazu(1) + 1
        elseif((x_m(i) >= -2.0d0) .and. (x_m(i) < -1.0d0)) then
            kazu(2) = kazu(2) + 1
        elseif((x_m(i) >= -1.0d0) .and. (x_m(i) < 0.0d0)) then
            kazu(3) = kazu(3) + 1
        elseif((x_m(i) >= 0.0d0) .and. (x_m(i) < 1.0d0)) then
            kazu(4) = kazu(4) + 1
        elseif((x_m(i) >= 1.0d0) .and. (x_m(i) < 2.0d0)) then
            kazu(5) = kazu(5) + 1
        elseif((x_m(i) >= 2.0d0) .and. (x_m(i) < 3.0d0)) then
            kazu(6) = kazu(6) + 1
        elseif((x_m(i) > 3.0d0) .or. (x_m(i) < -3.0d0)) then
            other = other + 1
        endif
    enddo

    open(10,file="boxmuller.d")
    do i = 1, 6
        ! write(10,*) i, x_m(i), y_m(i)
        write(10,*) i, kazu(i)
    end do
    close(10)
    saidai = 0.0d0
    saisyo = 100000.0d0
    do i = 1, 10000
        if(x_m(i) > saidai) then
            saidai = x_m(i)
        elseif(x_m(i) < saisyo) then
            saisyo = x_m(i)
        endif
    enddo
    write(*,*) saidai, saisyo

    open(11,file="seiki.d")
    write(11,*) "other      -1 < x_m < 1      -2 < x_m < 2      -3 < x_m < 3"
    write(11,*) other, (real(kazu(3)+kazu(4))) / 10000.0d0, (real(kazu(2)+kazu(3)+kazu(4)+kazu(5))) / 10000.0d0, (real(kazu(1)+kazu(2)+kazu(3)+kazu(4)+kazu(5)+kazu(6))) / 10000.0d0
    close(11)
end program main