program main
    use omp_lib
    implicit none

    integer,parameter :: xmax = 100
    integer,parameter :: ymax = 100
    integer,parameter :: zmax = 100
    real(8) f(xmax, ymax, zmax)
    integer xi, yi, zi
    real(8) time1, time2
    integer num

    f(:,:,:) = 0.0d0

    call cpu_time(time1)

    !$omp parallel
    do zi = 1, zmax
        do yi = 1, ymax
            do xi = 1, xmax
                f(xi, yi, zi) = 100.0d0 * dble(xi) + 10.0d0 * dble(yi) + dble(zi)
            enddo
        enddo
    enddo
    !$omp end parallel

    call cpu_time(time2)
    write(*,*) time2 - time1

end program main