program main
    implicit none
    real(8) step(118), kinetic_energy(118)
    real(8) t
    real(8),parameter:: D = 80.0d0
    real(8),parameter:: umax = 0.05d0
    integer i

    open(10, file="./kinetic_re200_2.d")
    do i = 1, 118
        read(10,*) step(i), kinetic_energy(i)
    enddo
    close(10)

    open(11, file="./kinetic_re200_3.d")
    do i = 1, 118
        write(11,*) step(i)*umax/D, kinetic_energy(i)
    enddo
    close(11)

end program main