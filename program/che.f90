program main
implicit none
    real(8) che
    real(8) phi
    real(8),parameter:: a = 9.0d0/49.0d0
    real(8),parameter:: b = 2.0d0/21.0d0
    real(8),parameter:: T = 0.55d0
    integer i
    real(8),parameter:: bin = 0.0001d0 

    open(10,file="che.d")
    do i = 0, 10000000
        write(*,*) i
        phi = dble(i) * bin
        che = T*log(phi/(1.0d0-b*phi)) + T/(1.0d0-b*phi) - 2.0d0*a*phi
        write(10, *) phi, che
    enddo
    close(10)
end program main