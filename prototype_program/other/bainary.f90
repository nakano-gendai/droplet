program main
implicit none
real(8) buf(0:5,0:5)
integer i,j
! real(8) buf(5) 

do j=0,5
    do i=0,5
        buf(i,j) = 10.0d0*dble(i)+dble(j)
    enddo
enddo
open(10,file="buf.bin",form='UNFORMATTED')
do j=0,5
    do i=0,5
        ! write(10) buf(i,j)
        read(10) buf(i,j)
    enddo
enddo
close(10)
do j=0,5
    do i=0,5
        write(*,"(8es24.16)") dble(i),dble(j),buf(i,j)
    enddo
enddo


end program main