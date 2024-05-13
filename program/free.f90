program main
implicit none
real(8) f_max, f_min
real(8) f(1:16, 1:16, 1:64)
integer xi, yi, zi, i, j, k

open(10,file="v_15.d")
do zi=1,64
    do yi=1,16
        do xi=1,16
            read(10,*) i, j, k, f(xi,yi,zi)
        enddo
    enddo
enddo
close(10)

f_max = 0.0d0
f_min = 100.0d0
do zi=1,64
    do yi=1,16
        do xi=1,16
            if(f(xi,yi,zi) > f_max) then
                f_max = f(xi,yi,zi)
            endif
            if(f(xi,yi,zi) < f_min) then
                f_min = f(xi,yi,zi)
            endif
        enddo
    enddo
enddo
write(*,*) "max = ", f_max
write(*,*) "min = ", f_min

end program main