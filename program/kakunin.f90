program main
implicit none
integer,parameter:: xmax = 255 !ｘ方向格子数
integer,parameter:: ymax = 255 !ｙ方向格子数
integer,parameter:: zmax = 255 !ｚ方向格子数
integer i, j, k, xi, yi, zi, Nxx, Nyy, Nzz, step
real(8) err(1:15,0:xmax,0:ymax,0:zmax)

open(10,file="./err.d")
do zi=0,zmax
    do yi=0,ymax
        do xi=0,xmax
            do i=1,15
                read(10,*) Nxx, Nyy, Nzz, err(i,xi,yi,zi)
            enddo
        enddo
    enddo
enddo
close(10)

open(11,file="./err2.d")
do zi=0,zmax
    do yi=0,ymax
        do xi=0,xmax
            do i=1,15
                if(err(i,xi,yi,zi) > 0.0d0) then
                    write(11,*) xi,yi,zi,err(i,xi,yi,zi)
                endif
            enddo
        enddo
    enddo
enddo
close(11)
end program main