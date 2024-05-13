program main
    implicit none
    integer,parameter:: xmax = 255 !ｘ方向格子数
    integer,parameter:: ymax = 255 !ｙ方向格子数
    integer,parameter:: zmax = 255 !ｚ方向格子数
    integer i, xi, yi, zi

    real(8) gosa(1:15,0:xmax,0:ymax,0:zmax), s
    real(8) max, min

    open(10,file="./gosa3.d")
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                do i = 1, 15
                    read(10,*)  s, gosa(i,xi,yi,zi)
                enddo
            enddo
        enddo
    enddo
    close(10)

    open(11,file="./max.d")
    write(11,*) "simulation start."
    close(11)
    max = 0.0d0
    min = 1.0d0
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                do i = 1, 15
                    if(gosa(i,xi,yi,zi) > max) then
                        max = gosa(i,xi,yi,zi)
                    endif
                    if(gosa(i,xi,yi,zi) < min) then
                        min = gosa(i,xi,yi,zi)
                    endif
                enddo
            enddo
        enddo
        open(11,file="./max.d",action="write",position="append")
        write(11,*) zi
        close(11)
    enddo

    open(11,file="./max.d",action="write",position="append")
    write(11,*) "max = ", max
    write(11,*) "min = ", min
    close(11)

end program main
    