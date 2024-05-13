program main
    integer,parameter:: xmax = 128 !ｘ方向格子数（０から数える）
    integer,parameter:: ymax = 64 !ｙ方向格子数（０から数える）
    integer,parameter:: zmax = 128 !ｚ方向格子数（０から数える）
    real(8),parameter:: phi1 = 2.211d0
    real(8),parameter:: phi2 = 4.895d0

    real(8) phi(-1:xmax+1,-1:ymax+1,-1:zmax+1)
    integer label(0:xmax,0:ymax,0:zmax)
    integer klass(1000)
    integer tmp(5)
    integer Itmp(5)
    integer labeltmp
    integer xi, yi, zi, i, j, k, a, b
    integer Vlabel, Vphi

    open(10,file="/data/n/n517/par/fg/collect/0_5000.bin",form="unformatted")
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                read(10)  phi(xi,yi,zi)
            enddo
        enddo
    enddo
    close(10)
    !phiの周期境界
    do zi=0,zmax
        do yi=-1,ymax+1
            do xi=0,xmax
                phi(xi,-1,zi) = phi(xi,ymax,zi)
                phi(xi,ymax+1,zi) = phi(xi,0,zi)
                phi(-1,yi,zi) = phi(xmax,yi,zi)
                phi(xmax+1,yi,zi) = phi(0,yi,zi)
            enddo
        enddo
    enddo
    do yi=-1,ymax+1
        do xi=-1,xmax+1
            phi(xi,yi,-1) = phi(xi,yi,zmax)
            phi(xi,yi,zmax+1) = phi(xi,yi,0)
        enddo
    enddo

    k = 1
    label(:,:,:) = 0
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                if(phi(xi,yi,zi) >= (phi1+phi2)/2.0d0) then
                    tmp(1) = label(xi-1,yi,zi)
                    tmp(2) = label(xi,yi-1,zi)
                    tmp(3) = label(xi,yi,zi-1)
                    tmp(4) = label(xi,yi+1,zi-1)
                    tmp(5) = label(xi,yi+2,zi-1)
                    labeltmp = xmax*ymax*zmax
                    Itmp(:) = 0
                    a = 0
                    b = 0
                    do i = 1, 5
                        if((tmp(i) > 0) .and. (labeltmp > tmp(i))) then
                            labeltmp = tmp(i)
                            Itmp(i) = i
                        endif
                    enddo

                    if(Itmp(1)+Itmp(2)+Itmp(3)+Itmp(4)+Itmp(5) == 0) then
                        label(xi,yi,zi) = k
                        klass(k) = k
                        k = k + 1
                    else
                        label(xi,yi,zi) = labeltmp

                        if(label(xi-1,yi,zi) > 0) then
                            klass(label(xi-1,yi,zi)) = label(xi,yi,zi)

                            a = label(xi-1,yi,zi)
                        1   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 1
                            else
                                label(xi-1,yi,zi) = a
                            endif

                        else if(label(xi,yi-1,zi) > 0) then
                            klass(label(xi,yi-1,zi)) = label(xi,yi,zi)

                            a = label(xi,yi-1,zi)
                        2   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 2
                            else
                                label(xi,yi-1,zi) = a
                            endif

                        else if(label(xi,yi,zi-1) > 0) then
                            klass(label(xi,yi,zi-1)) = label(xi,yi,zi)

                            a = label(xi,yi,zi-1)
                        3   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 3
                            else
                                label(xi,yi,zi-1) = a
                            endif

                        else if(label(xi,yi+1,zi-1) > 0) then
                            klass(label(xi,yi+1,zi-1)) = label(xi,yi,zi)

                            a = label(xi,yi+1,zi-1)
                        4   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 4
                            else
                                label(xi,yi+1,zi-1) = a
                            endif

                        else if(label(xi,yi+2,zi-1) > 0) then
                            klass(label(xi,yi+2,zi-1)) = label(xi,yi,zi)

                            a = label(xi,yi+2,zi-1)
                        5   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 5
                            else
                                label(xi,yi+2,zi-1) = a
                            endif
                        endif
                    endif
                endif
            enddo
        enddo
    enddo

    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                if(phi(xi,yi,zi) >= (phi1+phi2)/2.0d0) then
                    a = label(xi,yi,zi)
                7   b = a
                    a = klass(a)
                    if(a /= b) then
                        go to 7
                    else
                        label(xi,yi,zi) = a
                    endif
                endif
            enddo
        enddo
    enddo

    Vlabel = 0
    Vphi = 0
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                if(phi(xi,yi,zi) >= (phi1+phi2) / 2.0d0) then
                    Vphi = Vphi + 1
                endif
                if(label(xi,yi,zi) == 1) then
                    Vlabel = Vlabel +1
                endif
            enddo
        enddo
    enddo

    open(11,file="./label.d")
    zi = 52
    do yi = 0, ymax
        do xi = 0, xmax
            write(11,*) xi, yi, label(xi,yi,zi)
        enddo
        write(11,*)
    enddo
    close(11)

    open(12,file="./vol.d")
    write(12,*) Vlabel, Vphi
    close(12)
end program main