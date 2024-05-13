program main
    real(8),parameter:: ds = 1.0d0 !格子間隔（lattice unit）
    integer,parameter:: xmax = 255 !ｘ方向格子数（０から数える）
    integer,parameter:: ymax = 255 !ｙ方向格子数（０から数える）
    integer,parameter:: zmax = 255 !ｚ方向格子数（０から数える）
    real(8),parameter:: phi1 = 2.211d0
    real(8),parameter:: phi2 = 4.895d0

    real(8) phi(-1:xmax+1,-1:ymax+1,-1:zmax+1)
    integer label(0:xmax,0:ymax,0:zmax)
    ! integer label2(0:xmax,0:ymax,0:zmax)
    integer klass(1:100000)
    ! integer klass2(100)
    integer tmp(1:7)
    integer labeltmp
    integer xi, yi, zi, i, k, a, b, x_up, y_up, z_up, x_down, y_down, z_down
    real(8) Vlabel(1:100000), radius(1:100000)

    character(*),parameter :: datadir = "/data/sht/nakanog/taylor_re12000_drop8_2/collect/"
    integer step, step_num, dropnum
    character(8) file_num
    real(8) dummy1, dummy2, dummy3, dummy4

    !初期液滴中心
    real(8), parameter:: xc = 0.5d0*dble(xmax)	
    real(8), parameter:: yc = 0.5d0*dble(ymax)
    real(8), parameter:: zc = 0.5d0*dble(zmax)

step_num = 19
! DO step = 5000, 920000, 1000
    step = 900000
    step_num = step_num + 1
    if((step > 99) .and. (step < 1000)) then
        write(file_num, "(i3)") step
    elseif((step > 999) .and. (step < 10000)) then
        write(file_num,"(i4)") step
    elseif((step > 9999) .and. (step < 100000)) then
        write(file_num,"(i5)") step
    elseif((step > 99999) .and. (step < 1000000)) then
        write(file_num,"(i6)") step
    elseif((step > 999999) .and. (step < 10000000)) then
        write(file_num,"(i7)") step
    endif
    open(step_num, file=datadir//trim(file_num)//".bin", form="unformatted")
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                read(step_num)  phi(xi,yi,zi)
            enddo
        enddo
    enddo
    close(step_num)

    ! do zi=0,zmax
    !     do yi=0,ymax
    !         do xi=0,xmax
    !             phi(xi,yi,zi) = phi1
                
    !             if((dble(xi)*ds-0.5d0*xc)**2+(dble(yi)*ds-yc)**2+(dble(zi)*ds-zc)**2 <= (70.d0)**2) then
    !                 phi(xi,yi,zi) = phi2
    !             endif
    !             if((dble(xi)*ds-1.5d0*xc)**2+(dble(yi)*ds-yc)**2+(dble(zi)*ds-zc)**2 <= (70.d0)**2) then
    !                 phi(xi,yi,zi) = phi2
    !             endif
    !             ! if((dble(xi)*ds-xc)**2+(dble(yi)*ds-yc)**2+(dble(zi)*ds-zc)**2 <= (0.5d0*20.d0)**2) then
    !             !     phi(xi,yi,zi) = phi2
    !             ! endif
    !         enddo
    !     enddo
    ! enddo

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
    klass(:) = 0
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                if(phi(xi,yi,zi) >= (phi1+phi2)/2.0d0) then
                    x_up = xi + 1
                    x_down = xi - 1
                    y_up = yi + 1
                    y_down = yi - 1
                    z_up = zi + 1
                    z_down = zi - 1
                    if(xi == 0) then
                        x_down = xmax
                    else if(xi == xmax) then
                        x_up = 0
                    endif

                    if(yi == 0) then
                        y_down = ymax
                    else if(yi == ymax) then
                        y_up = 0
                    endif

                    if(zi == 0) then
                        z_down = zmax
                    else if(zi == zmax) then
                        z_up = 0
                    endif

                    tmp(:) = 0
                    tmp(1) = label(x_down,yi,zi)
                    tmp(2) = label(xi,y_down,zi)
                    tmp(3) = label(xi,yi,z_down)
                    tmp(4) = label(xi,y_up,z_down)
                    tmp(5) = label(x_up,yi,z_down)
                    tmp(6) = label(xi,y_down,z_down)
                    tmp(7) = label(x_down,yi,z_down)
                    labeltmp = xmax*ymax*zmax
                    a = 0
                    b = 0
                    do i = 1, 7
                        if((tmp(i) > 0) .and. (labeltmp > tmp(i))) then
                            labeltmp = tmp(i)
                        endif
                    enddo

                    if(tmp(1)+tmp(2)+tmp(3)+tmp(4)+tmp(5)+tmp(6)+tmp(7) == 0) then
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
                        endif

                        if(label(xi,yi-1,zi) > 0) then
                            klass(label(xi,yi-1,zi)) = label(xi,yi,zi)

                            a = label(xi,yi-1,zi)
                        2   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 2
                            else
                                label(xi,yi-1,zi) = a
                            endif
                        endif

                        if(label(xi,yi,zi-1) > 0) then
                            klass(label(xi,yi,zi-1)) = label(xi,yi,zi)

                            a = label(xi,yi,zi-1)
                        3   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 3
                            else
                                label(xi,yi,zi-1) = a
                            endif
                        endif

                        if(label(xi,yi+1,zi-1) > 0) then
                            klass(label(xi,yi+1,zi-1)) = label(xi,yi,zi)

                            a = label(xi,yi+1,zi-1)
                        4   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 4
                            else
                                label(xi,yi+1,zi-1) = a
                            endif
                        endif
                        if(label(xi+1,yi,zi-1) > 0) then
                            klass(label(xi+1,yi,zi-1)) = label(xi,yi,zi)

                            a = label(xi+1,yi,zi-1)
                        5   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 5
                            else
                                label(xi+1,yi,zi-1) = a
                            endif
                        endif
                        if(label(xi,yi-1,zi-1) > 0) then
                            klass(label(xi,yi-1,zi-1)) = label(xi,yi,zi)

                            a = label(xi,yi-1,zi-1)
                        6   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 6
                            else
                                label(xi,yi-1,zi-1) = a
                            endif
                        endif
                        if(label(xi-1,yi,zi-1) > 0) then
                            klass(label(xi-1,yi,zi-1)) = label(xi,yi,zi)

                            a = label(xi-1,yi,zi-1)
                        7   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 7
                            else
                                label(xi-1,yi,zi-1) = a
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
                8   b = a
                    a = klass(a)
                    if(a /= b) then
                        go to 8
                    else
                        label(xi,yi,zi) = a
                    endif
                endif
            enddo
        enddo
    enddo


    Vlabel(:) = 0.0d0
    DO i = 1, 10000
        do zi = 0, zmax
            do yi = 0, ymax
                do xi = 0, xmax
                    if(label(xi,yi,zi) == i) then
                        Vlabel(i) = Vlabel(i) + 1.0d0
                    endif
                enddo
            enddo
        enddo
        radius(i) = (3.0d0*Vlabel(i)/(4.0d0*3.14d0))**(1.0d0/3.0d0)
    ENDDO

    dropnum = 0
    do i = 1, 10000
        if(Vlabel(i) > 0.0d0) then
            dropnum = dropnum + 1
        endif
    enddo

    ! if(step == 5000) then
    !     open(10,file="./drop_num.d")
    !     write(10,*) dble(step), dble(dropnum)
    !     close(10)
    ! else
    !     open(10,file="./drop_num.d",action="write",position="append")
    !     write(10,*) dble(step), dble(dropnum)
    !     close(10)
    ! endif

    open(11,file="./radius.d")
    do i = 1, 100000
        if(radius(i) > 0.0d0) then
            write(11,*) radius(i)
        endif
    enddo
    close(11)
! ENDDO
end program main