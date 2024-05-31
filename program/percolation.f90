program main
    real(8),parameter:: ds = 1.0d0 !格子間隔（lattice unit）
    integer,parameter:: xmax = 63 !ｘ方向格子数（０から数える）
    integer,parameter:: ymax = 63 !ｙ方向格子数（０から数える）
    integer,parameter:: zmax = 63 !ｚ方向格子数（０から数える）
    real(8),parameter:: phi1 = 2.638d-1 !連続相のオーダーパラメータ
    real(8),parameter:: phi2 = 4.031d-1 !分散相のオーダーパラメータ

    real(8) phi(-1:xmax+1,-1:ymax+1,-1:zmax+1)
    integer label(0:xmax,0:ymax,0:zmax)
    integer klass(1:1000)
    integer tmp(1:7)
    integer labeltmp
    integer xi, yi, zi, i, k, a, b, x_up, y_up, z_up, x_down, y_down, z_down
    real(8) Vlabel(1:1000), radius(1:1000)
    real(8) Vall

    character(*),parameter :: datadir = "/data/sht/nakanog/DNS_turbulence_256_IHT/case1/collect/"
    character(*),parameter :: datadir2 = "/data/sht/nakanog/DNS_turbulence_256_IHT/case1/collect/radius/"
    integer step, dropnum
    character(8) file_num
    real(8) time1, time2

    real(8),parameter:: D = 32.0d0

    !初期液滴中心
    real(8), parameter:: xc = 0.5d0*dble(xmax+1)	
    real(8), parameter:: yc = 0.5d0*dble(ymax+1)
    real(8), parameter:: zc = 0.5d0*dble(zmax+1)

    open(30,file="./process.d")
    write(30,*) "percolation start!!!!"
    close(30)
! DO step = 5000, 100000, 1000
!     call cpu_time(time1)
!     if((step > 99) .and. (step < 1000)) then
!         write(file_num, "(i3)") step
!     elseif((step > 999) .and. (step < 10000)) then
!         write(file_num,"(i4)") step
!     elseif((step > 9999) .and. (step < 100000)) then
!         write(file_num,"(i5)") step
!     elseif((step > 99999) .and. (step < 1000000)) then
!         write(file_num,"(i6)") step
!     elseif((step > 999999) .and. (step < 10000000)) then
!         write(file_num,"(i7)") step
!     elseif((step > 9999999) .and. (step < 100000000)) then
!         write(file_num,"(i8)") step
!     endif
!     open(20, file=datadir//trim(file_num)//".bin", form="unformatted")
!     do zi = 0, zmax
!         do yi = 0, ymax
!             do xi = 0, xmax
!                 read(20)  phi(xi,yi,zi)
!             enddo
!         enddo
!     enddo
!     close(20)

    !phiの周期境界
    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                if((dble(xi)*ds-xc)**2 + (dble(yi)*ds-yc)**2 + (dble(zi)*ds-0.0d0)**2 <= (0.5d0*D)**2) then
                    phi(xi,yi,zi) = phi2
                else if((dble(xi)*ds-xc)**2 + (dble(yi)*ds-yc)**2 + (dble(zi)*ds-dble(zmax))**2 <= (0.5d0*D)**2) then
                    phi(xi,yi,zi) = phi2
                else
                    phi(xi,yi,zi) = phi1
                endif
            enddo
        enddo
    enddo

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
    open(30,file="./process.d",action="write",position="append")
    write(30,*) "input is OK!!!", step
    close(30)

    k = 1
    label(:,:,:) = 0
    klass(:) = 0

    ! zi = zmax
    ! do yi = 0, ymax
    !     do xi = 0, xmax
    !         if(phi(xi,yi,zi) >= (phi1+phi2)/2.0d0) then
    !             x_up = xi + 1
    !             x_down = xi - 1
    !             y_up = yi + 1
    !             y_down = yi - 1
    !             z_up = zi + 1
    !             z_down = zi - 1
    !             if(xi == 0) then
    !                 x_down = xmax
    !             else if(xi == xmax) then
    !                 x_up = 0
    !             endif

    !             if(yi == 0) then
    !                 y_down = ymax
    !             else if(yi == ymax) then
    !                 y_up = 0
    !             endif

    !             if(zi == 0) then
    !                 z_down = zmax
    !             else if(zi == zmax) then
    !                 z_up = 0
    !             endif

    !             tmp(:) = 0
    !             tmp(1) = label(x_down,yi,zi)
    !             tmp(2) = label(xi,y_down,zi)
    !             tmp(3) = label(xi,yi,z_down)
    !             tmp(4) = label(xi,y_up,z_down)
    !             tmp(5) = label(x_up,yi,z_down)
    !             tmp(6) = label(xi,y_down,z_down)
    !             tmp(7) = label(x_down,yi,z_down)
    !             labeltmp = xmax*ymax*zmax
    !             a = 0
    !             b = 0
    !             do i = 1, 7
    !                 if((tmp(i) > 0) .and. (labeltmp > tmp(i))) then
    !                     labeltmp = tmp(i)
    !                 endif
    !             enddo

    !             if(tmp(1)+tmp(2)+tmp(3)+tmp(4)+tmp(5)+tmp(6)+tmp(7) == 0) then
    !                 label(xi,yi,zi) = k
    !                 klass(k) = k
    !                 k = k + 1
    !             else
    !                 label(xi,yi,zi) = labeltmp

    !                 if(label(x_down,yi,zi) > 0) then
    !                     klass(label(x_down,yi,zi)) = label(xi,yi,zi)

    !                     a = label(x_down,yi,zi)
    !                 11   b = a
    !                     a = klass(a)
    !                     if(a /= b)  then 
    !                         go to 11
    !                     else
    !                         label(x_down,yi,zi) = a
    !                     endif
    !                 endif

    !                 if(label(xi,y_down,zi) > 0) then
    !                     klass(label(xi,y_down,zi)) = label(xi,yi,zi)

    !                     a = label(xi,y_down,zi)
    !                 21   b = a
    !                     a = klass(a)
    !                     if(a /= b)  then 
    !                         go to 21
    !                     else
    !                         label(xi,y_down,zi) = a
    !                     endif
    !                 endif

    !                 if(label(xi,yi,z_down) > 0) then
    !                     klass(label(xi,yi,z_down)) = label(xi,yi,zi)

    !                     a = label(xi,yi,z_down)
    !                 31   b = a
    !                     a = klass(a)
    !                     if(a /= b)  then 
    !                         go to 31
    !                     else
    !                         label(xi,yi,z_down) = a
    !                     endif
    !                 endif

    !                 if(label(xi,y_up,z_down) > 0) then
    !                     klass(label(xi,y_up,z_down)) = label(xi,yi,zi)

    !                     a = label(xi,y_up,z_down)
    !                 41   b = a
    !                     a = klass(a)
    !                     if(a /= b)  then 
    !                         go to 41
    !                     else
    !                         label(xi,y_up,z_down) = a
    !                     endif
    !                 endif

    !                 if(label(x_up,yi,z_down) > 0) then
    !                     klass(label(x_up,yi,z_down)) = label(xi,yi,zi)

    !                     a = label(x_up,yi,z_down)
    !                 51   b = a
    !                     a = klass(a)
    !                     if(a /= b)  then 
    !                         go to 51
    !                     else
    !                         label(x_up,yi,z_down) = a
    !                     endif
    !                 endif

    !                 if(label(xi,y_down,z_down) > 0) then
    !                     klass(label(xi,y_down,z_down)) = label(xi,yi,zi)

    !                     a = label(xi,y_down,z_down)
    !                 61   b = a
    !                     a = klass(a)
    !                     if(a /= b)  then 
    !                         go to 61
    !                     else
    !                         label(xi,y_down,z_down) = a
    !                     endif
    !                 endif

    !                 if(label(x_down,yi,z_down) > 0) then
    !                     klass(label(x_down,yi,z_down)) = label(xi,yi,zi)

    !                     a = label(x_down,yi,z_down)
    !                 71   b = a
    !                     a = klass(a)
    !                     if(a /= b)  then 
    !                         go to 71
    !                     else
    !                         label(x_down,yi,z_down) = a
    !                     endif
    !                 endif

    !             endif
    !         endif
    !     enddo
    ! enddo


    Vall = 0.0d0
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                if(phi(xi,yi,zi) >= (phi1+phi2)/2.0d0) then
                    Vall = Vall + 1.0d0 !全液滴の体積

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

                        if(label(x_down,yi,zi) > 0) then
                            klass(label(x_down,yi,zi)) = label(xi,yi,zi)

                            a = label(x_down,yi,zi)
                        1   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 1
                            else
                                label(x_down,yi,zi) = a
                            endif
                        endif

                        if(label(xi,y_down,zi) > 0) then
                            klass(label(xi,y_down,zi)) = label(xi,yi,zi)

                            a = label(xi,y_down,zi)
                        2   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 2
                            else
                                label(xi,y_down,zi) = a
                            endif
                        endif

                        if(label(xi,yi,z_down) > 0) then
                            klass(label(xi,yi,z_down)) = label(xi,yi,zi)

                            a = label(xi,yi,z_down)
                        3   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 3
                            else
                                label(xi,yi,z_down) = a
                            endif
                        endif

                        if(label(xi,y_up,z_down) > 0) then
                            klass(label(xi,y_up,z_down)) = label(xi,yi,zi)

                            a = label(xi,y_up,z_down)
                        4   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 4
                            else
                                label(xi,y_up,z_down) = a
                            endif
                        endif

                        if(label(x_up,yi,z_down) > 0) then
                            klass(label(x_up,yi,z_down)) = label(xi,yi,zi)

                            a = label(x_up,yi,z_down)
                        5   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 5
                            else
                                label(x_up,yi,z_down) = a
                            endif
                        endif

                        if(label(xi,y_down,z_down) > 0) then
                            klass(label(xi,y_down,z_down)) = label(xi,yi,zi)

                            a = label(xi,y_down,z_down)
                        6   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 6
                            else
                                label(xi,y_down,z_down) = a
                            endif
                        endif

                        if(label(x_down,yi,z_down) > 0) then
                            klass(label(x_down,yi,z_down)) = label(xi,yi,zi)

                            a = label(x_down,yi,z_down)
                        7   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 7
                            else
                                label(x_down,yi,z_down) = a
                            endif
                        endif

                    endif
                endif
            enddo
        enddo
    enddo

    open(30,file="./process.d",action="write",position="append")
    write(30,*) "1th is OK!!!", step
    close(30)

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

    open(30,file="./process.d",action="write",position="append")
    write(30,*) "2th is OK!!!", step
    close(30)


    Vlabel(:) = 0.0d0
    dropnum = 0
    DO i = 1, 1000
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
        if(Vlabel(i) > 0.0d0) then
            dropnum = dropnum + 1
        endif
    ENDDO
    ! open(30,file="./process.d",action="write",position="append")
    ! write(30,*) "radius is cal!!!", step
    ! close(30)

    ! dropnum = 0
    ! do i = 1, 10000
    !     if(Vlabel(i) > 0.0d0) then
    !         dropnum = dropnum + 1
    !     endif
    ! enddo
    open(30,file="./process.d",action="write",position="append")
    write(30,*) "dropnum is cal!!!", step
    close(30)

    ! if(step == 5000) then
        open(10,file="./drop_num.d")
        write(10,*) dble(step), dble(dropnum)
        close(10)
    ! else
        ! open(10,file="./drop_num.d",action="write",position="append")
        ! write(10,*) dble(step), dble(dropnum)
        ! close(10)
    ! endif

    ! open(11,file=datadir2//"radius_"//trim(file_num)//".d")
    ! do i = 1, 10000
    !     if(radius(i) > 0.0d0) then
    !         write(11,*) radius(i)
    !     endif
    ! enddo
    ! close(11)

    if(step == 5000) then
        open(12,file="./vall.d")
        write(12,*) dble(step), Vall
        close(12)
    else
        open(12,file="./vall.d",action="write",position="append")
        write(12,*) dble(step), Vall
        close(12)
    endif
    call cpu_time(time2)
    open(30,file="./process.d",action="write",position="append")
    write(30,*) "time is!!!", time2-time1
    close(30)
! ENDDO
end program main
