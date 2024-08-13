program main
    real(8),parameter:: ds = 1.0d0 !格子間隔（lattice unit）
    integer,parameter:: xmax = 255 !ｘ方向格子数（０から数える）
    integer,parameter:: ymax = 255 !ｙ方向格子数（０から数える）
    integer,parameter:: zmax = 255 !ｚ方向格子数（０から数える）
    real(8),parameter:: phi1 = 2.638d-1 !連続相のオーダーパラメータ
    real(8),parameter:: phi2 = 4.031d-1 !分散相のオーダーパラメータ
    integer,parameter:: case_initial_num = 49 !最初のケース番号
    integer,parameter:: case_end_num = 50!最後のケース番号

    real(8) phi(-1:xmax+1,-1:ymax+1,-1:zmax+1)
    integer label(0:xmax,0:ymax,0:zmax)
    integer klass(1:2000)
    integer tmp(1:7)
    integer labeltmp
    integer xi, yi, zi, i, k, a, b, x_up, y_up, z_up, x_down, y_down, z_down
    real(8) Vlabel(1:2000), radius(1:2000)
    real(8) Vall

    character(*),parameter :: datadir = "/data/sht/nakanog/IHT_drop_d70_we20/"
    character(*),parameter :: datadir2 = "/data/sht/nakanog/IHT_drop_d70_we20/drop_num/"
    integer step, dropnum, case_num
    character(8) file_num, file_num2
    real(8) time1, time2

    real(8),parameter:: D = 70.0d0
    real(8),parameter:: epsilon = 9.15d-10 !エネルギー散逸率
    real(8),parameter:: eddytime = epsilon**(-1.0d0/3.0d0)*D**(2.0d0/3.0d0)

    !初期液滴中心
    real(8), parameter:: xc = 0.5d0*dble(xmax+1)	
    real(8), parameter:: yc = 0.5d0*dble(ymax+1)
    real(8), parameter:: zc = 0.5d0*dble(zmax+1)

    open(30,file="./process.d")
    write(30,*) "percolation start!!!!"
    close(30)
DO case_num = case_initial_num, case_end_num
DO step = 5000, 50000, 1000
    ! step = 48000
    call cpu_time(time1)

    if((case_num > 0) .and. (case_num < 10)) then
        write(file_num2,"(i1)") case_num
    elseif((case_num > 9) .and. (case_num < 100)) then
        write(file_num2,"(i2)") case_num
    elseif((case_num > 99) .and. (case_num < 1000)) then
        write(file_num2,"(i3)") case_num
    endif

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
    elseif((step > 9999999) .and. (step < 100000000)) then
        write(file_num,"(i8)") step
    endif

    open(20, file=datadir//trim(file_num2)//"_"//trim(file_num)//".bin", form="unformatted")
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                read(20)  phi(xi,yi,zi)
            enddo
        enddo
    enddo
    close(20)

    open(30,file="./process.d",action="write",position="append")
    write(30,*) "input is OK!!!", step
    close(30)

    k = 1
    label(:,:,:) = 0
    klass(:) = 0
    Vall = 0.0d0

    do yi = 0, ymax
        do zi = 0, zmax
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
                        11   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 11
                            else
                                label(x_down,yi,zi) = a
                            endif
                        endif

                        if(label(xi,y_down,zi) > 0) then
                            klass(label(xi,y_down,zi)) = label(xi,yi,zi)

                            a = label(xi,y_down,zi)
                        12   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 12
                            else
                                label(xi,y_down,zi) = a
                            endif
                        endif

                        if(label(xi,yi,z_down) > 0) then
                            klass(label(xi,yi,z_down)) = label(xi,yi,zi)

                            a = label(xi,yi,z_down)
                        13   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 13
                            else
                                label(xi,yi,z_down) = a
                            endif
                        endif

                        if(label(xi,y_up,z_down) > 0) then
                            klass(label(xi,y_up,z_down)) = label(xi,yi,zi)

                            a = label(xi,y_up,z_down)
                        14   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 14
                            else
                                label(xi,y_up,z_down) = a
                            endif
                        endif

                        if(label(x_up,yi,z_down) > 0) then
                            klass(label(x_up,yi,z_down)) = label(xi,yi,zi)

                            a = label(x_up,yi,z_down)
                        15   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 15
                            else
                                label(x_up,yi,z_down) = a
                            endif
                        endif

                        if(label(xi,y_down,z_down) > 0) then
                            klass(label(xi,y_down,z_down)) = label(xi,yi,zi)

                            a = label(xi,y_down,z_down)
                        16   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 16
                            else
                                label(xi,y_down,z_down) = a
                            endif
                        endif

                        if(label(x_down,yi,z_down) > 0) then
                            klass(label(x_down,yi,z_down)) = label(xi,yi,zi)

                            a = label(x_down,yi,z_down)
                        17   b = a
                            a = klass(a)
                            if(a /= b)  then 
                                go to 17
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
    DO i = 1, 2000
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

    if(step == 5000) then
        open(10,file=datadir2//trim(file_num2)//"_num.d")
        write(10,*) dble(step), (dble(step)-5000.0d0)/eddytime, dble(dropnum)
        close(10)
    else
        open(10,file=datadir2//trim(file_num2)//"_num.d",action="write",position="append")
        write(10,*) dble(step), (dble(step)-5000.0d0)/eddytime, dble(dropnum)
        close(10)
    endif

    ! open(11,file=datadir2//"radius_"//trim(file_num)//".d")
    ! do i = 1, 10000
    !     if(radius(i) > 0.0d0) then
    !         write(11,*) radius(i)
    !     endif
    ! enddo
    ! close(11)

    ! if(step == 5000) then
    !     open(12,file="./vall.d")
    !     write(12,*) dble(step), Vall
    !     close(12)
    ! else
    !     open(12,file="./vall.d",action="write",position="append")
    !     write(12,*) dble(step), Vall
    !     close(12)
    ! endif
    call cpu_time(time2)
    open(30,file="./process.d",action="write",position="append")
    write(30,*) "time is!!!", time2-time1
    close(30)
ENDDO
ENDDO
end program main
