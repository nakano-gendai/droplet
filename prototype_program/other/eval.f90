program main
    implicit none
    integer, parameter :: n = 129 !データの行数
    integer, parameter :: m = 129
    integer,parameter:: xmax = 128 !ｘ方向格子数
    integer,parameter:: ymax = 64
    integer,parameter:: zmax = 128 !ｚ方向格子数
    real(8),parameter:: phi1 = 2.211d0
    real(8),parameter:: phi2 = 4.895d0
    integer,parameter:: num = 300000
    real(8),parameter:: uw = 0.035d0
    real(8) x, z
    real(8) u1(0:xmax, 0:ymax, 0:zmax)
    real(8) u2(0:xmax, 0:ymax, 0:zmax)
    real(8) u3(0:xmax, 0:ymax, 0:zmax)
    real(8) p(0:xmax, 0:ymax, 0:zmax)
    real(8) phi(0:xmax, 0:ymax, 0:zmax)
    real(8) u_m(0:xmax, 0:ymax, 0:zmax)
    integer xi, yi, zi, i, j, k, step

    !重心計算用変数
    real(8) xg, yg, zg
    real(8) phi_tmp
    !線形補間用変数
    real(8),parameter :: xx = 0.01d0
    real(8),parameter :: zz = 0.01d0
    real(8) xj, zj, phi_hokann
    integer x0, x1, x2, z0, z1, z2
    integer x_tmp, z_tmp
    real(8),parameter:: pi = acos(-1.0d0)
    real(8) distance

    !液滴の変形
    real(8) L1, B1, Dd

    character(8) file_num
    character(*),parameter :: dir = "/home/nakano/anime/data_test20/"

    open(10,file="eval20.d")
DO step = 100, num, 100

    yi = ymax / 2
    !ファイル読み込み
    u1(:,:,:) = 0.0d0
    u2(:,:,:) = 0.0d0
    u3(:,:,:) = 0.0d0
    u_m(:,:,:) = 0.0d0
    p(:,:,:) = 0.0d0
    phi(:,:,:) = 0.0d0
    if((step > 99) .and. (step < 1000)) then
        write(file_num, "(i3)") step
    elseif((step > 999) .and. (step < 10000)) then
        write(file_num,"(i4)") step
    elseif((step > 9999) .and. (step < 100000)) then
        write(file_num,"(i5)") step
    elseif((step > 99999) .and. (step < 1000000)) then
        write(file_num,"(i6)") step
    endif
    open(20, file=dir//"0_"trim(file_num)//".d")
    do zi = 0, m-1
        do xi = 0, n-1
            read(20,"(3es16.8)") x, z, phi(xi,yi,zi)
        enddo
        read(20,"()")
    enddo
    close(20)
    
    !==========================重心座標を求める===============================================
    xg = 0.0d0
    yg = 0.0d0
    zg = 0.0d0
    phi_tmp = 0.0d0
    
    do zi=0,zmax
        ! do yi=0,ymax
            do xi=0,xmax
                xg = xg + dble(xi) * (phi(xi,yi,zi) - phi1)
                yg = yg + dble(yi) * (phi(xi,yi,zi) - phi1)
                zg = zg + dble(zi) * (phi(xi,yi,zi) - phi1)
                phi_tmp = phi_tmp + (phi(xi,yi,zi) - phi1)
            enddo
        ! enddo
    enddo
    xg = xg / phi_tmp
    yg = yg / phi_tmp
    zg = zg / phi_tmp
    !========================線形補間========================================================
    L1 = 0.0d0
    B1 = 2.0d0 * dble(xmax)
    do zi = 0, m-2
        do xi = 0, n-2
            !座標設定
            x0 = int(xi) - 1
            x1 = int(xi)
            x2 = int(xi) + 1
            z0 = int(zi) - 1
            z1 = int(zi)
            z2 = int(zi) + 1
            !横
            k = 0
            do 
                k = k + 1
                xj = dble(xi) + dble(k) * xx
                zj = dble(zi)
                phi_hokann = ((dble(x2)-xj)*(dble(z2)-zj)*phi(x1,yi,z1) + (xj-dble(x1))*(dble(z2)-zj)*phi(x2,yi,z1) &
                                + (dble(x2)-xj)*(zj-dble(z1))*phi(x1,yi,z2) + (xj-dble(x1))*(zj-dble(z1))*phi(x2,yi,z2)) &
                                / ((dble(x2)-dble(x1))*(dble(z2)-dble(z1)))
                if((phi_hokann >= 0.5d0*(phi1+phi2)) .and. (0.5d0*(phi1+phi2)+0.1d0 >= phi_hokann )) then
                    distance = sqrt((xj-xg)**2 + (zj-zg)**2) * 2.0d0
                    if(L1 <= distance) then
                        L1 = distance
                    elseif(B1 >= distance) then
                        B1 = distance
                    else
                        L1 = L1
                        B1 = B1
                    endif
                endif
                if(xj > dble(x2)) exit
            enddo
            !縦
            k = 0
            do 
                k = k + 1
                xj = dble(xi)
                zj = dble(zi) + dble(k) * zz
                phi_hokann = ((dble(x2)-xj)*(dble(z2)-zj)*phi(x1,yi,z1) + (xj-dble(x1))*(dble(z2)-zj)*phi(x2,yi,z1) &
                                + (dble(x2)-xj)*(zj-dble(z1))*phi(x1,yi,z2) + (xj-dble(x1))*(zj-dble(z1))*phi(x2,yi,z2)) &
                                / ((dble(x2)-dble(x1))*(dble(z2)-dble(z1)))
                if((phi_hokann >= 0.5d0*(phi1+phi2)) .and. (0.5d0*(phi1+phi2)+0.1d0 >= phi_hokann )) then
                    distance = sqrt((xj-xg)**2 + (zj-zg)**2) * 2.0d0
                    if(L1 <= distance) then
                        L1 = distance
                    elseif(B1 >= distance) then
                        B1 = distance
                    else
                        L1 = L1
                        B1 = B1
                    endif
                endif
                if(zj > dble(z2)) exit
            enddo
            !右斜め
            k = 0
            do 
                k = k + 1
                xj = dble(xi) + dble(k) * xx
                zj = dble(zi) + dble(k) * zz
                phi_hokann = ((dble(x2)-xj)*(dble(z2)-zj)*phi(x1,yi,z1) + (xj-dble(x1))*(dble(z2)-zj)*phi(x2,yi,z1) &
                                + (dble(x2)-xj)*(zj-dble(z1))*phi(x1,yi,z2) + (xj-dble(x1))*(zj-dble(z1))*phi(x2,yi,z2)) &
                                / ((dble(x2)-dble(x1))*(dble(z2)-dble(z1)))
                if((phi_hokann >= 0.5d0*(phi1+phi2)) .and. (0.5d0*(phi1+phi2)+0.1d0 >= phi_hokann )) then
                    distance = sqrt((xj-xg)**2 + (zj-zg)**2) * 2.0d0
                    if(L1 <= distance) then
                        L1 = distance
                    elseif(B1 >= distance) then
                        B1 = distance
                    else
                        L1 = L1
                        B1 = B1
                    endif
                endif
                if(zj > dble(z2)) exit
            enddo
        enddo
    enddo

    do zi = 0, m-2
        do xi = 1, n-1
            !座標設定
            x0 = int(xi) - 1
            x1 = int(xi)
            x2 = int(xi) + 1
            z0 = int(zi) - 1
            z1 = int(zi)
            z2 = int(zi) + 1
            !左斜め
            k = 0
            do 
                k = k + 1
                xj = dble(xi) - dble(k) * xx
                zj = dble(zi) + dble(k) * zz
                phi_hokann = ((dble(x1)-xj)*(dble(z2)-zj)*phi(x0,yi,z1) + (xj-dble(x0))*(dble(z2)-zj)*phi(x1,yi,z1) &
                                + (dble(x1)-xj)*(zj-dble(z1))*phi(x0,yi,z2) + (xj-dble(x0))*(zj-dble(z1))*phi(x1,yi,z2)) &
                                / ((dble(x1)-dble(x0))*(dble(z2)-dble(z1)))
                if((phi_hokann >= 0.5d0*(phi1+phi2)) .and. (0.5d0*(phi1+phi2)+0.1d0 >= phi_hokann )) then
                    distance = sqrt((xj-xg)**2 + (zj-zg)**2) * 2.0d0
                    if(L1 <= distance) then
                        L1 = distance
                    elseif(B1 >= distance) then
                        B1 = distance
                    else
                        L1 = L1
                        B1 = B1
                    endif
                endif
                if(zj > dble(z2)) exit
            enddo
        enddo
    enddo

    !変形度
    Dd = abs(L1-B1) / (L1+B1)
    !出力
    write(10,"(4es16.8)") (dble(step)-dble(5000))*2.0d0*uw/dble(128),L1,B1,Dd
    write(*,*) "step =", step
ENDDO
    close(10)
end program main

