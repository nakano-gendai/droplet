module globals
    integer, parameter :: n = 300 !データの行数
    integer, parameter :: m = 72
    integer,parameter:: xmax = 299 !ｘ方向格子数
    integer,parameter:: ymax = 64
    integer,parameter:: zmax = 71 !ｚ方向格子数
    real(8),parameter:: phi1 = 2.211d0
    real(8),parameter:: phi2 = 4.895d0
    real(8) phi_kaimen 
    integer,parameter:: num = 140000
    real(8),parameter:: uw = 0.02d0
    character(8) file_num
    character(*),parameter :: dir = "/home/nakano/anime/test/"
    ! character(*),parameter :: dir = "C:\Users\genda\Documents\data_re0.2ca0.5_ini\"

contains
    subroutine senkei(xx1,xx2,zz1,zz2,xxg,zzg,L1,B1,phi)
        integer, intent(in) :: xx1,xx2,zz1,zz2
        real(8), intent(in) :: xxg,zzg
        real(8), intent(inout) :: L1, B1
        real(8), intent(in) :: phi(0:xmax, 0:ymax, 0:zmax)
        real(8) distance
        real(8) alpha, l
        real(8) xj, zj
        
        

        l = sqrt((dble(xx2)-dble(xx1))**2 +(dble(zz2)-dble(zz1))**2)
        alpha = l * (phi_kaimen - phi(xx1,ymax/2,zz1)) / (phi(xx2,ymax/2,zz2) - phi(xx1,ymax/2,zz1))
        if((alpha >= 0.0d0) .and. (alpha <= l)) then
            xj = dble(xx1) + alpha/l * (dble(xx2) - dble(xx1))
            zj = dble(zz1) + alpha/l * (dble(zz2) - dble(zz1))
            distance = sqrt((xj-xxg)**2 + (zj-zzg)**2)
            if(L1 <= distance) then
                L1 = distance
            elseif(B1 >= distance) then
                B1 = distance
            else
                L1 = L1
                B1 = B1
            endif
        endif
    end subroutine senkei
end module globals

program main
use globals
    implicit none
    real(8) x, z
    real(8) u1(0:xmax, 0:ymax, 0:zmax)
    real(8) u2(0:xmax, 0:ymax, 0:zmax)
    real(8) u3(0:xmax, 0:ymax, 0:zmax)
    real(8) p(0:xmax, 0:ymax, 0:zmax)
    real(8) phi(0:xmax, 0:ymax, 0:zmax)
    real(8) u_m(0:xmax, 0:ymax, 0:zmax)
    integer xi, yi, zi, i, j, k, step
    real(8) phi_max, phi_min

    !重心計算用変数
    real(8) xg, yg, zg
    real(8) phi_tmp
    !線形補間用変数
    integer x0, x1, x2, z0, z1, z2
    real(8),parameter:: pi = acos(-1.0d0)
    

    !液滴の変形
    real(8) L1, B1, Dd

    open(10,file="re0.2ca0.5ini.d")
DO step = 100, num, 100

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
    open(20, file=dir//trim(file_num)//".d")
    yi = ymax/2
    do zi = 0, m-1
        do xi = 0, n-1
            ! read(20,"(8es16.8)") x, z, u1(xi,yi,zi), u2(xi,yi,zi), u3(xi,yi,zi), u_m(xi,yi,zi),p(xi,yi,zi), phi(xi,yi,zi)
            read(20,"(3es16.8)") x, z, phi(xi,yi,zi)
        enddo
        read(20,"()")
    enddo
    close(20)

    phi_kaimen = 0.5d0 * (phi1 + phi2)

    !==========================重心座標を求める===============================================
    xg = 0.0d0
    yg = 0.0d0
    zg = 0.0d0
    phi_tmp = 0.0d0
    yi = ymax/2
    do zi=0,zmax
        ! do yi=0,ymax
            do xi=0,(xmax+1)/2
                if(phi(xi,yi,zi) >= phi_kaimen) then
                xg = xg + dble(xi)
                yg = yg + dble(yi) 
                zg = zg + dble(zi) 
                phi_tmp = phi_tmp + 1.0d0
                endif
            enddo
        ! enddo
    enddo

    xg = xg / phi_tmp 
    yg = yg / phi_tmp
    zg = zg / phi_tmp
    !========================線形補間========================================================
    L1 = 0.0d0
    B1 = 2.0d0 * dble(xmax)
    
    do zi = 1, m-2
        do xi = (xmax+1)/2, n-2
        ! do xi=0,(xmax+1)/2
            !座標設定
            x0 = int(xi) - 1
            x1 = int(xi)
            x2 = int(xi) + 1
            z0 = int(zi) - 1
            z1 = int(zi)
            z2 = int(zi) + 1
            
            !右
            call senkei(x1,x2,z1,z1,xg,zg,L1,B1,phi)
            !左
            call senkei(x1,x0,z1,z1,xg,zg,L1,B1,phi)
            !上
            call senkei(x1,x1,z1,z2,xg,zg,L1,B1,phi)
            !下
            call senkei(x1,x1,z1,z0,xg,zg,L1,B1,phi)
            !右斜め上
            call senkei(x1,x2,z1,z2,xg,zg,L1,B1,phi)
            !右斜め下
            call senkei(x1,x2,z1,z0,xg,zg,L1,B1,phi)
            !左斜め上
            call senkei(x1,x0,z1,z2,xg,zg,L1,B1,phi)
            !左斜め下
            call senkei(x1,x0,z1,z0,xg,zg,L1,B1,phi)
        enddo
    enddo


    !変形度
    Dd = abs(L1-B1) / (L1+B1)
    !出力
    write(10,"(6es16.8)") (dble(step)-dble(5000))*2.0d0*uw/dble(zmax),L1,B1,Dd,xg,zg

    write(*,*) "step =", step
ENDDO
    close(10)
end program main

