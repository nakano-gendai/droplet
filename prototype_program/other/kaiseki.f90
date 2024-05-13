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
    real(8),parameter:: uw = 0.021d0
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
    !慣性行列の固有値・固有ベクトル計算用変数
    real(8) xx_wa, yy_wa, zz_wa, xy_wa, yz_wa, zx_wa
    real(8) inertia_matrix(1:3,1:3), inertia_matrix_tmp(1:3,1:3)
    real(8) principal_moment_1, principal_moment_2, principal_moment_3
    real(8) inertia_principal_axis_1(1:3), inertia_principal_axis_2(1:3), inertia_principal_axis_3(1:3)
    real(8) w(1:3,1:3), w_t(1:3,1:3), t0(1:3,1:3), t_tmp(1:3,1:3)
    real(8) theta, mx
    integer pp, qq, check
    !液滴の変形
    real(8) L1, B1, Dd
    !線形補間用変数
    real(8), parameter :: weight = 0.01d0
    real(8) xj, zj, phi_hokann
    integer x1, x2, z1, z2

    real(8),parameter :: zz = 0.01d0
    real(8),parameter :: xx = 0.01d0
    integer x_tmp, z_tmp
    real(8),parameter:: pi = acos(-1.0d0)
    real(8) distance

    character(8) file_num
    character(*),parameter :: dir = "/home/nakano/anime/data_test9/"

    open(10,file="test9.d")
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
    open(20, file=dir//trim(file_num)//".d")
    do zi = 0, m-1
        do xi = 0, n-1
            read(20,"(8es16.8)") x, z, u1(xi,yi,zi), u2(xi,yi,zi), u3(xi,yi,zi), u_m(xi,yi,zi), p(xi,yi,zi), phi(xi,yi,zi)
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
    !==========================慣性行列を求める===============================================
    ! xx_wa = 0.0d0
    ! yy_wa = 0.0d0
    ! zz_wa = 0.0d0
    ! xy_wa = 0.0d0
    ! yz_wa = 0.0d0
    ! zx_wa = 0.0d0
    ! do zi=0,zmax
    !     ! do yi=0,ymax
    !         do xi=0,xmax
    !             xx_wa = xx_wa + (phi(xi,yi,zi)-phi1)*(dble(xi)-xg)*(dble(xi)-xg)
    !             yy_wa = yy_wa + (phi(xi,yi,zi)-phi1)*(dble(yi)-yg)*(dble(yi)-yg)
    !             zz_wa = zz_wa + (phi(xi,yi,zi)-phi1)*(dble(zi)-zg)*(dble(zi)-zg)
    !             xy_wa = xy_wa + (phi(xi,yi,zi)-phi1)*(dble(xi)-xg)*(dble(yi)-yg)
    !             yz_wa = yz_wa + (phi(xi,yi,zi)-phi1)*(dble(yi)-yg)*(dble(zi)-zg)
    !             zx_wa = zx_wa + (phi(xi,yi,zi)-phi1)*(dble(zi)-zg)*(dble(xi)-xg)
    !         enddo
    !     ! enddo
    ! enddo
    
    ! !慣性行列
    ! inertia_matrix(1,1) = yy_wa + zz_wa
    ! inertia_matrix(1,2) = -xy_wa
    ! inertia_matrix(1,3) = -zx_wa
    ! inertia_matrix(2,1) = -xy_wa
    ! inertia_matrix(2,2) = xx_wa + zz_wa
    ! inertia_matrix(2,3) = -yz_wa
    ! inertia_matrix(3,1) = -zx_wa
    ! inertia_matrix(3,2) = -yz_wa
    ! inertia_matrix(3,3) = xx_wa + yy_wa

    ! !==============固有値, 固有ベクトルを求める（ヤコビ法）=======================================
    ! do j=1,3
    !     do i=1,3
    !         t0(i,j) = 0.0d0
    !         !対角成分だけ1.0を代入
    !         if(i == j) then
    !             t0(i,j) = 1.0d0  
    !         endif
    !     enddo
    ! enddo

    ! do 
    !     !非対角成分の中から絶対値の最大成分を探す
    !     mx = abs(inertia_matrix(1,2))
    !     pp = 1
    !     qq = 2
    !     do j = 2,3
    !         do i = 1,j-1
    !             if(abs(inertia_matrix(i,j)) > mx) then
    !                 mx = abs(inertia_matrix(i,j))
    !                 pp = i
    !                 qq = j
    !             endif
    !         enddo
    !     enddo
    !     !thetaを求める
    !     theta = 0.5d0 * atan(2.0d0 * inertia_matrix(pp,qq) / (inertia_matrix(pp,pp) - inertia_matrix(qq,qq)))
    !     !（転置）直交行列を求める
    !     do j=1,3
    !         do i=1,3
    !             w(i,j) = 0.0d0
    !             !対角成分だけ1.0を代入
    !             if(i == j) then
    !                 w(i,j) = 1.0d0  
    !             endif
    !         enddo
    !     enddo
    !     w(pp,pp) = cos(theta)
    !     w(pp,qq) = -sin(theta)
    !     w(qq,pp) = sin(theta)
    !     w(qq,qq) = cos(theta)
    !     do j=1,3
    !         do i=1,3
    !             w_t(i,j) = w(j,i)
    !         enddo
    !     enddo
    !     !新しい行列aを計算する
    !     inertia_matrix_tmp(:,:) = 0.0d0
    !     inertia_matrix_tmp(:,:) = matmul(w_t, inertia_matrix)
    !     inertia_matrix(:,:) = matmul(inertia_matrix_tmp, w)
    !     !並行して固有ベクトルを計算するための行列も計算する
    !     t_tmp(:,:) = matmul(t0, w)
    !     t0(:,:) = t_tmp(:,:)
    !     !非対角成分がすべて0.0に近づいたかを調べて、収束判定する
    !     check = 0
    !     do j=1,3
    !         do i=j+1,3
    !             if(abs(inertia_matrix(i,j)) < 1.0d-6) then
    !                 check = check + 1
    !             endif
    !         enddo
    !     enddo
    !     if(check == (3 * 3 - 3) / 2) exit
    ! enddo
    ! !固有値
    ! principal_moment_1 = inertia_matrix(1,1)
    ! principal_moment_2 = inertia_matrix(2,2)
    ! principal_moment_3 = inertia_matrix(3,3)
    ! !固有ベクトル
    ! inertia_principal_axis_1(1:3) = t0(1:3,1)
    ! inertia_principal_axis_2(1:3) = t0(1:3,2)
    ! inertia_principal_axis_3(1:3) = t0(1:3,3)

    ! ! write(*,*) principal_moment_1,principal_moment_2,principal_moment_3
    ! ! 長辺 / 短辺
    ! L1 = 2.0d0 * sqrt(principal_moment_3 / phi_tmp)
    ! B1 = 2.0d0 * sqrt(principal_moment_1 / phi_tmp)


    !=================================線形補間================================================
        ! 短辺を計算
        ! k = 0
        ! do
        !     k = k + 1
        !     xj = inertia_principal_axis_3(1) * dble(k) * weight
        !     zj = inertia_principal_axis_3(3) * dble(k) * weight
        !     xj = xj + xg
        !     zj = zj + zg
        !     x1 = int(xj)
        !     x2 = int(xj) + 1
        !     z1 = int(zj)
        !     z2 = int(zj) + 1
        !     !線形補間
        !     phi_hokann = ((dble(x2)-xj)*(dble(z2)-zj)*phi(x1,yi,z1) + (xj-dble(x1))*(dble(z2)-zj)*phi(x2,yi,z1) &
        !                     + (dble(x2)-xj)*(zj-dble(z1))*phi(x1,yi,z2) + (xj-dble(x1))*(zj-dble(z1))*phi(x2,yi,z2)) &
        !                     / ((dble(x2)-dble(x1))*(dble(z2)-dble(z1)))
        !     if(phi_hokann >= 0.5d0*(phi1+phi2)) then
        !         B1 = sqrt((xj-xg)**2 + (zj-zg)**2) * 2.0d0
        !     endif
        !     if(phi_hokann < 0.5d0*(phi1+phi2)) exit
        ! enddo
        ! !長辺を計算
        ! k = 0
        ! do
        !     k = k + 1
        !     xj = inertia_principal_axis_1(1) * dble(k) * weight
        !     zj = inertia_principal_axis_1(3) * dble(k) * weight
        !     xj = xj + xg
        !     zj = zj + zg
        !     x1 = int(xj)
        !     x2 = int(xj) + 1
        !     z1 = int(zj)
        !     z2 = int(zj) + 1
        !     !線形補間
        !     phi_hokann = ((dble(x2)-xj)*(dble(z2)-zj)*phi(x1,yi,z1) + (xj-dble(x1))*(dble(z2)-zj)*phi(x2,yi,z1) &
        !                     + (dble(x2)-xj)*(zj-dble(z1))*phi(x1,yi,z2) + (xj-dble(x1))*(zj-dble(z1))*phi(x2,yi,z2)) &
        !                     / ((dble(x2)-dble(x1))*(dble(z2)-dble(z1)))
        !     if(phi_hokann >= 0.5d0*(phi1+phi2)) then
        !         L1 = sqrt((xj-xg)**2 + (zj-zg)**2) * 2.0d0
        !     endif
        !     if(phi_hokann < 0.5d0*(phi1+phi2)) exit
        ! enddo

    ! ! 回転角を計算
    !     distance = 0.0d0
    !     yi = (ymax)/2
    !     do zi=0,zmax
    !         do xi=0,xmax
    !             if(phi(xi,yi,zi) >= 0.5d0*(phi1+phi2)) then
    !                 if( distance <= sqrt((dble(xi)-xg)**2 + (dble(zi)-zg)**2) ) then
    !                     distance = sqrt((dble(xi)-xg)**2 + (dble(zi)-zg)**2)
    !                     x_tmp = xi
    !                     z_tmp = zi
    !                 endif
    !             endif
    !         enddo
    !     enddo
    !     theta = acos(abs((dble(x_tmp) - xg) / distance))
    !     if(((dble(x_tmp)-xg) > 0.0d0).and.((dble(z_tmp)-zg) > 0.0d0)) then !第一象限
    !         theta = theta
    !     else if(((dble(x_tmp)-xg) < 0.0d0).and.((dble(z_tmp)-zg) > 0.0d0)) then !第二象限
    !         theta = pi - theta
    !     else if(((dble(x_tmp)-xg) < 0.0d0).and.((dble(z_tmp)-zg) < 0.0d0)) then !第三象限
    !         theta = pi + theta
    !     else if(((dble(x_tmp)-xg) > 0.0d0).and.((dble(z_tmp)-zg) < 0.0d0)) then !第四象限
    !         theta = 2.0d0 * pi - theta 
    !     endif
    !     !短辺を計算
    !     k = 0
    !     do
    !         k = k + 1
    !         zj = (zz*dble(k)) / (tan(theta)*sin(theta) + cos(theta))
    !         xj = -zj*tan(theta)
    !         xj = xj + xg
    !         zj = zj + zg
    !         x1 = int(xj)
    !         x2 = int(xj) + 1
    !         z1 = int(zj)
    !         z2 = int(zj) + 1
    !         !線形補間
    !         phi_hokann = ((dble(x2)-xj)*(dble(z2)-zj)*phi(x1,yi,z1) + (xj-dble(x1))*(dble(z2)-zj)*phi(x2,yi,z1) &
    !                         + (dble(x2)-xj)*(zj-dble(z1))*phi(x1,yi,z2) + (xj-dble(x1))*(zj-dble(z1))*phi(x2,yi,z2)) &
    !                         / ((dble(x2)-dble(x1))*(dble(z2)-dble(z1)))
    !         if(phi_hokann >= 0.5d0*(phi1+phi2)) then
    !             B1 = sqrt((xj-xg)**2 + (zj-zg)**2) * 2.0d0
    !         endif
    !         if(phi_hokann < 0.5d0*(phi1+phi2)) exit
    !     enddo
    !     !長辺を計算
    !     k = 0
    !     do
    !         k = k + 1
    !         xj = (xx*dble(k)) / (tan(theta)*sin(theta) + cos(theta))
    !         zj = xj*tan(theta)
    !         xj = xj + xg
    !         zj = zj + zg
    !         x1 = int(xj)
    !         x2 = int(xj) + 1
    !         z1 = int(zj)
    !         z2 = int(zj) + 1
    !         !線形補間
    !         phi_hokann = ((dble(x2)-xj)*(dble(z2)-zj)*phi(x1,yi,z1) + (xj-dble(x1))*(dble(z2)-zj)*phi(x2,yi,z1) &
    !                         + (dble(x2)-xj)*(zj-dble(z1))*phi(x1,yi,z2) + (xj-dble(x1))*(zj-dble(z1))*phi(x2,yi,z2)) &
    !                         / ((dble(x2)-dble(x1))*(dble(z2)-dble(z1)))
    !         if(phi_hokann >= 0.5d0*(phi1+phi2)) then
    !             L1 = sqrt((xj-xg)**2 + (zj-zg)**2) * 2.0d0
    !         endif
    !         if(phi_hokann < 0.5d0*(phi1+phi2)) exit
    !     enddo

    !変形度
    Dd = abs(L1-B1) / (L1+B1)

    write(10,"(4es16.8)") (dble(step)-dble(5000))*2.0d0*uw/dble(128),L1,B1,Dd
    ! write(10,"(4es16.8)") dble(step),L1,B1,Dd
    
ENDDO
    close(10)
end program main

