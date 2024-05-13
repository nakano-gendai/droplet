module globals
    integer,parameter:: xmax = 999 !ｘ方向格子数
    integer,parameter:: ymax = 220
    integer,parameter:: zmax = 220 !ｚ方向格子数
    integer,parameter:: Nx = 100 !ｘ方向の並列数
    integer,parameter:: Nz = 17 !ｚ方向の並列数
    integer,parameter:: x_procs = (xmax+1) / Nx
    integer,parameter:: z_procs = (zmax+1) / Nz
    integer,parameter:: new_procs = Nx * Nz
    real,parameter:: phi1 = 2.211d0
    real,parameter:: phi2 = 4.895d0
    real,parameter:: phi_kaimen = 0.5d0 * (phi1 + phi2)
    integer,parameter:: num = 100000
    real,parameter:: uw = 0.08d0
    character(8) file_num, file_num2, dropnum
    integer set

    character(*),parameter :: datadir = "/data/n/n517/newtest/re100we5/"

contains
    subroutine senkei(xx1,xx2,zz1,zz2,xxg,zzg,Ltmp,Btmp,phi)
        integer, intent(in) :: xx1,xx2,zz1,zz2
        real, intent(in) :: xxg,zzg
        real, intent(inout) :: Ltmp, Btmp
        real, intent(in) :: phi(0:xmax, 0:zmax)
        real distance
        real alpha, l
        real xj, zj
        
        l = sqrt((dble(xx2)-dble(xx1))**2 +(dble(zz2)-dble(zz1))**2)
        alpha = l * (phi_kaimen - phi(xx1,zz1)) / (phi(xx2,zz2) - phi(xx1,zz1))
        if((alpha >= 0.0d0) .and. (alpha <= l)) then
            xj = dble(xx1) + alpha/l * (dble(xx2) - dble(xx1))
            zj = dble(zz1) + alpha/l * (dble(zz2) - dble(zz1))
            distance = sqrt((xj-xxg)**2 + (zj-zzg)**2)
            if(Ltmp <= distance) then
                Ltmp = distance
            elseif(Btmp >= distance) then
                Btmp = distance
            else
                Ltmp = Ltmp
                Btmp = Btmp
            endif
        endif
    end subroutine senkei

end module globals

program main
use globals
    implicit none
    real x, y, z, s
    real phiout(1:x_procs,1:z_procs,0:new_procs-1)
    real phi(0:xmax, 0:zmax)
    integer i, k, xi, yi, zi, Nxx, Nzz, step

    !重心計算用変数
    real xg, yg, zg
    real phi_tmp
    !線形補間用変数
    integer x0, x1, x2, z0, z1, z2
    real Ltmp, Btmp
    !判定用変数
    integer frag, fragformer, knum, inum, dropnumtmp
    integer xstart(100), xend(100)
    real L1(100), B1(100), Dd(100)

DO step = 4000, num, 2000
    if((step > 99) .and. (step < 1000)) then
        write(file_num2, "(i3)") step
    elseif((step > 999) .and. (step < 10000)) then
        write(file_num2,"(i4)") step
    elseif((step > 9999) .and. (step < 100000)) then
        write(file_num2,"(i5)") step
    elseif((step > 99999) .and. (step < 1000000)) then
        write(file_num2,"(i6)") step
    endif
    do i=0,new_procs-1
        if((i >= 0) .and. (i < 10)) then
            write(file_num, "(i1)") i
        elseif((i > 9) .and. (i < 100)) then
            write(file_num,"(i2)") i
        elseif((i > 99) .and. (i < 1000)) then
            write(file_num,"(i3)") i
        elseif((i > 999) .and. (i < 10000)) then
            write(file_num,"(i4)") i
        endif

        set = 20 + i
        open(set, file=datadir//"0_"//trim(file_num)//"_"//trim(file_num2)//".bin", form="unformatted")
        do zi = 1, z_procs
            do xi = 1, x_procs
                read(set) phiout(xi,zi,i), x, y, z, s
            enddo
            read(set)
        enddo
        close(set)
    enddo
    k = 0
    do Nzz=0,Nz-1
        do Nxx=0,Nx-1
            do zi=1,z_procs
                do xi=1,x_procs
                    phi((xi-1)+Nxx*x_procs,(zi-1)+Nzz*z_procs) = phiout(xi,zi,k)
                enddo
            enddo
            k = k + 1
        enddo
    enddo
    !========================液滴の数と測定範囲を調べる================================
    frag = 0
    fragformer = 0
    inum = 1
    knum = 0 !液滴の数
    do xi = 0, xmax
        do zi = 0, zmax
            if(phi(xi,zi) > phi_kaimen) then
                frag = frag + 1
            endif
        enddo
        if((frag > 0) .and. (fragformer == 0)) then
            xstart(inum) = xi - 1
            knum = knum + 1
        elseif((frag == 0) .and. (fragformer > 0)) then
            xend(inum) = xi + 1
            inum = inum + 1
        endif
        fragformer = frag
        frag = 0
    enddo

    DO inum = 1, knum
        !==============それぞれの液滴の重心座標を調べる=====================
        xg = 0.0d0
        zg = 0.0d0
        phi_tmp = 0.0d0
        do zi = 0, zmax
            do xi = xstart(inum), xend(inum)
                if(phi(xi,zi) >= phi_kaimen) then
                    xg = xg + dble(xi) 
                    zg = zg + dble(zi) 
                    phi_tmp = phi_tmp + 1.0d0
                endif
            enddo
        enddo
        xg = xg / phi_tmp 
        zg = zg / phi_tmp
        !=============長軸短軸半径を線形補間で調べる=======================
        Ltmp = 0.0d0
        Btmp = 2.0d0 * dble(xmax)
        do zi = 1, zmax-1
            do xi = xstart(inum)+1, xend(inum)-1
                !座標設定
                x0 = xi - 1
                x1 = xi
                x2 = xi + 1
                z0 = zi - 1
                z1 = zi
                z2 = zi + 1
                !右
                call senkei(x1,x2,z1,z1,xg,zg,Ltmp,Btmp,phi)
                !左
                call senkei(x1,x0,z1,z1,xg,zg,Ltmp,Btmp,phi)
                !上
                call senkei(x1,x1,z1,z2,xg,zg,Ltmp,Btmp,phi)
                !下
                call senkei(x1,x1,z1,z0,xg,zg,Ltmp,Btmp,phi)
                !右斜め上
                call senkei(x1,x2,z1,z2,xg,zg,Ltmp,Btmp,phi)
                !右斜め下
                call senkei(x1,x2,z1,z0,xg,zg,Ltmp,Btmp,phi)
                !左斜め上
                call senkei(x1,x0,z1,z2,xg,zg,Ltmp,Btmp,phi)
                !左斜め下
                call senkei(x1,x0,z1,z0,xg,zg,Ltmp,Btmp,phi)
            enddo
        enddo
        L1(inum) = Ltmp
        B1(inum) = Btmp
        !変形度
        Dd(inum) = abs(L1(inum)-B1(inum)) / (L1(inum)+B1(inum))
    ENDDO
    !===================出力================================================
    dropnumtmp = knum + 2
    write(dropnum, "(i1)") dropnumtmp
    if(step == 4000) then
        open(10,file="./L.d")
        open(11,file="./B.d")
        open(12,file="./D.d")
        write(10,"("//trim(dropnum)//"es16.8)") real(step),(real(step)-real(5000))*2.0d0*uw/real(zmax),L1(1:knum)
        write(11,"("//trim(dropnum)//"es16.8)") real(step),(real(step)-real(5000))*2.0d0*uw/real(zmax),B1(1:knum)
        write(12,"("//trim(dropnum)//"es16.8)") real(step),(real(step)-real(5000))*2.0d0*uw/real(zmax),Dd(1:knum)
        close(10)
        close(11)
        close(12)
    else
        open(10,file="./L.d",action="write",position="append")
        open(11,file="./B.d",action="write",position="append")
        open(12,file="./D.d",action="write",position="append")
        write(10,"("//trim(dropnum)//"es16.8)") real(step),(real(step)-real(5000))*2.0d0*uw/real(zmax),L1(1:knum)
        write(11,"("//trim(dropnum)//"es16.8)") real(step),(real(step)-real(5000))*2.0d0*uw/real(zmax),B1(1:knum)
        write(12,"("//trim(dropnum)//"es16.8)") real(step),(real(step)-real(5000))*2.0d0*uw/real(zmax),Dd(1:knum)
        close(10)
        close(11)
        close(12)
    endif
ENDDO
end program main