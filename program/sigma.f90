module globals
    integer,parameter:: xmax = 63 !ｘ方向格子数
    integer,parameter:: ymax = 63 !ｙ方向格子数
    integer,parameter:: zmax = 63 !ｚ方向格子数
    integer,parameter:: Nx = 8 !ｘ方向の並列数
    integer,parameter:: Ny = 8 !ｚ方向の並列数
    integer,parameter:: x_procs = (xmax+1) / Nx
    integer,parameter:: y_procs = (zmax+1) / Ny
    integer,parameter:: new_procs = Nx * Ny
    ! real(8),parameter:: kappag = (2.156d-4) / (2.59d-3)
    real(8),parameter:: kappag = (2.156d-4/1.7039d0)**(1.0d0/0.9991d0)  !界面張力を決めるパラメータ
    real(8) kappagtmp
    ! real(8),parameter:: phi1 = 2.638d-1 !連続相のオーダーパラメータ
    ! real(8),parameter:: phi2 = 4.031d-1 !分散相のオーダーパラメータ
    real(8),parameter:: phi1 = 2.211d0 !連続相のオーダーパラメータ
    real(8),parameter:: phi2 = 4.895d0 !分散相のオーダーパラメータ
    real(8) phi_kaimen 
    integer,parameter:: step = 20000
contains

    subroutine par(cx,cy,cz,cr)
        integer, intent(in):: cx(15),cy(15),cz(15)
        real(8), intent(out):: cr(1:3,1:15)
        do i=1,15
            cr(1,i) = dble(cx(i))
            cr(2,i) = dble(cy(i))
            cr(3,i) = dble(cz(i))
        enddo
    end subroutine par

    subroutine senkei(xx1,xx2,zz1,zz2,xxg,zzg,L1,B1,phi)
        integer, intent(in) :: xx1,xx2,zz1,zz2
        real(8), intent(in) :: xxg,zzg
        real(8), intent(inout) :: L1, B1
        real(8), intent(in) :: phi(-1:xmax+1, -1:ymax+1, -1:zmax+1)
        real(8) distance
        real(8) alpha, l
        real(8) xj, zj
        
        l = sqrt((dble(xx2)-dble(xx1))**2 +(dble(zz2)-dble(zz1))**2)
        alpha = l * (phi_kaimen - phi(xx1,ymax/2,zz1)) / (phi(xx2,ymax/2,zz2) - phi(xx1,(ymax+1)/2,zz1))
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

    !粒子速度（整数）
    integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
    integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
    integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)
    real(8):: cr(1:3, 1:15)  !粒子速度（実数）
    real(8) sigma
    integer i, k, j, xi, yi, zi, Nxx, Nyy, alpha
    character :: filename*200
    character(8) file_num, file_num2

    real(8),parameter:: pi = acos(-1.0d0)
    !ティレクトリー読み込み
    character(*),parameter :: datadir = "/data/sht/nakanog/droplet_test2/"

    real(8) phi(-1:xmax+1,-1:ymax+1,-1:zmax+1)
    real(8) grad_phi(1:3,-1:xmax+1,-1:ymax+1,-1:zmax+1)
    real(8) phiout(1:x_procs,1:y_procs,1:zmax+1,0:new_procs-1)
    real(8) dummy1, dummy2, dummy3, dummy4
    real(8) Vol
    real(8) p(-1:xmax+1,-1:ymax+1,-1:zmax+1)
    real(8) pout(1:x_procs,1:y_procs,1:zmax+1,0:new_procs-1)
    real(8) p_true(-1:xmax+1,-1:ymax+1,-1:zmax+1)

    call par(cx,cy,cz,cr)

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

        if((step > 99) .and. (step < 1000)) then
            write(file_num2, "(i3)") step
        elseif((step > 999) .and. (step < 10000)) then
            write(file_num2,"(i4)") step
        elseif((step > 9999) .and. (step < 100000)) then
            write(file_num2,"(i5)") step
        elseif((step > 99999) .and. (step < 1000000)) then
            write(file_num2,"(i6)") step
        endif

        open(9, file=datadir//"0_"//trim(file_num)//"_"//trim(file_num2)//".bin", form="unformatted")
        do zi = 1, zmax+1
            do yi = 1, y_procs
                do xi = 1, x_procs
                    ! read(9) phiout(xi,yi,zi,i), dummy1, dummy2, dummy3, dummy4
                    read(9) phiout(xi,yi,zi,i), dummy1, dummy2, dummy3, pout(xi,yi,zi,i)
                enddo
            enddo
        enddo
        close(9)
    enddo
    k = 0
    do Nxx=0,Nx-1
        do Nyy=0,Ny-1
            do zi=1,zmax+1
                do yi=1,y_procs
                    do xi=1,x_procs
                        phi((xi-1)+Nxx*x_procs,(yi-1)+Nyy*y_procs,zi-1) = phiout(xi,yi,zi,k)
                        p((xi-1)+Nxx*x_procs,(yi-1)+Nyy*y_procs,zi-1) = pout(xi,yi,zi,k)
                    enddo
                enddo
            enddo
            k = k + 1
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
    !grad_phiの計算
    Vol = 0.0d0
    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                do alpha=1,3
                    grad_phi(alpha,xi,yi,zi) = 0.0d0
                    do i = 2,15
                        grad_phi(alpha,xi,yi,zi) = grad_phi(alpha,xi,yi,zi) + cr(alpha,i)*phi(xi+cx(i),yi+cy(i),zi+cz(i))
                    enddo
                    grad_phi(alpha,xi,yi,zi) = grad_phi(alpha,xi,yi,zi) / (10.0d0)
                enddo
                if((phi1+phi2)/2.0d0 <= phi(xi,yi,zi)) then
                    Vol = Vol + 1.0d0
                endif
            enddo
        enddo
    enddo
    !表面張力の計算
    ! yi = (ymax+1) / 2
    ! zi = (zmax+1) / 2
    ! sigma = 0.0d0
    ! do xi=0,(xmax+1)/2
    !     sigma = sigma + grad_phi(1,xi,yi,zi)**2
    ! enddo
    ! sigma = sigma * kappag

    ! open(10,file="./sigma_kg.d")
    ! write(10,*) sigma, kappag, Vol/17157.0d0
    ! close(10)

    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                p_true(xi,yi,zi) = p(xi,yi,zi) + 2.0d0/3.0d0*kappag*(grad_phi(1,xi,yi,zi)**2+grad_phi(2,xi,yi,zi)**2+grad_phi(3,xi,yi,zi)**2)
            enddo
        enddo
    enddo

    open(11,file="./p_3dmap_2.d")
    zi = (zmax+1) / 2
    do yi=0,ymax
        do xi=0,xmax
            write(11,*) xi, yi, p(xi,yi,zi), p_true(xi,yi,zi)
        enddo
        write(11,*)
    enddo
    close(11)

    open(12,file="./p_xmap_2.d")
    zi = (zmax+1) / 2
    yi = (ymax+1) / 2
    do xi=0,xmax
        write(12,*) xi, p(xi,yi,zi), p_true(xi,yi,zi)
    enddo
    close(12)

    open(13,file="./phi_2.d")
    zi = (zmax+1) / 2
    yi = (ymax+1) / 2
    do xi=0,xmax
        write(13,*) xi, phi(xi,yi,zi)
    enddo
    close(13)

end program main