module globals
contains
    subroutine mk_dirs(outdir)
        implicit none
        character(len=*), intent(in) :: outdir
        character(len=256) command
        write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
        call system(command)
    end subroutine mk_dirs
end module globals

program main
use globals
    implicit none
    integer,parameter:: xmax = 63 !ｘ方向格子数
    integer,parameter:: ymax = 63 !ｙ方向格子数
    integer,parameter:: zmax = 63 !ｚ方向格子数
    integer,parameter:: Nx = 4 !ｘ方向の並列数
    integer,parameter:: Ny = 4 !ｘ方向の並列数
    integer,parameter:: Nz = 4 !ｚ方向の並列数
    integer,parameter:: x_procs = (xmax+1) / Nx
    integer,parameter:: y_procs = (ymax+1) / Ny
    integer,parameter:: z_procs = (zmax+1) / Nz
    integer,parameter:: new_procs = Nx * Nz
    integer i, j, k, xi, yi, zi, Nxx, Nzz, step, sumnum, beta, alpha, Nyy
    real x, y, z, s
    character(8) file_num, file_num2
    character :: filename*200
    integer set
    real(8),parameter :: ds = 1.0d0

    !ティレクトリー読み込み
    character(*),parameter :: datadir = "/data/sht/nakanog/debug/fg/"
    !ディレクトリ作成
    character(*),parameter :: datadir2 = "/data/sht/nakanog/debug2/fg/"

    real(8) f_fft(1:15,0:xmax,0:ymax,0:zmax)
    real(8) fout_fft(1:15,1:x_procs,1:y_procs,1:zmax+1,0:new_procs-1)

    real(8) f_ori(1:15,0:xmax,0:ymax,0:zmax)
    real(8) fout_ori(1:15,1:x_procs,1:ymax+1,1:z_procs,0:new_procs-1)

    integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
    integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
    integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)
    real(8):: cr(1:3, 1:15)  !粒子速度（実数）

    real(8) wa
    real(8),parameter:: pi = acos(-1.0d0) !円周率



    do i=1,15
        cr(1,i) = dble(cx(i))
        cr(2,i) = dble(cy(i))
        cr(3,i) = dble(cz(i))
    enddo

    step = 100
    if((step > 99) .and. (step < 1000)) then
        write(file_num2, "(i3)") step
    elseif((step > 999) .and. (step < 10000)) then
        write(file_num2,"(i4)") step
    elseif((step > 9999) .and. (step < 100000)) then
        write(file_num2,"(i5)") step
    elseif((step > 99999) .and. (step < 1000000)) then
        write(file_num2,"(i6)") step
    elseif((step > 999999) .and. (step < 10000000)) then
        write(file_num2,"(i7)") step
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
        open(set, file=datadir//"0_"//trim(file_num)//"_"//trim(file_num2)//"_fg.bin", form="unformatted")
        do zi = 1, zmax+1
            do yi = 1, y_procs
                do xi = 1, x_procs
                    do j = 1, 15
                        read(set)  fout_fft(j,xi,yi,zi,i)
                    enddo
                enddo
            enddo
        enddo
        close(set)
    enddo

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
        open(set, file=datadir2//"0_"//trim(file_num)//"_"//trim(file_num2)//"_fg.bin", form="unformatted")
        do zi = 1, z_procs
            do yi = 1, ymax+1
                do xi = 1, x_procs
                    do j = 1, 15
                        read(set)  fout_ori(j,xi,yi,zi,i)
                    enddo
                enddo
            enddo
        enddo
        close(set)
    enddo

    k = 0
    do Nzz=0,Nz-1
        do Nxx=0,Nx-1
            do zi=1,z_procs
                do yi=1,ymax+1
                    do xi=1,x_procs
                        do j=1,15
                            f_ori(j,(xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*z_procs) = fout_ori(j,xi,yi,zi,k)
                        enddo
                    enddo
                enddo
            enddo
            k = k + 1
        enddo
    enddo

    k = 0
    do Nyy=0,Ny-1
        do Nxx=0,Nx-1
            do zi=1,zmax+1
                do yi=1,y_procs
                    do xi=1,x_procs
                        do j=1,15
                            f_fft(j,(xi-1)+Nxx*x_procs,(yi-1)+Nyy*y_procs,(zi-1)) = fout_fft(j,xi,yi,zi,k)
                        enddo
                    enddo
                enddo
            enddo
            k = k + 1
        enddo
    enddo

    wa = 0.0d0
    open(100,file="./error.d")
    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                do j=1,15
                    wa = wa + dble(f_ori(j,xi,yi,zi)-f_fft(j,xi,yi,zi)) / dble(f_ori(j,xi,yi,zi))
                    write(100,*) dble(f_ori(j,xi,yi,zi)-f_fft(j,xi,yi,zi)) / dble(f_ori(j,xi,yi,zi))
                enddo
            enddo
        enddo
    enddo
    write(100,*) "wa = ", wa
    close(100)

end program main

        ! do zi=0,zmax
        !     do yi=0,ymax
        !         do xi=0,xmax
        !             do beta=1,3
        !                 grad_u(1,beta,xi,yi,zi) = 0.0d0
        !                 grad_u(2,beta,xi,yi,zi) = 0.0d0
        !                 grad_u(3,beta,xi,yi,zi) = 0.0d0
        !                 do i=2,15
        !                     grad_u(1,beta,xi,yi,zi) = grad_u(1,beta,xi,yi,zi)+cr(beta,i)*u1(xi+cx(i),yi+cy(i),zi+cz(i))
        !                     grad_u(2,beta,xi,yi,zi) = grad_u(2,beta,xi,yi,zi)+cr(beta,i)*u2(xi+cx(i),yi+cy(i),zi+cz(i))
        !                     grad_u(3,beta,xi,yi,zi) = grad_u(3,beta,xi,yi,zi)+cr(beta,i)*u3(xi+cx(i),yi+cy(i),zi+cz(i))
        !                 enddo
        !                 grad_u(1,beta,xi,yi,zi) = grad_u(1,beta,xi,yi,zi)/(10.0d0*ds)
        !                 grad_u(2,beta,xi,yi,zi) = grad_u(2,beta,xi,yi,zi)/(10.0d0*ds)
        !                 grad_u(3,beta,xi,yi,zi) = grad_u(3,beta,xi,yi,zi)/(10.0d0*ds)
        !             enddo
        !         enddo
        !     enddo
        ! enddo

        ! do zi=0,zmax
        !     do yi=0,ymax
        !         do xi=0,xmax
        !             do beta=1,3
        !                 do alpha=1,3
        !                     hizumi(alpha,beta,xi,yi,zi) = grad_u(alpha,beta,xi,yi,zi) + grad_u(beta,alpha,xi,yi,zi)
        !                     uzudo(alpha,beta,xi,yi,zi) = grad_u(alpha,beta,xi,yi,zi) - grad_u(beta,alpha,xi,yi,zi)
        !                 enddo
        !             enddo
        !         enddo
        !     enddo
        ! enddo

        ! do zi=0,zmax
        !     do yi=0,ymax
        !         do xi=0,xmax
        !             qti(xi,yi,zi) = (uzudo(1,1,xi,yi,zi)*uzudo(2,2,xi,yi,zi)*uzudo(3,3,xi,yi,zi) &
        !             +uzudo(1,2,xi,yi,zi)*uzudo(2,3,xi,yi,zi)*uzudo(3,1,xi,yi,zi) &
        !             +uzudo(1,3,xi,yi,zi)*uzudo(2,1,xi,yi,zi)*uzudo(3,2,xi,yi,zi) &
        !             -uzudo(1,3,xi,yi,zi)*uzudo(2,2,xi,yi,zi)*uzudo(3,1,xi,yi,zi) &
        !             -uzudo(1,1,xi,yi,zi)*uzudo(2,3,xi,yi,zi)*uzudo(3,2,xi,yi,zi) &
        !             -uzudo(1,2,xi,yi,zi)*uzudo(2,1,xi,yi,zi)*uzudo(3,3,xi,yi,zi))**2 &
        !             - (hizumi(1,1,xi,yi,zi)*hizumi(2,2,xi,yi,zi)*hizumi(3,3,xi,yi,zi) &
        !             +hizumi(1,2,xi,yi,zi)*hizumi(2,3,xi,yi,zi)*hizumi(3,1,xi,yi,zi) &
        !             +hizumi(1,3,xi,yi,zi)*hizumi(2,1,xi,yi,zi)*hizumi(3,2,xi,yi,zi) &
        !             -hizumi(1,3,xi,yi,zi)*hizumi(2,2,xi,yi,zi)*hizumi(3,1,xi,yi,zi) &
        !             -hizumi(1,1,xi,yi,zi)*hizumi(2,3,xi,yi,zi)*hizumi(3,2,xi,yi,zi) &
        !             -hizumi(1,2,xi,yi,zi)*hizumi(2,1,xi,yi,zi)*hizumi(3,3,xi,yi,zi))**2

        !             qti(xi,yi,zi) = 1.0d0/8.0d0 * qti(xi,yi,zi)
        !         enddo
        !     enddo
        ! enddo
