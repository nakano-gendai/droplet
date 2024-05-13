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

    real(8),parameter:: ds = 1.0d0
    integer,parameter:: xmax = 255 !ｘ方向格子数
    integer,parameter:: ymax = 255 !ｙ方向格子数
    integer,parameter:: zmax = 255 !ｚ方向格子数
    integer,parameter:: Nx = 8 !ｘ方向の並列数
    integer,parameter:: Nz = 16 !ｚ方向の並列数
    integer,parameter:: x_procs = (xmax+1) / Nx
    integer,parameter:: z_procs = (zmax+1) / Nz
    integer,parameter:: new_procs = Nx * Nz

    real(8),parameter:: delta = 10.0d0*ds

    integer i, j, k, xi, yi, zi, Nxx, Nzz, step, beta

    real(8) u1(-1:xmax+1, -1:ymax+1, -1:zmax+1), u2(-1:xmax+1, -1:ymax+1, -1:zmax+1), u3(-1:xmax+1, -1:ymax+1, -1:zmax+1), p(-1:xmax+1, -1:ymax+1, -1:zmax+1)
    real(8) u1_f(0:xmax, 0:ymax, 0:zmax), u2_f(0:xmax, 0:ymax, 0:zmax), u3_f(0:xmax, 0:ymax, 0:zmax), p_f(0:xmax, 0:ymax, 0:zmax)
    real(8) u1out(1:x_procs, 1:ymax+1, 1:z_procs, 0:new_procs-1), u2out(1:x_procs, 1:ymax+1, 1:z_procs, 0:new_procs-1), u3out(1:x_procs, 1:ymax+1, 1:z_procs, 0:new_procs-1), pout(1:x_procs, 1:ymax+1, 1:z_procs, 0:new_procs-1)
    real(8) dummy1, dummy2

    !粒子速度
    integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
    integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
    integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)
    real(8):: cr(1:3, 1:15)

    real(8),parameter:: pi = acos(-1.0d0) !円周率

    !ティレクトリー読み込み
    ! character(*),parameter :: datadir = "/data/sht/nakanog/taylor_re12000_drop14/"
    character(*),parameter :: datadir = "/data/sht/nakanog/taylor_re12000_ran/"
    !ディレクトリ作成
    character(*),parameter :: datadir2 = "/data/sht/nakanog/filter_test/"

    character(8) file_num, file_num2
    character :: filename*200
    integer set

    call mk_dirs(datadir2)
    open(10,file="./log.d")
    close(10)
!===============データの読み込み=========================================
    step = 1200000
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
        open(set, file=datadir//"0_"//trim(file_num)//"_"//trim(file_num2)//".bin", form="unformatted")
        do zi = 1, z_procs
            do yi = 1, ymax
                do xi = 1, x_procs
                    read(set) u1out(xi,yi,zi,i), u2out(xi,yi,zi,i), u3out(xi,yi,zi,i), pout(xi,yi,zi,i), dummy1
                enddo
            enddo
        enddo
        close(set)
    enddo
    k = 0
    do Nzz=0,Nz-1
        do Nxx=0,Nx-1
            do zi=1,z_procs
                do yi=1,ymax
                    do xi=1,x_procs
                        u1((xi)+Nxx*x_procs,yi,(zi)+Nzz*z_procs) = u1out(xi,yi,zi,k)
                        u2((xi)+Nxx*x_procs,yi,(zi)+Nzz*z_procs) = u2out(xi,yi,zi,k)
                        u3((xi)+Nxx*x_procs,yi,(zi)+Nzz*z_procs) = u3out(xi,yi,zi,k)
                        p((xi)+Nxx*x_procs,yi,(zi)+Nzz*z_procs) = pout(xi,yi,zi,k)
                    enddo
                enddo
            enddo
            k = k + 1
        enddo
    enddo
    open(10,file="./log.d",action="write",position="append")
    write(10,*) "Input is OK!"
    close(10)
!============================周期境界条件=============================================
    do zi=0,zmax
        do yi=-1,ymax+1
            do xi=0,xmax
                u1(xi,-1,zi) = u1(xi,ymax,zi)
                u1(xi,ymax+1,zi) = u1(xi,0,zi)
                u1(-1,yi,zi) = u1(xmax,yi,zi)
                u1(xmax+1,yi,zi) = u1(0,yi,zi)

                u2(xi,-1,zi) = u2(xi,ymax,zi)
                u2(xi,ymax+1,zi) = u2(xi,0,zi)
                u2(-1,yi,zi) = u2(xmax,yi,zi)
                u2(xmax+1,yi,zi) = u2(0,yi,zi)

                u3(xi,-1,zi) = u3(xi,ymax,zi)
                u3(xi,ymax+1,zi) = u3(xi,0,zi)
                u3(-1,yi,zi) = u3(xmax,yi,zi)
                u3(xmax+1,yi,zi) = u3(0,yi,zi)

                p(xi,-1,zi) = p(xi,ymax,zi)
                p(xi,ymax+1,zi) = p(xi,0,zi)
                p(-1,yi,zi) = p(xmax,yi,zi)
                p(xmax+1,yi,zi) = p(0,yi,zi)
            enddo
        enddo
    enddo
    do yi=-1,ymax+1
        do xi=-1,xmax+1
            u1(xi,yi,-1) = u1(xi,yi,zmax)
            u1(xi,yi,zmax+1) = u1(xi,yi,0)

            u2(xi,yi,-1) = u2(xi,yi,zmax)
            u2(xi,yi,zmax+1) = u2(xi,yi,0)

            u3(xi,yi,-1) = u3(xi,yi,zmax)
            u3(xi,yi,zmax+1) = u3(xi,yi,0)

            p(xi,yi,-1) = p(xi,yi,zmax)
            p(xi,yi,zmax+1) = p(xi,yi,0)
        enddo
    enddo
    open(10,file="./log.d",action="write",position="append")
    write(10,*) "B.C. is OK!"
    close(10)
!=============================フィルター==============================================
    !X方向
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                u1_f(xi,yi,zi) = u1(xi,yi,zi) &
                                    + delta * delta / 24.0d0 &
                                    * ( u1(xi-1,yi,zi) - 2.0d0 * u1(xi,yi,zi) + u1(xi+1,yi,zi) ) &
                                    / ds**2.0d0
                u2_f(xi,yi,zi) = u2(xi,yi,zi) &
                                    + delta * delta / 24.0d0 &
                                    * ( u2(xi-1,yi,zi) - 2.0d0 * u2(xi,yi,zi) + u2(xi+1,yi,zi) ) &
                                    / ds**2.0d0
                u3_f(xi,yi,zi) = u3(xi,yi,zi) &
                                    + delta * delta / 24.0d0 &
                                    * ( u3(xi-1,yi,zi) - 2.0d0 * u3(xi,yi,zi) + u3(xi+1,yi,zi) ) &
                                    / ds**2.0d0
            enddo
        enddo
    enddo
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                u1(xi,yi,zi) = u1_f(xi,yi,zi)
                u2(xi,yi,zi) = u2_f(xi,yi,zi)
                u3(xi,yi,zi) = u3_f(xi,yi,zi)
            enddo
        enddo
    enddo
    open(10,file="./log.d",action="write",position="append")
    write(10,*) "X filter is OK!"
    close(10)
    !Y方向
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                u1_f(xi,yi,zi) = u1(xi,yi,zi) &
                                    + delta * delta / 24.0d0 &
                                    * ( u1(xi,yi-1,zi) - 2.0d0 * u1(xi,yi,zi) + u1(xi,yi+1,zi) ) &
                                    / ds**2.0d0
                u2_f(xi,yi,zi) = u2(xi,yi,zi) &
                                    + delta * delta / 24.0d0 &
                                    * ( u2(xi,yi-1,zi) - 2.0d0 * u2(xi,yi,zi) + u2(xi,yi+1,zi) ) &
                                    / ds**2.0d0
                u3_f(xi,yi,zi) = u3(xi,yi,zi) &
                                    + delta * delta / 24.0d0 &
                                    * ( u3(xi,yi-1,zi) - 2.0d0 * u3(xi,yi,zi) + u3(xi,yi+1,zi) ) &
                                    / ds**2.0d0
            enddo
        enddo
    enddo
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                u1(xi,yi,zi) = u1_f(xi,yi,zi)
                u2(xi,yi,zi) = u2_f(xi,yi,zi)
                u3(xi,yi,zi) = u3_f(xi,yi,zi)
            enddo
        enddo
    enddo
    open(10,file="./log.d",action="write",position="append")
    write(10,*) "Y filter is OK!"
    close(10)
    !Z方向
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                u1_f(xi,yi,zi) = u1(xi,yi,zi) &
                                    + delta * delta / 24.0d0 &
                                    * ( u1(xi,yi,zi-1) - 2.0d0 * u1(xi,yi,zi) + u1(xi,yi,zi+1) ) &
                                    / ds**2.0d0
                u2_f(xi,yi,zi) = u2(xi,yi,zi) &
                                    + delta * delta / 24.0d0 &
                                    * ( u2(xi,yi,zi-1) - 2.0d0 * u2(xi,yi,zi) + u2(xi,yi,zi+1) ) &
                                    / ds**2.0d0
                u3_f(xi,yi,zi) = u3(xi,yi,zi) &
                                    + delta * delta / 24.0d0 &
                                    * ( u3(xi,yi,zi-1) - 2.0d0 * u3(xi,yi,zi) + u3(xi,yi,zi+1) ) &
                                    / ds**2.0d0
            enddo
        enddo
    enddo
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                u1(xi,yi,zi) = u1_f(xi,yi,zi)
                u2(xi,yi,zi) = u2_f(xi,yi,zi)
                u3(xi,yi,zi) = u3_f(xi,yi,zi)
            enddo
        enddo
    enddo
    open(10,file="./log.d",action="write",position="append")
    write(10,*) "Z filter is OK!"
    close(10)
    
    open(11, file="filter.bin", form='unformatted',status='replace')
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                write(11) u1(xi,yi,zi), u2(xi,yi,zi), u3(xi,yi,zi)
            enddo
        enddo
    enddo
    close(11)
end program main