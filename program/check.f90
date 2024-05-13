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
    integer,parameter:: xmax = 511 !ｘ方向格子数
    integer,parameter:: ymax = 511 !ｙ方向格子数
    integer,parameter:: zmax = 511 !ｚ方向格子数
    integer,parameter:: Nx = 32 !ｘ方向の並列数
    integer,parameter:: Ny = 64 !ｘ方向の並列数
    integer,parameter:: Nz = 64 !ｚ方向の並列数
    integer,parameter:: x_procs = (xmax+1) / Nx
    integer,parameter:: y_procs = (ymax+1) / Ny
    integer,parameter:: z_procs = (zmax+1) / Nz
    integer,parameter:: procs_xz = Nx * Nz
    integer,parameter:: procs_xy = Nx * Ny
    integer i, j, k, xi, yi, zi, Nxx, Nyy, Nzz, step
    character(8) file_num, file_num2
    character :: filename*200
    integer set
    real(8),parameter :: ds = 1.0d0

    !ティレクトリー読み込み
    character(*),parameter :: datadir = "/data/lng/nakanog/taylor_512_re100000_large/fg/"
    !ディレクトリ作成
    character(*),parameter :: datadir2 = "/data/lng/nakanog/taylor_512_re100000_large_xy/fg/"

    real(8) f_xz(0:xmax,0:ymax,0:zmax)
    real(8) fout_xz(1:15,1:x_procs,1:ymax+1,1:z_procs,0:procs_xz-1)

    real(8) f_xy(0:xmax,0:ymax,0:zmax)
    real(8) fout_xy(1:15,1:x_procs,1:y_procs,1:zmax+1,0:procs_xy-1)

    real(8) gosa(0:xmax,0:ymax,0:zmax)

    integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
    integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
    integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)
    real(8):: cr(1:3, 1:15)  !粒子速度（実数）

    real(8) wa
    real(8),parameter:: pi = acos(-1.0d0) !円周率
    real(8) dummy

    call mk_dirs(datadir2)
    open(15,file="./time.d")
    write(15,*) "time"
    close(15)

    do i=1,15
        cr(1,i) = dble(cx(i))
        cr(2,i) = dble(cy(i))
        cr(3,i) = dble(cz(i))
    enddo

    !xz並列のデータ読み込み
    step = 10000000
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
    elseif((step > 9999999) .and. (step < 100000000)) then
        write(file_num2,"(i8)") step
    endif
    do i=0,procs_xz-1
        if((i >= 0) .and. (i < 10)) then
            write(file_num, "(i1)") i
        elseif((i > 9) .and. (i < 100)) then
            write(file_num,"(i2)") i
        elseif((i > 99) .and. (i < 1000)) then
            write(file_num,"(i3)") i
        elseif((i > 999) .and. (i < 10000)) then
            write(file_num,"(i4)") i
        endif

        set = 20
        open(set, file=datadir//"0_"//trim(file_num)//"_"//trim(file_num2)//"_fg.bin", form="unformatted")
        do zi = 1, z_procs
            do yi = 1, ymax+1
                do xi = 1, x_procs
                    do j = 1, 15
                        read(set)  fout_xz(j,xi,yi,zi,i)
                    enddo
                enddo
            enddo
        enddo
        close(set)

        open(15,file="./time.d", action="write",position="append")
        write(15,*) "Read step XZ ", step, i
        close(15)
    enddo
    !xz並列のグローバル変数
    k = 0
    do Nzz=0,Nz-1
        do Nxx=0,Nx-1
            do zi=1,z_procs
                do yi=1,ymax+1
                    do xi=1,x_procs
                        j = 1
                        f_xz((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*z_procs) = fout_xz(j,xi,yi,zi,k)
                    enddo
                enddo
            enddo
            k = k + 1
        enddo
    enddo

    step = 10000000
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
    elseif((step > 9999999) .and. (step < 100000000)) then
        write(file_num2,"(i8)") step
    endif
    !xy並列のデータ読み込み
    do i=0,procs_xy-1
        if((i >= 0) .and. (i < 10)) then
            write(file_num, "(i1)") i
        elseif((i > 9) .and. (i < 100)) then
            write(file_num,"(i2)") i
        elseif((i > 99) .and. (i < 1000)) then
            write(file_num,"(i3)") i
        elseif((i > 999) .and. (i < 10000)) then
            write(file_num,"(i4)") i
        endif

        set = 21
        open(set, file=datadir2//"0_"//trim(file_num)//"_"//trim(file_num2)//".bin", form="unformatted")
        do zi = 1, zmax+1
            do yi = 1, y_procs
                do xi = 1, x_procs
                    do j = 1, 15
                        read(set) fout_xy(j,xi,yi,zi,i)
                    enddo
                enddo
            enddo
        enddo
        close(set)

        open(15,file="./time.d", action="write",position="append")
        write(15,*) "Read step XY ", step, i
        close(15)
    enddo
    !xy並列グローバル変数
    k = 0
    do Nxx=0,Nx-1
        do Nyy=0,Ny-1
            do zi=1,zmax+1
                do yi=1,y_procs
                    do xi=1,x_procs
                        j = 1
                        f_xy((xi-1)+Nxx*x_procs,(yi-1)+Nyy*y_procs,(zi-1)) = fout_xy(j,xi,yi,zi,k)
                    enddo
                enddo
            enddo
            k = k + 1
        enddo
    enddo

    wa = 0.0d0
    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                do j=1,15
                    gosa(xi,yi,zi) = abs(f_xz(xi,yi,zi)-f_xy(xi,yi,zi)) / abs(f_xz(xi,yi,zi))
                    wa = wa + gosa(xi,yi,zi)
                enddo
            enddo
        enddo
    enddo

    open(10,file="err.d")
    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                write(10,*) xi,yi,zi,gosa(xi,yi,zi)
            enddo
        enddo
    enddo
    write(10,*) wa
    close(10)

end program main