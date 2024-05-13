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
    integer,parameter:: xmax = 255 !ｘ方向格子数
    integer,parameter:: ymax = 255 !ｙ方向格子数
    integer,parameter:: zmax = 255 !ｚ方向格子数
    integer,parameter:: Nx = 8 !ｘ方向の並列数
    integer,parameter:: Ny = 16 !ｘ方向の並列数
    integer,parameter:: Nz = 16 !ｚ方向の並列数
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
    character(*),parameter :: datadir = "/data/sht/nakanog/DNS_turbulence_256_re80000/fg/"
    !ディレクトリ作成
    character(*),parameter :: datadir2 = "/data/sht/nakanog/DNS_turbulence_256_re80000/xy/fg/"

    real(8) f_xz(1:15,0:xmax,0:ymax,0:zmax)
    real(8) fout_xz(1:15,1:x_procs,1:ymax+1,1:z_procs,0:procs_xz-1)

    real(8) f_xy(1:15,0:xmax,0:ymax,0:zmax)
    real(8) fout_xy(1:15,1:x_procs,1:y_procs,1:zmax+1,0:procs_xy-1)

    real(8) err(1:15,0:xmax,0:ymax,0:zmax)

    integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
    integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
    integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)
    real(8):: cr(1:3, 1:15)  !粒子速度（実数）

    real(8) wa
    real(8),parameter:: pi = acos(-1.0d0) !円周率

    call mk_dirs(datadir2)

    do i=1,15
        cr(1,i) = dble(cx(i))
        cr(2,i) = dble(cy(i))
        cr(3,i) = dble(cz(i))
    enddo

    !xz並列のデータ読み込み
    step = 2400000
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

        open(10, file=datadir//"0_"//trim(file_num)//"_"//trim(file_num2)//"_fg.bin", form="unformatted")
        do zi = 1, z_procs
            do yi = 1, ymax+1
                do xi = 1, x_procs
                    do j = 1, 15
                        read(10)  fout_xz(j,xi,yi,zi,i)
                    enddo
                enddo
            enddo
        enddo
        close(10)
    enddo

    !xz並列のグローバル変数
    k = 0
    do Nzz=0,Nz-1
        do Nxx=0,Nx-1
            do zi=1,z_procs
                do yi=1,ymax+1
                    do xi=1,x_procs
                        do j=1,15
                            f_xz(j,(xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*z_procs) = fout_xz(j,xi,yi,zi,k)
                        enddo
                    enddo
                enddo
            enddo
            k = k + 1
        enddo
    enddo

    fout_xy(:,:,:,:,:) = 0.0d0
    k = 0
    do Nyy=0,Ny-1
        do Nxx=0,Nx-1
            do zi=1,zmax+1
                do yi=1,y_procs
                    do xi=1,x_procs
                        do j=1,15
                            fout_xy(j,xi,yi,zi,k) = f_xz(j,(xi-1)+Nyy*x_procs,(yi-1)+Nxx*y_procs,(zi-1))
                        enddo
                    enddo
                enddo
            enddo
            k = k + 1
        enddo
    enddo

    !出力
    do i = 0, procs_xy-1
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
        open(set, file=datadir2//'0_'//trim(adjustl(file_num))//'_2400000_fg.bin', form="unformatted")
        do zi=1,zmax+1
            do yi=1,y_procs
                do xi=1,x_procs
                    do j=1,15
                        write(set) fout_xy(j,xi,yi,zi,i)
                    enddo
                enddo
            enddo
        enddo
        close(set)
    enddo
end program main