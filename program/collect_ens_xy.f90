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
    integer,parameter:: Ny = 64 !ｚ方向の並列数
    integer,parameter:: x_procs = (xmax+1) / Nx
    integer,parameter:: y_procs = (ymax+1) / Ny
    integer,parameter:: new_procs = Nx * Ny
    integer i, j, k, xi, yi, zi, Nxx, Nyy, step, beta
    real(8) x, y, z, dummy
    character :: filename*200
    character(8) file_num, file_num2
    !ティレクトリー読み込み
    character(*),parameter :: datadir = "/data/sht/nakanog/taylor_512_drop_movie/medium/"
    !ディレクトリ作成
    character(*),parameter :: datadir2 = "/data/sht/nakanog/taylor_512_drop_movie/medium/collect/"
    real(8) phi(0:xmax,0:ymax,0:zmax)
    real(8) phiout(1:x_procs,1:y_procs,1:zmax+1,0:new_procs-1)
    real(8) s, ave
    real(8) a, b, c, d

    real(8) time1, time2
    integer set

    call mk_dirs(datadir2)

!============phiをまとめて出力するプログラム=============================================
    do step=5000, 1000000, 5000
        ! step = 400000
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

            set = 20
            open(set, file=datadir//"0_"//trim(file_num)//"_"//trim(file_num2)//".bin", form="unformatted")
            do zi = 1, zmax+1
                do yi = 1, y_procs
                    do xi = 1, x_procs
                        read(set) phiout(xi,yi,zi,i)
                    enddo
                enddo
            enddo
            close(set)
        enddo
        k = 0
        do Nxx=0,Nx-1
            do Nyy=0,Ny-1
                do zi=1,zmax+1
                    do yi=1,y_procs
                        do xi=1,x_procs
                            phi((xi-1)+Nxx*x_procs,(yi-1)+Nyy*y_procs,zi-1) = phiout(xi,yi,zi,k)
                        enddo
                    enddo
                enddo
                k = k + 1
            enddo
        enddo
        
        write(filename,*) step !i->filename 変換
        filename=datadir2//trim(adjustl(filename))//'.bin' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
        print *, filename !表示してみる
        open(11, file=filename, form="unformatted", status='replace') 
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    write(11) phi(xi,yi,zi)
                enddo
            enddo
        enddo
        close(11)

        ! write(filename,*) step !i->filename 変換
        ! filename=datadir2//trim(adjustl(filename))//'.d' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
        ! print *, filename !表示してみる
        ! open(11, file=filename, form="formatted", status='replace') 
        ! yi = (ymax+1)/2
        ! do zi=0,zmax
        !     do xi=0,xmax
        !         write(11,*) xi, zi, phi(xi,yi,zi)
        !     enddo
        !     write(11,*)
        ! enddo
        ! close(11)

    enddo
end program main