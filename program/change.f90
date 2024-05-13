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
    integer,parameter:: Nz = 64 !ｚ方向の並列数
    integer,parameter:: x_procs = (xmax+1) / Nx
    integer,parameter:: z_procs = (zmax+1) / Nz
    integer,parameter:: new_procs = Nx * Nz
    integer i, j, k, xi, yi, zi, Nxx, Nzz, step, beta, eq
    real(8) x, y, z, dummy
    character :: filename*200
    character(8) file_num, file_num2
    !ティレクトリー読み込み
    character(*),parameter :: datadir = "/data/lng/nakanog/taylor_512_re100000_large/fg/"
    !ディレクトリ作成
    character(*),parameter :: datadir2 = "/data/lng/nakanog/taylor_512_re100000_large_xy/fg/"

    real(8) phi(1:15,0:xmax,0:ymax,0:zmax)
    real(8) phiout(1:15,1:x_procs,1:ymax+1,1:z_procs,0:new_procs-1)
    real(8) phiout2(1:15,1:x_procs,1:z_procs,1:zmax+1,0:new_procs-1)
    ! real(8) u1(0:xmax,0:ymax,0:zmax), u2(0:xmax,0:ymax,0:zmax), u3(0:xmax,0:ymax,0:zmax)
    ! real(8) u1out(1:x_procs,1:ymax+1,1:z_procs,0:new_procs-1), u2out(1:x_procs,1:ymax+1,1:z_procs,0:new_procs-1), u3out(1:x_procs,1:ymax+1,1:z_procs,0:new_procs-1)
    real(8) s, ave
    real(8) a, b, c, d

    real(8) time1, time2

    call mk_dirs(datadir2)

    open(10,file="./time.d")
    write(10,*) "time"
    close(10)

!============phiをまとめて出力するプログラム=============================================
    ! DO step=5000, 100000000, 10000
        step = 10000000
        !データの読み込み(x-z並列)
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

            open(20, file=datadir//"0_"//trim(file_num)//"_"//trim(file_num2)//"_fg.bin", form="unformatted")
            do zi = 1, z_procs
                do yi = 1, ymax+1
                    do xi = 1, x_procs
                        ! read(20) u1out(xi,yi,zi,i), u2out(xi,yi,zi,i), u3out(xi,yi,zi,i), b
                        do eq = 1, 15
                            read(20) phiout(eq,xi,yi,zi,i)
                        enddo
                    enddo
                enddo
            enddo
            close(20)

            open(10,file="./time.d", action="write",position="append")
            write(10,*) "Read step", step, i
            close(10)
        enddo

        !並列化した変数を1つの変数にまとめる(x-z並列)
        k = 0
        do Nzz=0,Nz-1
            do Nxx=0,Nx-1
                do zi=1,z_procs
                    do yi=1,ymax+1
                        do xi=1,x_procs
                            ! u1((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*z_procs) = u1out(xi,yi,zi,k)
                            ! u2((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*z_procs) = u2out(xi,yi,zi,k)
                            ! u3((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*z_procs) = u3out(xi,yi,zi,k)
                            do eq=1,15
                                phi(eq,(xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*z_procs) = phiout(eq,xi,yi,zi,k)
                            enddo
                        enddo
                    enddo
                enddo
                k = k + 1
            enddo
        enddo

        ! u1out(:,:,:,:) = 0.0d0
        ! u2out(:,:,:,:) = 0.0d0
        ! u3out(:,:,:,:) = 0.0d0
        ! phiout(:,:,:,:,:) = 0.0d0

        !並列化の方向を変換する(x-z→x-y)
        k = 0
        do Nxx=0,Nx-1
            do Nzz=0,Nz-1
                do zi=1,zmax+1
                    do yi=1,z_procs
                        do xi=1,x_procs
                            ! u1out(xi,yi,zi,k) = u1((xi-1)+Nxx*x_procs,(yi-1)+Nzz*z_procs,zi-1)
                            ! u2out(xi,yi,zi,k) = u2((xi-1)+Nxx*x_procs,(yi-1)+Nzz*z_procs,zi-1)
                            ! u3out(xi,yi,zi,k) = u3((xi-1)+Nxx*x_procs,(yi-1)+Nzz*z_procs,zi-1)
                            do eq=1,15
                                phiout2(eq,xi,yi,zi,k) = phi(eq,(xi-1)+Nxx*x_procs,(yi-1)+Nzz*z_procs,zi-1)
                            enddo
                        enddo
                    enddo
                enddo
                k = k + 1
            enddo
        enddo

        !各プロセスで出力
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

            open(30, file=datadir2//"0_"//trim(file_num)//"_"//trim(file_num2)//".bin", form="unformatted")
            do zi = 1, zmax+1
                do yi = 1, z_procs
                    do xi = 1, x_procs
                        ! write(30) u1out(xi,yi,zi,i), u2out(xi,yi,zi,i), u3out(xi,yi,zi,i)
                        do eq = 1, 15
                            write(30) phiout2(eq,xi,yi,zi,i)
                        enddo
                    enddo
                enddo
            enddo
            close(30)

            open(10,file="./time.d", action="write",position="append")
            write(10,*) "OUTPUT step", step, i
            close(10)
        enddo
    ! ENDDO
end program main
