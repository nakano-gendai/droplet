
program main
use globals
    implicit none
    integer,parameter:: xmax = 255 !ｘ方向格子数
    integer,parameter:: ymax = 255 !ｙ方向格子数
    integer,parameter:: zmax = 255 !ｚ方向格子数
    integer i, xi, yi, zi, n
    character(8) file_num, file_num2
    character :: filename*200
    !ティレクトリー読み込み
    character(*),parameter :: datadir = "/data/sht/nakanog/omp/fg/collect/"
    !ティレクトリー読み込み
    character(*),parameter :: datadir2 = "/data/sht/nakanog/mpi/fg/collect/"

    real(8) f1(1:15,0:xmax,0:ymax,0:zmax)
    real(8) f2(1:15,0:xmax,0:ymax,0:zmax)
    real(8) gosa(1:15,0:xmax,0:ymax,0:zmax)
    real(8) wa


    f1(:,:,:,:) = 0.0d0
    f2(:,:,:,:) = 0.0d0

    filename=datadir//"omp_1000.bin"
    open(10,file=filename, form="unformatted")
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                do i = 1, 15
                    read(10)  f1(i,xi,yi,zi)
                enddo
            enddo
        enddo
    enddo
    close(10)

    filename=datadir2//"mpi_1000.bin"
    open(11,file=filename, form="unformatted")
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                do i = 1, 15
                    read(11)  f2(i,xi,yi,zi)
                enddo
            enddo
        enddo
    enddo
    close(11)

    wa = 0.0d0
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                do i = 1, 15
                    gosa(i,xi,yi,zi) = (abs(f1(i,xi,yi,zi) - f2(i,xi,yi,zi))) / f2(i,xi,yi,zi)

                    wa = wa + gosa(i,xi,yi,zi)
                enddo
            enddo
        enddo
    enddo

    open(12,file="./gosa3.d")
    n = 0
    do zi = 0, zmax
        do yi = 0, ymax
            do xi = 0, xmax
                do i = 1, 15
                    write(12,*) n, gosa(i,xi,yi,zi)
                    n = n + 1
                enddo
            enddo
        enddo
    enddo
    write(12,*) wa
    close(12)

end program main
