module globals
    include "mpif.h"
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
    integer,parameter:: xmax = 999 !ｘ方向格子数
    integer,parameter:: ymax = 290 !ｙ方向格子数
    integer,parameter:: zmax = 290 !ｚ方向格子数
    integer,parameter:: Nx = 25 !ｘ方向の並列数
    integer,parameter:: Nz = 97 !ｚ方向の並列数
    integer,parameter:: x_procs = (xmax+1) / Nx
    integer,parameter:: z_procs = (zmax+1) / Nz
    integer,parameter:: new_procs = Nx * Nz
    integer i, k, xi, yi, zi, Nxx, Nzz, step, kk, j, n, start
    real x, y, z, s
    character(8) file_num, file_num2
    character :: filename*200
    character :: filename3*200
    integer tags, tagr, recv_rank

    integer ierr, comm_procs, comm_rank
    integer req1s, req1r, sta1s, sta1r
    !ティレクトリー読み込み
    character(*),parameter :: datadir = "/data/n/n517/re500we7_box_d-3/fg/"
    !ディレクトリ作成
    character(*),parameter :: datadir2 = "/data/n/n517/re500we7_box_d-3/fg/phi/"

    real(8) f(1:15,1:x_procs,1:ymax+1,1:z_procs)
    real(8) phi(1:x_procs,1:ymax+1,1:z_procs)
    real(8) phiout(1:x_procs*Nx,1:ymax+1,1:z_procs)
    real(8) tmp(1:x_procs,1:ymax+1,1:z_procs)
    real(8) tmp1(0:Nx-1,1:x_procs,1:ymax+1,1:z_procs)
    real(8) dummy

    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, comm_procs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, comm_rank, ierr)

    call mk_dirs(datadir2)

    write(filename3,*) comm_rank
    filename=datadir//'0_'//trim(adjustl(filename3))//'_90000_fg.bin' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
    print *, filename !表示してみる
    open(103, file=filename, form="unformatted")
    do zi=1,z_procs
        do yi=1,ymax+1
            do xi=1,x_procs
                do i=1,15
                    read(103) f(i,xi,yi,zi), dummy
                enddo
            enddo
        enddo
    enddo
    close(103)

    do zi=1,z_procs
        do yi=1,ymax+1
            do xi=1,x_procs
                phi(xi,yi,zi) = 0.0d0
                do i=1,15
                    phi(xi,yi,zi) = phi(xi,yi,zi) + f(i,xi,yi,zi)
                enddo
            enddo
        enddo
    enddo
            
    open(10,file=datadir2//"phi.bin",form='unformatted')
    close(10)
    start = 0
DO j=0,Nz-1
    do i=0,Nx-1
        kk = i
        if(start == 0) then
            kk = kk + 1
            start = 1
            do zi=1,z_procs
                do yi=1,ymax+1
                    do xi=1,x_procs
                        tmp1(0,xi,yi,zi) = phi(xi,yi,zi)
                    enddo
                enddo
            enddo
        endif
        if(comm_rank ==  kk + (Nx * j)) then
            tags = kk + (Nx * j)
            call MPI_Isend(phi(1,1,1),x_procs*z_procs*(ymax+1),MPI_REAL8,0,tags,MPI_COMM_WORLD,req1s,ierr)
            call MPI_Wait(req1s,sta1s,ierr)
        endif
        if(comm_rank == 0) then
            tagr = kk + (Nx * j)
            recv_rank = kk + (Nx * j)
            call MPI_Irecv(tmp(1,1,1),x_procs*z_procs*(ymax+1),MPI_REAL8,recv_rank,tagr,MPI_COMM_WORLD,req1r,ierr)
            call MPI_Wait(req1r,sta1r,ierr)
            write(*,*) j
            do zi=1,z_procs
                do yi=1,ymax+1
                    do xi=1,x_procs
                        tmp1(recv_rank,xi,yi,zi) = tmp(xi,yi,zi)
                    enddo
                enddo
            enddo
        endif
    enddo

    if(comm_rank == 0) then
        k = 0
        do Nxx=0,Nx-1
            do zi=1,z_procs
                do yi=1,ymax+1
                    do xi=1,x_procs
                        phiout((xi)+Nxx*x_procs,yi,(zi)) = tmp1(k,xi,yi,zi)
                    enddo
                enddo
            enddo
            k = k + 1
        enddo

        open(10,file=datadir2//"phi.bin",form='unformatted',action="write",position="append")
        do zi = 1, z_procs
            do xi = 1, x_procs*Nx
                do yi = 1, ymax+1
                    write(10) phiout(xi,yi,zi)
                enddo
            enddo
        enddo
        close(10)
        
    endif

ENDDO
    call MPI_Finalize(ierr)
end program main