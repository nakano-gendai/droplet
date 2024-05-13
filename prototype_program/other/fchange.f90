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
    integer i, k, xi, yi, zi, Nxx, Nzz, step, kk, j
    real x, y, z, s
    character(8) file_num, file_num2
    character :: filename*200
    character :: filename3*200

    integer ierr, comm_procs, comm_rank
    integer req1s, req1r, sta1s, sta1r
    !ティレクトリー読み込み
    character(*),parameter :: datadir = "/data/n/n517/re500we7_box_d-3/fg/"
    !ディレクトリ作成
    character(*),parameter :: datadir2 = "/data/n/n517/re500we7_box_d-3/fg/phi/"

    real(8) f(1:15,1:x_procs,1:ymax+1,1:z_procs)
    real(8) phi(1:x_procs,1:ymax+1,1:z_procs)
    real(8) tmp(1:x_procs,1:ymax+1,1:z_procs)
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
            
    ! filename=datadir2//"2_"//trim(adjustl(filename3))//'_100000_phi.bin' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
    ! print *, filename !表示してみる
    ! open(104, file=filename, form="unformatted", status='replace') 
    ! do zi=1,z_procs
    !     do yi=1,ymax+1
    !         do xi=1,x_procs
    !             write(104) phi(xi,yi,zi)
    !         enddo
    !     enddo
    ! enddo
    ! close(104)

    if(comm_rank == 0) then
        do zi = 1, z_procs
            do yi = 1, ymax+1
                do xi = 1, x_procs
                    tmp(xi,yi,zi) = phi(xi,yi,zi)
                enddo
            enddo
        enddo
    endif
    open(10,file=datadir2//"phi.bin",form='unformatted')
    close(10)

DO j = 0, comm_procs-1, Nx
    yi = 0
    zi = 0
    do k = 1, z_procs
        zi = k
        do kk = 1, ymax+1
            yi = kk
            do i = 0, Nx-1
                if(comm_rank == j) then
                    open(10,file=datadir2//"phi.bin",form='unformatted',action="write",position="append")
                    write(*,*) j, xi, yi, zi !表示してみる
                    do xi = 1, x_procs
                        write(10) tmp(xi,yi,zi)
                    enddo
                    close(10)
                endif
                if(mod(i+1,Nx) /= 0) then
                    if(comm_rank == i+1+j) then
                        call MPI_Isend(phi(1,1,1),x_procs*(ymax+1)*z_procs,MPI_REAL8,j,1,MPI_COMM_WORLD,req1s,ierr)
                        call MPI_Wait(req1s,sta1s,ierr)
                    endif
                    if(comm_rank == j) then
                        call MPI_Irecv(tmp(1,1,1),x_procs*(ymax+1)*z_procs,MPI_REAL8,i+1+j,1,MPI_COMM_WORLD,req1r,ierr)
                        call MPI_Wait(req1r,sta1r,ierr)
                    endif
                elseif((mod(i+1,Nx) == 0) .and. (j /= comm_procs-1)) then
                    if(comm_rank == j) then
                        do zi = 1, z_procs
                            do yi = 1, ymax+1
                                do xi = 1, x_procs
                                    tmp(xi,yi,zi) = phi(xi,yi,zi)
                                enddo
                            enddo
                        enddo
                    endif
                endif
            enddo
            ! if(kk == ymax+1) then
            !     yi = 1
            ! endif
        enddo
        ! if(k == z_procs) then
        !     zi = 1
        ! endif
    enddo
ENDDO

    call MPI_Finalize(ierr)
end program main
