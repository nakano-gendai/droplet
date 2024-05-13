program main
    implicit none
    include "mpif.h"
    integer,parameter:: xmax = 999 !ｘ方向格子数
    integer,parameter:: ymax = 220 !ｙ方向格子数
    integer,parameter:: zmax = 220 !ｚ方向格子数
    integer i, j, k, n, xi, yi, zi
    !MPI用変数
    integer ierr,comm_procs,comm_rank
    integer Nx,Ny !xi、yi方向の分割数
    integer x_procs,y_procs
    integer key_x, group_x, key_y, group_y
    integer newy_comm_world,newy_procs,newy_rank
    integer newx_comm_world,newx_procs,newx_rank
    integer rectangle_type,rectangle_type_f,xtype,xtype_f,output_type
    integer next_rank_x,former_rank_x,next_rank_y,former_rank_y
    integer req1s,req1r,req2s,req2r,req3s,req3r,req4s,req4r,req5s,req5r,req6s,req6r,req7s,req7r,req8s,req8r,req9s,req9r,req10s,req10r
    integer, dimension(MPI_STATUS_SIZE) :: sta1s,sta1r,sta2s,sta2r,sta3s,sta3r,sta4s,sta4r,sta5s,sta5r,sta6s,sta6r,sta7s,sta7r,sta8s,sta8r,sta9s,sta9r,sta10s,sta10r
    integer Nyy,Nxx
    real(8), allocatable:: tmp1(:,:,:,:),tmp2(:,:,:,:),tmp3(:,:,:,:),tmp4(:,:,:,:),tmp5(:,:,:,:),tmpf(:,:,:,:,:),tmpff(:,:,:,:,:),tmp(:,:,:),tmp_f(:,:,:,:),tmp_nu(:,:,:,:)
    integer tags, tagr, recv_rank
    character(2) chmyrank

    real(8) u1(-1:xmax+1,-1:ymax+1,-1:zmax+1)
    real(8),allocatable :: u1_procs(:,:,:)
    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, comm_procs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, comm_rank, ierr)

    !配列のallocate(xi:0とx_procs+1がのりしろ)(yi:0とymax+2がのりしろ)(zi:0とy_procs+1がのりしろ)
    Nx = 100 !ｘ方向の並列数（ただし，Nx/=comm_procs）
    Ny = comm_procs / Nx !ｙ方向の並列数
    x_procs = (xmax+1) / Nx
    y_procs = (zmax+1) / Ny
    !以下はのりしろ有りの変数
    allocate(f_procs(0:x_procs+1,0:ymax+2,0:y_procs+1))
    allocate(tmp(1:x_procs,1:ymax+1,1:y_procs))
    allocate(tmp1(0:comm_procs-1,1:x_procs,1:ymax+1,1:y_procs))
    
    !初期化
    u1_procs(:,:,:) = 0.0d0
    !u1の初期条件
    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                u1(xi,yi,zi) = 100.0d0*dble(xi)+10.0d0*dble(yi)+dble(zi)
            enddo
        enddo
    enddo
    n = 0
    do Nyy=0,Ny-1
        do Nxx=0,Nx-1
            do zi=1,y_procs
                do yi=1,ymax+1
                    do xi=1,x_procs
                        tmp1(n,xi,yi,zi) = u1((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs))
                    enddo
                enddo
            enddo
            n = n + 1
        enddo
    enddo
    do n=0,comm_procs-1
        if(n == comm_rank) then
            do zi=1,y_procs
                do yi=1,ymax+1
                    do xi=1,x_procs
                        u1_procs(xi,yi,zi) = tmp1(n,xi,yi,zi)
                    enddo
                enddo
            enddo
        endif
    enddo
    write(*,*) "0"

    !u1
    tmp(:,:,:) = 0.0d0
    do zi=1,y_procs
        do yi=1,ymax+1
            do xi=1,x_procs
                tmp(xi,yi,zi) = u1_procs(xi,yi,zi)
            enddo
        enddo
    enddo
    write(*,*) "1"
    do i=1,comm_procs-1
        if(comm_rank == i) then
            tags = i
            call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,tags,MPI_COMM_WORLD,req1s,ierr)
            call MPI_Wait(req1s,sta1s,ierr)
            if(comm_rank == 0) then
                tagr = i
                recv_rank = i
                call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,recv_rank,tagr,MPI_COMM_WORLD,req1r,ierr)
                call MPI_Wait(req1r,sta1r,ierr)
                write(*,*) "recv",comm_rank
                do zi=1,y_procs
                    do yi=1,ymax+1
                        do xi=1,x_procs
                            tmp1(recv_rank,xi,yi,zi) = tmp(xi,yi,zi)
                        enddo
                    enddo
                enddo
            endif
        endif
    enddo
    if(comm_rank == 0) then
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp1(0,xi,yi,zi) = u1_procs(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    !物理量をまとめる
    if(comm_rank == 0) then
        k = 0
        do Nyy=0,Ny-1
            do Nxx=0,Nx-1
                do zi=1,y_procs
                    do yi=1,ymax+1
                        do xi=1,x_procs
                            u1((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*y_procs) = tmp1(k,xi,yi,zi)
                        enddo
                    enddo
                enddo
                k = k + 1
            enddo
        enddo
    endif
    write(*,*) "finish",comm_rank
    if(comm_rank==0) then
        open(10,file="u1.d")
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    write(10,*) xi,yi,zi,u1(xi,yi,zi)-(dble(xi)*100.0d0+dble(yi)*10.0d0+dble(zi))
                enddo
            enddo
        enddo
        close(10)
    endif

    call MPI_Finalize(ierr)
end program main