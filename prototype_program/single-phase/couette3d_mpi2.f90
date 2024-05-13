module globals
    include "mpif.h"
    real(8),parameter:: ds = 1.0d0 !格子間隔（lattice unit）
    integer,parameter:: xmax = 3 !ｘ方向格子数
    integer,parameter:: ymax = 3 !ｙ方向格子数
    integer,parameter:: zmax = 99 !ｚ方向格子数
    integer,parameter:: step = 100000000
    integer i, n, xi, yi, zi
    real(8) err,ut,errs

    !粒子速度（整数）
    integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
    integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
    integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)
    real(8):: cr(1:3, 1:15)  !粒子速度（実数）
    
    !パラメータ
    real(8),parameter:: rho = 1.0d0 !液体の密度
    real(8),parameter:: nu = 0.01d0 !液体の動粘度
    real(8),parameter:: U(3) = (/0.01d0, 0.0d0, 0.0d0/) !壁の速度
    real(8),parameter:: dp = 0.0d0
    real(8) c
    real(8),parameter:: tau = 0.5d0 + 3.0d0 * nu / ds
    real(8),parameter:: ep = 1.0d0 / tau
    real(8),parameter:: E(15) = (/ 2.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, &
                                1.0d0/9.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, &
                                1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0 /)
    !MPI用変数
    integer ierr,comm_procs,comm_rank
    integer Nx,Ny !xi、yi方向の分割数
    integer x_procs,y_procs
    integer key_x, group_x, key_y, group_y
    integer newy_comm_world,newy_procs,newy_rank
    integer newx_comm_world,newx_procs,newx_rank
    integer rectangle_type,rectangle_type_f,xtype,xtype_f
    character(2) chmyrank
    integer next_rank_x,former_rank_x,next_rank_y,former_rank_y
    integer req1s,req1r,req2s,req2r,req3s,req3r,req4s,req4r
    integer, dimension(MPI_STATUS_SIZE) :: sta1s,sta1r,sta2s,sta2r,sta3s,sta3r,sta4s,sta4r
    contains
    subroutine parv(cx,cy,cz,cr)
        integer, intent(in):: cx(15),cy(15),cz(15)
        real(8), intent(out):: cr(1:3,1:15)
        do i=1,15
            cr(1,i) = dble(cx(i))
            cr(2,i) = dble(cy(i))
            cr(3,i) = dble(cz(i))
        enddo
    end subroutine parv

    subroutine initial(p_procs,u1_procs,u2_procs,u3_procs,f_procs,fnext_procs)
        real(8), allocatable :: f_procs(:,:,:,:),fnext_procs(:,:,:,:)
        real(8), allocatable :: u1_procs(:,:,:),u2_procs(:,:,:),u3_procs(:,:,:),p_procs(:,:,:)

        !配列のallocate(xi:0とx_procs+1がのりしろ)(yi:0とy_procs+1がのりしろ)
        Nx = 2 !ｘ方向の並列数（ただし，Nx/=comm_procs）
        Ny = comm_procs / Nx !ｙ方向の並列数
        x_procs = (xmax+1) / Nx
        y_procs = (zmax+1) / Ny
        allocate(p_procs(0:x_procs+1,0:ymax+2,0:y_procs+1))
        allocate(u1_procs(0:x_procs+1,0:ymax+2,0:y_procs+1))
        allocate(u2_procs(0:x_procs+1,0:ymax+2,0:y_procs+1))
        allocate(u3_procs(0:x_procs+1,0:ymax+2,0:y_procs+1))
        allocate(f_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1))
        allocate(fnext_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1))
        !初期化
        p_procs(:,:,:) = 0.0d0
        u1_procs(:,:,:) = 0.0d0
        u2_procs(:,:,:) = 0.0d0
        u3_procs(:,:,:) = 0.0d0
        f_procs(:,:,:,:) = 0.0d0
        fnext_procs(:,:,:,:) = 0.0d0
        do zi=0,y_procs+1
            do yi=0,ymax+2
                do xi=0,x_procs+1
                    p_procs(xi,yi,zi) = rho / 3.0d0
                    u1_procs(xi,yi,zi) = 0.0d0
                    u2_procs(xi,yi,zi) = 0.0d0
                    u3_procs(xi,yi,zi) = 0.0d0
                    do i =1,15
                        f_procs(i,xi,yi,zi) = E(i)*(3.0d0*p_procs(xi,yi,zi)+3.0d0*(cr(1,i)*u1_procs(xi,yi,zi)+cr(2,i)*u2_procs(xi,yi,zi) &
                        +cr(3,i)*u3_procs(xi,yi,zi))+4.5d0*(cr(1,i)*u1_procs(xi,yi,zi)+cr(2,i)*u2_procs(xi,yi,zi)+cr(3,i)*u3_procs(xi,yi,zi))**2 &
                        -1.5d0*(u1_procs(xi,yi,zi)**2+u2_procs(xi,yi,zi)**2+u3_procs(xi,yi,zi)**2))
                    enddo
                enddo
            enddo
        enddo
        !z方向にコミュニケータを分割する
        key_y = comm_rank
        group_y = comm_rank / Nx
        call MPI_Comm_Split(MPI_COMM_WORLD,group_y,key_y,newy_comm_world,ierr)
        call MPI_Comm_Size(newy_comm_world,newy_procs,ierr)
        call MPI_Comm_Rank(newy_comm_world,newy_rank,ierr)

        !x方向にコミュニケータを分割する
        key_x = comm_rank
        group_x = mod(comm_rank,Nx)
        call MPI_Comm_Split(MPI_COMM_WORLD,group_x,key_x,newx_comm_world,ierr)
        call MPI_Comm_Size(newx_comm_world,newx_procs,ierr)
        call MPI_Comm_Rank(newx_comm_world,newx_rank,ierr)

        !のりしろ境界の通信先設定
        next_rank_x = newx_rank + 1
        former_rank_x = newx_rank - 1
        if(newx_rank == 0) then
            former_rank_x = Ny - 1
        else if(newx_rank == Ny - 1) then
            next_rank_x = 0
        endif
        next_rank_y = newy_rank + 1
        former_rank_y = newy_rank - 1
        if(newy_rank == 0) then
            former_rank_y = Nx - 1
        else if(newy_rank == Nx - 1) then
            next_rank_y = 0
        endif
        !X世界での受け渡しをする際の型作成
        call MPI_Type_Vector(ymax+3,x_procs,x_procs+2,MPI_REAL8,xtype,ierr)
        call MPI_Type_Commit(xtype,ierr)
        call MPI_Type_Vector(ymax+3,15*x_procs,15*(x_procs+2),MPI_REAL8,xtype_f,ierr)
        call MPI_Type_Commit(xtype_f,ierr)

        !Y世界での受け渡しをする際の型作成
        call MPI_Type_Vector((y_procs+2)*(ymax+3),1,x_procs+2,MPI_REAL8,rectangle_type,ierr)
        call MPI_Type_Commit(rectangle_type,ierr)
        call MPI_Type_Vector((y_procs+2)*(ymax+3),15,15*(x_procs+2),MPI_REAL8,rectangle_type_f,ierr)
        call MPI_Type_Commit(rectangle_type_f,ierr)
    end subroutine initial

    subroutine glue(var)
        real(8),intent(inout) :: var(0:x_procs+1,0:ymax+2,0:y_procs+1)
        do zi=1,y_procs
            do xi=1,x_procs
                var(xi,0,zi) = var(xi,ymax+1,zi)
                var(xi,ymax+2,zi) = var(xi,1,zi)
            enddo
        enddo
        !X世界でののりしろ通信
        call MPI_Isend(var(1,0,y_procs),1,xtype,next_rank_x,1,newx_comm_world,req1s,ierr)
        call MPI_Irecv(var(1,0,0),1,xtype,former_rank_x,1,newx_comm_world,req1r,ierr)

        call MPI_Isend(var(1,0,1),1,xtype,former_rank_x,2,newx_comm_world,req2s,ierr)
        call MPI_Irecv(var(1,0,y_procs+1),1,xtype,next_rank_x,2,newx_comm_world,req2r,ierr)

        call MPI_Wait(req1s,sta1s,ierr)
        call MPI_Wait(req1r,sta1r,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
        call MPI_Wait(req2r,sta2r,ierr)
        !Y世界でののりしろ通信
        call MPI_Isend(var(x_procs,0,0),1,rectangle_type,next_rank_y,1,newy_comm_world,req1s,ierr)
        call MPI_Irecv(var(0,0,0),1,rectangle_type,former_rank_y,1,newy_comm_world,req1r,ierr)

        call MPI_Isend(var(1,0,0),1,rectangle_type,former_rank_y,2,newy_comm_world,req2s,ierr)
        call MPI_Irecv(var(x_procs+1,0,0),1,rectangle_type,next_rank_y,2,newy_comm_world,req2r,ierr)

        call MPI_Wait(req1s,sta1s,ierr)
        call MPI_Wait(req1r,sta1r,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
        call MPI_Wait(req2r,sta2r,ierr)
    end subroutine glue

    subroutine glue_f(var_f)
        real(8),intent(inout) :: var_f(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        do zi=1,y_procs
            do xi=1,x_procs
                do i=1,15
                    var_f(i,xi,0,zi) = var_f(i,xi,ymax+1,zi)
                    var_f(i,xi,ymax+2,zi) = var_f(i,xi,1,zi)
                enddo
            enddo
        enddo
        !X世界でののりしろ通信
        call MPI_Isend(var_f(1,1,0,y_procs),1,xtype_f,next_rank_x,1,newx_comm_world,req1s,ierr)
        call MPI_Irecv(var_f(1,1,0,0),1,xtype_f,former_rank_x,1,newx_comm_world,req1r,ierr)

        call MPI_Isend(var_f(1,1,0,1),1,xtype_f,former_rank_x,2,newx_comm_world,req2s,ierr)
        call MPI_Irecv(var_f(1,1,0,y_procs+1),1,xtype_f,next_rank_x,2,newx_comm_world,req2r,ierr)

        call MPI_Wait(req1s,sta1s,ierr)
        call MPI_Wait(req1r,sta1r,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
        call MPI_Wait(req2r,sta2r,ierr)

        !Y世界でののりしろ通信
        call MPI_Isend(var_f(1,x_procs,0,0),1,rectangle_type_f,next_rank_y,1,newy_comm_world,req1s,ierr)
        call MPI_Irecv(var_f(1,0,0,0),1,rectangle_type_f,former_rank_y,1,newy_comm_world,req1r,ierr)

        call MPI_Isend(var_f(1,1,0,0),1,rectangle_type_f,former_rank_y,2,newy_comm_world,req2s,ierr)
        call MPI_Irecv(var_f(1,x_procs+1,0,0),1,rectangle_type_f,next_rank_y,2,newy_comm_world,req2r,ierr)

        call MPI_Wait(req1s,sta1s,ierr)
        call MPI_Wait(req1r,sta1r,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
        call MPI_Wait(req2r,sta2r,ierr)
    end subroutine glue_f

    subroutine MPI_boundary(p_procs,u1_procs,u2_procs,u3_procs,f_procs)
        real(8),intent(inout) :: f_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(inout) :: u1_procs(0:x_procs+1,0:ymax+2,0:y_procs+1),u2_procs(0:x_procs+1,0:ymax+2,0:y_procs+1), &
                                u3_procs(0:x_procs+1,0:ymax+2,0:y_procs+1),p_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        call glue(u1_procs)
        call glue(u2_procs)
        call glue(u3_procs)
        call glue(p_procs)
        call glue_f(f_procs)
        ! write(chmyrank, '(i2.2)') comm_rank   !myrankを文字型変数に格納
        ! open(60, file = './tasikame_'//chmyrank//'.d')    !ランク番号が付いたファイルを作成
        ! write(60, *) 'Myrank :', comm_rank
        ! do zi=0,y_procs+1
        !     do yi=0,ymax+2
        !         do xi=0,x_procs+1
        !                 write(60, *) xi,zi,u1_procs(xi,yi,zi)
        !         enddo
        !     enddo
        ! enddo
        ! close(60)
    endsubroutine MPI_boundary

    subroutine bulk_area(p_procs,u1_procs,u2_procs,u3_procs,f_procs,fnext_procs)
        real(8),intent(in) :: f_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(in) :: u1_procs(0:x_procs+1,0:ymax+2,0:y_procs+1),u2_procs(0:x_procs+1,0:ymax+2,0:y_procs+1), &
                                u3_procs(0:x_procs+1,0:ymax+2,0:y_procs+1),p_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(out) :: fnext_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        fnext_procs(i,xi,yi,zi) = f_procs(i,xi-cx(i),yi-cy(i),zi-cz(i)) - ep * (f_procs(i,xi-cx(i),yi-cy(i),zi-cz(i)) &
                            - feq(i,p_procs(xi-cx(i),yi-cy(i),zi-cz(i)),u1_procs(xi-cx(i),yi-cy(i),zi-cz(i)),u2_procs(xi-cx(i),yi-cy(i),zi-cz(i)),u3_procs(xi-cx(i),yi-cy(i),zi-cz(i))))
                    enddo
                enddo
            enddo
        enddo
    end subroutine bulk_area

    subroutine periodic_pressure(fnext_procs)
        real(8),intent(inout) :: fnext_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        call glue_f(fnext_procs) !のりしろ境界
        do zi=1,y_procs
            do yi=1,ymax+1
                c = dp - (fnext_procs(1,1,yi,zi)-fnext_procs(1,0,yi,zi) + fnext_procs(3,1,yi,zi)-fnext_procs(3,0,yi,zi) &
                        + fnext_procs(4,1,yi,zi)-fnext_procs(4,0,yi,zi) + fnext_procs(6,1,yi,zi)-fnext_procs(6,0,yi,zi) &
                        + fnext_procs(7,1,yi,zi)-fnext_procs(7,0,yi,zi)) / 3.0d0
                if(newy_rank == 0) then
                    !左境界
                    fnext_procs(2,1,yi,zi) = fnext_procs(2,0,yi,zi) + c
                    fnext_procs(8,1,yi,zi) = fnext_procs(8,0,yi,zi) + 0.125d0 * c
                    fnext_procs(10,1,yi,zi) = fnext_procs(10,0,yi,zi) + 0.125d0 * c
                    fnext_procs(11,1,yi,zi) = fnext_procs(11,0,yi,zi) + 0.125d0 * c
                    fnext_procs(13,1,yi,zi) = fnext_procs(13,0,yi,zi) + 0.125d0 * c
                else if(newy_rank == newy_procs-1) then
                    !右境界
                    fnext_procs(5,x_procs,yi,zi) = fnext_procs(5,x_procs+1,yi,zi) - c
                    fnext_procs(9,x_procs,yi,zi) = fnext_procs(9,x_procs+1,yi,zi) - 0.125d0 * c
                    fnext_procs(12,x_procs,yi,zi) = fnext_procs(12,x_procs+1,yi,zi) - 0.125d0 * c
                    fnext_procs(14,x_procs,yi,zi) = fnext_procs(14,x_procs+1,yi,zi) - 0.125d0 * c
                    fnext_procs(15,x_procs,yi,zi) = fnext_procs(15,x_procs+1,yi,zi) - 0.125d0 * c
                endif
            enddo
        enddo
    end subroutine periodic_pressure

    subroutine bounce_back(cr,fnext_procs)
        real(8),intent(in) :: cr(1:3,1:15)
        real(8),intent(inout) :: fnext_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        !下の壁
        if(group_y == 0) then
            zi = 1
            do yi=1,ymax+1
                do xi=1,x_procs
                    fnext_procs(4,xi,yi,zi) = fnext_procs(7,xi,yi,zi) - 6.0d0 * E(7) * ((-U(1))*cr(1,7)+(-U(2))*cr(2,7)+(-U(3))*cr(3,7))
                    fnext_procs(8,xi,yi,zi) = fnext_procs(12,xi,yi,zi) - 6.0d0 * E(12) * ((-U(1))*cr(1,12)+(-U(2))*cr(2,12)+(-U(3))*cr(3,12))
                    fnext_procs(9,xi,yi,zi) = fnext_procs(13,xi,yi,zi) - 6.0d0 * E(13) * ((-U(1))*cr(1,13)+(-U(2))*cr(2,13)+(-U(3))*cr(3,13))
                    fnext_procs(10,xi,yi,zi) = fnext_procs(14,xi,yi,zi) - 6.0d0 * E(14) * ((-U(1))*cr(1,14)+(-U(2))*cr(2,14)+(-U(3))*cr(3,14))
                    fnext_procs(15,xi,yi,zi) = fnext_procs(11,xi,yi,zi) - 6.0d0 * E(11) * ((-U(1))*cr(1,11)+(-U(2))*cr(2,11)+(-U(3))*cr(3,11))
                enddo
            enddo
        !上の壁
        else if(group_y == Ny-1) then
            zi = y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    fnext_procs(7,xi,yi,zi) = fnext_procs(4,xi,yi,zi) - 6.0d0 * E(4) * ((U(1))*cr(1,4)+(U(2))*cr(2,4)+(U(3))*cr(3,4))
                    fnext_procs(11,xi,yi,zi) = fnext_procs(15,xi,yi,zi) - 6.0d0 * E(15) * ((U(1))*cr(1,15)+(U(2))*cr(2,15)+(U(3))*cr(3,15))
                    fnext_procs(12,xi,yi,zi) = fnext_procs(8,xi,yi,zi) - 6.0d0 * E(8) * ((U(1))*cr(1,8)+(U(2))*cr(2,8)+(U(3))*cr(3,8))
                    fnext_procs(13,xi,yi,zi) = fnext_procs(9,xi,yi,zi) - 6.0d0 * E(9) * ((U(1))*cr(1,9)+(U(2))*cr(2,9)+(U(3))*cr(3,9))
                    fnext_procs(14,xi,yi,zi) = fnext_procs(10,xi,yi,zi) - 6.0d0 * E(10) * ((U(1))*cr(1,10)+(U(2))*cr(2,10)+(U(3))*cr(3,10))
                enddo
            enddo
        endif
    end subroutine bounce_back

    subroutine renew(f_procs,fnext_procs)
        real(8),intent(in) :: fnext_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(out) :: f_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        f_procs(i,xi,yi,zi) = fnext_procs(i,xi,yi,zi)
                    enddo
                enddo
            enddo
        enddo
    end subroutine renew

    subroutine reset(p_procs,u1_procs,u2_procs,u3_procs)
        real(8),intent(inout) :: u1_procs(0:x_procs+1,0:ymax+2,0:y_procs+1),u2_procs(0:x_procs+1,0:ymax+2,0:y_procs+1), &
                                u3_procs(0:x_procs+1,0:ymax+2,0:y_procs+1),p_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    p_procs(xi,yi,zi) = 0.0d0
                    u1_procs(xi,yi,zi) = 0.0d0
                    u2_procs(xi,yi,zi) = 0.0d0
                    u3_procs(xi,yi,zi) = 0.0d0
                enddo
            enddo
        enddo
    end subroutine reset

    subroutine pressure_cal(f_procs,p_procs)
        real(8),intent(in) :: f_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(out) :: p_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i =1,15
                        p_procs(xi,yi,zi) = p_procs(xi,yi,zi) + f_procs(i,xi,yi,zi)
                    enddo
                    p_procs(xi,yi,zi) = p_procs(xi,yi,zi) / 3.0d0
                enddo
            enddo
        enddo
    end subroutine pressure_cal

    subroutine velocity_cal(f_procs,u1_procs,u2_procs,u3_procs)
        real(8),intent(in) :: f_procs(1:15,0:x_procs+1,0:ymax+2,0:y_procs+1)
        real(8),intent(out) :: u1_procs(0:x_procs+1,0:ymax+2,0:y_procs+1),u2_procs(0:x_procs+1,0:ymax+2,0:y_procs+1), &
                                u3_procs(0:x_procs+1,0:ymax+2,0:y_procs+1)
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i =1,15
                        u1_procs(xi,yi,zi) = u1_procs(xi,yi,zi) + f_procs(i,xi,yi,zi) * cr(1,i)
                        u2_procs(xi,yi,zi) = u2_procs(xi,yi,zi) + f_procs(i,xi,yi,zi) * cr(2,i)
                        u3_procs(xi,yi,zi) = u3_procs(xi,yi,zi) + f_procs(i,xi,yi,zi) * cr(3,i)
                    enddo
                enddo
            enddo
        enddo
    end subroutine velocity_cal

    function feq(i,p, u1, u2, u3) result(feqr)
        integer,intent(in) :: i
        real(8),intent(in) :: p,u1,u2,u3
        real(8) feqr
        feqr = E(i)*(3.0d0*p + 3.0d0*(u1*dble(cx(i)) + u2*dble(cy(i)) + u3*dble(cz(i))) &
                    + 4.5d0*(u1*dble(cx(i)) + u2*dble(cy(i)) + u3*dble(cz(i)))**2 -1.5d0*(u1**2 + u2**2 + u3**2))
    end function feq
    
end module globals
    
program main
use globals
    implicit none
    ! real(8) f(1:15, 0:xmax, 0:ymax, 0:zmax)
    ! real(8) p(0:xmax, 0:ymax, 0:zmax)  !圧力
    ! real(8) u1(0:xmax, 0:ymax, 0:zmax), u2(0:xmax, 0:ymax, 0:zmax), u3(0:xmax, 0:ymax, 0:zmax) !流速
    real(8),allocatable :: f_procs(:,:,:,:), fnext_procs(:,:,:,:) 
    real(8),allocatable :: u1_procs(:,:,:),u2_procs(:,:,:),u3_procs(:,:,:),p_procs(:,:,:)
    real(8) time1,time2
    call cpu_time(time1)

    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, comm_procs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, comm_rank, ierr)

    call parv(cx,cy,cz,cr)
    !初期値の設定
    call initial(p_procs,u1_procs,u2_procs,u3_procs,f_procs,fnext_procs)
    write(*,*) "initial is OK!"
!===========時間発展=====================
    DO n=1,step
        !のりしろ境界の通信
        call MPI_boundary(p_procs,u1_procs,u2_procs,u3_procs,f_procs)
        write(*,*) "MPI_boundary is OK!"
        !次の時刻の速度分布関数fnextを計算
        call bulk_area(p_procs,u1_procs,u2_procs,u3_procs,f_procs,fnext_procs)
        write(*,*) "bulk_area is OK!"
        !左右 周期境界条件（dp=0）
        call periodic_pressure(fnext_procs)
        write(*,*) "periodic_pressure is OK!"
        !壁
        call bounce_back(cr,fnext_procs)
        write(*,*) "bounceback is OK!"
        ! 速度分布関数の更新
        call renew(f_procs,fnext_procs)
        write(*,*) "renew is OK!"
        !初期化
        call reset(p_procs,u1_procs,u2_procs,u3_procs)
        write(*,*) "reset is OK!"
        !圧力の計算
        call pressure_cal(f_procs,p_procs)
        write(*,*) "pressure is OK!"
        !流速の計算
        call velocity_cal(f_procs,u1_procs,u2_procs,u3_procs)
        write(*,*) "velocity is OK!"
        write(*,*) "step = ", n
    ENDDO
    write(chmyrank, '(i2.2)') comm_rank   !myrankを文字型変数に格納
    open(60, file = './mpi2_3d_'//chmyrank//'.d')    !ランク番号が付いたファイルを作成
    write(60, *) 'Myrank :', comm_rank
    do zi=1,y_procs
        do yi=1,ymax+1
            do xi=1,x_procs
                write(60, *) xi,yi,zi,u1_procs(xi,yi,zi)
            enddo
        enddo
    enddo
    close(60)
    write(*,*) "output is OK!"

    call MPI_Finalize(ierr)

    call cpu_time(time2)
    write(*,*) time2-time1
end program main
