module globals
    include "mpif.h"
    real(8),parameter:: ds = 1.0d0 !格子間隔（lattice unit）
    integer,parameter:: xmax = 39 !ｘ方向格子数
    integer,parameter:: ymax = 1 !ｙ方向格子数
    integer,parameter:: zmax = 39 !ｚ方向格子数
    integer,parameter:: step = 1
    integer i, n, xi, yi, zi
    real(8) err_u1,err_up_u1,err_down_u1,err_u2,err_up_u2,err_down_u2, &
            err_p,err_up_p,err_down_p
    real(8) ut1,ut2,pt
    character(2) chmyrank
    !MPI用変数
    integer ierr,comm_procs,comm_rank
    integer Nx,Ny !xi、yi方向の分割数
    integer x_procs,y_procs
    integer key_x, group_x, key_y, group_y
    integer newy_comm_world,newy_procs,newy_rank
    integer newx_comm_world,newx_procs,newx_rank
    integer rectangle_type,rectangle_type_f,xtype,xtype_f,output_type
    integer next_rank_x,former_rank_x,next_rank_y,former_rank_y
    integer req1s,req1r,req2s,req2r,req3s,req3r,req4s,req4r
    integer, dimension(MPI_STATUS_SIZE) :: sta1s,sta1r,sta2s,sta2r,sta3s,sta3r,sta4s,sta4r
    integer Nyy,Nxx
    real(8), allocatable:: tmp1(:,:,:,:),tmp2(:,:,:,:),tmp3(:,:,:,:),tmp4(:,:,:,:),tmpf(:,:,:,:,:),tmp(:,:,:),tmp_f(:,:,:,:)

    !粒子速度（整数）
    integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
    integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
    integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)
    real(8):: cr(1:3, 1:15)  !粒子速度（実数）
    
    !パラメータ
    real(8),parameter:: rho = 1.0d0 !液体の密度
    real(8),parameter:: tau = 0.8d0
    real(8),parameter:: nu = (tau - 0.5d0)/3.0d0 !液体の動粘度
    real(8),parameter:: u0 = 0.01d0 !流速の初期値
    real(8),parameter:: p0 = rho / 3.0d0 !圧力の初期値
    real(8),parameter:: ep = 1.0d0 / tau
    real(8),parameter:: pi = acos(-1.0d0)
    real(8),parameter:: k1 = 2.0d0*pi/dble(xmax+1)
    real(8),parameter:: k2 = 2.0d0*pi/dble(ymax+1)
    real(8),parameter:: k3 = 2.0d0*pi/dble(zmax+1)
    real(8),parameter:: E(15) = (/ 2.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, &
                                1.0d0/9.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, &
                                1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0 /)
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

    subroutine initial(p,u1,u2,u3,f,p_procs,u1_procs,u2_procs,u3_procs,f_procs,fnext_procs)
        real(8), intent(inout) :: p(0:xmax,0:ymax,0:zmax),u1(0:xmax,0:ymax,0:zmax),&
                                u2(0:xmax,0:ymax,0:zmax),u3(0:xmax,0:ymax,0:zmax)
        real(8),intent(inout) :: f(1:15,0:xmax,0:ymax,0:zmax)
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
        allocate(tmp1(0:comm_procs-1,1:x_procs,1:ymax+1,1:y_procs))
        allocate(tmp2(0:comm_procs-1,1:x_procs,1:ymax+1,1:y_procs))
        allocate(tmp3(0:comm_procs-1,1:x_procs,1:ymax+1,1:y_procs))
        allocate(tmp4(0:comm_procs-1,1:x_procs,1:ymax+1,1:y_procs))
        allocate(tmpf(0:comm_procs-1,1:15,1:x_procs,1:ymax+1,1:y_procs))
        !初期化
        p_procs(:,:,:) = 0.0d0
        u1_procs(:,:,:) = 0.0d0
        u2_procs(:,:,:) = 0.0d0
        u3_procs(:,:,:) = 0.0d0
        f_procs(:,:,:,:) = 0.0d0
        fnext_procs(:,:,:,:) = 0.0d0
        if(comm_rank == 0) then
            open(94,file="para.d")
            write(94,"(5es23.16)") k1,k2,k3,pi,p0
            close(94)
        endif
        if(comm_rank==0) then
        open(95,file="ini_u.d")
        open(96,file="ini.d")
        endif
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    p(xi,yi,zi) = p0 - rho * u0**2 / 4.0d0 * (cos(2.0d0*k1*dble(xi)) + &
                         k1**2 / k3**2 * cos(2.0d0*k3*dble(zi)))
                    u1(xi,yi,zi) = -u0*cos(k1*dble(xi))*sin(k3*dble(zi))
                    u2(xi,yi,zi) = 0.0d0
                    u3(xi,yi,zi) = k1/k3*u0*sin(k1*dble(xi))*cos(k3*dble(zi))
                    if(comm_rank == 0) then
                        write(95,"(7es23.16)") dble(xi),dble(yi),dble(zi),u1(xi,yi,zi),u2(xi,yi,zi),&
                                    u3(xi,yi,zi),p(xi,yi,zi)
                    endif
                    do i =1,15
                        f(i,xi,yi,zi) = E(i)*(3.0d0*p(xi,yi,zi)+3.0d0*(cr(1,i)*u1(xi,yi,zi)+cr(2,i)*u2(xi,yi,zi) &
                        +cr(3,i)*u3(xi,yi,zi))+4.5d0*(cr(1,i)*u1(xi,yi,zi)+cr(2,i)*u2(xi,yi,zi)+cr(3,i)*u3(xi,yi,zi))**2 &
                        -1.5d0*(u1(xi,yi,zi)**2+u2(xi,yi,zi)**2+u3(xi,yi,zi)**2))
                        
                        if(comm_rank == 0) then
                            write(96,"(5es23.16)") dble(xi),dble(yi),dble(zi),dble(i),f(i,xi,yi,zi)
                        endif
                    enddo
                enddo
            enddo
        enddo
        if(comm_rank==0) then
        close(95)
        close(96)
        endif
        n = 0
        do Nyy=0,Ny-1
            do Nxx=0,Nx-1
                do zi=1,y_procs
                    do yi=1,ymax+1
                        do xi=1,x_procs
                            tmp1(n,xi,yi,zi) = u1((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs))
                            tmp2(n,xi,yi,zi) = u2((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs))
                            tmp3(n,xi,yi,zi) = u3((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs))
                            tmp4(n,xi,yi,zi) = p((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs))
                            do i=1,15
                                tmpf(n,i,xi,yi,zi) = f(i,(xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs))
                            enddo
                        enddo
                    enddo
                enddo
                n = n + 1
            enddo
        enddo
        do n=0,comm_procs-1
            if(n == comm_rank)then
                do zi=1,y_procs
                    do yi=1,ymax+1
                        do xi=1,x_procs
                            u1_procs(xi,yi,zi) = tmp1(n,xi,yi,zi)
                            u2_procs(xi,yi,zi) = tmp2(n,xi,yi,zi)
                            u3_procs(xi,yi,zi) = tmp3(n,xi,yi,zi)
                            p_procs(xi,yi,zi) = tmp4(n,xi,yi,zi)
                            do i=1,15
                                f_procs(i,xi,yi,zi) = tmpf(n,i,xi,yi,zi)
                            enddo
                        enddo
                    enddo
                enddo
            endif
        enddo
        ! write(chmyrank, '(i2.2)') comm_rank   !myrankを文字型変数に格納
        ! open(60, file = './test_'//chmyrank//'.d')    !ランク番号が付いたファイルを作成
        ! write(60, *) 'Myrank :', comm_rank
        ! do zi=1,y_procs
        !     do yi=1,ymax+1
        !         do xi=1,x_procs
        !             write(60, *) xi,yi,zi,u1_procs(xi,yi,zi)
        !         enddo
        !     enddo
        ! enddo
        ! close(60)
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
        
        open(100,file="tasi.d")
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i =1,15
                        u1_procs(xi,yi,zi) = u1_procs(xi,yi,zi) + f_procs(i,xi,yi,zi) * dble(cx(i))
                        u2_procs(xi,yi,zi) = u2_procs(xi,yi,zi) + f_procs(i,xi,yi,zi) * dble(cy(i))
                        u3_procs(xi,yi,zi) = u3_procs(xi,yi,zi) + f_procs(i,xi,yi,zi) * dble(cz(i))
                    enddo
                enddo
            enddo
        enddo
        close(100)
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
    real(8) f(1:15, 0:xmax, 0:ymax, 0:zmax)
    real(8) p(0:xmax, 0:ymax, 0:zmax)  !圧力
    real(8) u1(0:xmax, 0:ymax, 0:zmax), u2(0:xmax, 0:ymax, 0:zmax), u3(0:xmax, 0:ymax, 0:zmax) !流速
    real(8),allocatable :: f_procs(:,:,:,:), fnext_procs(:,:,:,:) 
    real(8),allocatable :: u1_procs(:,:,:),u2_procs(:,:,:),u3_procs(:,:,:),p_procs(:,:,:)
    real(8) time1,time2
    call cpu_time(time1)

    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, comm_procs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, comm_rank, ierr)

    call parv(cx,cy,cz,cr)
    !初期値の設定
    call initial(p,u1,u2,u3,f,p_procs,u1_procs,u2_procs,u3_procs,f_procs,fnext_procs)
!===========時間発展=====================
    DO n=1,step
        !のりしろ境界の通信
        call MPI_boundary(p_procs,u1_procs,u2_procs,u3_procs,f_procs)
        !次の時刻の速度分布関数fnextを計算
        call bulk_area(p_procs,u1_procs,u2_procs,u3_procs,f_procs,fnext_procs)
        ! 速度分布関数の更新
        call renew(f_procs,fnext_procs)
        !初期化
        call reset(p_procs,u1_procs,u2_procs,u3_procs)
        !圧力の計算
        call pressure_cal(f_procs,p_procs)
        !流速の計算
        call velocity_cal(f_procs,u1_procs,u2_procs,u3_procs)
        ! write(*,*) "step = ", n
    ENDDO

    allocate(tmp(1:x_procs,1:ymax+1,1:y_procs))
    allocate(tmp_f(1:15,1:x_procs,1:ymax+1,1:y_procs))
    tmp1(:,:,:,:) = 0.0d0
    tmp2(:,:,:,:) = 0.0d0
    tmp3(:,:,:,:) = 0.0d0
    tmp4(:,:,:,:) = 0.0d0
    tmp(:,:,:) = 0.0d0
    !u1をまとめる
    do zi=1,y_procs
        do yi=1,ymax+1
            do xi=1,x_procs
                tmp(xi,yi,zi) = u1_procs(xi,yi,zi)
            enddo
        enddo
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

    if(comm_rank == 1) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,1,MPI_COMM_WORLD,req1s,ierr)
        call MPI_Wait(req1s,sta1s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,1,1,MPI_COMM_WORLD,req1r,ierr)
        call MPI_Wait(req1r,sta1r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp1(1,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 2) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,2,MPI_COMM_WORLD,req2s,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,2,2,MPI_COMM_WORLD,req2r,ierr)
        call MPI_Wait(req2r,sta2r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp1(2,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 3) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,3,MPI_COMM_WORLD,req3s,ierr)
        call MPI_Wait(req3s,sta3s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,3,3,MPI_COMM_WORLD,req3r,ierr)
        call MPI_Wait(req3r,sta3r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp1(3,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif
    !u2をまとめる
    tmp(:,:,:) = 0.0d0
    do zi=1,y_procs
        do yi=1,ymax+1
            do xi=1,x_procs
                tmp(xi,yi,zi) = u2_procs(xi,yi,zi)
            enddo
        enddo
    enddo
    if(comm_rank == 0) then
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp2(0,xi,yi,zi) = u2_procs(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 1) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,1,MPI_COMM_WORLD,req1s,ierr)
        call MPI_Wait(req1s,sta1s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,1,1,MPI_COMM_WORLD,req1r,ierr)
        call MPI_Wait(req1r,sta1r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp2(1,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 2) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,2,MPI_COMM_WORLD,req2s,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,2,2,MPI_COMM_WORLD,req2r,ierr)
        call MPI_Wait(req2r,sta2r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp2(2,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 3) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,3,MPI_COMM_WORLD,req3s,ierr)
        call MPI_Wait(req3s,sta3s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,3,3,MPI_COMM_WORLD,req3r,ierr)
        call MPI_Wait(req3r,sta3r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp2(3,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif
    !u3
    tmp(:,:,:) = 0.0d0
    do zi=1,y_procs
        do yi=1,ymax+1
            do xi=1,x_procs
                tmp(xi,yi,zi) = u3_procs(xi,yi,zi)
            enddo
        enddo
    enddo
    if(comm_rank == 0) then
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp3(0,xi,yi,zi) = u3_procs(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 1) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,1,MPI_COMM_WORLD,req1s,ierr)
        call MPI_Wait(req1s,sta1s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,1,1,MPI_COMM_WORLD,req1r,ierr)
        call MPI_Wait(req1r,sta1r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp3(1,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 2) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,2,MPI_COMM_WORLD,req2s,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,2,2,MPI_COMM_WORLD,req2r,ierr)
        call MPI_Wait(req2r,sta2r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp3(2,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 3) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,3,MPI_COMM_WORLD,req3s,ierr)
        call MPI_Wait(req3s,sta3s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,3,3,MPI_COMM_WORLD,req3r,ierr)
        call MPI_Wait(req3r,sta3r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp3(3,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif
    !p
    tmp(:,:,:) = 0.0d0
    do zi=1,y_procs
        do yi=1,ymax+1
            do xi=1,x_procs
                tmp(xi,yi,zi) = p_procs(xi,yi,zi)
            enddo
        enddo
    enddo
    if(comm_rank == 0) then
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp4(0,xi,yi,zi) = p_procs(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 1) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,1,MPI_COMM_WORLD,req1s,ierr)
        call MPI_Wait(req1s,sta1s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,1,1,MPI_COMM_WORLD,req1r,ierr)
        call MPI_Wait(req1r,sta1r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp4(1,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 2) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,2,MPI_COMM_WORLD,req2s,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,2,2,MPI_COMM_WORLD,req2r,ierr)
        call MPI_Wait(req2r,sta2r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp4(2,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 3) then
        call MPI_Isend(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,0,3,MPI_COMM_WORLD,req3s,ierr)
        call MPI_Wait(req3s,sta3s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp(1,1,1),x_procs*y_procs*(ymax+1),MPI_REAL8,3,3,MPI_COMM_WORLD,req3r,ierr)
        call MPI_Wait(req3r,sta3r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    tmp4(3,xi,yi,zi) = tmp(xi,yi,zi)
                enddo
            enddo
        enddo
    endif
    !fをまとめる
    tmp_f(:,:,:,:) = 0.0d0
    tmpf(:,:,:,:,:) = 0.0d0
    do zi=1,y_procs
        do yi=1,ymax+1
            do xi=1,x_procs
                do i=1,15
                    tmp_f(i,xi,yi,zi) = f_procs(i,xi,yi,zi)
                enddo
            enddo
        enddo
    enddo
    if(comm_rank == 0) then
        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        tmpf(0,i,xi,yi,zi) = f_procs(i,xi,yi,zi)
                    enddo
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 1) then
        call MPI_Isend(tmp_f(1,1,1,1),15*x_procs*y_procs*(ymax+1),MPI_REAL8,0,1,MPI_COMM_WORLD,req1s,ierr)
        call MPI_Wait(req1s,sta1s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp_f(1,1,1,1),15*x_procs*y_procs*(ymax+1),MPI_REAL8,1,1,MPI_COMM_WORLD,req1r,ierr)
        call MPI_Wait(req1r,sta1r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        tmpf(1,i,xi,yi,zi) = tmp_f(i,xi,yi,zi)
                    enddo
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 2) then
        call MPI_Isend(tmp_f(1,1,1,1),15*x_procs*y_procs*(ymax+1),MPI_REAL8,0,2,MPI_COMM_WORLD,req2s,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp_f(1,1,1,1),15*x_procs*y_procs*(ymax+1),MPI_REAL8,2,2,MPI_COMM_WORLD,req2r,ierr)
        call MPI_Wait(req2r,sta2r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        tmpf(2,i,xi,yi,zi) = tmp_f(i,xi,yi,zi)
                    enddo
                enddo
            enddo
        enddo
    endif

    if(comm_rank == 3) then
        call MPI_Isend(tmp_f(1,1,1,1),15*x_procs*y_procs*(ymax+1),MPI_REAL8,0,3,MPI_COMM_WORLD,req3s,ierr)
        call MPI_Wait(req3s,sta3s,ierr)
    endif
    if(comm_rank == 0) then
        call MPI_Irecv(tmp_f(1,1,1,1),15*x_procs*y_procs*(ymax+1),MPI_REAL8,3,3,MPI_COMM_WORLD,req3r,ierr)
        call MPI_Wait(req3r,sta3r,ierr)

        do zi=1,y_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        tmpf(3,i,xi,yi,zi) = tmp_f(i,xi,yi,zi)
                    enddo
                enddo
            enddo
        enddo
    endif

    ! if(comm_rank == 4) then
    !     call MPI_Isend(tmp_f(1,1,1,1),15*x_procs*y_procs*(ymax+1),MPI_REAL8,0,4,MPI_COMM_WORLD,req3s,ierr)
    !     call MPI_Wait(req3s,sta3s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp_f(1,1,1,1),15*x_procs*y_procs*(ymax+1),MPI_REAL8,4,4,MPI_COMM_WORLD,req3r,ierr)
    !     call MPI_Wait(req3r,sta3r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 do i=1,15
    !                     tmpf(4,i,xi,yi,zi) = tmp_f(i,xi,yi,zi)
    !                 enddo
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 5) then
    !     call MPI_Isend(tmp_f(1,1,1,1),15*x_procs*y_procs*(ymax+1),MPI_REAL8,0,5,MPI_COMM_WORLD,req3s,ierr)
    !     call MPI_Wait(req3s,sta3s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp_f(1,1,1,1),15*x_procs*y_procs*(ymax+1),MPI_REAL8,5,5,MPI_COMM_WORLD,req3r,ierr)
    !     call MPI_Wait(req3r,sta3r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 do i=1,15
    !                     tmpf(5,i,xi,yi,zi) = tmp_f(i,xi,yi,zi)
    !                 enddo
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 6) then
    !     call MPI_Isend(tmp_f(1,1,1,1),15*x_procs*y_procs*(ymax+1),MPI_REAL8,0,6,MPI_COMM_WORLD,req3s,ierr)
    !     call MPI_Wait(req3s,sta3s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp_f(1,1,1,1),15*x_procs*y_procs*(ymax+1),MPI_REAL8,6,6,MPI_COMM_WORLD,req3r,ierr)
    !     call MPI_Wait(req3r,sta3r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 do i=1,15
    !                     tmpf(6,i,xi,yi,zi) = tmp_f(i,xi,yi,zi)
    !                 enddo
    !             enddo
    !         enddo
    !     enddo
    ! endif

    ! if(comm_rank == 7) then
    !     call MPI_Isend(tmp_f(1,1,1,1),15*x_procs*y_procs*(ymax+1),MPI_REAL8,0,7,MPI_COMM_WORLD,req3s,ierr)
    !     call MPI_Wait(req3s,sta3s,ierr)
    ! endif
    ! if(comm_rank == 0) then
    !     call MPI_Irecv(tmp_f(1,1,1,1),15*x_procs*y_procs*(ymax+1),MPI_REAL8,7,7,MPI_COMM_WORLD,req3r,ierr)
    !     call MPI_Wait(req3r,sta3r,ierr)

    !     do zi=1,y_procs
    !         do yi=1,ymax+1
    !             do xi=1,x_procs
    !                 do i=1,15
    !                     tmpf(7,i,xi,yi,zi) = tmp_f(i,xi,yi,zi)
    !                 enddo
    !             enddo
    !         enddo
    !     enddo
    ! endif

    if(comm_rank == 0) then
        n = 0
        do Nyy=0,Ny-1
            do Nxx=0,Nx-1
                do zi=1,y_procs
                    do yi=1,ymax+1
                        do xi=1,x_procs
                            u1((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs)) = tmp1(n,xi,yi,zi)
                            u2((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs)) = tmp2(n,xi,yi,zi)
                            u3((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs)) = tmp3(n,xi,yi,zi)
                            p((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs)) = tmp4(n,xi,yi,zi)
                            do i=1,15
                                f(i,(xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nyy*(y_procs)) = tmpf(n,i,xi,yi,zi)
                            enddo
                        enddo
                    enddo
                enddo
                n = n + 1
            enddo
        enddo
        write(chmyrank, '(i2.2)') comm_rank   !myrankを文字型変数に格納
        open(61, file = './u1_'//chmyrank//'.d')    !ランク番号が付いたファイルを作成
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    write(61, "(7es23.16)") dble(xi),dble(yi),dble(zi),u1(xi,yi,zi),u2(xi,yi,zi),u3(xi,yi,zi),p(xi,yi,zi)
                enddo
            enddo
        enddo
        close(61)

        write(chmyrank, '(i2.2)') comm_rank   !myrankを文字型変数に格納
        open(71, file = 'data2.d')    !ランク番号が付いたファイルを作成
        do zi=0,zmax
            ! do yi=0,ymax
                do xi=0,xmax
                    ! do i=1,15
                        write(71,"(1es24.17)") u1(xi,0,zi)
                    ! enddo
                enddo
            ! enddo
        enddo
        close(71)
    endif

    call MPI_Finalize(ierr)

    call cpu_time(time2)
    open(72,file="time.d")
    write(72,*) time2-time1
    close(72)
end program main
