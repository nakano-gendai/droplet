! 【Multi-times simulations】
! このプログラムは、多数回の液滴分裂の数値シミュレーションを実行する。
! このプログラムは、界面を「Cahn-Hilliard 方程式」で表現している。
! 液滴を長時間保持できないが、系全体の保存性には優れている。
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module globals
    include "mpif.h"
    real(8),parameter:: ds = 1.0d0 !格子間隔（lattice unit）
    integer,parameter:: xmax = 3 !ｘ方向格子数（０から数える）
    integer,parameter:: ymax = 3 !ｙ方向格子数（０から数える）
    integer,parameter:: zmax = 3 !ｚ方向格子数（０から数える）
    integer,parameter:: Nxall = 3 !x方向の分割数（Nxall*Nzall=全体の計算数）
    integer,parameter:: Nzall = 1 !z方向の分割数（Nxall*Nzall=全体の計算数）
    integer,parameter:: xall = (xmax + 1) * Nxall !全体のx方向格子数
    integer,parameter:: zall = (zmax + 1) * Nzall !全体のz方向格子数
    integer,parameter:: step = 1 !計算時間step
    integer,parameter:: start = 5000 !壁を動かし始める時間step
    integer,parameter:: dirstep = 100 !ディレクトリーを作成するstep
    integer i, j, k, n, nall, nstar, xi, yi, zi, alpha, beta
    character :: filename*200
    character :: filename2*200
    character :: filename3*200
    ! character(*),parameter :: datadir = "/data/sht/nakanog/newinisim/"
    ! character(*),parameter :: datadir2 = "/data/sht/nakanog/newinisim/fg/"
    character(*),parameter :: datadir = "/data/n/n517/re500we7_box_d-4/"
    character(*),parameter :: datadir2 = "/data/n/n517/re500we7_box_d-4/fg/"
    character(*),parameter :: datadir3 = "/data/n/n517/ini_0.05_999_290_290/fg/"
    ! character(*),parameter :: datadir_ran = "/data/n/n517/rantest/"
    ! character(*),parameter :: datadir_ran = "/data/sht/nakanog/rantest/"
    character(*),parameter :: datadir_ran = "/data/2022/nakano/rantest/"

    !無次元数
    real(8),parameter:: We = 7.0d0 !粒子ウエーバー数
    real(8),parameter:: Re = 500.0d0 !粒子レイノルズ数
    real(8),parameter:: eta = 1.0d0 !粘度比(nu2/nu1)

    !支配パラメータ
    real(8) D !設置する液滴径
    real(8) uw !壁の移動速度
    real(8) nu1  !連続相の動粘度
    real(8) nu2  !分散相の動粘度
    real(8) sigma !界面張力
    real(8) kappag  !界面張力を決めるパラメータ
    real(8),parameter:: kappaf = 0.01d0*ds**2 !界面厚さを決めるパラメータ
    real(8),parameter:: H1 = dble(zmax) !代表長さ
    real(8),parameter:: phi1 = 2.211d0
    real(8),parameter:: phi2 = 4.895d0
    real(8),parameter:: a = 9.0d0/49.0d0
    real(8),parameter:: b = 2.0d0/21.0d0
    real(8),parameter:: T = 0.55d0
    real(8),parameter:: tauf = 0.7d0
    real(8),parameter:: Anu = 0.1d0

    !初期液滴中心
    real(8), parameter:: xc = 0.5d0*dble(xmax)	
    real(8), parameter:: yc = 0.5d0*dble(ymax)
    real(8), parameter:: zc = 0.5d0*dble(zmax)

    !粒子速度（整数）
    integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
    integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
    integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)
    real(8):: cr(1:3, 1:15)  !粒子速度（実数）
    real(8) krone(1:3,1:3) !クロネッカーのデルタ

    !パラメータ
    real(8),parameter:: E(15) = (/ 2.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, &
                                1.0d0/9.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, &
                                1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0 /)
    real(8),parameter:: H(15) = (/ 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
                                0.0d0, 0.0d0, 0.0d0 /)
    real(8),parameter:: F(15) = (/ -7.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, &
                                1.d0/24.d0, 1.d0/24.d0, 1.d0/24.d0, 1.d0/24.d0, &
                                1.d0/24.d0, 1.d0/24.d0, 1.d0/24.d0, 1.d0/24.d0 /)
    real(8),parameter:: pi = acos(-1.0d0)

    !その他変数
    real(8) phi_min, phi_max, min, max
    real(8) gtemp, ftemp
    real(8) dif
    integer grobalx, grobaly, grobalz
    real(8) dummy

    !MPI用変数
    integer ierr, comm_procs, comm_rank
    integer Nx, Nz !xi、yi方向の分割数
    integer x_procs, z_procs
    integer key_new, group_new, key_x, group_x, key_z, group_z
    integer new_comm_world, new_procs, new_rank
    integer newx_comm_world, newx_procs, newx_rank
    integer newz_comm_world, newz_procs, newz_rank
    integer rectangle_type,rectangle_type_f,xtype,xtype_f,output_type
    integer next_rank_x,former_rank_x,next_rank_z,former_rank_z
    integer req1s,req1r,req2s,req2r,req3s,req3r,req4s,req4r,req5s,req5r,req6s,req6r,req7s,req7r
    integer, dimension(MPI_STATUS_SIZE) :: sta1s,sta1r,sta2s,sta2r,sta3s,sta3r,sta4s,sta4r,sta5s,sta5r,sta6s,sta6r,sta7s,sta7r
    integer Nzz,Nxx
    integer tags, tagr, recv_rank
    character(2) chmyrank

    !乱数
    real(8), parameter :: ran = 1.0d-4

contains
    subroutine par(cx,cy,cz,cr,krone,U)
        integer, intent(in):: cx(15),cy(15),cz(15)
        real(8), intent(out):: cr(1:3,1:15), krone(1:3,1:3), U(1:3)
        do i=1,15
            cr(1,i) = dble(cx(i))
            cr(2,i) = dble(cy(i))
            cr(3,i) = dble(cz(i))
        enddo
        krone(1,1) = 1.d0; krone(1,2) = 0.d0; krone(1,3) = 0.d0
        krone(2,1) = 0.d0; krone(2,2) = 1.d0; krone(2,3) = 0.d0
        krone(3,1) = 0.d0; krone(3,2) = 0.d0; krone(3,3) = 1.d0
        U(1) = 0.0d0
        U(2) = 0.0d0
        U(3) = 0.0d0
    end subroutine par

    subroutine ini(g_procs,gnext_procs,f_procs,fnext_procs,feq_procs,geq_procs,phi_procs,p_procs,u1_procs,u2_procs,u3_procs,grad_phi_procs,gphi_procs,lap_phi_procs,p0_procs,grad_u_procs,eps,eps2)
        real(8),allocatable :: f_procs(:,:,:,:), fnext_procs(:,:,:,:), feq_procs(:,:,:,:)
        real(8),allocatable :: g_procs(:,:,:,:), gnext_procs(:,:,:,:), geq_procs(:,:,:,:)
        real(8),allocatable :: phi_procs(:,:,:), p_procs(:,:,:)
        real(8),allocatable :: u1_procs(:,:,:), u2_procs(:,:,:), u3_procs(:,:,:)
        real(8),allocatable :: p0_procs(:,:,:), lap_phi_procs(:,:,:), grad_phi_procs(:,:,:,:), gphi_procs(:,:,:,:,:), grad_u_procs(:,:,:,:,:),eps(:,:,:,:),eps2(:,:,:,:)

    !========================並列数・コミュニケータ分割・通信先設定========================
        !配列のallocate(xi:0とx_procs+1がのりしろ)(yi:0とymax+2がのりしろ)(zi:0とz_procs+1がのりしろ)
        Nx = 2 !x方向の並列数（ただし，Nx/=comm_procs）
        Nz = comm_procs / (Nx * Nxall * Nzall) !z方向の並列数
        x_procs = (xmax+1) / Nx
        z_procs = (zmax+1) / Nz

        !各条件でコミュニケータを分割する
        key_new = comm_rank
        group_new = comm_rank / (Nx * Nz)
        call MPI_Comm_Split(MPI_COMM_WORLD,group_new,key_new,new_comm_world,ierr)
        call MPI_Comm_Size(new_comm_world,new_procs,ierr)
        call MPI_Comm_Rank(new_comm_world,new_rank,ierr)
        !z方向にコミュニケータを分割する
        key_z = new_rank
        group_z = new_rank / Nx
        call MPI_Comm_Split(new_comm_world,group_z,key_z,newz_comm_world,ierr)
        call MPI_Comm_Size(newz_comm_world,newz_procs,ierr)
        call MPI_Comm_Rank(newz_comm_world,newz_rank,ierr)

        !x方向にコミュニケータを分割する
        key_x = new_rank
        group_x = mod(new_rank,Nx)
        call MPI_Comm_Split(new_comm_world,group_x,key_x,newx_comm_world,ierr)
        call MPI_Comm_Size(newx_comm_world,newx_procs,ierr)
        call MPI_Comm_Rank(newx_comm_world,newx_rank,ierr)

        !のりしろ境界の通信先設定
        next_rank_x = newx_rank + 1
        former_rank_x = newx_rank - 1
        if(newx_rank == 0) then
            former_rank_x = Nz - 1
        else if(newx_rank == Nz - 1) then
            next_rank_x = 0
        endif
        next_rank_z = newz_rank + 1
        former_rank_z = newz_rank - 1
        if(newz_rank == 0) then
            former_rank_z = Nx - 1
        else if(newz_rank == Nx - 1) then
            next_rank_z = 0
        endif

        !x世界での受け渡しをする際の型作成
        call MPI_Type_Vector(ymax+3,x_procs,x_procs+2,MPI_REAL8,xtype,ierr)
        call MPI_Type_Commit(xtype,ierr)
        call MPI_Type_Vector(ymax+3,15*x_procs,15*(x_procs+2),MPI_REAL8,xtype_f,ierr)
        call MPI_Type_Commit(xtype_f,ierr)

        !z世界での受け渡しをする際の型作成
        call MPI_Type_Vector((z_procs+2)*(ymax+3),1,x_procs+2,MPI_REAL8,rectangle_type,ierr)
        call MPI_Type_Commit(rectangle_type,ierr)
        call MPI_Type_Vector((z_procs+2)*(ymax+3),15,15*(x_procs+2),MPI_REAL8,rectangle_type_f,ierr)
        call MPI_Type_Commit(rectangle_type_f,ierr)
    !===============================================================================================================
        !以下はのりしろ有りの変数
        allocate(phi_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(p_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(f_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(fnext_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(feq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(g_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(gnext_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(geq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(p0_procs(1:x_procs,1:ymax+1,1:z_procs))
        allocate(lap_phi_procs(1:x_procs,1:ymax+1,1:z_procs))
        allocate(grad_phi_procs(1:3,1:x_procs,1:ymax+1,1:z_procs))
        allocate(gphi_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs))
        allocate(grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs))
        allocate(eps(1:15,1:x_procs,1:ymax+1,1:z_procs))
        allocate(eps2(1:15,1:x_procs,1:ymax+1,1:z_procs))

        !初期化
        phi_procs(:,:,:) = 0.0d0
        p_procs(:,:,:) = 0.0d0
        u1_procs(:,:,:) = 0.0d0
        u2_procs(:,:,:) = 0.0d0
        u3_procs(:,:,:) = 0.0d0
        feq_procs(:,:,:,:) = 0.0d0
        f_procs(:,:,:,:) = 0.0d0
        fnext_procs(:,:,:,:) = 0.0d0
        geq_procs(:,:,:,:) = 0.0d0
        g_procs(:,:,:,:) = 0.0d0
        gnext_procs(:,:,:,:) = 0.0d0
        p0_procs(:,:,:) = 0.0d0
        lap_phi_procs(:,:,:) = 0.0d0
        grad_phi_procs(:,:,:,:) = 0.0d0
        gphi_procs(:,:,:,:,:) = 0.0d0
        grad_u_procs(:,:,:,:,:) = 0.0d0
        eps(:,:,:,:) = 0.0d0
        eps2(:,:,:,:) = 0.0d0

        if(group_new == 0) then
            ! uw = 0.05d0
            ! nu1 = 0.1d0
            ! nu2 = 0.1d0
            ! kappag = 3.78d-3
            ! !液滴径をランダムに発生させる（後でやる）
            ! D = 0.0d0

            uw = 0.05d0
            D = 145.0d0
            nu1 = 2.0d0*uw*D*D / (4.0d0*Re*H1)
            nu2 = eta*nu1
            sigma = uw*uw*D*D*D / (2.0d0*We*H1*H1)
            kappag = (sigma/1.7039d0)**(1.0d0/0.9991d0)

        elseif(group_new == 1) then
            uw = 0.05d0
            D = 145.0d0
            nu1 = 2.0d0*uw*D*D / (4.0d0*Re*H1)
            nu2 = eta*nu1
            sigma = uw*uw*D*D*D / (2.0d0*We*H1*H1)
            kappag = (sigma/1.7039d0)**(1.0d0/0.9991d0)
        elseif(group_new == 2) then
            uw = 0.05d0
            D = 145.0d0
            nu1 = 2.0d0*uw*D*D / (4.0d0*Re*H1)
            nu2 = eta*nu1
            sigma = uw*uw*D*D*D / (2.0d0*We*H1*H1)
            kappag = (sigma/1.7039d0)**(1.0d0/0.9991d0)
        ! elseif(group_new == 3) then
        !     uw = 0.02d0
        !     nu1 = 0.4d0
        !     nu2 = 0.4d0
        !     kappag = 2.85d-3
        !     !液滴径をランダムに発生させる（後でやる）
        !     D = 32.0d0
        endif
        open(121,file="para.d")
        write(121,*) uw, D, nu1, nu2, sigma, kappag
        close(121)
        !===================================初期条件の設定================================================
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    phi_procs(xi,yi,zi) = phi1
                    grobalx = (xi-1) + group_x * x_procs
                    grobaly = yi-1
                    grobalz = (zi-1) + group_z * z_procs
                    dif = (dble(grobalx)*ds-xc)**2 + (dble(grobaly)*ds-yc)**2 + (dble(grobalz)*ds-zc)**2
                    if(dif <= (0.5d0*D)**2) then
                        phi_procs(xi,yi,zi) = phi2
                    endif
                enddo
            enddo
        enddo
        call glue(phi_procs)
        call grad_cal(grad_phi_procs,phi_procs)
        call lap_cal(lap_phi_procs,phi_procs)
        call gphi_cal(gphi_procs,grad_phi_procs)
        call p0_cal(p0_procs,phi_procs)
        call grad_u_cal(grad_u_procs,u1_procs,u2_procs,u3_procs)
        call feq_cal(gphi_procs,phi_procs,p0_procs,lap_phi_procs,grad_phi_procs,u1_procs,u2_procs,u3_procs,f_procs)
        call geq_cal(gphi_procs,p_procs,u1_procs,u2_procs,u3_procs,g_procs,cr,grad_u_procs)
    end subroutine ini

    subroutine ini_op(taug_procs,nu_procs,phi_procs,x_m,y_m)
        real(8),intent(in):: phi_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),allocatable :: taug_procs(:,:,:), nu_procs(:,:,:), x_m(:,:,:,:), y_m(:,:,:,:)

        !以下はのりしろ無しの変数
        allocate(taug_procs(0:x_procs+1,0:ymax+2,0:z_procs+1))
        allocate(nu_procs(1:x_procs,1:ymax+1,1:z_procs))
        allocate(x_m(1:15,1:x_procs,1:ymax+1,1:z_procs))
        allocate(y_m(1:15,1:x_procs,1:ymax+1,1:z_procs))

        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    nu_procs(xi,yi,zi) = (phi_procs(xi,yi,zi)-phi1)/(phi2-phi1)*(nu2-nu1) + nu1
                    taug_procs(xi,yi,zi) = 0.5d0 + 3.0d0*nu_procs(xi,yi,zi)/ds + 2.0d0/3.0d0*Anu
                enddo
            enddo
        enddo
    end subroutine ini_op

    subroutine glue(var)
        real(8),intent(inout) :: var(0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
            do xi=1,x_procs
                var(xi,0,zi) = var(xi,ymax+1,zi)
                var(xi,ymax+2,zi) = var(xi,1,zi)
            enddo
        enddo
        !X世界でののりしろ通信
        call MPI_Isend(var(1,0,z_procs),1,xtype,next_rank_x,1,newx_comm_world,req1s,ierr)
        call MPI_Irecv(var(1,0,0),1,xtype,former_rank_x,1,newx_comm_world,req1r,ierr)

        call MPI_Isend(var(1,0,1),1,xtype,former_rank_x,2,newx_comm_world,req2s,ierr)
        call MPI_Irecv(var(1,0,z_procs+1),1,xtype,next_rank_x,2,newx_comm_world,req2r,ierr)

        call MPI_Wait(req1s,sta1s,ierr)
        call MPI_Wait(req1r,sta1r,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
        call MPI_Wait(req2r,sta2r,ierr)
        !Y世界でののりしろ通信
        call MPI_Isend(var(x_procs,0,0),1,rectangle_type,next_rank_z,1,newz_comm_world,req1s,ierr)
        call MPI_Irecv(var(0,0,0),1,rectangle_type,former_rank_z,1,newz_comm_world,req1r,ierr)

        call MPI_Isend(var(1,0,0),1,rectangle_type,former_rank_z,2,newz_comm_world,req2s,ierr)
        call MPI_Irecv(var(x_procs+1,0,0),1,rectangle_type,next_rank_z,2,newz_comm_world,req2r,ierr)

        call MPI_Wait(req1s,sta1s,ierr)
        call MPI_Wait(req1r,sta1r,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
        call MPI_Wait(req2r,sta2r,ierr)
    end subroutine glue

    subroutine glue_f(var_f)
        real(8),intent(inout) :: var_f(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
            do xi=1,x_procs
                do i=1,15
                    var_f(i,xi,0,zi) = var_f(i,xi,ymax+1,zi)
                    var_f(i,xi,ymax+2,zi) = var_f(i,xi,1,zi)
                enddo
            enddo
        enddo
        !X世界でののりしろ通信
        call MPI_Isend(var_f(1,1,0,z_procs),1,xtype_f,next_rank_x,1,newx_comm_world,req1s,ierr)
        call MPI_Irecv(var_f(1,1,0,0),1,xtype_f,former_rank_x,1,newx_comm_world,req1r,ierr)

        call MPI_Isend(var_f(1,1,0,1),1,xtype_f,former_rank_x,2,newx_comm_world,req2s,ierr)
        call MPI_Irecv(var_f(1,1,0,z_procs+1),1,xtype_f,next_rank_x,2,newx_comm_world,req2r,ierr)

        call MPI_Wait(req1s,sta1s,ierr)
        call MPI_Wait(req1r,sta1r,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
        call MPI_Wait(req2r,sta2r,ierr)

        !Y世界でののりしろ通信
        call MPI_Isend(var_f(1,x_procs,0,0),1,rectangle_type_f,next_rank_z,1,newz_comm_world,req1s,ierr)
        call MPI_Irecv(var_f(1,0,0,0),1,rectangle_type_f,former_rank_z,1,newz_comm_world,req1r,ierr)

        call MPI_Isend(var_f(1,1,0,0),1,rectangle_type_f,former_rank_z,2,newz_comm_world,req2s,ierr)
        call MPI_Irecv(var_f(1,x_procs+1,0,0),1,rectangle_type_f,next_rank_z,2,newz_comm_world,req2r,ierr)

        call MPI_Wait(req1s,sta1s,ierr)
        call MPI_Wait(req1r,sta1r,ierr)
        call MPI_Wait(req2s,sta2s,ierr)
        call MPI_Wait(req2r,sta2r,ierr)
    end subroutine glue_f

    subroutine MPI_boundary(phi_procs,u1_procs,u2_procs,u3_procs,f_procs,g_procs,taug_procs)
        real(8),intent(inout) :: phi_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout) :: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout) :: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout) :: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout) :: f_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout) :: g_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout) :: taug_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        call glue(phi_procs)
        call glue(u1_procs)
        call glue(u2_procs)
        call glue(u3_procs)
        call glue(taug_procs)
        call glue_f(f_procs)
        call glue_f(g_procs)
    endsubroutine MPI_boundary

    subroutine MPI_boundary_fg(feq_procs,geq_procs)
        real(8),intent(inout) :: feq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout) :: geq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        call glue_f(feq_procs)
        call glue_f(geq_procs)
    endsubroutine MPI_boundary_fg

    subroutine feq_cal(gphi_procs,phi_procs,p0_procs,lap_phi_procs,grad_phi_procs,u1_procs,u2_procs,u3_procs,feq_procs)
        real(8),intent(in):: gphi_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: phi_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: p0_procs(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: lap_phi_procs(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: grad_phi_procs(1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(out):: feq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        gtemp = 0.0d0
                        do beta=1,3
                            do alpha=1,3
                                gtemp = gtemp + gphi_procs(alpha,beta,xi,yi,zi)*cr(alpha,i)*cr(beta,i)
                            enddo
                        enddo
                        feq_procs(i,xi,yi,zi) = H(i)*phi_procs(xi,yi,zi) &
                                        + F(i)*(p0_procs(xi,yi,zi)-kappaf*phi_procs(xi,yi,zi)*lap_phi_procs(xi,yi,zi) &
                                        -kappaf/6.0d0*(grad_phi_procs(1,xi,yi,zi)**2+grad_phi_procs(2,xi,yi,zi)**2+grad_phi_procs(3,xi,yi,zi)**2)) &
                                        + E(i)*phi_procs(xi,yi,zi)*3.0d0*(cr(1,i)*u1_procs(xi,yi,zi)+cr(2,i)*u2_procs(xi,yi,zi)+cr(3,i)*u3_procs(xi,yi,zi)) &
                                        + E(i)*kappaf*gtemp
                    enddo
                enddo
            enddo
        enddo
    end subroutine feq_cal

    subroutine geq_cal(gphi_procs,p_procs,u1_procs,u2_procs,u3_procs,geq_procs,cr,grad_u_procs)
        real(8),intent(in):: gphi_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: p_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(out):: geq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: cr(1:3,1:15)
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        gtemp = 0.0d0
                        ftemp = 0.0d0
                        do beta=1,3
                            do alpha=1,3
                                gtemp = gtemp + gphi_procs(alpha,beta,xi,yi,zi)*cr(alpha,i)*cr(beta,i)

                                ftemp = ftemp + (grad_u_procs(beta,alpha,xi,yi,zi)+grad_u_procs(alpha,beta,xi,yi,zi))*cr(alpha,i)*cr(beta,i)
                            enddo
                        enddo
                        geq_procs(i,xi,yi,zi) = E(i)*(3.0d0*p_procs(xi,yi,zi) + 3.0d0*(u1_procs(xi,yi,zi)*dble(cx(i)) &
                        + u2_procs(xi,yi,zi)*dble(cy(i)) + u3_procs(xi,yi,zi)*dble(cz(i))) &
                        + 4.5d0*(u1_procs(xi,yi,zi)*dble(cx(i)) + u2_procs(xi,yi,zi)*dble(cy(i)) + u3_procs(xi,yi,zi)*dble(cz(i)))**2 &
                        - 1.5d0*(u1_procs(xi,yi,zi)**2 + u2_procs(xi,yi,zi)**2 + u3_procs(xi,yi,zi)**2) &
                        + Anu*ds*ftemp) &
                        + E(i)*kappag*gtemp
                    enddo
                enddo
            enddo
        enddo
    end subroutine geq_cal

    subroutine f_cal(fnext_procs,f_procs,feq_procs)
        real(8),intent(inout):: fnext_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: f_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: feq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        fnext_procs(i,xi,yi,zi) = f_procs(i,xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                - (f_procs(i,xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                - feq_procs(i,xi-cx(i),yi-cy(i),zi-cz(i))) / tauf
                    enddo
                enddo
            enddo
        enddo
    end subroutine f_cal

    subroutine g_cal(gnext_procs,g_procs,geq_procs,taug_procs)
        real(8),intent(inout):: gnext_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: g_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: geq_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: taug_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        gnext_procs(i,xi,yi,zi) = g_procs(i,xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                - (g_procs(i,xi-cx(i),yi-cy(i),zi-cz(i)) &
                                                - geq_procs(i,xi-cx(i),yi-cy(i),zi-cz(i))) / taug_procs(xi-cx(i),yi-cy(i),zi-cz(i))
                    enddo
                enddo
            enddo
        enddo
    end subroutine g_cal

    subroutine bounce_back_LBM(fun,U)
        real(8),intent(inout):: fun(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: U(1:3)
        !安定したら壁を動かす
        if((n >= start)) then
            U(1) = uw
            U(2) = 0.0d0
            U(3) = 0.0d0
        endif
        if(group_z == 0) then
            zi = 1
            do yi=1,ymax+1
                do xi=1,x_procs
                    fun(4,xi,yi,zi) = fun(7,xi,yi,zi) - 6.0d0 * E(7) * ((-U(1))*cr(1,7)+(-U(2))*cr(2,7)+(-U(3))*cr(3,7))
                    fun(8,xi,yi,zi) = fun(12,xi,yi,zi) - 6.0d0 * E(12) * ((-U(1))*cr(1,12)+(-U(2))*cr(2,12)+(-U(3))*cr(3,12))
                    fun(9,xi,yi,zi) = fun(13,xi,yi,zi) - 6.0d0 * E(13) * ((-U(1))*cr(1,13)+(-U(2))*cr(2,13)+(-U(3))*cr(3,13))
                    fun(10,xi,yi,zi) = fun(14,xi,yi,zi) - 6.0d0 * E(14) * ((-U(1))*cr(1,14)+(-U(2))*cr(2,14)+(-U(3))*cr(3,14))
                    fun(15,xi,yi,zi) = fun(11,xi,yi,zi) - 6.0d0 * E(11) * ((-U(1))*cr(1,11)+(-U(2))*cr(2,11)+(-U(3))*cr(3,11))
                enddo
            enddo
        else if(group_z == Nz-1) then
            zi = z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    fun(7,xi,yi,zi) = fun(4,xi,yi,zi) - 6.0d0 * E(4) * ((U(1))*cr(1,4)+(U(2))*cr(2,4)+(U(3))*cr(3,4))
                    fun(11,xi,yi,zi) = fun(15,xi,yi,zi) - 6.0d0 * E(15) * ((U(1))*cr(1,15)+(U(2))*cr(2,15)+(U(3))*cr(3,15))
                    fun(12,xi,yi,zi) = fun(8,xi,yi,zi) - 6.0d0 * E(8) * ((U(1))*cr(1,8)+(U(2))*cr(2,8)+(U(3))*cr(3,8))
                    fun(13,xi,yi,zi) = fun(9,xi,yi,zi) - 6.0d0 * E(9) * ((U(1))*cr(1,9)+(U(2))*cr(2,9)+(U(3))*cr(3,9))
                    fun(14,xi,yi,zi) = fun(10,xi,yi,zi) - 6.0d0 * E(10) * ((U(1))*cr(1,10)+(U(2))*cr(2,10)+(U(3))*cr(3,10))
                enddo
            enddo
        endif
    end subroutine bounce_back_LBM

    subroutine physics(phi_procs,p_procs,u1_procs,u2_procs,u3_procs,f_procs,g_procs,nu_procs,taug_procs)
        real(8),intent(inout):: phi_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: p_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: f_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: g_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: nu_procs(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(inout):: taug_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    phi_procs(xi,yi,zi) = 0.0d0
                    p_procs(xi,yi,zi) = 0.0d0
                    u1_procs(xi,yi,zi) = 0.0d0
                    u2_procs(xi,yi,zi) = 0.0d0
                    u3_procs(xi,yi,zi) = 0.0d0
                    do i=1,15
                        phi_procs(xi,yi,zi) = phi_procs(xi,yi,zi) + f_procs(i,xi,yi,zi)
                        p_procs(xi,yi,zi) = p_procs(xi,yi,zi) + g_procs(i,xi,yi,zi) / 3.0d0
                        u1_procs(xi,yi,zi) = u1_procs(xi,yi,zi) + dble(cx(i))*g_procs(i,xi,yi,zi)
                        u2_procs(xi,yi,zi) = u2_procs(xi,yi,zi) + dble(cy(i))*g_procs(i,xi,yi,zi)
                        u3_procs(xi,yi,zi) = u3_procs(xi,yi,zi) + dble(cz(i))*g_procs(i,xi,yi,zi)
                    enddo
                enddo
            enddo
        enddo
        min = phi2
        max = phi1
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    if(phi_procs(xi,yi,zi) > max) then
                        max = phi_procs(xi,yi,zi)
                    elseif(phi_procs(xi,yi,zi) < min) then
                        min = phi_procs(xi,yi,zi)
                    endif
                enddo
            enddo
        enddo
        call MPI_Allreduce(max,phi_max,1,MPI_REAL8,MPI_MAX,new_comm_world,ierr)
        call MPI_Allreduce(min,phi_min,1,MPI_REAL8,MPI_MIN,new_comm_world,ierr)
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    nu_procs(xi,yi,zi) = (phi_procs(xi,yi,zi)-phi_min)/(phi_max-phi_min)*(nu2-nu1) + nu1
                    taug_procs(xi,yi,zi) = 0.5d0 + 3.0d0*nu_procs(xi,yi,zi)/ds + 2.0d0/3.0d0*Anu
                enddo
            enddo
        enddo
    end subroutine physics

    subroutine renew(fun_procs,funnext_procs)
        real(8),intent(in) :: funnext_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(out) :: fun_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        fun_procs(i,xi,yi,zi) = funnext_procs(i,xi,yi,zi)
                    enddo
                enddo
            enddo
        enddo
    end subroutine renew

    subroutine lap_cal(lap_f,fun)
        real(8),intent(out):: lap_f(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: fun(0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    lap_f(xi,yi,zi) = -14.0d0*fun(xi,yi,zi)
                    do i=2,15
                        lap_f(xi,yi,zi) = lap_f(xi,yi,zi) + fun(xi+cx(i),yi+cy(i),zi+cz(i))
                    enddo
                    lap_f(xi,yi,zi) = lap_f(xi,yi,zi) / (5.0d0*ds**2)
                enddo
            enddo
        enddo
    end subroutine lap_cal

    subroutine grad_cal(grad_f,fun)
        real(8),intent(out):: grad_f(1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: fun(0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do alpha=1,3
                        grad_f(alpha,xi,yi,zi) = 0.0d0
                        do i = 2,15
                            grad_f(alpha,xi,yi,zi) = grad_f(alpha,xi,yi,zi) + cr(alpha,i)*fun(xi+cx(i),yi+cy(i),zi+cz(i))
                        enddo
                        grad_f(alpha,xi,yi,zi) = grad_f(alpha,xi,yi,zi) / (10.0d0*ds)
                    enddo
                enddo
            enddo
        enddo
    end subroutine grad_cal

    subroutine p0_cal(p0_procs,phi_procs)
        real(8),intent(out):: p0_procs(1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: phi_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    p0_procs(xi,yi,zi) = phi_procs(xi,yi,zi)*T/(1.0d0-b*phi_procs(xi,yi,zi)) - a*phi_procs(xi,yi,zi)**2
                enddo
            enddo
        enddo
    end subroutine p0_cal

    subroutine gphi_cal(gphi_procs,grad_phi_procs)
        real(8),intent(out):: gphi_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: grad_phi_procs(1:3,1:x_procs,1:ymax+1,1:z_procs)
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do beta=1,3
                        do alpha=1,3
                            gphi_procs(alpha,beta,xi,yi,zi) = 4.5d0*grad_phi_procs(alpha,xi,yi,zi)*grad_phi_procs(beta,xi,yi,zi) &
                            -1.5d0*(grad_phi_procs(1,xi,yi,zi)**2 + grad_phi_procs(2,xi,yi,zi)**2 + grad_phi_procs(3,xi,yi,zi)**2) &
                            *krone(alpha,beta)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine gphi_cal

    subroutine grad_u_cal(grad_u_procs,u1_procs,u2_procs,u3_procs)
        real(8),intent(out):: grad_u_procs(1:3,1:3,1:x_procs,1:ymax+1,1:z_procs)
        real(8),intent(in):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(in):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do beta=1,3
                        grad_u_procs(1,beta,xi,yi,zi) = 0.0d0
                        grad_u_procs(2,beta,xi,yi,zi) = 0.0d0
                        grad_u_procs(3,beta,xi,yi,zi) = 0.0d0
                        do i=2,15
                            grad_u_procs(1,beta,xi,yi,zi) = grad_u_procs(1,beta,xi,yi,zi)+cr(beta,i)*u1_procs(xi+cx(i),yi+cy(i),zi+cz(i))
                            grad_u_procs(2,beta,xi,yi,zi) = grad_u_procs(2,beta,xi,yi,zi)+cr(beta,i)*u2_procs(xi+cx(i),yi+cy(i),zi+cz(i))
                            grad_u_procs(3,beta,xi,yi,zi) = grad_u_procs(3,beta,xi,yi,zi)+cr(beta,i)*u3_procs(xi+cx(i),yi+cy(i),zi+cz(i))
                        enddo
                        grad_u_procs(1,beta,xi,yi,zi) = grad_u_procs(1,beta,xi,yi,zi)/(10.0d0*ds)
                        grad_u_procs(2,beta,xi,yi,zi) = grad_u_procs(2,beta,xi,yi,zi)/(10.0d0*ds)
                        grad_u_procs(3,beta,xi,yi,zi) = grad_u_procs(3,beta,xi,yi,zi)/(10.0d0*ds)
                    enddo
                enddo
            enddo
        enddo
    end subroutine grad_u_cal

    subroutine output(phi_procs,u1_procs,u2_procs,u3_procs,p_procs)
        real(8),intent(inout):: phi_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: p_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: u1_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: u2_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: u3_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)

        write(filename,*) group_new !i->filename 変換
        write(filename2,*) n
        write(filename3,*) new_rank
        filename=datadir//trim(adjustl(filename))//'_'//trim(adjustl(filename3))//'_'//trim(adjustl(filename2))//'.bin' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
        print *, filename !表示してみる
        open(100,file=filename, form='unformatted',status='replace') 
        ! open(20,file=filename,status='replace') 
        yi = ymax / 2 
        do zi=1,z_procs
            ! do yi=1,ymax+1
                do xi=1,x_procs
                    ! write(20) phi_procs(xi,yi,zi),u1_procs(xi,yi,zi),u2_procs(xi,yi,zi),u3_procs(xi,yi,zi),p_procs(xi,yi,zi)
                    ! write(100) real(xi),real(zi),real(phi_procs(xi,yi,zi)),real(u1_procs(xi,yi,zi))
                    write(100) real(phi_procs(xi,yi,zi)),real(u1_procs(xi,yi,zi)),real(u2_procs(xi,yi,zi)),real(u3_procs(xi,yi,zi)),real(p_procs(xi,yi,zi))
                enddo
                write(100)
            ! enddo
        enddo
        close(100)
    end subroutine output

    subroutine outputphi(phi_procs)
        real(8),intent(inout):: phi_procs(0:x_procs+1,0:ymax+2,0:z_procs+1)

        write(filename,*) group_new !i->filename 変換
        write(filename2,*) n
        write(filename3,*) new_rank
        filename=datadir//trim(adjustl(filename))//'_'//trim(adjustl(filename3))//'_'//trim(adjustl(filename2))//'.bin' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
        print *, filename !表示してみる
        open(101,file=filename, form='unformatted',status='replace') 
        ! open(20,file=filename,status='replace') 
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    write(101) phi_procs(xi,yi,zi)
                enddo
            enddo
        enddo
        close(101)
    end subroutine outputphi

    subroutine outputfg(f_procs, g_procs)
        real(8),intent(inout):: f_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)
        real(8),intent(inout):: g_procs(1:15,0:x_procs+1,0:ymax+2,0:z_procs+1)

        write(filename,*) group_new !i->filename 変換
        write(filename2,*) n
        write(filename3,*) new_rank
        filename=datadir2//trim(adjustl(filename))//'_'//trim(adjustl(filename3))//'_'//trim(adjustl(filename2))//'_fg.bin' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
        print *, filename !表示してみる
        open(102,file=filename, form='unformatted',status='replace') 
        ! open(20,file=filename,status='replace') 
        ! yi = ymax / 2 
        do zi=1,z_procs
            do yi=1,ymax+1
                do xi=1,x_procs
                    do i=1,15
                        write(102) f_procs(i,xi,yi,zi), g_procs(i,xi,yi,zi)
                        ! write(102) g_procs(i,xi,yi,zi)
                    enddo
                enddo
            enddo
        enddo
        close(102)
    end subroutine outputfg
    
    !ディレクトリ作成
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
    real(8) U(1:3) !壁の移動速度

    real(8),allocatable :: f_procs(:,:,:,:), fnext_procs(:,:,:,:), feq_procs(:,:,:,:)
    real(8),allocatable :: g_procs(:,:,:,:), gnext_procs(:,:,:,:), geq_procs(:,:,:,:)
    real(8),allocatable :: phi_procs(:,:,:), p_procs(:,:,:)
    real(8),allocatable :: u1_procs(:,:,:), u2_procs(:,:,:), u3_procs(:,:,:)

    real(8),allocatable :: taug_procs(:,:,:), nu_procs(:,:,:)
    real(8),allocatable :: p0_procs(:,:,:), lap_phi_procs(:,:,:), grad_phi_procs(:,:,:,:), gphi_procs(:,:,:,:,:)
    real(8),allocatable :: grad_u_procs(:,:,:,:,:), eps(:,:,:,:), eps2(:,:,:,:), x_m(:,:,:,:), y_m(:,:,:,:)

    real(8) time1,time2

    integer seedsize
    integer, allocatable :: seeds(:), seeds2(:)
    ! character(8) date
    ! character(10) time
    ! character(5) zone
    ! integer values(8)
    ! integer time
    integer seed, seed2
    
    call cpu_time(time1)

    call MPI_Init(ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, comm_procs, ierr)
    call MPI_Comm_Rank(MPI_COMM_WORLD, comm_rank, ierr)

    ! call mk_dirs(datadir)
    ! call mk_dirs(datadir2)
    call mk_dirs(datadir_ran)
    call par(cx,cy,cz,cr,krone,U)
    call ini(g_procs,gnext_procs,f_procs,fnext_procs,feq_procs,geq_procs,phi_procs,p_procs,u1_procs,u2_procs,u3_procs,grad_phi_procs,gphi_procs,lap_phi_procs,p0_procs,grad_u_procs,eps,eps2)
    call ini_op(taug_procs,nu_procs,phi_procs,x_m,y_m)

    call random_seed(size=seedsize)
    allocate(seeds(1:seedsize))
    call system_clock(seed)
    seed = seed + comm_rank
    do i=1,seedsize
        seeds(i) = seed + i
    enddo
    call random_seed(put=seeds)
    call random_number(eps)

    allocate(seeds2(1:seedsize))
    call system_clock(seed2)
    seed2 = seed + seed2 + (2 * comm_rank)
    do i=1,seedsize
        seeds2(i) = seed2 + i
    enddo
    call random_seed(put=seeds2)
    call random_number(eps2)

    do zi=1,z_procs
        do yi=1,ymax+1
            do xi=1,x_procs
                do i=1,15
                    x_m(i,xi,yi,zi) = sqrt(-2.0d0*log(eps(i,xi,yi,zi)))*cos(2.0d0*pi*eps2(i,xi,yi,zi))
                    y_m(i,xi,yi,zi) = sqrt(-2.0d0*log(eps(i,xi,yi,zi)))*sin(2.0d0*pi*eps2(i,xi,yi,zi))
                enddo
            enddo
        enddo
    enddo

    write(filename,*) group_new !i->filename 変換
    write(filename2,*) n
    write(filename3,*) new_rank
    filename=datadir_ran//trim(adjustl(filename))//'_'//trim(adjustl(filename3))//'_'//trim(adjustl(filename2))//'_ran.d' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
    print *, filename !表示してみる
    open(150,file=filename, form='formatted',status='replace') 
    do zi=1,z_procs
        do yi=1,ymax+1
            do xi=1,x_procs
                do i=1,15
                    write(150,*) x_m(i,xi,yi,zi), y_m(i,xi,yi,zi)
                enddo
            enddo
        enddo
    enddo
    close(150)


! DO n=1,step
!     if(n == 5000) then
!         !一様乱数の発生
!         ! call system_clock(seed)
!         ! seed = seed + comm_rank
!         ! seeds(1) = seed
!         ! call random_seed(put=seeds)
!         ! call random_number(eps)
!         ! eps(:,:,:,:) = eps(:,:,:,:) - 0.5d0
!         ! eps(:,:,:,:) = eps(:,:,:,:) * 2.0d0
!         ! eps(:,:,:,:) = eps(:,:,:,:) * 1.0d-2

!         write(filename,*) group_new !i->filename 変換
!         write(filename3,*) new_rank
!         filename=datadir3//'0_'//trim(adjustl(filename3))//'_200000_fg.bin' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
!         print *, filename !表示してみる
!         open(103, file=filename, form="unformatted")
!         do zi=1,z_procs
!             do yi=1,ymax+1
!                 do xi=1,x_procs
!                     do i=1,15
!                         ! read(103) dummy, g_procs(i,xi,yi,zi)
!                         read(103) g_procs(i,xi,yi,zi)
!                         ! g_procs(i,xi,yi,zi) = g_procs(i,xi,yi,zi) + eps(i,xi,yi,zi)
!                         x_m(i,xi,yi,zi) = x_m(i,xi,yi,zi) * ran
!                         g_procs(i,xi,yi,zi) = g_procs(i,xi,yi,zi) + x_m(i,xi,yi,zi)
!                     enddo
!                 enddo
!             enddo
!         enddo
!         close(103)

!         do zi=1,z_procs
!             do yi=1,ymax+1
!                 do xi=1,x_procs
!                     p_procs(xi,yi,zi) = 0.0d0
!                     u1_procs(xi,yi,zi) = 0.0d0
!                     u2_procs(xi,yi,zi) = 0.0d0
!                     u3_procs(xi,yi,zi) = 0.0d0
!                     do i=1,15
!                         p_procs(xi,yi,zi) = p_procs(xi,yi,zi) + g_procs(i,xi,yi,zi) / 3.0d0
!                         u1_procs(xi,yi,zi) = u1_procs(xi,yi,zi) + dble(cx(i))*g_procs(i,xi,yi,zi)
!                         u2_procs(xi,yi,zi) = u2_procs(xi,yi,zi) + dble(cy(i))*g_procs(i,xi,yi,zi)
!                         u3_procs(xi,yi,zi) = u3_procs(xi,yi,zi) + dble(cz(i))*g_procs(i,xi,yi,zi)
!                     enddo
!                 enddo
!             enddo
!         enddo
!     endif
! !======================のりしろ境界==========================================================
!     call MPI_boundary(phi_procs,u1_procs,u2_procs,u3_procs,f_procs,g_procs,taug_procs)
! !=======================div/gradの計算============================================================
!     call lap_cal(lap_phi_procs,phi_procs)
!     call grad_cal(grad_phi_procs,phi_procs)
!     call p0_cal(p0_procs,phi_procs)
!     call gphi_cal(gphi_procs,grad_phi_procs)
!     call grad_u_cal(grad_u_procs,u1_procs,u2_procs,u3_procs)
! !=======================局所平衡分布関数の計算========================================================
!     call feq_cal(gphi_procs,phi_procs,p0_procs,lap_phi_procs,grad_phi_procs,u1_procs,u2_procs,u3_procs,feq_procs)
!     call geq_cal(gphi_procs,p_procs,u1_procs,u2_procs,u3_procs,geq_procs,cr,grad_u_procs)
!     call MPI_boundary_fg(feq_procs,geq_procs)
! !=======================平衡分布関数gの計算=============================================================
!     call f_cal(fnext_procs,f_procs,feq_procs)
!     call g_cal(gnext_procs,g_procs,geq_procs,taug_procs)
! !=======================境界条件=============================================================
!     call bounce_back_LBM(gnext_procs,U)
!     call renew(f_procs,fnext_procs)
!     call renew(g_procs,gnext_procs)
! !==========================物理量の計算==========================================================
!     call physics(phi_procs,p_procs,u1_procs,u2_procs,u3_procs,f_procs,g_procs,nu_procs,taug_procs)

!     ! if((mod(n,200)==0) .and. (n >= 5000)) then
!     !     call output(phi_procs,u1_procs,u2_procs,u3_procs,p_procs)
!     ! endif
!     if((mod(n,3000)==0)) then
!         call output(phi_procs,u1_procs,u2_procs,u3_procs,p_procs)
!     endif
!     ! if(n==4999) then
!     !     call outputphi(phi_procs)
!     ! endif
!     ! if((n == 4999).or.(mod(n,100000)==0)) then
!     !     call outputfg(f_procs,g_procs)
!     ! endif
!     if(mod(n,10000)==0) then
!         call outputfg(f_procs,g_procs)
!     endif
!     ! if(n == 300000) then
!     !     call outputfg(f_procs,g_procs)
!     ! endif
!     if(comm_rank == 0) then
!         call cpu_time(time2)
!         if(mod(n,100)==0) then
!             if(n == 1000) then
!                 open(10,file="./time.d")
!                 write(10,*) n, time2-time1, u1_procs(1,1,1)
!                 close(10)
!             else
!                 open(10,file="./time.d",action="write",position="append")
!                 write(10,*) n, time2-time1, u1_procs(1,1,1)
!                 close(10)
!             endif
!         endif
!     endif
! ENDDO
    call MPI_Finalize(ierr)
end program main