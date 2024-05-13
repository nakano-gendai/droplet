module subprogs
    use decomp_2d
    use decomp_2d_fft
    use decomp_2d_io
    use glassman
    implicit none
    include 'mpif.h'
    include 'fftw3.f'

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !-------------------------------------パラメータ設定-------------------------------------------!
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    double precision, parameter :: pi = 3.14159265358979d0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!-----------------------領域・格子点数の設定-------------------------!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    integer, parameter          :: NX = 64, NY = 64, NZ = 64   
    double precision, parameter :: Xmax = 2.0d0*pi, Ymax =  2.0d0*pi, Zmax =  2.0d0*pi
    double precision, parameter :: dX = Xmax / dble(NX), dY = Ymax / dble(NY), dZ = Zmax / dble(NZ) 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!------------------------------MPI/2decomp---------------------------!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, parameter                  :: MPI_Y_NUM = 8, MPI_Z_NUM = 8                          !!!Y方向・Z方向の並列数
    integer, parameter                  :: NY_procs = NY / MPI_Y_NUM                               !!!各ランクにおけるY方向の格子点数
    integer, parameter                  :: NZ_procs = NZ / MPI_Z_NUM                               !!!各ランクにおけるZ方向の格子点数
    integer, save                       :: NY_procs_min, NZ_procs_min, NY_procs_max, NZ_procs_max  
    integer, save                       :: next_ranky, former_ranky, next_rankz, former_rankz 
    integer, save                       :: ierr, procs, myrank, req_send, req_recv, myranky, myrankz
    integer, dimension(MPI_STATUS_SIZE) :: sta_send, sta_recv
    integer, save                       :: boundary_y                                               
    !!!↑MPI用, ↓2decomp用
    integer, save                       :: sta(1:3), end(1:3), size(1:3)                            !!!2decomp用変数

    double precision, save  :: X(-1:NX+1), Y(-1: NY_procs+1), Z(-1: NZ_procs+1)

    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!----------------------------MPIの設定-----------------------------------!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine MPI_settings

        myranky = 0
        myrankz = 0

        !!!Y方向・Z方向のMPI並列数の確認
        if(myrank == 0)then
            write(*,*) 'Y方向並列数 : ', MPI_Y_NUM
            write(*,*) 'Z方向並列数 : ', MPI_Z_NUM
        end if

        !!!Y方向とZ方向の並列番号
        myranky = myrank / MPI_Z_NUM
        myrankz = myrank - myranky * MPI_Z_NUM

        !!!各プロセスの最小と最大番号
        !!!外力点の補間時に使用（0〜N_procsまで）
            NZ_procs_min = myrankz*NZ_procs
            NZ_procs_max = myrankz*NZ_procs + NZ_procs
            NY_procs_min = myranky*NY_procs
            NY_procs_max = myranky*NY_procs + NY_procs 

        !!!Y方向プロセスの通信番号の取得
        if(myrank >= 0 .and. myrank <= MPI_Z_NUM - 1) then
            former_ranky = myrank + MPI_Z_NUM*(MPI_Y_NUM-1)
        else
            former_ranky = myrank - MPI_Z_NUM
        endif
        if(myrank >= procs - MPI_Z_NUM .and. myrank <= procs -1) then
            next_ranky = myrank - MPI_Z_NUM*(MPI_Y_NUM-1)
        else
            next_ranky = myrank + MPI_Z_NUM
        end if

        !!!Z方向プロセスの通信番号の取得
        if(mod(myrank,MPI_Z_NUM) == 0) then
            former_rankz = myrank + (MPI_Z_NUM - 1)
        else
            former_rankz = myrank - 1
        endif
        if(mod(myrank,MPI_Z_NUM) == MPI_Z_NUM - 1) then
            next_rankz = myrank - (MPI_Z_NUM - 1)
        else
            next_rankz = myrank + 1
        end if

        !!!通信端末の確認
        ! write(*,*)'myrank, next_ranky, former_ranky, next_rankz, former_rankz',myrank, next_ranky, former_ranky, next_rankz, former_rankz

        !!!通信用データタイプ型の作成
        call MPI_Type_vector(NZ_procs+2, NX+2, (NX+2)*(NY_procs+2), MPI_DOUBLE_PRECISION, boundary_y, ierr)
        call MPI_Type_Commit(boundary_y, ierr)

    end subroutine MPI_settings

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!----------------------------格子・外力点の初期設定----------------------!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine lattice_init_settings
        integer i, j, k

        !!!各ランクの計算領域を設定
        do i = -1, NX + 1
            X(i) = dble(i) * dX
        enddo
        do j = -1, NY_procs + 1
            Y(j) = myranky * NY_procs * dY + dble(j) * dY
        enddo
        do k = -1, NZ_procs + 1 
            Z(k) = myrankz * NZ_procs * dZ + dble(k) * dZ
        enddo

    end subroutine lattice_init_settings

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!----------------------------物理量の境界条件----------------------------!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine boundary(in)
        double precision, intent(inout) :: in(0:NX+1, 0:NY_procs+1,0:NZ_procs+1)

        !!!X方向の周期境界
        in(NX+1, 1:NY_procs, 1:NZ_procs) = in(1, 1:NY_procs, 1:NZ_procs)
        in(0, 1:NY_procs, 1:NZ_procs) = in(NX, 1:NY_procs, 1:NZ_procs)

        !!!Z方向の通信(1 to NZ_procs+1)
        call MPI_Isend(in(0,0,1), (NX+2)*(NY_procs+2), MPI_DOUBLE_PRECISION,former_rankz, 0, MPI_COMM_WORLD, req_send, ierr)
        call MPI_Irecv(in(0,0,NZ_procs+1), (NX+2)*(NY_procs+2),MPI_DOUBLE_PRECISION, next_rankz, 0, MPI_COMM_WORLD, req_recv, ierr)
        call MPI_Wait(req_send, sta_send, ierr)
        call MPI_Wait(req_recv, sta_recv, ierr)

        !!!Z方向の通信(NZ_procs to 0)
        call MPI_Isend(in(0,0,NZ_procs), (NX+2)*(NY_procs+2),MPI_DOUBLE_PRECISION, next_rankz, 0, MPI_COMM_WORLD, req_send, ierr)
        call MPI_Irecv(in(0,0,0), (NX+2)*(NY_procs+2), MPI_DOUBLE_PRECISION,former_rankz, 0, MPI_COMM_WORLD, req_recv, ierr)
        call MPI_Wait(req_send, sta_send, ierr)
        call MPI_Wait(req_recv, sta_recv, ierr)

        !!!Y方向の通信(1 to NY_procs+1)
        call MPI_Isend(in(0,1,0), 1, boundary_y, former_ranky, 0,MPI_COMM_WORLD, req_send, ierr)
        call MPI_Irecv(in(0,NY_procs+1,0), 1, boundary_y, next_ranky, 0,MPI_COMM_WORLD, req_recv, ierr)
        call MPI_Wait(req_send, sta_send, ierr)
        call MPI_Wait(req_recv, sta_recv, ierr)

        !!!Y方向の通信(NY_procs to 0)
        call MPI_Isend(in(0,NY_procs,0), 1, boundary_y, next_ranky, 0,MPI_COMM_WORLD, req_send, ierr)
        call MPI_Irecv(in(0,0,0), 1, boundary_y, former_ranky, 0,MPI_COMM_WORLD, req_recv, ierr)
        call MPI_Wait(req_send, sta_send, ierr)
        call MPI_Wait(req_recv, sta_recv, ierr)

    end subroutine boundary

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!----------------------------実数をフーリエ変換するサブルーチン----------------------------!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine fft_r2c(q, f_q)
        real(8), intent(inout)            :: q(1 : NX, 1 : NY_procs, 1 : NZ_procs)
        complex(kind(0d0)), intent(inout) :: f_q(sta(1) : end(1), sta(2) : end(2), sta(3) : end(3))
        call decomp_2d_fft_3d(q, f_q)
        f_q = f_q / (NX * NY * NZ)
    end subroutine fft_r2c

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!----------------------------------複素数を逆フーリエ変換----------------------------------!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine fft_c2r(q, f_q)
        real(8), intent(inout)            :: q(1 : NX, 1 : NY_procs, 1 : NZ_procs)
        complex(kind(0d0)), intent(inout) :: f_q(sta(1) : end(1), sta(2) : end(2), sta(3) : end(3))
        call decomp_2d_fft_3d(f_q, q)
    end subroutine fft_c2r

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!-------------------------N/2より大きかったら-Nを返す関数（FFT用）-------------------------!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function judge(i, N)
        integer :: i, N, judge
        if(i >= N / 2 + 1)then
            judge = N
        else
            judge = 0
        endif
    end function

end module

    program main
        use subprogs
        use decomp_2d
        use decomp_2d_fft
        use decomp_2d_io
        use glassman
        implicit none
        double precision    :: in(1 : NX, 1 : NY_procs, 1 : NZ_procs)
        complex(kind(0d0)), allocatable :: out(:,:,:)
        !///////////////////////////////////////////////////////////////
        !!!ポアソン方程式の右辺
        double precision Psi(1:NX, 1:NY_procs, 1:NZ_procs)
        !!!ポアソン方程式のPhi
        double precision Phi(1:NX, 1:NY_procs, 1:NZ_procs)
        !!!フーリエ変換された右辺
        complex(kind(0d0)),allocatable  :: f_Psi(:,:,:)
        !!!フーリエ変換された圧力ポテンシャル
        complex(kind(0d0)), allocatable :: f_Phi(:,:,:)
        double precision kx, ky, kz
        integer i,j,k
        double precision error_poisson, error_poisson_all
    
        !!!MPI開始
        call MPI_Init(ierr)
        call MPI_Comm_Size(MPI_COMM_WORLD, procs, ierr)
        call MPI_Comm_Rank(MPI_COMM_WORLD, myrank, ierr)

        !!!MPI変数の設定
        call MPI_settings
        !!!格子の設定
        call lattice_init_settings

        !!!2decompe開始
        call decomp_2d_init(NX, NY, NZ, MPI_Y_NUM, MPI_Z_NUM)
        call decomp_2d_fft_init(PHYSICAL_IN_X)
        call decomp_2d_fft_get_size(sta, end, size)

        allocate(out(sta(1) : end(1), sta(2) : end(2), sta(3) : end(3)))
        allocate(f_Psi(sta(1) : end(1), sta(2) : end(2), sta(3) : end(3)))
        allocate(f_Phi(sta(1) : end(1), sta(2) : end(2), sta(3) : end(3)))

        ! !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! !!!fft実行r2c(in：実数，out：複素数)
        ! call fft_r2c(in, out)
        ! !!!fft実行c2r(in：複素数，out：複素数)
        ! call fft_r2c(out, in)
        ! !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        !!!ポアソン方程式の右辺にPsi = sin2xsin2ysin2zを代入
        !!!dPhi/dx + dPhi/dy + dPhi/dz = Psi
        do k = 1, NZ_procs
        do j = 1, NY_procs
        do i = 1, NX
            Psi(i,j,k) = sin(2.0d0*(X(i)+X(i-1))/2.0d0) &
            *sin(2.0d0*(Y(j)+Y(j-1))/2.0d0)*sin(2.0d0*(Z(k)+Z(k-1))/2.0d0)
        end do
        end do
        end do

        !!!右辺のフーリエ変換（Psi ⇒ f_Psi）
        call fft_r2c(Psi, f_Psi)

        !!!ポアソン方程式を解く
        do k = sta(3), end(3)
        do j = sta(2), end(2)
        do i = sta(1), end(1)
            ! kx = (2.0d0*pi/Xmax)*((i - 1)-judge(i-1,NX))
            ! ky = (2.0d0*pi/Ymax)*((j - 1)-judge(j-1,NY))
            ! kz = (2.0d0*pi/Zmax)*((k - 1)-judge(k-1,NZ))
            kx = dble((2.0d0*pi/Xmax)*((i - 1)))
            ky = dble((2.0d0*pi/Ymax)*((j - 1)))
            kz = dble((2.0d0*pi/Zmax)*((k - 1)))
            if(kx == 0.0d0 .and. ky == 0.0d0 .and. kz == 0.0d0)then
                !!!!!基準点を与える
                f_Phi(i, j, k) = (0.0d0,0.0d0)
            else
                f_Phi(i, j, k) = f_Psi(i,j,k) / (-2.0d0 *  ((1.0d0 - cos(kx * dx)) / (dx ** 2) + &
                                                            (1.0d0 - cos(ky * dy)) / (dy ** 2) + &
                                                            (1.0d0 - cos(kz * dz)) / (dz ** 2)))
            endif
        enddo
        enddo
        enddo

        !!!フーリエ逆変換（f_Phi ⇒ Phi_2decomp）
        call fft_c2r(Phi,f_Phi)

        error_poisson = 0.0d0
        do k = 1, NZ_procs
        do j = 1, NY_procs
        do i = 1, NX
            error_poisson = error_poisson + (-(1.0d0/12.0d0)*sin(2.0d0*(X(i)+X(i-1))/2.0d0) &
            *sin(2.0d0*(Y(j)+Y(j-1))/2.0d0)*sin(2.0d0*(Z(k)+Z(k-1))/2.0d0) - Phi(i,j,k))**2
        end do
        end do
        end do
        call MPI_Allreduce(error_poisson, error_poisson_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        if(myrank==0)then
            write(*,*) 'error_poisson = ', sqrt(error_poisson_all/dble(NX*NZ*NY))
        end if   

        !!!2decompの終了
        call decomp_2d_fft_finalize
        call decomp_2d_finalize

        !!!MPIの終了
        call MPI_Finalize(ierr)

end program main
