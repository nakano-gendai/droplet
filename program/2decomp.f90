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
    double precision, parameter :: pi = 4*atan(1.0d0)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!-----------------------領域・格子点数の設定-------------------------!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    integer, parameter          :: NX = 64, NY = 64, NZ = 64   
    double precision, parameter :: Xmax = 2.0d0*pi, Ymax =  2.0d0*pi, Zmax =  2.0d0*pi
    double precision, parameter :: dX = Xmax / dble(NX), dY = Ymax / dble(NY), dZ = Zmax / dble(NZ) 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!------------------------------MPI/2decomp---------------------------!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, parameter                  :: MPI_Y_NUM = 2, MPI_Z_NUM = 2                          !!!Y方向・Z方向の並列数
    integer, parameter                  :: NY_procs = NY / MPI_Y_NUM                               !!!各ランクにおけるY方向の格子点数
    integer, parameter                  :: NZ_procs = NZ / MPI_Z_NUM                               !!!各ランクにおけるZ方向の格子点数
    integer, save                       :: NY_procs_min, NZ_procs_min, NY_procs_max, NZ_procs_max  
    integer, save                       :: next_ranky, former_ranky, next_rankz, former_rankz 
    integer, save                       :: ierr, procs, myrank, req_send, req_recv, myranky, myrankz
    integer, dimension(MPI_STATUS_SIZE) :: sta_send, sta_recv                        
    !!!↑MPI用, ↓2decomp用
    integer, save                       :: sta(1:3), end(1:3), size(1:3)                            !!!2decomp用変数

    double precision, save  :: X(0:NX-1), Y(0: NY_procs-1), Z(0: NZ_procs-1)

    character(2) chmyrank

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

    end subroutine MPI_settings

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!----------------------------格子・外力点の初期設定----------------------!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine lattice_init_settings
        integer i, j, k

        !!!各ランクの計算領域を設定
        do i = 0, NX - 1
            X(i) = dble(i) * dX
        enddo
        do j = 0, NY_procs - 1
            Y(j) = myranky * NY_procs * dY + dble(j) * dY
        enddo
        do k = 0, NZ_procs - 1 
            Z(k) = myrankz * NZ_procs * dZ + dble(k) * dZ
        enddo

    end subroutine lattice_init_settings

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!----------------------------実数をフーリエ変換するサブルーチン----------------------------!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine fft_r2c(q, q_hat)
        real(8), intent(inout)            :: q(0:NX-1, 0:NY_procs-1, 0:NZ_procs-1)
        complex(kind(0d0)), intent(inout) :: q_hat(sta(1) : end(1), sta(2) : end(2), sta(3) : end(3))
        call decomp_2d_fft_3d(q, q_hat)
        q_hat = q_hat / (NX * NY * NZ)
    end subroutine fft_r2c

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!----------------------------------複素数を逆フーリエ変換----------------------------------!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine fft_c2r(q_hat, q)
        real(8), intent(inout)            :: q(0 : NX-1, 0 : NY_procs-1, 0 : NZ_procs-1)
        complex(kind(0d0)), intent(inout) :: q_hat(sta(1) : end(1), sta(2) : end(2), sta(3) : end(3))
        call decomp_2d_fft_3d(q_hat, q)
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
        double precision    :: q(0:NX-1, 0:NY_procs-1, 0:NZ_procs-1), q_all(0:NX-1, 0:NY-1, 0:NZ-1)
        complex(kind(0d0)), allocatable :: q_hat(:,:,:)
        complex(8) q_hat_all(0:NX/2, 0:NY-1, 0:NZ-1)

        integer j1, j2, j3, k1, k2, k3, k(3)

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

        allocate(q_hat(sta(1)-1 : end(1)-1, sta(2)-1 : end(2)-1, sta(3)-1 : end(3)-1))

        do j3 = 0, NZ_procs - 1
            do j2 = 0, NY_procs - 1
                do j1 = 0, NX - 1
                    q(j1, j2, j3) = sin(X(j1))*sin(Y(j2))*sin(Z(j3))
                enddo
            enddo
        enddo

        call fft_r2c(q, q_hat)
        call fft_c2r(q_hat, q)
        call fft_r2c(q, q_hat)

        write(chmyrank, '(i2.2)') myrank   !myrankを文字型変数に格納
        open(60, file = './test_'//chmyrank//'.d')    !ランク番号が付いたファイルを作成
        do k3 = sta(3)-1, end(3)-1
            k(3) = k3 - judge(k3, NZ)
            do k2 = sta(2)-1, end(2)-1
                k(2) = k2 - judge(k2, NY)
                do k1 = sta(1)-1, end(1)-1
                    k(1) = k1
                    ! if (abs(q_hat(k1, k2, k3)) > 1.0d-10) then
                        write(60, *) k(1), k(2), k(3), q_hat(k1, k2, k3)
                    ! endif
                enddo
            enddo
        enddo
        close(60)

        !!!2decompの終了
        call decomp_2d_fft_finalize
        call decomp_2d_finalize

        !!!MPIの終了
        call MPI_Finalize(ierr)

end program main