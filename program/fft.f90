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
    include 'fftw3.f'

    real(8),parameter:: ds = 1.0d0
    integer,parameter:: xmax = 512 !ｘ方向格子数
    integer,parameter:: ymax = 512 !ｙ方向格子数
    integer,parameter:: zmax = 512 !ｚ方向格子数
    integer,parameter:: Nx = 32 !ｘ方向の並列数
    integer,parameter:: Nz = 64 !ｚ方向の並列数
    integer,parameter:: x_procs = (xmax) / Nx
    integer,parameter:: z_procs = (zmax) / Nz
    integer,parameter:: new_procs = Nx * Nz

    integer,parameter:: step_start = 5000000
    integer,parameter:: step_end = 5036000
    integer,parameter:: step_bin = 12000
    integer,parameter:: step_num = (step_end - step_start) / step_bin +1 

    double precision PLAN_2DR2C, PLAN_2DC2R, PLAN_1DC2CF, PLAN_1DC2CB
    integer i, j, k, xi, yi, zi, Nxx, Nzz, step, beta

    real(8) a(xmax, ymax)
    complex(kind(0d0)) a_freq(int(xmax/2)+1, ymax)
    complex(kind(0d0)) b(zmax), b_freq(zmax)

    real(8) u1(xmax, ymax, zmax), u2(xmax, ymax, zmax), u3(xmax, ymax, zmax), p(xmax, ymax, zmax)
    real(8) u1out(x_procs, ymax, z_procs, 0:new_procs-1), u2out(x_procs, ymax, z_procs, 0:new_procs-1), u3out(x_procs, ymax, z_procs, 0:new_procs-1), pout(x_procs, ymax, z_procs, 0:new_procs-1)
    real(8) dummy1, dummy2

    complex(kind(0d0)) u_freq(int(xmax/2)+1, ymax, zmax), u1_freq_xyz(int(xmax/2)+1, ymax, zmax), u2_freq_xyz(int(xmax/2)+1, ymax, zmax), u3_freq_xyz(int(xmax/2)+1, ymax, zmax)

    complex(kind(0d0)) u1a_freq_xyz(int(xmax/2)+1, ymax, zmax), u2a_freq_xyz(int(xmax/2)+1, ymax, zmax), u3a_freq_xyz(int(xmax/2)+1, ymax, zmax)
    complex(kind(0d0)) u1b_freq_xyz(int(xmax/2)+1, ymax, zmax), u2b_freq_xyz(int(xmax/2)+1, ymax, zmax), u3b_freq_xyz(int(xmax/2)+1, ymax, zmax)
    complex(kind(0d0)) u1c_freq_xyz(int(xmax/2)+1, ymax, zmax), u2c_freq_xyz(int(xmax/2)+1, ymax, zmax), u3c_freq_xyz(int(xmax/2)+1, ymax, zmax)
    real(8) u1a(xmax, ymax, zmax), u2a(xmax, ymax, zmax), u3a(xmax, ymax, zmax)
    real(8) u1b(xmax, ymax, zmax), u2b(xmax, ymax, zmax), u3b(xmax, ymax, zmax)
    real(8) u1c(xmax, ymax, zmax), u2c(xmax, ymax, zmax), u3c(xmax, ymax, zmax)
    real(8) u1a_2(0:xmax+1, 0:ymax+1, 0:zmax+1), u2a_2(0:xmax+1, 0:ymax+1, 0:zmax+1), u3a_2(0:xmax+1, 0:ymax+1, 0:zmax+1)

    real(8),allocatable::  energy_spectrum(:)
    real(8) kx, ky, kz
    real(8) dkx, dky, dkz, dk
    real(8) kmax_abs, k_abs
    integer k_index

    real(8) kinetic(xmax, ymax, zmax)
    real(8) ene_sum, kine_sum

    !エンストロフィー変数
    real(8) grad_u(1:3,1:3,1:xmax,1:ymax,1:zmax)
    real(8) vorticity(1:3,1:xmax,1:ymax,1:zmax)
    real(8) enstrophy(1:xmax,1:ymax,1:zmax)
    real(8) ave_enstrophy

    !粒子速度
    integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
    integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
    integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)
    real(8):: cr(1:3, 1:15)

    real(8),parameter:: pi = acos(-1.0d0) !円周率

    !ティレクトリー読み込み
    character(*),parameter :: datadir = "/data/lng/nakanog/taylor_512_re100000_large/"
    ! character(*),parameter :: datadir = "/data/sht/nakanog/taylor_512_dsd/"
    !ディレクトリ作成
    character(*),parameter :: datadir2 = "/data/lng/nakanog/taylor_512_re100000_large/plot/"

    character(8) file_num, file_num2
    character :: filename*200
    integer set
    real(8) time1, time2

    ! call mk_dirs(datadir2)
    open(10,file="./procs.d")
    close(10)
    call cpu_time(time1)

!===============データの読み込み=========================================
DO step = step_start, step_end, step_bin
    ! step = step_start
    u1out(:,:,:,:) = 0.0d0
    u2out(:,:,:,:) = 0.0d0
    u3out(:,:,:,:) = 0.0d0
    pout(:,:,:,:) = 0.0d0
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
        call cpu_time(time2)
        open(10,file="./procs.d",action="write",position="append")
        write(10,*) i, time2-time1
        close(10)

        if((i >= 0) .and. (i < 10)) then
            write(file_num, "(i1)") i
        elseif((i > 9) .and. (i < 100)) then
            write(file_num,"(i2)") i
        elseif((i > 99) .and. (i < 1000)) then
            write(file_num,"(i3)") i
        elseif((i > 999) .and. (i < 10000)) then
            write(file_num,"(i4)") i
        endif

        open(20, file=datadir//"0_"//trim(file_num)//"_"//trim(file_num2)//".bin", form="unformatted")
        do zi = 1, z_procs
            do yi = 1, ymax
                do xi = 1, x_procs
                    ! read(20) dummy1, u1out(xi,yi,zi,i), u2out(xi,yi,zi,i), u3out(xi,yi,zi,i), pout(xi,yi,zi,i)
                    read(20) u1out(xi,yi,zi,i), u2out(xi,yi,zi,i), u3out(xi,yi,zi,i), pout(xi,yi,zi,i)
                enddo
            enddo
        enddo
        close(20)
    enddo
    k = 0
    do Nzz=0,Nz-1
        do Nxx=0,Nx-1
            do zi=1,z_procs
                do yi=1,ymax
                    do xi=1,x_procs
                        u1((xi)+Nxx*x_procs,yi,(zi)+Nzz*z_procs) = dble(u1out(xi,yi,zi,k))
                        u2((xi)+Nxx*x_procs,yi,(zi)+Nzz*z_procs) = dble(u2out(xi,yi,zi,k))
                        u3((xi)+Nxx*x_procs,yi,(zi)+Nzz*z_procs) = dble(u3out(xi,yi,zi,k))
                        p((xi)+Nxx*x_procs,yi,(zi)+Nzz*z_procs) = dble(pout(xi,yi,zi,k))
                    enddo
                enddo
            enddo
            k = k + 1
        enddo
    enddo

    kine_sum = 0.0d0
    do zi=1,zmax
        do yi=1,ymax
            do xi=1,xmax
                kinetic(xi,yi,zi) = 0.5d0 * (u1(xi,yi,zi)*u1(xi,yi,zi) + u2(xi,yi,zi)*u2(xi,yi,zi) + u3(xi,yi,zi)*u3(xi,yi,zi))
                kine_sum = kine_sum + kinetic(xi,yi,zi)
            enddo
        enddo
    enddo
    kine_sum = kine_sum / dble(xmax*ymax*zmax)

!==============================フーリエ変換=========================================
    call dfftw_plan_dft_r2c_2d(PLAN_2DR2C, xmax, ymax, a(:, :), a_freq(:, :), FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_2d(PLAN_2DC2R, xmax, ymax, a_freq(:, :), a(:, :), FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(PLAN_1DC2CF, zmax, b(:), b_freq(:), FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(PLAN_1DC2CB, zmax, b_freq(:), b(:), FFTW_BACKWARD, FFTW_ESTIMATE) 
    !初期化
    u1_freq_xyz(:,:,:) = (0.0d0, 0.0d0)
    u2_freq_xyz(:,:,:) = (0.0d0, 0.0d0)
    u3_freq_xyz(:,:,:) = (0.0d0, 0.0d0)

    u_freq(:,:,:) = (0.0d0, 0.0d0)
    !x, y方向に２次元fft
    do zi=1,zmax
        call dfftw_execute_dft_r2c(PLAN_2DR2C, u1(:, :, zi), u_freq(:, :, zi))
        u_freq(:, :, zi) = u_freq(:, :, zi) / dble(xmax*ymax)
    enddo

    !z方向に１次元fft
    do yi=1,ymax
        do xi=1,int(xmax/2)+1
            call dfftw_execute_dft(PLAN_1DC2CF, u_freq(xi, yi, :), u1_freq_xyz(xi, yi, :))
            u1_freq_xyz(xi, yi, :) = u1_freq_xyz(xi, yi, :) / dble(zmax)
        enddo
    enddo

    u_freq(:,:,:) = (0.0d0, 0.0d0)
    !x, y方向に２次元fft
    do zi=1,zmax
        call dfftw_execute_dft_r2c(PLAN_2DR2C, u2(:, :, zi), u_freq(:, :, zi))
        u_freq(:, :, zi) = u_freq(:, :, zi) / dble(xmax*ymax)
    enddo

    !z方向に１次元fft
    do yi=1,ymax
        do xi=1,int(xmax/2)+1
            call dfftw_execute_dft(PLAN_1DC2CF, u_freq(xi, yi, :), u2_freq_xyz(xi, yi, :))
            u2_freq_xyz(xi, yi, :) = u2_freq_xyz(xi, yi, :) / dble(zmax)
        enddo
    enddo

    u_freq(:,:,:) = (0.0d0, 0.0d0)
    !x, y方向に２次元fft
    do zi=1,zmax
        call dfftw_execute_dft_r2c(PLAN_2DR2C, u3(:, :, zi), u_freq(:, :, zi))
        u_freq(:, :, zi) = u_freq(:, :, zi) / dble(xmax*ymax)
    enddo

    !z方向に１次元fft
    do yi=1,ymax
        do xi=1,int(xmax/2)+1
            call dfftw_execute_dft(PLAN_1DC2CF, u_freq(xi, yi, :), u3_freq_xyz(xi, yi, :))
            u3_freq_xyz(xi, yi, :) = u3_freq_xyz(xi, yi, :) / dble(zmax)
        enddo
    enddo

    call dfftw_destroy_plan(PLAN_2DR2C)
    call dfftw_destroy_plan(PLAN_2DC2R)
    call dfftw_destroy_plan(PLAN_1DC2CF)
    call dfftw_destroy_plan(PLAN_1DC2CB)

!==========================エネルギースペクトル============================================
    if(step == step_start) then
        dkx = 2.0d0*pi/dble(xmax)
        dky = 2.0d0*pi/dble(ymax)
        dkz = 2.0d0*pi/dble(zmax)
        dk = ( (dkx)**2.0d0 + (dky)**2.0d0 + (dkz)**2.0d0 )**0.5d0
        kmax_abs = ( (dkx*dble(xmax)/2.0d0)**2 + (dky*dble(ymax)/2.0d0)**2 + (dkz*dble(zmax)/2.0d0)**2 )**0.5d0
        allocate(energy_spectrum( floor(kmax_abs / dk) + 1))
        energy_spectrum(:) = 0.0d0
    endif

    do zi=1,zmax
        do yi=1,ymax
            do xi=1,int(xmax/2)+1
                kx = dble(xi) - 1.0d0
                ky = sign(1.0d0, dble(ymax+3)/2.0d0 - dble(yi)) * (-abs(dble(ymax)/2.0d0 + 1.0d0 - dble(yi)) + dble(ymax)/2.0d0)
                kz = sign(1.0d0, dble(zmax+3)/2.0d0 - dble(zi)) * (-abs(dble(zmax)/2.0d0 + 1.0d0 - dble(zi)) + dble(zmax)/2.0d0)
                k_abs = ( kx**2.0d0 + ky**2.0d0 + kz**2.0d0 )**0.5d0
                k_index = int( k_abs ) + 1

                energy_spectrum(k_index) = energy_spectrum(k_index) + 0.5d0 * (abs(u1_freq_xyz(xi,yi,zi))**2.0d0 + abs(u2_freq_xyz(xi,yi,zi))**2.0d0 + abs(u3_freq_xyz(xi,yi,zi))**2.0d0) * min(2.0d0, dble(xi)) / dkx
            enddo
        enddo
    enddo

    open(10,file="./procs.d",action="write",position="append")
    write(10,*) step, "is finished!!!!!"
    close(10)

    ! open(11,file="energy_spectrum_DNS_1.d")
    ! do i = 1, size(energy_spectrum)
    !     write(11, *) (dble(i-1) + 0.5d0)*dk / (0.5d0*dk), energy_spectrum(i), energy_spectrum(i) * (((dble(i-1) + 0.5d0)*dk)**(5.0d0/3.0d0))
    ! enddo
    ! close(11)

    ! ene_sum = 0.0d0
    ! do i= 1, size(energy_spectrum)
    !     ene_sum = ene_sum + energy_spectrum(i)*dk
    ! enddo

    ene_sum = 0.0d0
    do zi=1,zmax
        do yi=1,ymax
            do xi=1,int(xmax/2)+1
                ene_sum = ene_sum + (abs(u1_freq_xyz(xi,yi,zi))**2 + abs(u2_freq_xyz(xi,yi,zi))**2 + abs(u3_freq_xyz(xi,yi,zi))**2)/2.0d0 * min(2.0d0, dble(xi))
            enddo
        enddo
    enddo

    open(12,file="kakunin.d")
    write(12,*) ene_sum, kine_sum, abs(ene_sum - kine_sum)/abs(kine_sum)
    close(12)
ENDDO



    open(11,file="enesupe_single_2.d")
    do i = 1, size(energy_spectrum)
        write(11, *) dble(i)-1.0d0, energy_spectrum(i) / real(step_num)
        ! write(11, *) dble(i)-1.0d0, energy_spectrum(i)
    enddo
    close(11)
end program main

! !================================逆フーリエ変換==========================================
!     u1_freq(:, :, :) = 0.0d0
!     reu1(:, :, :) = 0.0d0
!     !z方向に１次元逆fft
!     do yi=1,ymax
!         do xi=1,int(xmax/2)+1
!             call dfftw_execute_dft(PLAN_1DC2CB, u1_freq_xyz(xi, yi, :), u1_freq(xi, yi, :))
!         enddo
!     enddo

!     !x, y方向に２次元逆fft
!     do zi=1,zmax
!         call dfftw_execute_dft_c2r(PLAN_2DC2R, u1_freq(:, :, zi), reu1(:, :, zi))
!     enddo
