module globals
    ! integer, parameter :: n = 300 !データの行数
    ! integer, parameter :: l = 105
    ! integer, parameter :: m = 105
    integer,parameter:: xmax = 9 !ｘ方向格子数
    integer,parameter:: ymax = 9
    integer,parameter:: zmax = 9 !ｚ方向格子数
    ! real(8),parameter:: phi1 = 2.211d0
    ! real(8),parameter:: phi2 = 4.895d0
    ! real(8),parameter:: phi_kaimen = 0.5d0 * (phi1 + phi2)
    ! integer,parameter:: num = 200000 !ファイルの数
    ! real(8),parameter:: uw = 0.02d0
    ! real(8),parameter:: kappag = 7.69d-3
    ! real(8),parameter:: ds = 1.0d0
    integer xi, yi, zi
    integer step

    character(8) file_num
    character :: filename*200
    !ディレクトリ読み込み
    character(*),parameter :: dir = "./"
    ! character(*),parameter :: dir = "/data/2022/nakano/rechange_caconst/data_re100ca0.25/"
    !ディレクトリ作成
    character(*),parameter :: datadir = "./"

contains

    subroutine input(phi,step)
        real(8),intent(inout) :: phi(0:xmax, 0:ymax, 0:zmax)
        integer,intent(in) :: step
        !ファイル読み込み
        phi(:,:,:) = 0.0d0
        ! if((step > 99) .and. (step < 1000)) then
        !     write(file_num, "(i3)") step
        ! elseif((step > 999) .and. (step < 10000)) then
        !     write(file_num,"(i4)") step
        ! elseif((step > 9999) .and. (step < 100000)) then
        !     write(file_num,"(i5)") step
        ! elseif((step > 99999) .and. (step < 1000000)) then
        !     write(file_num,"(i6)") step
        ! endif
        write(file_num,"(i1)") step
        
        open(20, file=dir//"0_"//trim(file_num)//".d")
        do zi = 0, zmax
            do yi = 0, ymax
                do xi = 0, xmax
                    read(20,"(4es24.16)") x, y, z, phi(xi,yi,zi)
                enddo
            enddo
        enddo
        close(20)
        y = x*z*y
    end subroutine input

    subroutine output(phi, step)
        real(8),intent(inout):: phi(0:xmax,0:ymax,0:zmax)
        integer,intent(in) :: step

        write(filename,*) step !i->filename 変換
        filename=datadir//"gnu_"//trim(adjustl(filename))//'.d' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
        print *, filename !表示してみる
        open(30,file=filename, status='replace') 

        do zi=0,zmax
            do xi=0,xmax
                write(30,"(3es24.16)") dble(xi),dble(zi),phi(xi,(ymax+1)/2,zi)
            enddo
            write(30,*)
        enddo
        close(30)
    end subroutine output

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
    real(8) x, y, z
    real(8) phi(0:xmax, 0:ymax, 0:zmax)
    ! real(8) phi_ini(0:xmax, 0:ymax, 0:zmax)
    ! real(8) sum
    step = 0
    x =0.0d0
    y=0.0d0
    z=0.0d0
    y = x*z*y

    ! open(20, file=dir//"2_0.d")
    ! do zi = 0, zmax
    !     do yi = 0, ymax
    !         do xi = 0, xmax
    !             read(20,"(4es24.16)") x, y, z, phi(xi,yi,zi)
    !         enddo
    !     enddo
    ! enddo
    ! close(20)

    ! open(21, file="ini_02.d")
    ! do zi = 0, zmax
    !     do yi = 0, ymax
    !         do xi = 0, xmax
    !             read(21,"(4es24.16)") x, y, z, phi_ini(xi,yi,zi)
    !         enddo
    !     enddo
    ! enddo
    ! close(21)
    

    ! sum = 0.0d0
    ! open(22,file="dif.d")
    ! do zi=0,zmax
    !     do yi=0,ymax
    !         do xi=0,xmax
    !             write(22,*) phi(xi,yi,zi)-phi_ini(xi,yi,zi)
    !             sum = sum + (phi(xi,yi,zi)-phi_ini(xi,yi,zi))
    !         enddo
    !     enddo
    ! enddo
    ! write(22,*) "sum = ", sum
    ! close(22)


    call mk_dirs(datadir)
! do step = 100, num, 100
    call input(phi,step)
    call output(phi, step)
    write(*,*) "step = ", step
    
! enddo
end program main

