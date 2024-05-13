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
    integer,parameter:: xmax = 160 !ｘ方向格子数
    integer,parameter:: ymax = 160 !ｙ方向格子数
    integer,parameter:: zmax = 160 !ｚ方向格子数
    integer,parameter:: Nx = 7 !ｘ方向の並列数
    integer,parameter:: Nz = 7 !ｚ方向の並列数
    integer,parameter:: x_procs = (xmax+1) / Nx
    integer,parameter:: z_procs = (zmax+1) / Nz
    integer,parameter:: new_procs = Nx * Nz
    integer i, k, xi, yi, zi, Nxx, Nzz, step
    real x, y, z, s
    character(8) file_num, file_num2
    character :: filename*200
    integer set
    !ティレクトリー読み込み
    character(*),parameter :: datadir = "/data/sht/nakanog/taylor_re1000we3.5/"
    !ディレクトリ作成
    character(*),parameter :: datadir2 = "/data/sht/nakanog/taylor_re1000we3.5/gnuvec/"
    real phi(0:xmax,0:zmax)
    real phiout(1:x_procs,1:z_procs,0:new_procs-1)
    real u1(0:xmax,0:zmax)
    real u1out(1:x_procs,1:z_procs,0:new_procs-1)
    real u3(0:xmax,0:zmax)
    real u3out(1:x_procs,1:z_procs,0:new_procs-1)

    call mk_dirs(datadir2)
    ! do step=1000, 50000, 1000
    step = 10000
        if((step > 99) .and. (step < 1000)) then
            write(file_num2, "(i3)") step
        elseif((step > 999) .and. (step < 10000)) then
            write(file_num2,"(i4)") step
        elseif((step > 9999) .and. (step < 100000)) then
            write(file_num2,"(i5)") step
        elseif((step > 99999) .and. (step < 1000000)) then
            write(file_num2,"(i6)") step
        endif
        do i=0,new_procs-1
            if((i >= 0) .and. (i < 10)) then
                write(file_num, "(i1)") i
            elseif((i > 9) .and. (i < 100)) then
                write(file_num,"(i2)") i
            elseif((i > 99) .and. (i < 1000)) then
                write(file_num,"(i3)") i
            elseif((i > 999) .and. (i < 10000)) then
                write(file_num,"(i4)") i
            endif

            set = 20 + i
            open(set, file=datadir//"0_"//trim(file_num)//"_"//trim(file_num2)//".bin", form="unformatted")
            do zi = 1, z_procs
                do xi = 1, x_procs
                    read(set)  phiout(xi,zi,i), u1out(xi,zi,i), y, u3out(xi,zi,i), s
                enddo
                read(set)
            enddo
            close(set)
        enddo
        k = 0
        do Nzz=0,Nz-1
            do Nxx=0,Nx-1
                do zi=1,z_procs
                    do xi=1,x_procs
                        phi((xi-1)+Nxx*x_procs,(zi-1)+Nzz*z_procs) = phiout(xi,zi,k)
                        u1((xi-1)+Nxx*x_procs,(zi-1)+Nzz*z_procs) = u1out(xi,zi,k)
                        u3((xi-1)+Nxx*x_procs,(zi-1)+Nzz*z_procs) = u3out(xi,zi,k)
                    enddo
                enddo
                k = k + 1
            enddo
        enddo
            
        write(filename,*) step !i->filename 変換
        filename=datadir2//"0_"//trim(adjustl(filename))//'_v.d' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
        print *, filename !表示してみる
        ! open(11, file=filename, form="unformatted", status='replace') 
        open(11, file=filename, status='replace') 
        do zi=0,zmax,5
            do xi=0,xmax,5
                write(11,*) dble(xi), dble(zi), 1000.0d0*u1(xi,zi), 1000.0d0*u3(xi,zi)
            enddo
            ! write(11,*)
        enddo
        close(11)
    ! enddo

    ! open(12,file="couette_r65.d")
    ! do zi = 0,zmax
    !     write(12,*) real(zi), phi(64,zi)
    ! enddo
    ! close(12)

    ! step = 30000
    ! write(filename,*) step !i->filename 変換
    ! filename=datadir2//trim(adjustl(filename))//'.bin' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
    ! print *, filename !表示してみる
    ! open(12,file=filename, form="unformatted") 
    ! do zi = 0, zmax
    !     do yi = 0, ymax
    !         do xi = 0, xmax
    !             read(12) x,y,z,phi2(xi,yi,zi),u12(xi,yi,zi),u22(xi,yi,zi),u32(xi,yi,zi),p2(xi,yi,zi)
    !         enddo
    !     enddo
    ! enddo
    ! close(12)
    ! wa1 = 0.0d0
    ! dif1 = 0.0d0
    ! wa2 = 0.0d0
    ! dif2 = 0.0d0
    ! wa3 = 0.0d0
    ! dif3 = 0.0d0
    ! wa4 = 0.0d0
    ! dif4 = 0.0d0
    ! wa5 = 0.0d0
    ! dif5 = 0.0d0
    ! do zi=0,zmax
    !     do yi=0,ymax
    !         do xi=0,xmax
    !             dif1 = phi(xi,yi,zi) - phi2(xi,yi,zi)
    !             wa1 = wa1 + dif1
    !             dif2 = u1(xi,yi,zi) - u12(xi,yi,zi)
    !             wa2 = wa2 + dif2
    !             dif3 = u2(xi,yi,zi) - u22(xi,yi,zi)
    !             wa3 = wa3 + dif3
    !             dif4 = u3(xi,yi,zi) - u32(xi,yi,zi)
    !             wa4 = wa4 + dif4
    !             dif5 = p(xi,yi,zi) - p2(xi,yi,zi)
    !             wa5 = wa5 + dif5
    !         enddo
    !     enddo
    ! enddo
    ! open(14,file="ca0.4_30000.d")
    ! do zi=0,zmax
    !     do yi=0,ymax
    !         do xi=0,xmax
    !             write(14,"(8es24.16)") p(xi,yi,zi), p2(xi,yi,zi), wa5
    !         enddo
    !     enddo
    ! enddo
    ! close(14)
    ! open(13,file="wa0.4_60000.d")
    ! write(13,*) "phi  u1  u2  u3  p"
    ! write(13,*) wa1, wa2, wa3, wa4, wa5
    ! close(13)

end program main
