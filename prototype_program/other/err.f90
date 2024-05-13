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
        real(8),parameter:: ds = 1.0d0
        integer,parameter:: xmax = 255 !ｘ方向格子数
        integer,parameter:: ymax = 255 !ｙ方向格子数
        integer,parameter:: zmax = 255 !ｚ方向格子数
        integer,parameter:: Nx = 8 !ｘ方向の並列数
        integer,parameter:: Nz = 16 !ｚ方向の並列数
        integer,parameter:: x_procs = (xmax+1) / Nx
        integer,parameter:: z_procs = (zmax+1) / Nz
        integer,parameter:: new_procs = Nx * Nz
    
        integer i, j, k, xi, yi, zi, Nxx, Nzz, step, beta
    
        real(8) u1(0:xmax, 0:ymax, 0:zmax), u2(0:xmax, 0:ymax, 0:zmax), u3(0:xmax, 0:ymax, 0:zmax), p(0:xmax, 0:ymax, 0:zmax)
        real(8) u1out(x_procs, ymax+1, z_procs, 0:new_procs-1), u2out(x_procs, ymax+1, z_procs, 0:new_procs-1), u3out(x_procs, ymax+1, z_procs, 0:new_procs-1), pout(x_procs, ymax+1, z_procs, 0:new_procs-1)
        real(8) dummy1, dummy2
        real(8) u1_th(0:xmax, 0:ymax, 0:zmax), u2_th(0:xmax, 0:ymax, 0:zmax), u3_th(0:xmax, 0:ymax, 0:zmax), p_th(0:xmax, 0:ymax, 0:zmax)
        real(8) u1_err(0:xmax, 0:ymax, 0:zmax), u2_err(0:xmax, 0:ymax, 0:zmax), u3_err(0:xmax, 0:ymax, 0:zmax), p_err(0:xmax, 0:ymax, 0:zmax)
        real(8) u_err(0:xmax, 0:ymax, 0:zmax)
    
        !粒子速度
        integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
        integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
        integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)
        real(8):: cr(1:3, 1:15)
    
        real(8),parameter:: pi = acos(-1.0d0) !円周率
    
        !ティレクトリー読み込み
        character(*),parameter :: datadir = "/data/sht/nakanog/taylor_laminar_re10/"
    
        character(8) file_num, file_num2
        character :: filename*200
        integer set
        real(8) time1, time2
    
        real(8) errx, erru, errz, errp
        real(8) u_abs, u_abs_th
        real(8),parameter:: D_vortex = 127.5d0 !渦の大きさ
        real(8),parameter:: kw = pi/D_vortex !波数
        real(8),parameter:: umax = 0.001d0 !最大流速
        real(8),parameter:: Re = 10.0d0 !レイノルズ数
    
        open(10,file="./procs.d",status="replace")
        close(10)
        call cpu_time(time1)
        open(11,file="./err_sum.d",status="replace")
        close(11)
    
    !===============データの読み込み=========================================
    ! DO step = 4000, 660000, 4000
        step = 660000
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
    
            set = 20 + i
            open(set, file=datadir//"0_"//trim(file_num)//"_"//trim(file_num2)//".bin", form="unformatted")
            do zi = 1, z_procs
                do yi = 1, ymax+1
                    do xi = 1, x_procs
                        read(set) u1out(xi,yi,zi,i), u2out(xi,yi,zi,i), u3out(xi,yi,zi,i), pout(xi,yi,zi,i), dummy1
                    enddo
                enddo
            enddo
            close(set)
        enddo
        k = 0
        do Nzz=0,Nz-1
            do Nxx=0,Nx-1
                do zi=1,z_procs
                    do yi=1,ymax+1
                        do xi=1,x_procs
                            u1((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*z_procs) = u1out(xi,yi,zi,k)
                            u2((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*z_procs) = u2out(xi,yi,zi,k)
                            u3((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*z_procs) = u3out(xi,yi,zi,k)
                            p((xi-1)+Nxx*x_procs,yi-1,(zi-1)+Nzz*z_procs) = pout(xi,yi,zi,k)
                        enddo
                    enddo
                enddo
                k = k + 1
            enddo
        enddo

    !     do zi = 0, zmax
    !         do yi = 0, ymax
    !             do xi = 0, xmax
    !                 u1_th(xi,yi,zi) = umax*sin(kw*dble(xi))*cos(kw*dble(zi))
    !                 u2_th(xi,yi,zi) = 0.0d0
    !                 u3_th(xi,yi,zi) = -umax*cos(kw*dble(xi))*sin(kw*dble(zi))
    !                 p_th(xi,yi,zi) = 1.0d0/3.0d0 + 1.0d0/4.0d0*umax*umax*(cos(2.0d0*kw*dble(xi))+cos(2.0d0*kw*dble(zi)))

    !                 if(u1_th(xi,yi,zi) == 0.0d0) then
    !                     u1_th(xi,yi,zi) = 1.0d-16
    !                 endif
    !                 if(u3_th(xi,yi,zi) == 0.0d0) then
    !                     u3_th(xi,yi,zi) = 1.0d-16
    !                 endif
    !                 if(p_th(xi,yi,zi) == 0.0d0) then
    !                     p_th(xi,yi,zi) = 1.0d-16
    !                 endif
    !             enddo
    !         enddo
    !     enddo

    
    !     errx = 0.0d0
    !     errz = 0.0d0
    !     errp = 0.0d0
    !     erru = 0.0d0
    !     do zi = 0, zmax
    !         do yi = 0, ymax
    !             do xi = 0, xmax
    !                 errx = errx + abs( u1(xi,yi,zi) - u1_th(xi,yi,zi) )
    !                 errz = errz + abs( u3(xi,yi,zi) - u3_th(xi,yi,zi) ) 
    !                 errp = errp + abs( p(xi,yi,zi) - p_th(xi,yi,zi) ) 

    !                 u_abs_th = (u1_th(xi,yi,zi)*u1_th(xi,yi,zi) + u3_th(xi,yi,zi)*u3_th(xi,yi,zi))**0.5d0
    !                 u_abs = (u1(xi,yi,zi)*u1(xi,yi,zi) + u2(xi,yi,zi)*u2(xi,yi,zi)+ u3(xi,yi,zi)*u3(xi,yi,zi))**0.5d0
    !                 erru = erru + abs( u_abs - u_abs_th )
    !             enddo
    !         enddo
    !     enddo
    !     errx = errx / (dble(xmax+1)*dble(ymax+1)*dble(zmax+1))
    !     errz = errz / (dble(xmax+1)*dble(ymax+1)*dble(zmax+1))
    !     errp = errp / (dble(xmax+1)*dble(ymax+1)*dble(zmax+1))
    !     erru = erru / (dble(xmax+1)*dble(ymax+1)*dble(zmax+1))
    !     open(11,file="./err_sum.d",action="write",position="append")
    !     write(11,*) step, errx, errz, errp, erru
    !     close(11)
    ! ENDDO
    ! do zi = 0, zmax
    !     do yi = 0, ymax
    !         do xi = 0, xmax
    !             u1_err(xi,yi,zi) = abs(u1(xi,yi,zi) - u1_th(xi,yi,zi))
    !             u3_err(xi,yi,zi) = abs(u3(xi,yi,zi) - u3_th(xi,yi,zi))
    !             p_err(xi,yi,zi) = abs(p(xi,yi,zi) - p_th(xi,yi,zi))

    !             u_abs_th = (u1_th(xi,yi,zi)*u1_th(xi,yi,zi) + u3_th(xi,yi,zi)*u3_th(xi,yi,zi))**0.5d0
    !             u_abs = (u1(xi,yi,zi)*u1(xi,yi,zi) + u2(xi,yi,zi)*u2(xi,yi,zi)+ u3(xi,yi,zi)*u3(xi,yi,zi))**0.5d0
    !             u_err(xi,yi,zi) = abs(u_abs - u_abs_th)
    !         enddo
    !     enddo
    ! enddo
    open(12,file="./u.d",status="replace")
    do zi = 0, zmax
        yi = (ymax+1) / 2
        ! do yi = 0, ymax
            do xi = 0, xmax
                write(12,*) xi, zi, 1000.0d0*u1(xi,yi,zi), 1000.0d0*u3(xi,yi,zi)
            enddo
        ! enddo
    enddo
    close(12)
    end program main