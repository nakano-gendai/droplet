module globals
    integer, parameter :: n = 129 !データの行数
    integer, parameter :: m = 129
    integer,parameter:: xmax = 128 !ｘ方向格子数
    integer,parameter:: ymax = 64
    integer,parameter:: zmax = 128 !ｚ方向格子数
    real(8),parameter:: phi1 = 2.211d0
    real(8),parameter:: phi2 = 4.895d0
    real(8) phi_kaimen 
    integer,parameter:: num = 500000
    real(8),parameter:: uw = 0.1d0
    character(8) file_num
    character :: filename*200
    integer xi, yi, zi, i, j, k, step
    ! character(*),parameter :: dir = "/home/nakano/anime/eta=1.0/Re=20.0/data_eta1re20ca0.075/"
    character(*),parameter :: dir = "C:\Users\genda\Documents\data_eta=0.1_ca=0.4\"
    character(*),parameter :: datadir = "C:\Users\genda\Documents\tauxz_re0.2eta0.1ca0.4\"

contains
    subroutine senkei(xx1,xx2,zz1,zz2,xxg,zzg,L1,B1,phi)
        integer, intent(in) :: xx1,xx2,zz1,zz2
        real(8), intent(in) :: xxg,zzg
        real(8), intent(inout) :: L1, B1
        real(8), intent(in) :: phi(0:xmax, 0:ymax, 0:zmax)
        real(8) distance
        real(8) alpha, l
        real(8) xj, zj
        
        

        l = sqrt((dble(xx2)-dble(xx1))**2 +(dble(zz2)-dble(zz1))**2)
        alpha = l * (phi_kaimen - phi(xx1,ymax/2,zz1)) / (phi(xx2,ymax/2,zz2) - phi(xx1,ymax/2,zz1))
        if((alpha >= 0.0d0) .and. (alpha <= l)) then
            xj = dble(xx1) + alpha/l * (dble(xx2) - dble(xx1))
            zj = dble(zz1) + alpha/l * (dble(zz2) - dble(zz1))
            distance = sqrt((xj-xxg)**2 + (zj-zzg)**2)
            if(L1 <= distance) then
                L1 = distance
            elseif(B1 >= distance) then
                B1 = distance
            else
                L1 = L1
                B1 = B1
            endif
        endif
    end subroutine senkei

    subroutine output(tauxz)
        real(8),intent(inout):: tauxz(0:xmax,0:ymax,0:zmax)

        write(filename,*) step !i->filename 変換
        filename=//trim(adjustl(filename))//'.d' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
        print *, filename !表示してみる
        open(20,file=filename, status='replace') 

        yi = ymax /2 
        do zi=0,zmax
            ! do yi=0,ymax
                do xi=0,xmax
                    ! write(20,"(8es16.8)") dble(xi),dble(zi),u1(xi,yi,zi),u2(xi,yi,zi),u3(xi,yi,zi),sqrt(u1(xi,yi,zi)**2+u2(xi,yi,zi)**2+u3(xi,yi,zi)**2),p(xi,yi,zi),phi(xi,yi,zi)
                    write(20,"(3es16.8)") dble(xi),dble(zi),tauxz(xi,yi,zi)
                enddo
                write(20,*)
            ! enddo
            ! write(20,*)
        enddo
        close(20)
    end subroutine output

end module globals

program main
use globals
    implicit none
    real(8) x, z
    real(8) u1(-1:xmax+1, -1:ymax+1, -1:zmax+1)
    real(8) u2(-1:xmax+1, -1:ymax+1, -1:zmax+1)
    real(8) u3(-1:xmax+1, -1:ymax+1, -1:zmax+1)
    real(8) grad_u1(1:3, 0:xmax, 0:ymax, 0:zmax)
    real(8) grad_u2(1:3, 0:xmax, 0:ymax, 0:zmax)
    real(8) grad_u3(1:3, 0:xmax, 0:ymax, 0:zmax)
    real(8) tauxz(0:xmax, 0:ymax, 0:zmax)
    real(8) p(0:xmax, 0:ymax, 0:zmax)
    real(8) phi(0:xmax, 0:ymax, 0:zmax)
    real(8) u_m(0:xmax, 0:ymax, 0:zmax)
    
    real(8) phi_max, phi_min

    !重心計算用変数
    real(8) xg, yg, zg
    real(8) phi_tmp
    !線形補間用変数
    integer x0, x1, x2, z0, z1, z2
    real(8),parameter:: pi = acos(-1.0d0)
    

    !液滴の変形
    real(8) L1, B1, Dd


DO step = 100, num, 100

    !ファイル読み込み
    u1(:,:,:) = 0.0d0
    u2(:,:,:) = 0.0d0
    u3(:,:,:) = 0.0d0
    u_m(:,:,:) = 0.0d0
    p(:,:,:) = 0.0d0
    phi(:,:,:) = 0.0d0
    if((step > 99) .and. (step < 1000)) then
        write(file_num, "(i3)") step
    elseif((step > 999) .and. (step < 10000)) then
        write(file_num,"(i4)") step
    elseif((step > 9999) .and. (step < 100000)) then
        write(file_num,"(i5)") step
    elseif((step > 99999) .and. (step < 1000000)) then
        write(file_num,"(i6)") step
    endif
    open(20, file=dir//trim(file_num)//".d")
    yi = ymax/2
    do zi = 0, m-1
        do xi = 0, n-1
            ! read(20,"(8es16.8)") x, z, u1(xi,yi,zi), u2(xi,yi,zi), u3(xi,yi,zi), u_m(xi,yi,zi),p(xi,yi,zi), phi(xi,yi,zi)
            read(20,"(7es16.8)") x, z, u1(xi,yi,zi), u2(xi,yi,zi), u3(xi,yi,zi), p(xi,yi,zi), phi(xi,yi,zi)
        enddo
        read(20,"()")
    enddo
    close(20)

    !周期境界
    do zi=0,zmax
        do yi=-1,ymax+1
            do xi=0,xmax
                u1(xi,-1,zi) = u1(xi,ymax,zi)
                u1(xi,ymax+1,zi) = u1(xi,0,zi)
                u1(-1,yi,zi) = u1(xmax,yi,zi)
                u1(xmax+1,yi,zi) = u1(0,yi,zi)

                u2(xi,-1,zi) = u2(xi,ymax,zi)
                u2(xi,ymax+1,zi) = u2(xi,0,zi)
                u2(-1,yi,zi) = u2(xmax,yi,zi)
                u2(xmax+1,yi,zi) = u2(0,yi,zi)

                u3(xi,-1,zi) = u3(xi,ymax,zi)
                u3(xi,ymax+1,zi) = u3(xi,0,zi)
                u3(-1,yi,zi) = u3(xmax,yi,zi)
                u3(xmax+1,yi,zi) = u3(0,yi,zi)
            enddo
        enddo
    enddo
    do yi=-1,ymax+1
        do xi=-1,xmax+1
            u1(xi,yi,-1) = u1(xi,yi,zmax)
            u1(xi,yi,zmax+1) = u1(xi,yi,0)

            u2(xi,yi,-1) = u2(xi,yi,zmax)
            u2(xi,yi,zmax+1) = u2(xi,yi,0)

            u3(xi,yi,-1) = u3(xi,yi,zmax)
            u3(xi,yi,zmax+1) = u3(xi,yi,0)
        enddo
    enddo
    yi = ymax/2
    !grad計算
    do zi=0,zmax
        do xi=0,xmax
            grad_u1(3,xi,yi,zi) = (u1(xi,yi,zi+1) - u1(xi,yi,zi-1))/10.0d0
            grad_u3(1,xi,yi,zi) = (u3(xi+1,yi,zi) - u3(xi-1,yi,zi))/10.0d0
        enddo
    enddo
    call output(tauxz)
    write(*,*) "step =", step
ENDDO
end program main

