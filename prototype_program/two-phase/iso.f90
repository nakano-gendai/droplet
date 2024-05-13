module globals
    real(8),parameter:: ds = 1.0d0 !格子間隔（lattice unit）
    integer,parameter:: xmax = 30 !ｘ方向格子数
    integer,parameter:: ymax = 30 !ｙ方向格子数
    integer,parameter:: zmax = 30 !ｚ方向格子数
    integer,parameter:: step = 100000
    integer i, n, xi, yi, zi, alpha, beta
    character :: filename*20
    character(*),parameter :: datadir = "./data_p1/"

    !粒子速度（整数）
    integer,parameter:: cx(15) = (/0, 1, 0,  0, -1,  0,  0,  1, -1,  1,  1, -1,  1, -1, -1/)
    integer,parameter:: cy(15) = (/0, 0, 1,  0,  0, -1,  0,  1,  1, -1,  1, -1, -1,  1, -1/)
    integer,parameter:: cz(15) = (/0, 0, 0,  1,  0,  0, -1,  1,  1,  1, -1, -1, -1, -1,  1/)
    real(8):: cr(1:3, 1:15)  !粒子速度（実数）
    real(8) krone(1:3,1:3) !クロネッカーのデルタ

    
    !パラメータ
    real(8),parameter:: nu1 = 0.1d0 !液体の動粘度
    real(8),parameter:: nu2 = 0.1d0 !液滴の動粘度
    real(8),parameter:: pi = acos(-1.0d0)
    real(8),parameter:: E(15) = (/ 2.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, &
                                1.0d0/9.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0, &
                                1.0d0/72.0d0, 1.0d0/72.0d0, 1.0d0/72.0d0 /)
    real(8),parameter:: H(15) = (/ 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
                                0.0d0, 0.0d0, 0.0d0 /)
    real(8),parameter:: F(15) = (/ -7.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, &
                                1.d0/24.d0, 1.d0/24.d0, 1.d0/24.d0, 1.d0/24.d0, &
                                1.d0/24.d0, 1.d0/24.d0, 1.d0/24.d0, 1.d0/24.d0 /)
    real(8),parameter:: D = 14.0d0 !設置する液滴のサイズ
    !初期液滴中心
    real(8), parameter:: xc = 0.5d0*dble(xmax)	
    real(8), parameter:: yc = 0.5d0*dble(ymax)
    real(8), parameter:: zc = 0.5d0*dble(zmax)
    !phi計算用パラメータ/変数
    real(8),parameter:: phi1 = 2.638d-1
    real(8),parameter:: phi2 = 4.031d-1
    real(8),parameter:: a = 1.0d0
    real(8),parameter:: b = 1.0d0
    real(8),parameter:: T = 2.93d-1
    real(8),parameter:: kf = 0.01d0*ds**2
    real(8),parameter:: C = 0.0d0
    real(8),parameter:: M = (0.5d0-C/3.0d0)*ds !モビリティ
    real(8),dimension(0:xmax,0:ymax,0:zmax):: p0
    real(8),dimension(0:xmax,0:ymax,0:zmax):: lap_phi
    real(8),dimension(1:3, 0:xmax,0:ymax,0:zmax):: grad_phi
    real(8),dimension(1:3, 1:3, 0:xmax,0:ymax,0:zmax):: gphi
    real(8),dimension(1:3, 1:3, 0:xmax,0:ymax,0:zmax):: pcap
    real(8),dimension(1:3, 0:xmax,0:ymax,0:zmax):: div_pcap
    real(8) gtemp
    
    !速度・圧力計算用パラメータ/変数
    real(8),parameter:: kg = 0.06d0
    real(8) ftemp
    real(8),dimension(1:3,1:3,0:xmax,0:ymax,0:zmax):: grad_u
    
    !その他パラメータ
    real(8) phi_min, phi_max, p_in, p_out, veff, reff, sigma
    real(8) dp_cal, dp_th
    integer,parameter:: start = 5000
    integer,parameter:: stop = 90000

contains
    subroutine par(cx,cy,cz,cr,krone)
        integer, intent(in):: cx(15),cy(15),cz(15)
        real(8), intent(out):: cr(1:3,1:15), krone(1:3,1:3)
        do i=1,15
            cr(1,i) = dble(cx(i))
            cr(2,i) = dble(cy(i))
            cr(3,i) = dble(cz(i))
        enddo
        krone(1,1) = 1.d0; krone(1,2) = 0.d0; krone(1,3) = 0.d0
        krone(2,1) = 0.d0; krone(2,2) = 1.d0; krone(2,3) = 0.d0
        krone(3,1) = 0.d0; krone(3,2) = 0.d0; krone(3,3) = 1.d0
    end subroutine par

    !ディレクトリ作成
    subroutine mk_dirs(outdir)
        implicit none
        character(len=*), intent(in) :: outdir
        character(len=256) command
        write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
        call system(command)
    end subroutine mk_dirs

    function geq_cal(i,p, u1, u2, u3, ftemp) result(geq)
        integer,intent(in) :: i
        real(8),intent(in) :: p,u1,u2,u3,ftemp
        real(8) geq
        geq = E(i)*(3.0d0*p + 3.0d0*(u1*dble(cx(i)) + u2*dble(cy(i)) + u3*dble(cz(i))) &
                + 4.5d0*(u1*dble(cx(i)) + u2*dble(cy(i)) + u3*dble(cz(i)))**2 -1.5d0*(u1**2 + u2**2 + u3**2)) &
                + E(i)*kg*ftemp
    end function geq_cal
    
    
end module globals

program main
use globals
    implicit none
    real(8) feq(1:15, -1:xmax+1, -1:ymax+1, -1:zmax+1) !局所平衡分布関数（phi計算用）
    real(8) geq(1:15, -1:xmax+1, -1:ymax+1, -1:zmax+1) !局所平衡分布関数（速度・圧力計算用）
    real(8) phi(-1:xmax+1, -1:ymax+1, -1:zmax+1) !秩序パラメータ
    real(8) u1(-1:xmax+1, -1:ymax+1, -1:zmax+1), u2(-1:xmax+1, -1:ymax+1, -1:zmax+1),&
            u3(-1:xmax+1, -1:ymax+1, -1:zmax+1) !流速
    real(8) p(-1:xmax+1, -1:ymax+1, -1:zmax+1) !圧力
    real(8) Anu(0:xmax,0:ymax,0:zmax)
    real(8) nu(0:xmax,0:ymax,0:zmax)
    ! call mk_dirs(datadir)
    call par(cx,cy,cz,cr,krone)
    !初期条件
    do zi=0,zmax
        do yi=0,ymax
            do xi=0,xmax
                phi(xi,yi,zi) = phi1
                u1(xi,yi,zi) = 0.0d0
                u2(xi,yi,zi) = 0.0d0
                u3(xi,yi,zi) = 0.0d0
                p(xi,yi,zi) = 0.0d0

                if((dble(xi)*ds-xc)**2+(dble(yi)*ds-yc)**2+(dble(zi)*ds-zc)**2 <= (0.5d0*D)**2) then
                    phi(xi,yi,zi) = phi2
                endif

                nu(xi,yi,zi) = 0.5d0*(nu2-nu1)*(dsin(pi*(phi(xi,yi,zi)-0.5d0*(phi2+phi1))/(phi2-phi1)) + 1.d0) + nu1
                Anu(xi,yi,zi) = 4.5d0*(1.0d0/6.0d0-nu(xi,yi,zi)/ds)
            enddo
        enddo
    enddo

    !時間発展
    DO n=1,step
        !のりしろ境界
        do zi=0,zmax
            do yi=-1,ymax+1
                do xi=0,xmax
                    phi(xi,-1,zi) = phi(xi,ymax,zi)
                    u1(xi,-1,zi) = u1(xi,ymax,zi)
                    u2(xi,-1,zi) = u2(xi,ymax,zi)
                    u3(xi,-1,zi) = u3(xi,ymax,zi)
                    p(xi,-1,zi) = p(xi,ymax,zi)
                    
                    phi(xi,ymax+1,zi) = phi(xi,0,zi)
                    u1(xi,ymax+1,zi) = u1(xi,0,zi)
                    u2(xi,ymax+1,zi) = u2(xi,0,zi)
                    u3(xi,ymax+1,zi) = u3(xi,0,zi)
                    p(xi,ymax+1,zi) = p(xi,0,zi)
                    
                    phi(-1,yi,zi) = phi(xmax,yi,zi)
                    u1(-1,yi,zi) = u1(xmax,yi,zi)
                    u2(-1,yi,zi) = u2(xmax,yi,zi)
                    u3(-1,yi,zi) = u3(xmax,yi,zi)
                    p(-1,yi,zi) = p(xmax,yi,zi)
                    
                    phi(xmax+1,yi,zi) = phi(0,yi,zi)
                    u1(xmax+1,yi,zi) = u1(0,yi,zi)
                    u2(xmax+1,yi,zi) = u2(0,yi,zi)
                    u3(xmax+1,yi,zi) = u3(0,yi,zi)
                    p(xmax+1,yi,zi) = p(0,yi,zi)
                    
                enddo
            enddo
        enddo
        do yi=-1,ymax+1
            do xi=-1,xmax+1
                phi(xi,yi,-1) = phi(xi,yi,zmax)
                u1(xi,yi,-1) = u1(xi,yi,zmax)
                u2(xi,yi,-1) = u2(xi,yi,zmax)
                u3(xi,yi,-1) = u3(xi,yi,zmax)
                p(xi,yi,-1) = p(xi,yi,zmax)
                
                phi(xi,yi,zmax+1) = phi(xi,yi,0)
                u1(xi,yi,zmax+1) = u1(xi,yi,0)
                u2(xi,yi,zmax+1) = u2(xi,yi,0)
                u3(xi,yi,zmax+1) = u3(xi,yi,0)
                p(xi,yi,zmax+1) = p(xi,yi,0)
                
            enddo
        enddo
!=======================phiの計算============================================================
        !lap_phiの計算
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    lap_phi(xi,yi,zi) = -14.0d0*phi(xi,yi,zi)
                    do i=2,15
                        lap_phi(xi,yi,zi) = lap_phi(xi,yi,zi) +phi(xi+cx(i),yi+cy(i),zi+cz(i))
                    enddo
                    lap_phi(xi,yi,zi) = lap_phi(xi,yi,zi) / (5.0d0*ds**2)
                enddo
            enddo
        enddo
        !grad_phiの計算
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    do alpha=1,3
                        grad_phi(alpha,xi,yi,zi) = 0.0d0
                        do i = 2,15
                            grad_phi(alpha,xi,yi,zi) = grad_phi(alpha,xi,yi,zi) + cr(alpha,i)*phi(xi+cx(i),yi+cy(i),zi+cz(i))
                        enddo
                        grad_phi(alpha,xi,yi,zi) = grad_phi(alpha,xi,yi,zi) / (10.0d0*ds)
                    enddo
                enddo
            enddo
        enddo
        !p0の計算
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    p0(xi,yi,zi) = phi(xi,yi,zi)*T/(1.0d0-b*phi(xi,yi,zi)) - a*phi(xi,yi,zi)**2
                enddo
            enddo
        enddo
        !gphiの計算計算
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    do beta=1,3
                        do alpha=1,3
                            gphi(alpha,beta,xi,yi,zi) = 4.5d0*grad_phi(alpha,xi,yi,zi)*grad_phi(beta,xi,yi,zi)&
                            -1.5d0*(grad_phi(1,xi,yi,zi)**2+grad_phi(2,xi,yi,zi)**2+grad_phi(3,xi,yi,zi)**2)&
                            *krone(alpha,beta)
                        enddo
                    enddo
                enddo
            enddo
        enddo
        !pcapの計算
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    do beta=1,3
                        do alpha=1,3
                            pcap(alpha,beta,xi,yi,zi) = (p0(xi,yi,zi)-kf*phi(xi,yi,zi)*lap_phi(xi,yi,zi)&
                                    -0.5d0*kf*(grad_phi(1,xi,yi,zi)**2+grad_phi(2,xi,yi,zi)**2&
                                    +grad_phi(3,xi,yi,zi)**2))*krone(alpha,beta)&
                                    +kf*grad_phi(alpha,xi,yi,zi)*grad_phi(beta,xi,yi,zi)
                        enddo
                    enddo
                enddo
            enddo
        enddo
        !div_pcapの計算
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    do alpha=1,3
                        div_pcap(alpha,xi,yi,zi) = 0.0d0
                        do i=2,15
                            do beta=1,3
                                div_pcap(alpha,xi,yi,zi) = div_pcap(alpha,xi,yi,zi) + cr(beta,i)* &
                                                            pcap(alpha,beta,xi+cx(i),yi+cy(i),zi+cz(i))
                            enddo
                        enddo
                        div_pcap(alpha,xi,yi,zi) = div_pcap(alpha,xi,yi,zi) / (10.0d0*ds)
                    enddo
                enddo
            enddo
        enddo
        !feqの計算
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    do i=1,15
                        gtemp = 0.0d0
                        do beta=1,3
                            do alpha=1,3
                                gtemp = gtemp + gphi(alpha,beta,xi,yi,zi)*cr(alpha,i)*cr(beta,i)
                            enddo
                        enddo
                        feq(i,xi,yi,zi) = H(i)*phi(xi,yi,zi) &
                                        + F(i)*(p0(xi,yi,zi)-kf*phi(xi,yi,zi)*lap_phi(xi,yi,zi) &
                                        -kf/6.0d0*(grad_phi(1,xi,yi,zi)**2+grad_phi(2,xi,yi,zi)**2+ &
                                                grad_phi(3,xi,yi,zi)**2)) &
                                        + 3.0d0*E(i)*phi(xi,yi,zi)*(cr(1,i)*u1(xi,yi,zi)+ &
                                                cr(2,i)*u2(xi,yi,zi)+cr(3,i)*u3(xi,yi,zi)) &
                                        + E(i)*kf*gtemp &
                                        + E(i)*C*(div_pcap(1,xi,yi,zi)*cr(1,i)+ &
                                                div_pcap(2,xi,yi,zi)*cr(2,i)+div_pcap(3,xi,yi,zi)*cr(3,i))*ds
                    enddo
                enddo
            enddo
        enddo
        !のりしろ境界
        do zi=0,zmax
            do yi=-1,ymax+1
                do xi=0,xmax
                    do i=1,15
                        feq(i,xi,-1,zi) = feq(i,xi,ymax,zi)
                        feq(i,xi,ymax+1,zi) = feq(i,xi,0,zi)
                        feq(i,-1,yi,zi) = feq(i,xmax,yi,zi)
                        feq(i,xmax+1,yi,zi) = feq(i,0,yi,zi)
                    enddo
                enddo
            enddo
        enddo
        do yi=-1,ymax+1
            do xi=-1,xmax+1
                do i=1,15
                    feq(i,xi,yi,-1) = feq(i,xi,yi,zmax)
                    feq(i,xi,yi,zmax+1) = feq(i,xi,yi,0)
                enddo
            enddo
        enddo
        !phiの時間発展
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    phi(xi,yi,zi) = 0.0d0
                    do i=1,15
                        phi(xi,yi,zi) = phi(xi,yi,zi) + feq(i,xi-cx(i),yi-cy(i),zi-cz(i))
                    enddo
                enddo
            enddo
        enddo
!=====================速度・圧力の計算===============================================================
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    do beta=1,3
                        grad_u(1,beta,xi,yi,zi)=0.0d0
                        grad_u(2,beta,xi,yi,zi)=0.0d0
                        grad_u(3,beta,xi,yi,zi)=0.0d0
                        do i=2,15
                            grad_u(1,beta,xi,yi,zi)=grad_u(1,beta,xi,yi,zi)+cr(beta,i)*u1(xi+cx(i),yi+cy(i),zi+cz(i))
                            grad_u(2,beta,xi,yi,zi)=grad_u(2,beta,xi,yi,zi)+cr(beta,i)*u2(xi+cx(i),yi+cy(i),zi+cz(i))
                            grad_u(3,beta,xi,yi,zi)=grad_u(3,beta,xi,yi,zi)+cr(beta,i)*u3(xi+cx(i),yi+cy(i),zi+cz(i))
                        enddo
                        grad_u(1,beta,xi,yi,zi)=grad_u(1,beta,xi,yi,zi)/(10.0d0*ds)
                        grad_u(2,beta,xi,yi,zi)=grad_u(2,beta,xi,yi,zi)/(10.0d0*ds)
                        grad_u(3,beta,xi,yi,zi)=grad_u(3,beta,xi,yi,zi)/(10.0d0*ds)
                    enddo
                enddo
            enddo
        enddo
        !geqの計算
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    do i=1,15
                        gtemp = 0.0d0
                        ftemp = 0.0d0
                        do beta=1,3
                            do alpha=1,3
                                gtemp = gtemp + gphi(alpha,beta,xi,yi,zi)*cr(alpha,i)*cr(beta,i)

                                ftemp = ftemp + (grad_u(beta,alpha,xi,yi,zi)+grad_u(alpha,beta,xi,yi,zi))*cr(alpha,i)*cr(beta,i)
                            enddo
                        enddo
                        geq(i,xi,yi,zi) = E(i)*(3.0d0*p(xi,yi,zi) + 3.0d0*(u1(xi,yi,zi)*dble(cx(i)) &
                        + u2(xi,yi,zi)*dble(cy(i)) + u3(xi,yi,zi)*dble(cz(i))) &
                        + 4.5d0*(u1(xi,yi,zi)*dble(cx(i)) + u2(xi,yi,zi)*dble(cy(i)) + u3(xi,yi,zi)*dble(cz(i)))**2 &
                        -1.5d0*(u1(xi,yi,zi)**2 + u2(xi,yi,zi)**2 + u3(xi,yi,zi)**2) &
                        +Anu(xi,yi,zi)*ds*ftemp) &
                        + E(i)*kg*gtemp
                    enddo
                enddo
            enddo
        enddo
        !のりしろ境界
        do zi=0,zmax
            do yi=-1,ymax+1
                do xi=0,xmax
                    do i=1,15
                        geq(i,xi,-1,zi) = geq(i,xi,ymax,zi)
                        geq(i,xi,ymax+1,zi) = geq(i,xi,0,zi)
                        geq(i,-1,yi,zi) = geq(i,xmax,yi,zi)
                        geq(i,xmax+1,yi,zi) = geq(i,0,yi,zi)
                    enddo
                enddo
            enddo
        enddo
        do yi=-1,ymax+1
            do xi=-1,xmax+1
                do i=1,15
                    geq(i,xi,yi,-1) = geq(i,xi,yi,zmax)
                    geq(i,xi,yi,zmax+1) = geq(i,xi,yi,0)
                enddo
            enddo
        enddo

        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    p(xi,yi,zi) = 0.0d0
                    do i=1,15
                        p(xi,yi,zi) = p(xi,yi,zi) + geq(i,xi-cx(i),yi-cy(i),zi-cz(i)) / 3.0d0
                    enddo
                enddo
            enddo
        enddo
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    u1(xi,yi,zi) = 0.0d0
                    u2(xi,yi,zi) = 0.0d0
                    u3(xi,yi,zi) = 0.0d0
                    do i=1,15
                        u1(xi,yi,zi) = u1(xi,yi,zi) + dble(cx(i))*geq(i,xi-cx(i),yi-cy(i),zi-cz(i))
                        u2(xi,yi,zi) = u2(xi,yi,zi) + dble(cy(i))*geq(i,xi-cx(i),yi-cy(i),zi-cz(i))
                        u3(xi,yi,zi) = u3(xi,yi,zi) + dble(cz(i))*geq(i,xi-cx(i),yi-cy(i),zi-cz(i))
                    enddo
                enddo
            enddo
        enddo

        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    nu(xi,yi,zi) = 0.5d0*(nu2-nu1)*(dsin(pi*(phi(xi,yi,zi)-0.5d0*(phi2+phi1))/(phi2-phi1)) + 1.d0) + nu1
                    Anu(xi,yi,zi) = 4.5d0*(1.0d0/6.0d0-nu(xi,yi,zi)/ds)
                enddo
            enddo
        enddo

        !圧力差を求める phi1 = 2.638d-1 phi2 = 4.031d-1
        ! phi_min = phi1
        ! phi_max = phi2
        veff = 0.0d0
        do zi=0,zmax
            do yi=0,ymax
                do xi=0,xmax
                    if(phi(xi,yi,zi) >= 0.5d0*(phi1+phi2)) then
                        veff = veff + ds**3
                    end if
        !             if(phi_min > phi(xi,yi,zi))then
        !                 phi_min = phi(xi,yi,zi)
        !                 p_out = p(xi,yi,zi)
        !             endif

        !             if(phi_max < phi(xi,yi,zi))then
        !                 phi_max = phi(xi,yi,zi)
        !                 p_in = p(xi,yi,zi)
        !             endif
                enddo
            enddo
        enddo
        reff = (3.d0/4.d0*veff / pi)**(1.d0/3.d0)
        p_in = p(xmax/2,ymax/2,zmax/2)
        p_out = p(0,ymax/2,zmax/2)

        dp_cal = p_in - p_out
        sigma = dp_cal*reff/2.0d0
        write(*,*) n,dp_cal,veff,reff,sigma,kg


        
        ! write(*,*) "step =",n
        if ( mod(n+1,100)==0 ) then
            

            write(filename,*) n !i->filename 変換

            filename='step='//trim(adjustl(filename))//'.d' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける

            print *, filename !表示してみる
            open(20,file=filename, status='replace') !使う
            ! do xi=0,xmax
            !     write(20,*) xi,p(xi,ymax/2,zmax/2),u1(xi,ymax/2,zmax/2),u2(xi,ymax/2,zmax/2), &
            !     u3(xi,ymax/2,zmax/2),phi(xi,ymax/2,zmax/2),nu(xi,ymax/2,zmax/2)
            ! enddo
            zi = zmax/2
            do xi=0,xmax
                do yi=0,ymax
                    write(20,*) xi,yi,zi,u1(xi,yi,zi),u2(xi,yi,zi),u3(xi,yi,zi) &
                    ,sqrt(u1(xi,yi,zi)**2+u2(xi,yi,zi)**2+u3(xi,yi,zi)**2)
                enddo
                write(20,*)
            enddo
            close(20)
        endif
    ENDDO

    


    open(10,file="phi.txt")
    do xi=0,xmax
        write(10,*) xi,phi(xi,ymax/2,zmax/2)
    enddo
    close(10)
    

    
end program main