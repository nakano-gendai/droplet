!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   2次元の格子ボルツマン法のプログラム。                         !!!!
!!!!   格子気体モデルは2次元9速度。                                  !!!!
!!!!   計算領域は［0:xmax, 0:ymax］で定義している。                  !!!!
!!!!   なので、格子点数は (xmax+1)*(ymax+1) 個となる。               !!!!
!!!!   初期条件はテーラーグリーン渦の速度場および圧力場としている。  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module globals
        real(8),parameter:: ds = 1.0d0 !格子間隔（lattice unit）
        integer,parameter:: xmax = 39 !ｘ方向の1辺の長さ
        integer,parameter:: ymax = 39 !ｙ方向の1辺の長さ
        integer,parameter:: step = (xmax+1)*(ymax+1) !計算時間
    
        !粒子速度（整数）
        integer,parameter:: cx(9) = (/0, 1, 0, -1, 0, 1, -1, -1, 1/)
        integer,parameter:: cy(9) = (/0, 0, 1, 0, -1, 1, 1, -1, -1/)

        !重み付け関数
        real(8),parameter:: E(9) = (/ 4.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/36.0d0, &
                            1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0 /)

        !格子ボルツマン方程式のパラメータ
        real(8),parameter:: tau = 0.8d0 !緩和時間
        real(8),parameter:: nu = (tau - 0.5d0)/3.0d0 !液体の動粘性係数
        real(8),parameter:: ep = 1.0d0 / tau !緩和時間の逆数

        !物理量のパラメータ
        real(8),parameter:: rho = 1.0d0 !液体の密度
        real(8),parameter:: u0 = 0.01d0 !流速の初期値（音速よりも十分に小さくする。感覚的には0.1より小さく。）
        real(8),parameter:: p0 = rho / 3.0d0 !圧力の初期値
        real(8),parameter:: pi = acos(-1.0d0)　!円周率π（3.14...）
        real(8),parameter:: k1 = 2.0d0*pi/dble(xmax+1) !初期条件で使う波数。
        real(8),parameter:: k2 = 2.0d0*pi/dble(ymax+1) !初期条件で使う波数。

        !インデックス
        integer i, n, xi, yi

        !エラーとかを評価する変数（outputで使う。）
        real(8) err_u1,err_up_u1,err_down_u1,err_u2,err_up_u2,err_down_u2, &
                err_p,err_up_p,err_down_p
        real(8) ut1,ut2,pt
    contains
        subroutine parv(cr)
            real(8), intent(out):: cr(1:2,1:9)
                do i=1,9
                    cr(1,i) = dble(cx(i))
                    cr(2,i) = dble(cy(i))
                enddo
        end subroutine parv
    
        subroutine initial(p,u1,u2,f)
            real(8), intent(inout):: p(-1:xmax+1, -1:ymax+1),u1(-1:xmax+1, -1:ymax+1),&
                                u2(-1:xmax+1, -1:ymax+1),f(1:9, -1:xmax+1, -1:ymax+1)
            do yi=0,ymax
                do xi=0,xmax
                    p(xi,yi) = p0 -rho/4.0d0*(u0**2*sin(2.0d0*k1*dble(xi))+ &
                                k1**2/k2**2*u0**2*cos(2.0d0*k2*dble(yi)))
                    u1(xi,yi) = -u0*cos(k1*dble(xi))*sin(k2*dble(yi))
                    u2(xi,yi) = (k1/k2)*u0*sin(k1*dble(xi))*cos(k2*dble(yi))
                    do i =1,9
                        f(i,xi,yi) = E(i)*(3.0d0*p(xi,yi)+3.0d0*(dble(cx(i))*u1(xi,yi)+dble(cy(i))*u2(xi,yi)) &
                                +4.5d0*(dble(cx(i))*u1(xi,yi)+dble(cy(i))*u2(xi,yi))**2 -1.5d0*(u1(xi,yi)**2+u2(xi,yi)**2))
                    enddo
                enddo
            enddo
        end subroutine initial

        subroutine periodic(u1,u2,p,f,feq)
            real(8), intent(inout):: p(-1:xmax+1, -1:ymax+1),u1(-1:xmax+1, -1:ymax+1),u2(-1:xmax+1, -1:ymax+1),f(1:9, -1:xmax+1, -1:ymax+1),feq(1:9, -1:xmax+1, -1:ymax+1)
            do yi=-1,ymax+1
                do xi=0,xmax
                    u1(xi,-1) = u1(xi,ymax)
                    u2(xi,-1) = u2(xi,ymax)
                    p(xi,-1) = p(xi,ymax)
                    u1(xi,ymax+1) = u1(xi,0)
                    u2(xi,ymax+1) = u2(xi,0)
                    p(xi,ymax+1) = p(xi,0)
                    u1(-1,yi) = u1(xmax,yi)
                    u2(-1,yi) = u2(xmax,yi)
                    p(-1,yi) = p(xmax,yi)
                    u1(xmax+1,yi) = u1(0,yi)
                    u2(xmax+1,yi) = u2(0,yi)
                    p(xmax+1,yi) = p(0,yi)
                    do i=1,9
                        f(i,xi,-1) = f(i,xi,ymax)
                        f(i,xi,ymax+1) = f(i,xi,0)
                        f(i,-1,yi) = f(i,xmax,yi)
                        f(i,xmax+1,yi) = f(i,0,yi)

                        feq(i,xi,-1) = feq(i,xi,ymax)
                        feq(i,xi,ymax+1) = feq(i,xi,0)
                        feq(i,-1,yi) = feq(i,xmax,yi)
                        feq(i,xmax+1,yi) = feq(i,0,yi)
                    enddo
                enddo
            enddo
        end subroutine periodic

        subroutine renew(f,fnext)
            real(8),intent(in) :: fnext(1:9,-1:xmax+1,-1:ymax+1)
            real(8),intent(out) :: f(1:9,-1:xmax+1,-1:ymax+1)
            do yi=0,ymax
                do xi=0,xmax
                    do i=1,9
                        f(i,xi,yi) = fnext(i,xi,yi)
                    enddo
                enddo
            enddo
        end subroutine renew

        subroutine feq_cal(feq,u1,u2,p)
            real(8),intent(out) :: feq(1:9, -1:xmax+1, -1:ymax+1)
            real(8),intent(inout) :: u1(-1:xmax+1, -1:ymax+1),u2(-1:xmax+1, -1:ymax+1),p(-1:xmax+1, -1:ymax+1)
            do yi=0,ymax
                do xi=0,xmax
                    do i=1,9
                        feq(i,xi,yi) = E(i)*(3.0d0*p(xi,yi) + 3.0d0*(u1(xi,yi)*dble(cx(i)) + u2(xi,yi)*dble(cy(i))) &
                                    + 4.5d0*(u1(xi,yi)*dble(cx(i)) + u2(xi,yi)*dble(cy(i)))**2 -1.5d0*(u1(xi,yi)**2 + u2(xi,yi)**2))
                    enddo
                enddo
            enddo
        end subroutine feq_cal

        subroutine f_cal(fnext,f,feq)
            real(8),intent(out) :: fnext(1:9, -1:xmax+1, -1:ymax+1)
            real(8),intent(inout) :: f(1:9, -1:xmax+1, -1:ymax+1),feq(1:9, -1:xmax+1, -1:ymax+1)
            do yi=0,ymax
                do xi=0,xmax
                    do i=1,9
                        fnext(i,xi,yi) = f(i,xi-cx(i),yi-cy(i)) &
                                        - ep * ( f(i,xi-cx(i),yi-cy(i)) - feq(i,xi-cx(i),yi-cy(i)) )
                    enddo
                enddo
            enddo
        end subroutine f_cal

        subroutine physics_cal(f,u1,u2,p)
            real(8),intent(in) :: f(1:9,-1:xmax+1,-1:ymax+1)
            real(8),intent(out) :: u1(-1:xmax+1, -1:ymax+1),u2(-1:xmax+1, -1:ymax+1),p(-1:xmax+1, -1:ymax+1)
            do yi=0,ymax
                do xi=0,xmax
                    p(xi,yi) = 0.0d0
                    u1(xi,yi) = 0.0d0
                    u2(xi,yi) = 0.0d0
                    do i =1,9
                        p(xi,yi) = p(xi,yi) + f(i,xi,yi)
                        u1(xi,yi) = u1(xi,yi) + f(i,xi,yi) * dble(cx(i))
                        u2(xi,yi) = u2(xi,yi) + f(i,xi,yi) * dble(cy(i))
                    enddo
                    p(xi,yi) = p(xi,yi) / 3.0d0
                enddo
            enddo
        end subroutine physics_cal

        subroutine output(u1,u2,p)
            real(8),intent(in) :: u1(-1:xmax+1, -1:ymax+1), u2(-1:xmax+1, -1:ymax+1)
            real(8),intent(in) :: p(-1:xmax+1, -1:ymax+1)
            open(10,file = "output01.d")
            do xi=0,xmax
                do yi=0,ymax
                    write(10,*) xi,yi,u1(xi,yi),u2(xi,yi),p(xi,yi)
                enddo
            enddo
            close(10)

            open(11,file = "output02.d")
            n = step
            do xi=0,xmax
                do yi=0,ymax
                    ut1 = -u0*exp(-nu*dble(n)*(k1**2+k2**2))*cos(k1*dble(xi))*sin(k2*dble(yi))
                    ut2 = k1/k2*u0*exp(-nu*n*(k1**2+k2**2))*sin(k1*dble(xi))*cos(k2*dble(yi))
                    pt = p0 -rho/4.0d0*exp(-2.0d0*nu*n*(k1**2+k2**2))*(u0**2*sin(2.0d0*k1*dble(xi))+ &
                    k1**2/k2**2*u0**2*cos(2.0d0*k2*dble(yi)))
                    write(11,*) xi,yi,ut1,ut2,pt
                enddo
            enddo
            close(11)

            open(12,file = "output03.d")
            err_u1 = 0.0d0
            err_up_u1 = 0.0d0
            err_down_u1 = 0.0d0
            err_u2 = 0.0d0
            err_up_u2 = 0.0d0
            err_down_u2 = 0.0d0
            err_p = 0.0d0
            err_up_p = 0.0d0
            err_down_p = 0.0d0
            n = step
            do xi=0,xmax
                do yi=0,ymax
                    ut1 = -u0*exp(-nu*dble(n)*(k1**2+k2**2))*cos(k1*dble(xi))*sin(k2*dble(yi))
                    ut2 = k1/k2*u0*exp(-nu*dble(n)*(k1**2+k2**2))*sin(k1*dble(xi))*cos(k2*dble(yi))
                    pt = p0 -rho/4.0d0*exp(-2.0d0*nu*dble(n)*(k1**2+k2**2))*(u0**2*sin(2.0d0*k1*dble(xi))+ &
                    k1**2/k2**2*u0**2*cos(2.0d0*k2*dble(yi)))
                    err_up_u1 = err_up_u1 + abs(ut1 - u1(xi,yi))
                    err_down_u1 = err_down_u1 + abs(ut1)
                    err_up_u2 = err_up_u2 + abs(ut2 - u2(xi,yi))
                    err_down_u2 = err_down_u2 + abs(ut2)
                    err_up_p = err_up_p + abs(pt - p(xi,yi))
                    err_down_p = err_down_p + abs(pt)
                enddo
            enddo
            err_u1 = (err_up_u1 / err_down_u1)
            err_u2 = (err_up_u2 / err_down_u2)
            err_p = (err_up_p / err_down_p)
            write(12,*) err_u1,err_u2,err_p
            close(12)
        end subroutine output
    
end module globals
    
program main
use globals
    implicit none
    real(8) f(1:9, -1:xmax+1, -1:ymax+1), fnext(1:9, -1:xmax+1, -1:ymax+1) !速度分布関数
    real(8) feq(1:9, -1:xmax+1, -1:ymax+1) !局所平衡分布関数
    real(8) p(-1:xmax+1, -1:ymax+1)  !圧力
    real(8) u1(-1:xmax+1, -1:ymax+1), u2(-1:xmax+1, -1:ymax+1) !流速
    real(8):: cr(1:2, 1:9)  !粒子速度（実数）
    !粒子速度を倍精度へと変換（必要があればコメントアウトをはずす。）
    ! call parv(cr)
    !初期値の設定
    call initial(p,u1,u2,f)
!===========時間発展=====================
    DO n=1,step
        !局所平衡分布関数の計算
        call feq_cal(feq,u1,u2,p)
        !周期境界条件
        call periodic(u1,u2,p,f,feq)
        !次の時刻の速度分布関数fnextを計算
        call f_cal(fnext,f,feq)
        !速度分布関数の更新
        call renew(f,fnext)
        !流速の計算
        call physics_cal(f,u1,u2,p)
        !何ステップ計算できているかの確認（標準出力）
        write(*,*) "step = ", n
    ENDDO
!=============出力======================
    call output(u1,u2,p)
end program main
