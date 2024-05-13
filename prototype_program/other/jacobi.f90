program main
    implicit none
    integer, parameter :: n = 3
    real(8) a(n,n), w(n,n), w_t(n,n), a_tmp(n,n), t0(n,n), t_tmp(n,n)
    real(8) n1(n), n2(n), n3(n)
    real(8) mx, theta
    integer i,j,k
    integer p,q
    integer check, number
    !計算したい行列
    a(1,1:n) = (/3.0d0, 1.0d0, 1.0d0/)
    a(2,1:n) = (/1.0d0, 2.0d0, 0.0d0/)
    a(3,1:n) = (/1.0d0, 0.0d0, 2.0d0/)

    do j=1,n
        do i=1,n
            t0(i,j) = 0.0d0
            !対角成分だけ1.0を代入
            if(i == j) then
                t0(i,j) = 1.0d0  
            endif
        enddo
    enddo

    number = 0
    do 
        number = number + 1
        !非対角成分の中から絶対値の最大成分を探す
        mx = abs(a(1,2))
        p = 1
        q = 2
        do j = 2,n
            do i = 1,j-1
                if(abs(a(i,j)) > mx) then
                    mx = abs(a(i,j))
                    p = i
                    q = j
                endif
            enddo
        enddo
        !thetaを求める
        theta = 0.5d0 * atan(2.0d0 * a(p,q) / (a(p,p) - a(q,q)))
        !（転置）直交行列を求める
        do j=1,n
            do i=1,n
                w(i,j) = 0.0d0
                !対角成分だけ1.0を代入
                if(i == j) then
                    w(i,j) = 1.0d0  
                endif
            enddo
        enddo
        w(p,p) = cos(theta)
        w(p,q) = -sin(theta)
        w(q,p) = sin(theta)
        w(q,q) = cos(theta)
        do j=1,n
            do i=1,n
                w_t(i,j) = w(j,i)
            enddo
        enddo
        !新しい行列aを計算する
        a_tmp(:,:) = 0.0d0
        a_tmp(:,:) = matmul(w_t, a)
        a(:,:) = matmul(a_tmp, w)
        write(*,*) number
        do j=1,n
            write(*,*) a(1:3,j)
        enddo
        write(*,*)
        !並行して固有ベクトルを計算するための行列も計算する
        t_tmp(:,:) = matmul(t0, w)
        t0(:,:) = t_tmp(:,:)
        !非対角成分がすべて0.0に近づいたかを調べて、収束判定する
        check = 0
        do j=1,n
            do i=j+1,n
                if(abs(a(i,j)) < 1.0d-6) then
                    check = check + 1
                endif
            enddo
        enddo
        if(check == (n * n - n) / 2) exit
    enddo
    n1(1:3) = t0(1:3,1)
    n2(1:3) = t0(1:3,2)
    n3(1:3) = t0(1:3,3)
    
    do i=1,n
        write(*,*) a(i,i)
    enddo
    write(*,*)
    do i=1,n
        write(*,*) t0(i,1:n)
    enddo
    write(*,*)
    write(*,*)
    write(*,*) n1(1:3)
    write(*,*)
    write(*,*) n2(1:3)
    write(*,*)
    write(*,*) n3(1:3)


end program main