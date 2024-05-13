module subprogs
    implicit none
contains
    subroutine gauss_seidel(a, b, x, n, itrmax, er0)
        integer, intent(in) :: n, itrmax
        real(8), intent(in) :: a(n, n), b(n), er0
        real(8), intent(out) :: x(n)
        real(8) s, er, rd(n), r(n)
        integer i, itr
        do i = 1, n
            if(a(i,i) == 0.0d0) stop "a(i,i) == 0.0d0"
            rd(i) = 1.0d0 / a(i,i)
        enddo
        x(1:n) = 0.0d0
        do itr = 1, itrmax
            do i = 1, n
                s = dot_product(a(i,1:i-1), x(1:i-1))
                s = s + dot_product(a(i,i+1:n), x(i+1:n))
                x(i) = rd(i) * (b(i) - s)
            enddo
            r(1:n) = b(1:n) -matmul(a, x)
            er = dot_product(r,r)
            write(*,*) "itr =", itr, "err =", er
            if(er <= er0) then
                write(*,*) "# converged #"
                exit
            endif
        enddo
    end subroutine gauss_seidel
end module subprogs

program main
    use subprogs
    implicit none
    integer, parameter :: n = 3, itrmax = 100
    real(8) a(n,n), b(n), x(n) !a(n,n)は狭義対角優位ギョウレツ
    real(8) :: er0 = 1.0d-6
    a(1,1:n) = (/-1.0d0, 1.0d0, 1.0d0/)
    a(2,1:n) = (/0.0d0, 0.0d0, 0.0d0/)
    a(3,1:n) = (/0.0d0, 0.0d0, 0.0d0/)
    b(1:n) = (/0.0d0, 0.0d0, 0.0d0/)
    call gauss_seidel(a, b, x, n, itrmax, er0)
    write(*,*) x(1:n)
end program main