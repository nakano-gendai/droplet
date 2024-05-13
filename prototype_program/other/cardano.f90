program main
    implicit none
    real(8) e1,e2,e3,e4
    real(8) m1,m2,m3
    real(8) n1,n2
    real(8) discriminant
    real(8) principal_moment_1, principal_moment_2, principal_moment_3
    e1 = 1.0d0
    e2 = 2.0d0
    e3 = 2.0d0
    e4 = 1.0d0
    
    m1 = e2 / e1
    m2 = e3 / e1
    m3 = e4 / e1

    n1 = m2 - m1*m1/3.0d0
    n2 = 2.0d0*m1*m1*m1/27.0d0 - m1*m2/3.0d0 + m3

    discriminant = n2*n2/4.0d0 + n1*n1*n1/27.0d0
    if(discriminant < 0.0d0) then
        theta = atan2(sqrt(-discriminant), -n2/2.0d0)
        principal_moment_1 = 2.0d0 * sqrt(-n1/3.0d0) * cos(theta / 3.0d0) - m1 / 3.0d0
        principal_moment_2 = 2.0d0 * sqrt(-n1/3.0d0) * cos((theta + 2.0d0*pi) / 3.0d0) - m1 / 3.0d0
        principal_moment_3 = 2.0d0 * sqrt(-n1/3.0d0) * cos((theta + 4.0d0*pi) / 3.0d0) - m1 / 3.0d0

        principal_moment_1 = 1.0d0 / sqrt(principal_moment_1)
        principal_moment_2 = 1.0d0 / sqrt(principal_moment_2)
        principal_moment_3 = 1.0d0 / sqrt(principal_moment_3)
    endif

end program main