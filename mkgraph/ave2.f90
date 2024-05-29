program main
implicit none
    integer,parameter:: xmax = 256 !ｘ方向格子数
    integer,parameter:: ymax = 256 !ｙ方向格子数
    integer,parameter:: zmax = 256 !ｚ方向格子数
    real(8),allocatable ::  f(:), g(:), t(:), ff(:), gg(:), k(:), E(:)
    real(8) dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7
    integer i
    real(8) wa_f, ave_f, wa_g, ave_g, wa_ff, ave_ff, wa_gg, ave_gg
    integer,parameter:: num = 180
    integer,parameter:: start = 1
    real(8),parameter:: pi = acos(-1.0d0) !円周率
    real(8),parameter:: nu = 0.001d0

    real(8) integral_scale, integral_top, integral_bottom

    allocate(t(num))
    allocate(f(num))
    allocate(g(num))
    allocate(ff(num))
    allocate(gg(num))
    allocate(k(xmax/2+1))
    allocate(E(xmax/2+1))

    ! open(10,file="kolmogorov.d")
    ! do i = 1, num
    !     read(10,*) t(i), dummy4, f(i), g(i), dummy1
    ! enddo
    ! close(10)

    ! open(11,file="taylor_para.d")
    ! do i = 1, num
    !     read(11,*) dummy1, dummy2, dummy3,  ff(i), gg(i)
    ! enddo
    ! close(11)

    open(12,file="enesupe_phi_30we1.4_4.d")
    do i = 1, xmax/2+1
        read(12,*) k(i), E(i)
    enddo
    close(12)

    ! wa_f = 0.0d0
    ! wa_g = 0.0d0
    ! wa_ff = 0.0d0
    ! wa_gg = 0.0d0
    ! do i = start, num
    !     wa_f = wa_f + f(i)
    !     wa_g = wa_g + g(i)
    !     wa_ff = wa_ff + ff(i)
    !     wa_gg = wa_gg + gg(i)
    ! enddo
    ! ave_f = wa_f / dble(num - start + 1)
    ! ave_g = wa_g / dble(num - start + 1)
    ! ave_ff = wa_ff / dble(num - start + 1)
    ! ave_gg = wa_gg / dble(num - start + 1)

    ! !積分長
    ! integral_top = 0.0d0
    ! integral_bottom = 0.0d0
    ! integral_scale = 0.0d0
    ! do i = 2, xmax/2+1
    !     integral_top = integral_top + E(i)/(dble(i)-1.0d0)
    !     integral_bottom = integral_bottom + E(i)
    ! enddo
    ! integral_scale = 3.0d0 / 4.0d0 * pi * integral_top / integral_bottom
    ! integral_scale = integral_scale / pi * dble(xmax) / 2.0d0

    ! write(*,*) "Re = ", ave_ff
    ! write(*,*) "u_rms = ", ave_gg

    ! write(*,*) "epsilon = ", ave_g
    ! write(*,*) "kolmogorov scale = ", ave_f

    ! write(*,*) "integral scale = ", integral_scale

    open(20,file = "enesupe_phi_30we1.4_5.d")
    do i = 1, xmax/2+1
        write(20,"(2es16.8)") k(i)*0.93d0/127.5d0, E(i)
    enddo
    close(20)

end program
