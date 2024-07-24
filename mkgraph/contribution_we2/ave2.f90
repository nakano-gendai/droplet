program main
implicit none
    integer,parameter:: xmax = 256 !ｘ方向格子数
    integer,parameter:: ymax = 256 !ｙ方向格子数
    integer,parameter:: zmax = 256 !ｚ方向格子数
    real(8),allocatable ::  f(:), g(:), t(:), ff(:), gg(:), k(:), E(:)
    real(8) dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7
    integer i
    real(8) wa_f, ave_f, wa_g, ave_g, wa_ff, ave_ff, wa_gg, ave_gg
    integer,parameter:: num = 71
    integer,parameter:: start = 113
    real(8),parameter:: pi = acos(-1.0d0) !円周率
    real(8),parameter:: nu = 0.001d0

    real(8) integral_scale, integral_top, integral_bottom
    integer kazu

    real(8) f1(71), f2(71), f3(71), f4(71), f5(71), f6(71), f7(71), f8(71), f9(71), f10(71)
    real(8) g1(71), g2(71), g3(71), g4(71), g5(71), g6(71), g7(71), g8(71), g9(71), g10(71)

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

    ! open(12,file="enesupe_phi_30we1.4_4.d")
    ! do i = 1, xmax/2+1
    !     read(12,*) k(i), E(i)
    ! enddo
    ! close(12)

    ! open(13,file="energy_interface.d")
    ! do i = 1, 33
    !     read(13,*) t(i), f(i)
    ! enddo
    ! close(13)

    ! open(14,file="contribution_m.d")
    ! do i = 1, 33
    !     read(14,*) t(i), g(i)
    ! enddo
    ! close(14)

    ! open(15,file="contribution_s.d")
    ! do i = 1, 33
    !     read(15,*) t(i), ff(i)
    ! enddo
    ! close(15)


    ! open(20,file = "energy_interface.d")
    ! do i = 1, 33
    !     write(20,"(4es16.8)") t(i)-5000.0d0, f(i)
    ! enddo
    ! close(20)

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

    ! open(20,file = "enesupe_phi_30we1.4_5.d")
    ! do i = 1, xmax/2+1
    !     write(20,"(2es16.8)") k(i)*0.93d0/127.5d0, E(i)
    ! enddo
    ! close(20)



    ! open(10,file="d70we1.4_break_time.d")
    ! do i = 1, num
    !     read(10,*) dummy4, f(i)
    ! enddo
    ! close(10)

    ! wa_f = 0.0d0
    ! kazu = 0
    ! do i = 1, num
    !     if(f(i) > 0) then
    !         wa_f = wa_f + f(i)
    !         kazu = kazu + 1
    !     endif
    ! enddo
    ! wa_f = wa_f / dble(kazu)
    ! write(*,*) wa_f

    ! do i = 1, num
    !     wa_g = wa_f + (f(i) - wa_f)**2.0d0
    ! enddo
    ! wa_g = wa_g / dble(kazu)
    ! wa_g = sqrt(wa_g)
    ! write(*,*) wa_g


    g1(:) = 0.0d0
    g2(:) = 0.0d0
    g3(:) = 0.0d0
    g4(:) = 0.0d0
    g5(:) = 0.0d0
    g6(:) = 0.0d0
    g7(:) = 0.0d0
    g8(:) = 0.0d0
    g9(:) = 0.0d0
    g10(:) = 0.0d0
    open(10,file="1_10_ave1_4.d")
    do i = 1, 71
        read(10,*) f(i), f1(i), f2(i), f3(i), f4(i), f5(i), f6(i), f7(i), f8(i), f9(i), f10(i)
    enddo
    close(10)

    do i = 1, 71
        g1(i) = g1(i) + f1(i)
        g2(i) = g2(i) + f2(i)
        g3(i) = g3(i) + f3(i)
        g4(i) = g4(i) + f4(i)
        g5(i) = g5(i) + f5(i)
        g6(i) = g6(i) + f6(i)
        g7(i) = g7(i) + f7(i)
        g8(i) = g8(i) + f8(i)
        g9(i) = g9(i) + f9(i)
        g10(i) = g10(i) + f10(i)
    enddo

    open(11,file="1_10_ave6_9.d")
    do i = 1, 71
        read(11,"(11es16.8)") f(i), f1(i), f2(i), f3(i), f4(i), f5(i), f6(i), f7(i), f8(i), f9(i), f10(i)
    enddo
    close(11)

    do i = 1, 71
        g1(i) = g1(i) + f1(i)
        g2(i) = g2(i) + f2(i)
        g3(i) = g3(i) + f3(i)
        g4(i) = g4(i) + f4(i)
        g5(i) = g5(i) + f5(i)
        g6(i) = g6(i) + f6(i)
        g7(i) = g7(i) + f7(i)
        g8(i) = g8(i) + f8(i)
        g9(i) = g9(i) + f9(i)
        g10(i) = g10(i) + f10(i)
    enddo

    open(12,file="1_10_ave11_24.d")
    do i = 1, 71
        read(12,"(11es16.8)") f(i), f1(i), f2(i), f3(i), f4(i), f5(i), f6(i), f7(i), f8(i), f9(i), f10(i)
    enddo
    close(12)

    do i = 1, 71
        g1(i) = g1(i) + f1(i)
        g2(i) = g2(i) + f2(i)
        g3(i) = g3(i) + f3(i)
        g4(i) = g4(i) + f4(i)
        g5(i) = g5(i) + f5(i)
        g6(i) = g6(i) + f6(i)
        g7(i) = g7(i) + f7(i)
        g8(i) = g8(i) + f8(i)
        g9(i) = g9(i) + f9(i)
        g10(i) = g10(i) + f10(i)
    enddo

    do i = 1, 71
        g1(i) = g1(i) / 22.0d0
        g2(i) = g2(i) / 22.0d0
        g3(i) = g3(i) / 22.0d0
        g4(i) = g4(i) / 22.0d0
        g5(i) = g5(i) / 22.0d0
        g6(i) = g6(i) / 22.0d0
        g7(i) = g7(i) / 22.0d0
        g8(i) = g8(i) / 22.0d0
        g9(i) = g9(i) / 22.0d0
        g10(i) = g10(i) / 22.0d0
    enddo

    do i = 1, 71
        f1(i) = g1(i)
        f2(i) = g2(i)
        f3(i) = g3(i)
        f4(i) = g4(i)
        f5(i) = g5(i)
        f6(i) = g6(i)
        f7(i) = g7(i)
        f8(i) = g8(i)
        f9(i) = g9(i)
        f10(i) = g10(i)
    enddo

    open(13,file="1_10_ave1_24_we2.d")
    do i = 1, 71
        write(13,"(11es16.8)") f(i), f1(i), f2(i), f3(i), f4(i), f5(i), f6(i), f7(i), f8(i), f9(i), f10(i)
    enddo
    close(13)

end program