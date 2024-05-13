program main
implicit none
    real(8) average_value, dis_ave
    real(8) dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7
    integer i
    integer, parameter:: data_num = 8
    integer, parameter:: data_end = 23
    real(8) f(data_num),g(data_num),fg(data_num)
    real(8) epsilon, Dv, taul, tauk, kolmogorov, taud, We, Dd
    real(8) tb1, tb2, tb3, tb4, tb5, tb6, tb7, tb8

    epsilon = 7.19d-11
    Dv = 127.5d0 !最大渦のスケール
    kolmogorov = 1.02d0 !コルモゴロフスケール

    !液滴直径
    Dd = 30.0d0

    !分裂時間
    tb1 = 30000.0d0
    tb2 = 35000.0d0
    tb3 = 32000.0d0
    tb4 = 31000.0d0
    tb5 = 31000.0d0
    tb6 = 31000.0d0
    tb7 = 30000.0d0
    tb8 = 29000.0d0

    taul = (epsilon)**(-1.0d0/3.0d0)*Dv**(2.0d0/3.0d0)
    tauk = (epsilon)**(-1.0d0/3.0d0)*kolmogorov**(2.0d0/3.0d0)
    taud = (epsilon)**(-1.0d0/3.0d0)*Dd**(2.0d0/3.0d0)


    open(113,file="D30.d")
    write(113,*) "## We   tb/taul   tb/tauk   tb/taud" 
    write(113,"(4es16.8)") 1.5d0, tb1/taul, tb1/tauk, tb1/taud
    write(113,"(4es16.8)") 3.0d0, tb2/taul, tb2/tauk, tb2/taud
    write(113,"(4es16.8)") 5.0d0, tb3/taul, tb3/tauk, tb3/taud
    write(113,"(4es16.8)") 10.0d0, tb4/taul, tb4/tauk, tb4/taud
    write(113,"(4es16.8)") 15.0d0, tb5/taul, tb5/tauk, tb5/taud
    write(113,"(4es16.8)") 20.0d0, tb6/taul, tb6/tauk, tb6/taud
    write(113,"(4es16.8)") 25.0d0, tb7/taul, tb7/tauk, tb7/taud
    write(113,"(4es16.8)") 30.0d0, tb8/taul, tb8/tauk, tb8/taud
    close(113)


    ! open(112,file="D70_2.d")
    ! write(112,"(5es16.8)") We, tb/taul, tb/tauk, tb/taud, tb/tauh
    ! close(112)

    ! open(10, file="taylor_para.d")
    ! do i = 1, data_num
    !     read(10,*) dummy1, dummy2, f(i), g(i), dummy3
    ! enddo
    ! close(10)

    ! open(11, file="taylor_re.d")
    ! do i = 1, data_num
    !     write(11,*) f(i)/286085.8d0, g(i)
    ! enddo
    ! close(11)

    ! average_value = 0.0d0
    ! do i = 1, data_num
    !     average_value = average_value + f(i)
    ! enddo
    ! average_value = average_value / dble(data_num - 252 + 1)
    ! write(*,*) "lambda = ", average_value

    ! average_value = 0.0d0
    ! do i = 1, data_num
    !     average_value = average_value + g(i)
    ! enddo
    ! average_value = average_value / dble(data_num)
    ! write(*,*) "Re = ", average_value

    ! do i = 1, data_num
    !     fg(i) = g(i) - average_value
    ! enddo

    ! dis_ave = 0.0d0
    ! do i = 1, data_num
    !     dis_ave = dis_ave + fg(i)
    ! enddo
    ! dis_ave = dis_ave / dble(data_num)
    ! write(*,*) "dis ave = ", dis_ave

    ! open(11, file="epsilon_dis.d")
    ! do i = 1, data_num
    !     write(11,*) f(i)*0.1d0/255.5d0, fg(i)
    ! enddo
    ! close(11)

    ! average_value = 0.0d0
    ! do i = 1, data_num
    !     average_value = average_value + fg(i)
    ! enddo
    ! average_value = average_value / data_num
    ! write(*,*) average_value
    
end program