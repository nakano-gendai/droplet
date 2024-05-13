program main
implicit none
    real(8) average_value
    real(8) dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7
    integer i
    integer, parameter:: data_num = 316
    real(8) f(data_num),g(data_num),fg(data_num)


    open(10, file="taylor_para.d")
    ! open(10, file="kolmogorov.d")
    do i = 1, data_num
        read(10,*) dummy1, dummy2,  f(i), g(i), dummy3
    enddo
    close(10)

    average_value = 0.0d0
    do i = 1, data_num
        average_value = average_value + f(i)
    enddo
    average_value = average_value / data_num
    write(*,*) average_value

    average_value = 0.0d0
    do i = 1, data_num
        average_value = average_value + g(i)
    enddo
    average_value = average_value / data_num
    write(*,*) average_value

    ! average_value = 0.0d0
    ! do i = 1, data_num
    !     average_value = average_value + fg(i)
    ! enddo
    ! average_value = average_value / data_num
    ! write(*,*) average_value
    
end program