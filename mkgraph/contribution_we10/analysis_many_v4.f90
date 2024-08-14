program main
implicit none
    integer,parameter:: xmax = 255 !ｘ方向格子数（０から数える）
    integer,parameter:: ymax = 255 !ｙ方向格子数（０から数える）
    integer,parameter:: zmax = 255 !ｚ方向格子数（０から数える）
    integer,parameter:: step_start = 5000
    integer,parameter:: step_end = 50000
    integer,parameter:: step_bin = 1000
    integer,parameter:: step_num2 = (step_end - step_start) / step_bin + 1 
    integer,parameter:: case_initial_num = 1 !最初のケース番号
    integer,parameter:: case_end_num = 50 !最後のケース番号
    character(*),parameter :: datadir_input = "./"
    character(*),parameter :: datadir_output = "./kikakuka3/"
    character(*),parameter :: datadir_output2 = "./sum_each_time/"

    real(8),parameter:: D = 70.0d0
    real(8),parameter:: epsilon = 9.15d-10 !エネルギー散逸率
    real(8),parameter:: eddytime = epsilon**(-1.0d0/3.0d0)*D**(2.0d0/3.0d0)
    real(8),parameter:: k_d = 0.5d0*dble(xmax) / D

    real(8) contribution_each_scale(step_num2,(xmax+1)/2+1)
    real(8) contribution_each_scale_time_sum((xmax+1)/2+1)
    real(8) contribution_each_scale_wave_sum(step_num2)

    real(8) contribution_each_scale_sum(step_num2,(xmax+1)/2 + 1)

    integer case_num, step_num
    character :: filename*200
    character :: filename_input*200
    character :: filename_output*200
    integer k_analysis, k_step

    real(8) time(1:50)
    integer time_num(1:50)
    real(8) wa
    integer dummy_int

    open(9,file="./d70we10_break_time.d")
    do case_num = 1, 50
        read(9,*) dummy_int, time(case_num)
    enddo
    close(9)

    do case_num = 1, 50
        time(case_num) = time(case_num)*eddytime
        time(case_num) = time(case_num) / 1000.0d0
        time_num(case_num) = int(time(case_num))
        if(time_num(case_num)==0) then
            time_num(case_num) = step_num2
        endif
    enddo

    DO case_num = case_initial_num, case_end_num
        contribution_each_scale(:,:) = 0.0d0

        write(filename,*) case_num
        filename_input = datadir_input//trim(adjustl(filename))//'_contribution.bin'
        open(110, file=filename_input, form="unformatted")
        do k_analysis = 1, (xmax+1) / 2
            do k_step = 1, step_num2
                read(110) contribution_each_scale(k_step, k_analysis)
            enddo
        enddo
        close(110)

        contribution_each_scale_time_sum(:) = 0.0d0
        contribution_each_scale_wave_sum(:) = 0.0d0
        do k_analysis = 1, (xmax+1) / 2
            do k_step = 1, time_num(case_num)
                contribution_each_scale_time_sum(k_analysis) = contribution_each_scale_time_sum(k_analysis) + contribution_each_scale(k_step, k_analysis)
            enddo
        enddo

        wa = 0.0d0
        do k_analysis = 1, (xmax+1) / 2
            contribution_each_scale_time_sum(k_analysis) = contribution_each_scale_time_sum(k_analysis) / dble(time_num(case_num))
            wa = wa + contribution_each_scale_time_sum(k_analysis)
        enddo

        do k_analysis = 1, (xmax+1) / 2
            do k_step = 1, step_num2
                contribution_each_scale_wave_sum(k_step) = contribution_each_scale_wave_sum(k_step) + contribution_each_scale(k_step, k_analysis)
            enddo
        enddo

        ! write(filename,*) case_num
        ! filename_output = datadir_output//trim(adjustl(filename))//'.d'
        ! open(120, file=filename_output, form="formatted")
        ! do k_analysis = 1, (xmax+1) / 2
        !     write(120,*) dble(k_analysis), dble(k_analysis) / k_d, contribution_each_scale_time_sum(k_analysis) / wa
        ! enddo
        ! close(120)

        ! write(filename,*) case_num
        ! filename_output = datadir_output2//trim(adjustl(filename))//'_wavesum_eachtime.d'
        ! open(130, file=filename_output, form="formatted")
        ! do k_step = 1, step_num2
        !     write(130,*) 1000.0d0*dble(k_step)-1000.0d0, (1000.0d0*dble(k_step)-1000.0d0)/eddytime , contribution_each_scale_wave_sum(k_step)
        ! enddo
        ! close(130)

        contribution_each_scale_sum(:,:) = 0.0d0
        do k_analysis = 1, (xmax+1) / 2
            do k_step = 1, step_num2
                contribution_each_scale_sum(k_step, k_analysis) = contribution_each_scale_sum(k_step, k_analysis) + contribution_each_scale(k_step, k_analysis)
            enddo
        enddo

        k_analysis = 1
        write(filename,*) case_num
        filename_output = datadir_output//trim(adjustl(filename))//'_4.d'
        open(20,file=filename_output,status='replace')
        do k_step = 1, step_num2, 4
            write(20,"(11es16.8)") (1000.0d0*dble(k_step)-1000.0d0)/eddytime, contribution_each_scale_sum(k_step, k_analysis)/wa, contribution_each_scale_sum(k_step, k_analysis+1)/wa, contribution_each_scale_sum(k_step, k_analysis+2)/wa, contribution_each_scale_sum(k_step, k_analysis+3)/wa, contribution_each_scale_sum(k_step, k_analysis+4)/wa, contribution_each_scale_sum(k_step, k_analysis+5)/wa, contribution_each_scale_sum(k_step, k_analysis+6)/wa, contribution_each_scale_sum(k_step, k_analysis+7)/wa, contribution_each_scale_sum(k_step, k_analysis+8)/wa, contribution_each_scale_sum(k_step, k_analysis+9)/wa
        enddo
        close(20) 

    ENDDO 

end program main
