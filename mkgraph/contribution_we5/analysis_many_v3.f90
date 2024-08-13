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
    ! character(*),parameter :: datadir_input = "/data/sht/nakanog/IHT_drop_d70_we2/contribution/"
    ! character(*),parameter :: datadir_output = "/data/sht/nakanog/IHT_drop_d70_we2/contribution/ave/"
    character(*),parameter :: datadir_input = "./ave_to_breaktime/"
    character(*),parameter :: datadir_input2 = "./sum_each_time/"
    character(*),parameter :: datadir_input3 = "./"
    character(*),parameter :: datadir_output = "./ensemble_ave_to_breaktime/"
    character(*),parameter :: datadir_output2 = "./kikakuka/"
    character(*),parameter :: datadir_output3 = "./kikakuka2/"

    real(8),parameter:: D = 70.0d0
    real(8),parameter:: epsilon = 9.15d-10 !エネルギー散逸率
    real(8),parameter:: eddytime = epsilon**(-1.0d0/3.0d0)*D**(2.0d0/3.0d0)
    real(8),parameter:: k_d = 0.5d0*dble(xmax) / D

    real(8),parameter:: We = 5.0d0
    real(8),parameter:: sigma = epsilon**(2.0d0/3.0d0)*D**(5.0d0/3.0d0) / We
    real(8),parameter:: Dh = 0.725*sigma**(3.0d0/5.0d0)*epsilon**(-2.0d0/5.0d0)
    real(8),parameter:: k_h = 0.5d0*dble(xmax) / Dh

    real(8) contribution_each_scale(step_num2,(xmax+1)/2 + 1)
    real(8) contribution_each_scale_sum(step_num2,(xmax+1)/2 + 1)
    real(8) contribution_each_scale_time_sum((xmax+1)/2 + 1)
    real(8) contribution_each_scale_wave_sum(step_num2)
    real(8) ensemble((xmax+1)/2 + 1)

    integer case_num, step_num
    character :: filename*200
    character :: filename_input*200
    character :: filename_output*200
    integer k_analysis, k_step

    real(8) time(1:case_end_num)
    integer time_num(1:case_end_num)
    real(8) dummy1, dummy2

    !!!!!!!!!!!!!!!!!!!!!!!!!分裂時間までのアンサンブル平均!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ensemble(:) = 0.0d0
    DO case_num = case_initial_num, case_end_num
        contribution_each_scale_time_sum(:) = 0.0d0

        write(filename,*) case_num
        filename_input = datadir_input//trim(adjustl(filename))//'.d'
        open(9, file=filename_input)
        do k_analysis = 1, (xmax+1) / 2
            read(9,*) dummy1, dummy2, contribution_each_scale_time_sum(k_analysis)
        enddo
        close(9)

        do k_analysis = 1, (xmax+1) / 2
            ensemble(k_analysis) = ensemble(k_analysis) + contribution_each_scale_time_sum(k_analysis)
            write(*,*) contribution_each_scale_time_sum(k_analysis)
        enddo
    ENDDO 

    filename_output = datadir_output//'d70we5_ensemble.d'
    open(120, file=filename_output, form="formatted")
    do k_analysis = 1, (xmax+1) / 2
        write(120,*) dble(k_analysis) / k_d, dble(k_analysis) / k_h, ensemble(k_analysis) / dble(case_end_num)
    enddo
    close(120)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!それぞれの時刻で全寄与で無次元化!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DO case_num = case_initial_num, case_end_num
    !     contribution_each_scale(:,:) = 0.0d0

    !     write(filename,*) case_num
    !     filename_input = datadir_input3//trim(adjustl(filename))//'_contribution.bin'
    !     open(110, file=filename_input, form="unformatted")
    !     do k_analysis = 1, (xmax+1) / 2
    !         do k_step = 1, step_num2
    !             read(110) contribution_each_scale(k_step, k_analysis)
    !         enddo
    !     enddo
    !     close(110)

    !     write(filename,*) case_num
    !     filename_input = datadir_input2//trim(adjustl(filename))//'_wavesum_eachtime.d'
    !     open(111, file=filename_input, form="formatted")
    !         do k_step = 1, step_num2
    !             read(111,*) dummy1, dummy2, contribution_each_scale_wave_sum(k_step)
    !         enddo
    !     close(111)

        ! k_analysis = 1
        ! write(filename,*) case_num
        ! filename_output = datadir_output2//trim(adjustl(filename))//'_4.d'
        ! open(20,file=filename_output,status='replace')
        ! do k_step = 1, step_num2, 4
        !     write(20,"(11es16.8)") (1000.0d0*dble(k_step)-1000.0d0)/eddytime, contribution_each_scale(k_step, k_analysis)/contribution_each_scale_wave_sum(k_step), contribution_each_scale(k_step, k_analysis+1)/contribution_each_scale_wave_sum(k_step), contribution_each_scale(k_step, k_analysis+2)/contribution_each_scale_wave_sum(k_step), contribution_each_scale(k_step, k_analysis+3)/contribution_each_scale_wave_sum(k_step), contribution_each_scale(k_step, k_analysis+4)/contribution_each_scale_wave_sum(k_step), contribution_each_scale(k_step, k_analysis+5)/contribution_each_scale_wave_sum(k_step), contribution_each_scale(k_step, k_analysis+6)/contribution_each_scale_wave_sum(k_step), contribution_each_scale(k_step, k_analysis+7)/contribution_each_scale_wave_sum(k_step), contribution_each_scale(k_step, k_analysis+8)/contribution_each_scale_wave_sum(k_step)
        ! enddo
        ! close(20) 

    !     k_analysis = 1
    !     write(filename,*) case_num
    !     filename_output = datadir_output3//trim(adjustl(filename))//'_4.d'
    !     open(20,file=filename_output,status='replace')
    !     do k_step = 1, step_num2, 4
    !         write(20,"(11es16.8)") (1000.0d0*dble(k_step)-1000.0d0)/eddytime, contribution_each_scale(k_step, k_analysis)/abs(contribution_each_scale_wave_sum(k_step)), contribution_each_scale(k_step, k_analysis+1)/abs(contribution_each_scale_wave_sum(k_step)), contribution_each_scale(k_step, k_analysis+2)/abs(contribution_each_scale_wave_sum(k_step)), contribution_each_scale(k_step, k_analysis+3)/abs(contribution_each_scale_wave_sum(k_step)), contribution_each_scale(k_step, k_analysis+4)/abs(contribution_each_scale_wave_sum(k_step)), contribution_each_scale(k_step, k_analysis+5)/abs(contribution_each_scale_wave_sum(k_step)), contribution_each_scale(k_step, k_analysis+6)/abs(contribution_each_scale_wave_sum(k_step)), contribution_each_scale(k_step, k_analysis+7)/abs(contribution_each_scale_wave_sum(k_step)), contribution_each_scale(k_step, k_analysis+8)/abs(contribution_each_scale_wave_sum(k_step))
    !     enddo
    !     close(20) 

    ! ENDDO 

end program main
