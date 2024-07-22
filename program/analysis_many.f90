program main
implicit none
    integer,parameter:: xmax = 255 !ｘ方向格子数（０から数える）
    integer,parameter:: ymax = 255 !ｙ方向格子数（０から数える）
    integer,parameter:: zmax = 255 !ｚ方向格子数（０から数える）
    integer,parameter:: step_start = 5000
    integer,parameter:: step_end = 75000
    integer,parameter:: step_bin = 1000
    integer,parameter:: step_num2 = (step_end - step_start) / step_bin + 1 
    integer,parameter:: case_initial_num = 11 !最初のケース番号
    integer,parameter:: case_end_num = 24 !最後のケース番号
    character(*),parameter :: datadir_input = "/data/sht/nakanog/IHT_drop_d70_we2/contribution/"
    character(*),parameter :: datadir_output = "/data/sht/nakanog/IHT_drop_d70_we2/contribution/ave/"

    real(8),parameter:: D = 70.0d0
    real(8),parameter:: epsilon = 9.15d-10 !エネルギー散逸率
    real(8),parameter:: eddytime = epsilon**(-1.0d0/3.0d0)*D**(2.0d0/3.0d0)

    real(8) contribution_each_scale(step_num2,(xmax+1)/2 + 1)
    real(8) contribution_each_scale_sum(step_num2,(xmax+1)/2 + 1)

    integer case_num, step_num
    character :: filename*200
    character :: filename_input*200
    character :: filename_output*200
    integer k_analysis, k_step

    contribution_each_scale_sum(:,:) = 0.0d0
    DO case_num = case_initial_num, case_end_num
        write(filename,*) case_num
        filename_input = datadir_input//trim(adjustl(filename))//'_contribution.bin'
        open(110, file=filename_input, form="unformatted")
        do k_analysis = 1, (xmax+1) / 2
            do k_step = 1, step_num2
                read(110) contribution_each_scale(k_step, k_analysis)
            enddo
        enddo
        close(110)

        do k_analysis = 1, (xmax+1) / 2
            do k_step = 1, step_num2
                contribution_each_scale_sum(k_step, k_analysis) = contribution_each_scale_sum(k_step, k_analysis) + contribution_each_scale(k_step, k_analysis)
            enddo
        enddo
    ENDDO

    do k_analysis = 1, (xmax+1) / 2
        do k_step = 1, step_num2
            ! contribution_each_scale_sum(k_step, k_analysis) = contribution_each_scale_sum(k_step, k_analysis) / dble(case_end_num-case_initial_num+1)
            contribution_each_scale_sum(k_step, k_analysis) = contribution_each_scale_sum(k_step, k_analysis)
        enddo
    enddo

    k_analysis = 1
    filename_output = datadir_output//'1_10_ave11_24.d'
    open(20,file=filename_output,status='replace')
    do k_step = 1, step_num2
        write(20,"(11es16.8)") (1000.0d0*dble(k_step)-1000.0d0)/eddytime, contribution_each_scale_sum(k_step, k_analysis), contribution_each_scale_sum(k_step, k_analysis+1), contribution_each_scale_sum(k_step, k_analysis+2), contribution_each_scale_sum(k_step, k_analysis+3), contribution_each_scale_sum(k_step, k_analysis+4), contribution_each_scale_sum(k_step, k_analysis+5), contribution_each_scale_sum(k_step, k_analysis+6), contribution_each_scale_sum(k_step, k_analysis+7), contribution_each_scale_sum(k_step, k_analysis+8), contribution_each_scale_sum(k_step, k_analysis+9)
    enddo
    close(20)  

end program main