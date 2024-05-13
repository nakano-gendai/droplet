program main
    implicit none

    real(8), parameter:: bin = 0.2d0 !度数分布の刻み幅
    real(8),parameter :: rhoc = 1.0d0 !連続相の密度
    real(8),parameter :: sigma = 9.0d-6 !液滴の界面張力係数
    real(8),parameter :: energy_dissipation = 2.5d-12 !乱流場のエネルギー散逸率
    real(8),parameter :: integral_scale = 138.0d0 !積分長
    real(8),parameter :: kolmogorov_scale = 1.135d0 !乱流場のKolmogorovスケール
    real(8),parameter :: Dh = 0.725d0*(rhoc/sigma)**(-3.0d0/5.0d0)*energy_dissipation**(-2.0d0/5.0d0) !Hinzeスケール
    real(8),parameter :: Di = 0.585d0*(sigma/rhoc)*energy_dissipation**(-2.0d0/3.0d0)*integral_scale**(-2.0d0/3.0d0)
    real(8),parameter :: Dk = 0.585d0*(sigma/rhoc)*energy_dissipation**(-2.0d0/3.0d0)*kolmogorov_scale**(-2.0d0/3.0d0)
    real(8),parameter :: Dd = 255.5d0 !液滴の初期直径

    integer,parameter :: flag = 1 !0：Kolmogorovスケールで無次元化　1：Hinzeスケールで無次元化  2:無次元化なし
    integer,parameter :: file_start = 8585000
    integer,parameter :: file_end = 13255000
    integer,parameter :: file_bin = 10000
    integer,parameter :: file_kazu = (file_end - file_start) / file_bin + 1
    integer,parameter :: kazu_max = 300 !全ステップ数における最大液滴数（見積もり）
    integer,parameter :: data_index = kazu_max * file_kazu !サンプルする液滴の全データ数（見積もり）
    character(*),parameter :: datadir = "/data/sht/nakanog/taylor_512_dsd/collect/radius/"

    real(8) radius(data_index),  tmp(file_kazu, kazu_max)
    integer,allocatable :: kazu(:)
    real(8),allocatable :: frequency(:)

    integer i, j, kk, step, file_step
    real(8) max_radius
    integer kazu_index, sample_number
    character :: filename*200

    open(21,file = "./cal.d")
    file_step = 0
DO step = file_start, file_end, file_bin
    file_step = file_step + 1
    write(21,*) file_step
    write(filename,*) step !i->filename 変換
    filename=datadir//"radius_"//trim(adjustl(filename))//'.d' !adjustlで左寄せにしてからtrimで末尾の空白除去，拡張子等をくっつける
    open(10,file=filename)
    do i = 1, kazu_max
        read(10,*,end=999)  tmp(file_step, i)
    enddo
999 close(10)
ENDDO
    close(21)

    open(22,file="tasikame.d")
    radius(:) = 0.0d0
    do j = 1, file_kazu
        do i = 1, kazu_max
            radius(kazu_max * (j-1) + i) = tmp(j, i)
            write(22,*) tmp(j, i)
        enddo
    enddo
    close(22)

!========数分布を計算=============================
    do i = 1, data_index
        if(flag == 0) then
            radius(i) = 2.0d0*radius(i) / kolmogorov_scale ! Kolmogorovスケールで無次元化
        else if(flag == 1) then
            radius(i) = 2.0d0*radius(i) / Dh ! Hinzeスケールで無次元化
        else if(flag == 2) then
            radius(i) = 2.0d0*radius(i) ! 無次元化なし
        endif
    enddo

    do i = 1, data_index
        if(radius(i) > Dd) then
            radius(i) = 0.0d0
        endif
    enddo

    max_radius = radius(1)
    do i = 2, data_index
        if(radius(i) > max_radius) then
            max_radius = radius(i)
        endif
    enddo
    kazu_index = int(max_radius / bin) + 1
    allocate(kazu(kazu_index))
    allocate(frequency(kazu_index))

    kazu(:) = 0
    do i = 1, data_index
        if(radius(i) > 0) then
            kk = int(radius(i) / bin) + 1
            kazu(kk) = kazu(kk) + 1
        endif
    enddo

!=========度数分布を計算==========================
    sample_number = 0
    do i = 1, kazu_index
        sample_number = sample_number + kazu(i)
    enddo

    do i = 1, kazu_index
        frequency(i) = dble(kazu(i)) / dble(sample_number)
    enddo

!==================出力==========================
    open(20,file="./dsd1.d")
    do i = 1, kazu_index
        write(20,"(4es16.8)") bin*dble(i)-0.5d0*bin, frequency(i), dble(kazu(i)), dble(sample_number)
    enddo
    close(20)
end program main