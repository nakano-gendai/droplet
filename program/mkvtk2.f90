module globals
    integer,parameter :: nx = 512
    integer,parameter :: ny = 512
    integer,parameter :: nz = 512
    real(8),parameter :: dx = 1.0d0
    real(8),parameter :: dy = 1.0d0
    real(8),parameter :: dz = 1.0d0
    !読み込みディレクトリー
    character(*),parameter :: dir = "/data/sht/nakanog/taylor_512_drop_movie/small/collect/"
    ! character(*),parameter :: dir = "./"
    !出力ディレクトリー
    character(*),parameter :: datadir = "/data/sht/nakanog/taylor_512_drop_movie/small/collect/vtk/"
    ! character(*),parameter :: datadir = "./case17_eddy/"
    integer,parameter :: num = 100 !読み込むファイルの数
    integer,parameter :: num_vtk = 30
    character(8) file_num
    character :: filename*200
    integer i, j, k
contains
    subroutine read_files(scalar, step)
        implicit none
        real(8), intent(inout) :: scalar( 0:nx+1, 0:ny+1, 0:nz+1 )
        integer, intent(in) :: step
        real(8) x, y, z
        scalar(:,:,:) = 1.0d-10
        ! scalar(:,:,:) = 2.0d0
        if((step >= 0) .and. (step < 10)) then
            write(file_num, "(i1)") step
        elseif((step > 9) .and. (step < 100)) then
            write(file_num, "(i2)") step
        elseif((step > 99) .and. (step < 1000)) then
            write(file_num, "(i3)") step
        elseif((step > 999) .and. (step < 10000)) then
            write(file_num,"(i4)") step
        elseif((step > 9999) .and. (step < 100000)) then
            write(file_num,"(i5)") step
        elseif((step > 99999) .and. (step < 1000000)) then
            write(file_num,"(i6)") step
        elseif((step > 999999) .and. (step < 10000000)) then
            write(file_num,"(i7)") step
        endif
        ! open(20, file=dir//"small_2.bin", form='unformatted')
        open(20, file=dir//trim(file_num)//".bin", form='unformatted')
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    read(20) scalar(i,j,k)
                enddo
            enddo
        enddo
        close(20)
    end subroutine read_files

    subroutine output_scalar(scalar, stepnum)
        implicit none
        real(8), intent(in) :: scalar( 0:nx+1, 0:ny+1, 0:nz+1 )
        integer, intent(in) :: stepnum
        character(len=100) :: write_dataname
        character(len=120) :: buffer
        character :: lf
        lf = char(10)   ! Line feed character
        !------------------------------------------------------------------------------------------
        !   # Write vtk File #
        !------------------------------------------------------------------------------------------
        write_dataname='       '
        write(write_dataname,'(1i7)') stepnum
        do j=1,7
            if(write_dataname(j:j)==' ') write_dataname(j:j)='0'
        end do
        write_dataname = datadir//trim(write_dataname)//".vtk"
        print *, write_dataname
        ! write(write_dataname,"(a)") //trim(file_num)//".vtk"
        open(unit=num_vtk,file=trim(write_dataname), status="replace", form="unformatted", action="write", access="stream", convert="big_endian")  ! num_vtkは装置番号, endianの指定が必要
        write(buffer,"(a)") '# vtk DataFile Version 3.0'//lf
        write(num_vtk) trim(buffer)
        write(buffer,"(a)") 'hoge'//lf  ! hogeはなんでもいい
        write(num_vtk) trim(buffer)
        write(buffer,"(a)") 'BINARY'//lf
        write(num_vtk) trim(buffer)
        write(buffer,"(a)") 'DATASET STRUCTURED_POINTS'//lf
        write(num_vtk) trim(buffer)
        write(buffer,"(a,3(1x,i4))") 'DIMENSIONS',nx+2,ny+2,nz+2
        write(num_vtk) trim(buffer)
        write(buffer,"(a,3(1x,i3))") lf//'ORIGIN',0,0,0
        write(num_vtk) trim(buffer)
        write(buffer,"(a,3(1x,f16.10))") lf//'SPACING',dx,dy,dz
        write(num_vtk) trim(buffer)
        write(buffer,"(a,i10)") lf//'POINT_DATA ', (nx+2)*(ny+2)*(nz+2)
        write(num_vtk) trim(buffer)
        write(buffer,"(a)") lf//'SCALARS '//'hoge'//' double'//lf  ! hogeはなんでもいい
        write(num_vtk) trim(buffer)
        write(buffer,"(a)") 'LOOKUP_TABLE default'//lf
        write(num_vtk) trim(buffer)

        do k = 0, nz+1
            do j = 0, ny+1
                do i = 0, nx+1
                    write(num_vtk) scalar(i,j,k)
                enddo
            enddo
        enddo
        close(num_vtk)
    end subroutine output_scalar

    !ディレクトリ作成
    subroutine mk_dirs(outdir)
        implicit none
        character(len=*), intent(in) :: outdir
        character(len=256) command
        write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
        call system(command)
    end subroutine mk_dirs
end module globals

program main
use globals
    implicit none
    real(8) scalar( 0:nx+1, 0:ny+1, 0:nz+1 )
    integer step, stepnum
    call mk_dirs(datadir)
    stepnum = 1
    DO step = 5000, 100000, 5000
    ! step = 54000
        call read_files(scalar, step)
        write(*,*) "read"
        call output_scalar(scalar, stepnum)
        stepnum = stepnum + 1
    ENDDO
end program main
