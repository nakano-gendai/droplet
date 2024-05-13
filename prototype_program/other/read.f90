program main
    real(8),allocatable :: data1(:)
    real(8),allocatable :: data2(:)
    ! real(8),allocatable :: data1(:,:)
    ! real(8),allocatable :: data2(:,:)
    integer i, xmax, ymax, zmax, N, xi, yi, zi
    xmax = 299
    ymax = 104
    zmax = 104
    N = (xmax+1)*(ymax+1)*(zmax+1)
    allocate(data1(N))
    allocate(data2(N))
    ! allocate(data1(ymax+1,zmax+1))
    ! allocate(data2(ymax+1,zmax+1))
    !data1の読み込み
    open(10,file="10000.d")
    ! do zi=1,zmax+1
    !     do yi=1,ymax+1
    !         read(10,*) data1(yi,zi)
    !     enddo
    ! enddo
    do i=1,N
        read(10,*) data1(i)
    enddo
    close(10)
    !data2の読み込み
    open(11,file="10000.d")
    ! do zi=1,zmax+1
    !     do yi=1,ymax+1
    !         read(11,*) data2(yi,zi)
    !     enddo
    ! enddo
    do i=1,N
        read(11,*) data2(i)
    enddo
    close(11)
    open(20,file="re0.2ca0.1.d")
    ! do zi=1,zmax+1
    !     do yi=1,ymax+1
    !         ! write(20,*) yi, zi, abs(data2(yi,zi)-data1(yi,zi))
    !         if((abs(data1(yi,zi))<=1d-15).or.(abs(data2(yi,zi))<=1d-15)) then
    !             write(20,*) yi, zi, 0.0d0
    !         else
    !             write(20,*) yi, zi, abs(abs(data2(yi,zi))-abs(data1(yi,zi)))/abs(data1(yi,zi))
    !         endif
    !     enddo
    !     write(20,*) " "
    ! enddo
    do i=1,N
        write(20,*) i, abs(data2(i)-data1(i))/abs(data1(i))
    enddo
    close(20)
end program main