!****************************************************
! 正規乱数生成（by ボックス＝ミューラー法）
! : 0 〜 20 の正規乱数をボックス＝ミューラー法により
!   計算し、ヒストグラムを出力する。
!
!   date          name            version
!   2018.12.03    mk-mode.com     1.00 新規作成
!
! Copyright(C) 2018 mk-mode.com All Rights Reserved.
!****************************************************
!
module const
  implicit none
  ! SP: 単精度(4), DP: 倍精度(8)
  integer,     parameter :: SP = kind(1.0)
  integer(SP), parameter :: DP = selected_real_kind(2 * precision(1.0_SP))
  integer(SP), parameter :: M  = 10                     ! 平均
  real(DP),    parameter :: S  = 2.5_DP                 ! 標準偏差
  integer(SP), parameter :: N  = 10000                  ! 発生させる乱数の個数
  real(DP),    parameter :: PI = 4.0_DP * atan(1.0_DP)  ! 円周率
  real(DP),    parameter :: SC = N / 100.0_DP           ! ヒストグラム用スケール
end module const

module box_muller
  use const, only : SP, DP, S, M, PI
  implicit none
  private
  public :: rnd_seed, rnd

contains
  ! 乱数の種の設定
  subroutine rnd_seed
    implicit none
    integer(SP) :: seed_size, clock
    integer(SP), allocatable :: seed(:)

    call system_clock(clock)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = clock
    call random_seed(put=seed)
    deallocate(seed)
  end subroutine rnd_seed

  ! 正規乱数生成
  !
  ! :param(out) integer(4) r(2)
  subroutine rnd(r)
    integer(SP), intent(out) :: r(2)
    real(DP) :: r_u(2)  ! [0, 1) の一様乱数 2 個

    call random_number(r_u(1))
    call random_number(r_u(2))
    r(1) = int(S * sqrt(-2 * log(r_u(1))) * cos(2 * PI * r_u(2)) + M)
    r(2) = int(S * sqrt(-2 * log(r_u(1))) * sin(2 * PI * r_u(2)) + M)
  end subroutine rnd
end module box_muller

program rndnum_bm
  use const, only : SP, N, M, SC
  use box_muller
  implicit none
  integer(SP) :: i, j, r(2), hist(0:M * 2)

  ! 乱数の種の設定
  ! (正規乱数生成に使用する一様乱数の種)
  call rnd_seed

  ! 正規乱数の生成
  hist(:) = 0
  do i = 1 , N
    call rnd(r)
    hist(r) = hist(r) + 1
  end do

  ! 結果出力
  do i = 0, M * 2
    write (*, '(I3, ":", I4, " | ", A)') &
      & i, hist(i), repeat("*", int(hist(i) / SC))
  end do
end program rndnum_bm