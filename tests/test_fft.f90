!-------------------------------------------------------------------------------
! 月行九道 Phase 2: FFT 验证程序
!
! 5 项独立验证:
!   T1 — Delta 函数测试 (N=8)
!   T2 — 单频正弦峰值检测 (N=64)
!   T3 — 与暴力 DFT 对比 (N=8)
!   T4 — 正向+逆向往返测试 (N=16)
!   T5 — 非 2 的幂输入检查 (N=7)
!
! 结果输出到 data/test_fft.log
!-------------------------------------------------------------------------------
program test_fft
    use mod_precision, only: dp
    use mod_fft, only: fft
    implicit none

    logical :: all_pass

    ! ---- 输出文件 ----
    integer, parameter :: FUNIT_LOG = 30

    ! ---- 时间字符串 ----
    character(len=8)  :: date_str
    character(len=10) :: time_str

    ! ---- 初始化日志 ----
    open(unit=FUNIT_LOG, file='data/test_fft.log', &
         access='sequential', form='formatted', status='replace')

    call date_and_time(date=date_str, time=time_str)
    write(FUNIT_LOG, '(a)') '# Lunar Nine Paths Phase 2 — FFT Module Verification'
    write(FUNIT_LOG, '(a)') &
        '# Generated: ' // date_str(1:4) // '-' // date_str(5:6) // '-' // &
        date_str(7:8) // ' ' // time_str(1:2) // ':' // time_str(3:4) // ':' // time_str(5:10)
    write(FUNIT_LOG, '(a)') '#'

    all_pass = .true.

    ! ==================================================================
    ! T1: Delta 函数
    ! ==================================================================
    write(FUNIT_LOG, '(a)') '=== T1: Delta function (N=8) ==='
    call test_t1(all_pass, FUNIT_LOG)

    ! ==================================================================
    ! T2: 单频正弦
    ! ==================================================================
    write(FUNIT_LOG, '(a)') ''
    write(FUNIT_LOG, '(a)') '=== T2: Single sinusoid peak detection (N=64) ==='
    call test_t2(all_pass, FUNIT_LOG)

    ! ==================================================================
    ! T3: 与暴力 DFT 对比
    ! ==================================================================
    write(FUNIT_LOG, '(a)') ''
    write(FUNIT_LOG, '(a)') '=== T3: Comparison with direct DFT (N=8) ==='
    call test_t3(all_pass, FUNIT_LOG)

    ! ==================================================================
    ! T4: 往返测试
    ! ==================================================================
    write(FUNIT_LOG, '(a)') ''
    write(FUNIT_LOG, '(a)') '=== T4: Round-trip FFT + IFFT (N=16) ==='
    call test_t4(all_pass, FUNIT_LOG)

    ! ==================================================================
    ! T5: 非 2 的幂
    ! ==================================================================
    write(FUNIT_LOG, '(a)') ''
    write(FUNIT_LOG, '(a)') '=== T5: Non-power-of-2 input (N=7) ==='
    call test_t5(all_pass, FUNIT_LOG)

    ! ==================================================================
    ! 汇总
    ! ==================================================================
    write(FUNIT_LOG, '(a)') ''
    if (all_pass) then
        write(FUNIT_LOG, '(a)') 'ALL TESTS PASSED'
    else
        write(FUNIT_LOG, '(a)') 'SOME TESTS FAILED'
    end if

    close(FUNIT_LOG)

    if (all_pass) then
        write(*, '(a)') 'All 5 FFT tests PASSED. See data/test_fft.log'
    else
        write(*, '(a)') 'Some FFT tests FAILED. See data/test_fft.log'
        stop 1
    end if

contains

    ! ==================================================================
    ! T1: Delta 函数测试
    ! 输入: x[0]=1, x[1..7]=0
    ! 预期: FFT 后所有元素 = (1, 0)
    ! ==================================================================
    subroutine test_t1(all_pass, logunit)
        logical, intent(inout) :: all_pass
        integer, intent(in)    :: logunit

        integer,  parameter :: NN = 8
        complex(dp) :: x(NN)
        integer :: info, i
        real(dp) :: max_diff

        x = (0.0_dp, 0.0_dp)
        x(1) = (1.0_dp, 0.0_dp)  ! Fortran 1-indexed: x(1) = delta at index 0

        call fft(x, info)

        if (info /= 0) then
            write(logunit, '(a, i0)') '  ERROR: fft returned info = ', info
            write(logunit, '(a)') '  T1: FAIL'
            all_pass = .false.
            return
        end if

        max_diff = 0.0_dp
        do i = 1, NN
            max_diff = max(max_diff, abs(x(i) - (1.0_dp, 0.0_dp)))
        end do

        write(logunit, '(a, es12.4)') '  Max deviation from (1,0): ', max_diff

        if (max_diff < 1.0e-12_dp) then
            write(logunit, '(a)') '  T1: PASS'
        else
            write(logunit, '(a)') '  T1: FAIL'
            all_pass = .false.
        end if
    end subroutine test_t1

    ! ==================================================================
    ! T2: 单频正弦峰值检测
    ! x[n] = sin(2π · 8 · n / 64), n=0..63
    ! 预期: 峰值在 k=8 和 k=56，幅度 = 32
    ! ==================================================================
    subroutine test_t2(all_pass, logunit)
        logical, intent(inout) :: all_pass
        integer, intent(in)    :: logunit

        integer,  parameter :: NN = 64
        real(dp), parameter :: FREQ = 8.0_dp
        complex(dp) :: x(NN)
        integer :: info, n, k1, k2
        integer :: i
        real(dp) :: mag(NN), max_mag, second_mag
        integer :: max_idx, second_idx
        integer :: tol_passed

        ! 构造正弦信号
        do n = 0, NN - 1
            x(n + 1) = cmplx(sin(2.0_dp * 4.0_dp * atan(1.0_dp) * FREQ &
                                  * real(n, dp) / real(NN, dp)), &
                              0.0_dp, kind=dp)
        end do

        call fft(x, info)

        if (info /= 0) then
            write(logunit, '(a, i0)') '  ERROR: fft returned info = ', info
            write(logunit, '(a)') '  T2: FAIL'
            all_pass = .false.
            return
        end if

        ! 计算幅度谱
        do i = 1, NN
            mag(i) = abs(x(i))
        end do

        ! 找最大和第二大峰值
        max_idx = 1
        do i = 2, NN
            if (mag(i) > mag(max_idx)) then
                second_idx = max_idx
                max_idx = i
            else if (mag(i) > mag(second_idx) .or. second_idx == max_idx) then
                second_idx = i
            end if
        end do
        ! Ensure second is actually second (not same as max)
        if (max_idx == second_idx) second_idx = 1
        do i = 1, NN
            if (i /= max_idx .and. mag(i) > mag(second_idx)) second_idx = i
        end do

        ! Fortran 1-indexed: FFT bin k (0-based) is at position k+1
        k1 = max_idx - 1   ! convert to 0-based
        k2 = second_idx - 1

        write(logunit, '(a, i0, a, f10.6)') '  Peak 1 at k=', k1, ', |X|=', mag(max_idx)
        write(logunit, '(a, i0, a, f10.6)') '  Peak 2 at k=', k2, ', |X|=', mag(second_idx)

        tol_passed = 0

        ! 检查峰值在 k=8 和 k=56 (N-8)
        if ((k1 == 8 .and. k2 == 56) .or. (k1 == 56 .and. k2 == 8)) then
            tol_passed = tol_passed + 1
            write(logunit, '(a)') '  Peak positions correct.'
        else
            write(logunit, '(a)') '  Peak positions WRONG (expected 8 and 56).'
        end if

        ! 检查幅度 = 32 ± 1e-10
        if (abs(mag(max_idx) - 32.0_dp) < 1.0e-10_dp .and. &
            abs(mag(second_idx) - 32.0_dp) < 1.0e-10_dp) then
            tol_passed = tol_passed + 1
            write(logunit, '(a)') '  Peak magnitudes correct (32.0).'
        else
            write(logunit, '(a)') '  Peak magnitudes WRONG (expected 32.0).'
        end if

        if (tol_passed == 2) then
            write(logunit, '(a)') '  T2: PASS'
        else
            write(logunit, '(a)') '  T2: FAIL'
            all_pass = .false.
        end if
    end subroutine test_t2

    ! ==================================================================
    ! T3: 与暴力 O(N²) DFT 对比
    ! 随机复数输入, N=8
    ! ==================================================================
    subroutine test_t3(all_pass, logunit)
        logical, intent(inout) :: all_pass
        integer, intent(in)    :: logunit

        integer,  parameter :: NN = 8
        complex(dp) :: x(NN), x_save(NN), x_dft(NN)
        integer :: info, k, n
        real(dp) :: angle, max_diff
        complex(dp) :: sum_val

        ! 生成确定性的复数列（非随机，保证可重复）
        do n = 0, NN - 1
            x(n + 1) = cmplx(real(n + 1, dp), real(2 * n, dp), kind=dp)
        end do
        x_save = x

        ! 暴力 O(N²) DFT
        do k = 0, NN - 1
            sum_val = (0.0_dp, 0.0_dp)
            do n = 0, NN - 1
                angle = -2.0_dp * 4.0_dp * atan(1.0_dp) * real(k * n, dp) / real(NN, dp)
                sum_val = sum_val + x_save(n + 1) * cmplx(cos(angle), sin(angle), kind=dp)
            end do
            x_dft(k + 1) = sum_val
        end do

        ! FFT
        call fft(x, info)

        if (info /= 0) then
            write(logunit, '(a, i0)') '  ERROR: fft returned info = ', info
            write(logunit, '(a)') '  T3: FAIL'
            all_pass = .false.
            return
        end if

        ! 对比
        max_diff = 0.0_dp
        do k = 1, NN
            max_diff = max(max_diff, abs(x(k) - x_dft(k)))
        end do

        write(logunit, '(a, es12.4)') '  Max |FFT - DFT|: ', max_diff

        if (max_diff < 1.0e-12_dp) then
            write(logunit, '(a)') '  T3: PASS'
        else
            write(logunit, '(a)') '  T3: FAIL'
            all_pass = .false.
        end if
    end subroutine test_t3

    ! ==================================================================
    ! T4: 正向+逆向往返测试
    ! 输入: x[n] = sin(2πn/16) + i·cos(2πn/16), N=16
    ! 预期: FFT 后 IFFT 恢复原值
    ! ==================================================================
    subroutine test_t4(all_pass, logunit)
        logical, intent(inout) :: all_pass
        integer, intent(in)    :: logunit

        integer,  parameter :: NN = 16
        complex(dp) :: x(NN), original(NN)
        integer :: info, i
        real(dp) :: max_diff

        ! 构造复正弦
        do i = 0, NN - 1
            original(i + 1) = cmplx(sin(2.0_dp * 4.0_dp * atan(1.0_dp) &
                                         * real(i, dp) / real(NN, dp)), &
                                     cos(2.0_dp * 4.0_dp * atan(1.0_dp) &
                                         * real(i, dp) / real(NN, dp)), &
                                     kind=dp)
        end do

        x = original

        ! 正向 FFT
        call fft(x, info)
        if (info /= 0) then
            write(logunit, '(a, i0)') '  ERROR: forward fft returned info = ', info
            write(logunit, '(a)') '  T4: FAIL'
            all_pass = .false.
            return
        end if

        ! 逆变换
        call fft(x, info, inverse=.true.)
        if (info /= 0) then
            write(logunit, '(a, i0)') '  ERROR: inverse fft returned info = ', info
            write(logunit, '(a)') '  T4: FAIL'
            all_pass = .false.
            return
        end if

        ! 对比
        max_diff = 0.0_dp
        do i = 1, NN
            max_diff = max(max_diff, abs(x(i) - original(i)))
        end do

        write(logunit, '(a, es12.4)') '  Max |original - roundtrip|: ', max_diff

        if (max_diff < 1.0e-12_dp) then
            write(logunit, '(a)') '  T4: PASS'
        else
            write(logunit, '(a)') '  T4: FAIL'
            all_pass = .false.
        end if
    end subroutine test_t4

    ! ==================================================================
    ! T5: 非 2 的幂输入检查
    ! N=7, 应返回 info=1 且输入不变
    ! ==================================================================
    subroutine test_t5(all_pass, logunit)
        logical, intent(inout) :: all_pass
        integer, intent(in)    :: logunit

        integer,  parameter :: NN = 7
        complex(dp) :: x(NN), x_save(NN)
        integer :: info, i
        real(dp) :: max_diff

        x = (1.0_dp, 0.0_dp)
        x_save = x

        call fft(x, info)

        write(logunit, '(a, i0)') '  info = ', info

        ! 检查输入未变
        max_diff = 0.0_dp
        do i = 1, NN
            max_diff = max(max_diff, abs(x(i) - x_save(i)))
        end do
        write(logunit, '(a, es12.4)') '  Max change to input: ', max_diff

        if (info == 1 .and. max_diff == 0.0_dp) then
            write(logunit, '(a)') '  T5: PASS'
        else
            write(logunit, '(a)') '  T5: FAIL'
            all_pass = .false.
        end if
    end subroutine test_t5

end program test_fft
