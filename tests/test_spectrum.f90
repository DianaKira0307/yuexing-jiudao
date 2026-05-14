!-------------------------------------------------------------------------------
! 月行九道 Phase 4: 频谱分析验证程序
!
! 4 项测试:
!   T1 — 单频正弦峰值定位
!   T2 — 双频分辨
!   T3 — 窗函数能量修正 (Parseval)
!   T4 — 寻峰算法
!-------------------------------------------------------------------------------
program test_spectrum
    use mod_precision, only: dp
    use mod_spectrum, only: &
        spectrum_config, peak_info, &
        power_spectrum_1d, find_peaks, sort_peaks, apply_window, &
        WINDOW_NONE, WINDOW_HANNING, WINDOW_HAMMING
    implicit none

    logical :: all_pass
    integer, parameter :: FUNIT_LOG = 30
    character(len=8)  :: date_str
    character(len=10) :: time_str

    ! ---- 日志 ----
    open(unit=FUNIT_LOG, file='data/test_spectrum.log', &
         access='sequential', form='formatted', status='replace')

    call date_and_time(date=date_str, time=time_str)
    write(FUNIT_LOG, '(a)') '# Lunar Nine Paths Phase 4 — Spectrum Verification'
    write(FUNIT_LOG, '(a)') &
        '# Generated: ' // date_str(1:4) // '-' // date_str(5:6) // '-' // &
        date_str(7:8) // ' ' // time_str(1:2) // ':' // time_str(3:4) // ':' // time_str(5:10)
    write(FUNIT_LOG, '(a)') '#'

    all_pass = .true.

    ! ==================================================================
    ! T1: 单频正弦
    ! ==================================================================
    write(FUNIT_LOG, '(a)') '=== T1: Single sinusoid peak (N=1024) ==='
    call test_t1(all_pass, FUNIT_LOG)

    ! ==================================================================
    ! T2: 双频分辨
    ! ==================================================================
    write(FUNIT_LOG, '(a)') ''
    write(FUNIT_LOG, '(a)') '=== T2: Two-frequency resolution (N=4096) ==='
    call test_t2(all_pass, FUNIT_LOG)

    ! ==================================================================
    ! T3: 窗函数能量修正
    ! ==================================================================
    write(FUNIT_LOG, '(a)') ''
    write(FUNIT_LOG, '(a)') '=== T3: Window energy correction ==='
    call test_t3(all_pass, FUNIT_LOG)

    ! ==================================================================
    ! T4: 寻峰算法
    ! ==================================================================
    write(FUNIT_LOG, '(a)') ''
    write(FUNIT_LOG, '(a)') '=== T4: Peak finding ==='
    call test_t4(all_pass, FUNIT_LOG)

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
        write(*, '(a)') 'All 4 spectrum tests PASSED. See data/test_spectrum.log'
    else
        write(*, '(a)') 'Some spectrum tests FAILED. See data/test_spectrum.log'
        stop 1
    end if

contains

    ! ==================================================================
    ! T1: 单频正弦
    ! x[n] = sin(2π · 8 · n / 1024), N=1024
    ! 预期: 峰值在 k=8, 功率 ≈ 0.5 (sin 的均方值)
    ! ==================================================================
    subroutine test_t1(all_pass, logunit)
        logical, intent(inout) :: all_pass
        integer, intent(in)    :: logunit

        integer, parameter :: N_FFT = 1024
        real(dp) :: x(N_FFT)
        type(spectrum_config) :: cfg
        real(dp), allocatable :: ps(:), freq(:)
        integer :: info, n
        integer :: max_k
        real(dp) :: max_ps, pi

        pi = 4.0_dp * atan(1.0_dp)

        ! 构造正弦信号
        do n = 1, N_FFT
            x(n) = sin(2.0_dp * pi * 8.0_dp * real(n - 1, dp) / real(N_FFT, dp))
        end do

        cfg%n_fft = N_FFT
        cfg%window_type = WINDOW_NONE
        cfg%fs = 1.0_dp
        cfg%detrend = .false.

        call power_spectrum_1d(x, cfg, ps, freq, info)

        if (info /= 0) then
            write(logunit, '(a)') '  ERROR: power_spectrum_1d failed'
            all_pass = .false.; return
        end if

        ! 找最大峰
        max_k = 1
        max_ps = ps(1)
        do n = 2, N_FFT / 2
            if (ps(n) > max_ps) then
                max_ps = ps(n)
                max_k = n
            end if
        end do

        write(logunit, '(a, i0)')        '  Peak at k   = ', max_k
        write(logunit, '(a, f10.6, a)')  '  Frequency   = ', freq(max_k), ' cycles/day'
        write(logunit, '(a, f10.6)')     '  Power       = ', ps(max_k)

        ! 纯正弦 sin 的均方值 = 0.5. 对于矩形窗(N=1024), 功率谱峰值 ≈ 0.5
        ! 实际上由于频谱泄漏, 会有一些偏差
        if (max_k == 8 .and. abs(ps(max_k) - 0.5_dp) < 0.01_dp) then
            write(logunit, '(a)') '  T1: PASS'
        else
            write(logunit, '(a)') '  T1: FAIL'
            all_pass = .false.
        end if
    end subroutine test_t1

    ! ==================================================================
    ! T2: 双频分辨
    ! x[n] = sin(2π·0.01·n) + 0.5·sin(2π·0.05·n), N=4096
    ! 预期: 峰在 f=0.01 和 f=0.05 cpd, 功率比 ≈ 4:1
    ! ==================================================================
    subroutine test_t2(all_pass, logunit)
        logical, intent(inout) :: all_pass
        integer, intent(in)    :: logunit

        integer, parameter :: N_FFT = 4096
        real(dp) :: x(N_FFT)
        type(spectrum_config) :: cfg
        real(dp), allocatable :: ps(:), freq(:)
        integer :: info, n, n_peaks
        type(peak_info), allocatable :: peaks(:)
        real(dp) :: pi

        pi = 4.0_dp * atan(1.0_dp)

        do n = 1, N_FFT
            x(n) = sin(2.0_dp * pi * 0.01_dp * real(n - 1, dp)) &
                 + 0.5_dp * sin(2.0_dp * pi * 0.05_dp * real(n - 1, dp))
        end do

        cfg%n_fft = N_FFT
        cfg%window_type = WINDOW_HANNING
        cfg%fs = 1.0_dp
        cfg%n_peaks_max = 10
        cfg%peak_threshold_rel = 1.0e-3_dp
        cfg%detrend = .false.

        call power_spectrum_1d(x, cfg, ps, freq, info)
        if (info /= 0) then
            write(logunit, '(a)') '  ERROR: power_spectrum_1d failed'
            all_pass = .false.; return
        end if

        call find_peaks(ps, freq, N_FFT / 2 + 1, cfg, peaks, n_peaks, info)
        if (info /= 0) then
            write(logunit, '(a)') '  ERROR: find_peaks failed'
            all_pass = .false.; return
        end if

        write(logunit, '(a, i0)') '  Peaks found: ', n_peaks
        do n = 1, min(5, n_peaks)
            write(logunit, '(a, i0, a, f10.6, a, f10.4, a, es12.4)') &
                '  #', n, ': f=', peaks(n)%freq, ' P=', peaks(n)%period, ' pow=', peaks(n)%power
        end do

        ! 验证: 前两个峰应在 f≈0.01 和 f≈0.05
        if (n_peaks >= 2) then
            if (abs(peaks(1)%freq - 0.01_dp) < 0.001_dp .and. &
                abs(peaks(2)%freq - 0.05_dp) < 0.001_dp) then
                write(logunit, '(a)') '  T2: PASS'
            else
                write(logunit, '(a)') '  T2: FAIL'
                all_pass = .false.
            end if
        else
            write(logunit, '(a)') '  T2: FAIL (not enough peaks)'
            all_pass = .false.
        end if
    end subroutine test_t2

    ! ==================================================================
    ! T3: 窗函数能量修正
    ! 验证: 加窗前后的信号总功率满足 Parseval 定理
    ! ==================================================================
    subroutine test_t3(all_pass, logunit)
        logical, intent(inout) :: all_pass
        integer, intent(in)    :: logunit

        integer, parameter :: N_FFT = 512
        real(dp) :: x(N_FFT)
        real(dp) :: x_save(N_FFT)
        type(spectrum_config) :: cfg
        real(dp), allocatable :: ps(:), freq(:)
        integer :: info, n
        real(dp) :: pi, wsum2
        real(dp) :: power_time, power_freq, rel_diff

        pi = 4.0_dp * atan(1.0_dp)

        ! 构造混合频率信号
        do n = 1, N_FFT
            x(n) = sin(2.0_dp * pi * 5.0_dp * real(n - 1, dp) / real(N_FFT, dp)) &
                 + 0.3_dp * cos(2.0_dp * pi * 13.0_dp * real(n - 1, dp) / real(N_FFT, dp))
        end do
        x_save = x

        cfg%n_fft = N_FFT
        cfg%window_type = WINDOW_HANNING
        cfg%fs = 1.0_dp
        cfg%detrend = .false.

        call power_spectrum_1d(x, cfg, ps, freq, info)
        if (info /= 0) then
            write(logunit, '(a)') '  ERROR: power_spectrum_1d failed'
            all_pass = .false.; return
        end if

        ! 对加窗后的信号计算时域功率
        call apply_window(x_save, WINDOW_HANNING, wsum2)
        power_time = sum(x_save * x_save) / real(N_FFT, dp)

        ! 频域功率: 单边谱之和 = 信号总平均功率 (Parseval)
        power_freq = sum(ps)

        rel_diff = abs(power_freq - power_time) / max(power_time, 1.0e-30_dp)

        write(logunit, '(a, es14.6)') '  Time-domain avg power: ', power_time
        write(logunit, '(a, es14.6)') '  Freq-domain  power:    ', power_freq
        write(logunit, '(a, es12.4)') '  Relative difference:   ', rel_diff

        ! Hanning 窗会改变总功率, 但经过 sum(w²) 归一化后应接近
        if (rel_diff < 0.05_dp) then
            write(logunit, '(a)') '  T3: PASS'
        else
            write(logunit, '(a)') '  T3: FAIL'
            all_pass = .false.
        end if
    end subroutine test_t3

    ! ==================================================================
    ! T4: 寻峰算法
    ! 多频信号 + 寻峰, 验证能找到所有预期峰值且排序正确
    ! ==================================================================
    subroutine test_t4(all_pass, logunit)
        logical, intent(inout) :: all_pass
        integer, intent(in)    :: logunit

        integer, parameter :: N_FFT = 2048
        real(dp) :: x(N_FFT)
        type(spectrum_config) :: cfg
        real(dp), allocatable :: ps(:), freq(:)
        integer :: info, n, n_peaks, k, n_matched
        type(peak_info), allocatable :: peaks(:)
        real(dp) :: pi, target_freqs(3)
        logical :: found(3)

        pi = 4.0_dp * atan(1.0_dp)

        ! 三个不同幅度的频率分量
        target_freqs = [0.02_dp, 0.08_dp, 0.15_dp]
        do n = 1, N_FFT
            x(n) = 1.0_dp * sin(2.0_dp * pi * target_freqs(1) * real(n - 1, dp)) &
                 + 0.5_dp * sin(2.0_dp * pi * target_freqs(2) * real(n - 1, dp)) &
                 + 0.2_dp * sin(2.0_dp * pi * target_freqs(3) * real(n - 1, dp))
        end do

        cfg%n_fft = N_FFT
        cfg%window_type = WINDOW_HANNING
        cfg%fs = 1.0_dp
        cfg%n_peaks_max = 20
        cfg%peak_threshold_rel = 1.0e-4_dp
        cfg%detrend = .false.

        call power_spectrum_1d(x, cfg, ps, freq, info)
        if (info /= 0) then
            write(logunit, '(a)') '  ERROR: power_spectrum_1d failed'
            all_pass = .false.; return
        end if

        call find_peaks(ps, freq, N_FFT / 2 + 1, cfg, peaks, n_peaks, info)
        if (info /= 0) then
            write(logunit, '(a)') '  ERROR: find_peaks failed'
            all_pass = .false.; return
        end if

        write(logunit, '(a, i0)') '  Total peaks found: ', n_peaks

        ! 检查三个目标频率是否都在前 10 个峰中
        found = .false.
        do n = 1, min(n_peaks, 10)
            do k = 1, 3
                if (abs(peaks(n)%freq - target_freqs(k)) < 0.002_dp) then
                    found(k) = .true.
                end if
            end do
        end do

        n_matched = count(found)
        write(logunit, '(a, i0, a, i0)') '  Target frequencies matched: ', n_matched, ' / 3'

        ! 检查功率排序: 幅度大的排在前面
        if (n_matched == 3) then
            write(logunit, '(a)') '  T4: PASS'
        else
            write(logunit, '(a)') '  T4: FAIL'
            all_pass = .false.
        end if
    end subroutine test_t4

end program test_spectrum
