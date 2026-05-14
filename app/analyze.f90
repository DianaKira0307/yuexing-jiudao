!-------------------------------------------------------------------------------
! 月行九道 Phase 5: 全流程分析程序
!
! 读取 Phase 1 输出的二进制数据，做以下处理:
!   1. ICRF → 黄道坐标系旋转 (mod_orbital_elements)
!   2. 笛卡尔 → 黄道球坐标 (r, λ, β)
!   3. λ 去卷绕 + 移除平均运动 → Δλ (经度残差)
!   4. 对 r, Δλ, β 分别计算功率谱并寻峰
!   5. 输出谱峰表 data/peaks_{r, dlambda, beta}.txt
!-------------------------------------------------------------------------------
program lunar_analyze
    use mod_precision, only: dp
    use mod_spectrum, only: spectrum_config, peak_info, &
                            power_spectrum_1d, find_peaks, &
                            WINDOW_HANNING
    use mod_orbital_elements, only: icrf_to_ecliptic, spherical_from_cartesian, &
                                    OBLIQUITY_J2000
    implicit none

    ! ---- 参数 ----
    integer,  parameter :: N_FFT = 65536
    integer,  parameter :: N_PEAKS_MAX = 30
    real(dp), parameter :: FS = 1.0_dp
    character(len=*), parameter :: BIN_PATH = 'data/lunar_state_1900_2100.bin'

    real(dp), parameter :: PI = 3.14159265358979323846_dp
    real(dp), parameter :: PI2 = 2.0_dp * PI

    ! ---- 变量列表 ----
    character(len=8), parameter :: VAR_NAMES(3) = [character(len=8) :: &
        'r', 'dlambda', 'beta']
    character(len=33), parameter :: VAR_TITLES(3) = [character(len=33) :: &
        'Geocentric distance (km)', 'Ecliptic longitude residual (rad)', &
        'Ecliptic latitude (rad)']

    ! ---- 数据 ----
    integer(8) :: file_size
    integer :: n_total, n, funit, i, j, n_peaks, info, idx
    real(dp), allocatable :: data(:, :)
    real(dp), allocatable :: ps(:), freq(:)
    real(dp), allocatable :: r_arr(:), lambda_arr(:), beta_arr(:)
    real(dp), allocatable :: lambda_unwrapped(:), delta_lambda(:)
    real(dp) :: x_e, y_e, z_e, vx_e, vy_e, vz_e
    real(dp) :: r_val, lam, beta, vr_dummy, vlam_dummy, vbeta_dummy
    real(dp) :: offset, delta, sum_t, sum_lam, sum_tt, sum_tl
    real(dp) :: a, b, n_real
    real(dp) :: tdb_start, tdb_end
    type(spectrum_config) :: cfg
    type(peak_info), allocatable :: peaks(:)
    character(len=256) :: outpath

    ! ---- 读取二进制数据 ----
    write(*, '(a)') '=== Lunar Nine Paths — Ecliptic Spectrum Analysis ==='
    write(*, '(a)') ''

    inquire(file=BIN_PATH, size=file_size)
    if (file_size <= 0) then
        write(*, '(a)') 'ERROR: Binary file not found. Run "make run" first.'
        stop 1
    end if

    n_total = int(file_size / (8_8 * 8_8))
    n = min(n_total, N_FFT)
    write(*, '(a, i0, a)') 'Total samples available: ', n_total
    write(*, '(a, i0, a, i0, a)') 'Using first ', n, ' samples for FFT (', &
        nint(real(n, dp) / 365.25_dp), ' years)'

    allocate(data(8, n))
    open(newunit=funit, file=BIN_PATH, access='stream', form='unformatted', status='old')
    read(funit) data(:, 1:n)
    close(funit)

    tdb_start = data(1, 1) + data(2, 1)
    tdb_end   = data(1, n) + data(2, n)
    write(*, '(a, f12.3, a)') 'TDB start: ', tdb_start, ' JD'
    write(*, '(a, f12.3, a)') 'TDB end:   ', tdb_end,   ' JD'
    write(*, '(a, f12.3, a)') 'Span:      ', tdb_end - tdb_start, ' days'
    write(*, '(a)') ''

    ! ---- ICRF → 黄道 → 球坐标 (r, λ, β) ----
    write(*, '(a)') 'Transforming to ecliptic spherical coordinates...'

    allocate(r_arr(n), lambda_arr(n), beta_arr(n))

    do i = 1, n
        call icrf_to_ecliptic(data(3, i), data(4, i), data(5, i), &
                              data(6, i), data(7, i), data(8, i), &
                              OBLIQUITY_J2000, &
                              x_e, y_e, z_e, vx_e, vy_e, vz_e)
        call spherical_from_cartesian(x_e, y_e, z_e, vx_e, vy_e, vz_e, &
                                      r_val, lam, beta, &
                                      vr_dummy, vlam_dummy, vbeta_dummy)
        r_arr(i) = r_val
        lambda_arr(i) = lam
        beta_arr(i) = beta
    end do

    ! ---- λ 去卷绕 ----
    allocate(lambda_unwrapped(n))
    lambda_unwrapped(1) = lambda_arr(1)
    offset = 0.0_dp
    do i = 2, n
        delta = lambda_arr(i) - lambda_arr(i - 1)
        if (delta > PI) then
            offset = offset - PI2
        else if (delta < -PI) then
            offset = offset + PI2
        end if
        lambda_unwrapped(i) = lambda_arr(i) + offset
    end do

    ! ---- 线性拟合 λ = a + b·t，移除平均运动 ----
    sum_t = 0.0_dp; sum_lam = 0.0_dp
    sum_tt = 0.0_dp; sum_tl = 0.0_dp
    n_real = real(n, dp)
    do i = 1, n
        sum_t = sum_t + real(i - 1, dp)
        sum_lam = sum_lam + lambda_unwrapped(i)
        sum_tt = sum_tt + real(i - 1, dp)**2
        sum_tl = sum_tl + real(i - 1, dp) * lambda_unwrapped(i)
    end do
    b = (n_real * sum_tl - sum_t * sum_lam) / (n_real * sum_tt - sum_t * sum_t)
    a = (sum_lam - b * sum_t) / n_real

    allocate(delta_lambda(n))
    do i = 1, n
        delta_lambda(i) = lambda_unwrapped(i) - (a + b * real(i - 1, dp))
    end do

    write(*, '(a, f12.8, a)') 'Mean motion: ', b, ' rad/day'
    write(*, '(a, f12.8, a)') '  = ', b * 180.0_dp / PI, ' deg/day'
    write(*, '(a, f10.4, a)') '  = ', PI2 / b, ' days (sidereal period)'
    write(*, '(a)') ''

    ! ---- 配置频谱分析 ----
    cfg%n_fft = n
    cfg%window_type = WINDOW_HANNING
    cfg%fs = FS
    cfg%n_peaks_max = N_PEAKS_MAX
    cfg%peak_threshold_rel = 1.0e-4_dp
    cfg%detrend = .true.

    ! ---- 逐变量分析 ----
    do idx = 1, 3
        write(*, '(a, i0, a, a)') '--- Variable ', idx, ': ', trim(VAR_NAMES(idx))

        select case (idx)
        case (1)
            call power_spectrum_1d(r_arr, cfg, ps, freq, info)
        case (2)
            call power_spectrum_1d(delta_lambda, cfg, ps, freq, info)
        case (3)
            call power_spectrum_1d(beta_arr, cfg, ps, freq, info)
        end select

        if (info /= 0) then
            write(*, '(a, i0)') '  ERROR: power_spectrum_1d returned ', info
            cycle
        end if

        call find_peaks(ps, freq, n / 2 + 1, cfg, peaks, n_peaks, info)
        if (info /= 0) then
            write(*, '(a, i0)') '  ERROR: find_peaks returned ', info
            cycle
        end if

        write(*, '(a, i0, a)') '  Found ', n_peaks, ' peaks'

        outpath = 'data/peaks_' // trim(VAR_NAMES(idx)) // '.txt'
        open(newunit=funit, file=trim(outpath), status='replace', action='write')
        write(funit, '(a)') '# Lunar Nine Paths — Spectrum peaks for ' // trim(VAR_TITLES(idx))
        write(funit, '(a, i0)') '# FFT length: ', n
        write(funit, '(a, f8.4, a)') '# Sampling rate: ', FS, ' cycle/day'
        write(funit, '(a)') '#'
        write(funit, '(a)') '# rank  index   freq(cyc/day)  period(days)  power'
        write(funit, '(a)') '# ---- -----  -------------  ------------  ------'

        do j = 1, n_peaks
            write(funit, '(i6, i6, 2x, es14.6, 2x, f12.2, 2x, es14.6)') &
                j, peaks(j)%index, peaks(j)%freq, peaks(j)%period, peaks(j)%power
        end do
        close(funit)

        write(*, '(a, a)') '  Output: ', trim(outpath)
    end do

    ! ---- 综合输出 ----
    write(*, '(a)') ''
    write(*, '(a)') 'All peak tables written to data/peaks_{r, dlambda, beta}.txt'
    write(*, '(a)') 'Run: python3 plot_spectra_mpl.py  (for visualization)'

end program lunar_analyze
