!-------------------------------------------------------------------------------
! 月行九道 Phase 4: 频谱分析模块
! 模块: mod_spectrum
!
! 功能: 加窗、功率谱、找峰、峰值排序
!-------------------------------------------------------------------------------

module mod_spectrum
    use mod_precision, only: dp
    use mod_fft, only: fft
    implicit none
    private

    public :: apply_window, power_spectrum_1d, find_peaks, sort_peaks
    public :: spectrum_config, peak_info

    ! 窗类型常量
    integer, parameter, public :: WINDOW_NONE    = 0
    integer, parameter, public :: WINDOW_HANNING = 1
    integer, parameter, public :: WINDOW_HAMMING = 2

    real(dp), parameter :: PI = 3.14159265358979323846_dp
    real(dp), parameter :: PI2 = 2.0_dp * PI

    !-------------------------------------------------------------------
    ! 频谱分析配置
    !-------------------------------------------------------------------
    type :: spectrum_config
        integer  :: n_fft = 65536       ! FFT 长度 (2 的幂)
        integer  :: window_type = 1     ! 0=none, 1=hanning, 2=hamming
        real(dp) :: fs = 1.0_dp         ! 采样率 (cycles/day)
        integer  :: n_peaks_max = 50    ! 输出最大峰数
        real(dp) :: peak_threshold_rel  = 1.0e-3_dp  ! 相对阈值 (相对最大峰)
        logical  :: detrend = .true.    ! 去均值
    end type

    !-------------------------------------------------------------------
    ! 峰值信息
    !-------------------------------------------------------------------
    type :: peak_info
        integer  :: index               ! FFT bin 索引 (0-based)
        real(dp) :: freq                ! 频率 (cycles/day)
        real(dp) :: period              ! 周期 (days)
        real(dp) :: power               ! 功率谱值
    end type

contains

    !-------------------------------------------------------------------
    ! apply_window: 对信号 x 加窗
    !
    ! 参数:
    !   x           (inout) — 信号序列
    !   window_type (in)    — 窗类型 (WINDOW_NONE/HANNING/HAMMING)
    !
    ! 输出:
    !   wsum2 — 窗能量和 sum(w²), 用于后续归一化
    !-------------------------------------------------------------------
    subroutine apply_window(x, window_type, wsum2)
        real(dp), intent(inout) :: x(:)
        integer,  intent(in)    :: window_type
        real(dp), intent(out), optional :: wsum2

        integer  :: n, i
        real(dp) :: w, s2

        n = size(x)
        s2 = 0.0_dp

        select case (window_type)
        case (WINDOW_NONE)
            do i = 1, n
                w = 1.0_dp
                x(i) = x(i) * w
                s2 = s2 + w * w
            end do

        case (WINDOW_HANNING)
            do i = 1, n
                w = 0.5_dp * (1.0_dp - cos(PI2 * real(i - 1, dp) / real(n, dp)))
                x(i) = x(i) * w
                s2 = s2 + w * w
            end do

        case (WINDOW_HAMMING)
            do i = 1, n
                w = 0.54_dp - 0.46_dp * cos(PI2 * real(i - 1, dp) / real(n, dp))
                x(i) = x(i) * w
                s2 = s2 + w * w
            end do

        case default
            do i = 1, n
                w = 1.0_dp
                x(i) = x(i) * w
                s2 = s2 + w * w
            end do
        end select

        if (present(wsum2)) wsum2 = s2
    end subroutine apply_window

    !-------------------------------------------------------------------
    ! power_spectrum_1d: 计算单边功率谱
    !
    ! 1. 截断/补零到 FFT 长度 (n_fft)
    ! 2. 去均值
    ! 3. 加窗
    ! 4. FFT
    ! 5. 单边功率谱 P[k] = 2·|X[k]|² / sum(w²) (k=1..N/2-1)
    !    P[0] = |X[0]|² / sum(w²), P[N/2] = |X[N/2]|² / sum(w²)
    !
    ! 输出:
    !   ps(0:n/2)   — 功率谱值
    !   freq(0:n/2) — 对应频率 (cycles/day)
    !   info        — 返回码
    !-------------------------------------------------------------------
    subroutine power_spectrum_1d(x, config, ps, freq, info)
        real(dp),            intent(in)  :: x(:)
        type(spectrum_config), intent(in)  :: config
        real(dp), allocatable, intent(out) :: ps(:), freq(:)
        integer,             intent(out) :: info

        integer  :: n_signal, n_fft, n_half, i
        real(dp), allocatable :: buf(:)
        complex(dp), allocatable :: cbuf(:)
        real(dp) :: mean
        integer  :: fft_info

        info = 0
        n_signal = size(x)
        n_fft    = config%n_fft
        n_half   = n_fft / 2

        ! ---- 1. 准备 FFT 缓冲区 (截断或补零) ----
        allocate(buf(n_fft))
        buf = 0.0_dp

        if (n_signal <= n_fft) then
            buf(1:n_signal) = x(1:n_signal)
        else
            buf(1:n_fft) = x(1:n_fft)
        end if

        ! ---- 2. 去均值 ----
        if (config%detrend) then
            mean = sum(buf) / real(n_fft, dp)
            buf = buf - mean
        end if

        ! ---- 3. 加窗 ----
        call apply_window(buf, config%window_type)

        ! ---- 4. FFT (打包为复数) ----
        allocate(cbuf(n_fft))
        cbuf = cmplx(buf, 0.0_dp, kind=dp)
        call fft(cbuf, fft_info)
        if (fft_info /= 0) then
            info = fft_info
            return
        end if

        ! ---- 5. 单边功率谱 ----
        allocate(ps(0:n_half))
        allocate(freq(0:n_half))

        do i = 0, n_half
            freq(i) = real(i, dp) * config%fs / real(n_fft, dp)
        end do

        ! 归一化: P[k] = |X[k]|² / N²  (加窗后信号的功率谱)
        ! Parseval 验证: sum(P_two_sided) = sum(xw²) / N = 平均功率
        ! 其中 xw = x · w 为加窗后的信号
        do i = 0, n_half
            ps(i) = abs(cbuf(i + 1)) ** 2 / (real(n_fft, dp) ** 2)
        end do
        ! 单边谱: 正频率加倍 (DC 和 Nyquist 除外)
        do i = 1, n_half - 1
            ps(i) = 2.0_dp * ps(i)
        end do

    end subroutine power_spectrum_1d

    !-------------------------------------------------------------------
    ! find_peaks: 在功率谱中寻找局部峰值
    !
    ! 从 k=1 到 n_freq-2 扫描, 局部极大值需大于两侧相邻 bin.
    ! 过滤: peaks(:,power) > max(power) * config%peak_threshold_rel
    !
    ! 输出:
    !   peaks  — 峰值数组, 已按功率降序排列
    !   n_found — 实际找到的峰数
    !-------------------------------------------------------------------
    subroutine find_peaks(ps, freq, n_freq, config, peaks, n_found, info)
        real(dp),            intent(in)  :: ps(0:), freq(0:)
        integer,             intent(in)  :: n_freq
        type(spectrum_config), intent(in)  :: config
        type(peak_info), allocatable, intent(out) :: peaks(:)
        integer,             intent(out) :: n_found
        integer,             intent(out) :: info

        integer  :: k, n_candidates, n_max
        real(dp) :: threshold
        type(peak_info), allocatable :: candidates(:)

        info = 0
        n_max = config%n_peaks_max
        if (n_max < 1) n_max = 1

        ! 第一阶段: 收集所有局部峰值
        allocate(candidates(n_freq / 2 + 1))
        n_candidates = 0

        ! DC (k=0) 和 Nyquist (k=n_freq-1) 不计入
        do k = 1, n_freq - 2
            if (ps(k) > ps(k - 1) .and. ps(k) > ps(k + 1)) then
                n_candidates = n_candidates + 1
                candidates(n_candidates)%index  = k
                candidates(n_candidates)%freq   = freq(k)
                if (freq(k) > 0.0_dp) then
                    candidates(n_candidates)%period = 1.0_dp / freq(k)
                else
                    candidates(n_candidates)%period = huge(1.0_dp)
                end if
                candidates(n_candidates)%power  = ps(k)
            end if
        end do

        if (n_candidates == 0) then
            allocate(peaks(0))
            n_found = 0
            return
        end if

        ! 过滤: 相对阈值
        threshold = -huge(1.0_dp)
        do k = 1, n_candidates
            if (candidates(k)%power > threshold) threshold = candidates(k)%power
        end do
        threshold = threshold * config%peak_threshold_rel

        ! 排序 (按功率降序)
        call sort_peaks(candidates, n_candidates)

        ! 截取前 n_max 个
        n_found = min(n_candidates, n_max)
        allocate(peaks(n_found))
        peaks(1:n_found) = candidates(1:n_found)

    end subroutine find_peaks

    !-------------------------------------------------------------------
    ! sort_peaks: 将峰值数组按功率降序排列 (插入排序)
    !
    ! 参数:
    !   peaks — 待排序数组
    !   n     — 有效元素数
    !-------------------------------------------------------------------
    subroutine sort_peaks(peaks, n)
        type(peak_info), intent(inout) :: peaks(:)
        integer,         intent(in)    :: n

        integer        :: i, j
        type(peak_info) :: key

        do i = 2, n
            key = peaks(i)
            j = i - 1
            do while (j >= 1 .and. peaks(j)%power < key%power)
                peaks(j + 1) = peaks(j)
                j = j - 1
            end do
            peaks(j + 1) = key
        end do
    end subroutine sort_peaks

end module mod_spectrum
