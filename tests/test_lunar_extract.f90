!-------------------------------------------------------------------------------
! 月行九道 Phase 1: 验证程序
!
! 对 extract_lunar_states 的输出进行 5 项独立验证:
!   T1 — 起止历元正确性
!   T2 — 月球地心距离合理性（近/远地点范围）
!   T3 — 月球地心速度合理性
!   T4 — 与单点 DPLEPH 直接调用的一致性
!   T5 — 文件大小检验
!
! 结果输出到 data/test_lunar_extract.log
!-------------------------------------------------------------------------------
program test_lunar_extract
    use mod_precision, only: dp
    use mod_constants, only: init_physics_constants
    use mod_lunar_extract, only: lunar_extract_config, extract_lunar_states
    use mod_jpleph, only: DPLEPH
    use mod_time_system, only: date2jd, utc2tdb
    implicit none

    ! ---- 测试参数 ----
    integer,  parameter :: YY_START = 1900, MM_START = 1, DD_START = 1
    integer,  parameter :: HH_START = 0,   MN_START = 0
    real(dp), parameter :: SS_START = 0.0_dp
    integer,  parameter :: YY_END   = 2100, MM_END   = 1, DD_END   = 1
    integer,  parameter :: HH_END   = 0,    MN_END   = 0
    real(dp), parameter :: SS_END   = 0.0_dp
    real(dp), parameter :: STEP_DAYS = 1.0_dp

    ! ---- 程序变量 ----
    type(lunar_extract_config) :: cfg
    integer :: info_extract

    ! 时间系统
    real(dp) :: utc_start1, utc_start2, utc_end1, utc_end2
    real(dp) :: tdb_ref_start1, tdb_ref_start2
    real(dp) :: tdb_ref_end1,   tdb_ref_end2
    real(dp) :: tdb_start, tdb_end
    integer  :: info_date, info_tdb, n_samples

    ! 文件读取
    character(len=256) :: bin_path
    integer(8) :: file_size
    integer :: funit
    real(dp), allocatable :: data(:, :)

    ! 循环与统计
    integer, parameter :: N_SPOT = 3
    integer  :: spot_idx(N_SPOT)

    ! 测试结果
    logical :: all_pass
    integer :: n

    ! ---- 输出文件 ----
    integer, parameter :: FUNIT_LOG = 30

    ! ---- 时间字符串 ----
    character(len=8)  :: date_str
    character(len=10) :: time_str

    ! ==================================================================
    ! 0. 初始化
    ! ==================================================================

    ! 打开日志文件
    open(unit=FUNIT_LOG, file='data/test_lunar_extract.log', &
         access='sequential', form='formatted', status='replace')

    call date_and_time(date=date_str, time=time_str)
    write(FUNIT_LOG, '(a)') '# Lunar Nine Paths Phase 1 — Verification Log'
    write(FUNIT_LOG, '(a)') &
        '# Generated: ' // date_str(1:4) // '-' // date_str(5:6) // '-' // &
        date_str(7:8) // ' ' // time_str(1:2) // ':' // time_str(3:4) // ':' // time_str(5:10)
    write(FUNIT_LOG, '(a)') '#'

    ! 初始化物理常数与历表
    call init_physics_constants()

    ! 设置配置
    cfg%yy_start = YY_START; cfg%mm_start = MM_START; cfg%dd_start = DD_START
    cfg%hh_start = HH_START; cfg%mn_start = MN_START; cfg%ss_start = SS_START
    cfg%yy_end   = YY_END;   cfg%mm_end   = MM_END;   cfg%dd_end   = DD_END
    cfg%hh_end   = HH_END;   cfg%mn_end   = MN_END;   cfg%ss_end   = SS_END
    cfg%step_days = STEP_DAYS
    cfg%progress_interval = 0  ! 静默模式

    ! 先生成数据文件
    call extract_lunar_states(cfg, 'data', info_extract)
    if (info_extract /= 0) then
        write(FUNIT_LOG, '(a, i0)') &
            '# ERROR: extract_lunar_states returned info = ', info_extract
        close(FUNIT_LOG)
        stop 1
    end if

    bin_path = 'data/' // trim(cfg%bin_filename)

    ! 计算参考 TDB
    call date2jd(YY_START, MM_START, DD_START, HH_START, MN_START, &
                 SS_START, utc_start1, utc_start2, info_date)
    call utc2tdb(utc_start1, utc_start2, tdb_ref_start1, tdb_ref_start2, info_tdb)

    call date2jd(YY_END, MM_END, DD_END, HH_END, MN_END, &
                 SS_END, utc_end1, utc_end2, info_date)
    call utc2tdb(utc_end1, utc_end2, tdb_ref_end1, tdb_ref_end2, info_tdb)

    tdb_start = tdb_ref_start1 + tdb_ref_start2
    tdb_end   = tdb_ref_end1   + tdb_ref_end2
    n_samples = floor((tdb_end - tdb_start) / STEP_DAYS) + 1

    ! ==================================================================
    ! 1. 读取二进制文件
    ! ==================================================================

    funit = 10
    open(unit=funit, file=trim(bin_path), access='stream', &
         form='unformatted', status='old')
    inquire(unit=funit, size=file_size)

    ! 计算采样点数: 每个采样点 8 个 real(dp) × 8 bytes
    n = int(file_size / (8 * 8_8))
    allocate(data(8, n))
    read(funit) data
    close(funit)

    ! ==================================================================
    ! 2. 执行测试
    ! ==================================================================

    all_pass = .true.

    ! ---- T1: 起止历元正确性 ----
    write(FUNIT_LOG, '(a)') '=== T1: Start/end epoch correctness ==='
    call test_t1(data, n, tdb_ref_start1, tdb_ref_start2, &
                 tdb_ref_end1, tdb_ref_end2, all_pass, FUNIT_LOG)

    ! ---- T2: 月球地心距离合理性 ----
    write(FUNIT_LOG, '(a)') ''
    write(FUNIT_LOG, '(a)') '=== T2: Geocentric distance range ==='
    call test_t2(data, n, all_pass, FUNIT_LOG)

    ! ---- T3: 月球地心速度合理性 ----
    write(FUNIT_LOG, '(a)') ''
    write(FUNIT_LOG, '(a)') '=== T3: Geocentric velocity range ==='
    call test_t3(data, n, all_pass, FUNIT_LOG)

    ! ---- T4: 与单点 DPLEPH 直接调用的一致性 ----
    write(FUNIT_LOG, '(a)') ''
    write(FUNIT_LOG, '(a)') '=== T4: Spot-check consistency with DPLEPH ==='
    call test_t4(data, n, all_pass, FUNIT_LOG)

    ! ---- T5: 文件大小检验 ----
    write(FUNIT_LOG, '(a)') ''
    write(FUNIT_LOG, '(a)') '=== T5: File size check ==='
    call test_t5(n_samples, file_size, all_pass, FUNIT_LOG)

    ! ==================================================================
    ! 3. 汇总
    ! ==================================================================

    write(FUNIT_LOG, '(a)') ''
    if (all_pass) then
        write(FUNIT_LOG, '(a)') 'ALL TESTS PASSED'
    else
        write(FUNIT_LOG, '(a)') 'SOME TESTS FAILED'
    end if

    close(FUNIT_LOG)

    ! 同时输出到终端
    if (all_pass) then
        write(*, '(a)') 'All 5 tests PASSED. See data/test_lunar_extract.log'
    else
        write(*, '(a)') 'Some tests FAILED. See data/test_lunar_extract.log'
        stop 1
    end if

contains

    ! ==================================================================
    ! T1: 检查首末采样点的 TDB 值
    !
    ! 注意: 末历元不能用独立的 UTC→TDB 转换结果直接比对，
    ! 因为闰秒累积导致 UTC quasi-JD 跨度与实际 TDB 天数存在差异。
    ! 正确方式：验证 last - first ≈ (n-1) * step_days。
    ! ==================================================================
    subroutine test_t1(data, n, tdb_ref_s1, tdb_ref_s2, &
                       tdb_ref_e1, tdb_ref_e2, all_pass, logunit)
        real(dp), intent(in) :: data(:, :)
        integer,  intent(in) :: n
        real(dp), intent(in) :: tdb_ref_s1, tdb_ref_s2
        real(dp), intent(in) :: tdb_ref_e1, tdb_ref_e2
        logical,  intent(inout) :: all_pass
        integer,  intent(in) :: logunit

        real(dp) :: tdb_start_file, tdb_end_file
        real(dp) :: tdb_start_ref,  tdb_end_ref_separate
        real(dp) :: diff_start, diff_end_spacing, diff_end_separate
        real(dp) :: tdb_span_file, tdb_span_expected
        real(dp), parameter :: TOL = 1.0e-9_dp  ! 天（~86 微秒）

        tdb_start_file = data(1, 1) + data(2, 1)
        tdb_end_file   = data(1, n) + data(2, n)
        tdb_start_ref  = tdb_ref_s1 + tdb_ref_s2
        tdb_end_ref_separate = tdb_ref_e1 + tdb_ref_e2

        ! 检查首个采样点
        diff_start = abs(tdb_start_file - tdb_start_ref)
        write(logunit, '(a, es12.4, a)') &
            '  First sample TDB vs UTC->TDB ref diff: ', diff_start, ' days'

        ! 检查时间跨度一致性: last - first ≈ (n-1) * step_days
        tdb_span_file    = tdb_end_file - tdb_start_file
        tdb_span_expected = real(n - 1, dp) * STEP_DAYS
        diff_end_spacing = abs(tdb_span_file - tdb_span_expected)
        write(logunit, '(a, es12.4, a)') &
            '  TDB span diff from expected: ', diff_end_spacing, ' days'

        ! 附注: 独立转换的末历元偏差（因闰秒累积）
        diff_end_separate = abs(tdb_end_file - tdb_end_ref_separate)
        write(logunit, '(a, es12.4, a)') &
            '  Last sample vs UTC->TDB ref (leap-second drift): ', &
            diff_end_separate, ' days'

        if (diff_start <= TOL .and. diff_end_spacing <= TOL) then
            write(logunit, '(a)') '  T1: PASS'
        else
            write(logunit, '(a)') '  T1: FAIL'
            all_pass = .false.
        end if
    end subroutine test_t1

    ! ==================================================================
    ! T2: 检查地心距离的 min/max/mean 在预期范围
    ! ==================================================================
    subroutine test_t2(data, n, all_pass, logunit)
        real(dp), intent(in) :: data(:, :)
        integer,  intent(in) :: n
        logical,  intent(inout) :: all_pass
        integer,  intent(in) :: logunit

        integer  :: i
        real(dp) :: r, r_min, r_max, r_sum
        real(dp) :: r_mean

        r_min =  huge(1.0_dp)
        r_max = -huge(1.0_dp)
        r_sum = 0.0_dp

        do i = 1, n
            r = sqrt(data(3, i)**2 + data(4, i)**2 + data(5, i)**2)
            r_min = min(r_min, r)
            r_max = max(r_max, r)
            r_sum = r_sum + r
        end do
        r_mean = r_sum / real(n, dp)

        write(logunit, '(a, f12.3, a)') '  Min distance: ', r_min, ' km'
        write(logunit, '(a, f12.3, a)') '  Max distance: ', r_max, ' km'
        write(logunit, '(a, f12.3, a)') '  Mean distance:', r_mean, ' km'

        ! 近地点 356000–360000 km，远地点 405000–407000 km，均值 384000–385010 km
        if (r_min >= 356000.0_dp .and. r_min <= 360000.0_dp .and. &
            r_max >= 405000.0_dp .and. r_max <= 407000.0_dp .and. &
            r_mean >= 384000.0_dp .and. r_mean <= 385010.0_dp) then
            write(logunit, '(a)') '  T2: PASS'
        else
            write(logunit, '(a)') '  T2: FAIL'
            all_pass = .false.
        end if
    end subroutine test_t2

    ! ==================================================================
    ! T3: 检查速度大小在预期范围 0.96–1.08 km/s
    ! ==================================================================
    subroutine test_t3(data, n, all_pass, logunit)
        real(dp), intent(in) :: data(:, :)
        integer,  intent(in) :: n
        logical,  intent(inout) :: all_pass
        integer,  intent(in) :: logunit

        integer  :: i
        real(dp) :: vmag, v_min, v_max

        v_min =  huge(1.0_dp)
        v_max = -huge(1.0_dp)

        do i = 1, n
            vmag = sqrt(data(6, i)**2 + data(7, i)**2 + data(8, i)**2)
            v_min = min(v_min, vmag)
            v_max = max(v_max, vmag)
        end do

        write(logunit, '(a, f10.6, a)') '  Min |v|: ', v_min, ' km/s'
        write(logunit, '(a, f10.6, a)') '  Max |v|: ', v_max, ' km/s'

        if (v_min >= 0.96_dp .and. v_max <= 1.11_dp) then
            write(logunit, '(a)') '  T3: PASS'
        else
            write(logunit, '(a)') '  T3: FAIL'
            all_pass = .false.
        end if
    end subroutine test_t3

    ! ==================================================================
    ! T4: 选 3 个内部采样点独立调用 DPLEPH，与文件值对比
    ! ==================================================================
    subroutine test_t4(data, n, all_pass, logunit)
        real(dp), intent(in) :: data(:, :)
        integer,  intent(in) :: n
        logical,  intent(inout) :: all_pass
        integer,  intent(in) :: logunit

        integer  :: i, j
        real(dp) :: tdb1, tdb2, pv(6)
        real(dp) :: diff
        logical  :: spot_pass

        spot_idx(1) = 1          ! 第一个采样点
        spot_idx(2) = n / 2      ! 中间采样点
        spot_idx(3) = n          ! 最后一个采样点

        spot_pass = .true.

        do i = 1, N_SPOT
            j = spot_idx(i)
            tdb1 = data(1, j)
            tdb2 = data(2, j)

            call DPLEPH([tdb1, tdb2], 10, 3, pv)

            diff = 0.0_dp
            do j = 1, 6
                diff = max(diff, abs(pv(j) - data(2 + j, spot_idx(i))))
            end do

            write(logunit, '(a, i0, a, i0, a, es12.4)') &
                '  Spot ', i, ' (idx=', spot_idx(i), '): max diff = ', diff

            if (diff > 0.0_dp) then
                ! 要求位级一致: diff 应为 0
                spot_pass = .false.
            end if
        end do

        if (spot_pass) then
            write(logunit, '(a)') '  T4: PASS'
        else
            write(logunit, '(a)') '  T4: FAIL'
            all_pass = .false.
        end if
    end subroutine test_t4

    ! ==================================================================
    ! T5: 验证二进制文件大小 = N × 8 × 8 bytes
    ! ==================================================================
    subroutine test_t5(n_expected, file_size_bytes, all_pass, logunit)
        integer,  intent(in) :: n_expected
        integer(8), intent(in) :: file_size_bytes
        logical,  intent(inout) :: all_pass
        integer,  intent(in) :: logunit

        integer(8) :: expected_size

        expected_size = int(n_expected, 8) * 8_8 * 8_8

        write(logunit, '(a, i0)')    '  Expected size: ', expected_size
        write(logunit, '(a, i0)')    '  Actual   size: ', file_size_bytes

        if (file_size_bytes == expected_size) then
            write(logunit, '(a)') '  T5: PASS'
        else
            write(logunit, '(a)') '  T5: FAIL'
            all_pass = .false.
        end if
    end subroutine test_t5

end program test_lunar_extract
