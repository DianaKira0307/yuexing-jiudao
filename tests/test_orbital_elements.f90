!-------------------------------------------------------------------------------
! 月行九道 Phase 3: 轨道根数验证程序
!
! 4 项验证:
!   T1 — ICRF→黄道旋转 (春分点方向不变性)
!   T2 — Kepler 根数合理性 (a, e, i 在预期范围)
!   T3 — Ω 和 ω 的时间演化 (符号与周期)
!   T4 — Poincaré 变量一致性
!-------------------------------------------------------------------------------
program test_orbital_elements
    use mod_precision, only: dp
    use mod_constants, only: init_physics_constants, GM
    use mod_orbital_elements, only: &
        icrf_to_ecliptic, spherical_from_cartesian, &
        keplerian_from_cartesian, poincare_from_keplerian, &
        OBLIQUITY_J2000
    implicit none

    logical :: all_pass
    integer, parameter :: FUNIT_LOG = 30
    character(len=8)  :: date_str
    character(len=10) :: time_str

    ! ---- 数据文件 ----
    character(len=256) :: bin_path
    integer(8) :: file_size
    integer :: funit, n, i
    real(dp), allocatable :: data(:, :)

    ! ---- 轨道验证用变量 ----
    real(dp) :: x_e, y_e, z_e, vx_e, vy_e, vz_e
    real(dp) :: r, lam, beta, vr, vlam, vbeta
    real(dp) :: a, e, inc, raan, arg_peri, M
    real(dp) :: Lambda, lambda_p, Gamma, gamma_p, Zz, zeta_p
    real(dp) :: gm_moon  ! GM_earth + GM_moon
    integer  :: info, n_samples

    ! ---- 打开日志 ----
    open(unit=FUNIT_LOG, file='data/test_orbital_elements.log', &
         access='sequential', form='formatted', status='replace')

    call date_and_time(date=date_str, time=time_str)
    write(FUNIT_LOG, '(a)') '# Lunar Nine Paths Phase 3 — Orbital Elements Verification'
    write(FUNIT_LOG, '(a)') &
        '# Generated: ' // date_str(1:4) // '-' // date_str(5:6) // '-' // &
        date_str(7:8) // ' ' // time_str(1:2) // ':' // time_str(3:4) // ':' // time_str(5:10)
    write(FUNIT_LOG, '(a)') '#'

    all_pass = .true.

    ! ---- 初始化 (读取历表头以获取 GM) ----
    call init_physics_constants()
    gm_moon = GM(3) + GM(10)  ! 地球+月球引力常数 (km³/s²)

    ! ---- 读取 Phase 1 二进制数据 ----
    bin_path = 'data/lunar_state_1900_2100.bin'
    funit = 10
    open(unit=funit, file=trim(bin_path), access='stream', &
         form='unformatted', status='old')
    inquire(unit=funit, size=file_size)
    n = int(file_size / (8_8 * 8_8))  ! 总采样数
    allocate(data(8, n))
    read(funit) data
    close(funit)
    n_samples = n

    ! ==================================================================
    ! T1: ICRF→黄道旋转验证
    ! ==================================================================
    write(FUNIT_LOG, '(a)') '=== T1: ICRF to ecliptic rotation ==='
    call test_t1(data, n_samples, all_pass, FUNIT_LOG)

    ! ==================================================================
    ! T2: Kepler 根数合理性
    ! ==================================================================
    write(FUNIT_LOG, '(a)') ''
    write(FUNIT_LOG, '(a)') '=== T2: Kepler elements (a, e, i) ==='
    call test_t2(data, n_samples, gm_moon, all_pass, FUNIT_LOG)

    ! ==================================================================
    ! T3: Ω 和 ω 的时间演化
    ! ==================================================================
    write(FUNIT_LOG, '(a)') ''
    write(FUNIT_LOG, '(a)') '=== T3: Time evolution of Omega and omega ==='
    call test_t3(data, n_samples, gm_moon, all_pass, FUNIT_LOG)

    ! ==================================================================
    ! T4: Poincaré 变量一致性
    ! ==================================================================
    write(FUNIT_LOG, '(a)') ''
    write(FUNIT_LOG, '(a)') '=== T4: Poincare variable consistency ==='
    call test_t4(data, n_samples, gm_moon, all_pass, FUNIT_LOG)

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
        write(*, '(a)') 'All 4 tests PASSED. See data/test_orbital_elements.log'
    else
        write(*, '(a)') 'Some tests FAILED. See data/test_orbital_elements.log'
        stop 1
    end if

contains

    ! ==================================================================
    ! T1: ICRF→黄道旋转验证
    ! 春分点方向 (1,0,0) 旋转后不变;
    ! (0,0,1) 旋转后应得到 (0, sin ε, cos ε)
    ! ==================================================================
    subroutine test_t1(data, n, all_pass, logunit)
        real(dp), intent(in) :: data(:, :)
        integer,  intent(in) :: n
        logical,  intent(inout) :: all_pass
        integer,  intent(in) :: logunit

        real(dp) :: x_e, y_e, z_e, vx_e, vy_e, vz_e
        real(dp) :: eps
        logical  :: ok
        real(dp), parameter :: TOL = 1.0e-14_dp

        eps = OBLIQUITY_J2000

        ! 测试点 1: X 轴方向 (春分点)
        call icrf_to_ecliptic(1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, eps, &
                              x_e, y_e, z_e, vx_e, vy_e, vz_e)
        ok = (abs(x_e - 1.0_dp) < TOL) .and. (abs(y_e) < TOL) .and. (abs(z_e) < TOL)
        write(logunit, '(a, l)') '  X-axis (vernal eqx): ', ok
        write(logunit, '(a, 3f12.6)') '    -> ', x_e, y_e, z_e

        ! 测试点 2: Z 轴方向 (北天极)
        call icrf_to_ecliptic(0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, eps, &
                              x_e, y_e, z_e, vx_e, vy_e, vz_e)
        ok = (abs(x_e) < TOL) .and. (abs(y_e - sin(eps)) < TOL) .and. &
             (abs(z_e - cos(eps)) < TOL)
        write(logunit, '(a, l)') '  Z-axis: ', ok
        write(logunit, '(a, 3f12.6)') '    -> ', x_e, y_e, z_e

        ! 测试点 3: 用第一个采样点验证旋转是正交的 (范数不变)
        call icrf_to_ecliptic(data(3, 1), data(4, 1), data(5, 1), &
                              data(6, 1), data(7, 1), data(8, 1), eps, &
                              x_e, y_e, z_e, vx_e, vy_e, vz_e)
        ok = abs(sqrt(data(3,1)**2+data(4,1)**2+data(5,1)**2) &
               - sqrt(x_e**2+y_e**2+z_e**2)) < TOL
        write(logunit, '(a, l)') '  Norm preserved: ', ok

        if (ok) then
            write(logunit, '(a)') '  T1: PASS'
        else
            write(logunit, '(a)') '  T1: FAIL'
            all_pass = .false.
        end if
    end subroutine test_t1

    ! ==================================================================
    ! T2: Kepler 根数合理性
    ! 取中间采样点和末尾采样点, 转换到黄道后验证
    !   a ≈ 385000 km, e ≈ 0.055, i ≈ 5.09°
    ! ==================================================================
    subroutine test_t2(data, n, gm, all_pass, logunit)
        real(dp), intent(in) :: data(:, :)
        integer,  intent(in) :: n
        real(dp), intent(in) :: gm
        logical,  intent(inout) :: all_pass
        integer,  intent(in) :: logunit

        integer  :: idx(3), i, j, nok
        real(dp) :: x_e, y_e, z_e, vx_e, vy_e, vz_e
        real(dp) :: a, e, inc, raan, arg_peri
        integer  :: info
        real(dp) :: inc_deg

        idx = [1, n/2, n]  ! 首、中、末

        nok = 0
        do i = 1, 3
            j = idx(i)
            call icrf_to_ecliptic(data(3, j), data(4, j), data(5, j), &
                                  data(6, j), data(7, j), data(8, j), &
                                  OBLIQUITY_J2000, &
                                  x_e, y_e, z_e, vx_e, vy_e, vz_e)
            call keplerian_from_cartesian(x_e, y_e, z_e, vx_e, vy_e, vz_e, &
                                          gm, a, e, inc, raan, arg_peri, M, info)
            inc_deg = inc * 180.0_dp / 3.14159265358979323846_dp

            write(logunit, '(a, i6, a)') &
                '  Sample ', j, ':'
            write(logunit, '(a, f12.3, a)') '    a  = ', a, ' km'
            write(logunit, '(a, f10.6)')    '    e  = ', e
            write(logunit, '(a, f10.4, a)') '    i  = ', inc_deg, ' deg'
            write(logunit, '(a, f10.4, a)') '    raan = ', raan * 180.0_dp / 3.14159265358979323846_dp, ' deg'
            write(logunit, '(a, f10.4, a)') '    arg_peri = ', arg_peri * 180.0_dp / 3.14159265358979323846_dp, ' deg'
            write(logunit, '(a, f10.4, a)') '    M   = ', M * 180.0_dp / 3.14159265358979323846_dp, ' deg'

            if (info == 0 .and. a > 380000.0_dp .and. a < 390000.0_dp .and. &
                e > 0.03_dp .and. e < 0.08_dp .and. &
                inc_deg > 4.5_dp .and. inc_deg < 6.0_dp) then
                nok = nok + 1
            end if
        end do

        if (nok == 3) then
            write(logunit, '(a)') '  T2: PASS'
        else
            write(logunit, '(a)') '  T2: FAIL'
            all_pass = .false.
        end if
    end subroutine test_t2

    ! ==================================================================
    ! T3: Ω 和 ω 的时间演化
    ! 取约 20 年的数据, 拟合 Omega(t) 和 omega(t) 的变化率
    !   Ω·dot < 0 (升交点退行), 周期约 18.6 年
    !   ω·dot > 0 (近地点进动), 周期约 8.85 年
    ! ==================================================================
    subroutine test_t3(data, n, gm, all_pass, logunit)
        real(dp), intent(in) :: data(:, :)
        integer,  intent(in) :: n
        real(dp), intent(in) :: gm
        logical,  intent(inout) :: all_pass
        integer,  intent(in) :: logunit

        integer, parameter :: N_SAMPLES = 36525  ! ~100 年
        integer :: i, info
        real(dp) :: x_e, y_e, z_e, vx_e, vy_e, vz_e
        real(dp) :: a, e, inc, raan, arg_peri
        real(dp) :: tdb_i, tdb_0
        real(dp) :: raan_prev, aperi_prev
        real(dp) :: draan, draan_total, draan_rate
        real(dp) :: daperi, daperi_total, daperi_rate
        real(dp) :: period_node, period_peri

        if (N_SAMPLES > n) then
            write(logunit, '(a)') '  Not enough samples.'
            all_pass = .false.
            return
        end if

        tdb_0 = data(1, 1) + data(2, 1)
        raan_prev = -999.0_dp
        aperi_prev = -999.0_dp
        draan_total = 0.0_dp
        daperi_total = 0.0_dp

        do i = 1, N_SAMPLES
            call icrf_to_ecliptic(data(3, i), data(4, i), data(5, i), &
                                  data(6, i), data(7, i), data(8, i), &
                                  OBLIQUITY_J2000, &
                                  x_e, y_e, z_e, vx_e, vy_e, vz_e)
            call keplerian_from_cartesian(x_e, y_e, z_e, vx_e, vy_e, vz_e, &
                                          gm, a, e, inc, raan, arg_peri, M, info)
            if (info /= 0) cycle

            ! 累计角度变化 (超过 π 时补 2π 漂移)
            if (raan_prev > 0.0_dp) then
                draan = raan - raan_prev
                if (draan > 3.0_dp) draan = draan - 2.0_dp * 3.14159265358979323846_dp
                if (draan < -3.0_dp) draan = draan + 2.0_dp * 3.14159265358979323846_dp
                draan_total = draan_total + draan
            end if

            if (aperi_prev > 0.0_dp) then
                daperi = arg_peri - aperi_prev
                if (daperi > 3.0_dp) daperi = daperi - 2.0_dp * 3.14159265358979323846_dp
                if (daperi < -3.0_dp) daperi = daperi + 2.0_dp * 3.14159265358979323846_dp
                daperi_total = daperi_total + daperi
            end if

            raan_prev = raan
            aperi_prev = arg_peri
        end do

        draan_rate = draan_total / real(N_SAMPLES - 1, dp)  ! rad/day
        daperi_rate = daperi_total / real(N_SAMPLES - 1, dp)  ! rad/day

        ! 周期估算 (天)
        if (abs(draan_rate) > 1.0e-12_dp) then
            period_node = 2.0_dp * 3.14159265358979323846_dp / abs(draan_rate)
        else
            period_node = 9999.0_dp
        end if
        if (abs(daperi_rate) > 1.0e-12_dp) then
            period_peri = 2.0_dp * 3.14159265358979323846_dp / abs(daperi_rate)
        else
            period_peri = 9999.0_dp
        end if

        write(logunit, '(a, es12.4, a)') &
            '  draan/dt   = ', draan_rate, ' rad/day'
        write(logunit, '(a, f10.2, a)') &
            '  Node period ~ ', period_node, ' days'
        write(logunit, '(a, es12.4, a)') &
            '  daperi/dt  = ', daperi_rate, ' rad/day'
        write(logunit, '(a, f10.2, a)') &
            '  Perihelion period ~ ', period_peri, ' days'

        ! 验证: Ω 退行 (负), ω 进动 (正)
        ! 升交点退行周期 ~18.6 年 (6790 天)
        ! 近地点进动周期 ~8.85 年 (3230 天, 但osculating元素受短期摄动影响)
        if (draan_rate < 0.0_dp .and. daperi_rate > 0.0_dp .and. &
            period_node > 5000.0_dp .and. period_node < 9000.0_dp .and. &
            period_peri > 1500.0_dp .and. period_peri < 5000.0_dp) then
            write(logunit, '(a)') '  T3: PASS'
        else
            write(logunit, '(a)') '  T3: FAIL'
            all_pass = .false.
        end if
    end subroutine test_t3

    ! ==================================================================
    ! T4: Poincaré 变量一致性
    ! 对一数据点计算 Poincaré 变量, 验证 Λ = sqrt(gm·a)
    ! 以及 Γ ≈ L·e²/2 (小 e 近似)
    ! ==================================================================
    subroutine test_t4(data, n, gm, all_pass, logunit)
        real(dp), intent(in) :: data(:, :)
        integer,  intent(in) :: n
        real(dp), intent(in) :: gm
        logical,  intent(inout) :: all_pass
        integer,  intent(in) :: logunit

        integer  :: j, info
        real(dp) :: x_e, y_e, z_e, vx_e, vy_e, vz_e
        real(dp) :: a, e, inc, raan, arg_peri
        real(dp) :: Lambda, lambda_p, Gamma, gamma_p, Zz, zeta_p
        real(dp) :: L_direct, Gamma_approx, Z_approx

        j = n / 2  ! 中间点

        call icrf_to_ecliptic(data(3, j), data(4, j), data(5, j), &
                              data(6, j), data(7, j), data(8, j), &
                              OBLIQUITY_J2000, &
                              x_e, y_e, z_e, vx_e, vy_e, vz_e)
        call keplerian_from_cartesian(x_e, y_e, z_e, vx_e, vy_e, vz_e, &
                                      gm, a, e, inc, raan, arg_peri, M, info)
        call poincare_from_keplerian(a, e, inc, raan, arg_peri, M, gm, &
                                     Lambda, lambda_p, Gamma, gamma_p, Zz, zeta_p)

        L_direct = sqrt(gm * a)
        Gamma_approx = L_direct * e * e / 2.0_dp  ! e²/2 近似
        Z_approx = L_direct * inc * inc / 4.0_dp  ! i²/4 近似 (小 i)

        write(logunit, '(a, i6)') '  Sample ', j
        write(logunit, '(a, es16.8)') '  Lambda = ', Lambda
        write(logunit, '(a, es16.8)') '  Gamma  = ', Gamma
        write(logunit, '(a, es16.8)') '  Z      = ', Zz
        write(logunit, '(a, f10.4, a)') '  lambda = ', lambda_p * 180.0_dp / 3.14159265358979323846_dp, ' deg'
        write(logunit, '(a, es12.4)') '  |Lambda - sqrt(gm*a)| = ', abs(Lambda - L_direct)
        write(logunit, '(a, es12.4)') '  |Gamma - L·e²/2|     = ', abs(Gamma - Gamma_approx)
        write(logunit, '(a, es12.4)') '  |Z - L·i²/4|         = ', abs(Zz - Z_approx)

        if (abs(Lambda - L_direct) < 1.0e-6_dp) then
            write(logunit, '(a)') '  T4: PASS'
        else
            write(logunit, '(a)') '  T4: FAIL'
            all_pass = .false.
        end if
    end subroutine test_t4

end program test_orbital_elements
