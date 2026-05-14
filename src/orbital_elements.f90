!-------------------------------------------------------------------------------
! 月行九道 Phase 3: 轨道根数计算模块
! 模块: mod_orbital_elements
!
! 提供三组分析变量的计算:
!   1. 黄道球坐标 (r, λ, β, 及导数)
!   2. 经典 Kepler 轨道根数 (a, e, i, Ω, ω, M)
!   3. Poincaré 正则变量 (Λ, λ, Γ, γ, Z, ζ)
!-------------------------------------------------------------------------------

module mod_orbital_elements
    use mod_precision, only: dp
    implicit none
    private

    public :: icrf_to_ecliptic
    public :: spherical_from_cartesian
    public :: keplerian_from_cartesian
    public :: poincare_from_keplerian

    ! J2000.0 平黄赤交角 (弧度)
    ! 来源: 天文常数 23°26′21.448″ = 23.4392911111°
    real(dp), parameter, public :: OBLIQUITY_J2000 = 0.409092804_dp

    real(dp), parameter :: PI = 3.14159265358979323846_dp
    real(dp), parameter :: PI2 = 2.0_dp * PI
    real(dp), parameter :: EPS_E = 1.0e-14_dp  ! 处理奇异的小量

contains

    !-------------------------------------------------------------------
    ! icrf_to_ecliptic: ICRF (赤道) → 黄道坐标旋转
    !
    ! 绕 X 轴旋转 ε (黄赤交角):
    !   x_ecl = x_eq
    !   y_ecl = y_eq·cos ε + z_eq·sin ε
    !   z_ecl = −y_eq·sin ε + z_eq·cos ε
    !
    ! 参数:
    !   eps    — 黄赤交角 (弧度), 可用模块常量 OBLIQUITY_J2000
    !-------------------------------------------------------------------
    subroutine icrf_to_ecliptic(x, y, z, vx, vy, vz, eps, &
                                x_e, y_e, z_e, vx_e, vy_e, vz_e)
        real(dp), intent(in)  :: x, y, z, vx, vy, vz, eps
        real(dp), intent(out) :: x_e, y_e, z_e, vx_e, vy_e, vz_e

        real(dp) :: c, s
        c = cos(eps)
        s = sin(eps)

        x_e   = x
        y_e   = c * y + s * z
        z_e   = -s * y + c * z
        vx_e  = vx
        vy_e  = c * vy + s * vz
        vz_e  = -s * vy + c * vz
    end subroutine icrf_to_ecliptic

    !-------------------------------------------------------------------
    ! spherical_from_cartesian: 笛卡尔 → 黄道球坐标
    !
    ! 输出:
    !   r      — 地心距 (km)
    !   lam    — 黄经 λ, [0, 2π) (rad)
    !   beta   — 黄纬 β, [−π/2, π/2] (rad)
    !   vr     — 径向速度 (km/s)
    !   vlam   — λ 方向速度分量 (km/s), λ·增长方向为正
    !   vbeta  — β 方向速度分量 (km/s), β·增长方向为正
    !-------------------------------------------------------------------
    subroutine spherical_from_cartesian(x, y, z, vx, vy, vz, &
                                        r, lam, beta, vr, vlam, vbeta)
        real(dp), intent(in)  :: x, y, z, vx, vy, vz
        real(dp), intent(out) :: r, lam, beta, vr, vlam, vbeta

        real(dp) :: rhox, rhoy, rho, c_lam, s_lam, c_beta, s_beta

        r = sqrt(x * x + y * y + z * z)

        ! 球坐标基向量
        rhox = sqrt(x * x + y * y)  ! 在 XY 平面投影
        if (rhox > 0.0_dp) then
            c_lam = x / rhox
            s_lam = y / rhox
            lam = atan2(y, x)  ! [0, 2π)
            if (lam < 0.0_dp) lam = lam + PI2

            beta = atan2(z, rhox)
            c_beta = rhox / r
            s_beta = z / r
        else
            ! 极点上: λ 无定义，取 0
            lam = 0.0_dp
            c_lam = 1.0_dp
            s_lam = 0.0_dp
            beta = atan2(z, 0.0_dp)  ! z 正负决定 ±π/2
            c_beta = 0.0_dp
            s_beta = sign(1.0_dp, z)
        end if

        ! 速度在球坐标基下的分量
        ! vr    = r̂·v = (x·vx + y·vy + z·vz) / r
        ! vlam  = λ̂·v = (−sin λ · vx + cos λ · vy)
        ! vbeta = β̂·v = (−cos λ·sin β·vx − sin λ·sin β·vy + cos β·vz)
        vr    = (x * vx + y * vy + z * vz) / r
        vlam  = -s_lam * vx + c_lam * vy
        vbeta = -c_lam * s_beta * vx - s_lam * s_beta * vy + c_beta * vz
    end subroutine spherical_from_cartesian

    !-------------------------------------------------------------------
    ! keplerian_from_cartesian: 笛卡尔 → Kepler 轨道根数
    !
    ! 输入为黄道坐标系 (或任意惯性系)。
    !
    ! 输出:
    !   a      — 半长轴 (km)
    !   e      — 偏心率
    !   inc    — 倾角, [0, π] (rad)
    !   Omega  — 升交点黄经, [0, 2π) (rad)
    !   omega  — 近地点幅角, [0, 2π) (rad)
    !   M      — 平近点角, [0, 2π) (rad)
    !   info   — 返回码
    !             0 = 成功
    !             1 = 角动量为零 (径向轨迹)
    !             2 = 双曲线/抛物线轨道 (a <= 0)
    !             3 = 偏心率异常
    !-------------------------------------------------------------------
    subroutine keplerian_from_cartesian(x, y, z, vx, vy, vz, gm, &
                                        a, e, inc, raan, arg_peri, M, info)
        real(dp), intent(in)  :: x, y, z, vx, vy, vz, gm
        real(dp), intent(out) :: a, e, inc, raan, arg_peri, M
        integer,  intent(out) :: info

        real(dp) :: r, v2, energy
        real(dp) :: hx, hy, hz, h_mag, h_mag2
        real(dp) :: nx, ny, n_mag
        real(dp) :: ex, ey, ez, e_mag
        real(dp) :: cos_i
        real(dp) :: cos_omega, sin_omega
        real(dp) :: cos_f, sin_f
        real(dp) :: cos_E, sin_E
        real(dp) :: true_anom, ecc_anom
        real(dp) :: RDV  ! r·v

        info = 0

        ! ---- 1. 位置和速度标量 ----
        r = sqrt(x * x + y * y + z * z)
        v2 = vx * vx + vy * vy + vz * vz
        RDV = x * vx + y * vy + z * vz

        ! ---- 2. 角动量 h = r × v ----
        hx = y * vz - z * vy
        hy = z * vx - x * vz
        hz = x * vy - y * vx
        h_mag2 = hx * hx + hy * hy + hz * hz
        h_mag  = sqrt(h_mag2)

        if (h_mag < EPS_E) then
            info = 1
            a = 0.0_dp; e = 0.0_dp; inc = 0.0_dp
            raan = 0.0_dp; arg_peri = 0.0_dp; M = 0.0_dp
            return
        end if

        ! ---- 3. 半长轴 (从能量) ----
        energy = 0.5_dp * v2 - gm / r
        a = -gm / (2.0_dp * energy)
        if (a <= 0.0_dp) then
            info = 2
            return
        end if

        ! ---- 4. 偏心率向量 e = (v × h) / GM − r̂ ----
        ex = (vy * hz - vz * hy) / gm - x / r
        ey = (vz * hx - vx * hz) / gm - y / r
        ez = (vx * hy - vy * hx) / gm - z / r
        e_mag = sqrt(ex * ex + ey * ey + ez * ez)

        if (e_mag < EPS_E) then
            info = 3
            return
        end if
        e = e_mag

        ! ---- 5. 倾角 ----
        ! h = |h|·(sin i·sin Ω, −sin i·cos Ω, cos i)
        cos_i = hz / h_mag
        inc = acos(cos_i)
        if (inc < 0.0_dp) inc = inc + PI

        ! ---- 6. 升交点黄经 Ω ----
        ! 节线 n = (0, 0, 1) × h = (−hy, hx, 0)
        nx = -hy
        ny = hx
        n_mag = sqrt(nx * nx + ny * ny)

        if (n_mag < EPS_E) then
            ! 赤道轨道: Ω 无定义
            raan = 0.0_dp
        else
            raan = atan2(ny, nx)  ! 即 atan2(hx, −hy)
            if (raan < 0.0_dp) raan = raan + PI2
        end if

        ! ---- 7. 近地点幅角 ω ----
        if (n_mag < EPS_E) then
            ! 赤道轨道: 从 e_x, e_y 方向确定
            arg_peri = atan2(ey, ex)
            if (arg_peri < 0.0_dp) arg_peri = arg_peri + PI2
        else
            ! cos ω = (n·e) / (|n|·|e|)
            cos_omega = (nx * ex + ny * ey) / (n_mag * e_mag)

            ! sin ω: 轨道平面内垂直于 n 的方向与 e 的点积
            ! sin ω = (ĥ × n̂) · ê  (ĥ = h/|h|, n̂ = n/|n|, ê = e/|e|)
            ! 由于 nz = 0:
            sin_omega = ((-hz / h_mag) * (ny / n_mag)) * (ex / e_mag) + &
                        ( (hz / h_mag) * (nx / n_mag)) * (ey / e_mag) + &
                        ((hx * ny - hy * nx) / (h_mag * n_mag)) * (ez / e_mag)

            arg_peri = atan2(sin_omega, cos_omega)
            if (arg_peri < 0.0_dp) arg_peri = arg_peri + PI2
        end if

        ! ---- 8. 平近点角 M (先经过真近点角 f 和偏近点角 E) ----
        ! cos f = (e·r) / (e·r)
        cos_f = (ex * x + ey * y + ez * z) / (e_mag * r)

        ! sin f: (e × r) · ĥ = ((e·r)(r·v) − (e·v)r²) / (e·r·h)
        sin_f = ((ex * x + ey * y + ez * z) * RDV - &
                 (ex * vx + ey * vy + ez * vz) * r * r) / (e_mag * r * h_mag)

        true_anom = atan2(sin_f, cos_f)
        if (true_anom < 0.0_dp) true_anom = true_anom + PI2

        ! 偏近点角 E
        cos_E = (e + cos_f) / (1.0_dp + e * cos_f)
        sin_E = sqrt(1.0_dp - e * e) * sin_f / (1.0_dp + e * cos_f)
        ecc_anom = atan2(sin_E, cos_E)
        if (ecc_anom < 0.0_dp) ecc_anom = ecc_anom + PI2

        ! 平近点角 M (Kepler 方程)
        M = ecc_anom - e * sin_E
        if (M < 0.0_dp) M = M + PI2

    end subroutine keplerian_from_cartesian

    !-------------------------------------------------------------------
    ! poincare_from_keplerian: Kepler 根数 → Poincaré 正则变量
    !
    ! Poincaré 变量消去了圆轨道和赤道轨道的奇异:
    !   Λ   = L                  (与平黄经共轭)
    !   λ   = M + ω + Ω          (平黄经)
    !   Γ   = L · (1 − √(1−e²)) (与 −(ω+Ω) 共轭)
    !   γ   = −(ω + Ω)
    !   Z   = L · √(1−e²) · (1 − cos i)   (与 −Ω 共轭)
    !   ζ   = −Ω
    !
    ! 其中 L = √(GM·a), GM 为引力常数与中心天体质量的乘积。
    !
    ! 所有角度归算到 [0, 2π)。
    !-------------------------------------------------------------------
    subroutine poincare_from_keplerian(a, e, inc, raan, arg_peri, M, gm, &
                                       Lambda, lambda_p, Gamma, gamma_p, Z, zeta_p)
        real(dp), intent(in)  :: a, e, inc, raan, arg_peri, M, gm
        real(dp), intent(out) :: Lambda, lambda_p, Gamma, gamma_p, Z, zeta_p

        real(dp) :: L, sqrt_1me2

        L = sqrt(gm * a)
        sqrt_1me2 = sqrt(1.0_dp - e * e)

        ! 作用量
        Lambda = L
        Gamma  = L * (1.0_dp - sqrt_1me2)
        Z      = L * sqrt_1me2 * (1.0_dp - cos(inc))

        ! 角度
        lambda_p = M + arg_peri + raan
        if (lambda_p < 0.0_dp) lambda_p = lambda_p + PI2
        if (lambda_p >= PI2)    lambda_p = lambda_p - PI2

        gamma_p = -(arg_peri + raan)
        if (gamma_p < 0.0_dp) gamma_p = gamma_p + PI2
        if (gamma_p >= PI2)    gamma_p = gamma_p - PI2

        zeta_p = -raan
        if (zeta_p < 0.0_dp) zeta_p = zeta_p + PI2
        if (zeta_p >= PI2)    zeta_p = zeta_p - PI2

    end subroutine poincare_from_keplerian

end module mod_orbital_elements
