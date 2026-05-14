!-------------------------------------------------------------------------------
! 月行九道 — 黄道交点漂移
!
! 检测月球穿越黄道平面 (β=0) 的时刻，线性插值求对应黄经 λ。
! 交点的持续漂移即升交点退行，是"月行九道"最直接的观测体现。
!
! 输出: data/node_crossings.txt
!       列: day  lambda(deg)  type  (A=升交点 β:−→+, D=降交点 β:+→−)
!-------------------------------------------------------------------------------
program lunar_nodes
    use mod_precision,        only: dp
    use mod_orbital_elements, only: icrf_to_ecliptic, spherical_from_cartesian, &
                                    OBLIQUITY_J2000
    implicit none

    character(len=*), parameter :: BIN_PATH = 'data/lunar_state_1900_2100.bin'
    real(dp), parameter :: PI     = 3.14159265358979323846_dp
    real(dp), parameter :: PI2    = 2.0_dp * PI
    real(dp), parameter :: RAD2DEG = 180.0_dp / PI

    integer(8) :: file_size
    integer    :: n_total, funit, fout, i, n_cross
    real(dp), allocatable :: buf(:,:), beta(:), lam(:), t_days(:)
    real(dp) :: xe, ye, ze, vxe, vye, vze
    real(dp) :: lam_val, beta_val
    real(dp) :: r_val, vr, vlam, vbeta  ! required by interface, not used here
    real(dp) :: tdb0, frac, t_cross, lam_cross, delta_lam

    write(*, '(a)') '=== 月行九道 — Node Crossings ==='

    inquire(file=BIN_PATH, size=file_size)
    if (file_size <= 0) stop 'ERROR: data file not found. Run: make run'
    n_total = int(file_size / 64_8)
    write(*, '(a, i0, a)') 'Loading ', n_total, ' samples...'

    allocate(buf(8, n_total), beta(n_total), lam(n_total), t_days(n_total))
    open(newunit=funit, file=BIN_PATH, access='stream', form='unformatted', status='old')
    read(funit) buf
    close(funit)

    tdb0 = buf(1,1) + buf(2,1)
    do i = 1, n_total
        t_days(i) = (buf(1,i) + buf(2,i)) - tdb0
        call icrf_to_ecliptic(buf(3,i), buf(4,i), buf(5,i), &
                               buf(6,i), buf(7,i), buf(8,i), &
                               OBLIQUITY_J2000, xe, ye, ze, vxe, vye, vze)
        call spherical_from_cartesian(xe, ye, ze, vxe, vye, vze, &
                                       r_val, lam_val, beta_val, vr, vlam, vbeta)
        beta(i) = beta_val
        lam(i)  = lam_val   ! λ ∈ [0, 2π)
    end do
    deallocate(buf)

    open(newunit=fout, file='data/node_crossings.txt', status='replace', action='write')
    write(fout, '(a)') '# 月行九道 — Ecliptic node crossings (β = 0)'
    write(fout, '(a)') '# DE440 1900–2100, linear interpolation between 1-day steps'
    write(fout, '(a)') '#'
    write(fout, '(a)') '# day          lambda(deg)   type'
    write(fout, '(a)') '# ----------   -----------   ----  (A=ascending, D=descending)'

    n_cross = 0
    do i = 1, n_total - 1
        if (beta(i) * beta(i+1) >= 0.0_dp) cycle  ! 同号，无穿越

        ! 线性插值：β(i) + frac*(β(i+1)-β(i)) = 0
        frac = beta(i) / (beta(i) - beta(i+1))
        t_cross = t_days(i) + frac * (t_days(i+1) - t_days(i))

        ! λ 插值，处理 0/2π 边界（月球日行 ~13°，不存在真正跳变）
        delta_lam = lam(i+1) - lam(i)
        if (delta_lam > PI)  delta_lam = delta_lam - PI2
        if (delta_lam < -PI) delta_lam = delta_lam + PI2
        lam_cross = lam(i) + frac * delta_lam
        ! 规范到 [0, 360°)
        lam_cross = mod(lam_cross * RAD2DEG, 360.0_dp)
        if (lam_cross < 0.0_dp) lam_cross = lam_cross + 360.0_dp

        if (beta(i+1) > beta(i)) then
            write(fout, '(f12.4, 2x, f12.6, 2x, a1)') t_cross, lam_cross, 'A'
        else
            write(fout, '(f12.4, 2x, f12.6, 2x, a1)') t_cross, lam_cross, 'D'
        end if
        n_cross = n_cross + 1
    end do
    close(fout)

    write(*, '(a, i0, a)') 'Crossings: ', n_cross, '  →  data/node_crossings.txt'
end program lunar_nodes
