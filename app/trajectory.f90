!-------------------------------------------------------------------------------
! 月行九道 — 黄道轨迹 (λ, β)
!
! 输出月球在黄道球坐标中的完整轨迹，供 Python 分不同时间窗口绘制
! "月行九道"在天球上的几何图像。
!
! 输出: data/traj_ecliptic_j2000.txt  固定 J2000 黄道坐标系
!       data/traj_ecliptic_date.txt   SOFA 历元黄道坐标系
!       列: day  lambda(deg)  beta(deg)
!-------------------------------------------------------------------------------
program lunar_trajectory
    use mod_precision,        only: dp
    use mod_orbital_elements, only: icrf_to_ecliptic, icrf_to_ecliptic_date, &
                                    spherical_from_cartesian, &
                                    OBLIQUITY_J2000
    implicit none

    character(len=*), parameter :: BIN_PATH = 'data/lunar_state_1900_2100.bin'
    real(dp), parameter :: RAD2DEG = 180.0_dp / 3.14159265358979323846_dp

    integer(8) :: file_size
    integer    :: n_total, funit, fout_j2000, fout_date, i
    real(dp), allocatable :: buf(:,:)
    real(dp) :: xe, ye, ze, vxe, vye, vze
    real(dp) :: lam, beta_val
    real(dp) :: r_val, vr, vlam, vbeta  ! required by interface, not used here
    real(dp) :: tdb0, day_i

    write(*, '(a)') '=== 月行九道 — Ecliptic Trajectory ==='

    inquire(file=BIN_PATH, size=file_size)
    if (file_size <= 0) stop 'ERROR: data file not found. Run: make run'
    n_total = int(file_size / 64_8)
    write(*, '(a, i0, a)') 'Loading ', n_total, ' samples...'

    allocate(buf(8, n_total))
    open(newunit=funit, file=BIN_PATH, access='stream', form='unformatted', status='old')
    read(funit) buf
    close(funit)

    tdb0 = buf(1,1) + buf(2,1)

    open(newunit=fout_j2000, file='data/traj_ecliptic_j2000.txt', status='replace', action='write')
    call write_header(fout_j2000, 'Frame: fixed J2000 mean ecliptic approximation')

    open(newunit=fout_date, file='data/traj_ecliptic_date.txt', status='replace', action='write')
    call write_header(fout_date, 'Frame: SOFA IAU 2006 mean ecliptic/equinox of date')

    do i = 1, n_total
        day_i = (buf(1,i) + buf(2,i)) - tdb0

        call icrf_to_ecliptic(buf(3,i), buf(4,i), buf(5,i), &
                               buf(6,i), buf(7,i), buf(8,i), &
                               OBLIQUITY_J2000, xe, ye, ze, vxe, vye, vze)
        call spherical_from_cartesian(xe, ye, ze, vxe, vye, vze, &
                                       r_val, lam, beta_val, vr, vlam, vbeta)
        write(fout_j2000, '(f12.4, 2x, f12.6, 2x, f12.6)') &
            day_i, lam * RAD2DEG, beta_val * RAD2DEG

        call icrf_to_ecliptic_date(buf(1,i), buf(2,i), &
                                   buf(3,i), buf(4,i), buf(5,i), &
                                   buf(6,i), buf(7,i), buf(8,i), &
                                   xe, ye, ze, vxe, vye, vze)
        call spherical_from_cartesian(xe, ye, ze, vxe, vye, vze, &
                                       r_val, lam, beta_val, vr, vlam, vbeta)
        write(fout_date, '(f12.4, 2x, f12.6, 2x, f12.6)') &
            day_i, lam * RAD2DEG, beta_val * RAD2DEG
    end do
    close(fout_j2000)
    close(fout_date)

    deallocate(buf)
    write(*, '(a, i0, a)') 'Written ', n_total, ' points  →  data/traj_ecliptic_j2000.txt'
    write(*, '(a, i0, a)') 'Written ', n_total, ' points  →  data/traj_ecliptic_date.txt'

contains

    subroutine write_header(fout, frame_note)
        integer, intent(in) :: fout
        character(len=*), intent(in) :: frame_note

        write(fout, '(a)') '# 月行九道 — Ecliptic trajectory (λ, β)'
        write(fout, '(a)') '# DE440 1900–2100, 1-day steps'
        write(fout, '(a)') '# ' // frame_note
        write(fout, '(a)') '#'
        write(fout, '(a)') '# day          lambda(deg)   beta(deg)'
    end subroutine write_header
end program lunar_trajectory
