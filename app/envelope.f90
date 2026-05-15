!-------------------------------------------------------------------------------
! 月行九道 — β 振幅包络
!
! 读取地心状态矢量，转换为黄道球坐标，提取 β(t) 的局部极值。
! 包络的缓变幅度直接反映升交点退行对月球南北范围的调制。
!
! 输出: data/envelope_beta_j2000.txt  固定 J2000 黄道坐标系
!       data/envelope_beta_date.txt   SOFA 历元黄道坐标系
!       列: day  beta(deg)  type  (1=极大, -1=极小)
!-------------------------------------------------------------------------------
program lunar_envelope
    use mod_precision,        only: dp
    use mod_orbital_elements, only: icrf_to_ecliptic, icrf_to_ecliptic_date, &
                                    spherical_from_cartesian, &
                                    OBLIQUITY_J2000
    implicit none

    character(len=*), parameter :: BIN_PATH = 'data/lunar_state_1900_2100.bin'
    real(dp), parameter :: RAD2DEG = 180.0_dp / 3.14159265358979323846_dp

    integer(8) :: file_size
    integer    :: n_total, funit, i
    real(dp), allocatable :: buf(:,:), beta_j2000(:), beta_date(:), t_days(:)
    real(dp) :: xe, ye, ze, vxe, vye, vze
    real(dp) :: beta_val
    real(dp) :: r_val, lam, vr, vlam, vbeta  ! required by interface, not used here
    real(dp) :: tdb0

    write(*, '(a)') '=== 月行九道 — β Envelope ==='

    inquire(file=BIN_PATH, size=file_size)
    if (file_size <= 0) stop 'ERROR: data file not found. Run: make run'
    n_total = int(file_size / 64_8)
    write(*, '(a, i0, a)') 'Loading ', n_total, ' samples...'

    allocate(buf(8, n_total), beta_j2000(n_total), beta_date(n_total), t_days(n_total))

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
                                       r_val, lam, beta_val, vr, vlam, vbeta)
        beta_j2000(i) = beta_val * RAD2DEG

        call icrf_to_ecliptic_date(buf(1,i), buf(2,i), &
                                   buf(3,i), buf(4,i), buf(5,i), &
                                   buf(6,i), buf(7,i), buf(8,i), &
                                   xe, ye, ze, vxe, vye, vze)
        call spherical_from_cartesian(xe, ye, ze, vxe, vye, vze, &
                                       r_val, lam, beta_val, vr, vlam, vbeta)
        beta_date(i) = beta_val * RAD2DEG
    end do
    deallocate(buf)

    call write_envelope_file('data/envelope_beta_j2000.txt', &
        'Frame: fixed J2000 mean ecliptic approximation', beta_j2000)
    call write_envelope_file('data/envelope_beta_date.txt', &
        'Frame: SOFA IAU 2006 mean ecliptic/equinox of date', beta_date)

contains

    subroutine write_envelope_file(path, frame_note, beta)
        character(len=*), intent(in) :: path, frame_note
        real(dp), intent(in) :: beta(:)

        integer :: fout, i, n_out

        open(newunit=fout, file=path, status='replace', action='write')
        write(fout, '(a)') '# 月行九道 — Ecliptic latitude β envelope (local extrema)'
        write(fout, '(a)') '# DE440 1900–2100, 1-day steps'
        write(fout, '(a)') '# ' // frame_note
        write(fout, '(a)') '#'
        write(fout, '(a)') '# day         beta(deg)    type'
        write(fout, '(a)') '# ----------  -----------  ----  (1=max, -1=min)'

        n_out = 0
        do i = 2, size(beta) - 1
            if (beta(i) > beta(i-1) .and. beta(i) > beta(i+1)) then
                write(fout, '(f12.4, 2x, f12.6, 2x, i2)') t_days(i), beta(i),  1
                n_out = n_out + 1
            else if (beta(i) < beta(i-1) .and. beta(i) < beta(i+1)) then
                write(fout, '(f12.4, 2x, f12.6, 2x, i2)') t_days(i), beta(i), -1
                n_out = n_out + 1
            end if
        end do
        close(fout)

        write(*, '(a, i0, a, a)') 'Extrema: ', n_out, '  →  ', path
    end subroutine write_envelope_file
end program lunar_envelope
