!-------------------------------------------------------------------------------
! 月行九道 Phase 1: 从 DE440 批量提取月球地心状态矢量
! 模块: mod_lunar_extract
!
! 本模块提供了 extract_lunar_states 子程序，根据用户配置从
! JPL DE440 历表中按等间隔时间步长查询月球相对地球的状态矢量，
! 并将结果写入二进制文件与元数据文件。
!-------------------------------------------------------------------------------

module mod_lunar_extract
    use mod_precision, only: dp
    use mod_jpleph, only: DPLEPH
    use mod_time_system, only: date2jd, utc2tdb
    implicit none
    private

    public :: extract_lunar_states
    public :: lunar_extract_config

    ! 配置派生类型，集中管理所有可调参数
    type :: lunar_extract_config
        ! UTC 起止时间（日历）
        integer  :: yy_start = 1900, mm_start = 1, dd_start = 1
        integer  :: hh_start = 0,    mn_start = 0
        real(dp) :: ss_start = 0.0_dp
        integer  :: yy_end   = 2100, mm_end   = 1, dd_end   = 1
        integer  :: hh_end   = 0,    mn_end   = 0
        real(dp) :: ss_end   = 0.0_dp
        ! 采样步长（天）
        real(dp) :: step_days = 1.0_dp
        ! 输出文件名（不含路径）
        character(len=256) :: bin_filename  = "lunar_state_1900_2100.bin"
        character(len=256) :: meta_filename = "lunar_state_1900_2100.meta"
        ! 进度报告间隔（每 N 个采样点输出一次），<=0 表示静默
        integer :: progress_interval = 5000
    end type lunar_extract_config

contains

    !-------------------------------------------------------------------
    ! extract_lunar_states: 按配置批量提取月球地心状态矢量并写出文件
    !
    ! 参数:
    !   cfg         (in)  : 配置结构体
    !   output_dir  (in)  : 输出目录（不含末尾斜杠）
    !   info        (out) : 返回码
    !                       0 = 成功
    !                       1 = UTC->TDB 转换失败
    !                       2 = 起止历元逻辑错误（end <= start）
    !                       3 = 输出目录不可写 / I/O 错误
    !                       4 = DPLEPH 查询失败（历元超出 DE440 覆盖范围等）
    !                       5 = I/O 错误（写入文件时）
    !
    ! 注意:
    !   - 调用前需确保 init_physics_constants() 已被主程序调用
    !     （该函数会设置 DPLEPH 的单位为 km + km/s 并读取历表头）
    !   - DPLEPH 在遇到超出历表范围的时间时会直接 stop 程序，
    !     此时无法通过 info 返回码捕获，info=4 为文档性标记
    !   - 每个采样点的 TDB 历元独立计算，避免浮点累积误差
    !-------------------------------------------------------------------
    subroutine extract_lunar_states(cfg, output_dir, info)
        implicit none

        type(lunar_extract_config), intent(in)  :: cfg
        character(len=*),           intent(in)  :: output_dir
        integer,                    intent(out) :: info

        ! 时间系统
        real(dp) :: utc_start1, utc_start2
        real(dp) :: utc_end1,   utc_end2
        real(dp) :: tdb_start1, tdb_start2
        real(dp) :: tdb_end1,   tdb_end2
        real(dp) :: tdb_start, tdb_end
        integer  :: info_date, info_tdb

        ! 历表与循环
        real(dp) :: tdb1, tdb2
        real(dp) :: pv(6)
        integer  :: i, n_samples

        ! 文件 I/O
        integer, parameter :: FUNIT_BIN  = 20
        integer, parameter :: FUNIT_META = 21
        character(len=256) :: bin_path, meta_path
        integer  :: io_status

        ! 时间戳
        character(len=8)  :: date_str
        character(len=10) :: time_str

        ! ---- 初始化 ----
        info = 0

        ! ---- Step 1: UTC 日历 → UTC quasi-JD ----
        call date2jd(cfg%yy_start, cfg%mm_start, cfg%dd_start, &
                     cfg%hh_start, cfg%mn_start, cfg%ss_start, &
                     utc_start1, utc_start2, info_date)
        if (info_date < 0) then
            info = 1
            return
        end if

        call date2jd(cfg%yy_end, cfg%mm_end, cfg%dd_end, &
                     cfg%hh_end, cfg%mn_end, cfg%ss_end, &
                     utc_end1, utc_end2, info_date)
        if (info_date < 0) then
            info = 1
            return
        end if

        ! ---- Step 2: UTC quasi-JD → TDB ----
        call utc2tdb(utc_start1, utc_start2, tdb_start1, tdb_start2, info_tdb)
        if (info_tdb < 0) then
            info = 1
            return
        end if

        call utc2tdb(utc_end1, utc_end2, tdb_end1, tdb_end2, info_tdb)
        if (info_tdb < 0) then
            info = 1
            return
        end if

        ! ---- Step 3: 计算总采样数 ----
        tdb_start = tdb_start1 + tdb_start2
        tdb_end   = tdb_end1   + tdb_end2

        if (tdb_end <= tdb_start) then
            info = 2
            return
        end if

        n_samples = floor((tdb_end - tdb_start) / cfg%step_days) + 1

        ! ---- Step 4: 构造输出路径 ----
        bin_path  = trim(output_dir) // '/' // trim(cfg%bin_filename)
        meta_path = trim(output_dir) // '/' // trim(cfg%meta_filename)

        ! ---- Step 5: 打开二进制文件（stream access, unformatted）----
        open(unit=FUNIT_BIN, file=trim(bin_path), access='stream', &
             form='unformatted', status='replace', iostat=io_status)
        if (io_status /= 0) then
            info = 3
            return
        end if

        ! ---- Step 6: 打开元数据文件 ----
        open(unit=FUNIT_META, file=trim(meta_path), access='sequential', &
             form='formatted', status='replace', iostat=io_status)
        if (io_status /= 0) then
            close(FUNIT_BIN)
            info = 3
            return
        end if

        ! ---- Step 7: 获取当前时间戳 ----
        call date_and_time(date=date_str, time=time_str)

        ! ---- Step 8: 写入元数据头 ----
        write(FUNIT_META, '(a)') '# Lunar geocentric state vectors from DE440'
        write(FUNIT_META, '(a)') &
            '# Generated: ' // date_str(1:4) // '-' // date_str(5:6) // '-' // &
            date_str(7:8) // ' ' // time_str(1:2) // ':' // time_str(3:4) // ':' // time_str(5:10)
        write(FUNIT_META, '(a)') '# Frame: ICRF (DE440 native, no rotation)'
        write(FUNIT_META, '(a)') '# Center: Earth (NTARG_CENT = 3)'
        write(FUNIT_META, '(a)') '# Target: Moon (NTARG = 10)'
        write(FUNIT_META, '(a)') '# Units: km, km/s'
        write(FUNIT_META, '(a)') '# Time system: TDB'
        write(FUNIT_META, '(a)') '#'
        write(FUNIT_META, '(a, i4, a, i2.2, a, i2.2, a, i2.2, a, i2.2, a, f6.3, a)') &
            '#   Start (UTC): ', cfg%yy_start, '-', cfg%mm_start, '-', cfg%dd_start, &
            ' ', cfg%hh_start, ':', cfg%mn_start, ':', cfg%ss_start
        write(FUNIT_META, '(a, i4, a, i2.2, a, i2.2, a, i2.2, a, i2.2, a, f6.3, a)') &
            '#   End   (UTC): ', cfg%yy_end, '-', cfg%mm_end, '-', cfg%dd_end, &
            ' ', cfg%hh_end, ':', cfg%mn_end, ':', cfg%ss_end
        write(FUNIT_META, '(a, f10.6)') '#   Step (days): ', cfg%step_days
        write(FUNIT_META, '(a, i0)')    '#   Number of samples: ', n_samples
        write(FUNIT_META, '(a)') '#'
        write(FUNIT_META, '(a)') &
            '# Binary file layout (8 columns, real(dp), stream access):'
        write(FUNIT_META, '(a)') &
            '#   col 1: tdb1   [day]   TDB JD first component'
        write(FUNIT_META, '(a)') &
            '#   col 2: tdb2   [day]   TDB JD second component'
        write(FUNIT_META, '(a)') &
            '#   col 3: x      [km]    geocentric position X (ICRF)'
        write(FUNIT_META, '(a)') &
            '#   col 4: y      [km]    geocentric position Y (ICRF)'
        write(FUNIT_META, '(a)') &
            '#   col 5: z      [km]    geocentric position Z (ICRF)'
        write(FUNIT_META, '(a)') &
            '#   col 6: vx     [km/s]  geocentric velocity X (ICRF)'
        write(FUNIT_META, '(a)') &
            '#   col 7: vy     [km/s]  geocentric velocity Y (ICRF)'
        write(FUNIT_META, '(a)') &
            '#   col 8: vz     [km/s]  geocentric velocity Z (ICRF)'

        ! ---- Step 9: 主循环：逐采样点查询并写入 ----
        tdb1 = tdb_start1  ! TDB 第一分量保持固定

        write(*, '(a, i0, a)') 'Extracting ', n_samples, ' lunar state vectors...'

        do i = 0, n_samples - 1
            ! 独立计算 TDB 第二分量，避免浮点累积误差
            tdb2 = tdb_start2 + real(i, dp) * cfg%step_days

            ! 调用 DPLEPH: NTARG=10(月球), NCENT=3(地球)
            call DPLEPH([tdb1, tdb2], 10, 3, pv)

            ! 写入二进制文件
            write(FUNIT_BIN, iostat=io_status) &
                tdb1, tdb2, pv(1), pv(2), pv(3), pv(4), pv(5), pv(6)
            if (io_status /= 0) then
                info = 5
                close(FUNIT_BIN)
                close(FUNIT_META)
                return
            end if

            ! 进度报告
            if (cfg%progress_interval > 0) then
                if (mod(i + 1, cfg%progress_interval) == 0) then
                    write(*, '(a, i0, a, i0, a, f6.1, a)') &
                        '  Progress: ', i + 1, ' / ', n_samples, &
                        '  (', 100.0_dp * real(i + 1, dp) / real(n_samples, dp), '%)'
                end if
            end if
        end do

        ! ---- Step 10: 关闭文件 ----
        close(FUNIT_BIN)
        close(FUNIT_META)

        write(*, '(a)') 'Extraction complete.'
        write(*, '(a, a)') '  Binary: ', trim(bin_path)
        write(*, '(a, a)') '  Meta:   ', trim(meta_path)

    end subroutine extract_lunar_states

end module mod_lunar_extract
