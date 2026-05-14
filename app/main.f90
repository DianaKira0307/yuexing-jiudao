!-------------------------------------------------------------------------------
! 月行九道 Phase 1: 主程序
!
! 调用 mod_lunar_extract 从 DE440 批量提取月球地心状态矢量
! 并输出到 data/ 目录下的二进制文件和元数据文件。
!-------------------------------------------------------------------------------
program extract_main
    use mod_precision, only: dp
    use mod_constants, only: init_physics_constants
    use mod_lunar_extract, only: lunar_extract_config, extract_lunar_states
    implicit none

    type(lunar_extract_config) :: cfg
    integer :: info

    ! 初始化物理常数与历表（设置 DPLEPH 单位为 km + km/s）
    call init_physics_constants()

    ! 使用默认配置（1900-01-01 至 2100-01-01，步长 1 天）
    ! 配置已在 lunar_extract_config 中设置了默认值

    ! 执行批量提取
    call extract_lunar_states(cfg, 'data', info)

    ! 根据返回码判定后续处理
    select case (info)
    case (0)
        write(*, '(a)') 'Program completed successfully.'
    case (1)
        write(*, '(a)') 'ERROR: UTC to TDB conversion failed.'
        stop 1
    case (2)
        write(*, '(a)') 'ERROR: End epoch <= start epoch.'
        stop 1
    case (3)
        write(*, '(a)') 'ERROR: Cannot open output directory or file.'
        stop 1
    case (4)
        write(*, '(a)') 'ERROR: DPLEPH query failed (epoch out of range).'
        stop 1
    case (5)
        write(*, '(a)') 'ERROR: I/O error during file writing.'
        stop 1
    case default
        write(*, '(a, i0)') 'ERROR: Unknown return code: ', info
        stop 1
    end select

end program extract_main
