!-------------------------------------------------------------------------------
!  本代码是轨道优化项目 E2M_Single_DSM 的一部分
!  基于 IAU SOFA 库实现了时间系统转换和坐标系统转换的工具函数
!-------------------------------------------------------------------------------

!===========================================================================
! mod_time_system: 时间系统模块
!===========================================================================
module mod_time_system
    use mod_precision, only : dp
    implicit none 

contains 
    !------------------------------------------------------------------
    ! 1. utc2tdb : 将 UTC 时间转换为 TDB 时间
    ! 输入:
    !   utc1, utc2 : UTC 时间。严格来说因为闰秒的存在，UTC 时间对应的 Julian Date 并不能作为良好的时间
    !                间隔定义，因此被称为 quasi-JD。
    ! 输出:
    !   tdb1, tdb2 : 双精度存储的 TDB 时间，分别表示整数部分和小数部分，单位: JD
    !   info       : 状态码，0 表示成功，非零表示失败
    ! 注意:
    !   - 该函数内部调用了 IAU SOFA 库的时间系统
    !   - 忽略了 TT 和 TDB 之间的差异，这在大多数轨道优化问题中是可以接受的
    !------------------------------------------------------------------
    subroutine utc2tdb(utc1, utc2, tdb1, tdb2, info)
        implicit none 

        ! 声明输入输出变量
        real(dp), intent(in) :: utc1, utc2  ! 双精度存储的 UTC 时间
        real(dp), intent(out) :: tdb1, tdb2  ! 双精度存储的 TDB 时间
        integer, intent(out) :: info  ! 状态码

        ! 声明局部变量
        real(dp) :: tai1, tai2  ! 双精度存储的 TAI 时间
        real(dp) :: tt1, tt2   ! 双精度存储的 TT 时间
        integer :: info_UTCTAI, info_TAITT, info_TTTDB ! 各个时间系统转换的状态码

        ! 初始化变量
        info = -1

        ! UTC 转 TAI
        call iau_UTCTAI(utc1, utc2, tai1, tai2, info_UTCTAI)

        ! TAI 转 TT
        call iau_TAITT(tai1, tai2, tt1, tt2, info_TAITT)

        ! TT 转 TDB
        ! 此处假设 TT 和 TDB 之间的差异可以忽略，直接将 TT 时间赋值给 TDB 时间
        call iau_TTTDB(tt1, tt2, 0.0_dp, tdb1, tdb2, info_TTTDB)

        ! 综合状态码
        info = info_UTCTAI + info_TAITT + info_TTTDB ! 只要有一个转换失败，info 就会是非零的

    end subroutine utc2tdb

    !------------------------------------------------------------------
    ! 2. date2jd : 将年、月、日、时、分、秒转换为 UTC 时间的 quasi-JD 表示
    ! 输入:
    !   yy, mm, dd, hh, mn : 年、月、日、时、分
    !   ss : 秒，允许有小数部分
    ! 输出:
    !   utc1, utc2 : 双精度存储的 UTC 时间，分别表示整数部分和小数部分，单位: JD
    !   info : 状态码
    !               +3 = both of next two
    !               +2 = time is after end of day (Note 5)
    !               +1 = dubious year (Note 6)
    !               0 = OK
    !               −1 = bad year
    !               −2 = bad month
    !               −3 = bad day (Note 3)
    !               −4 = bad minute
    !               −5 = bad second (<0)
    ! 注意:
    !   - 该函数内部调用了 IAU SOFA 库的时间系统
    !   - 由于闰秒的存在，UTC 时间对应的 Julian Date 并不能作为良好的时间间隔定义，因此被称为 quasi-JD
    !------------------------------------------------------------------
    subroutine  date2jd(yy, mm, dd, hh, mn, ss, utc1, utc2, info)
        implicit none 

        ! 声明输入输出变量
        integer, intent(in) :: yy, mm, dd, hh, mn ! 年月日时分
        real(dp), intent(in) :: ss  ! 秒
        real(dp), intent(out) :: utc1, utc2  ! 双精度存储的 UTC quasi-JD
        integer, intent(out) :: info  ! 状态码

        ! 声明局部变量
        character(len = 3) :: scale ! 时间系统标识符

        ! 设定转换所需的时间系统标识符，IAU SOFA 库使用 "UTC" 来表示协调世界时
        scale = "UTC"

        ! Date 编码为 quasi-JD
        call iau_DTF2D(scale, yy, mm, dd, hh, mn, ss, utc1, utc2, info)

    end subroutine date2jd

end module mod_time_system