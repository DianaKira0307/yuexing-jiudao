!-------------------------------------------------------------------------------
!  本代码是轨道优化项目 E2M_Single_DSM 的一部分
!  用于存储项目的底层工具函数
!-------------------------------------------------------------------------------

!===============================================================================
!  模块: mod_precision
!  项目中采用双精度数
!===============================================================================
module mod_precision
    implicit none 
    integer, parameter :: dp = selected_real_kind(15, 307)  ! 双精度
end module mod_precision

!===============================================================================
!  模块: mod_constants
!  包含: 
!  1. init_physics_constants: 从 JPLEPH 读取历表 header，填充物理常数
!  2. get_ephemeris_constant: 通用历表常数查询接口（按名字）
!
!  不包含了星历读取程序中的 CONST 函数
!  这里定义了项目中使用的转换因子，如角度与弧度之间的转换 
!===============================================================================
module mod_constants
    use mod_precision, only : dp
    use mod_jpleph, only : CONST, PL_UNITS
    implicit none
    private

    ! 数学常数
    real(dp), parameter, public :: PI = 3.1415926535897932384626433832795_dp  ! 圆周率
    real(dp), parameter, public :: PI2 = 2.0_dp * PI  ! 2π
    real(dp), parameter, public :: RAD2DEG = 180.0_dp / PI  ! 弧度转角度的转换因子
    real(dp), parameter, public :: DEG2RAD = PI / 180.0_dp  ! 角度转弧度的转换因子


    ! 物理常数
    ! 运行期初始化，调用 init_physics_constants 后填充
    real(dp), protected, public :: AU_km   = 0.0_dp  ! 1 天文单位的 km 值
    real(dp), protected, public :: TU_s    = 0.0_dp  ! 时间归一化单位 Time Unit (s)
    real(dp), protected, public :: VU_km_s = 0.0_dp  ! 速度归一化单位 Velocity Unit (km/s)

    !-------------------------------------------------------------------
    ! 天体引力常数数组 GM(i)，单位 km^3/s^2
    ! 下标 i 对应 DE 历表天体编号:
    !   1: 水星    2: 金星     3: 地球    4: 火星    5: 木星
    !   6: 土星    7: 天王星   8: 海王星  9: 冥王星  10: 月球
    !   11: 太阳
    !-------------------------------------------------------------------
    real(dp), protected, public :: GM(11) = 0.0_dp


    ! 内部缓存（存储 CONST 返回的所有历表常数，供 get_ephemeris_constant 查询）
    logical, private :: initialized = .false.
    character(len=6), allocatable, private :: cached_names(:)
    real(dp),         allocatable, private :: cached_values(:)
    integer,                       private :: cached_ncon = 0

    public :: init_physics_constants
    public :: get_ephemeris_constant

contains

    !-------------------------------------------------------------------
    ! init_physics_constants: 从 JPLEPH 读取历表 header，填充物理常数
    ! 
    ! 本子程序需要在首次使用星历或物理常数之前调用（推荐在主程序开头）
    ! 重复调用是安全的（幂等）
    !
    ! 功能:
    !   1. 调用 PL_UNITS 设置 PLEPH 输出单位为 km + km/s
    !   2. 调用 CONST 读取 DE 历表 header 中的全部常数
    !   3. 从下标提取 AU、各天体 GM，转换到 km + s 单位
    !   4. 计算派生归一化因子 TU_s, VU_km_s
    !   5. 缓存原始 names/values 数组供后续通用查询
    !
    ! 下标约定（基于老代码，对 DE430/DE432/DE435 等版本通用）:
    !   VAL(10) = AU (km)
    !   VAL(11) = EMRAT（地月质量比）
    !   VAL(12) ~ VAL(20) = GM1 ~ GM9 (水星~冥王星，AU^3/day^2)
    !     （注意: VAL(14) 实际是 GMB，即地月质心总 GM，需用 EMRAT 分解）
    !   VAL(21) = GMS (太阳, AU^3/day^2)
    !-------------------------------------------------------------------
    subroutine init_physics_constants()
        character(len=6) :: nam(1000)
        real(dp) :: val(1000), sss(3), unit ! sss 是 CONST 的输出参数，即历表起始和终止的 JED 和步长，此处不用
        integer :: nn

        if (initialized) return

        ! 设定 PLEPH 输出单位: km + km/s
        ! 分别是 km, /sec, AU 值采用历表定义
        call PL_UNITS(.false., .false., .false.)

        ! 读取历表 header 中全部常数
        call CONST(nam, val, sss, nn)

        ! 提取 1 AU 的 km 值（历表里已经是 km，不需要再乘 1000）
        AU_km = val(10)

        ! 单位换算因子（AU^3/day^2 → km^3/s^2）
        unit = AU_km**3 / 86400.0_dp**2

        ! 行星 GM（水星~冥王星，直接对应 val(12)~val(20), 除地球外）
        ! 此部分的 val 序号在 DE430/DE440 中都是固定的
        ! DE440 的常量定义在 22号及以后和 DE430 不同
        ! 常量定义参考文件内的 de440_constants.txt 中的注释
        GM(1) = val(12) * unit  ! 水星
        GM(2) = val(13) * unit  ! 金星
        GM(4) = val(15) * unit  ! 火星
        GM(5) = val(16) * unit  ! 木星
        GM(6) = val(17) * unit  ! 土星
        GM(7) = val(18) * unit  ! 天王星
        GM(8) = val(19) * unit  ! 海王星
        GM(9) = val(20) * unit  ! 冥王星

        ! 地球与月球 GM（由地月质心 GMB 和质量比 EMRAT 分解）
        ! 这部分常量在星历中没有直接给出，需要用以下方法进行计算：
        !   GMB = GM_earth + GM_moon,  EMRAT = M_earth / M_moon
        !   → GM_earth = GMB * EMRAT/(EMRAT+1)
        !   → GM_moon  = GMB * 1/(EMRAT+1)
        GM(3)  = val(14) * val(11) / (val(11) + 1.0_dp) * unit  ! 地球
        GM(10) = val(14)           / (val(11) + 1.0_dp) * unit  ! 月球

        ! 太阳 GM
        GM(11) = val(21) * unit

        ! 派生归一化因子（基于日心引力常数 GM(11)）
        TU_s    = sqrt(AU_km**3 / GM(11))  ! 时间单位 (s)
        VU_km_s = AU_km / TU_s             ! 速度单位 (km/s)

        ! 缓存原始 names/values 供通用查询接口使用
        cached_ncon = nn
        allocate(cached_names(nn), cached_values(nn))
        cached_names (1:nn) = nam(1:nn)
        cached_values(1:nn) = val(1:nn)

        initialized = .true.
    end subroutine init_physics_constants

    !-------------------------------------------------------------------
    ! get_ephemeris_constant: 通用历表常数查询接口（按名字）
    !
    ! 从 DE 历表 header 中取任意一个常数。需要在调用前已完成
    ! init_physics_constants 初始化，否则将自动触发初始化。
    !
    ! 注意事项:
    !   - NAME 必须用大写（DE 历表约定），例如 'AU', 'GMS', 'EMRAT'
    !   - NAME 长度不足 6 字符时会自动右侧补空格后与缓存名对比
    !   - 未找到时 value 返回 0，found 返回 .false.
    !
    ! 参数:
    !   NAME(in)  - 常数名称（大写字符串，长度 ≤ 6）
    !   value(out) - 返回常数值（历表原始单位，未做单位转换）
    !   found(optional) - 可选，查询成功与否
    !-------------------------------------------------------------------
    subroutine get_ephemeris_constant(NAME, value, found)
        character(len=*), intent(in)  :: NAME
        real(dp),         intent(out) :: value
        logical, optional, intent(out) :: found

        integer :: i
        logical :: ok

        ok = .false.
        value = 0.0_dp

        ! 若未初始化则自动初始化
        if (.not. initialized) call init_physics_constants()

        ! 在缓存里做线性查找
        do i = 1, cached_ncon
            if (trim(cached_names(i)) == trim(NAME)) then
                value = cached_values(i)
                ok = .true.
                exit
            end if
        end do

        ! 查询失败且用户未传 found 时发出警告
        if (.not. ok .and. .not. present(found)) then
            write(*, '(3a)') 'Warning [get_ephemeris_constant]: constant "', &
                            trim(NAME), '" not found.'
        end if

        if (present(found)) found = ok
    end subroutine get_ephemeris_constant

end module mod_constants