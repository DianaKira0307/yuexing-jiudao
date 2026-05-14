!-------------------------------------------------------------------------------
!  本代码是轨道优化项目 E2M_Single_DSM 的一部分
!  是脚本 JPLEPH_Module.f90 的封装版本，提供 Fortran 90 Module 接口
!-------------------------------------------------------------------------------


! ==============================================================================
! JPLEPH_Module.f90
!
! JPL DE 系列行星历表读取库 —— Fortran 90 Module 封装版本
!
! 本文件在 JPLEPH_FreeForm.f90（忠实翻译版）的基础上做最小化改动：
!   1. 外层用 `module jpleph` 包裹，模块级变量替代 entry 共享状态
!   2. 原 PLEPH 的四个 entry 拆为四个独立的 public 子程序：
!        PLEPH / DPLEPH / CONST / PL_UNITS
!      （Module 内部不允许 entry 语句，故必须拆分）
!   3. READHD、INTCHB、FSIZER1/2/3 设为 private
!   4. 内部算法、变量名、循环结构、I/O 语句全部保持与 FreeForm 版一致
!
! 未改动内容：
!   - 所有插值数学逻辑
!   - Chebyshev 递推、地月质量分配、XSCALE/VSCALE 机制
!   - 文件头读取格式、记录缓存机制
!   - 原版 Disclaimer
!
! 使用方式：
!   ! 主程序通过 use 语句导入
!   program main
!     use jpleph
!     implicit none
!     ...
!     call PLEPH(TDB, NTARG, NCENT, PV)
!   end program
!
!   ! 编译（module 文件必须先编译）
!   gfortran -c JPLEPH_Module.f90
!   gfortran main.f90 JPLEPH_Module.o -o myprog
!
! 天体编号：     天体编号：
!       1=Mercury   2=Venus   3=Earth   4=Mars SSB   5=Jupiter SSB
!       6=Saturn SSB   7=Uranus SSB   8=Neptune SSB   9=Pluto SSB
!       10=Moon   11=Sun   12=Solar System Barycenter   13=Earth-Moon Barycenter
!       14=Nutations (经度+倾斜)   15=Lunar Euler angles (phi,theta,psi)
!       16=Lunar angular velocity (omegax,omegay,omegaz)   17=TT-TDB
! ==============================================================================

module mod_jpleph
  implicit none
  private

  ! ---- 公开接口 ----
  public :: PLEPH       ! 单精度时刻查询
  public :: DPLEPH      ! 双精度时刻查询
  public :: CONST       ! 获取历表常数
  public :: PL_UNITS    ! 切换输出单位

  ! ---- Module 级持久变量（替代原 entry 的共享 save 状态）----
  integer, parameter :: NMAX = 3000          ! 系数缓冲区最大字数

  character*6      :: NAMS(NMAX)             ! 历表常数名
  double precision :: VALS(NMAX)             ! 历表常数值
  double precision :: SS(3)                  ! 历表时间范围与步长
  double precision :: DATA(NMAX)             ! 当前记录的 Chebyshev 系数缓存
  integer          :: NCON                   ! 常数总数
  double precision :: FACTE, FACTM           ! 地月质量分配系数
  double precision :: XSCALE, VSCALE         ! 位置/速度缩放因子
  integer          :: NRFILE, NCOEFF         ! 文件单元号、每条记录双精度系数个数
  integer          :: IPT(3,15)              ! Chebyshev 指针表

  double precision :: SECSPAN                ! 记录时间跨度（秒）
  double precision, parameter :: SECDAY = 86400.D0

  ! 默认单位开关（PL_UNITS 可在首次查询前修改）
  logical :: AU_KM   = .true.                ! 默认输出 AU
  logical :: DAY_SEC = .true.                ! 默认速度单位 /day
  logical :: IAU_AU  = .false.               ! 默认使用历表自带 AU 值

  logical :: FIRST = .true.                  ! 首次调用标志
  integer :: NRREC = 0                       ! 上次读取的记录号（用于缓存）

contains


! ==============================================================================
! PLEPH —— 单精度时刻查询入口
!
! 输入：
!   TDB    —— 儒略历书时刻（double precision）
!   NTARG  —— 目标天体编号 (1..17)
!   NCENT  —— 中心天体编号 (1..13，对 NTARG≥14 被忽略)
!
! 输出：
!   PV(6)  —— 位置+速度（单位由 PL_UNITS 控制）
!             章动/天平动/角速度/TT-TDB 时含义见 FreeForm 版注释
! ==============================================================================
subroutine PLEPH( TDB, NTARG, NCENT, PV )

  implicit none

  double precision, intent(in)  :: TDB
  integer,          intent(in)  :: NTARG, NCENT
  double precision, intent(out) :: PV(6)

  double precision :: T1, T2

  T1 = TDB
  T2 = 0.D0

  call PLEPH_CORE( T1, T2, NTARG, NCENT, PV )

  return
end subroutine PLEPH


! ==============================================================================
! DPLEPH —— 双精度时刻查询入口
!
! 典型用法：TDB2(1) = 整数日，TDB2(2) = 小数日。
! 在 32 天历表块内精度可达 ~3 纳秒。
! ==============================================================================
subroutine DPLEPH( TDB2, NTARG, NCENT, PV )

  implicit none

  double precision, intent(in)  :: TDB2(2)
  integer,          intent(in)  :: NTARG, NCENT
  double precision, intent(out) :: PV(6)

  double precision :: T1, T2

  T1 = TDB2(1)
  T2 = TDB2(2)

  call PLEPH_CORE( T1, T2, NTARG, NCENT, PV )

  return
end subroutine DPLEPH


! ==============================================================================
! CONST —— 获取历表常数
!
! 无输入参数。首次调用会触发 READHD 读文件头。
!
! 输出：
!   NAMS1(*)  —— 常数名（character*6）
!   VALS1(*)  —— 常数值（double precision）
!   SS1(3)    —— SS1(1)=起始JED, SS1(2)=终止JED, SS1(3)=单条记录覆盖天数
!   NCON1     —— 常数总数
! ==============================================================================
subroutine CONST( NAMS1, VALS1, SS1, NCON1 )

  implicit none

  character*6,      intent(out) :: NAMS1(*)
  double precision, intent(out) :: VALS1(*)
  double precision, intent(out) :: SS1(*)
  integer,          intent(out) :: NCON1

  integer :: I

  ! 首次调用时触发初始化（借助一次虚拟 PLEPH 流程进入 READHD）
  if (FIRST) then
    call READHD                              ! 直接读头（module 内部可直接调用）
    FIRST = .false.
    SECSPAN = SECDAY * SS(3)
    do I = 1, NMAX
      DATA(I) = -999999999999999.d0
    enddo
  endif

  NCON1 = NCON
  do I = 1, NCON
    NAMS1(I) = NAMS(I)
    VALS1(I) = VALS(I)
  enddo
  do I = 1, 3
    SS1(I) = SS(I)
  enddo

  return
end subroutine CONST


! ==============================================================================
! PL_UNITS —— 覆盖默认输出单位
!
! 必须在首次调用 PLEPH/DPLEPH/CONST 之前调用，
! 否则后续的 XSCALE/VSCALE 已固化。
!
! AU_KM1   = .true.  → 位置单位 AU;  .false. → km
! DAY_SEC1 = .true.  → 速度单位 /day; .false. → /sec
! IAU_AU1  = .true.  → 使用 IAU 2012 AU (149597870.700 km)
!          = .false. → 使用历表自带 AU 值
! 默认 (.true., .true., .false.)
! ==============================================================================
subroutine PL_UNITS( AU_KM1, DAY_SEC1, IAU_AU1 )

  implicit none

  logical, intent(in) :: AU_KM1, DAY_SEC1, IAU_AU1

  AU_KM   = AU_KM1
  DAY_SEC = DAY_SEC1
  IAU_AU  = IAU_AU1

  return
end subroutine PL_UNITS


! ==============================================================================
! PLEPH_CORE —— 核心查询流程（private，不对外）
!
! 本子程序是原 testeph1.f 中 PLEPH 主体（合并了 PLEPH 和 DPLEPH 两个入口后
! 的公共处理段）的直接迁移。调用方把 TDB 拆为 T1+T2 后进入本例程。
! 算法与 JPLEPH_FreeForm.f90 中 PLEPH 的对应段一字不差。
! ==============================================================================
subroutine PLEPH_CORE( T1, T2, NTARG, NCENT, PV )

  implicit none

  double precision, intent(in)  :: T1, T2
  integer,          intent(in)  :: NTARG, NCENT
  double precision, intent(out) :: PV(6)

  double precision :: TF, TF1, TF2, SUMT
  double precision :: PV1(6), PV2(6), PVM(6), PVB(6)
  integer          :: I
  integer          :: IPVA
  data IPVA / 2 /                            ! INTCHB 模式：2=位置+速度

  ! 首次调用触发文件头读取
  if (FIRST) then
    call READHD
    if (NCOEFF .gt. NMAX) stop 'Coefficient array too small in PLEPH'
    FIRST = .false.
    SECSPAN = SECDAY * SS(3)
    do I = 1, NMAX
      DATA(I) = -999999999999999.d0
    enddo
  endif

!-----------------------------------------------------------------------
!     检验 NTARG 和 NCENT 的合法范围
!-----------------------------------------------------------------------

  if (NTARG .le. 0)  stop 'invalid NTARG < 0 in PLEPH'
  if (NTARG .gt. 17) stop 'invalid NTARG >17 in PLEPH'
  if (NCENT .le. 0) then
    if (NTARG .lt. 14) stop 'invalid NCENT < 0 in PLEPH'
  endif
  if (NCENT .gt. 13) then
    if (NTARG .lt. 14) stop 'invalid NCENT > 13 in PLEPH'
  endif

  if (NTARG .le. 13 .and. NCENT .ge. 14) then
    stop 'invalid NCENT >13 in for body PLEPH'
  endif

  ! 同天体对自身：返回零向量
  if (NTARG .le. 13 .and. NTARG .eq. NCENT) then
    do i = 1, 6
      PV(I) = 0.d0
    enddo
    return
  endif

!-----------------------------------------------------------------------
!     确定所需记录号，若非当前缓存记录则重新读盘
!-----------------------------------------------------------------------

  SUMT = T1 + T2

  if (SUMT .lt. SS(1)) stop 'input time before file start in PLEPH'
  if (SUMT .gt. SS(2)) stop 'input time after  file start in PLEPH'

  I = INT( (SUMT - SS(1)) / SS(3) ) + 3
  if (SUMT .eq. SS(2)) I = I - 1

  if (I .ne. NRREC) then
    NRREC = I
    read(NRFILE, rec=NRREC, err=99) (DATA(I), I=1, NCOEFF)
  endif

  ! 计算归一化时间 TF ∈ [0,1]
  TF1 = dble(int(T1))
  TF2 = T1 - TF1
  TF  = (TF1 - DATA(1)) / SS(3)
  TF  = TF + ((TF2 + T2) / SS(3))

!-----------------------------------------------------------------------
!     辅助量（章动/天平动/角速度/TT-TDB），不需要 NCENT
!-----------------------------------------------------------------------

!     章动 (NTARG=14)
  if (NTARG .eq. 14) then
    if (IPT(1,12) .gt. 0 .and. IPT(2,12) * IPT(3,12) .ne. 0) then
      call INTCHB( DATA(IPT(1,12)), TF, SECSPAN,     &
                   IPT(2,12), 2, IPT(3,12), IPVA, PV )
      PV(3) = SECDAY * PV(3)
      PV(4) = SECDAY * PV(4)
      PV(5) = -99.D99
      PV(6) = -99.D99
      return
    else
      write(6,*) 'Nutations requested but not found on file'
      do I = 1, 6
        PV(I) = -99.d99
      enddo
    endif
    return
  endif

!     月球天平动 Euler 角 (NTARG=15)
  if (NTARG .eq. 15) then
    if (IPT(1,13) .gt. 0 .and. IPT(2,13) * IPT(3,13) .ne. 0) then
      call INTCHB( DATA(IPT(1,13)), TF, SECSPAN,     &
                   IPT(2,13), 3, IPT(3,13), IPVA, PV )
      PV(4) = SECDAY * PV(4)
      PV(5) = SECDAY * PV(5)
      PV(6) = SECDAY * PV(6)
      return
    else
      write(6,*) 'Mantle Euler angles requested but not on file'
      do I = 1, 6
        PV(I) = -99.d99
      enddo
    endif
    return
  endif

!     月球角速度 (NTARG=16)
  if (NTARG .eq. 16) then
    if (IPT(1,14) .gt. 0 .and. IPT(2,14) * IPT(3,14) .ne. 0) then
      call INTCHB( DATA(IPT(1,14)), TF, SECSPAN,     &
                   IPT(2,14), 3, IPT(3,14), IPVA, PV )
      PV(4) = SECDAY * PV(4)
      PV(5) = SECDAY * PV(5)
      PV(6) = SECDAY * PV(6)
      return
    else
      write(6,*) 'Mantle angular velocity requested but not on file'
      do I = 1, 6
        PV(I) = -99.d99
      enddo
    endif
    return
  endif

!     TT-TDB (NTARG=17)
  if (NTARG .eq. 17) then
    if (IPT(1,15) .gt. 0 .and. IPT(2,15) * IPT(3,15) .ne. 0) then
      call INTCHB( DATA(IPT(1,15)), TF, SECSPAN,     &
                   IPT(2,15), 1, IPT(3,15), IPVA, PV )
      PV(2) = SECDAY * PV(2)
      do I = 3, 6
        PV(I) = -99.d99
      enddo
      return
    else
      write(6,*) 'TT-TDB requested but not found on file'
      do I = 1, 6
        PV(I) = -99.d99
      enddo
    endif
    return
  endif

!-----------------------------------------------------------------------
!     天体查询
!-----------------------------------------------------------------------

  ! 地球 ↔ 月球的快速通道
  if (NTARG .eq. 10 .and. NCENT .eq. 3) then
    call INTCHB( DATA(IPT(1,10)), TF, SECSPAN,       &
                 IPT(2,10), 3, IPT(3,10), IPVA, PVM )
    do I = 1, 3
      PV(i)   = PVM(i)   * XSCALE
      PV(i+3) = PVM(i+3) * VSCALE
    enddo
    return
  endif

  if (NTARG .eq. 3 .and. NCENT .eq. 10) then
    call INTCHB( DATA(IPT(1,10)), TF, SECSPAN,       &
                 IPT(2,10), 3, IPT(3,10), IPVA, PVM )
    do I = 1, 3
      PV(i)   = -PVM(i)   * XSCALE
      PV(i+3) = -PVM(i+3) * VSCALE
    enddo
    return
  endif

  do I = 1, 6
    PV1(I) = 0.d0
    PV2(I) = 0.d0
  enddo

  ! 涉及地球或月球的通用情况
  if (NTARG .eq. 3 .or. NTARG .eq. 10 .or.          &
      NCENT .eq. 3 .or. NCENT .eq. 10) then

    call INTCHB( DATA(IPT(1,10)),                   &
                 TF, SECSPAN, IPT(2,10), 3, IPT(3,10), IPVA, PVM )
    call INTCHB( DATA(IPT(1,3)),                    &
                 TF, SECSPAN, IPT(2,3),  3, IPT(3,3),  IPVA, PVB )

    if (NTARG .eq. 3 .or. NTARG .eq. 10) then

      if (NCENT .lt. 12) then
        call INTCHB( DATA(IPT(1,NCENT)), TF, SECSPAN, &
                     IPT(2,NCENT), 3, IPT(3,NCENT), IPVA, PV2 )
      else if (NCENT .eq. 13) then
        call INTCHB( DATA(IPT(1,3)), TF, SECSPAN,    &
                     IPT(2,3), 3, IPT(3,3), IPVA, PV2 )
      endif

      if (NTARG .eq. 3) then
        do I = 1, 6
          PV(I) = PVB(I) - FACTE * PVM(I) - PV2(I)
        enddo
      else
        do I = 1, 6
          PV(I) = PVB(I) - FACTM * PVM(I) - PV2(I)
        enddo
      endif

    else

      if (NTARG .lt. 12) then
        call INTCHB( DATA(IPT(1,NTARG)), TF, SECSPAN, &
                     IPT(2,NTARG), 3, IPT(3,NTARG), IPVA, PV1 )
      else if (NTARG .eq. 13) then
        call INTCHB( DATA(IPT(1,3)), TF, SECSPAN,    &
                     IPT(2,3), 3, IPT(3,3), IPVA, PV1 )
      endif

      if (NCENT .eq. 3) then
        do I = 1, 6
          PV(I) = PV1(I) - (PVB(I) - FACTE * PVM(I))
        enddo
      else
        do I = 1, 6
          PV(I) = PV1(I) - (PVB(I) - FACTM * PVM(I))
        enddo
      endif
    endif

    do I = 1, 3
      PV(i)   = PV(i)   * XSCALE
      PV(i+3) = PV(i+3) * VSCALE
    enddo
    return
  endif

!-----------------------------------------------------------------------
!     不涉及地球也不涉及月球的简化情况
!-----------------------------------------------------------------------

  if (NTARG .lt. 12) then
    call INTCHB( DATA(IPT(1,NTARG)), TF, SECSPAN,   &
                 IPT(2,NTARG), 3, IPT(3,NTARG), IPVA, PV1 )
  else if (NTARG .eq. 13) then
    call INTCHB( DATA(IPT(1,3)), TF, SECSPAN,       &
                 IPT(2,3), 3, IPT(3,3), IPVA, PV1 )
  endif

  if (NCENT .lt. 12) then
    call INTCHB( DATA(IPT(1,NCENT)), TF, SECSPAN,   &
                 IPT(2,NCENT), 3, IPT(3,NCENT), IPVA, PV2 )
  else if (NCENT .eq. 13) then
    call INTCHB( DATA(IPT(1,3)), TF, SECSPAN,       &
                 IPT(2,3), 3, IPT(3,3), IPVA, PV2 )
  endif

  do I = 1, 3
    PV(i)   = (PV1(I)   - PV2(I))   * XSCALE
    PV(i+3) = (PV1(I+3) - PV2(I+3)) * VSCALE
  enddo
  return

!-----------------------------------------------------------------------
!     读盘错误处理
!-----------------------------------------------------------------------

99 continue

  write(6,*) 'Error reading ephemeris data record in PLEPH'
  write(6,*) 'T1,T2 = ', T1, T2
  write(6,*) 'Record not found = ', NRREC
  stop

end subroutine PLEPH_CORE


! ==============================================================================
! READHD —— 读取历表文件头（private）
!
! 与 JPLEPH_FreeForm.f90 的 READHD 逻辑一致，
! 但改为无参数——所有输出直接写入 module 级变量。
!
!$ Disclaimer: 见 FreeForm 版文件头。
! ==============================================================================
subroutine READHD

  implicit none

  integer, parameter :: OLDMAX = 400

  integer          :: NRECL
  integer          :: KSIZE
  integer          :: NUMDE
  data NRECL / 0 /
  data NUMDE / 0 /
  character*6      :: TTL(14,3)
  character*80     :: NAMFIL

  double precision :: AU, EMRAT
  double precision :: iau
  data iau / 149597870.700d0 /               ! IAU 2012 天文单位 (km)

  integer :: I, J, K, L, IRECSZ

! ************************************************************************
!     FSIZER 版本选择（三选一）
! ************************************************************************

!       call FSIZER1( NRECL, KSIZE, NRFILE, NAMFIL )
        call FSIZER2( NRECL, KSIZE, NRFILE, NAMFIL )
!       call FSIZER3( NRECL, KSIZE, NRFILE, NAMFIL )

! ************************************************************************

  if (NRECL .EQ. 0) write(*,*) '  ***** FSIZER IS NOT WORKING *****'

  IRECSZ = NRECL * KSIZE
  NCOEFF = KSIZE / 2

  open( NRFILE,                    &
        FILE   = NAMFIL,           &
        ACCESS = 'DIRECT',         &
        FORM   = 'UNFORMATTED',    &
        RECL   = IRECSZ,           &
        STATUS = 'OLD' )

  read(NRFILE, REC=1) TTL, (NAMS(K), K=1, OLDMAX), SS, NCON, AU, EMRAT

  if (NCON .le. OLDMAX) then

    read(NRFILE, REC=1) TTL, (NAMS(K), K=1, OLDMAX), SS, NCON, AU, EMRAT,  &
         ((IPT(I,J), I=1,3), J=1,12), NUMDE, (IPT(I,13), I=1,3),           &
                                             (IPT(I,14), I=1,3),           &
                                             (IPT(I,15), I=1,3)

    read(NRFILE, REC=2) (VALS(I), I=1, OLDMAX)

  else

    if (NCON .gt. NMAX) then
      write(*,*) 'Number of ephemeris constants too big in READHD'
      stop
    endif

  read(NRFILE, REC=1) TTL, (NAMS(K), K=1, OLDMAX), SS, NCON, AU, EMRAT,  &
         ((IPT(I,J), I=1,3), J=1,12), NUMDE, (IPT(I,13), I=1,3),           &
         (NAMS(L), L=OLDMAX+1, NCON),                                      &
                                             (IPT(I,14), I=1,3),           &
                                             (IPT(I,15), I=1,3)

    read(NRFILE, REC=2) (VALS(I), I=1, NCON)

  endif

  FACTE = 1.D0 / (1.D0 + EMRAT)
  FACTM = FACTE - 1.D0

  if (FACTE .eq. 0.D0) then
    write(*,*) 'Invalid value of EMRAT from file in READHD'
    stop
  endif

  if (AU_KM) then
    if (IAU_AU) then
      XSCALE = 1.d0 / IAU
    else
      XSCALE = 1.d0 / AU
    endif
  else
    XSCALE = 1.d0
  endif

  if (DAY_SEC) then
    VSCALE = XSCALE * 86400.d0
  else
    VSCALE = XSCALE
  endif

  if (NUMDE .eq. 0) stop 'DENUM not found by READHD in constants'

  write(6,*)
  write(6,'(a27,i3.3)') ' JPL planetary ephemeris DE', NUMDE
  write(6,*) 'Requested output units are :'
  if (AU_KM) then
    if (IAU_AU) then
      write(6,*) 'IAU au for distance'
      if (DAY_SEC) then
        write(6,*) 'IAU au/day for velocity'
      else
        write(6,*) 'IAU au/sec for velocity'
      endif
    else
      write(6,'(a2,i3.3,a16)') 'DE', NUMDE, ' au for distance'
      if (DAY_SEC) then
        write(6,'(a2,i3.3,a20)') 'DE', NUMDE, ' au/day for velocity'
      else
        write(6,'(a2,i3.3,a20)') 'DE', NUMDE, ' au/sec for velocity'
      endif
    endif
  else
    write(6,*) 'km for distance'
    if (DAY_SEC) then
      write(6,*) 'km/day for velocity'
    else
      write(6,*) 'km/sec for velocity'
    endif
  endif

  return
end subroutine READHD


! ==============================================================================
! FSIZER1 (private) —— VAX 版
! 见 FreeForm 版对应注释。
! ==============================================================================
subroutine FSIZER1( NRECL, KSIZE, NRFILE, NAMFIL )

  implicit none

  integer,      intent(out) :: NRECL, KSIZE
  integer,      intent(out) :: NRFILE
  character*80, intent(out) :: NAMFIL

  integer :: IRECSZ

!     用户参数（NRECL、NRFILE、NAMFIL）需根据环境设置
!     NRECL =
!     NRFILE = 12
!     NAMFIL = 'JPLEPH'

  IRECSZ = 0
  inquire(FILE=NAMFIL, RECL=IRECSZ)

  if (IRECSZ .LE. 0) then
    write(*,*) ' INQUIRE STATEMENT PROBABLY DID NOT WORK'
  endif

  KSIZE = IRECSZ / NRECL

  return
end subroutine FSIZER1


! ==============================================================================
! FSIZER2 (private) —— UNIX 版（当前激活）
!
! 以临时大 RECL 打开文件，读 IPT 指针表后计算真实 KSIZE，再关闭。
! 调用方（READHD）随后以正确 RECL 重新打开。
! ==============================================================================
subroutine FSIZER2( NRECL, KSIZE, NRFILE, NAMFIL )

  implicit none

  integer,      intent(out) :: NRECL, KSIZE
  integer,      intent(out) :: NRFILE
  character*80, intent(out) :: NAMFIL

  double precision :: SS_LOCAL(3), AU, EMRAT

  integer, parameter :: OLDMAX = 400
  integer, parameter :: NMAX_LOCAL = 1000

  integer :: IPT_LOCAL(3,15)
  data IPT_LOCAL / 45*0 /

  integer :: I, J, K, L, ND, KHI, KMX, MRECL, NCON_LOCAL, NUMDE

  character*6 :: TTL(14,3), CNAM(NMAX_LOCAL)

  NRECL  = 4
  NRFILE = 12
  NAMFIL = 'JPLEPH-DE440-1Bytes'

  MRECL = NRECL * 10000

  open( NRFILE,                     &
        FILE   = NAMFIL,            &
        ACCESS = 'DIRECT',          &
        FORM   = 'UNFORMATTED',     &
        RECL   = MRECL,             &
        STATUS = 'OLD' )

  read(NRFILE, REC=1) TTL, (CNAM(K), K=1, OLDMAX), SS_LOCAL, NCON_LOCAL

  if (NCON_LOCAL .le. OLDMAX) then

    read(NRFILE, REC=1) TTL, (CNAM(K), K=1, OLDMAX),                       &
         SS_LOCAL, NCON_LOCAL, AU, EMRAT,                                  &
         ((IPT_LOCAL(I,J), I=1,3), J=1,12), NUMDE,                         &
         (IPT_LOCAL(I,13), I=1,3),                                         &
         (IPT_LOCAL(I,14), I=1,3),                                         &
         (IPT_LOCAL(I,15), I=1,3)

  else

    if (NCON_LOCAL .gt. NMAX_LOCAL) then
      write(*,*) 'Number of ephemeris constants too big for FSIZER'
      stop
    endif

    read(NRFILE, REC=1) TTL, (CNAM(K), K=1, OLDMAX),                       &
         SS_LOCAL, NCON_LOCAL, AU, EMRAT,                                  &
         ((IPT_LOCAL(I,J), I=1,3), J=1,12), NUMDE,                         &
         (IPT_LOCAL(I,13), I=1,3),                                         &
         (CNAM(L), L=OLDMAX+1, NCON_LOCAL),                                &
         (IPT_LOCAL(I,14), I=1,3),                                         &
         (IPT_LOCAL(I,15), I=1,3)

  endif

  close(NRFILE)

  KMX = 0
  KHI = 0

  do I = 1, 15
    if (IPT_LOCAL(1,I) .ge. KMX) then
      KMX = IPT_LOCAL(1,I)
      KHI = I
    endif
  enddo

  ND = 3
  if (KHI .EQ. 12) ND = 2
  if (KHI .EQ. 15) ND = 1

  KSIZE = IPT_LOCAL(1,KHI) - 1 + ND * IPT_LOCAL(2,KHI) * IPT_LOCAL(3,KHI)
  KSIZE = 2 * KSIZE

  return
end subroutine FSIZER2


! ==============================================================================
! FSIZER3 (private) —— 硬编码版
! 需要用户手动设置 NRECL 和 KSIZE（见 FreeForm 版注释表）。
! ==============================================================================
subroutine FSIZER3( NRECL, KSIZE, NRFILE, NAMFIL )

  implicit none

  integer,      intent(out) :: NRECL, KSIZE
  integer,      intent(out) :: NRFILE
  character*80, intent(out) :: NAMFIL

!     NRECL =       ! 用户填写
  NRFILE = 12
  NAMFIL = 'JPLEPH'
!     KSIZE =       ! 用户填写（DE430/DE440 → 2036）

  return
end subroutine FSIZER3


! ==============================================================================
! INTCHB (private) —— Chebyshev 多项式插值
! 与 JPLEPH_FreeForm.f90 的 INTCHB 代码完全一致。
!
!$ Disclaimer: 见 FreeForm 版文件头。
! ==============================================================================
subroutine INTCHB( BUF, T, LINT, NCF, NCM, NSC, REQ, PV )

  implicit none

  integer          :: NCF, NCM, NSC, REQ
  double precision :: BUF(NCF,NCM,NSC), T, LINT
  double precision :: PV(*)

  double precision :: PC(18), VC(2:18), TEMP, TC, TTC, BMA
  integer          :: NS, L, NP, NV, NC, I, J

  data PC(1) / 1.D0 /, PC(2) / 2.D0 /, VC(2) / 1.D0 /

  save

  NS   = NSC
  TEMP = T * dble(NS)
  TC   = 2.d0 * (TEMP - DINT(TEMP)) - 1.d0
  L    = TEMP + 1
  if (L .gt. NS) then
    L  = L - 1
    TC = TC + 2.d0
  endif

  if (TC .ne. PC(2)) then
    NP = 3
    NV = 4
    PC(2) = TC
    TTC   = TC + TC
    VC(3) = TTC + TTC
  endif

  NC = NCF
  do NP = NP, NC
    PC(NP) = TTC * PC(NP-1) - PC(NP-2)
  enddo

  do I = 1, NCM
    TEMP = 0.D0
    do J = NC, 1, -1
      TEMP = TEMP + PC(J) * BUF(J,I,L)
    enddo
    PV(I) = TEMP
  enddo

  IF (REQ .LE. 1) RETURN

  do NV = NV, NC
    VC(NV) = TTC * VC(NV-1) + PC(NV-1) + PC(NV-1) - VC(NV-2)
  enddo

  BMA = DBLE(2*NS) / LINT
  do I = 1, NCM
    TEMP = 0.D0
    do J = NC, 2, -1
      TEMP = TEMP + VC(J) * BUF(J,I,L)
    enddo
    PV(I+NCM) = TEMP * BMA
  enddo

  return
end subroutine INTCHB


end module mod_jpleph
