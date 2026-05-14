!-------------------------------------------------------------------------------
! 月行九道 Phase 2: Radix-2 Cooley-Tukey 复数 FFT 模块
! 模块: mod_fft
!
! 实现按时间抽取 (DIT) 的基-2 Cooley-Tukey 快速傅里叶变换。
! 支持正向和逆向变换（通过 optional inverse 参数控制）。
! N 必须为 2 的幂。
!-------------------------------------------------------------------------------

module mod_fft
    use mod_precision, only: dp
    implicit none
    private

    public :: fft

contains

    !-------------------------------------------------------------------
    ! fft: 一维复数 FFT（DIT radix-2 Cooley-Tukey）
    !
    ! 参数:
    !   x       (inout) : 复数列向量，长度 N 必须为 2 的幂
    !                     输入为时域信号，输出为频谱
    !   info    (out)   : 返回码
    !                      0 = 成功
    !                      1 = N 不是 2 的幂
    !                      2 = N < 2
    !   inverse (in, optional) :
    !                     .true.  → 逆变换（含 1/N 缩放）
    !                     缺省     → 正变换
    !
    ! 算法:
    !   1. 检查输入合法性
    !   2. 复制到 0 基内部工作数组
    !   3. 位逆序置换
    !   4. 逐级蝶形运算（旋转因子即时计算）
    !   5. 逆变换时乘以 1/N
    !   6. 复制回 1 基输出数组
    !-------------------------------------------------------------------
    subroutine fft(x, info, inverse)
        complex(dp), intent(inout)        :: x(:)
        integer,     intent(out)          :: info
        logical,     intent(in), optional :: inverse

        integer  :: n, nbits, m, m_half, i, j, k
        real(dp) :: angle, pi
        complex(dp) :: w, t, u
        logical  :: inv

        ! 0-indexed work array
        complex(dp), allocatable :: work(:)

        n = size(x)

        ! ---- 检查输入 ----
        if (n < 2) then
            info = 2
            return
        end if

        if (.not. is_power_of_two(n)) then
            info = 1
            return
        end if

        info = 0

        ! ---- 确定变换方向 ----
        inv = .false.
        if (present(inverse)) inv = inverse

        ! ---- 复制到 0 基工作数组 ----
        allocate(work(0:n - 1))
        work(0:n - 1) = x(1:n)

        ! ---- 位逆序置换 ----
        nbits = log2_i(n)
        do i = 0, n - 1
            j = bit_reverse_i(i, nbits)
            if (j > i) then
                ! Swap work(i) and work(j) without temp variable
                t = work(i)
                work(i) = work(j)
                work(j) = t
            end if
        end do

        ! ---- Cooley-Tukey DIT 蝶形各 stage ----
        pi = 4.0_dp * atan(1.0_dp)
        m = 2
        do while (m <= n)
            m_half = m / 2

            do k = 0, m_half - 1
                ! 旋转因子: W = exp(∓2πi · k / m)
                ! 正变换用负角，逆变换用正角
                angle = -2.0_dp * pi * real(k, dp) / real(m, dp)
                if (inv) angle = -angle
                w = cmplx(cos(angle), sin(angle), kind=dp)

                ! 蝶形运算
                do i = k, n - 1, m
                    j = i + m_half
                    u = work(i)
                    t = w * work(j)
                    work(i) = u + t
                    work(j) = u - t
                end do
            end do

            m = m * 2
        end do

        ! ---- 逆变换缩放 ----
        if (inv) then
            work(0:n - 1) = work(0:n - 1) * (1.0_dp / real(n, dp))
        end if

        ! ---- 复制回 1 基输出数组 ----
        x(1:n) = work(0:n - 1)

    end subroutine fft

    !-------------------------------------------------------------------
    ! is_power_of_two: 判断正整数 n 是否为 2 的幂
    !-------------------------------------------------------------------
    pure logical function is_power_of_two(n) result(res)
        integer, intent(in) :: n
        res = (n > 0) .and. (iand(n, n - 1) == 0)
    end function is_power_of_two

    !-------------------------------------------------------------------
    ! bit_reverse_i: 将 x 的低 nbits 位按位逆序重排
    ! 例如: nbits=4, x=0b0011(3) → 0b1100(12)
    !-------------------------------------------------------------------
    pure integer function bit_reverse_i(x, nbits) result(res)
        integer, intent(in) :: x, nbits
        integer :: i
        res = 0
        do i = 0, nbits - 1
            if (btest(x, i)) then
                res = ibset(res, nbits - 1 - i)
            end if
        end do
    end function bit_reverse_i

    !-------------------------------------------------------------------
    ! log2_i: 计算整数以 2 为底的对数
    ! 输入必须为正整数，返回 floor(log2(n))
    !-------------------------------------------------------------------
    pure integer function log2_i(n) result(res)
        integer, intent(in) :: n
        integer :: val
        res = 0
        val = n
        do while (val > 1)
            val = ishft(val, -1)
            res = res + 1
        end do
    end function log2_i

end module mod_fft
