module philox_rng
    use iso_fortran_env, only: int32, int64, real64
    implicit none

    integer(int32), parameter :: PHILOX_M0 = int(z'D2511F53', int32)
    integer(int32), parameter :: PHILOX_M1 = int(z'CD9E8D57', int32)
    integer(int32), parameter :: PHILOX_W0 = int(z'9E3779B9', int32) ! Weyl constant 1
    integer(int32), parameter :: PHILOX_W1 = int(z'BB67AE85', int32) ! Weyl constant 2

    integer(int32), parameter :: PHILOX_KEY0 = int(z'6D2B79F5', int32)  ! Randomly chosen by ChatGPT
    integer(int32), parameter :: PHILOX_KEY1 = int(z'1B873593', int32)  ! Randomly chosen by ChatGPT

    contains
    pure subroutine mulhilo32( a, b, hi, lo )
        implicit none
        integer(int32), intent(in)  :: a, b
        integer(int32), intent(out) :: hi, lo
        integer(int64)              :: aa, bb, prod

        aa = iand( int(a, int64), int(z'FFFFFFFF', int64) )
        bb = iand( int(b, int64), int(z'FFFFFFFF', int64) )
        prod = aa * bb

        lo = int( iand( prod, int(z'FFFFFFFF', int64) ), int32 )
        hi = int( iand( ishft(prod, -32), int(z'FFFFFFFF', int64) ), int32 )
    end subroutine mulhilo32

    pure function philox_real( word ) result(u)
        implicit none
        integer(int32), intent(in) :: word
        real(real64)               :: u
        integer(int64)             :: uw

        uw = iand( int(word, int64), int(z'FFFFFFFF', int64) )
        u = real(uw, real64) * ( 1.0_real64 / 4294967296.0_real64 )  ! / 2^32
    end function philox_real

    ! --- core Philox4x32, 10 rounds ---
    ! ctr(4): 128-bit counter as four 32-bit words (ctr0..ctr3)
    ! key(2): 64-bit key as two 32-bit words (key0, key1)
    ! out(4): the four 32-bit pseudorandom outputs
    pure subroutine philox_uniform( ctr, out )
        implicit none
        integer(int32), intent(in)  :: ctr(4)
        real(real64),   intent(out) :: out(4)

        integer(int32) :: c0, c1, c2, c3
        integer(int32) :: k0, k1
        integer(int32) :: hi0, lo0, hi1, lo1
        integer        :: round

        c0 = ctr(1); c1 = ctr(2); c2 = ctr(3); c3 = ctr(4)
        k0 = PHILOX_KEY0; k1 = PHILOX_KEY1

        do round = 1, 10
            call mulhilo32( PHILOX_M0, c0, hi0, lo0 )
            call mulhilo32( PHILOX_M1, c2, hi1, lo1 )

            c0 = ieor( hi1, ieor( c1, k0 ) )
            c1 = lo1
            c2 = ieor( hi0, ieor( c3, k1 ) )
            c3 = lo0

            ! Weyl sequence key bump (skip on final round, matches reference)
            if (round < 10) then
                k0 = k0 + PHILOX_W0
                k1 = k1 + PHILOX_W1
            endif
        enddo

        out(1) = philox_real( c0 )
        out(2) = philox_real( c1 )
        out(3) = philox_real( c2 )
        out(4) = philox_real( c3 )

    end subroutine philox_uniform

    pure subroutine philox_normal_3seq( ctr1, ctr2, ctr3, ctr4, n1, n2, n3 )
        implicit none
        integer(int32), intent(in)  :: ctr1, ctr2, ctr3, ctr4
        real(real64), intent(out)   :: n1, n2, n3

        real(real64), parameter       :: TWOPI = 6.283185307179586
        real(real64), parameter       :: TINY  = 1.0e-35
        integer(int32)                :: ctr(4), key(2)
        real(real64)                  :: out(4)
        real(real64)                  :: u1, u2, u3, u4, n4
        real(real64)                  :: r1, r2, th1, th2

        ! counter uniquely identifies WHICH random draw this is —
        ! no shared mutable state, fully deterministic given these 4 integers
        ctr(1) = ctr1; ctr(2) = ctr2; ctr(3) = ctr3; ctr(4) = ctr4

        call philox_uniform( ctr, out )
        u1 = out(1); u2 = out(2); u3 = out(3); u4 = out(4)

        u1 = max( u1, TINY )   ! guard log(0)
        u3 = max( u3, TINY )   ! guard log(0)

        ! Two Box-Muller transform for Gaussian variates
        r1  = sqrt( -2.0 * log(u1) )
        th1 = TWOPI * u2

        r2  = sqrt( -2.0 * log(u3) )
        th2 = TWOPI * u4

        n1 = r1 * cos(th1)
        n2 = r1 * sin(th1)

        n3 = r2 * cos(th2)
        n4 = r2 * sin(th2)

    end subroutine philox_normal_3seq

end module philox_rng

module xmath
    use iso_fortran_env, only: real64
    implicit none
    real(real64), parameter :: pi = 3.1415926535897932384626433832795

    interface operator(.x.)
        module procedure vmdot
        module procedure mvdot
        module procedure mmdot
    end interface

    contains

    pure function cross(a, b) result (out)
        implicit none
        real(real64), intent(in), dimension(3) :: a, b
        real(real64), dimension(3)             :: out

        out(1) = a(2)*b(3) - a(3)*b(2)
        out(2) = a(3)*b(1) - a(1)*b(3)
        out(3) = a(1)*b(2) - a(2)*b(1)

    end function cross

    pure function determinant(M) result (out)
        implicit none
        real(real64), intent(in), dimension(3,3) :: M
        real(real64)                             :: out

        out = M(1, 1) * ( M(2, 2) * M(3, 3) - M(3, 2) * M(2, 3) ) &
            - M(1, 2) * ( M(2, 1) * M(3, 3) - M(3, 1) * M(2, 3) ) &
            + M(1, 3) * ( M(2, 1) * M(3, 2) - M(3, 1) * M(2, 2) )

    end function determinant

    pure function inverse(M) result (out)
        implicit none
        real(real64), intent(in), dimension(3,3) :: M
        real(real64), dimension(3,3)             :: out
        real(real64)                             :: det

        out(1, 1) = M(2, 2) * M(3, 3) - M(3, 2) * M(2, 3)
        out(1, 2) = M(1, 3) * M(3, 2) - M(1, 2) * M(3, 3)
        out(1, 3) = M(1, 2) * M(2, 3) - M(1, 3) * M(2, 2)

        out(2, 1) = M(2, 3) * M(3, 1) - M(2, 1) * M(3, 3)
        out(2, 2) = M(1, 1) * M(3, 3) - M(1, 3) * M(3, 1)    
        out(2, 3) = M(1, 3) * M(2, 1) - M(1, 1) * M(2, 3)

        out(3, 1) = M(2, 1) * M(3, 2) - M(2, 2) * M(3, 1)
        out(3, 2) = M(1, 2) * M(3, 1) - M(1, 1) * M(3, 2)
        out(3, 3) = M(1, 1) * M(2, 2) - M(1, 2) * M(2, 1)

        det = determinant(M)
        out = out / det

    end function inverse

    pure function as_matrix( w, x, y, z ) result (out)
        implicit none
        real(real64), intent(in)     :: w, x, y, z
        real(real64), dimension(3,3) :: out

        out(1, 1) = 1 - 2 * y ** 2 - 2 * z ** 2
        out(1, 2) = 2 * x * y - 2 * w * z
        out(1, 3) = 2 * w * y + 2 * x * z

        out(2, 1) = 2 * x * y + 2 * w * z
        out(2, 2) = 1 - 2 * x ** 2 - 2 * z ** 2
        out(2, 3) = 2 * y * z - 2 * w * x

        out(3, 1) = 2 * w * z - 2 * x * y
        out(3, 2) = 2 * y * z + 2 * w * x
        out(3, 3) = 1 - 2 * x ** 2 - 2 * y ** 2

    end function as_matrix

    pure function rotation_matrix( w, x, y, z ) result (rot_matrix)
        real(real64), intent(in)         :: w, x, y, z
        real(real64), dimension(3, 3)    :: rot_matrix

        rot_matrix(1, 1) = w ** 2 + x ** 2 - y ** 2 - z ** 2
        rot_matrix(1, 2) = 2.0 * ( x * y + w * z )
        rot_matrix(1, 3) = 2.0 * ( x * z - w * y )
        rot_matrix(2, 1) = 2.0 * ( x * y - w * z )
        rot_matrix(2, 2) = w ** 2 - x ** 2 + y ** 2 - z ** 2
        rot_matrix(2, 3) = 2.0 * ( y * z + w * x )
        rot_matrix(3, 1) = 2.0 * ( x * z + w * y )
        rot_matrix(3, 2) = 2.0 * ( y * z - w * x )
        rot_matrix(3, 3) = w ** 2 - x ** 2 - y ** 2 + z ** 2

    end function rotation_matrix

    pure function vvdot(a, b) result (out)
        implicit none
        real(real64), intent(in), dimension(3) :: a, b
        real(real64)                           :: out

        out = sum(a * b)

    end function vvdot

    pure function vmdot(a, b) result (out)
        implicit none
        real(real64), intent(in), dimension(3)    :: a
        real(real64), intent(in), dimension(3, 3) :: b
        real(real64), dimension(3)                :: out

        out(1) = a(1) * b(1, 1) + a(2) * b(2, 1) + a(3) * b(3, 1)
        out(2) = a(1) * b(1, 2) + a(2) * b(2, 2) + a(3) * b(3, 2)
        out(3) = a(1) * b(1, 3) + a(2) * b(2, 3) + a(3) * b(3, 3)

    end function vmdot

    pure function mvdot(a, b) result (out)
        implicit none
        real(real64), intent(in), dimension(3, 3) :: a
        real(real64), intent(in), dimension(3)    :: b
        real(real64), dimension(3)                :: out

        out(1) = a(1, 1) * b(1) + a(1, 2) * b(2) + a(1, 3) * b(3)
        out(2) = a(2, 1) * b(1) + a(2, 2) * b(2) + a(2, 3) * b(3)
        out(3) = a(3, 1) * b(1) + a(3, 2) * b(2) + a(3, 3) * b(3)

    end function mvdot

    pure function mmdot(a, b) result (out)
        implicit none
        real(real64), intent(in), dimension(3, 3) :: a, b
        real(real64), dimension(3, 3)             :: out

        out(1, 1) = a(1, 1) * b(1, 1) + a(1, 2) * b(2, 1) + a(1, 3) * b(3, 1)
        out(1, 2) = a(1, 1) * b(1, 2) + a(1, 2) * b(2, 2) + a(1, 3) * b(3, 2)
        out(1, 3) = a(1, 1) * b(1, 3) + a(1, 2) * b(2, 3) + a(1, 3) * b(3, 3)

        out(2, 1) = a(2, 1) * b(1, 1) + a(2, 2) * b(2, 1) + a(2, 3) * b(3, 1)
        out(2, 2) = a(2, 1) * b(1, 2) + a(2, 2) * b(2, 2) + a(2, 3) * b(3, 2)
        out(2, 3) = a(2, 1) * b(1, 3) + a(2, 2) * b(2, 3) + a(2, 3) * b(3, 3)

        out(3, 1) = a(3, 1) * b(1, 1) + a(3, 2) * b(2, 1) + a(3, 3) * b(3, 1)
        out(3, 2) = a(3, 1) * b(1, 2) + a(3, 2) * b(2, 2) + a(3, 3) * b(3, 2)
        out(3, 3) = a(3, 1) * b(1, 3) + a(3, 2) * b(2, 3) + a(3, 3) * b(3, 3)

    end function mmdot

    pure function matmul_atba(a, b) result (out)
        implicit none
        real(real64), intent(in), dimension(3, 3) :: a, b
        real(real64), dimension(3, 3)             :: out

        out(1,1) =   A(1, 1) * B(1, 1) * A(1, 1)   +   A(2, 1) * B(2, 2) * A(2, 1)   +   A(3, 1) * B(3, 3) * A(3, 1)
        out(1,2) =   A(1, 1) * B(1, 1) * A(1, 2)   +   A(2, 1) * B(2, 2) * A(2, 2)   +   A(3, 1) * B(3, 3) * A(3, 2)
        out(1,3) =   A(1, 1) * B(1, 1) * A(1, 3)   +   A(2, 1) * B(2, 2) * A(2, 3)   +   A(3, 1) * B(3, 3) * A(3, 3)
 
        out(2,1) =   out(1, 2)
        out(2,2) =   A(1, 2) * B(1, 1) * A(1, 2)   +   A(2, 2) * B(2, 2) * A(2, 2)   +   A(3, 2) * B(3, 3) * A(3, 2)
        out(2,3) =   A(1, 2) * B(1, 1) * A(1, 3)   +   A(2, 2) * B(2, 2) * A(2, 3)   +   A(3, 2) * B(3, 3) * A(3, 3) 

        out(3,1) =   out(1, 3)
        out(3,2) =   out(2, 3)
        out(3,3) =   A(1, 3) * B(1, 1) * A(1, 3)   +   A(2, 3) * B(2, 2) * A(2, 3)   +   A(3, 3) * B(3, 3) * A(3, 3) 
 
    end function matmul_atba

    pure function cholesky_solve(A, b) result (out)
        real(real64), intent(in)  :: A(3,3), b(3)
        real(real64)              :: out(3)
        real(real64)              :: inv_l11, inv_l22, inv_l33, l21, l31, l32
        real(real64)              :: y1, y2, y3

        inv_l11 = 1.0 / sqrt(A(1,1))
        l21     = A(2,1) * inv_l11
        l31     = A(3,1) * inv_l11

        inv_l22 = 1.0 / sqrt(A(2,2) - l21 * l21)
        l32     = (A(3,2) - l31 * l21) * inv_l22

        inv_l33 = 1.0 / sqrt(A(3,3) - l31 * l31 - l32 * l32)

        y1 = b(1) * inv_l11
        y2 = (b(2) - l21 * y1) * inv_l22
        y3 = (b(3) - l31 * y1 - l32 * y2) * inv_l33

        out(3) = y3 * inv_l33
        out(2) = (y2 - l32 * out(3)) * inv_l22
        out(1) = (y1 - l21 * out(2) - l31 * out(3)) * inv_l11

    end function cholesky_solve

    pure function quat_product(q1, q2) result (out)
        implicit none
        real(real64), dimension(4), intent(in) :: q1, q2
        real(real64), dimension(4)             :: out

        out(1) = q1(1) * q2(1) - q1(2) * q2(2) - q1(3) * q2(3) - q1(4) * q2(4)
        out(2) = q1(1) * q2(2) + q1(2) * q2(1) + q1(3) * q2(4) - q1(4) * q2(3)
        out(3) = q1(1) * q2(3) - q1(2) * q2(4) + q1(3) * q2(1) + q1(4) * q2(2)
        out(4) = q1(1) * q2(4) + q1(2) * q2(3) - q1(3) * q2(2) + q1(4) * q2(1)

    end function quat_product

end module xmath