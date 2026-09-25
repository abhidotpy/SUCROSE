module verlet
    use system
    use constraints
    implicit none

    contains
    subroutine vv_initial_step( DT )
        implicit none
        real(real64), intent(in) :: DT
        real(real64)             :: rx_pbc, ry_pbc, rz_pbc
        integer                  :: I


        do I = 1, N
                
            VX(I) = VX(I) + 0.5 * DT * ( FX(I) * inv_mass(I) )
            VY(I) = VY(I) + 0.5 * DT * ( FY(I) * inv_mass(I) )
            VZ(I) = VZ(I) + 0.5 * DT * ( FZ(I) * inv_mass(I) )

            RX(I) = RX(I) + DT * VX(I)
            RY(I) = RY(I) + DT * VY(I)
            RZ(I) = RZ(I) + DT * VZ(I)

            RX_PBC = RX(I) - ANINT( RX(I) / box(1) ) * box(1)
            RY_PBC = RY(I) - ANINT( RY(I) / box(2) ) * box(2)
            RZ_PBC = RZ(I) - ANINT( RZ(I) / box(3) ) * box(3)
            
            RX(I) = MERGE( RX_PBC, RX(I), periodic )
            RY(I) = MERGE( RY_PBC, RY(I), periodic )
            RZ(I) = MERGE( RZ_PBC, RZ(I), periodic )

        enddo

    end subroutine vv_initial_step

    subroutine vv_final_step( DT )
        implicit none
        real(real64), intent(in) :: DT
        integer                  :: I

        kin_eng = 0.0

        do I = 1, N

            VX(I) = VX(I) + 0.5 * DT * ( FX(I) * inv_mass(I) )
            VY(I) = VY(I) + 0.5 * DT * ( FY(I) * inv_mass(I) )
            VZ(I) = VZ(I) + 0.5 * DT * ( FZ(I) * inv_mass(I) )

            kin_eng = kin_eng + 0.5 * mass(I) * ( VX(I) * VX(I) + VY(I) * VY(I) + VZ(I) * VZ(I) )

        enddo

    end subroutine vv_final_step

    subroutine vvc_initial_step( DT )
        implicit none
        real(real64), intent(in) :: DT
        real(real64)             :: rx_pbc, ry_pbc, rz_pbc
        integer                  :: I

        do I = 1, N
                
            VX(I) = VX(I) + 0.5 * DT * ( FX(I) * inv_mass(I) )
            VY(I) = VY(I) + 0.5 * DT * ( FY(I) * inv_mass(I) )
            VZ(I) = VZ(I) + 0.5 * DT * ( FZ(I) * inv_mass(I) )

            RX_old(I) = RX(I)
            RY_old(I) = RY(I)
            RZ_old(I) = RZ(I)

            RX(I) = RX(I) + DT * VX(I)
            RY(I) = RY(I) + DT * VY(I)
            RZ(I) = RZ(I) + DT * VZ(I)

            RX_PBC = RX(I) - ANINT( RX(I) / box(1) ) * box(1)
            RY_PBC = RY(I) - ANINT( RY(I) / box(2) ) * box(2)
            RZ_PBC = RZ(I) - ANINT( RZ(I) / box(3) ) * box(3)
            
            RX(I) = MERGE( RX_PBC, RX(I), periodic )
            RY(I) = MERGE( RY_PBC, RY(I), periodic )
            RZ(I) = MERGE( RZ_PBC, RZ(I), periodic )

        enddo

        call apply_constraints_a( DT )

    end subroutine vvc_initial_step

    subroutine vvc_final_step( DT )
        implicit none
        real(real64), intent(in) :: DT
        integer                  :: I

        kin_eng = 0.0

        do I = 1, N

            VX(I) = VX(I) + 0.5 * DT * ( FX(I) * inv_mass(I) )
            VY(I) = VY(I) + 0.5 * DT * ( FY(I) * inv_mass(I) )
            VZ(I) = VZ(I) + 0.5 * DT * ( FZ(I) * inv_mass(I) )

            kin_eng = kin_eng + 0.5 * mass(I) * ( VX(I) * VX(I) + VY(I) * VY(I) + VZ(I) * VZ(I) )

        enddo

        call apply_constraints_b( DT )

    end subroutine vvc_final_step

end module verlet

module quaternion
    use xmath
    use system
    implicit none

    real(real64), private :: RM(3, 3)
    real(real64), private :: QDW, QDX, QDY, QDZ
    real(real64), private :: LX_body, LY_body, LZ_body
    real(real64), private :: WX_body, WY_body, WZ_body

    contains
    subroutine qq_initial_step( DT )
        implicit none
        real(real64), intent(in) :: DT
        real(real64)             :: QW_Old, QX_Old, QY_Old, QZ_Old
        real(real64)             :: QW_New, QX_New, QY_New, QZ_New
        real(real64)             :: Q_mag
        integer                  :: I, J

        do I = 1, N
            
            QW_New = QW(I); QX_New = QX(I)
            QY_New = QY(I); QZ_New = QZ(I)

            LX(I) = LX(I) + 0.5 * DT * TX(I)
            LY(I) = LY(I) + 0.5 * DT * TY(I)
            LZ(I) = LZ(I) + 0.5 * DT * TZ(I)

            do J = 1, 2
                QW_Old = QW_New; QX_Old = QX_New
                QY_Old = QY_New; QZ_Old = QZ_New
                
                RM(1, 1) = QW_old**2 + QX_old**2 - QY_old**2 - QZ_old**2
                RM(1, 2) = 2.0 * ( QX_old * QY_old + QW_old * QZ_old )
                RM(1, 3) = 2.0 * ( QX_old * QZ_old - QW_old * QY_old )
                RM(2, 1) = 2.0 * ( QX_old * QY_old - QW_old * QZ_old )
                RM(2, 2) = QW_old**2 - QX_old**2 + QY_old**2 - QZ_old**2
                RM(2, 3) = 2.0 * ( QY_old * QZ_old + QW_old * QX_old )
                RM(3, 1) = 2.0 * ( QX_old * QZ_old + QW_old * QY_old )
                RM(3, 2) = 2.0 * ( QY_old * QZ_old - QW_old * QX_old )
                RM(3, 3) = QW_old**2 - QX_old**2 - QY_old**2 + QZ_old**2

                WX_body = ( RM(1,1) * LX(I) + RM(1,2) * LY(I) + RM(1,3) * LZ(I) ) * Ixx_inv(I)
                WY_body = ( RM(2,1) * LX(I) + RM(2,2) * LY(I) + RM(2,3) * LZ(I) ) * Iyy_inv(I)
                WZ_body = ( RM(3,1) * LX(I) + RM(3,2) * LY(I) + RM(3,3) * LZ(I) ) * Izz_inv(I)

                QDW = 0.5 * ( - QX_Old * WX_body - QY_Old * WY_body - QZ_Old * WZ_body )
                QDX = 0.5 * (   QW_Old * WX_body + QY_Old * WZ_body - QZ_Old * WY_body )
                QDY = 0.5 * (   QW_Old * WY_body - QX_Old * WZ_body + QZ_Old * WX_body )
                QDZ = 0.5 * (   QW_Old * WZ_body + QX_Old * WY_body - QY_Old * WX_body )

                QW_New = QW(I) + 0.5 * DT * QDW
                QX_New = QX(I) + 0.5 * DT * QDX
                QY_New = QY(I) + 0.5 * DT * QDY
                QZ_New = QZ(I) + 0.5 * DT * QDZ
            
            enddo

            QW(I) = QW(I) + DT * QDW
            QX(I) = QX(I) + DT * QDX
            QY(I) = QY(I) + DT * QDY
            QZ(I) = QZ(I) + DT * QDZ

            Q_mag = SQRT( QW(I)**2 + QX(I)**2 + QY(I)**2 + QZ(I)**2 )
            QW(I) = QW(I) / Q_mag
            QX(I) = QX(I) / Q_mag
            QY(I) = QY(I) / Q_mag
            QZ(I) = QZ(I) / Q_mag

        enddo

    end subroutine qq_initial_step

    subroutine qq_final_step( DT )
        implicit none
        real(real64), intent(in) :: DT
        integer                  :: I

        kin_rot = 0.0

        do I = 1, N

            LX(I) = LX(I) + 0.5 * DT * TX(I)
            LY(I) = LY(I) + 0.5 * DT * TY(I)
            LZ(I) = LZ(I) + 0.5 * DT * TZ(I)

            RM(1, 1) = QW(I)**2 + QX(I)**2 - QY(I)**2 - QZ(I)**2
            RM(1, 2) = 2.0 * ( QX(I) * QY(I) + QW(I) * QZ(I) )
            RM(1, 3) = 2.0 * ( QX(I) * QZ(I) - QW(I) * QY(I) )
            RM(2, 1) = 2.0 * ( QX(I) * QY(I) - QW(I) * QZ(I) )
            RM(2, 2) = QW(I)**2 - QX(I)**2 + QY(I)**2 - QZ(I)**2
            RM(2, 3) = 2.0 * ( QY(I) * QZ(I) + QW(I) * QX(I) )
            RM(3, 1) = 2.0 * ( QX(I) * QZ(I) + QW(I) * QY(I) )
            RM(3, 2) = 2.0 * ( QY(I) * QZ(I) - QW(I) * QX(I) )
            RM(3, 3) = QW(I)**2 - QX(I)**2 - QY(I)**2 + QZ(I)**2

            LX_body = ( RM(1,1) * LX(I) + RM(1,2) * LY(I) + RM(1,3) * LZ(I) )
            LY_body = ( RM(2,1) * LX(I) + RM(2,2) * LY(I) + RM(2,3) * LZ(I) )
            LZ_body = ( RM(3,1) * LX(I) + RM(3,2) * LY(I) + RM(3,3) * LZ(I) )

            kin_rot = kin_rot + 0.5 * ( ( LX_body ** 2.0) * Ixx_inv(I) + &
                                        ( LY_body ** 2.0) * Iyy_inv(I) + &
                                        ( LZ_body ** 2.0) * Izz_inv(I) )

        enddo

    end subroutine qq_final_step

end module quaternion

module nose_hoover
    use system
    use constraints
    implicit none

    integer, parameter            :: N_th = 3
    real(real64), parameter       :: tau = 2.0
    real(real64)                  :: n_dof
    real(real64), dimension(N_th) :: HB = 1.0, eta = 0.0, p_eta = 0.0

    contains
    pure function polyval ( x, c ) RESULT ( f )
        implicit none
        real(real64),                intent(in) :: x ! argument
        real(real64), dimension(0:), intent(in) :: c ! given coefficients (ascending powers of x)
        real(real64)                            :: f ! Returns polynomial in ...
        integer                                 :: i, upper

        ! Uses Horner's rule

        upper = UBOUND(c,1)
        f = c(upper)

        do i = upper - 1, 0, -1
            f = f * x + c(i)
        enddo

    end function polyval

    pure function exprel ( x ) RESULT ( f )
        implicit none
        real(real64), intent(in) :: x ! Argument
        real(real64)             :: f ! Returns value of (exp(x)-1)/x

        ! At small x, we must guard against the ratio of imprecise small values.
        ! There are various ways of doing this.
        ! We follow some others and use the identity: (exp(x)-1)/x = exp(x/2)*[sinh(x/2)/(x/2)].
        ! For small x, sinh(x)/x = g0 + g1*x**2 + g2*x**4 + ...
        ! where the coefficient of x**(2n) is gn = 1/(2*n+1)!
        ! Alternatively, the exprel function is available in some math and scientific libraries.

        real(real64), dimension(0:4), parameter :: g = 1.0 / [1,6,120,5040,362880]
        real(real64),                 parameter :: tol = 0.01

        if ( abs(x) > tol ) then
            f = ( exp(x) - 1.0 ) / x
        else
            f = exp(x/2) * polyval ( (x/2)**2, g )
        endif

    end function exprel

    subroutine u4_propagator ( t, j_start, j_stop, j_stride )
        implicit none
        real(real64), intent(in)    :: t               
        integer, intent(in)         :: j_start, j_stop
        integer                     :: j, j_stride
        real(real64)                :: gj, x, c

        n_dof = 3 * ( N / Nrigid ) - 3
        do j = j_start, j_stop, j_stride

            if ( j == 1 ) then
                gj = SUM( mass * ( VX ** 2 + VY ** 2 + VZ ** 2 ) ) - n_dof * target_temp
            else
                gj = ( p_eta(j-1)**2 / HB(j-1) ) - target_temp
            endif

            if ( j == N_th ) then
                p_eta(j)  = p_eta(j) + t * gj
            else
                x = t * p_eta(j+1)/HB(j+1)
                c = exprel(-x) ! (1-exp(-x))/x, preserving accuracy for small x

                p_eta(j) = p_eta(j)*EXP(-x) + t * gj * c
            endif

        enddo

    end subroutine u4_propagator

    subroutine nht_initial_step( DT, coupling )
        implicit none
        real(real64),           intent(in) :: DT
        real(real64), optional, intent(in) :: coupling

        real(real64)    :: rx_pbc, ry_pbc, rz_pbc
        integer         :: I

        if ( .not. target_temp_set ) then
            write(*, "(1x, 'Error: Target temperature not set.')")
            stop
        endif

        n_dof = 3 * ( N / Nrigid ) - 3
        if (present(coupling)) then
            HB = coupling
        else
            HB = target_temp * tau ** 2.0
            HB(1) = n_dof * target_temp * tau ** 2.0
        endif

        do I = 1, N

            call u4_propagator( DT / 4.0, N_th, 1, -1)

            VX(I) = VX(I) * EXP( -0.5 * DT * p_eta(1) / HB(1) )
            VY(I) = VY(I) * EXP( -0.5 * DT * p_eta(1) / HB(1) )
            VZ(I) = VZ(I) * EXP( -0.5 * DT * p_eta(1) / HB(1) )

            eta = eta + 0.5 * DT * p_eta / HB

            CALL u4_propagator( DT / 4.0, 1, N_th, 1)

            VX(I) = VX(I) + 0.5 * DT * ( FX(I) * inv_mass(I) )
            VY(I) = VY(I) + 0.5 * DT * ( FY(I) * inv_mass(I) )
            VZ(I) = VZ(I) + 0.5 * DT * ( FZ(I) * inv_mass(I) )

            RX(I) = RX(I) + DT * VX(I)
            RY(I) = RY(I) + DT * VY(I)
            RZ(I) = RZ(I) + DT * VZ(I)

            RX_PBC = RX(I) - ANINT( RX(I) / box(1) ) * box(1)
            RY_PBC = RY(I) - ANINT( RY(I) / box(2) ) * box(2)
            RZ_PBC = RZ(I) - ANINT( RZ(I) / box(3) ) * box(3)

            RX(I) = MERGE( RX_PBC, RX(I), periodic )
            RY(I) = MERGE( RY_PBC, RY(I), periodic )
            RZ(I) = MERGE( RZ_PBC, RZ(I), periodic )

        enddo

    end subroutine nht_initial_step

    subroutine nht_final_step( DT, coupling )
        implicit none
        real(real64),           intent(in) :: DT
        real(real64), optional, intent(in) :: coupling
        integer                            :: I

        n_dof = 3 * ( N / Nrigid ) - 3
        if (present(coupling)) then
            HB = coupling
        else
            HB = target_temp * tau ** 2.0
            HB(1) = n_dof * target_temp * tau ** 2.0
        endif

        kin_eng = 0.0

        do I = 1, N

            VX(I) = VX(I) + 0.5 * DT * ( FX(I) * inv_mass(I) )
            VY(I) = VY(I) + 0.5 * DT * ( FY(I) * inv_mass(I) )
            VZ(I) = VZ(I) + 0.5 * DT * ( FZ(I) * inv_mass(I) )
            
            call u4_propagator( DT / 4.0, N_th, 1, -1)
            
            VX(I) = VX(I) * EXP( -0.5 * DT * p_eta(1) / HB(1) )
            VY(I) = VY(I) * EXP( -0.5 * DT * p_eta(1) / HB(1) )
            VZ(I) = VZ(I) * EXP( -0.5 * DT * p_eta(1) / HB(1) )
            
            eta = eta + 0.5 * DT * p_eta / HB
            
            CALL u4_propagator( DT / 4.0, 1, N_th, 1)

            kin_eng = kin_eng + 0.5 * mass(I) * ( VX(I) * VX(I) + VY(I) * VY(I) + VZ(I) * VZ(I) )

        enddo

    end subroutine nht_final_step

    subroutine nhtc_initial_step( DT, coupling )
        implicit none
        real(real64), intent(in) :: DT
        real(real64), optional, intent(in) :: coupling
        real(real64)    :: rx_pbc, ry_pbc, rz_pbc
        integer         :: I

        if ( .not. target_temp_set ) then
            write(*, "(1x, 'Error: Target temperature not set.')")
            stop
        endif

        n_dof = 3 * ( N / Nrigid ) - 3
        if (present(coupling)) then
            HB = coupling
        else
            HB = target_temp * tau ** 2.0
            HB(1) = n_dof * target_temp * tau ** 2.0
        endif

        do I = 1, N

            call u4_propagator( DT / 4.0, N_th, 1, -1)

            VX(I) = VX(I) * EXP( -0.5 * DT * p_eta(1) / HB(1) )
            VY(I) = VY(I) * EXP( -0.5 * DT * p_eta(1) / HB(1) )
            VZ(I) = VZ(I) * EXP( -0.5 * DT * p_eta(1) / HB(1) )

            eta = eta + 0.5 * DT * p_eta / HB

            CALL u4_propagator( DT / 4.0, 1, N_th, 1)

            VX(I) = VX(I) + 0.5 * DT * ( FX(I) * inv_mass(I) )
            VY(I) = VY(I) + 0.5 * DT * ( FY(I) * inv_mass(I) )
            VZ(I) = VZ(I) + 0.5 * DT * ( FZ(I) * inv_mass(I) )

            RX_old(I) = RX(I)
            RY_old(I) = RY(I)
            RZ_old(I) = RZ(I)

            RX(I) = RX(I) + DT * VX(I)
            RY(I) = RY(I) + DT * VY(I)
            RZ(I) = RZ(I) + DT * VZ(I)

            RX_PBC = RX(I) - ANINT( RX(I) / box(1) ) * box(1)
            RY_PBC = RY(I) - ANINT( RY(I) / box(2) ) * box(2)
            RZ_PBC = RZ(I) - ANINT( RZ(I) / box(3) ) * box(3)

            RX(I) = MERGE( RX_PBC, RX(I), periodic )
            RY(I) = MERGE( RY_PBC, RY(I), periodic )
            RZ(I) = MERGE( RZ_PBC, RZ(I), periodic )

        enddo

        call apply_constraints_a( DT )

    end subroutine nhtc_initial_step

    subroutine nhtc_final_step( DT, coupling )
        implicit none
        real(real64), intent(in) :: DT
        real(real64), optional, intent(in) :: coupling
        integer         :: I

        n_dof = 3 * ( N / Nrigid ) - 3
        if (present(coupling)) then
            HB = coupling
        else
            HB = target_temp * tau ** 2.0
            HB(1) = n_dof * target_temp * tau ** 2.0
        endif

        kin_eng = 0.0

        do I = 1, N

            VX(I) = VX(I) + 0.5 * DT * ( FX(I) * inv_mass(I) )
            VY(I) = VY(I) + 0.5 * DT * ( FY(I) * inv_mass(I) )
            VZ(I) = VZ(I) + 0.5 * DT * ( FZ(I) * inv_mass(I) )
            
            call u4_propagator( DT / 4.0, N_th, 1, -1)
            
            VX(I) = VX(I) * EXP( -0.5 * DT * p_eta(1) / HB(1) )
            VY(I) = VY(I) * EXP( -0.5 * DT * p_eta(1) / HB(1) )
            VZ(I) = VZ(I) * EXP( -0.5 * DT * p_eta(1) / HB(1) )
            
            eta = eta + 0.5 * DT * p_eta / HB
            
            CALL u4_propagator( DT / 4.0, 1, N_th, 1)

            kin_eng = kin_eng + 0.5 * mass(I) * ( VX(I) * VX(I) + VY(I) * VY(I) + VZ(I) * VZ(I) )

        enddo

        call apply_constraints_b( DT )

    end subroutine nhtc_final_step

end module nose_hoover

module langevin
    use philox_rng
    use system
    use constraints
    implicit none
    real(real64), private :: RM(3, 3)
    real(real64), private :: noise_1x, noise_1y, noise_1z, noise_2x, noise_2y, noise_2z
    real(real64), private :: QDW, QDX, QDY, QDZ
    real(real64), private :: LX_body, LY_body, LZ_body
    real(real64), private :: WX_body, WY_body, WZ_body
    real(real64)          :: kk1, kk2, damp = 1.0
    
    contains
    subroutine lgv_initial_step( DT, step, damping )
        implicit none
        real(real64), optional, intent(in)  :: damping
        real(real64), intent(in)            :: DT
        integer, intent(in)                 :: step

        real(real64)    :: RX_PBC, RY_PBC, RZ_PBC
        real(real64)    :: QW_Old, QX_Old, QY_Old, QZ_Old
        real(real64)    :: QW_New, QX_New, QY_New, QZ_New
        real(real64)    :: Q_mag
        integer         :: I, J, temp_int

        if ( .not. target_temp_set ) then
            write(*, "(1x, 'Error: Target temperature not set.')")
            stop
        endif

        if (present(damping)) then
            damp = damping
        else
            damp = 1.0
        endif

        kk1 = 1.0 - (damp * DT) / 2.0
        kk2 = 1.0 / (1.0 + (damp * DT) / 2.0)
        temp_int = int( target_temp * 1000.0 )
        noise_1x = 0.0; noise_1y = 0.0; noise_1z = 0.0
        noise_2x = 0.0; noise_2y = 0.0; noise_2z = 0.0

        do I = 1, N

            VX(I) = VX(I) + 0.5 * DT * ( FX(I) * inv_mass(I) )
            VY(I) = VY(I) + 0.5 * DT * ( FY(I) * inv_mass(I) )
            VZ(I) = VZ(I) + 0.5 * DT * ( FZ(I) * inv_mass(I) )

            RX(I) = RX(I) + 0.5 * DT * VX(I)
            RY(I) = RY(I) + 0.5 * DT * VY(I)
            RZ(I) = RZ(I) + 0.5 * DT * VZ(I)

            call philox_normal_3seq( I, step, temp_int, 0, noise_1x, noise_1y, noise_1z )
            call philox_normal_3seq( I, step, temp_int, 1, noise_2x, noise_2y, noise_2z )

            VX(I) = EXP( -damp * DT ) * VX(I) + sqrt( target_temp * inv_mass(I) * ( 1 - EXP( -2.0 * damp * DT ))) * noise_1x
            VY(I) = EXP( -damp * DT ) * VY(I) + sqrt( target_temp * inv_mass(I) * ( 1 - EXP( -2.0 * damp * DT ))) * noise_1y
            VZ(I) = EXP( -damp * DT ) * VZ(I) + sqrt( target_temp * inv_mass(I) * ( 1 - EXP( -2.0 * damp * DT ))) * noise_1z

            RX(I) = RX(I) + 0.5 * DT * VX(I)
            RY(I) = RY(I) + 0.5 * DT * VY(I)
            RZ(I) = RZ(I) + 0.5 * DT * VZ(I)

            RX_PBC = RX(I) - ANINT( RX(I) / box(1) ) * box(1)
            RY_PBC = RY(I) - ANINT( RY(I) / box(2) ) * box(2)
            RZ_PBC = RZ(I) - ANINT( RZ(I) / box(3) ) * box(3)

            RX(I) = MERGE( RX_PBC, RX(I), periodic )
            RY(I) = MERGE( RY_PBC, RY(I), periodic )
            RZ(I) = MERGE( RZ_PBC, RZ(I), periodic )

            QW_New = QW(I); QX_New = QX(I)
            QY_New = QY(I); QZ_New = QZ(I)

            LX(I) = LX(I) + 0.5 * DT * TX(I)
            LY(I) = LY(I) + 0.5 * DT * TY(I)
            LZ(I) = LZ(I) + 0.5 * DT * TZ(I)

            do J = 1, 2
                QW_Old = QW_New; QX_Old = QX_New
                QY_Old = QY_New; QZ_Old = QZ_New
                
                RM(1, 1) = QW_old**2 + QX_old**2 - QY_old**2 - QZ_old**2
                RM(1, 2) = 2.0 * ( QX_old * QY_old + QW_old * QZ_old )
                RM(1, 3) = 2.0 * ( QX_old * QZ_old - QW_old * QY_old )
                RM(2, 1) = 2.0 * ( QX_old * QY_old - QW_old * QZ_old )
                RM(2, 2) = QW_old**2 - QX_old**2 + QY_old**2 - QZ_old**2
                RM(2, 3) = 2.0 * ( QY_old * QZ_old + QW_old * QX_old )
                RM(3, 1) = 2.0 * ( QX_old * QZ_old + QW_old * QY_old )
                RM(3, 2) = 2.0 * ( QY_old * QZ_old - QW_old * QX_old )
                RM(3, 3) = QW_old**2 - QX_old**2 - QY_old**2 + QZ_old**2

                LX_body = ( RM(1,1) * LX(I) + RM(1,2) * LY(I) + RM(1,3) * LZ(I) )
                LY_body = ( RM(2,1) * LX(I) + RM(2,2) * LY(I) + RM(2,3) * LZ(I) )
                LZ_body = ( RM(3,1) * LX(I) + RM(3,2) * LY(I) + RM(3,3) * LZ(I) )

                WX_body = LX_body * Ixx_inv(I)
                WY_body = LY_body * Iyy_inv(I)
                WZ_body = LZ_body * Izz_inv(I)

                QDW = 0.5 * ( - QX_Old * WX_body - QY_Old * WY_body - QZ_Old * WZ_body )
                QDX = 0.5 * (   QW_Old * WX_body + QY_Old * WZ_body - QZ_Old * WY_body )
                QDY = 0.5 * (   QW_Old * WY_body - QX_Old * WZ_body + QZ_Old * WX_body )
                QDZ = 0.5 * (   QW_Old * WZ_body + QX_Old * WY_body - QY_Old * WX_body )

                QW_New = QW(I) + 0.5 * DT * QDW
                QX_New = QX(I) + 0.5 * DT * QDX
                QY_New = QY(I) + 0.5 * DT * QDY
                QZ_New = QZ(I) + 0.5 * DT * QDZ
            
            enddo

            LX_body = kk2 * kk1 * LX_body + kk2 * sqrt( 2.0 * damp * target_temp * DT * Ixx(I) ) * noise_2x
            LY_body = kk2 * kk1 * LY_body + kk2 * sqrt( 2.0 * damp * target_temp * DT * Iyy(I) ) * noise_2y
            LZ_body = kk2 * kk1 * LZ_body + kk2 * sqrt( 2.0 * damp * target_temp * DT * Izz(I) ) * noise_2z

            LX(I) = ( RM(1,1) * LX_body + RM(2,1) * LY_body + RM(3,1) * LZ_body )
            LY(I) = ( RM(1,2) * LX_body + RM(2,2) * LY_body + RM(3,2) * LZ_body )
            LZ(I) = ( RM(1,3) * LX_body + RM(2,3) * LY_body + RM(3,3) * LZ_body )

            QW(I) = QW(I) + DT * QDW
            QX(I) = QX(I) + DT * QDX
            QY(I) = QY(I) + DT * QDY
            QZ(I) = QZ(I) + DT * QDZ

            Q_mag = SQRT( QW(I)**2 + QX(I)**2 + QY(I)**2 + QZ(I)**2 )
            QW(I) = QW(I) / Q_mag
            QX(I) = QX(I) / Q_mag
            QY(I) = QY(I) / Q_mag
            QZ(I) = QZ(I) / Q_mag

        enddo

    end subroutine lgv_initial_step

    subroutine lgv_final_step( DT, damping )
        implicit none
        real(real64), optional, intent(in)   :: damping
        real(real64), intent(in)             :: DT
        integer                              :: I

        if (present(damping)) then
            damp = damping
        else
            damp = 1.0
        endif

        kin_eng = 0.0
        kin_rot = 0.0

        do I = 1, N

            VX(I) = VX(I) + 0.5 * DT * ( FX(I) * inv_mass(I) )
            VY(I) = VY(I) + 0.5 * DT * ( FY(I) * inv_mass(I) )
            VZ(I) = VZ(I) + 0.5 * DT * ( FZ(I) * inv_mass(I) )

            LX(I) = LX(I) + 0.5 * DT * TX(I)
            LY(I) = LY(I) + 0.5 * DT * TY(I)
            LZ(I) = LZ(I) + 0.5 * DT * TZ(I)

            RM(1, 1) = QW(I)**2 + QX(I)**2 - QY(I)**2 - QZ(I)**2
            RM(1, 2) = 2.0 * ( QX(I) * QY(I) + QW(I) * QZ(I) )
            RM(1, 3) = 2.0 * ( QX(I) * QZ(I) - QW(I) * QY(I) )
            RM(2, 1) = 2.0 * ( QX(I) * QY(I) - QW(I) * QZ(I) )
            RM(2, 2) = QW(I)**2 - QX(I)**2 + QY(I)**2 - QZ(I)**2
            RM(2, 3) = 2.0 * ( QY(I) * QZ(I) + QW(I) * QX(I) )
            RM(3, 1) = 2.0 * ( QX(I) * QZ(I) + QW(I) * QY(I) )
            RM(3, 2) = 2.0 * ( QY(I) * QZ(I) - QW(I) * QX(I) )
            RM(3, 3) = QW(I)**2 - QX(I)**2 - QY(I)**2 + QZ(I)**2

            LX_body = ( RM(1,1) * LX(I) + RM(1,2) * LY(I) + RM(1,3) * LZ(I) )
            LY_body = ( RM(2,1) * LX(I) + RM(2,2) * LY(I) + RM(2,3) * LZ(I) )
            LZ_body = ( RM(3,1) * LX(I) + RM(3,2) * LY(I) + RM(3,3) * LZ(I) )

            kin_eng = kin_eng + 0.5 * mass(I) * ( VX(I) * VX(I) + VY(I) * VY(I) + VZ(I) * VZ(I) )

            kin_rot = kin_rot + 0.5 * ( ( LX_body ** 2.0) * Ixx_inv(I) + &
                                        ( LY_body ** 2.0) * Iyy_inv(I) + &
                                        ( LZ_body ** 2.0) * Izz_inv(I) )

        enddo

    end subroutine lgv_final_step

    subroutine lgvc_initial_step( DT, step, damping )
        implicit none
        real(real64), optional, intent(in)   :: damping
        real(real64), intent(in)             :: DT
        integer, intent(in)                  :: step

        real(real64)    :: RX_PBC, RY_PBC, RZ_PBC
        real(real64)    :: QW_Old, QX_Old, QY_Old, QZ_Old
        real(real64)    :: QW_New, QX_New, QY_New, QZ_New
        real(real64)    :: Q_mag
        integer         :: I, J, temp_int

        if (target_temp .le. 1e-3) then
            write(*, "(1x, 'Error: Target temperature too small.')")
            stop
        endif

        if (present(damping)) then
            damp = damping
        else
            damp = 1.0
        endif

        kk1 = 1.0 - (damp * DT) / 2.0
        kk2 = 1.0 / (1.0 + (damp * DT) / 2.0)
        temp_int = int( target_temp * 1000.0 )
        noise_1x = 0.0; noise_1y = 0.0; noise_1z = 0.0
        noise_2x = 0.0; noise_2y = 0.0; noise_2z = 0.0

        do I = 1, N

            VX(I) = VX(I) + 0.5 * DT * ( FX(I) * inv_mass(I) )
            VY(I) = VY(I) + 0.5 * DT * ( FY(I) * inv_mass(I) )
            VZ(I) = VZ(I) + 0.5 * DT * ( FZ(I) * inv_mass(I) )

            RX_old(I) = RX(I)
            RY_old(I) = RY(I)
            RZ_old(I) = RZ(I)

            RX(I) = RX(I) + 0.5 * DT * VX(I)
            RY(I) = RY(I) + 0.5 * DT * VY(I)
            RZ(I) = RZ(I) + 0.5 * DT * VZ(I)

            call philox_normal_3seq( I, step, temp_int, 0, noise_1x, noise_1y, noise_1z )
            call philox_normal_3seq( I, step, temp_int, 1, noise_2x, noise_2y, noise_2z )

            VX(I) = EXP( -damp * DT ) * VX(I) + sqrt( target_temp * inv_mass(I) * ( 1 - EXP( -2.0 * damp * DT ))) * noise_1x
            VY(I) = EXP( -damp * DT ) * VY(I) + sqrt( target_temp * inv_mass(I) * ( 1 - EXP( -2.0 * damp * DT ))) * noise_1y
            VZ(I) = EXP( -damp * DT ) * VZ(I) + sqrt( target_temp * inv_mass(I) * ( 1 - EXP( -2.0 * damp * DT ))) * noise_1z

            RX(I) = RX(I) + 0.5 * DT * VX(I)
            RY(I) = RY(I) + 0.5 * DT * VY(I)
            RZ(I) = RZ(I) + 0.5 * DT * VZ(I)

            RX_PBC = RX(I) - ANINT( RX(I) / box(1) ) * box(1)
            RY_PBC = RY(I) - ANINT( RY(I) / box(2) ) * box(2)
            RZ_PBC = RZ(I) - ANINT( RZ(I) / box(3) ) * box(3)

            RX(I) = MERGE( RX_PBC, RX(I), periodic )
            RY(I) = MERGE( RY_PBC, RY(I), periodic )
            RZ(I) = MERGE( RZ_PBC, RZ(I), periodic )

            QW_New = QW(I); QX_New = QX(I)
            QY_New = QY(I); QZ_New = QZ(I)

            LX(I) = LX(I) + 0.5 * DT * TX(I)
            LY(I) = LY(I) + 0.5 * DT * TY(I)
            LZ(I) = LZ(I) + 0.5 * DT * TZ(I)

            do J = 1, 2
                QW_Old = QW_New; QX_Old = QX_New
                QY_Old = QY_New; QZ_Old = QZ_New
                
                RM(1, 1) = QW_old**2 + QX_old**2 - QY_old**2 - QZ_old**2
                RM(1, 2) = 2.0 * ( QX_old * QY_old + QW_old * QZ_old )
                RM(1, 3) = 2.0 * ( QX_old * QZ_old - QW_old * QY_old )
                RM(2, 1) = 2.0 * ( QX_old * QY_old - QW_old * QZ_old )
                RM(2, 2) = QW_old**2 - QX_old**2 + QY_old**2 - QZ_old**2
                RM(2, 3) = 2.0 * ( QY_old * QZ_old + QW_old * QX_old )
                RM(3, 1) = 2.0 * ( QX_old * QZ_old + QW_old * QY_old )
                RM(3, 2) = 2.0 * ( QY_old * QZ_old - QW_old * QX_old )
                RM(3, 3) = QW_old**2 - QX_old**2 - QY_old**2 + QZ_old**2

                LX_body = ( RM(1,1) * LX(I) + RM(1,2) * LY(I) + RM(1,3) * LZ(I) )
                LY_body = ( RM(2,1) * LX(I) + RM(2,2) * LY(I) + RM(2,3) * LZ(I) )
                LZ_body = ( RM(3,1) * LX(I) + RM(3,2) * LY(I) + RM(3,3) * LZ(I) )

                WX_body = LX_body * Ixx_inv(I)
                WY_body = LY_body * Iyy_inv(I)
                WZ_body = LZ_body * Izz_inv(I)

                QDW = 0.5 * ( - QX_Old * WX_body - QY_Old * WY_body - QZ_Old * WZ_body )
                QDX = 0.5 * (   QW_Old * WX_body + QY_Old * WZ_body - QZ_Old * WY_body )
                QDY = 0.5 * (   QW_Old * WY_body - QX_Old * WZ_body + QZ_Old * WX_body )
                QDZ = 0.5 * (   QW_Old * WZ_body + QX_Old * WY_body - QY_Old * WX_body )

                QW_New = QW(I) + 0.5 * DT * QDW
                QX_New = QX(I) + 0.5 * DT * QDX
                QY_New = QY(I) + 0.5 * DT * QDY
                QZ_New = QZ(I) + 0.5 * DT * QDZ
            
            enddo

            LX_body = kk2 * kk1 * LX_body + kk2 * sqrt( 2.0 * damp * target_temp * DT * Ixx(I) ) * noise_2x
            LY_body = kk2 * kk1 * LY_body + kk2 * sqrt( 2.0 * damp * target_temp * DT * Iyy(I) ) * noise_2y
            LZ_body = kk2 * kk1 * LZ_body + kk2 * sqrt( 2.0 * damp * target_temp * DT * Izz(I) ) * noise_2z

            LX(I) = ( RM(1,1) * LX_body + RM(2,1) * LY_body + RM(3,1) * LZ_body )
            LY(I) = ( RM(1,2) * LX_body + RM(2,2) * LY_body + RM(3,2) * LZ_body )
            LZ(I) = ( RM(1,3) * LX_body + RM(2,3) * LY_body + RM(3,3) * LZ_body )

            QW(I) = QW(I) + DT * QDW
            QX(I) = QX(I) + DT * QDX
            QY(I) = QY(I) + DT * QDY
            QZ(I) = QZ(I) + DT * QDZ

            Q_mag = SQRT( QW(I)**2 + QX(I)**2 + QY(I)**2 + QZ(I)**2 )
            QW(I) = QW(I) / Q_mag
            QX(I) = QX(I) / Q_mag
            QY(I) = QY(I) / Q_mag
            QZ(I) = QZ(I) / Q_mag

        enddo

        call apply_constraints_a( DT )

    end subroutine lgvc_initial_step

    subroutine lgvc_final_step( DT, damping )
        implicit none
        real(real64), optional, intent(in) :: damping
        real(real64), intent(in)           :: DT
        integer                            :: I

        if (present(damping)) then
            damp = damping
        else
            damp = 1.0
        endif

        kin_eng = 0.0
        kin_rot = 0.0

        do I = 1, N

            VX(I) = VX(I) + 0.5 * DT * ( FX(I) * inv_mass(I) )
            VY(I) = VY(I) + 0.5 * DT * ( FY(I) * inv_mass(I) )
            VZ(I) = VZ(I) + 0.5 * DT * ( FZ(I) * inv_mass(I) )

            LX(I) = LX(I) + 0.5 * DT * TX(I)
            LY(I) = LY(I) + 0.5 * DT * TY(I)
            LZ(I) = LZ(I) + 0.5 * DT * TZ(I)

            RM(1, 1) = QW(I)**2 + QX(I)**2 - QY(I)**2 - QZ(I)**2
            RM(1, 2) = 2.0 * ( QX(I) * QY(I) + QW(I) * QZ(I) )
            RM(1, 3) = 2.0 * ( QX(I) * QZ(I) - QW(I) * QY(I) )
            RM(2, 1) = 2.0 * ( QX(I) * QY(I) - QW(I) * QZ(I) )
            RM(2, 2) = QW(I)**2 - QX(I)**2 + QY(I)**2 - QZ(I)**2
            RM(2, 3) = 2.0 * ( QY(I) * QZ(I) + QW(I) * QX(I) )
            RM(3, 1) = 2.0 * ( QX(I) * QZ(I) + QW(I) * QY(I) )
            RM(3, 2) = 2.0 * ( QY(I) * QZ(I) - QW(I) * QX(I) )
            RM(3, 3) = QW(I)**2 - QX(I)**2 - QY(I)**2 + QZ(I)**2

            LX_body = ( RM(1,1) * LX(I) + RM(1,2) * LY(I) + RM(1,3) * LZ(I) )
            LY_body = ( RM(2,1) * LX(I) + RM(2,2) * LY(I) + RM(2,3) * LZ(I) )
            LZ_body = ( RM(3,1) * LX(I) + RM(3,2) * LY(I) + RM(3,3) * LZ(I) )

            kin_eng = kin_eng + 0.5 * mass(I) * ( VX(I) * VX(I) + VY(I) * VY(I) + VZ(I) * VZ(I) )

            kin_rot = kin_rot + 0.5 * ( ( LX_body ** 2.0) * Ixx_inv(I) + &
                                        ( LY_body ** 2.0) * Iyy_inv(I) + &
                                        ( LZ_body ** 2.0) * Izz_inv(I) )

        enddo

        call apply_constraints_b( DT )

    end subroutine lgvc_final_step

end module langevin

module monte_carlo_barostat
    use system
    implicit none

    contains
    subroutine mc_barostat( compute_potential, accept, volume_change )
        implicit none
        real(real64) :: V_old, V_new, volume, scale
        real(real64) :: U_old, U_new, dU
        real(real64) :: Lx_old, Ly_old, Lz_old
        real(real64) :: exponent, rand_u, beta
        real(real64) :: max_dV
        real(real64), optional, intent(in) :: volume_change
        logical, intent(out) :: accept

        interface
        real(real64) function compute_potential()
            import real64 
            end function compute_potential
        end interface
        
        if (target_pres .le. 1e-3) then
            write(*, "(1x, 'Error: Target pressure too small.')")
            stop
        endif

        if (present(volume_change)) then
            max_dV = volume_change
        else
            max_dV = 0.01
        endif

        beta = 1.0d0 / temp

        V_old  = product( box )
        Lx_old = box(1)
        Ly_old = box(2)
        Lz_old = box(3)

        RX_old = RX
        RY_old = RY
        RZ_old = RZ
        
        U_old = compute_potential()

        call random_number(rand_u)
        V_new = V_old + max_dV * (2.0d0 * rand_u - 1.0d0) * V_old

        if (V_new <= 0.0d0) then      ! unphysical — reject without counting
            return
        end if

        scale = (V_new / V_old)**(1.0d0 / 3.0d0)
        
        box(1) = Lx_old * scale
        box(2) = Ly_old * scale
        box(3) = Lz_old * scale

        RX = RX_old * scale
        RY = RY_old * scale
        RZ = RZ_old * scale

        U_new = compute_potential()
        dU = U_new - U_old

        exponent = - beta * dU - beta * target_pres * (V_new - V_old) + dble(N) * log(V_new / V_old)

        call random_number(rand_u)
        if (log(rand_u) < exponent) then
            accept = .true.
            volume = product( box )
            density = dble(N) / volume
        else

            box(1) = Lx_old
            box(2) = Ly_old
            box(3) = Lz_old

            RX = RX_old
            RY = RY_old
            RZ = RZ_old

            accept = .false.

        end if

    end subroutine mc_barostat

end module monte_carlo_barostat
    