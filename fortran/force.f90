module lennard_jones
    use system
    implicit none

    contains
    pure subroutine lj_calculate_forces( I, J, sigma, epsilon, cutoff, PE, FXIJ)
        implicit none
        integer, intent(in)         :: I, J
        real(real64), intent(in)    :: cutoff, sigma, epsilon
        real(real64), intent(out)   :: PE, FXIJ(3)

        real(real64)                :: RIJ(3), RIJ_pbc(3), rij_sq, rcut_sq
        real(real64)                :: sr2_lj, sr6_lj, sr12_lj
        real(real64)                :: sr2_cut, sr6_cut, sr12_cut
        real(real64)                :: pot, pot_cut, coeff
        logical                     :: cmask

        PE = 0.0; FXIJ = 0.0

        RIJ(1) = RX(I) - RX(J)
        RIJ(2) = RY(I) - RY(J)
        RIJ(3) = RZ(I) - RZ(J)
        RIJ_pbc = RIJ - ANINT( RIJ / box ) * box
        RIJ = MERGE( RIJ_pbc, RIJ, periodic )

        rij_sq  = RIJ(1) * RIJ(1) + RIJ(2) * RIJ(2) + RIJ(3) * RIJ(3)
        rcut_sq = cutoff * cutoff

        cmask = (rij_sq < rcut_sq)

        sr2_lj   = ( sigma * sigma ) / rij_sq
        sr6_lj   = sr2_lj * sr2_lj * sr2_lj
        sr12_lj  = sr6_lj * sr6_lj

        sr2_cut  = ( sigma * sigma ) / rcut_sq
        sr6_cut  = sr2_cut * sr2_cut * sr2_cut
        sr12_cut = sr6_cut * sr6_cut

        pot      = 4.0 * epsilon * ( sr12_lj - sr6_lj )
        pot_cut  = 4.0 * epsilon * ( sr12_cut - sr6_cut )
        coeff    = 24.0 * epsilon * ( 2.0 * sr12_lj - sr6_lj )

        PE       = ( pot - pot_cut )
        FXIJ     = ( coeff * RIJ / rij_sq )

        PE = merge( PE, 0.0, cmask )
        FXIJ = merge( FXIJ, 0.0, cmask )

    end subroutine lj_calculate_forces

    pure subroutine lj12_calculate_forces( I, J, sigma, epsilon, cutoff, PE, FXIJ )
        implicit none
        integer, intent(in)         :: I, J
        real(real64), intent(in)    :: cutoff, sigma, epsilon
        real(real64), intent(out)   :: PE, FXIJ(3)

        real(real64)                :: RIJ(3), rij_sq, rcut_sq
        real(real64)                :: sr2_lj, sr6_lj, sr12_lj
        real(real64)                :: sr2_cut, sr6_cut, sr12_cut
        real(real64)                :: pot, pot_cut, coeff

        PE = 0.0; FXIJ = 0.0

        RIJ(1) = RX(I) - RX(J)
        RIJ(2) = RY(I) - RY(J)
        RIJ(3) = RZ(I) - RZ(J)

        if (periodic) then
            RIJ = RIJ - ANINT( RIJ / box ) * box
        endif

        rij_sq  = RIJ(1) * RIJ(1) + RIJ(2) * RIJ(2) + RIJ(3) * RIJ(3)
        rcut_sq = cutoff * cutoff

        if (rij_sq < rcut_sq) then

            sr2_lj   = ( sigma * sigma ) / rij_sq
            sr6_lj   = sr2_lj * sr2_lj * sr2_lj
            sr12_lj  = sr6_lj * sr6_lj

            sr2_cut  = ( sigma * sigma ) / rcut_sq
            sr6_cut  = sr2_cut * sr2_cut * sr2_cut
            sr12_cut = sr6_cut * sr6_cut

            pot      = epsilon * ( sr12_lj - 2.0 * sr6_lj )
            pot_cut  = epsilon * ( sr12_cut - 2.0 * sr6_cut )
            coeff    = 12.0 * epsilon * ( sr12_lj - sr6_lj )

            PE       = pot - pot_cut
            FXIJ     = coeff * RIJ / rij_sq
        endif

    end subroutine lj12_calculate_forces

    pure subroutine ljw_calculate_forces( I, J, sigma, epsilon, cutoff, PE, FXIJ, wf )
        implicit none
        integer, intent(in)         :: I, J
        real(real64), intent(in)    :: cutoff, sigma, epsilon, wf
        real(real64), intent(out)   :: PE, FXIJ(3)

        real(real64)                :: RIJ(3), rij_sq, rij_mag
        real(real64)                :: sr_lj, sr2_lj, sr6_lj, sr12_lj
        real(real64)                :: sr_cut, sr2_cut, sr6_cut, sr12_cut
        real(real64)                :: pot, pot_cut, coeff

        PE = 0.0; FXIJ = 0.0

        RIJ(1) = RX(I) - RX(J)
        RIJ(2) = RY(I) - RY(J)
        RIJ(3) = RZ(I) - RZ(J)

        if (periodic) then
            RIJ = RIJ - ANINT( RIJ / box ) * box
        endif

        rij_sq  = RIJ(1) * RIJ(1) + RIJ(2) * RIJ(2) + RIJ(3) * RIJ(3)
        rij_mag = SQRT( rij_sq )

        if (rij_mag < (cutoff + wf)) then

            if (rij_mag < sigma) then
                sr_lj = sigma / rij_mag
            else if ((rij_mag .ge. sigma) .and. (rij_mag < (sigma + wf))) then
                sr_lj = 1.0
            else
                sr_lj = sigma / ( rij_mag - wf )
            endif

            sr2_lj   = sr_lj * sr_lj
            sr6_lj   = sr2_lj * sr2_lj * sr2_lj
            sr12_lj  = sr6_lj * sr6_lj

            sr_cut   = sigma / cutoff
            sr2_cut  = sr_cut * sr_cut
            sr6_cut  = sr2_cut * sr2_cut * sr2_cut
            sr12_cut = sr6_cut * sr6_cut

            pot      = epsilon * ( sr12_lj - 2.0 * sr6_lj )
            pot_cut  = epsilon * ( sr12_cut - 2.0 * sr6_cut )
            coeff    = 12.0 * epsilon * ( sr12_lj - sr6_lj ) * ( sr_lj / sigma )

            PE       = pot - pot_cut
            FXIJ     = coeff * RIJ / rij_mag
        endif

    end subroutine ljw_calculate_forces

end module lennard_jones

module wca
    use system
    implicit none

    contains
    pure subroutine wca_calculate_forces( I, J, sigma, epsilon, PE, FXIJ )
        implicit none
        integer, intent(in)         :: I, J
        real(real64), intent(in)    :: sigma, epsilon
        real(real64), intent(out)   :: PE, FXIJ(3)

        real(real64)                :: RIJ(3), rij_sq, sigma_sq
        real(real64)                :: sr2_lj, sr6_lj, sr12_lj
        real(real64)                :: coeff

        PE = 0.0; FXIJ = 0.0

        RIJ(1) = RX(I) - RX(J)
        RIJ(2) = RY(I) - RY(J)
        RIJ(3) = RZ(I) - RZ(J)

        if (periodic) then
            RIJ = RIJ - ANINT( RIJ / box ) * box
        endif

        rij_sq   = RIJ(1) * RIJ(1) + RIJ(2) * RIJ(2) + RIJ(3) * RIJ(3)
        sigma_sq = sigma * sigma

        if (rij_sq .le. sigma_sq) then
            
            sr2_lj  = sigma_sq / rij_sq
            sr6_lj  = sr2_lj * sr2_lj * sr2_lj
            sr12_lj = sr6_lj * sr6_lj

            PE      =  epsilon * ( sr12_lj - 2.0 * sr6_lj ) + epsilon
            coeff   = 12.0 * epsilon * ( sr12_lj - sr6_lj )
            FXIJ    = coeff * RIJ / rij_sq

        endif

    end subroutine wca_calculate_forces

end module wca

module morse
    use system
    implicit none

    contains
    pure subroutine morse_calculate_forces( I, J, dissoc, width, rmin, cutoff, PE, FXIJ )
        implicit none
        integer, intent(in)        :: I, J
        real(real64), intent(in)   :: dissoc, width, rmin, cutoff
        real(real64), intent(out)  :: PE, FXIJ(3)

        real(real64)               :: RIJ(3), rij_mag, rij_sq
        real(real64)               :: exp_term, exp_cut, pot, pot_cut

        PE = 0.0; FXIJ = 0.0

        RIJ(1) = RX(I) - RX(J)
        RIJ(2) = RY(I) - RY(J)
        RIJ(3) = RZ(I) - RZ(J)

        if (periodic) then
            RIJ = RIJ - ANINT( RIJ / box ) * box
        endif

        rij_sq  = RIJ(1) * RIJ(1) + RIJ(2) * RIJ(2) + RIJ(3) * RIJ(3)
        rij_mag = SQRT( rij_sq )

        if (rij_mag < cutoff) then

            exp_term = EXP( -width * ( rij_mag - rmin ) )
            exp_cut  = EXP( -width * ( cutoff - rmin ) )

            pot      = dissoc * ( ( 1.0 - exp_term ) * ( 1.0 - exp_term ) - 1.0 )
            pot_cut  = dissoc * ( ( 1.0 - exp_cut  ) * ( 1.0 - exp_cut  ) - 1.0 )

            PE       = pot - pot_cut
            FXIJ     = -2.0 * width * dissoc * (1.0 - exp_term) * exp_term * RIJ / rij_mag
        endif

    end subroutine morse_calculate_forces

end module morse

module gay_berne
    use system
    implicit none
    integer, parameter, private :: meu = 1, neu = 2
    
    type :: gb_params
        real(real64) :: ru1, ru2, uu
    end type
    
    type :: gb_pair
        real(real64), dimension(3)   :: U1, U2, rij
        real(real64)                 :: rij_sq
    end type

    contains
    pure function g_func( chie, param, pair ) result (res)
        implicit none
        real(real64), intent(in)    :: chie
        type(gb_params), intent(in) :: param
        type(gb_pair), intent(in)   :: pair
        real(real64)                :: term1, term2, res

        term1 = (param % ru1 + param % ru2)*(param % ru1 + param % ru2) / (1 + chie * param % uu)
        term2 = (param % ru1 - param % ru2)*(param % ru1 - param % ru2) / (1 - chie * param % uu)
        
        res = 1 - (chie / 2.0 / pair % rij_sq) * (term1 + term2)

    end function g_func

    pure function dG_dr( chie, param, pair ) result (res)
        implicit none
        real(real64), intent(in)    :: chie
        type(gb_params), intent(in) :: param
        type(gb_pair), intent(in)   :: pair
        real(real64), dimension(3)  :: term1, term2, res

        term1 = ((param % ru1 + param % ru2) / (1 + chie * param % uu)) * (pair % U1 + pair % U2)
        term2 = ((param % ru1 - param % ru2) / (1 - chie * param % uu)) * (pair % U1 - pair % U2)
        
        res = (-chie / pair % RIJ_SQ ) * (term1 + term2) + (2 * pair % RIJ / pair % RIJ_SQ) * (1 - g_func( chie, param, pair ))

    end function dG_dr

    pure function dG_du1( chie, param, pair ) result(res)
        implicit none
        real(real64), intent(in)    :: chie
        type(gb_params), intent(in) :: param
        type(gb_pair), intent(in)   :: pair
        real(real64)                :: term1, term2, res(3)
        
        term1 = ((param % ru1 + param % ru2) / (1 + chie * param % uu))
        term2 = ((param % ru1 - param % ru2) / (1 - chie * param % uu))
        
        res = (-chie * pair % RIJ / pair % rij_sq ) * (term1 + term2) + ((chie*chie * pair % U2) / (2 * pair % rij_sq)) * (term1*term1 - term2*term2)
        
    end function dG_du1

    pure function dG_du2( chie, param, pair ) result(res)
        implicit none
        real(real64), intent(in)    :: chie
        type(gb_params), intent(in) :: param
        type(gb_pair), intent(in)   :: pair
        real(real64)                :: term1, term2, res(3)
        
        term1 = ((param % ru1 + param % ru2) / (1 + chie * param % uu))
        term2 = ((param % ru1 - param % ru2) / (1 - chie * param % uu))
        
        res = (-chie * pair % RIJ / pair % RIJ_SQ ) * (term1 - term2) + ((chie*chie * pair % U1) / (2 * pair % rij_sq)) * (term1*term1 - term2*term2)
        
    end function dG_du2

    pure subroutine gb_calculate_forces( I, J, cutoff, PE, FXIJ, TXI, TXJ )
        implicit none
        integer, intent(in)                  :: I, J
        real(real64), intent(in)             :: cutoff
        real(real64), intent(out)            :: PE, FXIJ(3), TXI(3), TXJ(3)
        
        real(real64), dimension(3)           :: rij_hat
        real(real64), dimension(3)           :: dR_dr, de2_dr, dR_du1, dR_du2
        real(real64), dimension(3)           :: de1_du1, de2_du1, de1_du2, de2_du2
        real(real64), dimension(3)           :: FIJ, TI, TJ, TORQ1, TORQ2
        
        real(real64)                         :: muinv, chi, xhi
        real(real64)                         :: lpar, lprp, epar, eprp
        real(real64)                         :: gchi, gxhi, ginv_chi_3by2
        real(real64)                         :: eps1, eps2, sigma, sigma_0, sigma_0_inv
        real(real64)                         :: sr_gb, sr2_gb, sr6_gb, sr12_gb
        real(real64)                         :: pot, rij_mag, srx_fac, srd_fac
        real(real64)                         :: switch, xx, sx, sdx, factor
        type(gb_params)                      :: param
        type(gb_pair)                        :: pair

        pot = 0.0; FIJ = 0.0; TORQ1 = 0.0; TORQ2 = 0.0; sigma_0 = 1.0;
        
        pair % RIJ(1) = RX(I) - RX(J)
        pair % RIJ(2) = RY(I) - RY(J)
        pair % RIJ(3) = RZ(I) - RZ(J)

        if (periodic) then
            pair % RIJ = pair % RIJ - ANINT( pair % RIJ / box ) * box
        endif

        pair % rij_sq  = pair % RIJ(1) * pair % RIJ(1) + pair % RIJ(2) * pair % RIJ(2) + pair % RIJ(3) * pair % RIJ(3)
        rij_mag        = SQRT( pair % rij_sq )
        rij_hat        = pair % RIJ / rij_mag

        pair % U1(1) = 2.0 *  ( QX(I) * QZ(I) + QW(I) * QY(I) )
        pair % U1(2) = 2.0 *  ( QY(I) * QZ(I) - QW(I) * QX(I) )
        pair % U1(3) = QW(I)*QW(I) - QX(I)*QX(I) - QY(I)*QY(I) + QZ(I)*QZ(I)

        pair % U2(1) = 2.0 *  ( QX(J) * QZ(J) + QW(J) * QY(J) )
        pair % U2(2) = 2.0 *  ( QY(J) * QZ(J) - QW(J) * QX(J) )
        pair % U2(3) = QW(J)*QW(J) - QX(J)*QX(J) - QY(J)*QY(J) + QZ(J)*QZ(J)

        lpar = sig_z(I);  lprp = sig_x(I);
        epar = eps_z(I); eprp = eps_x(I);
        muinv = 1.0 / meu

        ! CALCULATE GAY BERNE PARAMETERS

        param % ru1 = pair % RIJ(1) * pair % U1(1) + pair % RIJ(2) * pair % U1(2) + pair % RIJ(3) * pair % U1(3)
        param % ru2 = pair % RIJ(1) * pair % U2(1) + pair % RIJ(2) * pair % U2(2) + pair % RIJ(3) * pair % U2(3)
        param % uu  = pair % U1(1)  * pair % U2(1) + pair % U1(2)  * pair % U2(2) + pair % U1(3)  * pair % U2(3)
        
        chi = ( lpar * lpar - lprp * lprp )     / ( lpar * lpar + lprp * lprp )
        xhi = ( eprp ** muinv - epar ** muinv ) / ( eprp ** muinv + epar ** muinv )
        gchi = g_func( chi, param, pair )
        gxhi = g_func( xhi, param, pair )

        eps1  = 1.0 / SQRT( 1 - ( chi * param % uu ) * ( chi * param % uu ) )
        eps2  = gxhi / eps_z(I)
        sigma = lprp / SQRT( gchi )
        sigma_0_inv = 1.0 / sigma_0
        ginv_chi_3by2 = sigma / gchi

        ! EVALUATE SMOOTHSTEP FUNCTION

        switch = 0.9 * ( sigma + cutoff - sigma_0 )
        factor = 1.0 / ( ( sigma + cutoff - sigma_0 ) - switch )
        xx = ( rij_mag - switch ) * factor
        sx = 1.0 - 6.0 * xx*xx*xx*xx*xx + 15.0 * xx*xx*xx*xx - 10.0 * xx*xx*xx
        sdx = -30.0 * xx*xx*xx*xx + 60.0 * xx*xx*xx - 30.0 * xx*xx

        ! EVALUATE GAY-BERNE POTENTIAL

        if ( rij_mag < ( sigma + cutoff - sigma_0 ) ) then

            sr_gb = sigma_0 / ( rij_mag - sigma + sigma_0 )
            sr2_gb = sr_gb * sr_gb
            sr6_gb = sr2_gb * sr2_gb * sr2_gb
            sr12_gb = sr6_gb * sr6_gb

            srx_fac = ( sr12_gb - 2.0 * sr6_gb )
            srd_fac = 12.0 * ( sr6_gb - sr12_gb ) * sr_gb

            dR_dr   = sigma_0_inv * ( rij_hat + 0.5 * dG_dr( chi, param, pair ) * ginv_chi_3by2 )
            de2_dr  = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_dr( xhi, param, pair )

            de1_du1 = (eps2 ** meu) * (neu * eps1 ** (neu + 2)) * (chi * chi * param % uu * pair % U2)
            de2_du1 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_du1( xhi, param, pair )
            dR_du1  = srd_fac * 0.5 * sigma_0_inv * dG_du1( chi, param, pair ) * ginv_chi_3by2

            de1_du2 = (eps2 ** meu) * (neu * eps1 ** (neu + 2)) * (chi * chi * param % uu * pair % U1)
            de2_du2 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_du2( xhi, param, pair )
            dR_du2  = srd_fac * 0.5 * sigma_0_inv * dG_du2( chi, param, pair ) * ginv_chi_3by2
        
            pot = (eps1 ** neu) * (eps2 ** meu) * srx_fac
            FIJ = (-1.0) * ((eps1 ** neu) * (eps2 ** meu) * srd_fac * dR_dr  + de2_dr * srx_fac)
            TI  = (-1.0) * (srx_fac * (de1_du1 + de2_du1) + ((eps1 ** neu) * (eps2 ** meu) * dR_du1))
            TJ  = (-1.0) * (srx_fac * (de1_du2 + de2_du2) + ((eps1 ** neu) * (eps2 ** meu) * dR_du2))
            TORQ1(1) = pair % U1(2) * TI(3) - pair % U1(3) * TI(2)
            TORQ1(2) = pair % U1(3) * TI(1) - pair % U1(1) * TI(3)
            TORQ1(3) = pair % U1(1) * TI(2) - pair % U1(2) * TI(1)
            TORQ2(1) = pair % U2(2) * TJ(3) - pair % U2(3) * TJ(2)
            TORQ2(2) = pair % U2(3) * TJ(1) - pair % U2(1) * TJ(3)
            TORQ2(3) = pair % U2(1) * TJ(2) - pair % U2(2) * TJ(1)

            if (( rij_mag >= switch ) .and. ( rij_mag <= (sigma + cutoff - sigma_0) )) then
                FIJ   = sx * FIJ   - sdx * pot * factor
                TORQ1 = sx * TORQ1 - sdx * pot * factor
                TORQ2 = sx * TORQ2 - sdx * pot * factor
                pot   = sx * pot
            endif

        endif

        PE = pot
        FXIJ = FIJ
        TXI = TORQ1
        TXJ = TORQ2

    end subroutine gb_calculate_forces

    pure subroutine gbr_calculate_forces( I, J, PE, FXIJ, TXI, TXJ )
        implicit none
        integer, intent(in)                  :: I, J
        real(real64), intent(out)            :: PE, FXIJ(3), TXI(3), TXJ(3)
        
        real(real64), dimension(3)           :: rij_hat
        real(real64), dimension(3)           :: dR_dr, de2_dr, dR_du1, dR_du2
        real(real64), dimension(3)           :: de1_du1, de2_du1, de1_du2, de2_du2
        real(real64), dimension(3)           :: FIJ, TI, TJ, TORQ1, TORQ2
        
        real(real64)                         :: muinv, chi, xhi
        real(real64)                         :: lpar, lprp, epar, eprp
        real(real64)                         :: gchi, gxhi, ginv_chi_3by2
        real(real64)                         :: eps1, eps2, sigma, sigma_0, sigma_0_inv
        real(real64)                         :: sr_gb, sr2_gb, sr6_gb, sr12_gb
        real(real64)                         :: pot, rij_mag, srx_fac, srd_fac
        type(gb_params)                      :: param
        type(gb_pair)                        :: pair

        pot = 0.0; FIJ = 0.0; TORQ1 = 0.0; TORQ2 = 0.0; sigma_0 = 1.0;

        pair % RIJ(1) = RX(I) - RX(J)
        pair % RIJ(2) = RY(I) - RY(J)
        pair % RIJ(3) = RZ(I) - RZ(J)

        if (periodic) then
            pair % RIJ = pair % RIJ - ANINT( pair % RIJ / box ) * box
        endif

        pair % rij_sq  = pair % RIJ(1) * pair % RIJ(1) + pair % RIJ(2) * pair % RIJ(2) + pair % RIJ(3) * pair % RIJ(3)
        rij_mag        = SQRT( pair % rij_sq )
        rij_hat        = pair % RIJ / rij_mag

        pair % U1(1) = 2.0 *  ( QX(I) * QZ(I) + QW(I) * QY(I) )
        pair % U1(2) = 2.0 *  ( QY(I) * QZ(I) - QW(I) * QX(I) )
        pair % U1(3) = QW(I)*QW(I) - QX(I)*QX(I) - QY(I)*QY(I) + QZ(I)*QZ(I)

        pair % U2(1) = 2.0 *  ( QX(J) * QZ(J) + QW(J) * QY(J) )
        pair % U2(2) = 2.0 *  ( QY(J) * QZ(J) - QW(J) * QX(J) )
        pair % U2(3) = QW(J)*QW(J) - QX(J)*QX(J) - QY(J)*QY(J) + QZ(J)*QZ(J)

        lpar = sig_z(I);  lprp = sig_x(I);
        epar = eps_z(I); eprp = eps_x(I);
        muinv = 1.0 / meu

        ! CALCULATE GAY BERNE PARAMETERS

        param % ru1 = pair % RIJ(1) * pair % U1(1) + pair % RIJ(2) * pair % U1(2) + pair % RIJ(3) * pair % U1(3)
        param % ru2 = pair % RIJ(1) * pair % U2(1) + pair % RIJ(2) * pair % U2(2) + pair % RIJ(3) * pair % U2(3)
        param % uu  = pair % U1(1)  * pair % U2(1) + pair % U1(2)  * pair % U2(2) + pair % U1(3)  * pair % U2(3)
        
        chi = ( lpar * lpar - lprp * lprp )     / ( lpar * lpar + lprp * lprp )
        xhi = ( eprp ** muinv - epar ** muinv ) / ( eprp ** muinv + epar ** muinv )
        gchi = g_func( chi, param, pair )
        gxhi = g_func( xhi, param, pair )

        eps1  = 1.0 / SQRT( 1 - ( chi * param % uu ) * ( chi * param % uu ) )
        eps2  = gxhi / eps_z(I)
        sigma = lprp / SQRT( gchi )
        sigma_0_inv = 1.0 / sigma_0
        ginv_chi_3by2 = sigma / gchi

        ! EVALUATE GAY-BERNE POTENTIAL

        if ( rij_mag <= sigma ) then

            sr_gb = sigma_0 / ( rij_mag - sigma + sigma_0 )
            sr2_gb = sr_gb * sr_gb
            sr6_gb = sr2_gb * sr2_gb * sr2_gb
            sr12_gb = sr6_gb * sr6_gb

            srx_fac = ( sr12_gb - 2.0 * sr6_gb ) + 1.0
            srd_fac = 12.0 * ( sr6_gb - sr12_gb ) * sr_gb

            dR_dr   = sigma_0_inv * ( rij_hat + 0.5 * dG_dr( chi, param, pair ) * ginv_chi_3by2 )
            de2_dr  = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_dr( xhi, param, pair )

            de1_du1 = (eps2 ** meu) * (neu * eps1 ** (neu + 2)) * (chi * chi * param % uu * pair % U2)
            de2_du1 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_du1( xhi, param, pair )
            dR_du1  = srd_fac * 0.5 * sigma_0_inv * dG_du1( chi, param, pair ) * ginv_chi_3by2

            de1_du2 = (eps2 ** meu) * (neu * eps1 ** (neu + 2)) * (chi * chi * param % uu * pair % U1)
            de2_du2 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_du2( xhi, param, pair )
            dR_du2  = srd_fac * 0.5 * sigma_0_inv * dG_du2( chi, param, pair ) * ginv_chi_3by2
        
            pot = (eps1 ** neu) * (eps2 ** meu) * srx_fac
            FIJ = (-1.0) * ((eps1 ** neu) * (eps2 ** meu) * srd_fac * dR_dr  + de2_dr * srx_fac)
            TI  = (-1.0) * (srx_fac * (de1_du1 + de2_du1) + ((eps1 ** neu) * (eps2 ** meu) * dR_du1))
            TJ  = (-1.0) * (srx_fac * (de1_du2 + de2_du2) + ((eps1 ** neu) * (eps2 ** meu) * dR_du2))
            TORQ1(1) = pair % U1(2) * TI(3) - pair % U1(3) * TI(2)
            TORQ1(2) = pair % U1(3) * TI(1) - pair % U1(1) * TI(3)
            TORQ1(3) = pair % U1(1) * TI(2) - pair % U1(2) * TI(1)
            TORQ2(1) = pair % U2(2) * TJ(3) - pair % U2(3) * TJ(2)
            TORQ2(2) = pair % U2(3) * TJ(1) - pair % U2(1) * TJ(3)
            TORQ2(3) = pair % U2(1) * TJ(2) - pair % U2(2) * TJ(1)

        endif

        PE = pot
        FXIJ = FIJ
        TXI = TORQ1
        TXJ = TORQ2

    end subroutine gbr_calculate_forces

    pure subroutine gbw_calculate_forces( I, J, cutoff, PE, FXIJ, TXI, TXJ, wf )
        implicit none
        integer, intent(in)                  :: I, J
        real(real64), intent(in)             :: cutoff, wf
        real(real64), intent(out)            :: PE, FXIJ(3), TXI(3), TXJ(3)
        
        real(real64), dimension(3)           :: rij_hat
        real(real64), dimension(3)           :: dR_dr, de2_dr, dR_du1, dR_du2
        real(real64), dimension(3)           :: de1_du1, de2_du1, de1_du2, de2_du2
        real(real64), dimension(3)           :: FIJ, TI, TJ, TORQ1, TORQ2
        
        real(real64)                         :: muinv, chi, xhi
        real(real64)                         :: lpar, lprp, epar, eprp
        real(real64)                         :: gchi, gxhi, ginv_chi_3by2
        real(real64)                         :: eps1, eps2, sigma, sigma_0, sigma_0_inv
        real(real64)                         :: sr_gb, sr2_gb, sr6_gb, sr12_gb
        real(real64)                         :: pot, rij_mag, srx_fac, srd_fac
        real(real64)                         :: switch, xx, sx, sdx, factor
        type(gb_params)                      :: param
        type(gb_pair)                        :: pair

       
        pot = 0.0; FIJ = 0.0; TORQ1 = 0.0; TORQ2 = 0.0; sigma_0 = 1.0;

        pair % RIJ(1) = RX(I) - RX(J)
        pair % RIJ(2) = RY(I) - RY(J)
        pair % RIJ(3) = RZ(I) - RZ(J)

        if (periodic) then
            pair % RIJ = pair % RIJ - ANINT( pair % RIJ / box ) * box
        endif

        pair % rij_sq  = pair % RIJ(1) * pair % RIJ(1) + pair % RIJ(2) * pair % RIJ(2) + pair % RIJ(3) * pair % RIJ(3)
        rij_mag        = SQRT( pair % rij_sq )
        rij_hat        = pair % RIJ / rij_mag

        pair % U1(1) = 2.0 *  ( QX(I) * QZ(I) + QW(I) * QY(I) )
        pair % U1(2) = 2.0 *  ( QY(I) * QZ(I) - QW(I) * QX(I) )
        pair % U1(3) = QW(I)*QW(I) - QX(I)*QX(I) - QY(I)*QY(I) + QZ(I)*QZ(I)

        pair % U2(1) = 2.0 *  ( QX(J) * QZ(J) + QW(J) * QY(J) )
        pair % U2(2) = 2.0 *  ( QY(J) * QZ(J) - QW(J) * QX(J) )
        pair % U2(3) = QW(J)*QW(J) - QX(J)*QX(J) - QY(J)*QY(J) + QZ(J)*QZ(J)

        lpar = sig_z(I);  lprp = sig_x(I);
        epar = eps_z(I); eprp = eps_x(I);
        muinv = 1.0 / meu

        ! CALCULATE GAY BERNE PARAMETERS

        param % ru1 = pair % RIJ(1) * pair % U1(1) + pair % RIJ(2) * pair % U1(2) + pair % RIJ(3) * pair % U1(3)
        param % ru2 = pair % RIJ(1) * pair % U2(1) + pair % RIJ(2) * pair % U2(2) + pair % RIJ(3) * pair % U2(3)
        param % uu  = pair % U1(1)  * pair % U2(1) + pair % U1(2)  * pair % U2(2) + pair % U1(3)  * pair % U2(3)
        
        chi = ( lpar * lpar - lprp * lprp )     / ( lpar * lpar + lprp * lprp )
        xhi = ( eprp ** muinv - epar ** muinv ) / ( eprp ** muinv + epar ** muinv )
        gchi = g_func( chi, param, pair )
        gxhi = g_func( xhi, param, pair )
        
        eps1  = 1.0 / SQRT( 1 - ( chi * param % uu ) * ( chi * param % uu ) )
        eps2  = gxhi / eps_z(I)
        sigma = lprp / SQRT( gchi )
        sigma_0_inv = 1.0 / sigma_0
        ginv_chi_3by2 = sigma / gchi

        ! EVALUATE SMOOTHSTEP FUNCTION

        switch = 0.9 * ( sigma + cutoff + wf - sigma_0 )
        factor = 1.0 / ( ( sigma + cutoff + wf - sigma_0 ) - switch )
        xx = ( rij_mag - switch ) * factor
        sx = 1.0 - 6.0 * xx*xx*xx*xx*xx + 15.0 * xx*xx*xx*xx - 10.0 * xx*xx*xx
        sdx = -30.0 * xx*xx*xx*xx + 60.0 * xx*xx*xx - 30.0 * xx*xx

        ! EVALUATE GAY-BERNE POTENTIAL

        if ( rij_mag < ( sigma + cutoff + wf - sigma_0 ) ) then

            if ( rij_mag < sigma ) then
                sr_gb = sigma_0 / ( rij_mag - sigma + sigma_0 )
            else if ((rij_mag .ge. sigma) .and. (rij_mag .le. (sigma + wf))) then
                sr_gb = 1.0
            else
                sr_gb = sigma_0 / ( rij_mag - sigma + sigma_0 - wf )
            endif

            sr2_gb = sr_gb * sr_gb
            sr6_gb = sr2_gb * sr2_gb * sr2_gb
            sr12_gb = sr6_gb * sr6_gb

            srx_fac = ( sr12_gb - 2.0 * sr6_gb )
            srd_fac = 12.0 * ( sr6_gb - sr12_gb ) * sr_gb

            dR_dr   = sigma_0_inv * ( rij_hat + 0.5 * dG_dr( chi, param, pair ) * ginv_chi_3by2 )
            de2_dr  = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_dr( xhi, param, pair )

            de1_du1 = (eps2 ** meu) * (neu * eps1 ** (neu + 2)) * (chi * chi * param % uu * pair % U2)
            de2_du1 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_du1( xhi, param, pair )
            dR_du1  = srd_fac * 0.5 * sigma_0_inv * dG_du1( chi, param, pair ) * ginv_chi_3by2

            de1_du2 = (eps2 ** meu) * (neu * eps1 ** (neu + 2)) * (chi * chi * param % uu * pair % U1)
            de2_du2 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_du2( xhi, param, pair )
            dR_du2  = srd_fac * 0.5 * sigma_0_inv * dG_du2( chi, param, pair ) * ginv_chi_3by2
        
            pot = (eps1 ** neu) * (eps2 ** meu) * srx_fac
            FIJ = (-1.0) * ((eps1 ** neu) * (eps2 ** meu) * srd_fac * dR_dr  + de2_dr * srx_fac)
            TI  = (-1.0) * (srx_fac * (de1_du1 + de2_du1) + ((eps1 ** neu) * (eps2 ** meu) * dR_du1))
            TJ  = (-1.0) * (srx_fac * (de1_du2 + de2_du2) + ((eps1 ** neu) * (eps2 ** meu) * dR_du2))
            TORQ1(1) = pair % U1(2) * TI(3) - pair % U1(3) * TI(2)
            TORQ1(2) = pair % U1(3) * TI(1) - pair % U1(1) * TI(3)
            TORQ1(3) = pair % U1(1) * TI(2) - pair % U1(2) * TI(1)
            TORQ2(1) = pair % U2(2) * TJ(3) - pair % U2(3) * TJ(2)
            TORQ2(2) = pair % U2(3) * TJ(1) - pair % U2(1) * TJ(3)
            TORQ2(3) = pair % U2(1) * TJ(2) - pair % U2(2) * TJ(1)

            if (( rij_mag >= switch ) .and. ( rij_mag <= (sigma + cutoff + wf - sigma_0) )) then
                FIJ   = sx * FIJ   - sdx * pot * factor
                TORQ1 = sx * TORQ1 - sdx * pot * factor
                TORQ2 = sx * TORQ2 - sdx * pot * factor
                pot   = sx * pot
            endif

        endif

        PE = pot
        FXIJ = FIJ
        TXI = TORQ1
        TXJ = TORQ2

    end subroutine gbw_calculate_forces

    pure subroutine gbc_calculate_forces( I, J, cutoff, PE, FXIJ, TXI, TXJ, chirality )
        implicit none
        integer, intent(in)                  :: I, J
        real(real64), intent(in)             :: cutoff, chirality
        real(real64), intent(out)            :: PE, FXIJ(3), TXI(3), TXJ(3)
        
        real(real64), dimension(3)           :: rij_hat, u2xrij, rijxu1, u1xu2
        real(real64), dimension(3)           :: dR_dr, de2_dr, dR_du1, dR_du2
        real(real64), dimension(3)           :: de1_du1, de2_du1, de1_du2, de2_du2
        real(real64), dimension(3)           :: duur_dr, duur_du1, duur_du2, duu_du1, duu_du2
        real(real64), dimension(3)           :: FIJ, TI, TJ, TORQ1, TORQ2
        real(real64), dimension(3)           :: FIJ_chiral, TI_chiral, TJ_chiral, TORQ1_chiral, TORQ2_chiral
        
        real(real64)                         :: muinv, chi, xhi, uur, uu
        real(real64)                         :: lpar, lprp, epar, eprp
        real(real64)                         :: gchi, gxhi, ginv_chi_3by2
        real(real64)                         :: eps1, eps2, sigma, sigma_0, sigma_0_inv
        real(real64)                         :: sr_gb, sr2_gb, sr6_gb, sr12_gb
        real(real64)                         :: pot, pot_chiral, rij_mag,  srx_fac, srd_fac
        real(real64)                         :: switch, xx, sx, sdx, factor
        type(gb_params)                      :: param
        type(gb_pair)                        :: pair
        

        pot = 0.0; FIJ = 0.0; TORQ1 = 0.0; TORQ2 = 0.0; sigma_0 = 1.0;
        pot_chiral = 0.0; FIJ_chiral = 0.0; TORQ1_chiral = 0.0; TORQ2_chiral = 0.0

        pair % RIJ(1) = RX(I) - RX(J)
        pair % RIJ(2) = RY(I) - RY(J)
        pair % RIJ(3) = RZ(I) - RZ(J)

        if (periodic) then
            pair % RIJ = pair % RIJ - ANINT( pair % RIJ / box ) * box
        endif

        pair % rij_sq  = pair % RIJ(1) * pair % RIJ(1) + pair % RIJ(2) * pair % RIJ(2) + pair % RIJ(3) * pair % RIJ(3)
        rij_mag        = SQRT( pair % rij_sq )
        rij_hat        = pair % RIJ / rij_mag

        pair % U1(1) = 2.0 *  ( QX(I) * QZ(I) + QW(I) * QY(I) )
        pair % U1(2) = 2.0 *  ( QY(I) * QZ(I) - QW(I) * QX(I) )
        pair % U1(3) = QW(I)*QW(I) - QX(I)*QX(I) - QY(I)*QY(I) + QZ(I)*QZ(I)

        pair % U2(1) = 2.0 *  ( QX(J) * QZ(J) + QW(J) * QY(J) )
        pair % U2(2) = 2.0 *  ( QY(J) * QZ(J) - QW(J) * QX(J) )
        pair % U2(3) = QW(J)*QW(J) - QX(J)*QX(J) - QY(J)*QY(J) + QZ(J)*QZ(J)

        lpar = sig_z(I);  lprp = sig_x(I);
        epar = eps_z(I); eprp = eps_x(I);
        muinv = 1.0 / meu

        ! CALCULATE GAY BERNE PARAMETERS

        param % ru1 = pair % RIJ(1) * pair % U1(1) + pair % RIJ(2) * pair % U1(2) + pair % RIJ(3) * pair % U1(3)
        param % ru2 = pair % RIJ(1) * pair % U2(1) + pair % RIJ(2) * pair % U2(2) + pair % RIJ(3) * pair % U2(3)
        param % uu  = pair % U1(1)  * pair % U2(1) + pair % U1(2)  * pair % U2(2) + pair % U1(3)  * pair % U2(3)
        
        chi = ( lpar * lpar - lprp * lprp )     / ( lpar * lpar + lprp * lprp )
        xhi = ( eprp ** muinv - epar ** muinv ) / ( eprp ** muinv + epar ** muinv )
        gchi = g_func( chi, param, pair )
        gxhi = g_func( xhi, param, pair )
        
        eps1  = 1.0 / SQRT( 1 - ( chi * param % uu ) * ( chi * param % uu ) )
        eps2  = gxhi / eps_z(I)
        sigma = lprp / SQRT( gchi )
        sigma_0_inv = 1.0 / sigma_0
        ginv_chi_3by2 = sigma / gchi

        ! EVALUATE SMOOTHSTEP FUNCTION

        switch = 0.9 * ( sigma + cutoff - sigma_0 )
        factor = 1.0 / ( ( sigma + cutoff - sigma_0 ) - switch )
        xx = ( rij_mag - switch ) * factor
        sx = 1.0 - 6.0 * xx*xx*xx*xx*xx + 15.0 * xx*xx*xx*xx - 10.0 * xx*xx*xx
        sdx = -30.0 * xx*xx*xx*xx + 60.0 * xx*xx*xx - 30.0 * xx*xx

        ! EVALUATE GAY-BERNE POTENTIAL

        if ( rij_mag < ( sigma + cutoff - sigma_0 ) ) then

            sr_gb = sigma_0 / ( rij_mag - sigma + sigma_0 )
            sr2_gb = sr_gb * sr_gb
            sr6_gb = sr2_gb * sr2_gb * sr2_gb
            sr12_gb = sr6_gb * sr6_gb

            srx_fac = ( sr12_gb - 2.0 * sr6_gb )
            srd_fac = 12.0 * ( sr6_gb - sr12_gb ) * sr_gb

            ! EVALUATE NORMAL GB POTENTIAL

            dR_dr   = (eps1 ** neu) * (eps2 ** meu) * sigma_0_inv * ( rij_hat + 0.5 * dG_dr( chi, param, pair ) * ginv_chi_3by2 )
            de2_dr  = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_dr( xhi, param, pair )

            de1_du1 = (eps2 ** meu) * (neu * eps1 ** (neu + 2)) * (chi * chi * param % uu * pair % U2)
            de2_du1 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_du1( xhi, param, pair )
            dR_du1  = (eps1 ** neu) * (eps2 ** meu) * 0.5 * sigma_0_inv * dG_du1( chi, param, pair ) * ginv_chi_3by2

            de1_du2 = (eps2 ** meu) * (neu * eps1 ** (neu + 2)) * (chi * chi * param % uu * pair % U1)
            de2_du2 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_du2( xhi, param, pair )
            dR_du2  = (eps1 ** neu) * (eps2 ** meu) * 0.5 * sigma_0_inv * dG_du2( chi, param, pair ) * ginv_chi_3by2
        
            pot = (eps1 ** neu) * (eps2 ** meu) * srx_fac
            FIJ = (-1.0) * ( srx_fac *  de2_dr             + srd_fac * dR_dr  )
            TI  = (-1.0) * ( srx_fac * (de1_du1 + de2_du1) + srd_fac * dR_du1 )
            TJ  = (-1.0) * ( srx_fac * (de1_du2 + de2_du2) + srd_fac * dR_du2 )
            TORQ1(1) = pair % U1(2) * TI(3) - pair % U1(3) * TI(2)
            TORQ1(2) = pair % U1(3) * TI(1) - pair % U1(1) * TI(3)
            TORQ1(3) = pair % U1(1) * TI(2) - pair % U1(2) * TI(1)
            TORQ2(1) = pair % U2(2) * TJ(3) - pair % U2(3) * TJ(2)
            TORQ2(2) = pair % U2(3) * TJ(1) - pair % U2(1) * TJ(3)
            TORQ2(3) = pair % U2(1) * TJ(2) - pair % U2(2) * TJ(1)

            ! EVALUATE CHIRAL TERMS

            srx_fac  = sr6_gb * sr_gb
            srd_fac  = -7.0 * sr6_gb * sr2_gb

            u1xu2(1) = pair % U1(2) * pair % U2(3) - pair % U1(3) * pair % U2(2)
            u1xu2(2) = pair % U1(3) * pair % U2(1) - pair % U1(1) * pair % U2(3)
            u1xu2(3) = pair % U1(1) * pair % U2(2) - pair % U1(2) * pair % U2(1) 

            u2xrij(1) = pair % U2(2) * rij_hat(3) - pair % U2(3) * rij_hat(2)
            u2xrij(2) = pair % U2(3) * rij_hat(1) - pair % U2(1) * rij_hat(3)
            u2xrij(3) = pair % U2(1) * rij_hat(2) - pair % U2(2) * rij_hat(1)

            rijxu1(1) = rij_hat(2) * pair % U1(3) - rij_hat(3) * pair % U1(2)
            rijxu1(2) = rij_hat(3) * pair % U1(1) - rij_hat(1) * pair % U1(3)
            rijxu1(3) = rij_hat(1) * pair % U1(2) - rij_hat(2) * pair % U1(1)

            uur      = u1xu2(1) * rij_hat(1) + u1xu2(2) * rij_hat(2) + u1xu2(3) * rij_hat(3)
            uu       = param % uu

            dR_dr   = dR_dr  * srd_fac * uur * uu
            de2_dr  = de2_dr * srx_fac * uur * uu
            duur_dr = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uu * ( u1xu2 - uur * rij_hat ) / rij_mag

            de1_du1  = de1_du1 * srx_fac * uur * uu
            de2_du1  = de2_du1 * srx_fac * uur * uu
            dR_du1   = dR_du1  * srd_fac * uur * uu
            duur_du1 = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uu * u2xrij
            duu_du1  = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uur * pair % U2

            de1_du2  = de1_du2 * srx_fac * uur * uu
            de2_du2  = de2_du2 * srx_fac * uur * uu
            dR_du2   = dR_du2  * srd_fac * uur * uu
            duur_du2 = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uu * rijxu1
            duu_du2  = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uur * pair % U1

            pot_chiral   = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uur * uu
            FIJ_chiral   = (-1.0) * (          de2_dr  + dR_dr  + duur_dr           )
            TI_chiral    = (-1.0) * (de1_du1 + de2_du1 + dR_du1 + duur_du1 + duu_du1)
            TJ_chiral    = (-1.0) * (de1_du2 + de2_du2 + dR_du2 + duur_du2 + duu_du2)
            TORQ1_chiral(1) = pair % U1(2) * TI_chiral(3) - pair % U1(3) * TI_chiral(2)
            TORQ1_chiral(2) = pair % U1(3) * TI_chiral(1) - pair % U1(1) * TI_chiral(3)
            TORQ1_chiral(3) = pair % U1(1) * TI_chiral(2) - pair % U1(2) * TI_chiral(1)
            TORQ2_chiral(1) = pair % U2(2) * TJ_chiral(3) - pair % U2(3) * TJ_chiral(2)
            TORQ2_chiral(2) = pair % U2(3) * TJ_chiral(1) - pair % U2(1) * TJ_chiral(3)
            TORQ2_chiral(3) = pair % U2(1) * TJ_chiral(2) - pair % U2(2) * TJ_chiral(1)

            pot   = pot   + chirality * pot_chiral
            FIJ   = FIJ   + chirality * FIJ_chiral
            TORQ1 = TORQ1 + chirality * TORQ1_chiral
            TORQ2 = TORQ2 + chirality * TORQ2_chiral

            if (( rij_mag >= switch ) .and. ( rij_mag <= (sigma + cutoff - sigma_0) ) ) then
                FIJ   = sx * FIJ   - sdx * pot * factor
                TORQ1 = sx * TORQ1 - sdx * pot * factor
                TORQ2 = sx * TORQ2 - sdx * pot * factor
                pot   = sx * pot
            endif

        endif

        PE = pot
        FXIJ = FIJ
        TXI = TORQ1
        TXJ = TORQ2

    end subroutine gbc_calculate_forces

    pure subroutine gbwc_calculate_forces( I, J, cutoff, PE, FXIJ, TXI, TXJ, wf, chirality )
        implicit none
        integer, intent(in)                  :: I, J
        real(real64), intent(in)             :: cutoff, chirality, wf
        real(real64), intent(out)            :: PE, FXIJ(3), TXI(3), TXJ(3)
        
        real(real64), dimension(3)           :: rij_hat, u2xrij, rijxu1, u1xu2
        real(real64), dimension(3)           :: dR_dr, de2_dr, dR_du1, dR_du2
        real(real64), dimension(3)           :: de1_du1, de2_du1, de1_du2, de2_du2
        real(real64), dimension(3)           :: duur_dr, duur_du1, duur_du2, duu_du1, duu_du2
        real(real64), dimension(3)           :: FIJ, TI, TJ, TORQ1, TORQ2
        real(real64), dimension(3)           :: FIJ_chiral, TI_chiral, TJ_chiral, TORQ1_chiral, TORQ2_chiral
        
        real(real64)                         :: muinv, chi, xhi, uur, uu
        real(real64)                         :: lpar, lprp, epar, eprp
        real(real64)                         :: gchi, gxhi, ginv_chi_3by2
        real(real64)                         :: eps1, eps2, sigma, sigma_0, sigma_0_inv
        real(real64)                         :: sr_gb, sr2_gb, sr6_gb, sr12_gb
        real(real64)                         :: pot, pot_chiral, rij_mag,  srx_fac, srd_fac
        real(real64)                         :: switch, xx, sx, sdx, factor
        type(gb_params)                      :: param
        type(gb_pair)                        :: pair
        

        pot = 0.0; FIJ = 0.0; TORQ1 = 0.0; TORQ2 = 0.0; sigma_0 = 1.0;
        pot_chiral = 0.0; FIJ_chiral = 0.0; TORQ1_chiral = 0.0; TORQ2_chiral = 0.0

        pair % RIJ(1) = RX(I) - RX(J)
        pair % RIJ(2) = RY(I) - RY(J)
        pair % RIJ(3) = RZ(I) - RZ(J)

        if (periodic) then
            pair % RIJ = pair % RIJ - ANINT( pair % RIJ / box ) * box
        endif

        pair % rij_sq  = pair % RIJ(1) * pair % RIJ(1) + pair % RIJ(2) * pair % RIJ(2) + pair % RIJ(3) * pair % RIJ(3)
        rij_mag        = SQRT( pair % rij_sq )
        rij_hat        = pair % RIJ / rij_mag

        pair % U1(1) = 2.0 *  ( QX(I) * QZ(I) + QW(I) * QY(I) )
        pair % U1(2) = 2.0 *  ( QY(I) * QZ(I) - QW(I) * QX(I) )
        pair % U1(3) = QW(I)*QW(I) - QX(I)*QX(I) - QY(I)*QY(I) + QZ(I)*QZ(I)

        pair % U2(1) = 2.0 *  ( QX(J) * QZ(J) + QW(J) * QY(J) )
        pair % U2(2) = 2.0 *  ( QY(J) * QZ(J) - QW(J) * QX(J) )
        pair % U2(3) = QW(J)*QW(J) - QX(J)*QX(J) - QY(J)*QY(J) + QZ(J)*QZ(J)

        lpar = sig_z(I);  lprp = sig_x(I);
        epar = eps_z(I); eprp = eps_x(I);
        muinv = 1.0 / meu

        ! CALCULATE GAY BERNE PARAMETERS

        param % ru1 = pair % RIJ(1) * pair % U1(1) + pair % RIJ(2) * pair % U1(2) + pair % RIJ(3) * pair % U1(3)
        param % ru2 = pair % RIJ(1) * pair % U2(1) + pair % RIJ(2) * pair % U2(2) + pair % RIJ(3) * pair % U2(3)
        param % uu  = pair % U1(1)  * pair % U2(1) + pair % U1(2)  * pair % U2(2) + pair % U1(3)  * pair % U2(3)
        
        chi = ( lpar * lpar - lprp * lprp )     / ( lpar * lpar + lprp * lprp )
        xhi = ( eprp ** muinv - epar ** muinv ) / ( eprp ** muinv + epar ** muinv )
        gchi = g_func( chi, param, pair )
        gxhi = g_func( xhi, param, pair )
        
        eps1  = 1.0 / SQRT( 1 - ( chi * param % uu ) * ( chi * param % uu ) )
        eps2  = gxhi / eps_z(I)
        sigma = lprp / SQRT( gchi )
        sigma_0_inv = 1.0 / sigma_0
        ginv_chi_3by2 = sigma / gchi

        ! EVALUATE SMOOTHSTEP FUNCTION

        switch = 0.9 * ( sigma + cutoff + wf - sigma_0 )
        factor = 1.0 / ( ( sigma + cutoff + wf - sigma_0 ) - switch )
        xx = ( rij_mag - switch ) * factor
        sx = 1.0 - 6.0 * xx*xx*xx*xx*xx + 15.0 * xx*xx*xx*xx - 10.0 * xx*xx*xx
        sdx = -30.0 * xx*xx*xx*xx + 60.0 * xx*xx*xx - 30.0 * xx*xx

        ! EVALUATE GAY-BERNE POTENTIAL

        if ( rij_mag < ( sigma + cutoff + wf - sigma_0 ) ) then

            if ( rij_mag < sigma ) then
                sr_gb = sigma_0 / ( rij_mag - sigma + sigma_0 )
            else if ((rij_mag .ge. sigma) .and. (rij_mag .le. (sigma + wf))) then
                sr_gb = 1.0
            else
                sr_gb = sigma_0 / ( rij_mag - sigma + sigma_0 - wf )
            endif

            sr2_gb = sr_gb * sr_gb
            sr6_gb = sr2_gb * sr2_gb * sr2_gb
            sr12_gb = sr6_gb * sr6_gb

            srx_fac = ( sr12_gb - 2.0 * sr6_gb )
            srd_fac = 12.0 * ( sr6_gb - sr12_gb ) * sr_gb

            ! EVALUATE NORMAL GB POTENTIAL

            dR_dr   = (eps1 ** neu) * (eps2 ** meu) * sigma_0_inv * ( rij_hat + 0.5 * dG_dr( chi, param, pair ) * ginv_chi_3by2 )
            de2_dr  = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_dr( xhi, param, pair )

            de1_du1 = (eps2 ** meu) * (neu * eps1 ** (neu + 2)) * (chi * chi * param % uu * pair % U2)
            de2_du1 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_du1( xhi, param, pair )
            dR_du1  = (eps1 ** neu) * (eps2 ** meu) * 0.5 * sigma_0_inv * dG_du1( chi, param, pair ) * ginv_chi_3by2

            de1_du2 = (eps2 ** meu) * (neu * eps1 ** (neu + 2)) * (chi * chi * param % uu * pair % U1)
            de2_du2 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_du2( xhi, param, pair )
            dR_du2  = (eps1 ** neu) * (eps2 ** meu) * 0.5 * sigma_0_inv * dG_du2( chi, param, pair ) * ginv_chi_3by2
        
            pot = (eps1 ** neu) * (eps2 ** meu) * srx_fac
            FIJ = (-1.0) * ( srx_fac *  de2_dr             + srd_fac * dR_dr  )
            TI  = (-1.0) * ( srx_fac * (de1_du1 + de2_du1) + srd_fac * dR_du1 )
            TJ  = (-1.0) * ( srx_fac * (de1_du2 + de2_du2) + srd_fac * dR_du2 )
            TORQ1(1) = pair % U1(2) * TI(3) - pair % U1(3) * TI(2)
            TORQ1(2) = pair % U1(3) * TI(1) - pair % U1(1) * TI(3)
            TORQ1(3) = pair % U1(1) * TI(2) - pair % U1(2) * TI(1)
            TORQ2(1) = pair % U2(2) * TJ(3) - pair % U2(3) * TJ(2)
            TORQ2(2) = pair % U2(3) * TJ(1) - pair % U2(1) * TJ(3)
            TORQ2(3) = pair % U2(1) * TJ(2) - pair % U2(2) * TJ(1)

            ! EVALUATE CHIRAL TERMS

            srx_fac = ( sr12_gb - (12.0 / 7.0) * sr6_gb * sr_gb)
            srd_fac = 12.0 * ( sr6_gb * sr2_gb - sr12_gb * sr_gb )

            u1xu2(1) = pair % U1(2) * pair % U2(3) - pair % U1(3) * pair % U2(2)
            u1xu2(2) = pair % U1(3) * pair % U2(1) - pair % U1(1) * pair % U2(3)
            u1xu2(3) = pair % U1(1) * pair % U2(2) - pair % U1(2) * pair % U2(1) 

            u2xrij(1) = pair % U2(2) * rij_hat(3) - pair % U2(3) * rij_hat(2)
            u2xrij(2) = pair % U2(3) * rij_hat(1) - pair % U2(1) * rij_hat(3)
            u2xrij(3) = pair % U2(1) * rij_hat(2) - pair % U2(2) * rij_hat(1)

            rijxu1(1) = rij_hat(2) * pair % U1(3) - rij_hat(3) * pair % U1(2)
            rijxu1(2) = rij_hat(3) * pair % U1(1) - rij_hat(1) * pair % U1(3)
            rijxu1(3) = rij_hat(1) * pair % U1(2) - rij_hat(2) * pair % U1(1)

            uur      = u1xu2(1) * rij_hat(1) + u1xu2(2) * rij_hat(2) + u1xu2(3) * rij_hat(3)
            uu       = param % uu

            dR_dr   = dR_dr  * srd_fac * uur * uu
            de2_dr  = de2_dr * srx_fac * uur * uu
            duur_dr = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uu * ( u1xu2 - uur * rij_hat ) / rij_mag

            de1_du1  = de1_du1 * srx_fac * uur * uu
            de2_du1  = de2_du1 * srx_fac * uur * uu
            dR_du1   = dR_du1  * srd_fac * uur * uu
            duur_du1 = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uu * u2xrij
            duu_du1  = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uur * pair % U2

            de1_du2  = de1_du2 * srx_fac * uur * uu
            de2_du2  = de2_du2 * srx_fac * uur * uu
            dR_du2   = dR_du2  * srd_fac * uur * uu
            duur_du2 = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uu * rijxu1
            duu_du2  = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uur * pair % U1

            pot_chiral   = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uur * uu
            FIJ_chiral   = (-1.0) * (          de2_dr  + dR_dr  + duur_dr           )
            TI_chiral    = (-1.0) * (de1_du1 + de2_du1 + dR_du1 + duur_du1 + duu_du1)
            TJ_chiral    = (-1.0) * (de1_du2 + de2_du2 + dR_du2 + duur_du2 + duu_du2)
            TORQ1_chiral(1) = pair % U1(2) * TI_chiral(3) - pair % U1(3) * TI_chiral(2)
            TORQ1_chiral(2) = pair % U1(3) * TI_chiral(1) - pair % U1(1) * TI_chiral(3)
            TORQ1_chiral(3) = pair % U1(1) * TI_chiral(2) - pair % U1(2) * TI_chiral(1)
            TORQ2_chiral(1) = pair % U2(2) * TJ_chiral(3) - pair % U2(3) * TJ_chiral(2)
            TORQ2_chiral(2) = pair % U2(3) * TJ_chiral(1) - pair % U2(1) * TJ_chiral(3)
            TORQ2_chiral(3) = pair % U2(1) * TJ_chiral(2) - pair % U2(2) * TJ_chiral(1)

            pot   = pot   + chirality * pot_chiral
            FIJ   = FIJ   + chirality * FIJ_chiral
            TORQ1 = TORQ1 + chirality * TORQ1_chiral
            TORQ2 = TORQ2 + chirality * TORQ2_chiral

            if (( rij_mag >= switch ) .and. ( rij_mag <= (sigma + cutoff + wf - sigma_0) ) ) then
                FIJ   = sx * FIJ   - sdx * pot * factor
                TORQ1 = sx * TORQ1 - sdx * pot * factor
                TORQ2 = sx * TORQ2 - sdx * pot * factor
                pot   = sx * pot
            endif

        endif

        PE = pot
        FXIJ = FIJ
        TXI = TORQ1
        TXJ = TORQ2

    end subroutine gbwc_calculate_forces

end module gay_berne

module ecp
    use xmath
    use system
    implicit none
    integer, parameter, private :: meu = 2, neu = 1

    type :: ecp_pair
        real(real64), dimension(3, 3) :: AE, BE, AR, BR
        real(real64), dimension(3)    :: RIJ
    end type

    contains
    pure real(real64) function brent(func, pair) result (res)
        implicit none
        type(ecp_pair), intent(in)  :: pair
        integer, parameter          :: max_iter = 10000
        real(real64), parameter     :: eps = 1e-8, psi = 0.5 * ( 3.0 - sqrt(5.0) )
        real(real64)                :: a, b, c, x, w, v, u
        real(real64)                :: fa, fb, fc, ft, fw, fv, fu
        real(real64)                :: deltax, atol, rtol
        real(real64)                :: tol1, tol2, xmid, tmp1, tmp2, rat, p, dx_temp
        integer                     :: iter

        interface
            pure real(real64) function func(xx, pp)
                import ecp_pair, real64
                real(real64), intent(in)   :: xx
                type(ecp_pair), intent(in) :: pp
            end function func
        end interface

        deltax = 0.0; atol = 1e-3; rtol = 1e-2; iter = 0
        a = 0.0; b = 1.0 
        c = (1 - psi) * a + psi * b
        fa = func(a, pair); fb = func(b, pair); fc = func(c, pair)
        v = c; w = v; x = w
        fv = fc; fw = fv; ft = fw

        do while (iter < max_iter)
            tol1 = rtol * abs(x) + atol
            tol2 = 2.0 * tol1
            xmid = 0.5 * (a + b)

            if (abs(x - xmid) < (tol2 - 0.5 * (b - a))) exit

            if (abs(deltax) <= tol1) then
                if (x >= xmid) then
                    deltax = a - x
                else
                    deltax = b - x
                endif
                rat = psi * deltax

            else
                tmp1 = (x - w) * (ft - fv)
                tmp2 = (x - v) * (ft - fw)
                p = (x - v) * tmp2 - (x - w) * tmp1
                tmp2 = 2.0 * (tmp2 - tmp1)

                if (tmp2 > 0.0) p = -p

                tmp2 = abs(tmp2)
                dx_temp = deltax
                deltax = rat
                
                if ((p > tmp2 * (a - x)) .and. (p < tmp2 * (b - x)) .and. (abs(p) < abs(0.5 * tmp2 * dx_temp))) then
                    rat = p * 1.0 / tmp2
                    u = x + rat
                    if ((u - a) < tol2 .or. (b - u) < tol2) then
                        if (xmid - x >= 0) then
                            rat = tol1
                        else
                            rat = -tol1
                        endif
                    endif
                else
                    if (x >= xmid) then
                        deltax = a - x 
                    else
                        deltax = b - x
                    endif
                    rat = psi * deltax
                endif
            endif

            if (abs(rat) < tol1) then
                if (rat >= 0) then
                    u = x + tol1
                else
                    u = x - tol1
                endif
            else
                u = x + rat

            endif

            fu = func(u, pair)

            if (fu > ft) then

                if (u < x) then
                    a = u
                else
                    b = u
                endif

                if ((fu <= fw) .or. (w == x)) then
                    v = w
                    w = u
                    fv = fw
                    fw = fu

                else if ((fu <= fv) .or. (v == x) .or. (v == w)) then
                    v = u
                    fv = fu
                endif

            else

                if (u >= x) then
                    a = x
                else
                    b = x
                endif

                v = w
                w = x
                x = u
                fv = fw
                fw = ft
                ft = fu

            endif
            iter = iter + 1
    
        enddo

        res = x
    end function brent

    pure real(real64) function optim_eps( L, pair )
        implicit none
        real(real64), intent(in)        :: L
        type(ecp_pair), intent(in)      :: pair
        real(real64), dimension(3, 3)   :: GM
        real(real64), dimension(3)       :: KM

        GM = ( 1 - L ) * pair % AE + L * pair % BE
        KM = cholesky_solve( GM, pair % rij )

        optim_eps = -L * ( 1 - L ) * ( pair % rij(1) * KM(1) + pair % rij(2) * KM(2) + pair % rij(3) * KM(3) ) 
        
    end function optim_eps

    pure real(real64) function optim_dist( L, pair )
        implicit none
        real(real64), intent(in)         :: L
        type(ecp_pair), intent(in)  :: pair
        real(real64), dimension(3, 3)    :: GM_R
        real(real64), dimension(3)       :: KM_R

        GM_R = ( 1 - L ) * pair % AR + L * pair % BR
        KM_R = cholesky_solve( GM_R, pair % rij )

        optim_dist = -L * ( 1 - L ) * ( pair % rij(1) * KM_R(1) + pair % rij(2) * KM_R(2) + pair % rij(3) * KM_R(3) )
        
    end function optim_dist

    pure subroutine ecp_calculate_forces( I, J, cutoff, PE, FXIJ, TXI, TXJ )
        implicit none
        integer, intent(in)                  :: I, J
        real(real64), intent(in)             :: cutoff
        real(real64), intent(out)            :: PE, FXIJ(3), TXI(3), TXJ(3)
        
        real(real64), dimension(3, 3)        :: U1, U2, S1, S2, E1, E2, EM
        real(real64), dimension(3, 3)        :: U1_Inv, U2_Inv, S1_sq, S2_sq
        real(real64), dimension(3, 3)        :: GE, GR, MM1, MM2, EAB
        real(real64), dimension(3)           :: MM1x, MM1y, MM1z, MM2x, MM2y, MM2z
        real(real64), dimension(3)           :: EM_X_MM1x, EM_X_MM1y, EM_X_MM1z, EM_X_MM2x, EM_X_MM2y, EM_X_MM2z
        real(real64), dimension(3)           :: KE, KR, KE_X_AE, KR_X_AR, KE_X_BE, KR_X_BR
        real(real64)                         :: lambda_E, lambda_R
        real(real64)                         :: rij_sq, rij_mag, rij_hat(3)
        
        real(real64), dimension(3)           :: dR_dr, de2_dr
        real(real64), dimension(3)           :: de1_du1x, de1_du1y, de1_du1z, de1_du1, de2_du1, dR_du1
        real(real64), dimension(3)           :: de1_du2x, de1_du2y, de1_du2z, de1_du2, de2_du2, dR_du2
        real(real64), dimension(3)           :: FIJ, TORQ1, TORQ2
        
        real(real64)                         :: muinv
        real(real64)                         :: eps1, eps2, phi, phi_inv_3by2, sigma, sigma_0, sigma_0_inv
        real(real64)                         :: e10, e20, sigma_t1, sigma_t2, sigma_t
        real(real64)                         :: sr_ecp, sr2_ecp, sr6_ecp, sr12_ecp
        real(real64)                         :: pot, efac, rfac, srx_fac, srd_fac
        real(real64)                         :: switch, xx, sx, sdx, factor
        type(ecp_pair)                       :: pair

        pot = 0.0; FIJ = 0.0; TORQ1 = 0.0; TORQ2 = 0.0; sigma_0 = 1.0; 
        
        pair % RIJ(1) = RX(I) - RX(J)
        pair % RIJ(2) = RY(I) - RY(J)
        pair % RIJ(3) = RZ(I) - RZ(J)

        if (periodic) then
            pair % RIJ = pair % RIJ - ANINT( pair % RIJ / box ) * box
        endif

        rij_sq  = pair % RIJ(1)*pair % RIJ(1) + pair % RIJ(2)*pair % RIJ(2) + pair % RIJ(3)*pair % RIJ(3)
        rij_mag = SQRT( rij_sq )
        rij_hat = pair % RIJ / rij_mag

        U1(1, 1) = QW(I)*QW(I) + QX(I)*QX(I) - QY(I)*QY(I) - QZ(I)*QZ(I)
        U1(1, 2) = 2.0 * ( QX(I) * QY(I) + QW(I) * QZ(I) )
        U1(1, 3) = 2.0 * ( QX(I) * QZ(I) - QW(I) * QY(I) )
        U1(2, 1) = 2.0 * ( QX(I) * QY(I) - QW(I) * QZ(I) )
        U1(2, 2) = QW(I)*QW(I) - QX(I)*QX(I) + QY(I)*QY(I) - QZ(I)*QZ(I)
        U1(2, 3) = 2.0 * ( QY(I) * QZ(I) + QW(I) * QX(I) )
        U1(3, 1) = 2.0 * ( QX(I) * QZ(I) + QW(I) * QY(I) )
        U1(3, 2) = 2.0 * ( QY(I) * QZ(I) - QW(I) * QX(I) )
        U1(3, 3) = QW(I)*QW(I) - QX(I)*QX(I) - QY(I)*QY(I) + QZ(I)*QZ(I)

        U2(1, 1) = QW(J)*QW(J) + QX(J)*QX(J) - QY(J)*QY(J) - QZ(J)*QZ(J)
        U2(1, 2) = 2.0 * ( QX(J) * QY(J) + QW(J) * QZ(J) )
        U2(1, 3) = 2.0 * ( QX(J) * QZ(J) - QW(J) * QY(J) )
        U2(2, 1) = 2.0 * ( QX(J) * QY(J) - QW(J) * QZ(J) )
        U2(2, 2) = QW(J)*QW(J) - QX(J)*QX(J) + QY(J)*QY(J) - QZ(J)*QZ(J)
        U2(2, 3) = 2.0 * ( QY(J) * QZ(J) + QW(J) * QX(J) )
        U2(3, 1) = 2.0 * ( QX(J) * QZ(J) + QW(J) * QY(J) )
        U2(3, 2) = 2.0 * ( QY(J) * QZ(J) - QW(J) * QX(J) )
        U2(3, 3) = QW(J)*QW(J) - QX(J)*QX(J) - QY(J)*QY(J) + QZ(J)*QZ(J)

        e10 = max( eps_x(I), eps_y(I), eps_z(I) )
        e20 = max( eps_x(J), eps_y(J), eps_z(J) )
        muinv = 1.0 / meu
        E1 = 0.0_8; E2 = 0.0_8; S1 = 0.0_8; S2 = 0.0_8

        E1(1, 1) = ( (e10 / eps_x(I)) ** muinv ) / 4.0
        E1(2, 2) = ( (e10 / eps_y(I)) ** muinv ) / 4.0
        E1(3, 3) = ( (e10 / eps_z(I)) ** muinv ) / 4.0

        E2(1, 1) = ( (e20 / eps_x(J)) ** muinv ) / 4.0
        E2(2, 2) = ( (e20 / eps_y(J)) ** muinv ) / 4.0
        E2(3, 3) = ( (e20 / eps_z(J)) ** muinv ) / 4.0

        S1(1, 1) = ( sig_x(I) / 2.0 )
        S1(2, 2) = ( sig_y(I) / 2.0 )
        S1(3, 3) = ( sig_z(I) / 2.0 )

        S2(1, 1) = ( sig_x(J) / 2.0 )
        S2(2, 2) = ( sig_y(J) / 2.0 )
        S2(3, 3) = ( sig_z(J) / 2.0 )

        U1_Inv = transpose(U1)
        U2_Inv = transpose(U2)
        S1_sq = S1*S1
        S2_sq = S2*S2
        
        MM1 = U1_Inv .x. S1
        MM2 = U2_Inv .x. S2
        MM1x = MM1(:, 1); MM1y = MM1(:, 2); MM1z = MM1(:, 3)
        MM2x = MM2(:, 1); MM2y = MM2(:, 2); MM2z = MM2(:, 3)
        pair % AR = matmul_atba(U1, S1_sq)
        pair % BR = matmul_atba(U2, S2_sq)
        pair % AE = matmul_atba(U1, E1)
        pair % BE = matmul_atba(U2, E2)
        
        sigma_t1 = (sig_x(I) * sig_y(I) + sig_z(I)*sig_z(I)) * sqrt( 2.0 * sig_x(I) * sig_y(I) ) / 8.0
        sigma_t2 = (sig_x(J) * sig_y(J) + sig_z(J)*sig_z(J)) * sqrt( 2.0 * sig_x(J) * sig_y(J) ) / 8.0
        sigma_t  = sqrt( sigma_t1 * sigma_t2 )
        
        EAB  = pair % AR + pair % BR
        EM   = inverse( EAB )
        eps1 = sigma_t * sqrt( determinant( EM ) )
        
        lambda_R = brent( optim_dist, pair )
        lambda_E = lambda_R

        GE   = ( 1 - lambda_E ) * pair % AE + lambda_E * pair % BE
        KE   = cholesky_solve( GE, pair % RIJ )
        eps2 = lambda_E * (1.0 - lambda_E) * ( rij_hat(1) * KE(1) + rij_hat(2) * KE(2) + rij_hat(3) * KE(3) ) / rij_mag
        
        GR    = ( 1 - lambda_R ) * pair % AR + lambda_R * pair % BR
        KR    = cholesky_solve( GR, pair % RIJ )
        phi   = lambda_R * (1.0 - lambda_R) * ( rij_hat(1) * KR(1) + rij_hat(2) * KR(2) + rij_hat(3) * KR(3)) / rij_mag
        sigma = 1.0 / sqrt( phi )
        efac  = 2.0 * lambda_E * ( 1.0 - lambda_E ) / rij_sq
        rfac  = 2.0 * lambda_R * ( 1.0 - lambda_R ) / rij_sq
        sigma_0_inv = 1.0 / sigma_0
        phi_inv_3by2 = sigma * sigma * sigma

        ! EVALUATE SMOOTHSTEP FUNCTION

        switch = 0.9 * ( sigma + cutoff - sigma_0 )
        factor = 1.0 / ( ( sigma + cutoff - sigma_0 ) - switch )
        xx = ( rij_mag - switch ) * factor
        sx = 1.0 - 6.0 * xx*xx*xx*xx*xx + 15.0 * xx*xx*xx*xx - 10.0 * xx*xx*xx
        sdx = -30.0 * xx*xx*xx*xx + 60.0 * xx*xx*xx - 30.0 * xx*xx

        ! EVALUATE ELLIPSOID CONTACT POTENTIAL

        if ( rij_mag < ( sigma + cutoff - sigma_0 ) ) then

            sr_ecp   = sigma_0 / ( rij_mag - sigma + sigma_0 )
            sr2_ecp  = sr_ecp * sr_ecp
            sr6_ecp  = sr2_ecp * sr2_ecp * sr2_ecp
            sr12_ecp = sr6_ecp * sr6_ecp

            srx_fac = ( sr12_ecp - 2.0 * sr6_ecp )
            srd_fac = 12.0 * ( sr6_ecp - sr12_ecp ) * sr_ecp

            dR_dr  = sigma_0_inv * ( rij_hat + 0.5 * rfac * ( KR - dot_product( rij_hat, KR ) * rij_hat ) * phi_inv_3by2 )
            de2_dr = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * efac * ( KE - dot_product( rij_hat, KE ) * rij_hat )

            EM_X_MM1x = EM .x. MM1x; EM_X_MM1y = EM .x. MM1y; EM_X_MM1z = EM .x. MM1z
            de1_du1x  = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM1x, EM_X_MM1x )
            de1_du1y  = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM1y, EM_X_MM1y )
            de1_du1z  = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM1z, EM_X_MM1z )

            KE_X_AE = KE .x. pair % AE; KR_X_AR = KR .x. pair % AR
            de1_du1 = (de1_du1x + de1_du1y + de1_du1z)
            de2_du1 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * srx_fac * (1 - lambda_E) * efac * cross( KE_X_AE, KE )
            dR_du1  = (eps1 ** neu) * (eps2 ** meu) * srd_fac * (0.5 * sigma_0_inv) * (1 - lambda_R) * rfac * cross( KR_X_AR, KR ) * phi_inv_3by2

            EM_X_MM2x = EM .x. MM2x; EM_X_MM2y = EM .x. MM2y; EM_X_MM2z = EM .x. MM2z
            de1_du2x  = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM2x, EM_X_MM2x )
            de1_du2y  = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM2y, EM_X_MM2y )
            de1_du2z  = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM2z, EM_X_MM2z )

            KE_X_BE = KE .x. pair % BE; KR_X_BR = KR .x. pair % BR
            de1_du2 = (de1_du2x + de1_du2y + de1_du2z)
            de2_du2 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * srx_fac * lambda_E * efac * cross( KE_X_BE, KE )
            dR_du2  = (eps1 ** neu) * (eps2 ** meu) * srd_fac * (0.5 * sigma_0_inv) * lambda_R * rfac * cross( KR_X_BR, KR ) * phi_inv_3by2
            
            pot   = (eps1 ** neu) * (eps2 ** meu) * srx_fac
            FIJ   = (-1.0) * ((eps1 ** neu) * (eps2 ** meu) * srd_fac * dR_dr  + srx_fac * de2_dr)
            TORQ1 = (de1_du1 + de2_du1 + dR_du1)
            TORQ2 = (de1_du2 + de2_du2 + dR_du2)

            if (( rij_mag >= switch ) .and. ( rij_mag <= (sigma + cutoff - sigma_0) )) then
                FIJ   = sx * FIJ   - sdx * pot * factor
                TORQ1 = sx * TORQ1 - sdx * pot * factor
                TORQ2 = sx * TORQ2 - sdx * pot * factor
                pot   = sx * pot
            endif

        endif
        
        PE = pot
        FXIJ = FIJ
        TXI = TORQ1
        TXJ = TORQ2
        
    end subroutine ecp_calculate_forces

    pure subroutine ecpr_calculate_forces( I, J, PE, FXIJ, TXI, TXJ )
        implicit none
        integer, intent(in)                  :: I, J
        real(real64), intent(out)            :: PE, FXIJ(3), TXI(3), TXJ(3)
        
        real(real64), dimension(3, 3)        :: U1, U2, S1, S2, E1, E2, EM
        real(real64), dimension(3, 3)        :: U1_Inv, U2_Inv, S1_sq, S2_sq
        real(real64), dimension(3, 3)        :: GE, GR, MM1, MM2, EAB
        real(real64), dimension(3)           :: MM1x, MM1y, MM1z, MM2x, MM2y, MM2z
        real(real64), dimension(3)           :: EM_X_MM1x, EM_X_MM1y, EM_X_MM1z, EM_X_MM2x, EM_X_MM2y, EM_X_MM2z
        real(real64), dimension(3)           :: KE, KR, KE_X_AE, KR_X_AR, KE_X_BE, KR_X_BR
        real(real64)                         :: lambda_E, lambda_R
        real(real64)                         :: rij_sq, rij_mag, rij_hat(3)
        
        real(real64), dimension(3)           :: dR_dr, de2_dr
        real(real64), dimension(3)           :: de1_du1x, de1_du1y, de1_du1z, de1_du1, de2_du1, dR_du1
        real(real64), dimension(3)           :: de1_du2x, de1_du2y, de1_du2z, de1_du2, de2_du2, dR_du2
        real(real64), dimension(3)           :: FIJ, TORQ1, TORQ2
        
        real(real64)                         :: muinv
        real(real64)                         :: eps1, eps2, phi, phi_inv_3by2, sigma, sigma_0, sigma_0_inv
        real(real64)                         :: e10, e20, sigma_t1, sigma_t2, sigma_t
        real(real64)                         :: sr_ecp, sr2_ecp, sr6_ecp, sr12_ecp
        real(real64)                         :: pot, efac, rfac, srx_fac, srd_fac
        type(ecp_pair)                       :: pair

        pot = 0.0; FIJ = 0.0; TORQ1 = 0.0; TORQ2 = 0.0; sigma_0 = 1.0; 
        
        pair % RIJ(1) = RX(I) - RX(J)
        pair % RIJ(2) = RY(I) - RY(J)
        pair % RIJ(3) = RZ(I) - RZ(J)

        if (periodic) then
            pair % RIJ = pair % RIJ - ANINT( pair % RIJ / box ) * box
        endif

        rij_sq = pair % RIJ(1)*pair % RIJ(1) + pair % RIJ(2)*pair % RIJ(2) + pair % RIJ(3)*pair % RIJ(3)
        rij_mag = SQRT( rij_sq )
        rij_hat = pair % RIJ / rij_mag

        U1(1, 1) = QW(I)*QW(I) + QX(I)*QX(I) - QY(I)*QY(I) - QZ(I)*QZ(I)
        U1(1, 2) = 2.0 * ( QX(I) * QY(I) + QW(I) * QZ(I) )
        U1(1, 3) = 2.0 * ( QX(I) * QZ(I) - QW(I) * QY(I) )
        U1(2, 1) = 2.0 * ( QX(I) * QY(I) - QW(I) * QZ(I) )
        U1(2, 2) = QW(I)*QW(I) - QX(I)*QX(I) + QY(I)*QY(I) - QZ(I)*QZ(I)
        U1(2, 3) = 2.0 * ( QY(I) * QZ(I) + QW(I) * QX(I) )
        U1(3, 1) = 2.0 * ( QX(I) * QZ(I) + QW(I) * QY(I) )
        U1(3, 2) = 2.0 * ( QY(I) * QZ(I) - QW(I) * QX(I) )
        U1(3, 3) = QW(I)*QW(I) - QX(I)*QX(I) - QY(I)*QY(I) + QZ(I)*QZ(I)

        U2(1, 1) = QW(J)*QW(J) + QX(J)*QX(J) - QY(J)*QY(J) - QZ(J)*QZ(J)
        U2(1, 2) = 2.0 * ( QX(J) * QY(J) + QW(J) * QZ(J) )
        U2(1, 3) = 2.0 * ( QX(J) * QZ(J) - QW(J) * QY(J) )
        U2(2, 1) = 2.0 * ( QX(J) * QY(J) - QW(J) * QZ(J) )
        U2(2, 2) = QW(J)*QW(J) - QX(J)*QX(J) + QY(J)*QY(J) - QZ(J)*QZ(J)
        U2(2, 3) = 2.0 * ( QY(J) * QZ(J) + QW(J) * QX(J) )
        U2(3, 1) = 2.0 * ( QX(J) * QZ(J) + QW(J) * QY(J) )
        U2(3, 2) = 2.0 * ( QY(J) * QZ(J) - QW(J) * QX(J) )
        U2(3, 3) = QW(J)*QW(J) - QX(J)*QX(J) - QY(J)*QY(J) + QZ(J)*QZ(J)

        e10 = max( eps_x(I), eps_y(I), eps_z(I) )
        e20 = max( eps_x(J), eps_y(J), eps_z(J) )
        muinv = 1.0 / meu
        E1 = 0.0_8; E2 = 0.0_8; S1 = 0.0_8; S2 = 0.0_8

        E1(1, 1) = ( (e10 / eps_x(I)) ** muinv ) / 4.0
        E1(2, 2) = ( (e10 / eps_y(I)) ** muinv ) / 4.0
        E1(3, 3) = ( (e10 / eps_z(I)) ** muinv ) / 4.0

        E2(1, 1) = ( (e20 / eps_x(J)) ** muinv ) / 4.0
        E2(2, 2) = ( (e20 / eps_y(J)) ** muinv ) / 4.0
        E2(3, 3) = ( (e20 / eps_z(J)) ** muinv ) / 4.0

        S1(1, 1) = ( sig_x(I) / 2.0 )
        S1(2, 2) = ( sig_y(I) / 2.0 )
        S1(3, 3) = ( sig_z(I) / 2.0 )

        S2(1, 1) = ( sig_x(J) / 2.0 )
        S2(2, 2) = ( sig_y(J) / 2.0 )
        S2(3, 3) = ( sig_z(J) / 2.0 )

        U1_Inv = transpose(U1)
        U2_Inv = transpose(U2)
        S1_sq = S1*S1
        S2_sq = S2*S2
        
        MM1 = U1_Inv .x. S1
        MM2 = U2_Inv .x. S2
        MM1x = MM1(:, 1); MM1y = MM1(:, 2); MM1z = MM1(:, 3)
        MM2x = MM2(:, 1); MM2y = MM2(:, 2); MM2z = MM2(:, 3)
        pair % AR = matmul_atba(U1, S1_sq)
        pair % BR = matmul_atba(U2, S2_sq)
        pair % AE = matmul_atba(U1, E1)
        pair % BE = matmul_atba(U2, E2)
        
        sigma_t1 = (sig_x(I) * sig_y(I) + sig_z(I)*sig_z(I)) * sqrt( 2.0 * sig_x(I) * sig_y(I) ) / 8.0
        sigma_t2 = (sig_x(J) * sig_y(J) + sig_z(J)*sig_z(J)) * sqrt( 2.0 * sig_x(J) * sig_y(J) ) / 8.0
        sigma_t = sqrt( sigma_t1 * sigma_t2 )
        
        EAB = pair % AR + pair % BR
        EM = inverse( EAB )
        eps1 = sigma_t * sqrt( determinant( EM ) )
        
        lambda_R = brent( optim_dist, pair )
        lambda_E = lambda_R

        GE = ( 1 - lambda_E ) * pair % AE + lambda_E * pair % BE
        KE = cholesky_solve( GE, pair % RIJ )
        eps2 = lambda_E * (1.0 - lambda_E) * ( rij_hat(1) * KE(1) + rij_hat(2) * KE(2) + rij_hat(3) * KE(3) ) / rij_mag
        
        GR    = ( 1 - lambda_R ) * pair % AR + lambda_R * pair % BR
        KR    = cholesky_solve( GR, pair % RIJ )
        phi   = lambda_R * (1.0 - lambda_R) * ( rij_hat(1) * KR(1) + rij_hat(2) * KR(2) + rij_hat(3) * KR(3)) / rij_mag
        sigma = 1.0 / sqrt( phi )
        efac  = 2.0 * lambda_E * ( 1.0 - lambda_E ) / rij_sq
        rfac  = 2.0 * lambda_R * ( 1.0 - lambda_R ) / rij_sq
        sigma_0_inv = 1.0 / sigma_0
        phi_inv_3by2 = sigma * sigma * sigma

        if ( rij_mag <= sigma ) then

            sr_ecp = sigma_0 / ( rij_mag - sigma + sigma_0 )
            sr2_ecp = sr_ecp * sr_ecp
            sr6_ecp = sr2_ecp * sr2_ecp * sr2_ecp
            sr12_ecp = sr6_ecp * sr6_ecp

            srx_fac = ( sr12_ecp - 2.0 * sr6_ecp ) + 1.0
            srd_fac = 12.0 * ( sr6_ecp - sr12_ecp ) * sr_ecp

            dR_dr  = sigma_0_inv * ( rij_hat + 0.5 * rfac * ( KR - dot_product( rij_hat, KR ) * rij_hat ) * phi_inv_3by2 )
            de2_dr = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * efac * ( KE - dot_product( rij_hat, KE ) * rij_hat )

            EM_X_MM1x = EM .x. MM1x; EM_X_MM1y = EM .x. MM1y; EM_X_MM1z = EM .x. MM1z
            de1_du1x  = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM1x, EM_X_MM1x )
            de1_du1y  = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM1y, EM_X_MM1y )
            de1_du1z  = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM1z, EM_X_MM1z )

            KE_X_AE = KE .x. pair % AE; KR_X_AR = KR .x. pair % AR
            de1_du1 = (de1_du1x + de1_du1y + de1_du1z)
            de2_du1 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * srx_fac * (1 - lambda_E) * efac * cross( KE_X_AE, KE )
            dR_du1  = (eps1 ** neu) * (eps2 ** meu) * srd_fac * (0.5 * sigma_0_inv) * (1 - lambda_R) * rfac * cross( KR_X_AR, KR ) * phi_inv_3by2

            EM_X_MM2x = EM .x. MM2x; EM_X_MM2y = EM .x. MM2y; EM_X_MM2z = EM .x. MM2z
            de1_du2x  = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM2x, EM_X_MM2x )
            de1_du2y  = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM2y, EM_X_MM2y )
            de1_du2z  = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM2z, EM_X_MM2z )

            KE_X_BE = KE .x. pair % BE; KR_X_BR = KR .x. pair % BR
            de1_du2 = (de1_du2x + de1_du2y + de1_du2z)
            de2_du2 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * srx_fac * lambda_E * efac * cross( KE_X_BE, KE )
            dR_du2  = (eps1 ** neu) * (eps2 ** meu) * srd_fac * (0.5 * sigma_0_inv) * lambda_R * rfac * cross( KR_X_BR, KR ) * phi_inv_3by2
            
            pot   = (eps1 ** neu) * (eps2 ** meu) * srx_fac
            FIJ   = (-1.0) * ((eps1 ** neu) * (eps2 ** meu) * srd_fac * dR_dr  + srx_fac * de2_dr)
            TORQ1 = (de1_du1 + de2_du1 + dR_du1)
            TORQ2 = (de1_du2 + de2_du2 + dR_du2)

        endif
        
        PE = pot
        FXIJ = FIJ
        TXI = TORQ1
        TXJ = TORQ2

    end subroutine ecpr_calculate_forces

    pure subroutine ecpw_calculate_forces( I, J, cutoff, PE, FXIJ, TXI, TXJ, wf )
        implicit none
        integer, intent(in)                  :: I, J
        real(real64), intent(in)             :: cutoff, wf
        real(real64), intent(out)            :: PE, FXIJ(3), TXI(3), TXJ(3)
        
        real(real64), dimension(3, 3)        :: U1, U2, S1, S2, E1, E2, EM
        real(real64), dimension(3, 3)        :: U1_Inv, U2_Inv, S1_sq, S2_sq
        real(real64), dimension(3, 3)        :: GE, GR, MM1, MM2, EAB
        real(real64), dimension(3)           :: MM1x, MM1y, MM1z, MM2x, MM2y, MM2z
        real(real64), dimension(3)           :: EM_X_MM1x, EM_X_MM1y, EM_X_MM1z, EM_X_MM2x, EM_X_MM2y, EM_X_MM2z
        real(real64), dimension(3)           :: KE, KR, KE_X_AE, KR_X_AR, KE_X_BE, KR_X_BR
        real(real64)                         :: lambda_E, lambda_R
        real(real64)                         :: rij_sq, rij_mag, rij_hat(3)
        
        real(real64), dimension(3)           :: dR_dr, de2_dr
        real(real64), dimension(3)           :: de1_du1x, de1_du1y, de1_du1z, de1_du1, de2_du1, dR_du1
        real(real64), dimension(3)           :: de1_du2x, de1_du2y, de1_du2z, de1_du2, de2_du2, dR_du2
        real(real64), dimension(3)           :: FIJ, TORQ1, TORQ2
        
        real(real64)                         :: muinv
        real(real64)                         :: eps1, eps2, phi, phi_inv_3by2, sigma, sigma_0, sigma_0_inv
        real(real64)                         :: e10, e20, sigma_t1, sigma_t2, sigma_t
        real(real64)                         :: sr_ecp, sr2_ecp, sr6_ecp, sr12_ecp
        real(real64)                         :: pot, efac, rfac, srx_fac, srd_fac
        real(real64)                         :: switch, xx, sx, sdx, factor
        type(ecp_pair)                       :: pair

        pot = 0.0; FIJ = 0.0; TORQ1 = 0.0; TORQ2 = 0.0; sigma_0 = 1.0; 
        
        pair % RIJ(1) = RX(I) - RX(J)
        pair % RIJ(2) = RY(I) - RY(J)
        pair % RIJ(3) = RZ(I) - RZ(J)

        if (periodic) then
            pair % RIJ = pair % RIJ - ANINT( pair % RIJ / box ) * box
        endif

        rij_sq = pair % RIJ(1)*pair % RIJ(1) + pair % RIJ(2)*pair % RIJ(2) + pair % RIJ(3)*pair % RIJ(3)
        rij_mag = SQRT( rij_sq )
        rij_hat = pair % RIJ / rij_mag

        U1(1, 1) = QW(I)*QW(I) + QX(I)*QX(I) - QY(I)*QY(I) - QZ(I)*QZ(I)
        U1(1, 2) = 2.0 * ( QX(I) * QY(I) + QW(I) * QZ(I) )
        U1(1, 3) = 2.0 * ( QX(I) * QZ(I) - QW(I) * QY(I) )
        U1(2, 1) = 2.0 * ( QX(I) * QY(I) - QW(I) * QZ(I) )
        U1(2, 2) = QW(I)*QW(I) - QX(I)*QX(I) + QY(I)*QY(I) - QZ(I)*QZ(I)
        U1(2, 3) = 2.0 * ( QY(I) * QZ(I) + QW(I) * QX(I) )
        U1(3, 1) = 2.0 * ( QX(I) * QZ(I) + QW(I) * QY(I) )
        U1(3, 2) = 2.0 * ( QY(I) * QZ(I) - QW(I) * QX(I) )
        U1(3, 3) = QW(I)*QW(I) - QX(I)*QX(I) - QY(I)*QY(I) + QZ(I)*QZ(I)

        U2(1, 1) = QW(J)*QW(J) + QX(J)*QX(J) - QY(J)*QY(J) - QZ(J)*QZ(J)
        U2(1, 2) = 2.0 * ( QX(J) * QY(J) + QW(J) * QZ(J) )
        U2(1, 3) = 2.0 * ( QX(J) * QZ(J) - QW(J) * QY(J) )
        U2(2, 1) = 2.0 * ( QX(J) * QY(J) - QW(J) * QZ(J) )
        U2(2, 2) = QW(J)*QW(J) - QX(J)*QX(J) + QY(J)*QY(J) - QZ(J)*QZ(J)
        U2(2, 3) = 2.0 * ( QY(J) * QZ(J) + QW(J) * QX(J) )
        U2(3, 1) = 2.0 * ( QX(J) * QZ(J) + QW(J) * QY(J) )
        U2(3, 2) = 2.0 * ( QY(J) * QZ(J) - QW(J) * QX(J) )
        U2(3, 3) = QW(J)*QW(J) - QX(J)*QX(J) - QY(J)*QY(J) + QZ(J)*QZ(J)

        e10 = max( eps_x(I), eps_y(I), eps_z(I) )
        e20 = max( eps_x(J), eps_y(J), eps_z(J) )
        muinv = 1.0 / meu
        E1 = 0.0_8; E2 = 0.0_8; S1 = 0.0_8; S2 = 0.0_8

        E1(1, 1) = ( (e10 / eps_x(I)) ** muinv ) / 4.0
        E1(2, 2) = ( (e10 / eps_y(I)) ** muinv ) / 4.0
        E1(3, 3) = ( (e10 / eps_z(I)) ** muinv ) / 4.0

        E2(1, 1) = ( (e20 / eps_x(J)) ** muinv ) / 4.0
        E2(2, 2) = ( (e20 / eps_y(J)) ** muinv ) / 4.0
        E2(3, 3) = ( (e20 / eps_z(J)) ** muinv ) / 4.0

        S1(1, 1) = ( sig_x(I) / 2.0 )
        S1(2, 2) = ( sig_y(I) / 2.0 )
        S1(3, 3) = ( sig_z(I) / 2.0 )

        S2(1, 1) = ( sig_x(J) / 2.0 )
        S2(2, 2) = ( sig_y(J) / 2.0 )
        S2(3, 3) = ( sig_z(J) / 2.0 )

        U1_Inv = transpose(U1)
        U2_Inv = transpose(U2)
        S1_sq = S1*S1
        S2_sq = S2*S2
        
        MM1 = U1_Inv .x. S1
        MM2 = U2_Inv .x. S2
        MM1x = MM1(:, 1); MM1y = MM1(:, 2); MM1z = MM1(:, 3)
        MM2x = MM2(:, 1); MM2y = MM2(:, 2); MM2z = MM2(:, 3)
        pair % AR = matmul_atba(U1, S1_sq)
        pair % BR = matmul_atba(U2, S2_sq)
        pair % AE = matmul_atba(U1, E1)
        pair % BE = matmul_atba(U2, E2)
        
        sigma_t1 = (sig_x(I) * sig_y(I) + sig_z(I)*sig_z(I)) * sqrt( 2.0 * sig_x(I) * sig_y(I) ) / 8.0
        sigma_t2 = (sig_x(J) * sig_y(J) + sig_z(J)*sig_z(J)) * sqrt( 2.0 * sig_x(J) * sig_y(J) ) / 8.0
        sigma_t = sqrt( sigma_t1 * sigma_t2 )
        
        EAB = pair % AR + pair % BR
        EM = inverse( EAB )
        eps1 = sigma_t * sqrt( determinant( EM ) )
        
        lambda_R = brent( optim_dist, pair )
        lambda_E = lambda_R

        GE = ( 1 - lambda_E ) * pair % AE + lambda_E * pair % BE
        KE = cholesky_solve( GE, pair % RIJ )
        eps2 = lambda_E * (1.0 - lambda_E) * ( rij_hat(1) * KE(1) + rij_hat(2) * KE(2) + rij_hat(3) * KE(3) ) / rij_mag
        
        GR    = ( 1 - lambda_R ) * pair % AR + lambda_R * pair % BR
        KR    = cholesky_solve( GR, pair % RIJ )
        phi   = lambda_R * (1.0 - lambda_R) * ( rij_hat(1) * KR(1) + rij_hat(2) * KR(2) + rij_hat(3) * KR(3)) / rij_mag
        sigma = 1.0 / sqrt( phi )
        efac  = 2.0 * lambda_E * ( 1.0 - lambda_E ) / rij_sq
        rfac  = 2.0 * lambda_R * ( 1.0 - lambda_R ) / rij_sq
        sigma_0_inv = 1.0 / sigma_0
        phi_inv_3by2 = sigma * sigma * sigma

        ! EVALUATE SMOOTHSTEP FUNCTION

        switch = 0.9 * ( sigma + cutoff + wf - sigma_0 )
        factor = 1.0 / ( ( sigma + cutoff + wf - sigma_0 ) - switch )
        xx = ( rij_mag - switch ) * factor
        sx = 1.0 - 6.0 * xx*xx*xx*xx*xx + 15.0 * xx*xx*xx*xx - 10.0 * xx*xx*xx
        sdx = -30.0 * xx*xx*xx*xx + 60.0 * xx*xx*xx - 30.0 * xx*xx

        ! EVALUATE ELLIPSOID CONTACT POTENTIAL

        if ( rij_mag < ( sigma + cutoff + wf - sigma_0 ) ) then

            if ( rij_mag < sigma ) then
                sr_ecp = sigma_0 / ( rij_mag - sigma + sigma_0 )
            else if ((rij_mag .ge. sigma) .and. (rij_mag .le. (sigma + wf))) then
                sr_ecp = 1.0
            else
                sr_ecp = sigma_0 / ( rij_mag - sigma + sigma_0 - wf )
            endif

            sr2_ecp = sr_ecp * sr_ecp
            sr6_ecp = sr2_ecp * sr2_ecp * sr2_ecp
            sr12_ecp = sr6_ecp * sr6_ecp

            srx_fac = ( sr12_ecp - 2.0 * sr6_ecp )
            srd_fac = 12.0 * ( sr6_ecp - sr12_ecp ) * sr_ecp

            dR_dr  = sigma_0_inv * ( rij_hat + 0.5 * rfac * ( KR - dot_product( rij_hat, KR ) * rij_hat ) * phi_inv_3by2 )
            de2_dr = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * efac * ( KE - dot_product( rij_hat, KE ) * rij_hat )

            EM_X_MM1x = EM .x. MM1x; EM_X_MM1y = EM .x. MM1y; EM_X_MM1z = EM .x. MM1z
            de1_du1x  = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM1x, EM_X_MM1x )
            de1_du1y  = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM1y, EM_X_MM1y )
            de1_du1z  = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM1z, EM_X_MM1z )

            KE_X_AE = KE .x. pair % AE; KR_X_AR = KR .x. pair % AR
            de1_du1 = (de1_du1x + de1_du1y + de1_du1z)
            de2_du1 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * srx_fac * (1 - lambda_E) * efac * cross( KE_X_AE, KE )
            dR_du1  = (eps1 ** neu) * (eps2 ** meu) * srd_fac * (0.5 * sigma_0_inv) * (1 - lambda_R) * rfac * cross( KR_X_AR, KR ) * phi_inv_3by2

            EM_X_MM2x = EM .x. MM2x; EM_X_MM2y = EM .x. MM2y; EM_X_MM2z = EM .x. MM2z
            de1_du2x  = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM2x, EM_X_MM2x )
            de1_du2y  = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM2y, EM_X_MM2y )
            de1_du2z  = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM2z, EM_X_MM2z )

            KE_X_BE = KE .x. pair % BE; KR_X_BR = KR .x. pair % BR
            de1_du2 = (de1_du2x + de1_du2y + de1_du2z)
            de2_du2 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * srx_fac * lambda_E * efac * cross( KE_X_BE, KE )
            dR_du2  = (eps1 ** neu) * (eps2 ** meu) * srd_fac * (0.5 * sigma_0_inv) * lambda_R * rfac * cross( KR_X_BR, KR ) * phi_inv_3by2
            
            pot   = (eps1 ** neu) * (eps2 ** meu) * srx_fac
            FIJ   = (-1.0) * ((eps1 ** neu) * (eps2 ** meu) * srd_fac * dR_dr  + srx_fac * de2_dr)
            TORQ1 = (de1_du1 + de2_du1 + dR_du1)
            TORQ2 = (de1_du2 + de2_du2 + dR_du2)

            if (( rij_mag .ge. switch ) .and. ( rij_mag <= (sigma + cutoff + wf - sigma_0) )) then
                FIJ   = sx * FIJ   - sdx * pot * factor
                TORQ1 = sx * TORQ1 - sdx * pot * factor
                TORQ2 = sx * TORQ2 - sdx * pot * factor
                pot   = sx * pot
            endif

        endif
        
        PE = pot
        FXIJ = FIJ
        TXI = TORQ1
        TXJ = TORQ2
    end subroutine ecpw_calculate_forces

    pure subroutine ecpc_calculate_forces( I, J, cutoff, PE, FXIJ, TXI, TXJ, chirality )
        implicit none
        integer, intent(in)                  :: I, J
        real(real64), intent(in)             :: cutoff, chirality
        real(real64), intent(out)            :: PE, FXIJ(3), TXI(3), TXJ(3)
        
        real(real64), dimension(3, 3)        :: U1, U2, S1, S2, E1, E2, EM
        real(real64), dimension(3, 3)        :: U1_Inv, U2_Inv, S1_sq, S2_sq
        real(real64), dimension(3, 3)        :: GE, GR, MM1, MM2, EAB
        real(real64), dimension(3)           :: MM1x, MM1y, MM1z, MM2x, MM2y, MM2z
        real(real64), dimension(3)           :: EM_X_MM1x, EM_X_MM1y, EM_X_MM1z, EM_X_MM2x, EM_X_MM2y, EM_X_MM2z
        real(real64), dimension(3)           :: KE, KR, KE_X_AE, KR_X_AR, KE_X_BE, KR_X_BR
        real(real64), dimension(3)           :: rij_hat, u1z, u2z, u2xrij, rijxu1, u1xu2
        real(real64)                         :: lambda_E, lambda_R
        real(real64)                         :: rij_sq, rij_mag
        
        real(real64), dimension(3)           :: dR_dr, de2_dr
        real(real64), dimension(3)           :: de1_du1x, de1_du1y, de1_du1z, de1_du1, de2_du1, dR_du1
        real(real64), dimension(3)           :: de1_du2x, de1_du2y, de1_du2z, de1_du2, de2_du2, dR_du2
        real(real64), dimension(3)           :: duur_dr, duur_du1z, duur_du2z, duu_du1z, duu_du2z
        real(real64), dimension(3)           :: FIJ, TORQ1, TORQ2, termu1_chiral, termu2_chiral
        real(real64), dimension(3)           :: FIJ_chiral, TORQ1_chiral, TORQ2_chiral, tu1_chiral, tu2_chiral
        
        real(real64)                         :: muinv, uur, uu
        real(real64)                         :: eps1, eps2, phi, phi_inv_3by2, sigma, sigma_0, sigma_0_inv
        real(real64)                         :: e10, e20, sigma_t1, sigma_t2, sigma_t
        real(real64)                         :: sr_ecp, sr2_ecp, sr6_ecp, sr12_ecp
        real(real64)                         :: pot, pot_chiral, efac, rfac, srx_fac, srd_fac
        real(real64)                         :: switch, xx, sx, sdx, factor
        type(ecp_pair)                       :: pair

        pot = 0.0; FIJ = 0.0; TORQ1 = 0.0; TORQ2 = 0.0; sigma_0 = 1.0;
        pot_chiral = 0.0; FIJ_chiral = 0.0; TORQ1_chiral = 0.0; TORQ2_chiral = 0.0 
        
        pair % RIJ(1) = RX(I) - RX(J)
        pair % RIJ(2) = RY(I) - RY(J)
        pair % RIJ(3) = RZ(I) - RZ(J)

        if (periodic) then
            pair % RIJ = pair % RIJ - ANINT( pair % RIJ / box ) * box
        endif

        rij_sq = pair % RIJ(1)*pair % RIJ(1) + pair % RIJ(2)*pair % RIJ(2) + pair % RIJ(3)*pair % RIJ(3)
        rij_mag = SQRT( rij_sq )
        rij_hat = pair % RIJ / rij_mag

        U1(1, 1) = QW(I)*QW(I) + QX(I)*QX(I) - QY(I)*QY(I) - QZ(I)*QZ(I)
        U1(1, 2) = 2.0 * ( QX(I) * QY(I) + QW(I) * QZ(I) )
        U1(1, 3) = 2.0 * ( QX(I) * QZ(I) - QW(I) * QY(I) )
        U1(2, 1) = 2.0 * ( QX(I) * QY(I) - QW(I) * QZ(I) )
        U1(2, 2) = QW(I)*QW(I) - QX(I)*QX(I) + QY(I)*QY(I) - QZ(I)*QZ(I)
        U1(2, 3) = 2.0 * ( QY(I) * QZ(I) + QW(I) * QX(I) )
        U1(3, 1) = 2.0 * ( QX(I) * QZ(I) + QW(I) * QY(I) )
        U1(3, 2) = 2.0 * ( QY(I) * QZ(I) - QW(I) * QX(I) )
        U1(3, 3) = QW(I)*QW(I) - QX(I)*QX(I) - QY(I)*QY(I) + QZ(I)*QZ(I)

        U2(1, 1) = QW(J)*QW(J) + QX(J)*QX(J) - QY(J)*QY(J) - QZ(J)*QZ(J)
        U2(1, 2) = 2.0 * ( QX(J) * QY(J) + QW(J) * QZ(J) )
        U2(1, 3) = 2.0 * ( QX(J) * QZ(J) - QW(J) * QY(J) )
        U2(2, 1) = 2.0 * ( QX(J) * QY(J) - QW(J) * QZ(J) )
        U2(2, 2) = QW(J)*QW(J) - QX(J)*QX(J) + QY(J)*QY(J) - QZ(J)*QZ(J)
        U2(2, 3) = 2.0 * ( QY(J) * QZ(J) + QW(J) * QX(J) )
        U2(3, 1) = 2.0 * ( QX(J) * QZ(J) + QW(J) * QY(J) )
        U2(3, 2) = 2.0 * ( QY(J) * QZ(J) - QW(J) * QX(J) )
        U2(3, 3) = QW(J)*QW(J) - QX(J)*QX(J) - QY(J)*QY(J) + QZ(J)*QZ(J)

        e10 = max( eps_x(I), eps_y(I), eps_z(I) )
        e20 = max( eps_x(J), eps_y(J), eps_z(J) )
        muinv = 1.0 / meu
        E1 = 0.0_8; E2 = 0.0_8; S1 = 0.0_8; S2 = 0.0_8

        E1(1, 1) = ( (e10 / eps_x(I)) ** muinv ) / 4.0
        E1(2, 2) = ( (e10 / eps_y(I)) ** muinv ) / 4.0
        E1(3, 3) = ( (e10 / eps_z(I)) ** muinv ) / 4.0

        E2(1, 1) = ( (e20 / eps_x(J)) ** muinv ) / 4.0
        E2(2, 2) = ( (e20 / eps_y(J)) ** muinv ) / 4.0
        E2(3, 3) = ( (e20 / eps_z(J)) ** muinv ) / 4.0

        S1(1, 1) = ( sig_x(I) / 2.0 )
        S1(2, 2) = ( sig_y(I) / 2.0 )
        S1(3, 3) = ( sig_z(I) / 2.0 )

        S2(1, 1) = ( sig_x(J) / 2.0 )
        S2(2, 2) = ( sig_y(J) / 2.0 )
        S2(3, 3) = ( sig_z(J) / 2.0 )

        U1_Inv = transpose(U1)
        U2_Inv = transpose(U2)
        S1_sq = S1*S1
        S2_sq = S2*S2
        
        MM1 = U1_Inv .x. S1
        MM2 = U2_Inv .x. S2
        MM1x = MM1(:, 1); MM1y = MM1(:, 2); MM1z = MM1(:, 3)
        MM2x = MM2(:, 1); MM2y = MM2(:, 2); MM2z = MM2(:, 3)
        pair % AR = matmul_atba(U1, S1_sq)
        pair % BR = matmul_atba(U2, S2_sq)
        pair % AE = matmul_atba(U1, E1)
        pair % BE = matmul_atba(U2, E2)
        
        sigma_t1 = (sig_x(I) * sig_y(I) + sig_z(I)*sig_z(I)) * sqrt( 2.0 * sig_x(I) * sig_y(I) ) / 8.0
        sigma_t2 = (sig_x(J) * sig_y(J) + sig_z(J)*sig_z(J)) * sqrt( 2.0 * sig_x(J) * sig_y(J) ) / 8.0
        sigma_t = sqrt( sigma_t1 * sigma_t2 )
        
        EAB = pair % AR + pair % BR
        EM = inverse( EAB )
        eps1 = sigma_t * sqrt( determinant( EM ) )
        
        lambda_R = brent( optim_dist, pair )
        lambda_E = lambda_R

        GE = ( 1 - lambda_E ) * pair % AE + lambda_E * pair % BE
        KE = cholesky_solve( GE, pair % RIJ )
        eps2 = lambda_E * (1.0 - lambda_E) * ( rij_hat(1) * KE(1) + rij_hat(2) * KE(2) + rij_hat(3) * KE(3) ) / rij_mag
        
        GR    = ( 1 - lambda_R ) * pair % AR + lambda_R * pair % BR
        KR    = cholesky_solve( GR, pair % RIJ )
        phi   = lambda_R * (1.0 - lambda_R) * ( rij_hat(1) * KR(1) + rij_hat(2) * KR(2) + rij_hat(3) * KR(3)) / rij_mag
        sigma = 1.0 / sqrt( phi )
        efac  = 2.0 * lambda_E * ( 1.0 - lambda_E ) / rij_sq
        rfac  = 2.0 * lambda_R * ( 1.0 - lambda_R ) / rij_sq
        sigma_0_inv = 1.0 / sigma_0
        phi_inv_3by2 = sigma * sigma * sigma

        ! EVALUATE SMOOTHSTEP FUNCTION

        switch = 0.9 * ( sigma + cutoff - sigma_0 )
        factor = 1.0 / ( ( sigma + cutoff - sigma_0 ) - switch )
        xx = ( rij_mag - switch ) * factor
        sx = 1.0 - 6.0 * xx*xx*xx*xx*xx + 15.0 * xx*xx*xx*xx - 10.0 * xx*xx*xx
        sdx = -30.0 * xx*xx*xx*xx + 60.0 * xx*xx*xx - 30.0 * xx*xx

        ! EVALUATE ELLIPSOID CONTACT POTENTIAL

        if ( rij_mag < ( sigma + cutoff - sigma_0 ) ) then

            sr_ecp = sigma_0 / ( rij_mag - sigma + sigma_0 )
            sr2_ecp = sr_ecp * sr_ecp
            sr6_ecp = sr2_ecp * sr2_ecp * sr2_ecp
            sr12_ecp = sr6_ecp * sr6_ecp

            srx_fac = ( sr12_ecp - 2.0 * sr6_ecp )
            srd_fac = 12.0 * ( sr6_ecp - sr12_ecp ) * sr_ecp

            dR_dr  = (eps1 ** neu) * (eps2 ** meu) * sigma_0_inv * ( rij_hat + 0.5 * rfac * ( KR - dot_product( rij_hat, KR ) * rij_hat ) * phi_inv_3by2 )
            de2_dr = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * efac * ( KE - dot_product( rij_hat, KE ) * rij_hat )

            EM_X_MM1x = EM .x. MM1x; EM_X_MM1y = EM .x. MM1y; EM_X_MM1z = EM .x. MM1z
            de1_du1x  = (eps2 ** meu) * (neu * eps1 ** neu) * cross( MM1x, EM_X_MM1x )
            de1_du1y  = (eps2 ** meu) * (neu * eps1 ** neu) * cross( MM1y, EM_X_MM1y )
            de1_du1z  = (eps2 ** meu) * (neu * eps1 ** neu) * cross( MM1z, EM_X_MM1z )

            KE_X_AE = KE .x. pair % AE; KR_X_AR = KR .x. pair % AR
            de1_du1 = (de1_du1x + de1_du1y + de1_du1z)
            de2_du1 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * (1 - lambda_E) * efac * cross( KE_X_AE, KE )
            dR_du1  = (eps1 ** neu) * (eps2 ** meu) * (0.5 * sigma_0_inv) * (1 - lambda_R) * rfac * cross( KR_X_AR, KR ) * phi_inv_3by2

            EM_X_MM2x = EM .x. MM2x; EM_X_MM2y = EM .x. MM2y; EM_X_MM2z = EM .x. MM2z
            de1_du2x  = (eps2 ** meu) * (neu * eps1 ** neu) * cross( MM2x, EM_X_MM2x )
            de1_du2y  = (eps2 ** meu) * (neu * eps1 ** neu) * cross( MM2y, EM_X_MM2y )
            de1_du2z  = (eps2 ** meu) * (neu * eps1 ** neu) * cross( MM2z, EM_X_MM2z )

            KE_X_BE = KE .x. pair % BE; KR_X_BR = KR .x. pair % BR
            de1_du2 = (de1_du2x + de1_du2y + de1_du2z)
            de2_du2 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * lambda_E * efac * cross( KE_X_BE, KE )
            dR_du2  = (eps1 ** neu) * (eps2 ** meu) * (0.5 * sigma_0_inv) * lambda_R * rfac * cross( KR_X_BR, KR ) * phi_inv_3by2
            
            pot   = (eps1 ** neu) * (eps2 ** meu) * srx_fac
            FIJ   = (-1.0) * ( srx_fac *  de2_dr             + srd_fac * dR_dr )
            TORQ1 =          ( srx_fac * (de1_du1 + de2_du1) + srd_fac * dR_du1)
            TORQ2 =          ( srx_fac * (de1_du2 + de2_du2) + srd_fac * dR_du2)

            ! EVALUATE CHIRAL TERMS

            srx_fac  = sr6_ecp * sr_ecp
            srd_fac  = -7.0 * sr6_ecp * sr2_ecp
            u1z(1) = U1(3, 1); u1z(2) = U1(3, 2); u1z(3) = U1(3, 3)
            u2z(1) = U2(3, 1); u2z(2) = U2(3, 2); u2z(3) = U2(3, 3)

            u1xu2(1) = U1Z(2) * U2Z(3) - U1Z(3) * U2Z(2)
            u1xu2(2) = U1Z(3) * U2Z(1) - U1Z(1) * U2Z(3)
            u1xu2(3) = U1Z(1) * U2Z(2) - U1Z(2) * U2Z(1) 

            u2xrij(1) = U2Z(2) * rij_hat(3) - U2Z(3) * rij_hat(2)
            u2xrij(2) = U2Z(3) * rij_hat(1) - U2Z(1) * rij_hat(3)
            u2xrij(3) = U2Z(1) * rij_hat(2) - U2Z(2) * rij_hat(1)

            rijxu1(1) = rij_hat(2) * U1Z(3) - rij_hat(3) * U1Z(2)
            rijxu1(2) = rij_hat(3) * U1Z(1) - rij_hat(1) * U1Z(3)
            rijxu1(3) = rij_hat(1) * U1Z(2) - rij_hat(2) * U1Z(1)

            uur      = u1xu2(1) * rij_hat(1) + u1xu2(2) * rij_hat(2) + u1xu2(3) * rij_hat(3)
            uu       = u1z(1) * u2z(1) + u1z(2) * u2z(2) + u1z(3) * u2z(3)

            dR_dr   = dR_dr  * srd_fac * uur * uu
            de2_dr  = de2_dr * srx_fac * uur * uu
            duur_dr = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uu * ( u1xu2 - uur * rij_hat ) / rij_mag

            de1_du1   = de1_du1 * srx_fac * uur * uu
            de2_du1   = de2_du1 * srx_fac * uur * uu
            dR_du1    = dR_du1  * srd_fac * uur * uu
            duur_du1z = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uu * u2xrij
            duu_du1z  = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uur * u2z

            de1_du2   = de1_du2 * srx_fac * uur * uu
            de2_du2   = de2_du2 * srx_fac * uur * uu
            dR_du2    = dR_du2  * srd_fac * uur * uu
            duur_du2z = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uu * rijxu1
            duu_du2z  = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uur * u1z

            termu1_chiral = ( duur_du1z + duu_du1z )
            termu2_chiral = ( duur_du2z + duu_du2z )

            tu1_chiral = cross( u1z, termu1_chiral )
            tu2_chiral = cross( u2z, termu2_chiral )

            pot_chiral   = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uur * uu
            FIJ_chiral   = (-1.0) * (de2_dr  + dR_dr  + duur_dr)
            TORQ1_chiral =          (de1_du1 + de2_du1 + dR_du1 - tu1_chiral)
            TORQ2_chiral =          (de1_du2 + de2_du2 + dR_du2 - tu2_chiral)

            pot   = pot   + chirality * pot_chiral
            FIJ   = FIJ   + chirality * FIJ_chiral
            TORQ1 = TORQ1 + chirality * TORQ1_chiral
            TORQ2 = TORQ2 + chirality * TORQ2_chiral

            if (( rij_mag >= switch ) .and. ( rij_mag <= (sigma + cutoff - sigma_0) )) then
                FIJ   = sx * FIJ   - sdx * pot * factor
                TORQ1 = sx * TORQ1 - sdx * pot * factor
                TORQ2 = sx * TORQ2 - sdx * pot * factor
                pot   = sx * pot
            endif

        endif
        
        PE = pot
        FXIJ = FIJ
        TXI = TORQ1
        TXJ = TORQ2
    end subroutine ecpc_calculate_forces

    pure subroutine ecpwc_calculate_forces( I, J, cutoff, PE, FXIJ, TXI, TXJ, wf, chirality )
        implicit none
        integer, intent(in)                  :: I, J
        real(real64), intent(in)             :: cutoff, chirality, wf
        real(real64), intent(out)            :: PE, FXIJ(3), TXI(3), TXJ(3)
        
        real(real64), dimension(3, 3)        :: U1, U2, S1, S2, E1, E2, EM
        real(real64), dimension(3, 3)        :: U1_Inv, U2_Inv, S1_sq, S2_sq
        real(real64), dimension(3, 3)        :: GE, GR, MM1, MM2, EAB
        real(real64), dimension(3)           :: MM1x, MM1y, MM1z, MM2x, MM2y, MM2z
        real(real64), dimension(3)           :: EM_X_MM1x, EM_X_MM1y, EM_X_MM1z, EM_X_MM2x, EM_X_MM2y, EM_X_MM2z
        real(real64), dimension(3)           :: KE, KR, KE_X_AE, KR_X_AR, KE_X_BE, KR_X_BR
        real(real64), dimension(3)           :: rij_hat, u1z, u2z, u2xrij, rijxu1, u1xu2
        real(real64)                         :: lambda_E, lambda_R
        real(real64)                         :: rij_sq, rij_mag
        
        real(real64), dimension(3)           :: dR_dr, de2_dr
        real(real64), dimension(3)           :: de1_du1x, de1_du1y, de1_du1z, de1_du1, de2_du1, dR_du1
        real(real64), dimension(3)           :: de1_du2x, de1_du2y, de1_du2z, de1_du2, de2_du2, dR_du2
        real(real64), dimension(3)           :: duur_dr, duur_du1z, duur_du2z, duu_du1z, duu_du2z
        real(real64), dimension(3)           :: FIJ, TORQ1, TORQ2, termu1_chiral, termu2_chiral
        real(real64), dimension(3)           :: FIJ_chiral, TORQ1_chiral, TORQ2_chiral, tu1_chiral, tu2_chiral
        
        real(real64)                         :: muinv, uur, uu
        real(real64)                         :: eps1, eps2, phi, phi_inv_3by2, sigma, sigma_0, sigma_0_inv
        real(real64)                         :: e10, e20, sigma_t1, sigma_t2, sigma_t
        real(real64)                         :: sr_ecp, sr2_ecp, sr6_ecp, sr12_ecp
        real(real64)                         :: pot, pot_chiral, efac, rfac, srx_fac, srd_fac
        real(real64)                         :: switch, xx, sx, sdx, factor
        type(ecp_pair)                       :: pair

        pot = 0.0; FIJ = 0.0; TORQ1 = 0.0; TORQ2 = 0.0; sigma_0 = 1.0;
        pot_chiral = 0.0; FIJ_chiral = 0.0; TORQ1_chiral = 0.0; TORQ2_chiral = 0.0 
        
        pair % RIJ(1) = RX(I) - RX(J)
        pair % RIJ(2) = RY(I) - RY(J)
        pair % RIJ(3) = RZ(I) - RZ(J)

        if (periodic) then
            pair % RIJ = pair % RIJ - ANINT( pair % RIJ / box ) * box
        endif

        rij_sq = pair % RIJ(1)*pair % RIJ(1) + pair % RIJ(2)*pair % RIJ(2) + pair % RIJ(3)*pair % RIJ(3)
        rij_mag = SQRT( rij_sq )
        rij_hat = pair % RIJ / rij_mag

        U1(1, 1) = QW(I)*QW(I) + QX(I)*QX(I) - QY(I)*QY(I) - QZ(I)*QZ(I)
        U1(1, 2) = 2.0 * ( QX(I) * QY(I) + QW(I) * QZ(I) )
        U1(1, 3) = 2.0 * ( QX(I) * QZ(I) - QW(I) * QY(I) )
        U1(2, 1) = 2.0 * ( QX(I) * QY(I) - QW(I) * QZ(I) )
        U1(2, 2) = QW(I)*QW(I) - QX(I)*QX(I) + QY(I)*QY(I) - QZ(I)*QZ(I)
        U1(2, 3) = 2.0 * ( QY(I) * QZ(I) + QW(I) * QX(I) )
        U1(3, 1) = 2.0 * ( QX(I) * QZ(I) + QW(I) * QY(I) )
        U1(3, 2) = 2.0 * ( QY(I) * QZ(I) - QW(I) * QX(I) )
        U1(3, 3) = QW(I)*QW(I) - QX(I)*QX(I) - QY(I)*QY(I) + QZ(I)*QZ(I)

        U2(1, 1) = QW(J)*QW(J) + QX(J)*QX(J) - QY(J)*QY(J) - QZ(J)*QZ(J)
        U2(1, 2) = 2.0 * ( QX(J) * QY(J) + QW(J) * QZ(J) )
        U2(1, 3) = 2.0 * ( QX(J) * QZ(J) - QW(J) * QY(J) )
        U2(2, 1) = 2.0 * ( QX(J) * QY(J) - QW(J) * QZ(J) )
        U2(2, 2) = QW(J)*QW(J) - QX(J)*QX(J) + QY(J)*QY(J) - QZ(J)*QZ(J)
        U2(2, 3) = 2.0 * ( QY(J) * QZ(J) + QW(J) * QX(J) )
        U2(3, 1) = 2.0 * ( QX(J) * QZ(J) + QW(J) * QY(J) )
        U2(3, 2) = 2.0 * ( QY(J) * QZ(J) - QW(J) * QX(J) )
        U2(3, 3) = QW(J)*QW(J) - QX(J)*QX(J) - QY(J)*QY(J) + QZ(J)*QZ(J)

        e10 = max( eps_x(I), eps_y(I), eps_z(I) )
        e20 = max( eps_x(J), eps_y(J), eps_z(J) )
        muinv = 1.0 / meu
        E1 = 0.0_8; E2 = 0.0_8; S1 = 0.0_8; S2 = 0.0_8

        E1(1, 1) = ( (e10 / eps_x(I)) ** muinv ) / 4.0
        E1(2, 2) = ( (e10 / eps_y(I)) ** muinv ) / 4.0
        E1(3, 3) = ( (e10 / eps_z(I)) ** muinv ) / 4.0

        E2(1, 1) = ( (e20 / eps_x(J)) ** muinv ) / 4.0
        E2(2, 2) = ( (e20 / eps_y(J)) ** muinv ) / 4.0
        E2(3, 3) = ( (e20 / eps_z(J)) ** muinv ) / 4.0

        S1(1, 1) = ( sig_x(I) / 2.0 )
        S1(2, 2) = ( sig_y(I) / 2.0 )
        S1(3, 3) = ( sig_z(I) / 2.0 )

        S2(1, 1) = ( sig_x(J) / 2.0 )
        S2(2, 2) = ( sig_y(J) / 2.0 )
        S2(3, 3) = ( sig_z(J) / 2.0 )

        U1_Inv = transpose(U1)
        U2_Inv = transpose(U2)
        S1_sq = S1*S1
        S2_sq = S2*S2
        
        MM1 = U1_Inv .x. S1
        MM2 = U2_Inv .x. S2
        MM1x = MM1(:, 1); MM1y = MM1(:, 2); MM1z = MM1(:, 3)
        MM2x = MM2(:, 1); MM2y = MM2(:, 2); MM2z = MM2(:, 3)
        pair % AR = matmul_atba(U1, S1_sq)
        pair % BR = matmul_atba(U2, S2_sq)
        pair % AE = matmul_atba(U1, E1)
        pair % BE = matmul_atba(U2, E2)
        
        sigma_t1 = (sig_x(I) * sig_y(I) + sig_z(I)*sig_z(I)) * sqrt( 2.0 * sig_x(I) * sig_y(I) ) / 8.0
        sigma_t2 = (sig_x(J) * sig_y(J) + sig_z(J)*sig_z(J)) * sqrt( 2.0 * sig_x(J) * sig_y(J) ) / 8.0
        sigma_t = sqrt( sigma_t1 * sigma_t2 )
        
        EAB = pair % AR + pair % BR
        EM = inverse( EAB )
        eps1 = sigma_t * sqrt( determinant( EM ) )
        
        lambda_R = brent( optim_dist, pair )
        lambda_E = lambda_R

        GE = ( 1 - lambda_E ) * pair % AE + lambda_E * pair % BE
        KE = cholesky_solve( GE, pair % RIJ )
        eps2 = lambda_E * (1.0 - lambda_E) * ( rij_hat(1) * KE(1) + rij_hat(2) * KE(2) + rij_hat(3) * KE(3) ) / rij_mag
        
        GR    = ( 1 - lambda_R ) * pair % AR + lambda_R * pair % BR
        KR    = cholesky_solve( GR, pair % RIJ )
        phi   = lambda_R * (1.0 - lambda_R) * ( rij_hat(1) * KR(1) + rij_hat(2) * KR(2) + rij_hat(3) * KR(3)) / rij_mag
        sigma = 1.0 / sqrt( phi )
        efac  = 2.0 * lambda_E * ( 1.0 - lambda_E ) / rij_sq
        rfac  = 2.0 * lambda_R * ( 1.0 - lambda_R ) / rij_sq
        sigma_0_inv = 1.0 / sigma_0
        phi_inv_3by2 = sigma * sigma * sigma

        ! EVALUATE SMOOTHSTEP FUNCTION

        switch = 0.9 * ( sigma + cutoff + wf - sigma_0 )
        factor = 1.0 / ( ( sigma + cutoff + wf - sigma_0 ) - switch )
        xx = ( rij_mag - switch ) * factor
        sx = 1.0 - 6.0 * xx*xx*xx*xx*xx + 15.0 * xx*xx*xx*xx - 10.0 * xx*xx*xx
        sdx = -30.0 * xx*xx*xx*xx + 60.0 * xx*xx*xx - 30.0 * xx*xx

        ! EVALUATE ELLIPSOID CONTACT POTENTIAL

        if ( rij_mag < ( sigma + cutoff + wf - sigma_0 ) ) then

            if ( rij_mag < sigma ) then
                sr_ecp = sigma_0 / ( rij_mag - sigma + sigma_0 )
            else if ((rij_mag .ge. sigma) .and. (rij_mag .le. (sigma + wf))) then
                sr_ecp = 1.0
            else
                sr_ecp = sigma_0 / ( rij_mag - sigma + sigma_0 - wf )
            endif

            sr2_ecp = sr_ecp * sr_ecp
            sr6_ecp = sr2_ecp * sr2_ecp * sr2_ecp
            sr12_ecp = sr6_ecp * sr6_ecp

            srx_fac = ( sr12_ecp - 2.0 * sr6_ecp )
            srd_fac = 12.0 * ( sr6_ecp - sr12_ecp ) * sr_ecp

            dR_dr  = (eps1 ** neu) * (eps2 ** meu) * sigma_0_inv * ( rij_hat + 0.5 * rfac * ( KR - dot_product( rij_hat, KR ) * rij_hat ) * phi_inv_3by2 )
            de2_dr = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * efac * ( KE - dot_product( rij_hat, KE ) * rij_hat )

            EM_X_MM1x = EM .x. MM1x; EM_X_MM1y = EM .x. MM1y; EM_X_MM1z = EM .x. MM1z
            de1_du1x  = (eps2 ** meu) * (neu * eps1 ** neu) * cross( MM1x, EM_X_MM1x )
            de1_du1y  = (eps2 ** meu) * (neu * eps1 ** neu) * cross( MM1y, EM_X_MM1y )
            de1_du1z  = (eps2 ** meu) * (neu * eps1 ** neu) * cross( MM1z, EM_X_MM1z )

            KE_X_AE = KE .x. pair % AE; KR_X_AR = KR .x. pair % AR
            de1_du1 = (de1_du1x + de1_du1y + de1_du1z)
            de2_du1 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * (1 - lambda_E) * efac * cross( KE_X_AE, KE )
            dR_du1  = (eps1 ** neu) * (eps2 ** meu) * (0.5 * sigma_0_inv) * (1 - lambda_R) * rfac * cross( KR_X_AR, KR ) * phi_inv_3by2

            EM_X_MM2x = EM .x. MM2x; EM_X_MM2y = EM .x. MM2y; EM_X_MM2z = EM .x. MM2z
            de1_du2x  = (eps2 ** meu) * (neu * eps1 ** neu) * cross( MM2x, EM_X_MM2x )
            de1_du2y  = (eps2 ** meu) * (neu * eps1 ** neu) * cross( MM2y, EM_X_MM2y )
            de1_du2z  = (eps2 ** meu) * (neu * eps1 ** neu) * cross( MM2z, EM_X_MM2z )

            KE_X_BE = KE .x. pair % BE; KR_X_BR = KR .x. pair % BR
            de1_du2 = (de1_du2x + de1_du2y + de1_du2z)
            de2_du2 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * lambda_E * efac * cross( KE_X_BE, KE )
            dR_du2  = (eps1 ** neu) * (eps2 ** meu) * (0.5 * sigma_0_inv) * lambda_R * rfac * cross( KR_X_BR, KR ) * phi_inv_3by2
            
            pot   = (eps1 ** neu) * (eps2 ** meu) * srx_fac
            FIJ   = (-1.0) * ( srx_fac *  de2_dr             + srd_fac * dR_dr )
            TORQ1 =          ( srx_fac * (de1_du1 + de2_du1) + srd_fac * dR_du1)
            TORQ2 =          ( srx_fac * (de1_du2 + de2_du2) + srd_fac * dR_du2)

            ! EVALUATE CHIRAL TERMS

            srx_fac = ( sr12_ecp - (12.0 / 7.0) * sr6_ecp * sr_ecp)
            srd_fac = 12.0 * ( sr6_ecp * sr2_ecp - sr12_ecp * sr_ecp )

            u1z(1) = U1(3, 1); u1z(2) = U1(3, 2); u1z(3) = U1(3, 3)
            u2z(1) = U2(3, 1); u2z(2) = U2(3, 2); u2z(3) = U2(3, 3)

            u1xu2(1) = U1Z(2) * U2Z(3) - U1Z(3) * U2Z(2)
            u1xu2(2) = U1Z(3) * U2Z(1) - U1Z(1) * U2Z(3)
            u1xu2(3) = U1Z(1) * U2Z(2) - U1Z(2) * U2Z(1) 

            u2xrij(1) = U2Z(2) * rij_hat(3) - U2Z(3) * rij_hat(2)
            u2xrij(2) = U2Z(3) * rij_hat(1) - U2Z(1) * rij_hat(3)
            u2xrij(3) = U2Z(1) * rij_hat(2) - U2Z(2) * rij_hat(1)

            rijxu1(1) = rij_hat(2) * U1Z(3) - rij_hat(3) * U1Z(2)
            rijxu1(2) = rij_hat(3) * U1Z(1) - rij_hat(1) * U1Z(3)
            rijxu1(3) = rij_hat(1) * U1Z(2) - rij_hat(2) * U1Z(1)

            uur      = u1xu2(1) * rij_hat(1) + u1xu2(2) * rij_hat(2) + u1xu2(3) * rij_hat(3)
            uu       = u1z(1) * u2z(1) + u1z(2) * u2z(2) + u1z(3) * u2z(3)

            dR_dr   = dR_dr  * srd_fac * uur * uu
            de2_dr  = de2_dr * srx_fac * uur * uu
            duur_dr = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uu * ( u1xu2 - uur * rij_hat ) / rij_mag

            de1_du1   = de1_du1 * srx_fac * uur * uu
            de2_du1   = de2_du1 * srx_fac * uur * uu
            dR_du1    = dR_du1  * srd_fac * uur * uu
            duur_du1z = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uu * u2xrij
            duu_du1z  = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uur * u2z

            de1_du2   = de1_du2 * srx_fac * uur * uu
            de2_du2   = de2_du2 * srx_fac * uur * uu
            dR_du2    = dR_du2  * srd_fac * uur * uu
            duur_du2z = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uu * rijxu1
            duu_du2z  = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uur * u1z

            termu1_chiral = ( duur_du1z + duu_du1z )
            termu2_chiral = ( duur_du2z + duu_du2z )

            tu1_chiral = cross( u1z, termu1_chiral )
            tu2_chiral = cross( u2z, termu2_chiral )

            pot_chiral   = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uur * uu
            FIJ_chiral   = (-1.0) * (de2_dr  + dR_dr  + duur_dr)
            TORQ1_chiral =          (de1_du1 + de2_du1 + dR_du1 - tu1_chiral)
            TORQ2_chiral =          (de1_du2 + de2_du2 + dR_du2 - tu2_chiral)

            pot   = pot   + chirality * pot_chiral
            FIJ   = FIJ   + chirality * FIJ_chiral
            TORQ1 = TORQ1 + chirality * TORQ1_chiral
            TORQ2 = TORQ2 + chirality * TORQ2_chiral

            if (( rij_mag >= switch ) .and. ( rij_mag <= (sigma + cutoff + wf - sigma_0) )) then
                FIJ   = sx * FIJ   - sdx * pot * factor
                TORQ1 = sx * TORQ1 - sdx * pot * factor
                TORQ2 = sx * TORQ2 - sdx * pot * factor
                pot   = sx * pot
            endif

        endif
        
        PE = pot
        FXIJ = FIJ
        TXI = TORQ1
        TXJ = TORQ2
    end subroutine ecpwc_calculate_forces

end module ecp

module bonded
    use xmath, only: PI
    use system
    implicit none

    contains
    pure subroutine bend_calculate_forces( I, J, K, PE, FXI, FXJ, FXK, coeff )
        implicit none
        integer, intent(in)         :: I, J, K
        real(real64), intent(in)    :: coeff
        real(real64), intent(out)   :: PE, FXI(3), FXJ(3), FXK(3)
        real(real64), dimension(3)  :: VI, VJ, VI_PBC, VJ_PBC
        real(real64)                :: CC11, CC12, CC22
        real(real64)                :: prefac, fac, fac1, fac2

        PE = 0.0; FXI = 0.0; FXJ = 0.0; FXK = 0.0

        VI(1) = RX(J) - RX(I)
        VI(2) = RY(J) - RY(I)
        VI(3) = RZ(J) - RZ(I)

        VJ(1) = RX(K) - RX(J)
        VJ(2) = RY(K) - RY(J)
        VJ(3) = RZ(K) - RZ(J)

        VI_PBC = VI - ANINT( VI / box ) * box
        VJ_PBC = VJ - ANINT( VJ / box ) * box

        VI = MERGE( VI_PBC, VI, periodic )
        VJ = MERGE( VJ_PBC, VJ, periodic )

        CC11 = VI(1) * VI(1) + VI(2) * VI(2) + VI(3) * VI(3)
        CC12 = VI(1) * VJ(1) + VI(2) * VJ(2) + VI(3) * VJ(3)
        CC22 = VJ(1) * VJ(1) + VJ(2) * VJ(2) + VJ(3) * VJ(3)

        prefac = 1.0 / SQRT( CC11 * CC22 )
        fac    = CC12
        fac1   = fac / CC22
        fac2   = fac / CC11

        PE  = coeff * (1.0 - prefac * fac)

        FXK = -coeff * prefac * ( fac1 * VJ - VI )
        FXJ =  coeff * prefac * ( fac1 * VJ - fac2 * VI + VJ - VI )
        FXI =  coeff * prefac * ( fac2 * VI - VJ )
        
    end subroutine bend_calculate_forces

    pure subroutine dih_calculate_forces( I, J, K, L, PE, FXI, FXJ, FXK, FXL, coeff )
        implicit none
        integer, intent(in)         :: I, J, K, L
        real(real64), intent(in)    :: coeff
        real(real64), intent(out)   :: PE, FXI(3), FXJ(3), FXK(3), FXL(3)
        real(real64), dimension(3)  :: VI, VJ, VK, VI_PBC, VJ_PBC, VK_PBC
        real(real64)                :: CC11, CC12, CC13, CC22, CC23, CC33
        real(real64)                :: DD12, DD23
        real(real64)                :: prefac, fac, fac1, fac2

        VI(1) = RX(J) - RX(I)
        VI(2) = RY(J) - RY(I)
        VI(3) = RZ(J) - RZ(I)

        VJ(1) = RX(K) - RX(J)
        VJ(2) = RY(K) - RY(J)
        VJ(3) = RZ(K) - RZ(J)

        VK(1) = RX(L) - RX(K)
        VK(2) = RY(L) - RY(K)
        VK(3) = RZ(L) - RZ(K)

        VI_PBC = VI - ANINT( VI / box ) * box
        VJ_PBC = VJ - ANINT( VJ / box ) * box
        VK_PBC = VK - ANINT( VK / box ) * box

        VI = MERGE( VI_PBC, VI, periodic )
        VJ = MERGE( VJ_PBC, VJ, periodic )
        VK = MERGE( VK_PBC, VK, periodic )

        CC11 = VI(1) * VI(1) + VI(2) * VI(2) + VI(3) * VI(3)
        CC12 = VI(1) * VJ(1) + VI(2) * VJ(2) + VI(3) * VJ(3)
        CC13 = VI(1) * VK(1) + VI(2) * VK(2) + VI(3) * VK(3)
        CC22 = VJ(1) * VJ(1) + VJ(2) * VJ(2) + VJ(3) * VJ(3)
        CC23 = VJ(1) * VK(1) + VJ(2) * VK(2) + VJ(3) * VK(3)
        CC33 = VK(1) * VK(1) + VK(2) * VK(2) + VK(3) * VK(3)

        DD12 = CC11 * CC22 - CC12 * CC12
        DD23 = CC22 * CC33 - CC23 * CC23

        prefac  = 1.0 / SQRT( DD23 * DD12 )
        fac     = CC23 * CC12 - CC13 * CC22
        fac1    = fac / DD23
        fac2    = fac / DD12

        PE  = coeff * (1.0 + prefac * fac)

        FXL = -coeff * prefac * ( CC12 * VJ - CC22 * VI - fac1 * ( CC22 * VK- CC23 * VJ ) )

        FXK = -coeff * prefac * ( CC12 * VK - CC12 * VJ + CC23 * VI + CC22 * VI - 2.0 * CC13 * VJ &
                        - fac2  * ( CC11 * VJ - CC12 * VI ) - fac1 * ( CC33 * VJ - CC22 * VK - CC23 * VK + CC23 * VJ ) )

        FXJ = -coeff * prefac * ( -CC12 * VK + CC23 * VJ - CC23 * VI - CC22 * VK + 2.0 * CC13 * VJ &
                        - fac2  * (  CC22 * VI - CC11 * VJ - CC12 * VJ + CC12 * VI ) - fac1 * ( -CC33 * VJ + CC23 * VK ) )

        FXI = -coeff * prefac * ( -CC23 * VJ + CC22 * VK - fac2 * ( -CC22 * VI + CC12 * VJ ) )

    end subroutine dih_calculate_forces

    pure subroutine cosine_bend_calculate_forces( I, J, K, PE, FXI, FXJ, FXK, coeff, angle )
        implicit none
        integer, intent(in)         :: I, J, K
        real(real64), intent(in)    :: coeff, angle
        real(real64), intent(out)   :: PE, FXI(3), FXJ(3), FXK(3)
        real(real64), dimension(3)  :: VI, VJ, VI_PBC, VJ_PBC
        real(real64)                :: CC11, CC12, CC22
        real(real64)                :: prefac, fac, fac1, fac2
        real(real64)                :: cos_theta, cos_theta0, theta0, strength

        PE = 0.0; FXI = 0.0; FXJ = 0.0; FXK = 0.0

        VI(1) = RX(J) - RX(I)
        VI(2) = RY(J) - RY(I)
        VI(3) = RZ(J) - RZ(I)

        VJ(1) = RX(K) - RX(J)
        VJ(2) = RY(K) - RY(J)
        VJ(3) = RZ(K) - RZ(J)

        VI_PBC = VI - ANINT( VI / box ) * box
        VJ_PBC = VJ - ANINT( VJ / box ) * box

        VI = MERGE( VI_PBC, VI, periodic )
        VJ = MERGE( VJ_PBC, VJ, periodic )

        CC11 = VI(1) * VI(1) + VI(2) * VI(2) + VI(3) * VI(3)
        CC12 = VI(1) * VJ(1) + VI(2) * VJ(2) + VI(3) * VJ(3)
        CC22 = VJ(1) * VJ(1) + VJ(2) * VJ(2) + VJ(3) * VJ(3)

        prefac = 1.0 / SQRT( CC11 * CC22 )
        fac    = CC12
        fac1   = fac / CC22
        fac2   = fac / CC11

        theta0     = PI - (PI / 180.0) * angle
        cos_theta0 = cos(theta0)
        cos_theta  = prefac * fac
        strength   = coeff * (cos_theta - cos_theta0)
        PE         = (coeff / 2.0) * (cos_theta - cos_theta0)*(cos_theta - cos_theta0)

        FXK = -strength * prefac * ( fac1 * VJ - VI )
        FXJ =  strength * prefac * ( fac1 * VJ - fac2 * VI + VJ - VI )
        FXI =  strength * prefac * ( fac2 * VI - VJ )
        
    end subroutine cosine_bend_calculate_forces

    pure subroutine cosine_dih_calculate_forces( I, J, K, L, PE, FXI, FXJ, FXK, FXL, coeff, angle )
        implicit none
        integer, intent(in)         :: I, J, K, L
        real(real64), intent(in)    :: coeff, angle
        real(real64), intent(out)   :: PE, FXI(3), FXJ(3), FXK(3), FXL(3)
        real(real64), dimension(3)  :: VI, VJ, VK, vixvj, VI_PBC, VJ_PBC, VK_PBC
        real(real64)                :: cos_phi, sin_phi, cos_phi0, phi0, sgn_phi
        real(real64)                :: CC11, CC12, CC13, CC22, CC23, CC33
        real(real64)                :: DD12, DD23
        real(real64)                :: prefac, fac, fac1, fac2, fac_sin, strength

        VI(1) = RX(J) - RX(I)
        VI(2) = RY(J) - RY(I)
        VI(3) = RZ(J) - RZ(I)

        VJ(1) = RX(K) - RX(J)
        VJ(2) = RY(K) - RY(J)
        VJ(3) = RZ(K) - RZ(J)

        VK(1) = RX(L) - RX(K)
        VK(2) = RY(L) - RY(K)
        VK(3) = RZ(L) - RZ(K)

        VI_PBC = VI - ANINT( VI / box ) * box
        VJ_PBC = VJ - ANINT( VJ / box ) * box
        VK_PBC = VK - ANINT( VK / box ) * box

        VI = MERGE( VI_PBC, VI, periodic )
        VJ = MERGE( VJ_PBC, VJ, periodic )
        VK = MERGE( VK_PBC, VK, periodic )

        CC11 = VI(1) * VI(1) + VI(2) * VI(2) + VI(3) * VI(3)
        CC12 = VI(1) * VJ(1) + VI(2) * VJ(2) + VI(3) * VJ(3)
        CC13 = VI(1) * VK(1) + VI(2) * VK(2) + VI(3) * VK(3)
        CC22 = VJ(1) * VJ(1) + VJ(2) * VJ(2) + VJ(3) * VJ(3)
        CC23 = VJ(1) * VK(1) + VJ(2) * VK(2) + VJ(3) * VK(3)
        CC33 = VK(1) * VK(1) + VK(2) * VK(2) + VK(3) * VK(3)

        DD12 = CC11 * CC22 - CC12 * CC12
        DD23 = CC22 * CC33 - CC23 * CC23

        prefac  = 1.0 / SQRT( DD23 * DD12 )
        fac     = CC23 * CC12 - CC13 * CC22
        fac1    = fac / DD23
        fac2    = fac / DD12

        vixvj(1) = vi(2) * vj(3) - vi(3) * vj(2)
        vixvj(2) = vi(3) * vj(1) - vi(1) * vj(3)
        vixvj(3) = vi(1) * vj(2) - vi(2) * vj(1)
        fac_sin = vixvj(1) * VK(1) + vixvj(2) * VK(2) + vixvj(3) * VK(3)

        phi0     = PI - (PI / 180.0) * angle
        cos_phi0 = cos(phi0)
        cos_phi  = -prefac * fac
        sin_phi  = prefac * fac_sin * sqrt( vj(1) * vj(1) + vj(2) * vj(2) + vj(3) * vj(3) )
        sgn_phi = atan2( sin_phi, cos_phi )

        strength = coeff * (cos(sgn_phi) - cos_phi0)
        PE       = (coeff / 2.0) * (cos(sgn_phi) - cos_phi0)*(cos(sgn_phi) - cos_phi0)

        FXL = -strength * prefac * ( CC12 * VJ - CC22 * VI - fac1 * ( CC22 * VK- CC23 * VJ ) )

        FXK = -strength * prefac * ( CC12 * VK - CC12 * VJ + CC23 * VI + CC22 * VI - 2.0 * CC13 * VJ &
                        - fac2  * ( CC11 * VJ - CC12 * VI ) - fac1 * ( CC33 * VJ - CC22 * VK - CC23 * VK + CC23 * VJ ) )

        FXJ = -strength * prefac * ( -CC12 * VK + CC23 * VJ - CC23 * VI - CC22 * VK + 2.0 * CC13 * VJ &
                        - fac2  * (  CC22 * VI - CC11 * VJ - CC12 * VJ + CC12 * VI ) - fac1 * ( -CC33 * VJ + CC23 * VK ) )

        FXI = -strength * prefac * ( -CC23 * VJ + CC22 * VK - fac2 * ( -CC22 * VI + CC12 * VJ ) )

    end subroutine cosine_dih_calculate_forces

    subroutine angle_bend_calculate_forces( I, J, K, PE, FXI, FXJ, FXK, coeff, angle )
        implicit none
        integer, intent(in)         :: I, J, K
        real(real64), intent(in)    :: coeff, angle
        real(real64), intent(out)   :: PE, FXI(3), FXJ(3), FXK(3)
        real(real64), dimension(3)  :: VI, VJ, VI_PBC, VJ_PBC
        real(real64)                :: CC11, CC12, CC22
        real(real64)                :: prefac, fac, fac1, fac2
        real(real64)                :: theta, theta_0, strength

        VI(1) = RX(J) - RX(I)
        VI(2) = RY(J) - RY(I)
        VI(3) = RZ(J) - RZ(I)

        VJ(1) = RX(K) - RX(J)
        VJ(2) = RY(K) - RY(J)
        VJ(3) = RZ(K) - RZ(J)

        VI_PBC = VI - ANINT( VI / box ) * box
        VJ_PBC = VJ - ANINT( VJ / box ) * box

        VI = MERGE( VI_PBC, VI, periodic )
        VJ = MERGE( VJ_PBC, VJ, periodic )

        CC11 = VI(1) * VI(1) + VI(2) * VI(2) + VI(3) * VI(3)
        CC12 = VI(1) * VJ(1) + VI(2) * VJ(2) + VI(3) * VJ(3)
        CC22 = VJ(1) * VJ(1) + VJ(2) * VJ(2) + VJ(3) * VJ(3)

        prefac = 1.0 / SQRT( CC11 * CC22 )
        fac    = CC12
        fac1   = fac / CC22
        fac2   = fac / CC11

        theta_0  = PI - (PI / 180.0) * angle
        theta    = acos(prefac * fac)
        strength = coeff * ((theta - theta_0)) / sin(theta)
        PE       = (coeff / 2.0) * (theta - theta_0)*(theta - theta_0)

        FXK = -strength * prefac * ( fac1 * VJ - VI )
        FXJ =  strength * prefac * ( fac1 * VJ - fac2 * VI + VJ - VI )
        FXI =  strength * prefac * ( fac2 * VI - VJ )
        
    end subroutine angle_bend_calculate_forces

    subroutine angle_dih_calculate_forces( I, J, K, L, PE, FXI, FXJ, FXK, FXL, coeff, angle )
        implicit none
        integer, intent(in)         :: I, J, K, L
        real(real64), intent(in)    :: coeff, angle
        real(real64), intent(out)   :: PE, FXI(3), FXJ(3), FXK(3), FXL(3)
        real(real64), dimension(3)  :: VI, VJ, VK, vixvj, VI_PBC, VJ_PBC, VK_PBC
        real(real64)                :: CC11, CC12, CC13, CC22, CC23, CC33
        real(real64)                :: DD12, DD23
        real(real64)                :: prefac, fac, fac1, fac2, fac_sin, strength
        real(real64)                :: phi0, cos_phi0, sin_phi0, cos_phi, sin_phi

        VI(1) = RX(J) - RX(I)
        VI(2) = RY(J) - RY(I)
        VI(3) = RZ(J) - RZ(I)

        VJ(1) = RX(K) - RX(J)
        VJ(2) = RY(K) - RY(J)
        VJ(3) = RZ(K) - RZ(J)

        VK(1) = RX(L) - RX(K)
        VK(2) = RY(L) - RY(K)
        VK(3) = RZ(L) - RZ(K)

        VI_PBC = VI - ANINT( VI / box ) * box
        VJ_PBC = VJ - ANINT( VJ / box ) * box
        VK_PBC = VK - ANINT( VK / box ) * box

        VI = MERGE( VI_PBC, VI, periodic )
        VJ = MERGE( VJ_PBC, VJ, periodic )
        VK = MERGE( VK_PBC, VK, periodic )

        CC11 = VI(1) * VI(1) + VI(2) * VI(2) + VI(3) * VI(3)
        CC12 = VI(1) * VJ(1) + VI(2) * VJ(2) + VI(3) * VJ(3)
        CC13 = VI(1) * VK(1) + VI(2) * VK(2) + VI(3) * VK(3)
        CC22 = VJ(1) * VJ(1) + VJ(2) * VJ(2) + VJ(3) * VJ(3)
        CC23 = VJ(1) * VK(1) + VJ(2) * VK(2) + VJ(3) * VK(3)
        CC33 = VK(1) * VK(1) + VK(2) * VK(2) + VK(3) * VK(3)

        DD12 = CC11 * CC22 - CC12 * CC12
        DD23 = CC22 * CC33 - CC23 * CC23

        prefac = 1.0 / SQRT( DD23 * DD12 )
        fac    = CC23 * CC12 - CC13 * CC22
        fac1   = fac / DD23
        fac2   = fac / DD12

        vixvj(1) = vi(2) * vj(3) - vi(3) * vj(2)
        vixvj(2) = vi(3) * vj(1) - vi(1) * vj(3)
        vixvj(3) = vi(1) * vj(2) - vi(2) * vj(1)
        fac_sin  = vixvj(1) * VK(1) + vixvj(2) * VK(2) + vixvj(3) * VK(3)

        phi0     = PI - (PI / 180.0) * angle
        cos_phi0 = cos(phi0)
        sin_phi0 = sin(phi0)

        cos_phi  = -prefac * fac
        sin_phi  = prefac * fac_sin * sqrt( vj(1) * vj(1) + vj(2) * vj(2) + vj(3) * vj(3) )

        strength = coeff * ( cos_phi0 - sin_phi0 * ( cos_phi / sin_phi ) )
        PE       = coeff * ( 1.0 - cos_phi * cos_phi0 - sin_phi * sin_phi0 )

        FXL = -strength * prefac * ( CC12 * VJ - CC22 * VI - fac1 * ( CC22 * VK- CC23 * VJ ) )

        FXK = -strength * prefac * ( CC12 * VK - CC12 * VJ + CC23 * VI + CC22 * VI - 2.0 * CC13 * VJ &
                        - fac2  * ( CC11 * VJ - CC12 * VI ) - fac1 * ( CC33 * VJ - CC22 * VK - CC23 * VK + CC23 * VJ ) )

        FXJ = -strength * prefac * ( -CC12 * VK + CC23 * VJ - CC23 * VI - CC22 * VK + 2.0 * CC13 * VJ &
                        - fac2  * (  CC22 * VI - CC11 * VJ - CC12 * VJ + CC12 * VI ) - fac1 * ( -CC33 * VJ + CC23 * VK ) )

        FXI = -strength * prefac * ( -CC23 * VJ + CC22 * VK - fac2 * ( -CC22 * VI + CC12 * VJ ) )

    end subroutine angle_dih_calculate_forces

end module bonded

module debye_huckel
    use xmath, only: PI
    use system
    implicit none

    contains
    pure subroutine dh_calculate_forces( I, J, cutoff, db_len, dielectric, PE, FXIJ )
        implicit none
        integer, intent(in)         :: I, J
        real(real64), intent(in)    :: cutoff, db_len, dielectric
        real(real64), intent(out)   :: PE, FXIJ(3)
        real(real64)                :: RIJ(3), RIJ_PBC(3), rij_sq, rij_mag
        logical                     :: cmask

        PE = 0.0; FXIJ = 0.0

        RIJ(1) = RX(I) - RX(J)
        RIJ(2) = RY(I) - RY(J)
        RIJ(3) = RZ(I) - RZ(J)

        RIJ_PBC = RIJ - ANINT( RIJ / box ) * box
        RIJ = MERGE( RIJ_PBC, RIJ, periodic )

        rij_sq  = RIJ(1)*RIJ(1) + RIJ(2)*RIJ(2) + RIJ(3)*RIJ(3)
        rij_mag = SQRT( rij_sq )

        cmask = (rij_mag < cutoff)

        PE    = ( charge(I) * charge(J) * exp(-rij_mag / db_len ) ) / (4.0 * PI * dielectric * rij_mag)
        FXIJ  = PE * ( 1.0 + rij_mag / db_len ) * RIJ / rij_sq

        PE = MERGE( PE, 0.0, cmask )
        FXIJ = MERGE( FXIJ, 0.0, cmask )

    end subroutine dh_calculate_forces

end module debye_huckel

module coulumbic
    use xmath, only: PI
    use system
    implicit none

    contains
    pure subroutine coulombic_calculate_forces( I, J, cutoff, dielectric, PE, FXIJ )
        implicit none
        integer, intent(in)         :: I, J
        real(real64), intent(in)    :: cutoff, dielectric
        real(real64), intent(out)   :: PE, FXIJ(3)
        real(real64)                :: RIJ(3), rij_sq, rij_mag
        real(real64)                :: switch, factor, xx, sx, sdx

        PE = 0.0; FXIJ = 0.0

        RIJ(1) = RX(I) - RX(J)
        RIJ(2) = RY(I) - RY(J)
        RIJ(3) = RZ(I) - RZ(J)
        
        if (periodic) then
            RIJ = RIJ - ANINT( RIJ / box ) * box
        endif

        rij_sq  = SUM( RIJ*RIJ )
        rij_mag = SQRT( rij_sq )

        switch = cutoff - 2.0    ! Apply switch function 2sigma_0 inside cutoff
        factor = 1.0 / ( cutoff - switch )
        xx = ( rij_mag - switch ) * factor
        sx = 1.0 - 6.0 * xx*xx*xx*xx*xx + 15.0 * xx*xx*xx*xx - 10.0 * xx*xx*xx
        sdx = -30.0 * xx*xx*xx*xx + 60.0 * xx*xx*xx - 30.0 * xx*xx

        if ( rij_mag < cutoff ) then
            PE    = charge(I) * charge(J) / (4.0 * PI * dielectric * rij_mag)
            FXIJ  = PE * RIJ / rij_sq

            if ((rij_mag .ge. switch)) then
                PE  = sx * PE
                FXIJ = sx * FXIJ - sdx * PE * factor
            endif

        endif

    end subroutine coulombic_calculate_forces

end module coulumbic

module dipole
    use xmath
    use system
    implicit none
    real(real64), dimension(:), allocatable :: DIPX, DIPY, DIPZ
    real(real64), dimension(:), allocatable :: RFLX, RFLY, RFLZ

    contains
    subroutine init_reaction_field
        implicit none

        allocate( DIPX(N), DIPY(N), DIPZ(N) )
        allocate( RFLX(N), RFLY(N), RFLZ(N) )
        DIPX = 0.0; DIPY = 0.0; DIPZ = 0.0
        RFLX = 0.0; RFLY = 0.0; RFLZ = 0.0

    end subroutine init_reaction_field

    subroutine dipole_calculate_forces( I, J, cutoff, dielectric, PE, FXIJ, TXI, TXJ )
        implicit none
        integer, intent(in)         :: I, J
        real(real64), intent(in)    :: cutoff, dielectric
        real(real64), intent(out)   :: PE, FXIJ(3), TXI(3), TXJ(3)
        real(real64), dimension(3)  :: RIJ, U1, U2
        real(real64), dimension(3)  :: term1, term2, term3
        real(real64)                :: RIJ_SQ, rij_cube, rij_pow5, RIJ_MAG
        real(real64)                :: dip_coeff, ci, cj, cij

        PE = 0.0; FXIJ = 0.0; TXI = 0.0; TXJ = 0.0
        dip_coeff = (2.0 * (dielectric - 1.0)) / (2.0 * dielectric + 1.0)
        dip_coeff = dip_coeff / ( cutoff*cutoff*cutoff )

        RIJ(1) = RX(I) - RX(J)
        RIJ(2) = RY(I) - RY(J)
        RIJ(3) = RZ(I) - RZ(J)

        if (periodic) then
            RIJ = RIJ - ANINT( RIJ / box ) * box
        endif

        RIJ_SQ  = SUM( RIJ * RIJ )
        RIJ_MAG = SQRT( RIJ_SQ )

        U1(1) = 2.0 *  ( QX(I) * QZ(I) + QW(I) * QY(I) )
        U1(2) = 2.0 *  ( QY(I) * QZ(I) - QW(I) * QX(I) )
        U1(3) = QW(I)*QW(I) - QX(I)*QX(I) - QY(I)*QY(I) + QZ(I)*QZ(I)

        U2(1) = 2.0 *  ( QX(J) * QZ(J) + QW(J) * QY(J) )
        U2(2) = 2.0 *  ( QY(J) * QZ(J) - QW(J) * QX(J) )
        U2(3) = QW(J)*QW(J) - QX(J)*QX(J) - QY(J)*QY(J) + QZ(J)*QZ(J)

        U1 = U1 * charge(I)
        U2 = U2 * charge(J)

        if ( RIJ_MAG < cutoff ) then

            ci  = U1(1) * RIJ(1) + U1(2) * RIJ(2) + U1(3) * RIJ(3)
            cj  = U2(1) * RIJ(1) + U2(2) * RIJ(2) + U2(3) * RIJ(3)
            cij = U1(1) * U2(1)  + U1(2) * U2(2)  + U1(3) * U2(3)

            term1 = (cij - (5.0 * ci * cj / rij_sq )) * RIJ + cj * U1 + ci * U2
            term2 = 3.0 * cj * ( RIJ / rij_sq ) - U2
            term3 = 3.0 * ci * ( RIJ / rij_sq ) - U1

            rij_cube = rij_sq * RIJ_mag
            rij_pow5   = rij_cube * rij_sq

            PE   = ( cij / rij_cube ) - ((3.0 * ci * cj) / rij_pow5 )
            FXIJ = ( 3.0 / rij_pow5 ) * term1
            TXI  = ( 1.0 / rij_cube ) * cross( U1, term2 )
            TXJ  = ( 1.0 / rij_cube ) * cross( U2, term3 )

            DIPX(I) = U1(1)
            DIPY(I) = U1(2)
            DIPZ(I) = U1(3)

            RFLX = RFLX + U2(1) * dip_coeff
            RFLY = RFLY + U2(2) * dip_coeff
            RFLZ = RFLZ + U2(3) * dip_coeff

        endif

    end subroutine dipole_calculate_forces

    subroutine eval_reaction_field()
        implicit none
        integer :: I

        do I = 1, N
            if ( charge(I) /= 0.0 ) then
                pot_eng = pot_eng + 0.5 * ( DIPX(I) * RFLX(I) + DIPY(I) * RFLY(I) + DIPZ(I) * RFLZ(I) )
                TX(I) = TX(I) - ( DIPY(I) * RFLZ(I) - DIPZ(I) * RFLY(I) )
                TY(I) = TY(I) - ( DIPZ(I) * RFLX(I) - DIPX(I) * RFLZ(I) )
                TZ(I) = TZ(I) - ( DIPX(I) * RFLY(I) - DIPY(I) * RFLX(I) )
            endif
        enddo

    end subroutine eval_reaction_field
end module dipole