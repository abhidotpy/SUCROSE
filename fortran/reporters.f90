module state_data_reporter
    use system
    implicit none
    character(len=256) :: file

    contains
    subroutine open_state_data_file( unit, filename )
        implicit none
        integer, intent(in)          :: unit
        character(len=*), intent(in) :: filename

        file = filename

        open(unit=unit, file=file, action='write')

        write(unit, "('Step', 2x, 'Potential_Energy', 2x, 'Kinetic_Energy', 2x, 'Total_Energy', 2x, 'Temperature', 2x, 'Pressure', 2x, 'Density')")

    end subroutine open_state_data_file

    subroutine report_state_data( unit, step )
        implicit none
        integer, intent(in) :: unit, step
        real(real64)        :: kin_tot, tot_eng
        real(real64)        :: pot_av, kin_av, rot_av

        pot_av = pot_eng / dble( N )
        kin_av = kin_eng / dble( N )
        rot_av = kin_rot / dble( N )
        
        kin_tot = kin_av + rot_av
        tot_eng = pot_av + kin_tot

        write(unit, "((I0, 2x), 100(F0.10, 2x))") &
        step, pot_av, kin_tot, tot_eng, temp, pressure, density
    
    end subroutine report_state_data

    subroutine close_state_data_file( unit )
        implicit none
        integer, intent(in) :: unit

        close(unit=unit)

    end subroutine close_state_data_file

end module state_data_reporter

module trajectory_reporter
    use system
    implicit none
    character(len=256) :: file

    contains
    subroutine open_trajectory_file( unit, filename )
        implicit none
        integer, intent(in)          :: unit
        character(len=*), intent(in) :: filename

        file = filename

        open(unit=unit, file=file, action='write')

    end subroutine open_trajectory_file

    subroutine report_trajectory( unit, step )
        implicit none
        integer, intent(in) :: unit, step
        real(real64)        :: UX, UY, UZ
        integer             :: I

        write(unit, '(A)') "ITEM: TIMESTEP"
        write(unit, *) step
        write(unit, '(A)') "ITEM: NUMBER OF ATOMS"
        write(unit, *) n
        write(unit, "(A)") "ITEM: BOX BOUNDS pp pp pp"
        write(unit, "(1x, 2(F0.3, 2x))") -box(1) / 2.0, box(1) / 2.0
        write(unit, "(1x, 2(F0.3, 2x))") -box(2) / 2.0, box(2) / 2.0
        write(unit, "(1x, 2(F0.3, 2x))") -box(3) / 2.0, box(3) / 2.0
        write(unit, "(A)") "ITEM: ATOMS id type shapex shapey shapez x y z vx vy vz fx fy fz quatw quati quatj quatk"

        do I = 1, N

            UX = 2.0 *  ( QX(I) * QZ(I) + QW(I) * QY(I) )
            UY = 2.0 *  ( QY(I) * QZ(I) - QW(I) * QX(I) )
            UZ = QW(I)**2 - QX(I)**2 - QY(I)**2 + QZ(I)**2

            write(unit, "(1x, 2(I0, 2x), 50(F0.3, 2x))") I, rtype(I), sig_x(I), sig_y(I), sig_z(I), RX(I), RY(I), RZ(I), UX, UY, UZ, FX(I), FY(I), FZ(I), &
            QW(I), QX(I), QY(I), QZ(I)
        enddo

    end subroutine report_trajectory

    subroutine close_trajectory_file( unit )
        implicit none
        integer, intent(in) :: unit

        close(unit=unit)

    end subroutine close_trajectory_file

end module trajectory_reporter

module config_reporter
    use system
    implicit none
    character(len=256) :: file

    contains
    subroutine open_confwriter_file( unit, filename )
        implicit none
        integer, intent(in) :: unit
        character(len=*), intent(in) :: filename

        file = filename

        open(unit=unit, file=file, action='write')

    end subroutine open_confwriter_file

    subroutine report_confwriter( unit )
        implicit none
        integer, intent(in) :: unit
        integer :: I

        do I = 1, N
            write(unit, "(50(F0.3, 2x))") RX(I), RY(I), RZ(I), VX(I), VY(I), VZ(I), LX(I), LY(I), LZ(I), QW(I), QX(I), QY(I), QZ(I)
        enddo

    end subroutine report_confwriter

    subroutine close_confwriter_file( unit )
        implicit none
        integer, intent(in) :: unit

        close(unit=unit)

    end subroutine close_confwriter_file
    
end module config_reporter

module checkpoint_reporter
    use system
    implicit none
    character(len=256) :: file
    integer            :: config_index = 0

    contains
    subroutine open_checkpoint_file( unit, filename )
        implicit none
        integer, intent(in) :: unit
        character(len=*), intent(in) :: filename

        file = filename

        open(unit=unit, file=file)
        config_index = 0

    end subroutine open_checkpoint_file

    subroutine report_checkpoint( unit )
        implicit none
        integer, intent(in) :: unit
        integer :: I

        write(unit, "('ATOMS ', I0)") N
        write(unit, *)
        write(unit, "('BONDS ', I0)") NC
        write(unit, *)
        write(unit, "('BOX ', 3(F0.10, 2x))") box(1), box(2), box(3)
        write(unit, *)
        write(unit, "('ATOM: ID TYPE MASS CHARGE SHAPEX SHAPEY SHAPEZ ESHAPEX ESHAPEY ESHAPEZ IXX IYY IZZ RX RY RZ VX VY VZ LX LY LZ QW QX QY QZ')")
        write(unit, *)
        do I = 1, N
            write(unit, "(1x, 2(I0, 2x), 50(F0.10, 2x))") I, rtype(I), mass(I), charge(I), sig_x(I), sig_y(I), sig_z(I), eps_x(I), eps_y(I), eps_z(I), &
                                          Ixx(I), Iyy(I), Izz(I), RX(I), RY(I), RZ(I), VX(I), VY(I), VZ(I), LX(I), LY(I), LZ(I), QW(I), QX(I), QY(I), QZ(I)
        enddo
        write(unit, *)
        write(unit, "('BOND: ID I J LEN')")
        write(unit, *)
        do I = 1, NC
            write(unit, "(1x, 3(I0, 2x), 50(F0.10, 2x))") I, BBI(I), BBJ(I), BB_len(I)
        enddo
        write(unit, *)
        write(unit, *) config_index

    end subroutine report_checkpoint

    subroutine load_checkpoint( unit )
        implicit none
        integer, intent(in) :: unit
        integer             :: I, J, na, nb, ia, ja
        real(real64)        :: la, bx, by, bz, ma
        character(len=256)  :: label

        read(unit, *) label, na
        call set_num_atoms(na)
        read(unit, *)
        read(unit, *) label, nb
        call set_num_constraints(nb)
        read(unit, *)
        read(unit, *) label, bx, by, bz
        call set_box( bx, by, bz )
        read(unit, *)
        read(unit, *) label
        read(unit, *)
        do I = 1, N
            read(unit, *) J, rtype(I), ma, charge(I), sig_x(I), sig_y(I), sig_z(I), eps_x(I), eps_y(I), eps_z(I), &                  
                            Ixx(I), Iyy(I), Izz(I), RX(I), RY(I), RZ(I), VX(I), VY(I), VZ(I), LX(I), LY(I), LZ(I), QW(I), QX(I), QY(I), QZ(I)
            call set_atom_mass(J, ma)
        enddo
        read(unit, *)
        read(unit, *) label
        read(unit, *)
        do I = 1, NC
            read(unit, *) J, ia, ja, la
            call add_constraint(J, ia, ja, la)
        enddo
        read(unit, *)
        read(unit, *) config_index
        
    end subroutine load_checkpoint

    subroutine close_checkpoint_file( unit )
        implicit none
        integer, intent(in) :: unit

        close(unit=unit)

    end subroutine close_checkpoint_file
    
end module checkpoint_reporter

module reporters
    use state_data_reporter
    use trajectory_reporter
    use config_reporter
    use checkpoint_reporter
end module reporters
