module system
    use iso_fortran_env, only: real64
    use xmath
    use philox_rng
    implicit none

    real(real64), allocatable, dimension(:) :: mass, charge, inv_mass
    real(real64), allocatable, dimension(:) :: sig_x, sig_y, sig_z
    real(real64), allocatable, dimension(:) :: eps_x, eps_y, eps_z
    real(real64), allocatable, dimension(:) :: Ixx, Iyy, Izz
    real(real64), allocatable, dimension(:) :: Ixx_inv, Iyy_inv, Izz_inv
    integer, allocatable, dimension(:)      :: rtype

    real(real64), allocatable, dimension(:) :: RX, RY, RZ
    real(real64), allocatable, dimension(:) :: VX, VY, VZ
    real(real64), allocatable, dimension(:) :: FX, FY, FZ

    real(real64), allocatable, dimension(:) :: QW, QX, QY, QZ
    real(real64), allocatable, dimension(:) :: LX, LY, LZ
    real(real64), allocatable, dimension(:) :: TX, TY, TZ

    real(real64), allocatable, dimension(:) :: RX_old, RY_old, RZ_old

    real(real64)    :: box(3)
    logical         :: periodic

    integer         :: N, NC, Nrigid, dof
    real(real64)    :: pot_eng, kin_eng, kin_rot, temp
    real(real64)    :: pressure, density, virial, virial_bonded
    real(real64)    :: target_temp, target_pres

    integer, dimension(:), allocatable      :: BBI, BBJ
    real(real64), dimension(:), allocatable :: BB_len
    
    logical :: num_atoms_set        = .FALSE.
    logical :: fcc_lattice_set      = .FALSE.
    logical :: remove_com_flag      = .FALSE.
    logical :: target_temp_set      = .FALSE.
    logical :: target_pres_set      = .FALSE.

    contains

    subroutine set_num_atoms( num_atoms, rigid_atoms )
        implicit none
        integer, intent(in) :: num_atoms
        integer, intent(in), optional :: rigid_atoms

        if ( present(rigid_atoms) ) then
            Nrigid = rigid_atoms
            N = rigid_atoms * num_atoms
            NC = 0
            dof = 6 * N
        else
            Nrigid = 1
            N = num_atoms
            NC = 0
            dof = 6 * N
        endif
        num_atoms_set = .TRUE.

        allocate( mass(N), inv_mass(N), charge(N), rtype(N) )
        allocate( sig_x(N), sig_y(N), sig_z(N) )
        allocate( eps_x(N), eps_y(N), eps_z(N) )
        allocate( Ixx(N), Iyy(N), Izz(N) )
        allocate( Ixx_inv(N), Iyy_inv(N), Izz_inv(N) )

        allocate( RX(N), RY(N), RZ(N) )
        allocate( VX(N), VY(N), VZ(N) )
        allocate( FX(N), FY(N), FZ(N) )
        allocate( LX(N), LY(N), LZ(N) )
        allocate( TX(N), TY(N), TZ(N) )
        allocate( QW(N), QX(N), QY(N), QZ(N) )

        mass = -1.0; inv_mass = -1.0; charge = 0.0; rtype = -1;
        sig_x = -1.0; sig_y = -1.0; sig_z = -1.0;
        eps_x = -1.0; eps_y = -1.0; eps_z = -1.0
        Ixx = -1.0; Iyy = -1.0; Izz = -1.0;
        Ixx_inv = -1.0; Iyy_inv = -1.0; Izz_inv = -1.0

        RX = 0.0; RY = 0.0; RZ = 0.0
        VX = 0.0; VY = 0.0; VZ = 0.0
        FX = 0.0; FY = 0.0; FZ = 0.0
        LX = 0.0; LY = 0.0; LZ = 0.0
        TX = 0.0; TY = 0.0; TZ = 0.0
        QW = 1.0; QX = 0.0; QY = 0.0; QZ = 0.0

    end subroutine set_num_atoms

    subroutine set_num_rigid_atoms( num_rigid )
        implicit none
        integer, intent(in) :: num_rigid

        Nrigid = num_rigid

    end subroutine set_num_rigid_atoms

    subroutine set_num_constraints( num_constraints )
        implicit none
        integer, intent(in) :: num_constraints

        NC = num_constraints
        
        allocate( RX_old(N), RY_old(N), RZ_old(N) )
        allocate( BBI(NC), BBJ(NC), BB_len(NC) )
        BBI = 0; BBJ = 0; BB_len = 0.0
        RX_old = 0.0; RY_old = 0.0; RZ_old = 0.0

    end subroutine set_num_constraints

    subroutine set_atom_type( id, atom_type )
        implicit none
        integer, intent(in) :: id
        integer, intent(in) :: atom_type

        if (num_atoms_set .eqv. .FALSE.) then
            write(*, "('Error: System has no atoms. Set number of atoms first.')")
            stop
        endif

        if ( id < 1 .or. id > N ) then
            write(*, "('Error: Atom index ', I0,' does not exist in system of ',I0,' atoms.')") id, N
            stop
        endif
        
        rtype(id) = atom_type

    end subroutine set_atom_type

    subroutine set_atom_mass( id, atom_mass )
        implicit none
        integer, intent(in)      :: id
        real(real64), intent(in) :: atom_mass

        if (num_atoms_set .eqv. .FALSE.) then
            write(*, "('Error: System has no atoms. Set number of atoms first.')")
			stop
        endif

        if ( id < 1 .or. id > N ) then
            write(*, "('Error: Atom index ', I0,' does not exist in system of ',I0,' atoms.')") id, N
            stop
        endif

        mass(id) = atom_mass

    end subroutine set_atom_mass

    subroutine set_atom_charge( id, atom_charge )
        implicit none
        integer, intent(in)      :: id
        real(real64), intent(in) :: atom_charge

        if (num_atoms_set .eqv. .FALSE.) then
            write(*, "('Error: System has no atoms. Set number of atoms first.')")
            stop
        endif

        if ( id < 1 .or. id > N ) then
            write(*, "('Error: Atom index ', I0,' does not exist in system of ',I0,' atoms.')") id, N
            stop
        endif
        
        charge(id) = atom_charge

    end subroutine set_atom_charge

    subroutine set_density( dens )
        implicit none
        real(real64), intent(in) :: dens
        real(real64)             :: vol, box_len

        if (dens == 0.0) then
            write(*, "('Error: Density cannot be zero.')")
            stop
        endif

        density = dens
        vol = DBLE( N / Nrigid ) / density
        box_len = vol ** (1.0/3.0)

        box = box_len

    end subroutine set_density

    subroutine set_box( length, width, height )
        implicit none
        real(real64), intent(in)           :: length
        real(real64), optional, intent(in) :: width, height

        if (length == 0.0 .or. width == 0.0 .or. height == 0.0) then
            write(*, "('Error: Box dimensions cannot be zero.')")
            stop
        endif

        if ( present(width) .and. present(height) ) then
            box(1) = length
            box(2) = width
            box(3) = height
        else
            box(1) = length
            box(2) = length
            box(3) = length
        endif

        density = DBLE( N / Nrigid ) / product( box )

    end subroutine set_box

    subroutine set_periodic_boundary
        implicit none

        periodic = .TRUE.
    end subroutine set_periodic_boundary

    subroutine set_target_temperature( tmp_val )
        implicit none
        real(real64), intent(in) :: tmp_val

        target_temp = tmp_val
        target_temp_set = .TRUE.

    end subroutine set_target_temperature

    subroutine set_target_pressure( pres_val )
        implicit none
        real(real64), intent(in) :: pres_val

        target_pres = pres_val

    end subroutine set_target_pressure

    subroutine set_atom_shape( id, length, width, height )
        implicit none
        integer, intent(in)      :: id
        real(real64), intent(in)  :: length, width, height

        if (num_atoms_set .eqv. .FALSE.) then
            write(*, "('Error: System has no atoms. Set number of atoms first.')")
            stop
        endif

        if ( id < 1 .or. id > N ) then
            write(*, "('Error: Atom index ', I0,' does not exist in system of ',I0,' atoms.')") id, N
            stop
        endif

        if (length == 0.0 .or. width == 0.0 .or. height == 0.0) then
            write(*, "('Error: Shape dimensions cannot be zero.')")
            stop
        endif

        sig_x(id) = length
        sig_y(id) = width
        sig_z(id) = height

    end subroutine set_atom_shape

    subroutine set_atom_energy( id, eng_x, eng_y, eng_z )
        implicit none
        integer, intent(in)      :: id
        real(real64), intent(in) :: eng_x, eng_y, eng_z

        if (num_atoms_set .eqv. .FALSE.) then
            write(*, "('Error: System has no atoms. Set number of atoms first.')")
            stop
        endif

        if ( id < 1 .or. id > N ) then
            write(*, "('Error: Atom index ', I0,' does not exist in system of ',I0,' atoms.')") id, N
            stop
        endif

        eps_x(id) = eng_x
        eps_y(id) = eng_y
        eps_z(id) = eng_z
        
    end subroutine set_atom_energy

    subroutine set_atom_position(id, pos_x, pos_y, pos_z)
        implicit none

        integer, intent(in)         :: id
        real(real64), intent(in)    :: pos_x, pos_y, pos_z

        if (num_atoms_set .eqv. .FALSE.) then
            write(*, "('Error: System has no atoms. Set number of atoms first.')")
            stop
        endif

        if ( id < 1 .or. id > N ) then
            write(*, "('Error: Atom index ', I0,' does not exist &
            in system of ',I0,' atoms.')") id, N
            stop
        endif
        
        RX(id) = pos_x
        RY(id) = pos_y
        RZ(id) = pos_z

    end subroutine set_atom_position

    subroutine set_atom_velocity(id, vel_x, vel_y, vel_z)
        implicit none

        integer, intent(in)         :: id
        real(real64), intent(in)    :: vel_x, vel_y, vel_z

        if (num_atoms_set .eqv. .FALSE.) then
            write(*, "('Error: System has no atoms. Set number of atoms first.')")
            stop
        endif

        if ( id < 1 .or. id > N ) then
            write(*, "('Error: Atom index ', I0,' does not exist in system of ',I0,' atoms.')") id, N
            stop
        endif
        
        VX(id) = vel_x
        VY(id) = vel_y
        VZ(id) = vel_z

    end subroutine set_atom_velocity

    subroutine set_atom_inertia( id, inx, iny, inz )
        implicit none
        integer, intent(in)         :: id
        real(real64), intent(in)    :: inx, iny, inz

        if (num_atoms_set .eqv. .FALSE.) then
            write(*, "('Error: System has no atoms. Set number of atoms first.')")
            stop
        endif

        if ( id < 1 .or. id > N ) then
            write(*, "('Error: Atom index ', I0,' does not exist in system of ',I0,' atoms.')") id, N
            stop
        endif

        IXX(id) = inx
        IYY(id) = iny
        IZZ(id) = inz

    end subroutine set_atom_inertia

    subroutine add_constraint(id, I, J, length)
        implicit none
        integer, intent(in)         :: id
        integer, intent(in)         :: I, J
        real(real64), intent(in)    :: length

        if ( id < 1 .or. id > NC ) then
            write(*, "('Error: Constraint index ', I0,' does not exist in system of ',I0,' constraints.')") id, NC
            stop
        endif

        BBI(id) = I
        BBJ(id) = J
        BB_len(id) = length

    end subroutine add_constraint

    subroutine set_velocity_to_temperature( rtemp, seed, remove_com_velocity )
        implicit none

        real(real64), intent(in)      :: rtemp
        integer, intent(in)           :: seed
        logical, optional, intent(in) :: remove_com_velocity

        real(real64)    :: factor
        real(real64)    :: vx_gen, vy_gen, vz_gen
        real(real64)    :: vx_com, vy_com, vz_com
        integer         :: IX, n_dof, temp_int

        if (present(remove_com_velocity)) then
            remove_com_flag = remove_com_velocity
        else
            remove_com_flag = .TRUE.
        endif

        if ( num_atoms_set .eqv. .FALSE. ) then
            write(*, "('Error: System has no atoms. Set number of atoms first.')")
            stop
        endif

        vx_gen = 0.0; vy_gen = 0.0; vz_gen = 0.0
        vx_com = 0.0; vy_com = 0.0; vz_com = 0.0
        temp_int = int( rtemp * 1000 )

        do IX = 1, N

            call philox_normal_3seq( ix, temp_int, seed, 42, vx_gen, vy_gen, vz_gen )

            VX(IX) = vx_gen
            VY(IX) = vy_gen
            VZ(IX) = vz_gen

        enddo

        n_dof = 6 * N - NC

        if (remove_com_flag) then

            vx_com = sum( vx ) / dble( N ) 
            vy_com = sum( vy ) / dble( N ) 
            vz_com = sum( vz ) / dble( N ) 

            VX = VX - vx_com
            VY = VY - vy_com
            VZ = VZ - vz_com

            n_dof = n_dof - 3

        endif

        factor = sqrt( dble(n_dof) * rtemp * Nrigid / sum( ( VX ** 2 + VY ** 2 + VZ ** 2 ) ) )

        VX = VX * factor
        VY = VY * factor
        VZ = VZ * factor

    end subroutine set_velocity_to_temperature

    subroutine initialize_forces
        implicit none

        FX = 0.0; FY = 0.0; FZ = 0.0
        TX = 0.0; TY = 0.0; TZ = 0.0

        pot_eng  = 0.0
        kin_eng  = 0.0
        kin_rot  = 0.0
        temp     = 0.0
        pressure = 0.0
        virial   = 0.0
        virial_bonded = 0.0
    end subroutine initialize_forces

    subroutine generate_fcc_lattice( unit_cells, dens, rigid_atoms )
        implicit none
        integer, intent(in)             :: unit_cells
        real(real64), intent(in)        :: dens
        integer, intent(in), optional   :: rigid_atoms
        integer                         :: num_atoms, mtemp, I, J, K, IREF
        real(real64)                    :: cell, half_cell, rroot3
        fcc_lattice_set = .TRUE.
        num_atoms = 4 * unit_cells ** 3

        if ( num_atoms_set .eqv. .FALSE. ) then
            if ( present(rigid_atoms) ) then
                call set_num_atoms( num_atoms, rigid_atoms )
            else
                call set_num_atoms( num_atoms )
            endif
        else if ( N < num_atoms ) then
            write(*, "('Error: FCC lattice of ', I0,' atoms &
                        &cannot be generated for system of ',I0,' atoms.')") num_atoms, N
            stop
        endif
            
        call set_density( dens )

        cell = box(1) / real(unit_cells, real64)
        half_cell = cell / 2.0
        rroot3 = 1.0 / sqrt(3.0)

        RX(1) = 0.0
        RY(1) = 0.0
        RZ(1) = 0.0
        QW(1) = sqrt( ( 1.0 + rroot3 ) / 2.0 )
        QX(1) = sqrt( ( 1.0 - rroot3 ) / 2.0 ) * (  rroot3 / sqrt( 1.0 - rroot3 ** 2.0 ) )
        QY(1) = sqrt( ( 1.0 - rroot3 ) / 2.0 ) * ( -rroot3 / sqrt( 1.0 - rroot3 ** 2.0 ) )
        QZ(1) = 0.0

        RX(2) = half_cell
        RY(2) = half_cell
        RZ(2) = 0.0
        QW(2) = sqrt( ( 1.0 - rroot3 ) / 2.0 )
        QX(2) = sqrt( ( 1.0 + rroot3 ) / 2.0 ) * ( -rroot3 / sqrt( 1.0 - rroot3 ** 2.0 ) )
        QY(2) = sqrt( ( 1.0 + rroot3 ) / 2.0 ) * ( -rroot3 / sqrt( 1.0 - rroot3 ** 2.0 ) )
        QZ(2) = 0.0

        RX(3) = 0.0
        RY(3) = half_cell
        RZ(3) = half_cell
        QW(3) = sqrt( ( 1.0 - rroot3 ) / 2.0 )
        QX(3) = sqrt( ( 1.0 + rroot3 ) / 2.0 ) * (  rroot3 / sqrt( 1.0 - rroot3 ** 2.0 ) )
        QY(3) = sqrt( ( 1.0 + rroot3 ) / 2.0 ) * (  rroot3 / sqrt( 1.0 - rroot3 ** 2.0 ) )
        QZ(3) = 0.0

        RX(4) = half_cell
        RY(4) = 0.0
        RZ(4) = half_cell
        QW(4) = sqrt( ( 1.0 + rroot3 ) / 2.0 )
        QX(4) = sqrt( ( 1.0 - rroot3 ) / 2.0 ) * ( -rroot3 / sqrt( 1.0 - rroot3 ** 2.0 ) )
        QY(4) = sqrt( ( 1.0 - rroot3 ) / 2.0 ) * (  rroot3 / sqrt( 1.0 - rroot3 ** 2.0 ) )
        QZ(4) = 0.0

        mtemp = 0
        do I = 1, unit_cells
            do J = 1, unit_cells
                do K = 1, unit_cells
                    do iref = 1, 4

                        RX( IREF + MTEMP ) = RX( IREF ) + cell * ( I - 1 )
                        RY( IREF + MTEMP ) = RY( IREF ) + cell * ( J - 1 )
                        RZ( IREF + MTEMP ) = RZ( IREF ) + cell * ( K - 1 )

                        QW( IREF + MTEMP ) = QW( IREF )
                        QX( IREF + MTEMP ) = QX( IREF )
                        QY( IREF + MTEMP ) = QY( IREF )
                        QZ( IREF + MTEMP ) = QZ( IREF )

                    enddo
                    mtemp = mtemp + 4
                enddo
            enddo
        enddo

        RX(1:num_atoms) = RX(1:num_atoms) - box(1) / 2.0
        RY(1:num_atoms) = RY(1:num_atoms) - box(2) / 2.0
        RZ(1:num_atoms) = RZ(1:num_atoms) - box(3) / 2.0

    end subroutine generate_fcc_lattice

    subroutine initialize_system( v )
        implicit none
        logical, optional :: v
        logical           :: verbose
        integer           :: IX, count_mass, count_charge, count_type
        integer           :: count_shape, count_energy, count_inertia
        real(real64)      :: temp_dens

        if ( .not. present(v) ) then
            verbose = .TRUE.
        else
            verbose = v
        endif

        if(verbose) write(*, *)
        if(verbose) write(*, "('This is SUCROSE v2.0')")

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (num_atoms_set) then
            if(verbose) write(*, "('Number of particles in system set to ', I0)") N
        else
            if(verbose) write(*, "('Error: Number of particles not set. Set number of &
            particles before initializing system.')")
            stop
        endif

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (any(box .le. 0)) then
            if(verbose) write(*, "('Error: Box dimensions not set. Set box dimensions &
            before initializing system.')")
            stop
        else
            if(verbose) write(*, "('Box dimensions set to ', F0.2, 3X, F0.2, 3X, F0.2)") box

            if (verbose) write(*, "('Density of system is ', F0.5)") density
        endif

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        dof = 6 * N - NC
        if (remove_com_flag) dof = dof - 3
        if (verbose) write(*, "('Degrees of freedom in system set to ', I0)") dof

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (periodic) then
            if(verbose) write(*, "('Periodic boundary conditions ENABLED')")
        else
            if(verbose) write(*, "('Periodic boundary conditions DISABLED')")
        endif

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (fcc_lattice_set) then
            if(verbose) write(*, "('FCC lattice generated with ', I0, ' particles')") N / Nrigid
        endif

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        count_type = 0; count_mass = 0; count_charge = 0
        count_shape = 0; count_energy = 0; count_inertia = 0
        do IX = 1, N

            if (rtype(IX) == -1 .or. rtype(IX) < 0) then
                rtype(IX) = 0
                count_type = count_type + 1
            endif

            if (charge(IX) == -1 .or. charge(IX) < 0) then
                charge(IX) = 0.0
                count_charge = count_charge + 1   
            endif
            
            if (sig_x(IX) == -1 .or. sig_x(IX) < 0) then
                sig_x(IX) = 1.0
                sig_y(IX) = 1.0
                sig_z(IX) = 1.0
                count_shape = count_shape + 1
            endif
            
            if (eps_x(IX) == -1 .or. eps_x(IX) < 0) then
                eps_x(IX) = 1.0
                eps_y(IX) = 1.0
                eps_z(IX) = 1.0
                count_energy = count_energy + 1
            endif
            
            if (mass(IX) == -1 .or. mass(IX) < 0) then
                mass(IX) = 1.0
                inv_mass(IX) = 1.0

                if (Ixx(IX) < 0 .or. Iyy(IX) < 0 .or. Izz(IX) < 0) then
                    Ixx(IX) = ( sig_y(IX) ** 2.0 + sig_z(IX) ** 2.0 ) / 20.0
                    Iyy(IX) = ( sig_x(IX) ** 2.0 + sig_z(IX) ** 2.0 ) / 20.0
                    Izz(IX) = ( sig_x(IX) ** 2.0 + sig_y(IX) ** 2.0 ) / 20.0

                    Ixx_inv(IX) = 1.0 / Ixx(IX)
                    Iyy_inv(IX) = 1.0 / Iyy(IX)
                    Izz_inv(IX) = 1.0 / Izz(IX)

                    count_inertia = count_inertia + 1

                else if (Ixx(IX) > 0 .and. Iyy(IX) > 0 .and. Izz(IX) > 0) then
                    Ixx_inv(IX) = 1.0 / Ixx(IX)
                    Iyy_inv(IX) = 1.0 / Iyy(IX)
                    Izz_inv(IX) = 1.0 / Izz(IX)
                
                else
                    Ixx_inv(IX) = 0.0
                    Iyy_inv(IX) = 0.0
                    Izz_inv(IX) = 0.0
                endif

                
                count_mass = count_mass + 1

            else if (mass(IX) == 0) then
                inv_mass(IX) = 0.0

                Ixx(IX) = 0.0
                Iyy(IX) = 0.0
                Izz(IX) = 0.0

                Ixx_inv(IX) = 0.0
                Iyy_inv(IX) = 0.0
                Izz_inv(IX) = 0.0

            else
                inv_mass(IX) = 1.0 / mass(IX)

                if (Ixx(IX) < 0 .or. Iyy(IX) < 0 .or. Izz(IX) < 0) then
                    Ixx(IX) = ( sig_y(IX) ** 2.0 + sig_z(IX) ** 2.0 ) / 20.0
                    Iyy(IX) = ( sig_x(IX) ** 2.0 + sig_z(IX) ** 2.0 ) / 20.0
                    Izz(IX) = ( sig_x(IX) ** 2.0 + sig_y(IX) ** 2.0 ) / 20.0

                    Ixx_inv(IX) = 1.0 / Ixx(IX)
                    Iyy_inv(IX) = 1.0 / Iyy(IX)
                    Izz_inv(IX) = 1.0 / Izz(IX)

                    count_inertia = count_inertia + 1
                
                else if (Ixx(IX) > 0 .and. Iyy(IX) > 0 .and. Izz(IX) > 0) then
                    Ixx_inv(IX) = 1.0 / Ixx(IX)
                    Iyy_inv(IX) = 1.0 / Iyy(IX)
                    Izz_inv(IX) = 1.0 / Izz(IX)
                
                else
                    Ixx_inv(IX) = 0.0
                    Iyy_inv(IX) = 0.0
                    Izz_inv(IX) = 0.0
                endif


            endif

        enddo

        if (verbose) then
            if (count_mass > 0) then
                write(*, "('Default mass 1.0 set for ', I0, ' atoms.')") count_mass
            endif
            if (count_charge > 0) then
                write(*, "('Default charge 0.0 set for ', I0, ' atoms.')") count_charge
            endif
            if (count_shape > 0) then
                write(*, "('Default shape (1.0, 1.0, 1.0) set for ', I0, ' atoms.')") count_shape
            endif
            if (count_energy > 0) then
                write(*, "('Default energy (1.0, 1.0, 1.0) set for ', I0, ' atoms.')") count_energy
            endif
            if (count_inertia > 0) then
                write(*, "('Default inertia (1.0, 1.0, 1.0) set for ', I0, ' atoms.')") count_inertia
            endif
        endif

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (Nrigid > 1) then
            target_temp = target_temp * Nrigid

            temp_dens = dble( N / Nrigid ) / product( box )

            if (abs( temp_dens - density ) > 0.01) then
                if(verbose) write(*, "('Density of system changed to ', F0.5, ' from ', F0.5)") temp_dens, density
                call set_density( temp_dens )
            endif

        endif
        if(verbose) write(*, "('Target temperature set to ', F0.5)") target_temp

        
        call initialize_forces

        if(verbose) write(*, "('System initialization complete.')")
        if(verbose) write(*, *)

    end subroutine initialize_system

    subroutine calculate_thermo()
        implicit none
        real(real64) :: volume, kin_tot

        kin_tot = kin_eng + kin_rot
        temp = 2.0 * kin_tot / dble(dof)

        if (NC > 0 .and. virial_bonded >= 0.0) then
            virial = SUM( RX * FX + RY * FY + RZ * FZ ) + virial_bonded
        else
            virial = SUM( RX * FX + RY * FY + RZ * FZ )
        endif

        volume = product( box )
        pressure = (density * Nrigid * temp) + (virial / (3.0 * volume) )

    end subroutine calculate_thermo

    subroutine copy_transform( from, to, xoffset, yoffset, zoffset, set_inertia )
        implicit none
        integer, intent(in)      :: from, to
        real(real64), intent(in)         :: xoffset, yoffset, zoffset
        real(real64)                     :: rot(3, 3), rot_offset(3)
        logical, optional        :: set_inertia

        rot(1, 1) = QW(from)**2 + QX(from)**2 - QY(from)**2 - QZ(from)**2
        rot(1, 2) = 2.0 * ( QX(from) * QY(from) - QW(from) * QZ(from) )
        rot(1, 3) = 2.0 * ( QX(from) * QZ(from) + QW(from) * QY(from) )
        rot(2, 1) = 2.0 * ( QX(from) * QY(from) + QW(from) * QZ(from) )
        rot(2, 2) = QW(from)**2 - QX(from)**2 + QY(from)**2 - QZ(from)**2
        rot(2, 3) = 2.0 * ( QY(from) * QZ(from) - QW(from) * QX(from) )
        rot(3, 1) = 2.0 * ( QX(from) * QZ(from) - QW(from) * QY(from) )
        rot(3, 2) = 2.0 * ( QY(from) * QZ(from) + QW(from) * QX(from) )
        rot(3, 3) = QW(from)**2 - QX(from)**2 - QY(from)**2 + QZ(from)**2

        rot_offset(1) = rot(1, 1) * xoffset + rot(1, 2) * yoffset + rot(1, 3) * zoffset
        rot_offset(2) = rot(2, 1) * xoffset + rot(2, 2) * yoffset + rot(2, 3) * zoffset
        rot_offset(3) = rot(3, 1) * xoffset + rot(3, 2) * yoffset + rot(3, 3) * zoffset

        RX(to) = RX(from) + rot_offset(1)
        RY(to) = RY(from) + rot_offset(2)
        RZ(to) = RZ(from) + rot_offset(3)
		
		QW(to) = QW(from)
        QX(to) = QX(from)
        QY(to) = QY(from)
        QZ(to) = QZ(from)

        if (present(set_inertia)) then

            if (set_inertia) then

                if (mass(to) <= 0.0) then

                    if (Ixx(from) > 0.0 .and. Iyy(from) > 0.0 .and. Izz(from) > 0.0 ) then

                        if (Ixx(to) > 0.0 .and. Iyy(to) > 0.0 .and. Izz(to) > 0.0 ) then

                            Ixx(from) = Ixx(from) + Ixx(to) + ( yoffset**2.0 + zoffset**2.0 )
                            Iyy(from) = Iyy(from) + Iyy(to) + ( xoffset**2.0 + zoffset**2.0 )
                            Izz(from) = Izz(from) + Izz(to) + ( xoffset**2.0 + yoffset**2.0 )

                        else if ( sig_x(to) > 0.0 .and. sig_y(to) > 0.0 .and. sig_z(to) > 0.0 ) then

                            Ixx(to) = ( sig_y(to) ** 2.0 + sig_z(to) ** 2.0 ) / 20.0
                            Iyy(to) = ( sig_x(to) ** 2.0 + sig_z(to) ** 2.0 ) / 20.0
                            Izz(to) = ( sig_x(to) ** 2.0 + sig_y(to) ** 2.0 ) / 20.0

                            Ixx(from) = Ixx(from) + Ixx(to) + ( yoffset**2.0 + zoffset**2.0 )
                            Iyy(from) = Iyy(from) + Iyy(to) + ( xoffset**2.0 + zoffset**2.0 )
                            Izz(from) = Izz(from) + Izz(to) + ( xoffset**2.0 + yoffset**2.0 )
                        
                        else
                            write(*, "('Error: Second particle ', I0, 'must have non-zero shape')") to
                            stop
                        endif

                    else if (sig_x(from) > 0.0 .and. sig_y(from) > 0.0 .and. sig_z(from) > 0.0 ) then

                        Ixx(from) = ( sig_y(from) ** 2.0 + sig_z(from) ** 2.0 ) / 20.0
                        Iyy(from) = ( sig_x(from) ** 2.0 + sig_z(from) ** 2.0 ) / 20.0
                        Izz(from) = ( sig_x(from) ** 2.0 + sig_y(from) ** 2.0 ) / 20.0

                        if (Ixx(to) > 0.0 .and. Iyy(to) > 0.0 .and. Izz(to) > 0.0 ) then

                            Ixx(from) = Ixx(from) + Ixx(to) + ( yoffset**2.0 + zoffset**2.0 )
                            Iyy(from) = Iyy(from) + Iyy(to) + ( xoffset**2.0 + zoffset**2.0 )
                            Izz(from) = Izz(from) + Izz(to) + ( xoffset**2.0 + yoffset**2.0 )

                        else if ( sig_x(to) > 0.0 .and. sig_y(to) > 0.0 .and. sig_z(to) > 0.0 ) then

                            Ixx(to) = ( sig_y(to) ** 2.0 + sig_z(to) ** 2.0 ) / 20.0
                            Iyy(to) = ( sig_x(to) ** 2.0 + sig_z(to) ** 2.0 ) / 20.0
                            Izz(to) = ( sig_x(to) ** 2.0 + sig_y(to) ** 2.0 ) / 20.0

                            Ixx(from) = Ixx(from) + Ixx(to) + ( yoffset**2.0 + zoffset**2.0 )
                            Iyy(from) = Iyy(from) + Iyy(to) + ( xoffset**2.0 + zoffset**2.0 )
                            Izz(from) = Izz(from) + Izz(to) + ( xoffset**2.0 + yoffset**2.0 )
                        
                        else
                            write(*, "('Error: Second particle ', I0, 'must have non-zero shape')") to
                            stop
                        endif

                    else
                        write(*, "('Error: First particle ', I0, 'must have non-zero shape')") from
                        stop
                    endif

                else
                    write(*, "('Error: Second particle ', I0, 'must have mass zero')") to
                    stop
                
                endif
            endif
        endif

    end subroutine copy_transform

    subroutine transfer_forces( from, to )
        implicit none
        integer, intent(in) :: from, to
        real(real64)                :: rij(3), rij_pbc(3)
        real(real64)                :: fij_from(3), tij_to(3)

        rij(1) = RX(from) - RX(to)
        rij(2) = RY(from) - RY(to)
        rij(3) = RZ(from) - RZ(to)

        rij_pbc = rij - ANINT( rij / box ) * box
        rij = MERGE( rij_pbc, rij, periodic )
        
        fij_from(1) = FX(from)
        fij_from(2) = FY(from)
        fij_from(3) = FZ(from)
        
        tij_to(1) = rij(2) * fij_from(3) - rij(3) * fij_from(2)
        tij_to(2) = rij(3) * fij_from(1) - rij(1) * fij_from(3)
        tij_to(3) = rij(1) * fij_from(2) - rij(2) * fij_from(1)

        FX(to) = FX(to) + FX(from)
        FY(to) = FY(to) + FY(from)
        FZ(to) = FZ(to) + FZ(from)

        TX(to) = TX(to) + TX(from) + tij_to(1)
        TY(to) = TY(to) + TY(from) + tij_to(2)
        TZ(to) = TZ(to) + TZ(from) + tij_to(3)

        FX(from) = 0.0
        FY(from) = 0.0
        FZ(from) = 0.0

        TX(from) = 0.0
        TY(from) = 0.0
        TZ(from) = 0.0

    end subroutine transfer_forces

end module system

module neighbour_list
    use system
    implicit none
    real(real64), allocatable, private :: x_last(:), y_last(:), z_last(:)
    real(real64), private              :: skin_width, r_skin_sq, global_cutoff

    integer, allocatable, public :: nlist(:), offset(:)
    integer, public              :: nl_size, nb_max

    contains
    subroutine init_nlist(skin, cutoff)
        implicit none
        real(real64), intent(in) :: skin, cutoff
        integer                  :: cap

		global_cutoff = cutoff
		
        if (any(global_cutoff > box / 2.0)) then
            write(*, "('A cutoff of ', F0.2, ' is too large &
                &for the box dimensions ', 3(F0.2, 2x))") global_cutoff, box
            stop
        endif

        skin_width = skin
        r_skin_sq  = skin ** 2

        if (.not. allocated(nlist)) then
            cap = 50 * N
            nl_size = cap
            allocate(nlist(cap))
        end if

        allocate( offset(N) )
        allocate( x_last(N), y_last(N), z_last(N) )

        offset = 1; nlist = 1;
        x_last = 0.0; y_last = 0.0; z_last = 0.0

        call build_init_nlist
		
    end subroutine init_nlist

    subroutine build_init_nlist()
        implicit none
        real(real64) :: RIJ(3), RIJ_PBC(3), rij_sq, rskin
        integer      :: I, J, K

        rskin = global_cutoff + skin_width
        r_skin_sq = rskin ** 2

        K = 0
        do I = 1, N-1

            offset(I) = K + 1
            do J = I+1, N
                
                RIJ(1) = RX(I) - RX(J)
                RIJ(2) = RY(I) - RY(J)
                RIJ(3) = RZ(I) - RZ(J)

                RIJ_PBC = RIJ - ANINT( RIJ / box ) * box
                RIJ = MERGE( RIJ_PBC, RIJ, periodic )

                rij_sq = SUM( RIJ * RIJ )
                
                if ( rij_sq < r_skin_sq ) then
                    K = K + 1
                    if (K > nl_size) call resize_nlist()
                    nlist(K) = J
                    
                end if
            end do

            nb_max = max(nb_max, K - offset(I) + 1)
        end do

        offset(N) = K + 1
        x_last = RX; y_last = RY; z_last = RZ

    end subroutine build_init_nlist

    subroutine build_nlist()
        implicit none
        real(real64) :: RIJ(3), RIJ_PBC(3), rij_sq, rskin
        integer      :: I, J, K

        rskin = global_cutoff + skin_width
        r_skin_sq = rskin ** 2

        K = 0
        do I = 1, N-1

            offset(I) = K + 1
            do J = I+1, N
                
                RIJ(1) = RX(I) - RX(J)
                RIJ(2) = RY(I) - RY(J)
                RIJ(3) = RZ(I) - RZ(J)

                RIJ_PBC = RIJ - ANINT( RIJ / box ) * box
                RIJ = MERGE( RIJ_PBC, RIJ, periodic )

                rij_sq = SUM( RIJ * RIJ )
                
                if ( rij_sq < r_skin_sq ) then
                    K = K + 1
                    if (K > nl_size) call resize_nlist()
                    nlist(K) = J
                    
                end if
            end do

            nb_max = max(nb_max, K - offset(I) + 1)
        end do

        offset(N) = K + 1
        x_last = RX; y_last = RY; z_last = RZ

    end subroutine build_nlist

    logical function needs_rebuild()
        implicit none
        real(real64)    :: dx, dy, dz
        real(real64)    :: dx_pbc, dy_pbc, dz_pbc
        real(real64)    :: max_disp, r_disp_sq
        integer         :: I

        max_disp = 0.0

        do I = 1, N

            dx = RX(I) - X_last(I)
            dy = RY(I) - Y_last(I)
            dz = RZ(I) - Z_last(I)

            dx_pbc = dx - ANINT( dx / box(1) ) * box(1)
            dy_pbc = dy - ANINT( dy / box(2) ) * box(2)
            dz_pbc = dz - ANINT( dz / box(3) ) * box(3)

            dx = MERGE( dx_pbc, dx, periodic )
            dy = MERGE( dy_pbc, dy, periodic )
            dz = MERGE( dz_pbc, dz, periodic )

            r_disp_sq = dx*dx + dy*dy + dz*dz
            max_disp = max( max_disp, r_disp_sq )

        enddo

        if (4.0 * max_disp > r_skin_sq) then
            needs_rebuild = .true.
            return
        end if

        needs_rebuild = .false.
    end function needs_rebuild

    subroutine resize_nlist()
        implicit none
        integer, allocatable :: tmp(:)
        integer              :: new_nl_size

        new_nl_size = ceiling( 2.0 * real(nl_size, real64) )
        allocate(tmp(new_nl_size))
        tmp(1:nl_size) = nlist
        call move_alloc(tmp, nlist)
        nl_size = new_nl_size
    end subroutine resize_nlist

end module neighbour_list

module constraints
    use system
    implicit none

    contains
    subroutine apply_constraints_a( DT )
        implicit none
        real(real64), intent(in)    :: DT
        real(real64)                :: RIJ(3), RIJ_PBC(3), RIJ_old(3), RIJ_old_pbc(3)
        real(real64)                :: L_ij, chi, vdot, r_tol, len
        integer                     :: I, J, IX, rattle, max_rattle = 100
        logical                     :: moved(N)

        L_ij = 0.0; rattle = 0; 
        r_tol = 1.0e-5; max_rattle = 100
        moved = .TRUE.

        do while (any(moved) .and. rattle < max_rattle)

            do IX = 1, NC
                I   = BBI(IX)
                J   = BBJ(IX)
                len = BB_len(IX)
                
                if ( moved(I) .or. moved(J) ) then

                    RIJ(1) = RX(I) - RX(J)
                    RIJ(2) = RY(I) - RY(J)
                    RIJ(3) = RZ(I) - RZ(J)

                    RIJ_pbc = RIJ - ANINT( RIJ / box ) * box
                    RIJ = MERGE( RIJ_pbc, RIJ, periodic )
                    
                    chi = len ** 2 - SUM( RIJ ** 2 )

                    if (abs(chi) > 2 * r_tol * len**2.0) then 

                        RIJ_old(1) = RX_old(I) - RX_old(J)
                        RIJ_old(2) = RY_old(I) - RY_old(J)
                        RIJ_old(3) = RZ_old(I) - RZ_old(J)

                        RIJ_old_pbc = RIJ_old - ANINT( RIJ_old / box ) * box
                        RIJ_old = MERGE( RIJ_old_pbc, RIJ_old, periodic )

                        vdot = RIJ(1) * RIJ_old(1) + RIJ(2) * RIJ_old(2) + RIJ(3) * RIJ_old(3)

                        if (vdot < r_tol * len ** 2.0) then
                            write(*, "(5F10.5)") chi, vdot, len, sum(RIJ**2), sum(RIJ_old**2)
                            write(*, "('Constraint failure between atoms ', I0, ' and ', I0)") I, J
                            stop
                        endif

                        L_ij = chi / ( 2.0 * vdot * ( inv_mass(I) + inv_mass(J) ) )

                        virial_bonded = virial_bonded + ( L_ij * len ** 2 ) / DT ** 2

                        RX(I) = RX(I) + ( L_ij * inv_mass(I) ) * RIJ_old(1)
                        RY(I) = RY(I) + ( L_ij * inv_mass(I) ) * RIJ_old(2)
                        RZ(I) = RZ(I) + ( L_ij * inv_mass(I) ) * RIJ_old(3)

                        RX(J) = RX(J) - ( L_ij * inv_mass(J) ) * RIJ_old(1)
                        RY(J) = RY(J) - ( L_ij * inv_mass(J) ) * RIJ_old(2)
                        RZ(J) = RZ(J) - ( L_ij * inv_mass(J) ) * RIJ_old(3)

                        VX(I) = VX(I) + ( L_ij * inv_mass(I) ) * RIJ_old(1) / DT
                        VY(I) = VY(I) + ( L_ij * inv_mass(I) ) * RIJ_old(2) / DT
                        VZ(I) = VZ(I) + ( L_ij * inv_mass(I) ) * RIJ_old(3) / DT

                        VX(J) = VX(J) - ( L_ij * inv_mass(J) ) * RIJ_old(1) / DT
                        VY(J) = VY(J) - ( L_ij * inv_mass(J) ) * RIJ_old(2) / DT
                        VZ(J) = VZ(J) - ( L_ij * inv_mass(J) ) * RIJ_old(3) / DT

                        moved(i) = .TRUE.
                        moved(j) = .TRUE.
                    else
                        moved(i) = .FALSE.
                        moved(j) = .FALSE.
                    endif

                endif

            enddo

            rattle = rattle + 1

        enddo

    end subroutine apply_constraints_a

    subroutine apply_constraints_b( DT )
        implicit none
        real(real64), intent(in)    :: DT
        real(real64)                :: RIJ(3), RIJ_pbc(3), VIJ(3)
        real(real64)                :: L_ij, vdot, r_tol, len
        integer                     :: I, J, IX, rattle, max_rattle = 100
        logical                     :: moved(N)

        L_ij = 0.0; rattle = 0; 
        r_tol = 1.0e-5; max_rattle = 100
        moved = .TRUE.

        do while (any(moved) .and. rattle < max_rattle)

            do IX = 1, NC
                I   = BBI(IX)
                J   = BBJ(IX)
                len = BB_len(IX)

                if( moved(I) .or. moved(J) ) then

                    RIJ(1) = RX(I) - RX(J)
                    RIJ(2) = RY(I) - RY(J)
                    RIJ(3) = RZ(I) - RZ(J)

                    RIJ_pbc = RIJ - ANINT( RIJ / box ) * box
                    RIJ = MERGE( RIJ_pbc, RIJ, periodic )
                    
                    VIJ(1) = VX(I) - VX(J)
                    VIJ(2) = VY(I) - VY(J)
                    VIJ(3) = VZ(I) - VZ(J)

                    vdot = RIJ(1) * VIJ(1) + RIJ(2) * VIJ(2) + RIJ(3) * VIJ(3)
                    L_ij = -vdot / ( ( inv_mass(I) + inv_mass(J) ) * len**2 )
                    
                    if (abs(L_ij) > r_tol) then

                        virial_bonded = virial_bonded + ( L_ij * len**2 ) / DT

                        VX(I) = VX(I) + ( L_ij * inv_mass(I) ) * RIJ(1)
                        VY(I) = VY(I) + ( L_ij * inv_mass(I) ) * RIJ(2)
                        VZ(I) = VZ(I) + ( L_ij * inv_mass(I) ) * RIJ(3)

                        VX(J) = VX(J) - ( L_ij * inv_mass(J) ) * RIJ(1)
                        VY(J) = VY(J) - ( L_ij * inv_mass(J) ) * RIJ(2)
                        VZ(J) = VZ(J) - ( L_ij * inv_mass(J) ) * RIJ(3)

                        moved(i) = .TRUE.
                        moved(j) = .TRUE.

                    else

                        moved(i) = .FALSE.
                        moved(j) = .FALSE.

                    endif

                endif
                
            enddo

            rattle = rattle + 1

        enddo

    end subroutine apply_constraints_b
    
end module constraints