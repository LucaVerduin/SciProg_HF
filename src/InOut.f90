! By Luca Verduin
! 03-2026

module inout
    use molecular_structure
    use ao_basis
    use atom_definition
    implicit none
    public getInput, ao_basis, molecule, defined_atoms, input_atoms, tolerance, max_cycles
    public generate_molecule, n_AO, n_occ
    public outfile, output_tofile, print_every
    private atom

    real(8) :: tolerance ! Tolerance criterium from input
    integer :: max_cycles ! max_cycles from input
    integer :: print_every ! Print every n cycles in SCF 
    integer :: n_AO, n_occ ! To be determined from system
    character(75) :: outfile ! Destination txt file for the output
    logical :: output_tofile ! Used to enable writing to output file (only if it is given in input file)

    ! Variable containing molecule data
    type(molecular_structure_t) :: molecule
    ! Variable containing the atomic orbital basis
    type(basis_set_info_t) :: ao_basis
    ! Variable containing definitions of atoms
    type(atom), allocatable, target :: defined_atoms(:)
    ! Variable containing atoms of the system
    type(atom_pointer), allocatable :: input_atoms(:)    

contains

subroutine getInput()
    character(75) :: filename
    character(300) :: line
    character(75) :: l_ang, exponents, temp_coordinates, dummy_char
    integer :: io, i, char, point_index
    integer :: n_defined_atoms, rewind_lines, index_atom, n_functions, n_atoms

    print *, "Enter input filename please."
    read *,filename
    filename = "inputs/"//trim(filename)

    io = 15                     ! Assign io value
    n_defined_atoms = 0         ! Used to count how much to allocate atom_definitions array
    n_atoms = 0                 ! Used to count how much to allocate input_atoms array
    index_atom = 0              ! internal counter to check which atom to set
    print_every = 5             ! default setting
    tolerance = 0.0001          ! default setting
    output_tofile = .false.     ! Logical to determine if output to file or terminal

    open(io, file=filename, status='old', action='read')

    do
        read(io, '(a)', end=10)line

        ! Lines starting with ! outside of "blocks" are ignored
        if (line(1:1) == "!") then
            continue
        end if

        ! If statement for defing atoms
        if (line(1:12)=="DEFINE ATOMS") then

            ! First find how many defined_atoms need to be allocated
            rewind_lines = 0 ! Count so instead of reopening just rewind # of lines
            do
                read(io, '(a)', end=101)line             
                rewind_lines = rewind_lines + 1
                if (line(1:13)=="/DEFINE ATOMS") then ! End with /DEFINE ATOMS
                    exit
                else if (line(1:1) == "!") then
                    continue
                else if (line(1:1) /= " ") then
                    n_defined_atoms = n_defined_atoms + 1 ! Non empty line count as atom definitions
                end if
            end do

            allocate(defined_atoms(n_defined_atoms)) ! Allocate array for defined atoms

            do i=1,rewind_lines ! Rewind back to first line after DEFINE ATOMS
                backspace(io)
            end do

            index_atom = 0 ! Count which atom we are assigning
            do i=1,rewind_lines ! Now now how many lines to read until /DEFINE ATOMS in principle
                read(io, '(a)')line
                if (line(1:13)=="/DEFINE ATOMS") then ! End with /DEFINE ATOMS
                    exit
                else if (line(1:1) == '!') then
                    continue
                else if (line(1:1)/=" ") then   ! Ignore lines starting with a space
                    index_atom = index_atom + 1 ! Increase index counter

                    ! Line should contain symbol, charge, n_electrons, array string for l, and array string for exponents.
                    read(line, *)defined_atoms(index_atom)%symbol, defined_atoms(index_atom)%charge, defined_atoms(index_atom)%n_electrons, l_ang, exponents

                    n_functions = 1 ! Count how many functions are given for this atom
                    do char=1,len(trim(l_ang))
                        if (l_ang(char:char)==" ") then
                            n_functions = n_functions + 1
                        end if
                    end do

                    ! Check if RHF even # electrons is obeyed
                    if (n_functions < defined_atoms(index_atom)%n_electrons/2) then
                        print '(a,1x,i3,2x,a)', "WARNING: amount of basis functions for atom", index_atom, "smaller than n_electrons/2"
                    end if

                    ! Allocate angular momenta and exponents for defined atom
                    allocate(defined_atoms(index_atom)%angular_momenta(n_functions))
                    allocate(defined_atoms(index_atom)%exponents(n_functions))
                    ! Read angular momenta and exponents for defined atom
                    read(l_ang, *)defined_atoms(index_atom)%angular_momenta
                    read(exponents, *)defined_atoms(index_atom)%exponents

                end if
            end do
        end if

        if (line(1:5) == "ATOMS") then
            rewind_lines = 0
            do
                ! First count how many atoms in system
                read(io, '(a)', end=102)line             
                    rewind_lines = rewind_lines + 1
                    if (line(1:6)=="/ATOMS") then ! End with /DEFINE ATOMS
                        exit
                    else if (line(1:1) == "!") then
                        continue
                    else if (line(1:1) /= " ") then ! Ignore empty lines
                        n_atoms = n_atoms + 1
                    end if
            end do

            ! Allocate room for atoms
            allocate(input_atoms(n_atoms))
            do i=1,rewind_lines ! Rewind back to first line after DEFINE ATOMS
                backspace(io)
            end do

            index_atom = 0 ! Count which atom is being read
            do
                read(io, '(a)')line
                if (line(1:6)=="/ATOMS") then ! If /ATOMS exit this loop
                    exit
                else if (line(1:1) == '!') then
                    continue
                else if (line(1:1) /= " ") then ! Ignore empty lines
                    index_atom = index_atom + 1
                    read(line, *)point_index, temp_coordinates ! Read temp_coordinates as text string
                    input_atoms(index_atom)%atom_type => defined_atoms(point_index) ! Point atom_type to the correct atom definition
                    read(temp_coordinates, *)input_atoms(index_atom)%coordinates ! Read from string the coordinates for the atom

                    ! Convert from A to Bohr
                    input_atoms(index_atom)%coordinates = input_atoms(index_atom)%coordinates*1.8897259886d0
                end if
            end do

        end if

        ! If statement for settings
        if (line(1:8) == "SETTINGS") then
            do
                read(io, '(a)', end=103)line
                if (line(1:9) == "/SETTINGS") then ! End this loop if /SETTINGS
                    exit
                else if (line(1:1) == "!") then
                    continue
                else if (line(1:9) == "tolerance") then ! Read tolerance
                    read(line, *)dummy_char, tolerance
                else if (line(1:9) == "maxcycles") then ! Read maxcycles
                    read(line, *)dummy_char, max_cycles
                else if (line(1:7) == "outfile") then   ! Read outputfile
                    output_tofile = .true.              ! Set output_tofile true
                    dummy_char = trim(line(8:))
                    read(dummy_char, *)outfile
                    outfile = "outputs/"//trim(outfile) ! Set output to ./outputs/ folder
                else if (line(1:10) == "printevery") then   ! Set printevery variable
                    read(line, *)dummy_char, print_every
                end if

            end do
        end if


    end do
    10 close(io)

    ! If output_tofile write set io to 10, else 6 for terminal
    if (output_tofile) then
        io = 10
        open(io, file=outfile, action='write')
    else
        io = 6
    end if

    ! Write input information
    write(io, '(a,/)')"INPUT --------"
    write(io, '(a)')"Atoms in system"

    do i=1,size(input_atoms)
        write(io, *)input_atoms(i)%atom_type%symbol, input_atoms(i)%coordinates/1.8897259886d0 ! Convert back to A (input) units
    end do

    write(io, *)""
    write(io, '(a,t17,e12.5)')"Tolerance",tolerance
    write(io, '(a,t18,i4)')"Max Cycles: ",max_cycles
    write(io, '(a,t16,2x,a)')"Output: ",outfile
    write(io, '(a,t15,i4)')"Print Every: ", print_every
    write(io, '(/,a)')"END INPUT ----"

    if (output_tofile) then
        close(io)
    end if

    return
    ! Errors for reaching end of file while reading a block (should not happen)
    101 stop "Error: block not ended with '/DEFINE ATOMS'"
    102 stop "Error: block not ended with '/ATOMS'"
    103 stop "Error: block not ended with '/SETTINGS'"

end subroutine

subroutine generate_molecule()
    real(8), allocatable :: charge(:), coords(:,:)
    real(8) :: n_electrons
    integer :: i, n_atoms

    ! get amount of atoms, allocate charge array and coords array
    n_atoms = size(input_atoms)
    allocate(charge(n_atoms))
    allocate(coords(3,n_atoms))

    ! set counter to 0
    n_electrons = 0

    ! loop over all input atoms, add charges to array and coords to arrays
    ! Also count # electrons
    do i=1,size(input_atoms)
        charge(i) = input_atoms(i)%atom_type%charge
        coords(:,i) = input_atoms(i)%coordinates
        n_electrons = n_electrons + input_atoms(i)%atom_type%n_electrons
    end do

    ! check if # electrons is even
    if (mod(n_electrons,2.0) == 0.0) then
        ! print '(a)', "Even amount of electrons"
    end if

    if (mod(n_electrons,2.0) /= 0.0) then
        ! print*, "WARNING amount of electrons in system is not even!"
        stop "Amount of electrons in system is not even!"
    end if

    ! set n_occ orbitals
    n_occ = nint(n_electrons/2.0)

    ! add atoms to molecule
    call add_atoms_to_molecule(molecule, charge, coords)

end subroutine

subroutine define_basis()
    integer :: atom, ifunc

    do atom = 1,size(input_atoms) ! Loop over all input atoms
        do ifunc = 1,size(input_atoms(atom)%atom_type%angular_momenta) ! Loop over all input angular momenta (all input functions in definition of %atom_Type)
            call add_shell_to_basis(ao_basis,&                                              ! Add to ao_basis
                                    input_atoms(atom)%atom_type%angular_momenta(ifunc),&    ! Angular momentum function
                                    input_atoms(atom)%coordinates,&                         ! At coordinates of atom
                                    input_atoms(atom)%atom_type%exponents(ifunc)&           ! With exponent
             )
        end do
    end do
    n_AO = ao_basis%nao

end subroutine

end module