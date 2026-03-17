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
    integer :: print_every ! Determine every how many cycles to print in scf
    integer :: n_AO, n_occ ! To be determined from system
    character(50) :: outfile ! Destination txt file for the output
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

subroutine getInput(filename)
    character(32), intent(in) :: filename
    character(300) :: line
    character(75) :: l_ang, exponents, temp_coordinates, dummy_char
    integer :: io, i, char, point_index
    integer :: n_defined_atoms, rewind_lines, index_atom, n_functions, n_atoms

    io = 15
    n_defined_atoms = 0
    n_atoms = 0
    index_atom = 0
    print_every = 5
    output_tofile = .false.

    open(io, file=filename, status='old', action='read')

    do
        read(io, '(a)', end=10)line
        if (line(1:1) == "!") then
            continue
        end if

        ! If statement for defing atoms
        if (line(1:12)=="DEFINE ATOMS") then

            ! First find how many defined_atoms need to be allocated
            rewind_lines = 0 ! Count so instead of reopening just rewind # of lines
            do
                read(io, '(a)')line             
                rewind_lines = rewind_lines + 1
                if (line(1:13)=="/DEFINE ATOMS") then ! End with /DEFINE ATOMS
                    exit
                else if (line(1:1) /= " ") then
                    n_defined_atoms = n_defined_atoms + 1
                end if
            end do

            allocate(defined_atoms(n_defined_atoms)) ! Allocate array for defined atoms

            do i=1,rewind_lines ! Rewind back to first line after DEFINE ATOMS
                backspace(io)
            end do

            index_atom = 0 ! Count which atom we are assigning
            do i=1,rewind_lines
                read(io, '(a)')line
                if (line(1:13)=="/DEFINE ATOMS") then ! End with /DEFINE ATOMS
                    exit
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

                    if (n_functions < defined_atoms(index_atom)%n_electrons/2) then
                        print '(a,1x,i3,2x,a)', "WARNING: amount of basis functions for atom", index_atom, "smaller than n_electrons/2"
                    end if

                    ! print *, i, n_functions
                    allocate(defined_atoms(index_atom)%angular_momenta(n_functions))
                    allocate(defined_atoms(index_atom)%exponents(n_functions))
                    read(l_ang, *)defined_atoms(index_atom)%angular_momenta
                    read(exponents, *)defined_atoms(index_atom)%exponents

                end if
            end do
        end if

        if (line(1:5) == "ATOMS") then
            rewind_lines = 0
            do
                read(io, '(a)')line             
                    rewind_lines = rewind_lines + 1
                    if (line(1:6)=="/ATOMS") then ! End with /DEFINE ATOMS
                        exit
                    else if (line(1:1) /= " ") then
                        n_atoms = n_atoms + 1
                    end if
            end do

            allocate(input_atoms(n_atoms))
            do i=1,rewind_lines ! Rewind back to first line after DEFINE ATOMS
                backspace(io)
            end do

            index_atom = 0
            do
                read(io, '(a)')line
                if (line(1:6)=="/ATOMS") then
                    exit
                else if (line(1:1) /= " ") then
                    index_atom = index_atom + 1
                    read(line, *)point_index, temp_coordinates
                    input_atoms(index_atom)%atom_type => defined_atoms(point_index)
                    read(temp_coordinates, *)input_atoms(index_atom)%coordinates
                end if
            end do

        end if

        if (line(1:8) == "SETTINGS") then
            do
                read(io, '(a)')line
                if (line(1:9) == "/SETTINGS") then
                    exit
                else if (line(1:9) == "tolerance") then
                    read(line, *)dummy_char, tolerance
                else if (line(1:9) == "maxcycles") then
                    read(line, *)dummy_char, max_cycles
                else if (line(1:7) == "outfile") then
                    output_tofile = .true.
                    dummy_char = trim(line(8:))
                    read(dummy_char, *)outfile
                else if (line(1:10) == "printevery") then
                    read(line, *)dummy_char, print_every
                end if

            end do
        end if


    end do
    10 close(io)

    if (output_tofile) then
        open(io, file=outfile, action='write')
            write(io, '(a,/)')"INPUT --------"
            write(io, '(a)')"Atoms in system"
            do i=1,size(input_atoms)
                write(io, *)input_atoms(i)%atom_type%symbol, input_atoms(i)%coordinates
            end do
            write(io, *)""
            write(io, '(a,t17,e12.5)')"Tolerance",tolerance
            write(io, '(a,t18,i4)')"Max Cycles: ",max_cycles
            write(io, '(a,t16,2x,a)')"Output: ",outfile
            write(io, '(a,t15,i4)')"Print Every: ", print_every
            write(io, '(/,a)')"END INPUT ----"
        close(io)
        print '(i4, a)', size(input_atoms), "  Atoms in system"
        print '(i4, a)', size(defined_atoms), "  Defined atom types"
    else
        print '(a,/)', "Atoms in system:"
        do i=1,size(input_atoms)
            print *, input_atoms(i)%atom_type%symbol, input_atoms(i)%coordinates
        end do
        print *, ""
        print '(a,e12.5)', "Tolerance: ",tolerance
        print '(a,i4)', "Max Cycles: ", max_cycles
        print '(a,t14,2x,a)', "Output: ",outfile
    end if


end subroutine

subroutine generate_molecule()
    real(8), allocatable :: charge(:), coords(:,:)
    real(8) :: n_electrons
    integer :: i, n_atoms

    n_atoms = size(input_atoms)
    allocate(charge(n_atoms))
    allocate(coords(3,n_atoms))
    n_electrons = 0

    do i=1,size(input_atoms)
        charge(i) = input_atoms(i)%atom_type%charge
        coords(:,i) = input_atoms(i)%coordinates
        n_electrons = n_electrons + input_atoms(i)%atom_type%n_electrons
    end do

    if (mod(n_electrons,2.0) == 0.0) then
        print '(/,a)', "Even amount of electrons"
    end if

    if (mod(n_electrons,2.0) /= 0.0) then
        print*, "WARNING amount of electrons in system is not even!"
        stop "Amount of electrons in system is not even!"
    end if

    n_occ = nint(n_electrons/2.0)

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