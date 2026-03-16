module inout
    use molecular_structure
    use ao_basis
    use atom_definition
    implicit none
    public getInput, ao_basis, molecule, defined_atoms, input_atoms, tolerance, max_cycles
    public generate_molecule, n_AO, n_occ
    private atom

    real(8) :: tolerance ! Tolerance criterium from input
    integer :: max_cycles ! max_cycles from input
    integer :: n_AO, n_occ ! To be determined from system

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
                end if
            end do
        end if


    end do
    10 close(io)

    do i=1,size(input_atoms)
        print *, input_atoms(i)%atom_type%symbol
        print *, input_atoms(i)%coordinates
    end do

    print *, tolerance, max_cycles

end subroutine

subroutine generate_molecule()

end subroutine

end module