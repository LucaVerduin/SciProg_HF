module inout
    use molecular_structure
    use ao_basis
    implicit none
    public getInput, ao_basis, molecule
    private atom

    type(molecular_structure_t) :: molecule
    ! Variable containing the atomic orbital basis
    type(basis_set_info_t) :: ao_basis

    type atom
        character(4) :: symbol
        real(8) :: charge
        real(8) :: n_electrons
        integer, allocatable :: angular_momenta(:)
        real(8), allocatable :: exponents(:)
    end type

contains

subroutine getInput(filename)
    character(32), intent(in) :: filename
    character(250) :: line
    character(20) :: l_ang
    character(20) :: exponents
    integer :: io, i, char
    integer :: n_defined_atoms, rewind_lines, index_atom, n_functions
    type(atom), allocatable :: defined_atoms(:)

    io = 15
    n_defined_atoms = 0
    index_atom = 0
    open(io, file=filename, status='old', action='read')

    do
        read(io, '(a)', end=10)line
    
        if (line(1:12)=="DEFINE ATOMS") then
            rewind_lines = 0
            do
                read(io, '(a)')line             
                rewind_lines = rewind_lines + 1
                if (line(1:13)=="/DEFINE ATOMS") then
                    exit
                else if (line(1:1) /= " ") then
                    n_defined_atoms = n_defined_atoms + 1
                end if
            end do

            allocate(defined_atoms(n_defined_atoms))

            do i=1,rewind_lines
                backspace(io)
            end do

            index_atom = 0
            do i=1,rewind_lines
                read(io, '(a)')line
                if (line(1:13)=="/DEFINE ATOMS") then
                    exit
                else if (line(1:1)/=" ") then
                    index_atom = index_atom + 1

                    read(line, *)defined_atoms(index_atom)%symbol, defined_atoms(index_atom)%charge, defined_atoms(index_atom)%n_electrons, l_ang, exponents


                    n_functions = 1
                    do char=1,len(trim(l_ang))
                        if (l_ang(char:char)==" ") then
                            n_functions = n_functions + 1
                        end if
                    end do

                    ! print *, i, n_functions
                    allocate(defined_atoms(index_atom)%angular_momenta(n_functions))
                    allocate(defined_atoms(index_atom)%exponents(n_functions))
                    read(l_ang, *)defined_atoms(index_atom)%angular_momenta
                    read(exponents, *)defined_atoms(index_atom)%exponents

                end if
            end do

        end if

        print *, n_defined_atoms, "Defined atoms"
        print *, defined_atoms(1)%symbol
        print *, defined_atoms(1)%charge
        print *, defined_atoms(1)%n_electrons
        print *, defined_atoms(1)%angular_momenta
        print *, defined_atoms(1)%exponents
    end do
    10 close(io)

end subroutine

end module