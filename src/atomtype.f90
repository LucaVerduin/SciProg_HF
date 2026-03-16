module atom_definition
    implicit none
    public atom, atom_pointer
    private

    type atom
        character(4) :: symbol
        real(8) :: charge
        real(8) :: n_electrons
        integer, allocatable :: angular_momenta(:)
        real(8), allocatable :: exponents(:)
    end type

    type atom_pointer
        type(atom), pointer :: atom_type
        real(8) :: coordinates(3)
    end type

contains


end module