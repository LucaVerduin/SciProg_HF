program HartreeFock

   ! Demonstration program that can be used as a starting point
   ! Lucas Visscher, March 2022

   use molecular_structure
   use ao_basis
   use compute_integrals
   use HartreeFock
   use diagonalization
   use InOut

     implicit none

     ! Variable naming as in the description of the exercise
     integer  :: n_AO, n_occ
     integer :: max_cycles
     real(8) :: tolerance
     character(32) :: filename
   
    filename = "testinput.txt"
    call getInput(filename)

     ! Definition of the molecule
     call define_molecule(molecule)

     ! Definition of the GTOs
     call define_basis(ao_basis)
     n_AO = ao_basis%nao
   
     ! Definition of the number of occupied orbitals
     n_occ = 3 ! hardwired for this demonstration program, should be set via input

    stop "stop here for generating input"

    call coreHamiltonian(n_AO, n_occ, molecule, ao_basis)

    tolerance = 0.0001d0
    max_cycles = 1000

    call SCFprocedure(n_AO, n_occ, max_cycles, tolerance)
   
     print*, "The Hartree-Fock energy:    ", E_HF

   end

   subroutine define_molecule(molecule)
     ! This routine should be improved such that an arbitrary molecule can be given as input
     ! the coordinates below are for a be-he dimer oriented along the x-axis with a bond length of 2 au
     use molecular_structure
     type(molecular_structure_t), intent(inout) :: molecule
     real(8) :: charge(2),coord(3,2)
     charge(1)   = 4.D0
     charge(2)   = 2.D0
     coord       = 0.D0
     coord(1,2)  = 2.D0
     call add_atoms_to_molecule(molecule,charge,coord)
   end subroutine

   subroutine define_basis(ao_basis)
    ! This routine can be extended to use better basis sets 
    ! The coordinates of the shell centers are the nuclear coordinates
    ! Think of a refactoring of define_molecule and define_basis to ensure consistency 
     use ao_basis
     type(basis_set_info_t), intent(inout) :: ao_basis
     type(basis_func_info_t) :: gto
     ! Be:  2 uncontracted s-funs:    l      coord          exp      
     call add_shell_to_basis(ao_basis,0,(/0.D0,0.D0,0.D0/),4.D0)
     call add_shell_to_basis(ao_basis,0,(/0.D0,0.D0,0.D0/),1.D0)
    ! Add 1 extra s function to Be for step 1
     call add_shell_to_basis(ao_basis,0,(/0.d0,0.d0,0.d0/),2.d0)
     ! He:  1 uncontracted s-fun:     l      coord          exp      
     call add_shell_to_basis(ao_basis,0,(/2.D0,0.D0,0.D0/),1.D0)
     ! Add 1 extra s function to H for step 1
     call add_shell_to_basis(ao_basis,0,(/2.d0,0.d0,0.d0/),2.d0)
   end subroutine


   
