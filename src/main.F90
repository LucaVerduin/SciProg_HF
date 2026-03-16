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
     character(32) :: filename
   
    filename = "CH4.txt"
    call getInput(filename)
    call generate_molecule()
    call define_basis()
    print '(3/,a)', ""

    call coreHamiltonian(n_AO, n_occ, molecule, ao_basis)

    call SCFprocedure(n_AO, n_occ, max_cycles, tolerance)
   
     print*, "The Hartree-Fock energy:    ", E_HF

   end
