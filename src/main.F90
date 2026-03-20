program HartreFock

   ! Demonstration program that can be used as a starting point
   ! Lucas Visscher, March 2022

   use molecular_structure
   use ao_basis
   use compute_integrals
   use HartreeFock
   use diagonalization
   use InOut

     implicit none
   
    ! filename = "CH4.txt"
    call getInput()
    call set_output(outfile, output_tofile)
    call generate_molecule()
    call define_basis()

    ! stop "STOP HERE FOR INPUT"

    call coreHamiltonian(n_AO, n_occ, molecule, ao_basis)

    ! stop "STOP HERE FOR CORE HAMILTONIAN"

    call SCFprocedure(n_AO, n_occ, max_cycles, tolerance, print_every)

    ! stop "STOP HERE FOR SCF"

   end
