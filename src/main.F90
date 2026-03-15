program HartreeFock

   ! Demonstration program that can be used as a starting point
   ! Lucas Visscher, March 2022

   use molecular_structure
   use ao_basis
   use compute_integrals
   use HartreeFock
   use diagonalization

     implicit none

     ! Variable containing the molecular structure
     type(molecular_structure_t) :: molecule
     ! Variable containing the atomic orbital basis
     type(basis_set_info_t) :: ao_basis

     ! Variable naming as in the description of the exercise
     integer  :: n_AO, n_occ
     integer  :: kappa, lambda, mu, nu
     integer :: max_cycles, icycle
     real(8)  :: E_HF
     real(8), allocatable :: H(:,:), F(:,:), V(:,:),T(:,:),S(:,:), C(:,:), eps(:), D(:,:), D_old(:,:)
     real(8) :: delta_D, tolerance

     ! The following large array can be eliminated when Fock matrix contruction is implemented
     real(8), allocatable :: ao_integrals (:,:,:,:)
   
     ! Definition of the molecule
     call define_molecule(molecule)

     ! Definition of the GTOs
     call define_basis(ao_basis)
     n_AO = ao_basis%nao
   
     ! Definition of the number of occupied orbitals
     n_occ = 3 ! hardwired for this demonstration program, should be set via input

     ! Compute the overlap matrix
     allocate (S(n_AO,n_AO))
     call   compute_1e_integrals ("OVL",ao_basis,ao_basis,S)

     ! Compute the kinetic matrix
     allocate (T(n_AO,n_AO))
     call   compute_1e_integrals ("KIN",ao_basis,ao_basis,T)

     ! Compute the potential matrix
     allocate (V(n_AO,n_AO))
     call   compute_1e_integrals ("POT",ao_basis,ao_basis,V,molecule)

     ! Compute the core Hamiltonian matrix (the potential is positive, we scale with -e = -1 to get to the potential energy matrix)
     allocate (F(n_AO,n_AO))
     allocate (H(n_AO,n_AO))
     H = T - V

     ! Diagonalize the Fock matrix
     allocate (C(n_AO,n_AO))
     allocate (eps(n_AO))
     call solve_genev (H,S,C,eps)
     print*, "Orbital energies for the core Hamiltonian:",eps

     ! Form the density matrix
     allocate (D(n_AO,n_AO))
     allocate (D_old(n_AO,n_AO))

     call calculateDensity(C,D,n_occ)

     allocate (ao_integrals(n_AO,n_AO,n_AO,n_AO))
     ! Compute all 2-electron integrals
     call generate_2int (ao_basis,ao_integrals)



    tolerance = 0.0001d0
    icycle = 0
    max_cycles = 1000

    ! SCF loop
    do
      
      ! Calculating Fock matrix from D
      do lambda = 1, n_ao
        do kappa = 1, n_ao
              F(kappa, lambda) = H(kappa, lambda) + sum( (2.d0*ao_integrals(:,:,kappa,lambda) - ao_integrals(:,lambda,kappa,:))*D )
        end do
      end do

      ! Saving copy of previous density matrix and calculating new density matrix
      D_old = D
      call solve_genev(F,S,C,eps)
      call calculateDensity(C,D,n_occ)

      ! Calculate energy of this itteration
      E_HF = sum(H*D)         ! Hcore part
      E_HF = E_HF + sum(F*D)  ! F part

      print '(a,2x,i4,2x,f14.8,a)', "Energy cycle ", icycle, E_HF, " Ha" 

      icycle = icycle + 1 ! Count # cycles

      ! Calculating delta_D and checking for convergence
      delta_D = sqrt(sum( (D_old - D)**2 ))
      if (delta_D < tolerance .or. icycle==max_cycles) then
        print '(a,t13,i4,t18,a)', "Converged in", icycle, "cycles"
        exit
      end if 

    end do
    ! End SCF loop
   
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


   
