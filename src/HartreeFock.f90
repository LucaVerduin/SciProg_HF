module HartreeFock
    use molecular_structure
    use ao_basis
    use compute_integrals
    use diagonalization
    implicit none
    private
    public SCFprocedure, coreHamiltonian, calculateDensity, H, F, V, T, S, C, eps, D, D_old, E_HF, delta_D, ao_integrals

real(8), allocatable :: H(:,:), F(:,:), V(:,:),T(:,:),S(:,:), C(:,:), eps(:), D(:,:), D_old(:,:)
real(8)  :: E_HF, delta_D
real(8), allocatable :: ao_integrals (:,:,:,:)

contains

subroutine coreHamiltonian(n_ao, n_occ, molecule, ao_basis)
    integer, intent(in) :: n_ao, n_occ
     type(molecular_structure_t), intent(in) :: molecule
     type(basis_set_info_t), intent(in) :: ao_basis
     integer :: kappa, lambda, mu, nu

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

end subroutine

subroutine SCFprocedure(n_AO, n_occ, max_cycles, tolerance)
    integer, intent(in) :: max_cycles, n_AO, n_occ
    real(8), intent(in) :: tolerance
    integer :: icycle, kappa, lambda, mu, nu

    icycle = 0

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
        print '(a,t13,i6,t20,a)', "Converged in", icycle, "cycles"
        exit
      end if 

    end do
    ! End SCF loop

end subroutine

subroutine calculateDensity(C, D, n_occ)
    real(8), intent(in) :: C(:,:)
    real(8), intent(out) :: D(:,:)
    integer, intent(in) :: n_occ
    integer :: kappa, lambda

    do lambda = 1, size(D,2)
        do kappa = 1, size(D,2)
            D(kappa,lambda) = sum(C(kappa,1:n_occ)*C(lambda,1:n_occ))
        end do
    end do

end subroutine   

end module