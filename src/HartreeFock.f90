module HartreeFock
    use molecular_structure
    use ao_basis
    use compute_integrals
    use diagonalization
    implicit none
    private write_tofile, output_file, print_every, io
    public SCFprocedure, coreHamiltonian, calculateDensity, H, F, V, T, S, C, eps, D, D_old, E_HF, delta_D, ao_integrals
    public set_output

real(8), allocatable :: H(:,:), F(:,:), V(:,:),T(:,:),S(:,:), C(:,:), eps(:), D(:,:), D_old(:,:)
real(8)  :: E_HF, delta_D
real(8), allocatable :: ao_integrals (:,:,:,:)
character(50) :: output_file, exit_status
logical :: write_tofile, is_converged
integer :: print_every
integer :: io = 10

contains

subroutine set_output(outfile, outtofile)
    character(50), intent(in) :: outfile
    logical, intent(in) :: outtofile
    
    write_tofile = outtofile
    output_file = outfile

end subroutine

subroutine coreHamiltonian(n_AO, n_occ, molecule, ao_basis)
    integer, intent(in) :: n_AO, n_occ
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
     
     open(io, file=output_file, status='old', access='append', action='write')
        write(io, '(2/,a,/)')"CALCULATION --------"
        write(io, '(a)')"Orbital energies for the core hamiltonian"
        write(io, '(8(f16.6))')eps
     close(io)

     ! Form the density matrix
     allocate (D(n_AO,n_AO))
     allocate (D_old(n_AO,n_AO))

     call calculateDensity(C,D,n_occ)

     allocate (ao_integrals(n_AO,n_AO,n_AO,n_AO))
     ! Compute all 2-electron integrals
     call generate_2int (ao_basis,ao_integrals)

end subroutine

subroutine SCFprocedure(n_AO, n_occ, max_cycles, tolerance, print_every)
    integer, intent(in) :: max_cycles, n_AO, n_occ
    real(8), intent(in) :: tolerance
    integer, intent(in) :: print_every
    integer :: icycle, kappa, lambda, mu, nu

    icycle = 0

    if (write_tofile) then
        open(io, file=output_file, status='old', access='append', action='write')
        write(io, '(/,a)')""
    end if

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

      ! Calculating delta_D
      delta_D = sqrt(sum( (D_old - D)**2 ))

      ! Checking for convergence
      if (delta_D < tolerance) then
        is_converged = .true.
        exit_status = "Converged"
      else if (icycle == max_cycles) then
        is_converged = .true.
        exit_status = "Max Iterations"
      end if
    
      if ( (mod(icycle, print_every)==0) .or. (icycle==max_cycles) ) then
        if (write_tofile) then
            write(io, '(a,2x,i4,2x,f14.8,a,f14.8)')"Energy cycle ", icycle, E_HF, " Ha       convergence: ", delta_D
        else
            print '(a,2x,i4,2x,f14.8,a,f14.8)', "Energy cycle ", icycle, E_HF, " Ha       convergence: ", delta_D 
        end if
      end if

    icycle = icycle + 1 ! Count # cycles

    if (is_converged) then
        if (write_tofile) then
            write(io, '(/,a)')"END CALCULATION -----"
            write(io, '(/,a,i4,a)')"Exited After ",icycle," SCF cycles"
            write(io, '(a,a)')"Exit status: ",exit_status
        else
            print '(/,a,i4,a)', "Exited After ",icycle," SCF cycles"
            print '(a,t15,a)',"Exit status: ", exit_status
        end if
        exit
    end if


    end do
    ! End SCF loop
    
    if (write_tofile) then
        close(io)
    end if

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