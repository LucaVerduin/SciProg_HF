module HartreeFock
    implicit none
    private
    public calculateDensity

contains

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