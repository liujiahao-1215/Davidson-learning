program main
  use numeric_kinds, only: dp
  use davidson, only: generalized_eigensolver, generate_diagonal_dominant
  use test_utils, only: compute_vector_generalized_eigenvalue, compute_vector_on_fly, write_matrix, write_vector

  implicit none

  real(dp), dimension(3) :: eigenvalues_DPR
  real(dp), dimension(50, 3) :: eigenvectors_DPR
  real(dp), dimension(50, 50) :: mtx, stx
  integer :: iter_i
  
  ! mtx = read_matrix("tests/matrix.txt", 100)
  mtx = generate_diagonal_dominant(50, 1d-3)

    ! Matrix to check the algorithm
  mtx = generate_diagonal_dominant(50, 1d-3)
  stx = generate_diagonal_dominant(50, 1d-3, 1d0)
  call write_matrix("matrix_free.txt", mtx)
  call write_matrix("stx_free.txt", stx)

    ! NOTE:
  ! compute_vector_on_fly and compute_vector_generalized_eigenvalue JUST READ THE VALUES FROM
  ! THE PREVIOUS CREATED FILES!!!!!
  
  ! call eigenvalue solver
  call generalized_eigensolver(compute_vector_on_fly, eigenvalues_DPR, eigenvectors_DPR, 3, "DPR", 1000, &
       1d-8, iter_i, 20, compute_vector_generalized_eigenvalue)

  call write_vector("eigenvalues_DPR_free.txt",eigenvalues_DPR)
  call write_matrix("eigenvectors_DPR_free.txt", eigenvectors_DPR)

end program main
