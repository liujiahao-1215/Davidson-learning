module array_utils

  use numeric_kinds, only: dp
  use lapack_wrapper, only: lapack_generalized_eigensolver, lapack_matmul, lapack_matrix_vector, &
    lapack_qr, lapack_solver, lapack_sort
  implicit none

  !> \private
  private
  !> \public
  public :: concatenate, diagonal,eye, generate_diagonal_dominant, norm, &
    generate_preconditioner

  interface norm
    !> compute the norm-2 of a vector
    procedure norm_r, norm_c
  end interface norm

contains

  pure function eye(m, n, alpha)
    !> Create a matrix with ones in the diagonal and zero everywhere else
    !> \param m: number of rows
    !> \param n: number of colums
    !> \param alpha: optional diagonal value
    !> \return matrix of size n x m
    integer, intent(in) :: n, m
    real(dp), dimension(m, n) :: eye
    real(dp), intent(in), optional :: alpha

    !local variable
    integer :: i, j
    real(dp) :: x

    ! check optional values
    x = 1.d0
    if (present(alpha)) x = alpha

    do i=1, m
      do j=1, n
        if (i /= j) then
          eye(i, j) = 0.d0
        else
          eye(i, i) = x
        end if
      end do
    end do

  end function eye

  pure function norm_c(vector)
    complex(dp), dimension(:), intent(in) :: vector
    real(dp) :: norm_c

    norm_c = sqrt(sum(abs(vector) ** 2.d0))

  end function norm_c

  pure function norm_r(vector)
    real(dp), dimension(:), intent(in) :: vector
    real(dp) :: norm_r

    norm_r = sqrt(sum(vector ** 2.d0))

  end function norm_r

  subroutine concatenate(arr, brr)

    !> Concatenate two matrices
    !> \param arr: first array
    !> \param brr: second array
    !> \return arr concatenate brr (overwrites arr)

    complex(dp), dimension(:, :), intent(inout), allocatable :: arr
    complex(dp), dimension(:, :), intent(in) :: brr
    complex(dp), dimension(:, :), allocatable :: tmp_array
    integer :: new_dim, dim_cols, dim_rows

    ! dimension
    dim_rows = size(arr, 1)
    dim_cols = size(arr, 2)

    ! Number of columns of the new matrix
    new_dim = dim_cols + size(brr, 2)

    ! move to temporal array
    allocate(tmp_array(dim_rows, new_dim))
    tmp_array(:, :dim_cols) = arr

    ! Move to new expanded matrix
    deallocate(arr)
    call move_alloc(tmp_array, arr)

    arr(:, dim_cols + 1:) = brr

  end subroutine concatenate

  function generate_diagonal_dominant(m, sparsity, diag_val) result(arr)
    !> Generate a diagonal dominant square matrix of dimension m
    !> \param m dimension of the matrix
    !> \param sparsity magnitude order of the off-diagonal values

    integer, intent(in) :: m ! size of the square matrix
    real(dp), optional :: diag_val
    integer :: i, j
    real(dp) :: sparsity
    complex(dp), dimension(m, m) :: arr

    ! Local variables
    real(dp), dimension(m, m) :: rea, img

    call random_number(rea)
    call random_number(img)

    arr = cmplx(rea, img, dp) * sparsity
    do j=1, m
      do i=1, m
        if (i > j) then
          arr(i, j) = conjg(arr(j, i))
        else if(i == j) then
          if (present(diag_val))then
            arr(i, i) = diag_val
          else
            arr(i, i) = i
          end if
        end if
      end do
    end do

  end function generate_diagonal_dominant

  function diagonal(matrix)
    !> return the diagonal of a matrix
    complex(dp), dimension(:, :), intent(in) :: matrix
    complex(dp), dimension(size(matrix, 1)) :: diagonal

    ! local variables
    integer :: i, m

    ! dimension of the matrix
    m = size(matrix, 1)

    do i=1,m
      diagonal(i) = matrix(i, i)
    end do

  end function diagonal

  function generate_preconditioner(diag, dim_sub) result(precond)
    !> \brief generates a diagonal preconditioner for `matrix`.
    !> \return diagonal matrix

    ! input variable
    complex(dp), dimension(:), intent(inout) :: diag
    integer, intent(in) :: dim_sub

    ! local variables
    complex(dp), dimension(size(diag), dim_sub) :: precond
    real(dp), dimension(size(diag)) :: r_diag
    logical(dp), dimension(size(diag)) :: mask
    integer :: i, k

    ! Initialize the norm of diag and mask for minloc
    r_diag = abs(diag)
    mask = .true.

    ! Fill matrix with zeros
    precond = 0.0_dp

    ! Add one depending on the order of the matrix diagonal
    do i=1, dim_sub
      k = minloc(r_diag, dim=1, mask=mask)
      mask(k) = .false.
      precond(k, i) = 1.d0
    end do

  end function generate_preconditioner

end module array_utils
