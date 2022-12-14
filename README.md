[![Build Status](https://travis-ci.org/NLESC-JCER/Fortran_Davidson.svg?branch=master)](https://travis-ci.org/NLESC-JCER/Fortran_Davidson) [![Pipeline](https://github.com/NLESC-JCER/Fortran_Davidson/workflows/FortranCI/badge.svg)](https://github.com/NLESC-JCER/Fortran_Davidson/actions) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3637964.svg)](https://doi.org/10.5281/zenodo.3637964)

Davidson Eigensolver
===================
This package contains a Modern Fortran implementation of the *Davidson diagonalization algorithms*.
Different schemas are available to compute the correction.

Available correction methods are:
 * **DPR**: Diagonal-Preconditioned-Residue
 * **GJD**: Generalized Jacobi Davidson


### Note:
The Davidson method is suitable for **diagonal-dominant hermitian matrices**, that are quite common
in certain scientific problems like [electronic structure](https://en.wikipedia.org/wiki/Electronic_structure). The Davidson method could be not practical
for other kinds of hermitian matrices.

Usage
-----
The following program call the `eigensolver` *subroutine* from the `davidson` module and computes
the lowest 3 eigenvalues and corresponding eigenvectors, using the *GJD* method with a tolerance
of `1e-8` and `100` maximum iteration.
```fortran
program main
  use numeric_kinds, only: dp
  use davidson, only: generalized_eigensolver, generate_diagonal_dominant
 
  implicit none

  integer, parameter :: dim = 50
  integer, parameter :: lowest = 3
  complex(dp), dimension(dim, dim) :: matrix, second_matrix
  real(dp), dimension(lowest) :: eigenvalues
  complex(dp), dimension(dim, lowest) :: eigenvectors
  real(dp) :: tolerance
  integer:: max_dim_subspace, max_iterations, lowest

  matrix = generate_diagonal_dominant(dim, 1d-4)
  second_matrix = generate_diagonal_dominant(dim, 1d-4, 1.0_dp)
  max_iterations = 1000
  max_dim_subspace = 20
  tolerance = 1d-8
  call generalized_eigensolver(matrix, eigenvalues, eigenvectors, lowest, "GJD", max_iterations, &
       tolerance, final_iterations, max_dim_subspace, second_matrix)
  print *, eigenvalues
  print *, eigenvectors

end program main
```
The helper  `generate_diagonal_dominant` function generates a diagonal dominant
matrix with entries to the diagonal close to row number `(i=1, number_of_rows)`
and random number of the order `1e-4` on the off-diagonal entries.

**Variables**:
 * `matrix` (*in*) matrix to diagonalize
 * `eigenvalues` (*out*) resulting eigenvalues
 * `eigenvectors` (*out*) resulting eigenvectors
 * `lowest`(*in*) number of eigenvalues to compute
 * `method`(*in*) Either "DPR" or "GJD"
 * `max_iterations`(*in*) maximum number of iterations
 * `tolerance`(*in*) Numerical tolerance for convergence
 * `final_iterations`(*output*) returns the number of iterations that were needed to converge
 * `max_dim_subspace`(*in*, *optional*) Dimension of the subspace of search
 * `second_matrix`(*in*, optional) Optional matrix to compute the generalized eigenvalue problem
 
### References:
 * [Davidson diagonalization method and its applications to electronic structure calculations](https://www.semanticscholar.org/paper/DAVIDSON-DIAGONALIZATION-METHOD-AND-ITS-APPLICATION-Liao/5811eaf768d1a006f505dfe24f329874a679ba59)
 * [Numerical Methods for Large Eigenvalue Problem](https://doi.org/10.1137/1.9781611970739)

Installation and Testing
------------------------

To compile execute:
```
cmake -H. -Bbuild && cmake --build build
```

To use another compiler (e.g. ifort):
```
cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=ifort && cmake --build build
```

To run the test:
```
cmake -H. -Bbuild -DENABLE_TEST=ON && cmake --build build
cd build && ctest -V
```

To Debug compile as:
```
cmake -H. -Bbuild  -DCMAKE_BUILD_TYPE=Debug && cmake --build build
```

Dependencies
------------
This packages assumes that you have installed the following packages:
 * A Fortran compiler >=  version 7.0 
 * [CMake](https://cmake.org/)
 * [Lapack](http://www.netlib.org/lapack/)
	
Optionally, If an [MKL](https://software.intel.com/en-us/mkl) library is available the package will try to find it.
