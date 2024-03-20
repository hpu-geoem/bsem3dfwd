// ----------------------------------------------------------------------------
//
// Copyright (C) 2023 by The authors
//
// This file is part of the BSEM3DFWD program.
//
// This code is licensed under the MIT license found in the
// LICENSE file in the root directory of this source tree.
//
// ----------------------------------------------------------------------------

#ifndef _EM_DEFS_H_
#define _EM_DEFS_H_ 1

#include <complex>
#include <tuple>

#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/petsc_block_sparse_matrix.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

typedef std::complex<double> Complex;

typedef dealii::Point<3> Point;
typedef dealii::Tensor<1, 3, double> VectorD;
typedef dealii::Tensor<1, 3, Complex> VectorZ;

typedef dealii::DoFHandler<3> DoFHandler;
typedef dealii::AffineConstraints<double> AffineConstaints;

typedef dealii::Triangulation<3> Triangulation;
typedef dealii::parallel::shared::Triangulation<3> SharedTriangulation;

typedef dealii::PETScWrappers::MPI::Vector PETScVector;
typedef dealii::PETScWrappers::MPI::BlockVector PETScBlockVector;
typedef dealii::PETScWrappers::MPI::SparseMatrix PETScSparseMatrix;
typedef dealii::PETScWrappers::MPI::BlockSparseMatrix PETScBlockSparseMatrix;

const double PI = 3.14159265358979323846;
const double MU = 4 * PI * 1E-7;
const Complex II = Complex(0.0, 1.0);

const double RTOD = 180.0 / PI;
const double DTOR = PI / 180.0;

const double EPS = 1E-6;

enum InnerPCType { Mixed = 0, AMS = 1, Direct = 2 };

enum DirectSolverType { MUMPS = 0, SUPERLU_DIST = 1, CPARDISO = 2 };

struct Transmitter {
  Point center;
  double azimuth, dip, length, current;
};

#define EM_ERR_USER -100000

#endif
