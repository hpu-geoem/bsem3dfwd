// ----------------------------------------------------------------------------
//
// Copyright (C) 2023 by Ce Qin <ce.qin@hpu.edu.cn>
//
// This file is part of the BSEM3DFWD program.
//
// This code is licensed under the MIT license found in the
// LICENSE file in the root directory of this source tree.
//
// ----------------------------------------------------------------------------

#ifndef _EM_CTX_H_
#define _EM_CTX_H_ 1

#include <deal.II/fe/fe_nedelec.h>

#include "em_defs.h"


/**
 * @brief The EMContext struct contains all the necessary data for the BSEM forward problem.
 */
struct EMContext {
  MPI_Comm world_comm;             // MPI communicator
  PetscInt world_size, world_rank; // MPI world size and rank

  std::vector<Point> rx;       // vector of receivers
  std::vector<double> freqs;   // vector of frequencies
  std::vector<Transmitter> tx; // vector of transmitters

  std::vector<Complex> rsp; // vector of responses

  std::vector<double> rho;                   // vector of resistivities
  std::shared_ptr<SharedTriangulation> mesh; // shared pointer to triangulation

  PetscInt order;                                    // order of finite element space
  std::shared_ptr<dealii::FE_Nedelec<3>> fe_nedelec; // shared pointer to Nedelec finite element
  std::shared_ptr<DoFHandler> dh;                    // shared pointer to degree of freedom handler
  std::shared_ptr<AffineConstaints> constraints;     // shared pointer to affine constraints

  std::shared_ptr<PETScBlockSparseMatrix> K;       // shared pointer to system matrix
  std::shared_ptr<PETScBlockVector> s, e, ghost_e; // shared pointers to solution and right-hand side vectors

  Mat A, B;         // PETSc matrices for preconditioners
  KSP A_ksp, B_ksp; // Krylov subspace solvers

  PetscViewer log; // PETSc viewer

  IS R;                          // Restrict operator
  Mat G, ND_Pi[3];               // discrete gradient and Pi operators
  std::vector<PetscReal> coords; // vertices coordinates

  char iprefix[256], oprefix[256];                      // input and output prefixes
  PetscBool use_ams;                                    // boolean for using ams solver
  PetscReal min_tx_cell_size, min_rx_cell_size, e_rtol; // minimum cell sizes and relative error tolerance

  // integers for maximum iterations, preconditioner type, solver type, threshold, and number of cell
  // refinements
  PetscInt K_max_it, inner_pc_type, direct_solver_type, pc_threshold, n_tx_cell_refinements, n_rx_cell_refinements,
      n_tx_divisions, n_global_refinements;

  PetscClassId EMCTX_ID; // PETSc class ID
  PetscLogEvent SetupSystem, AssembleMat, AssembleRHS, CreatePC, SolveLS, SaveMesh, CalculateRSP,
      SaveRSP; // PETSc log events correspond to each phase of the forward problem
};

PetscErrorCode create_context(EMContext *);
PetscErrorCode destroy_context(EMContext *);
PetscErrorCode process_options(EMContext *);

#endif
