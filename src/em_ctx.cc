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

#include <petsclog.h>

#include "em_ctx.h"
#include "em_utils.h"

/**
 * @brief Creates an EMContext object and initializes its members.
 *
 * @param ctx Pointer to the EMContext object to be created.
 * @return PetscErrorCode 0 on success, non-zero on failure.
 *
 * This function creates an EMContext object and initializes its members. It
 * sets up the MPI communicator, the triangulation, the finite element space,
 * and the PETSc matrices and vectors. It also registers PETSc log events and
 * opens a log file for output.
 */
PetscErrorCode create_context(EMContext *ctx) {
  PetscFunctionBegin;

  PetscCall(MPI_Comm_dup(MPI_COMM_WORLD, &ctx->world_comm));
  ctx->world_size = dealii::Utilities::MPI::n_mpi_processes(ctx->world_comm);
  ctx->world_rank = dealii::Utilities::MPI::this_mpi_process(ctx->world_comm);

  ctx->mesh.reset(new SharedTriangulation(ctx->world_comm, SharedTriangulation::none, false));

  ctx->fe_nedelec.reset(new dealii::FE_Nedelec<3>(ctx->order));
  ctx->dh.reset(new DoFHandler);
  ctx->constraints.reset(new AffineConstaints);
  ctx->K.reset(new PETScBlockSparseMatrix);
  ctx->s.reset(new PETScBlockVector);
  ctx->e.reset(new PETScBlockVector);
  ctx->ghost_e.reset(new PETScBlockVector);

  ctx->A = nullptr;
  ctx->B = nullptr;
  ctx->A_ksp = nullptr;
  ctx->B_ksp = nullptr;

  ctx->R = nullptr;
  ctx->G = nullptr;
  ctx->ND_Pi[0] = nullptr;
  ctx->ND_Pi[1] = nullptr;
  ctx->ND_Pi[2] = nullptr;

  PetscCall(PetscClassIdRegister("EMCTX", &ctx->EMCTX_ID));
  PetscCall(PetscLogEventRegister("SetupSystem", ctx->EMCTX_ID, &ctx->SetupSystem));
  PetscCall(PetscLogEventRegister("AssembleMat", ctx->EMCTX_ID, &ctx->AssembleMat));
  PetscCall(PetscLogEventRegister("AssembleRHS", ctx->EMCTX_ID, &ctx->AssembleRHS));
  PetscCall(PetscLogEventRegister("CreatePC", ctx->EMCTX_ID, &ctx->CreatePC));
  PetscCall(PetscLogEventRegister("SolveLS", ctx->EMCTX_ID, &ctx->SolveLS));
  PetscCall(PetscLogEventRegister("SaveMesh", ctx->EMCTX_ID, &ctx->SaveMesh));
  PetscCall(PetscLogEventRegister("CalculateRSP", ctx->EMCTX_ID, &ctx->CalculateRSP));
  PetscCall(PetscLogEventRegister("SaveRSP", ctx->EMCTX_ID, &ctx->SaveRSP));
  PetscCall(PetscLogDefaultBegin());

  PetscCall(PetscViewerASCIIOpen(ctx->world_comm, string_format("%s.log", ctx->oprefix).c_str(), &ctx->log));

  PetscFunctionReturn(0);
}

/**
 * @brief Destroys an EMContext object and frees its memory.
 *
 * @param ctx Pointer to the EMContext object to be destroyed.
 * @return PetscErrorCode 0 on success, non-zero on failure.
 *
 * This function destroys an EMContext object and frees its memory. It sets all
 * the PETSc matrices and vectors to nullptr, and frees the triangulation, the
 * finite element space, and the constraints. It also destroys the PETSc log
 * viewer and frees the MPI communicator.
 */
PetscErrorCode destroy_context(EMContext *ctx) {
  PetscFunctionBegin;

  ctx->K = nullptr;
  ctx->s = nullptr;
  ctx->e = nullptr;
  ctx->ghost_e = nullptr;

  ctx->constraints = nullptr;
  ctx->fe_nedelec = nullptr;

  ctx->mesh = nullptr;

  PetscCall(PetscLogView(ctx->log));
  PetscCall(PetscViewerDestroy(&ctx->log));

  PetscCall(MPI_Comm_free(&ctx->world_comm));

  PetscFunctionReturn(0);
}

/**
 * @brief Processes the command line options and sets the corresponding
 * EMContext members.
 *
 * @param ctx Pointer to the EMContext object.
 * @return PetscErrorCode 0 on success, non-zero on failure.
 *
 * This function processes the command line options and sets the corresponding
 * EMContext members. It uses the PETScOptions API to retrieve the values of the
 * options and sets the corresponding members of the EMContext object. The
 * options that are processed include the input and output file prefixes, the
 * type of inner preconditioner, the type of direct solver, the preconditioner
 * threshold, the order of the finite element space, the number of refinements
 * for the transmitter and receiver cells, the number of divisions for the
 * transmitter cells, the minimum size of the transmitter and receiver cells,
 * the maximum number of iterations for the Krylov solver, and the relative
 * tolerance for the electric field solver.
 *
 * The function returns 0 on success and a non-zero error code on failure.
 */
PetscErrorCode process_options(EMContext *ctx) {
  PetscBool flg;
  const char *InnerPCType[] = { "mixed", "ams", "direct" };
  const char *DirectSolverType[] = { "mumps", "superlu_dist", "cpardiso" };

  PetscFunctionBegin;

  PetscOptionsBegin(MPI_COMM_WORLD, "", "BSEM3DFWD options", "");

  flg = PETSC_FALSE;
  PetscCall(PetscOptionsGetString(nullptr, nullptr, "-iprefix", ctx->iprefix, sizeof(ctx->iprefix), &flg));
  if (!flg) {
    SETERRQ(ctx->world_comm, EM_ERR_USER, "Plase specify the prefix of the input file.");
  }

  flg = PETSC_FALSE;
  PetscCall(PetscOptionsGetString(nullptr, nullptr, "-oprefix", ctx->oprefix, sizeof(ctx->oprefix), &flg));
  if (!flg) {
    SETERRQ(ctx->world_comm, EM_ERR_USER, "Plase specify the prefix of the output file.");
  }

  ctx->use_ams = PETSC_TRUE;

  ctx->inner_pc_type = Mixed;
  PetscCall(PetscOptionsEList("-inner_pc_type", "", "", InnerPCType, sizeof(InnerPCType) / sizeof(InnerPCType[0]),
                              InnerPCType[ctx->inner_pc_type], &ctx->inner_pc_type, &flg));

  ctx->direct_solver_type = MUMPS;
  PetscCall(PetscOptionsEList("-direct_solver_type", "", "", DirectSolverType,
                              sizeof(DirectSolverType) / sizeof(DirectSolverType[0]), DirectSolverType[ctx->direct_solver_type],
                              &ctx->direct_solver_type, &flg));

  ctx->pc_threshold = 500000;
  PetscCall(PetscOptionsGetInt(nullptr, nullptr, "-pc_threshold", &ctx->pc_threshold, &flg));

  ctx->order = 0;
  PetscCall(PetscOptionsGetInt(nullptr, nullptr, "-order", &ctx->order, &flg));

  ctx->n_tx_cell_refinements = 0;
  PetscCall(PetscOptionsGetInt(nullptr, nullptr, "-n_tx_cell_refinements", &ctx->n_tx_cell_refinements, &flg));

  ctx->n_rx_cell_refinements = 0;
  PetscCall(PetscOptionsGetInt(nullptr, nullptr, "-n_rx_cell_refinements", &ctx->n_rx_cell_refinements, &flg));

  ctx->n_global_refinements = 0;
  PetscCall(PetscOptionsGetInt(nullptr, nullptr, "-n_global_refinements", &ctx->n_global_refinements, &flg));

  ctx->n_tx_divisions = 20;
  PetscCall(PetscOptionsGetInt(nullptr, nullptr, "-n_tx_divisions", &ctx->n_tx_divisions, &flg));

  ctx->min_tx_cell_size = -1.0;
  PetscCall(PetscOptionsGetReal(nullptr, nullptr, "-min_tx_cell_size", &ctx->min_tx_cell_size, &flg));

  ctx->min_rx_cell_size = -1.0;
  PetscCall(PetscOptionsGetReal(nullptr, nullptr, "-min_rx_cell_size", &ctx->min_rx_cell_size, &flg));

  ctx->K_max_it = 100;
  PetscCall(PetscOptionsGetInt(nullptr, nullptr, "-K_max_it", &ctx->K_max_it, &flg));

  ctx->e_rtol = 1.0E-8;
  PetscCall(PetscOptionsGetReal(nullptr, nullptr, "-e_rtol", &ctx->e_rtol, &flg));

  PetscOptionsEnd();

  PetscFunctionReturn(0);
}
