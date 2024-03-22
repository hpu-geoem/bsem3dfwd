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

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>

#include "em_ctx.h"
#include "em_la.h"
#include "em_utils.h"


/**
 * @brief Creates the restriction operator R for the Nedelec elements.
 *
 * This function creates the restriction operator R for the Nedelec elements
 * using the given context, affine constraints for the Nedelec elements, and
 * the locally owned Nedelec dofs. The resulting operator is stored in the
 * provided PETSc index set R.
 *
 * @param ctx Pointer to the EMContext object.
 * @param hnc_nedelec The affine constraints for the Nedelec elements.
 * @param locally_owned_nedelec_dofs The locally owned Nedelec dofs.
 * @param R The PETSc index set to store the resulting operator.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode create_ams_R(EMContext *ctx, const AffineConstaints &hnc_nedelec,
                            const dealii::IndexSet &locally_owned_nedelec_dofs, IS *R) {
  std::vector<PetscInt> local_master_nedelec_dofs;
  PetscInt i, n, beg, end, n_local_nedelec_dofs, n_local_master_nedelec_dofs;

  PetscFunctionBegin;

  n_local_nedelec_dofs = locally_owned_nedelec_dofs.n_elements();

  if (n_local_nedelec_dofs > 0) {
    beg = locally_owned_nedelec_dofs.nth_index_in_set(0);
    end = beg + n_local_nedelec_dofs;
  } else {
    beg = end = -1;
  }

  n_local_master_nedelec_dofs = n_local_nedelec_dofs;
  for (i = beg; i < end; ++i) {
    if (hnc_nedelec.is_constrained(i)) {
      --n_local_master_nedelec_dofs;
    }
  }

  local_master_nedelec_dofs.resize(n_local_master_nedelec_dofs);

  n = 0;
  for (i = beg; i < end; ++i) {
    if (!hnc_nedelec.is_constrained(i)) {
      local_master_nedelec_dofs[n++] = i;
    }
  }
  PetscCall(ISCreateGeneral(ctx->world_comm, n_local_master_nedelec_dofs, &local_master_nedelec_dofs[0], PETSC_COPY_VALUES, R));

  PetscFunctionReturn(0);
}

/**
 * @brief Creates the interpolation matrix P for the nodal elements.
 *
 * This function creates the interpolation matrix P for the nodal elements
 * using the given context, affine constraints for the nodal elements, and
 * the locally owned nodal dofs. The resulting operator is stored in the
 * provided PETSc matrix P.
 *
 * @param ctx Pointer to the EMContext object.
 * @param hnc_nodal The affine constraints for the nodal elements.
 * @param locally_owned_nodal_dofs The locally owned nodal dofs.
 * @param locally_relevant_nodal_dofs The locally relevant nodal dofs.
 * @param tn_to_mn The mapping from total dofs to master dofs.
 * @param P The PETSc matrix to store the resulting operator.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode create_ams_P(EMContext *ctx, const AffineConstaints &hnc_nodal, const dealii::IndexSet &locally_owned_nodal_dofs,
                            const dealii::IndexSet &locally_relevant_nodal_dofs, PETScVector &tn_to_mn, Mat *P) {
  PetscLayout pl;
  std::vector<PetscReal> vals;
  PETScVector total_dofs_to_master_dofs;
  std::vector<PetscInt> row_ptr, col_idx;
  const std::vector<std::pair<AffineConstaints::size_type, double>> *entries;
  PetscInt i, j, n, nnz, r_beg, r_end, c_beg, c_end, n_local_nodal_dofs, n_local_master_nodal_dofs;

  PetscFunctionBegin;

  n_local_nodal_dofs = locally_owned_nodal_dofs.n_elements();

  if (n_local_nodal_dofs > 0) {
    r_beg = locally_owned_nodal_dofs.nth_index_in_set(0);
    r_end = r_beg + n_local_nodal_dofs;
  } else {
    r_beg = r_end = -1;
  }

  n_local_master_nodal_dofs = n_local_nodal_dofs;
  for (i = r_beg; i < r_end; ++i) {
    if (hnc_nodal.is_constrained(i)) {
      --n_local_master_nodal_dofs;
    }
  }

  PetscCall(PetscLayoutCreate(ctx->world_comm, &pl));
  PetscCall(PetscLayoutSetLocalSize(pl, n_local_master_nodal_dofs));
  PetscCall(PetscLayoutSetUp(pl));
  PetscCall(PetscLayoutGetRange(pl, &c_beg, &c_end));
  PetscCall(PetscLayoutDestroy(&pl));

  total_dofs_to_master_dofs.reinit(locally_owned_nodal_dofs, ctx->world_comm);
  tn_to_mn.reinit(locally_owned_nodal_dofs, locally_relevant_nodal_dofs, ctx->world_comm);

  nnz = 0;
  n = c_beg;
  for (i = r_beg; i < r_end; ++i) {
    if (hnc_nodal.is_constrained(i)) {
      nnz += hnc_nodal.get_constraint_entries(i)->size();
      total_dofs_to_master_dofs[(dealii::types::global_dof_index)i] = -1;
    } else {
      total_dofs_to_master_dofs[(dealii::types::global_dof_index)i] = n++;
      nnz += 1;
    }
  }
  total_dofs_to_master_dofs.compress(dealii::VectorOperation::insert);
  tn_to_mn = total_dofs_to_master_dofs;

  row_ptr.resize(n_local_nodal_dofs + 1);
  col_idx.resize(nnz);
  vals.resize(nnz);

  nnz = 0;
  for (i = r_beg; i < r_end; ++i) {
    row_ptr[i - r_beg] = nnz;
    if (hnc_nodal.is_constrained(i)) {
      entries = hnc_nodal.get_constraint_entries(i);
      for (j = 0; j < (PetscInt)entries->size(); ++j) {
        col_idx[nnz] = (PetscInt)tn_to_mn[(*entries)[j].first];
        vals[nnz] = (*entries)[j].second;
        nnz += 1;
      }
    } else {
      col_idx[nnz] = (PetscInt)tn_to_mn[(dealii::types::global_dof_index)i];
      vals[nnz] = 1.0;
      nnz += 1;
    }
  }
  row_ptr[n_local_nodal_dofs] = nnz;

  PetscCall(MatCreateMPIAIJWithArrays(ctx->world_comm, n_local_nodal_dofs, n_local_master_nodal_dofs, PETSC_DECIDE, PETSC_DECIDE,
                                      &row_ptr[0], &col_idx[0], &vals[0], P));

  PetscFunctionReturn(0);
}

/**
 * @brief Creates the discrete gradient matrix G for the lowest order Nedelec element.
 *
 * This function creates the discrete gradient matrix G for the lowest order Nedelec elements
 * using the given context, nodal DoF handler, locally owned Nedelec dofs, and
 * locally owned nodal dofs. The resulting matrix is stored in the provided PETSc matrix G.
 *
 * @param ctx Pointer to the EMContext object.
 * @param dh_nodal The DoF handler for the nodal elements.
 * @param locally_owned_nedelec_dofs The locally owned Nedelec dofs.
 * @param locally_owned_nodal_dofs The locally owned nodal dofs.
 * @param G The PETSc matrix to store the resulting matrix.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode create_ams_G_lowest_order(EMContext *ctx, const DoFHandler &dh_nodal,
                                         const dealii::IndexSet &locally_owned_nedelec_dofs,
                                         const dealii::IndexSet &locally_owned_nodal_dofs, Mat *G) {
  std::vector<PetscReal> vals;
  std::vector<PetscInt> row_ptr, col_idx;
  Triangulation::active_cell_iterator cell;
  DoFHandler::active_cell_iterator cell_nodal, cell_nedelec;
  PetscInt i, n_nedelec_dofs, n_local_nedelec_dofs, n_nodal_dofs, n_local_nodal_dofs, r_beg, r_end, line_index;

  PetscFunctionBegin;

  n_nedelec_dofs = ctx->dh->n_dofs();
  n_local_nedelec_dofs = locally_owned_nedelec_dofs.n_elements();
  n_nodal_dofs = dh_nodal.n_dofs();
  n_local_nodal_dofs = locally_owned_nodal_dofs.n_elements();

  if (n_local_nedelec_dofs > 0) {
    r_beg = locally_owned_nedelec_dofs.nth_index_in_set(0);
    r_end = r_beg + n_local_nedelec_dofs;
  } else {
    r_beg = r_end = -1;
  }

  row_ptr.resize(n_local_nedelec_dofs + 1);
  col_idx.resize(n_local_nedelec_dofs * 2);
  vals.resize(n_local_nedelec_dofs * 2);

  std::fill(row_ptr.begin(), row_ptr.end(), -1);

  for (cell = ctx->mesh->begin_active(); cell != ctx->mesh->end(); ++cell) {
    if (!cell->is_locally_owned()) {
      continue;
    }

    cell_nedelec = DoFHandler::active_cell_iterator(*cell, ctx->dh.get());
    cell_nodal = DoFHandler::active_cell_iterator(*cell, &dh_nodal);

    for (i = 0; i < (PetscInt)dealii::GeometryInfo<3>::lines_per_cell; ++i) {
      line_index = cell_nedelec->line(i)->dof_index(0);

      if (line_index >= r_beg && line_index < r_end) {
        if (row_ptr[line_index - r_beg] < 0) {
          row_ptr[line_index - r_beg] = (line_index - r_beg) * 2;
          col_idx[(line_index - r_beg) * 2 + 0] =
              cell_nodal->vertex_dof_index(dealii::GeometryInfo<3>::line_to_cell_vertices(i, 0), 0);
          col_idx[(line_index - r_beg) * 2 + 1] =
              cell_nodal->vertex_dof_index(dealii::GeometryInfo<3>::line_to_cell_vertices(i, 1), 0);
          vals[(line_index - r_beg) * 2 + 0] = -1;
          vals[(line_index - r_beg) * 2 + 1] = 1;
        }
      }
    }
  }
  row_ptr[n_local_nedelec_dofs] = n_local_nedelec_dofs * 2;

  PetscCall(MatCreateMPIAIJWithArrays(ctx->world_comm, n_local_nedelec_dofs, n_local_nodal_dofs, n_nedelec_dofs, n_nodal_dofs,
                                      &row_ptr[0], &col_idx[0], &vals[0], G));

  PetscFunctionReturn(0);
}

/**
 * @brief Eliminates zero values from a matrix.
 *
 * This function eliminates all values in the given matrix `A` that are below
 * the given tolerance `tol`. The resulting matrix is stored in the provided
 * PETSc matrix `B`.
 *
 * @param A The PETSc matrix to eliminate zero values from.
 * @param tol The tolerance below which values are considered zero.
 * @param B The PETSc matrix to store the resulting matrix.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode eliminate_zero_values(Mat A, PetscReal tol, Mat *B) {
  PetscReal *row_vals;
  std::vector<PetscReal> vals;
  std::vector<PetscInt> row_ptr, col_idx;
  PetscInt i, j, r_beg, r_end, c_beg, c_end, nnz, ncols, *cols;

  PetscFunctionBegin;

  PetscCall(MatGetOwnershipRange(A, &r_beg, &r_end));
  PetscCall(MatGetOwnershipRangeColumn(A, &c_beg, &c_end));

  nnz = 0;
  for (i = r_beg; i < r_end; ++i) {
    PetscCall(MatGetRow(A, i, &ncols, nullptr, (const PetscReal **)&row_vals));
    for (j = 0; j < ncols; ++j) {
      if (std::abs(row_vals[j]) > tol) {
        ++nnz;
      }
    }
    PetscCall(MatRestoreRow(A, i, &ncols, nullptr, (const PetscReal **)&row_vals));
  }

  row_ptr.resize((r_end - r_beg) + 1);
  col_idx.resize(nnz);
  vals.resize(nnz);

  nnz = 0;
  for (i = r_beg; i < r_end; ++i) {
    row_ptr[i - r_beg] = nnz;
    PetscCall(MatGetRow(A, i, &ncols, (const PetscInt **)&cols, (const PetscReal **)&row_vals));
    for (j = 0; j < ncols; ++j) {
      if (std::abs(row_vals[j]) > tol) {
        col_idx[nnz] = cols[j];
        vals[nnz] = row_vals[j];
        ++nnz;
      }
    }
    PetscCall(MatRestoreRow(A, i, &ncols, (const PetscInt **)&cols, (const PetscReal **)&row_vals));
  }
  row_ptr[i - r_beg] = nnz;

  PetscCall(MatCreateMPIAIJWithArrays(PetscObjectComm((PetscObject)A), r_end - r_beg, c_end - c_beg, PETSC_DECIDE, PETSC_DECIDE,
                                      &row_ptr[0], &col_idx[0], &vals[0], B));

  PetscFunctionReturn(0);
}

/**
 * @brief Creates the discrete gradient matrix G and interpolation matrix Pi high-order elements.
 *
 * @param ctx Pointer to the EMContext object.
 * @param dh_nodal The nodal DoF handler.
 * @param locally_relevant_nedelec_dofs The locally relevant Nedelec DoFs.
 * @param ND_G The output matrix G.
 * @param ND_Pi The output matrices Pi_x, Pi_y, and Pi_z.
 *
 * @return PetscErrorCode The error code from the PETSc library.
 */
PetscErrorCode create_ams_G_Pi(EMContext *ctx, const DoFHandler &dh_nodal, const dealii::IndexSet &locally_relevant_nedelec_dofs,
                               Mat *ND_G, Mat ND_Pi[]) {
  size_t i, j, q;
  dealii::DynamicSparsityPattern dsp;
  PETScSparseMatrix G, Pi_x, Pi_y, Pi_z;
  Triangulation::active_cell_iterator cell;
  dealii::IndexSet locally_owned_nedelec_dofs;
  const dealii::FEValuesExtractors::Vector E(0);
  DoFHandler::active_cell_iterator cell_nodal, cell_nedelec;
  dealii::QGauss<3> quadrature_formula(ctx->fe_nedelec->degree + 1);
  std::vector<dealii::types::global_dof_index> nedelec_dof_indices, nodal_dof_indices;
  dealii::FEValues<3> fv_nodal(dh_nodal.get_fe(), quadrature_formula,
                               dealii::update_values | dealii::update_gradients | dealii::update_JxW_values);
  dealii::FEValues<3> fv_nedelec(*ctx->fe_nedelec, quadrature_formula, dealii::update_values | dealii::update_JxW_values);
  const size_t nedelec_dofs_per_cell = ctx->fe_nedelec->dofs_per_cell, nodal_dofs_per_cell = dh_nodal.get_fe().dofs_per_cell;
  dealii::FullMatrix<double> M(nedelec_dofs_per_cell, nedelec_dofs_per_cell), inv_M(nedelec_dofs_per_cell, nedelec_dofs_per_cell),
      D(nedelec_dofs_per_cell, nodal_dofs_per_cell), cell_G(nedelec_dofs_per_cell, nodal_dofs_per_cell),
      cell_Pi_x(nedelec_dofs_per_cell, nodal_dofs_per_cell), cell_Pi_y(nedelec_dofs_per_cell, nodal_dofs_per_cell),
      cell_Pi_z(nedelec_dofs_per_cell, nodal_dofs_per_cell);

  PetscFunctionBegin;

  nedelec_dof_indices.resize(nedelec_dofs_per_cell);
  nodal_dof_indices.resize(nodal_dofs_per_cell);

  dsp.reinit(ctx->dh->n_dofs(), dh_nodal.n_dofs(), locally_relevant_nedelec_dofs);
  for (cell = dh_nodal.get_triangulation().begin_active(); cell != dh_nodal.get_triangulation().end(); ++cell) {
    if (!cell->is_locally_owned()) {
      continue;
    }

    cell_nedelec = DoFHandler::active_cell_iterator(*cell, ctx->dh.get());
    cell_nodal = DoFHandler::active_cell_iterator(*cell, &dh_nodal);

    cell_nedelec->get_dof_indices(nedelec_dof_indices);
    cell_nodal->get_dof_indices(nodal_dof_indices);

    for (i = 0; i < nedelec_dofs_per_cell; ++i) {
      for (j = 0; j < nodal_dofs_per_cell; ++j) {
        dsp.add(nedelec_dof_indices[i], nodal_dof_indices[j]);
      }
    }
  }
  locally_owned_nedelec_dofs = ctx->dh->locally_owned_dofs();
  dealii::SparsityTools::distribute_sparsity_pattern(dsp, locally_owned_nedelec_dofs, ctx->world_comm,
                                                     locally_relevant_nedelec_dofs);
  dsp.compress();

  G.reinit(ctx->dh->locally_owned_dofs(), dh_nodal.locally_owned_dofs(), dsp, ctx->world_comm);
  Pi_x.reinit(ctx->dh->locally_owned_dofs(), dh_nodal.locally_owned_dofs(), dsp, ctx->world_comm);
  Pi_y.reinit(ctx->dh->locally_owned_dofs(), dh_nodal.locally_owned_dofs(), dsp, ctx->world_comm);
  Pi_z.reinit(ctx->dh->locally_owned_dofs(), dh_nodal.locally_owned_dofs(), dsp, ctx->world_comm);

  for (cell = dh_nodal.get_triangulation().begin_active(); cell != dh_nodal.get_triangulation().end(); ++cell) {
    if (!cell->is_locally_owned()) {
      continue;
    }

    cell_nedelec = DoFHandler::active_cell_iterator(*cell, ctx->dh.get());
    cell_nodal = DoFHandler::active_cell_iterator(*cell, &dh_nodal);

    cell_nedelec->get_dof_indices(nedelec_dof_indices);
    cell_nodal->get_dof_indices(nodal_dof_indices);

    fv_nedelec.reinit(cell_nedelec);
    fv_nodal.reinit(cell_nodal);

    M = 0.0;
    for (i = 0; i < nedelec_dofs_per_cell; ++i) {
      for (j = 0; j < nedelec_dofs_per_cell; ++j) {
        for (q = 0; q < quadrature_formula.size(); ++q) {
          M(i, j) += fv_nedelec[E].value(i, q) * fv_nedelec[E].value(j, q) * fv_nedelec.JxW(q);
        }
      }
    }
    inv_M.invert(M);

    cell_G = 0.0;
    cell_Pi_x = 0.0;
    cell_Pi_y = 0.0;
    cell_Pi_z = 0.0;
    for (i = 0; i < nedelec_dofs_per_cell; ++i) {
      for (j = 0; j < nodal_dofs_per_cell; ++j) {
        for (q = 0; q < quadrature_formula.size(); ++q) {
          cell_G(i, j) += fv_nedelec[E].value(i, q) * fv_nodal.shape_grad(j, q) * fv_nedelec.JxW(q);
          cell_Pi_x(i, j) += fv_nedelec[E].value(i, q)[0] * fv_nodal.shape_value(j, q) * fv_nedelec.JxW(q);
          cell_Pi_y(i, j) += fv_nedelec[E].value(i, q)[1] * fv_nodal.shape_value(j, q) * fv_nedelec.JxW(q);
          cell_Pi_z(i, j) += fv_nedelec[E].value(i, q)[2] * fv_nodal.shape_value(j, q) * fv_nedelec.JxW(q);
        }
      }
    }

    inv_M.mmult(D, cell_G);
    G.set(nedelec_dof_indices, nodal_dof_indices, D);

    inv_M.mmult(D, cell_Pi_x);
    Pi_x.set(nedelec_dof_indices, nodal_dof_indices, D);

    inv_M.mmult(D, cell_Pi_y);
    Pi_y.set(nedelec_dof_indices, nodal_dof_indices, D);

    inv_M.mmult(D, cell_Pi_z);
    Pi_z.set(nedelec_dof_indices, nodal_dof_indices, D);
  }

  G.compress(dealii::VectorOperation::insert);
  Pi_x.compress(dealii::VectorOperation::insert);
  Pi_y.compress(dealii::VectorOperation::insert);
  Pi_z.compress(dealii::VectorOperation::insert);

  PetscCall(eliminate_zero_values(G, EPS, ND_G));
  PetscCall(eliminate_zero_values(Pi_x, EPS, &ND_Pi[0]));
  PetscCall(eliminate_zero_values(Pi_y, EPS, &ND_Pi[1]));
  PetscCall(eliminate_zero_values(Pi_z, EPS, &ND_Pi[2]));

  PetscFunctionReturn(0);
}

/**
 * @brief Creates the vertex coordinates for the AMS context.
 *
 * This function creates the vertex coordinates for the AMS context using the given
 * nodal DoF handler, total-to-master DoF mapping, and interpolation matrix P.
 * The resulting vertex coordinates are stored in the EMContext object.
 *
 * @param ctx Pointer to the EMContext object.
 * @param dh_nodal The DoF handler for the nodal elements.
 * @param tn_to_mn The mapping from total dofs to master dofs.
 * @param P The interpolation matrix for the nodal elements.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode create_ams_vertex_coords(EMContext *ctx, const DoFHandler &dh_nodal, const PETScVector &tn_to_mn, Mat P) {
  Point vertex;
  PetscInt i, vidx, c_beg, c_end;
  DoFHandler::active_cell_iterator cell;

  PetscFunctionBegin;

  PetscCall(MatGetOwnershipRangeColumn(P, &c_beg, &c_end));

  ctx->coords.resize((c_end - c_beg) * 3);

  for (cell = dh_nodal.begin_active(); cell != dh_nodal.end(); ++cell) {
    if (!cell->is_locally_owned()) {
      continue;
    }

    for (i = 0; i < (PetscInt)dealii::GeometryInfo<3>::lines_per_cell; ++i) {
      vidx = (PetscInt)tn_to_mn[cell->vertex_dof_index(dealii::GeometryInfo<3>::line_to_cell_vertices(i, 0), 0)];
      if (vidx != -1 && vidx >= c_beg && vidx < c_end) {
        vertex = cell->vertex(dealii::GeometryInfo<3>::line_to_cell_vertices(i, 0));
        ctx->coords[(vidx - c_beg) * 3 + 0] = vertex[0];
        ctx->coords[(vidx - c_beg) * 3 + 1] = vertex[1];
        ctx->coords[(vidx - c_beg) * 3 + 2] = vertex[2];
      }

      vidx = (PetscInt)tn_to_mn[cell->vertex_dof_index(dealii::GeometryInfo<3>::line_to_cell_vertices(i, 1), 0)];
      if (vidx != -1 && vidx >= c_beg && vidx < c_end) {
        vertex = cell->vertex(dealii::GeometryInfo<3>::line_to_cell_vertices(i, 1));
        ctx->coords[(vidx - c_beg) * 3 + 0] = vertex[0];
        ctx->coords[(vidx - c_beg) * 3 + 1] = vertex[1];
        ctx->coords[(vidx - c_beg) * 3 + 2] = vertex[2];
      }
    }
  }

  PetscFunctionReturn(0);
}

/**
 * @brief Sets up the AMS preconditioner for the electromagnetic problem.
 *
 * This function sets up the AMS for the electromagnetic problem by creating the necessary
 * matrices and vectors, and initializing the Krylov subspace solver for the preconditioner.
 *
 * @param ctx Pointers to the EMContext object.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode setup_ams(EMContext *ctx) {
  PetscInt i;
  DoFHandler dh_nodal;
  PETScVector tn_to_mn;
  Mat P, G, Pi[3], RX, H;
  AffineConstaints hnc_nedelec, hnc_nodal;
  dealii::FE_Q<3> fe_nodal(ctx->fe_nedelec->degree);
  dealii::IndexSet locally_owned_nedelec_dofs, locally_relevant_nedelec_dofs, locally_owned_nodal_dofs,
      locally_relevant_nodal_dofs;

  PetscFunctionBegin;

  P = G = Pi[0] = Pi[1] = Pi[2] = RX = H = nullptr;

  locally_owned_nedelec_dofs = ctx->dh->locally_owned_dofs();
  dealii::DoFTools::extract_locally_relevant_dofs(*ctx->dh, locally_relevant_nedelec_dofs);

  hnc_nedelec.reinit(locally_relevant_nedelec_dofs);
  dealii::DoFTools::make_hanging_node_constraints(*ctx->dh, hnc_nedelec);
  hnc_nedelec.close();

  dh_nodal.reinit(ctx->dh->get_triangulation());
  dh_nodal.distribute_dofs(fe_nodal);

  locally_owned_nodal_dofs = dh_nodal.locally_owned_dofs();
  dealii::DoFTools::extract_locally_relevant_dofs(dh_nodal, locally_relevant_nodal_dofs);

  hnc_nodal.reinit(locally_relevant_nodal_dofs);
  dealii::DoFTools::make_hanging_node_constraints(dh_nodal, hnc_nodal);
  hnc_nodal.close();

  PetscCall(create_ams_R(ctx, hnc_nedelec, locally_owned_nedelec_dofs, &ctx->R));
  PetscCall(create_ams_P(ctx, hnc_nodal, locally_owned_nodal_dofs, locally_relevant_nodal_dofs, tn_to_mn, &P));

  if (ctx->order == 0) {
    PetscCall(create_ams_vertex_coords(ctx, dh_nodal, tn_to_mn, P));
    PetscCall(create_ams_G_lowest_order(ctx, dh_nodal, locally_owned_nedelec_dofs, locally_owned_nodal_dofs, &G));
  } else {
    PetscCall(create_ams_G_Pi(ctx, dh_nodal, locally_relevant_nedelec_dofs, &G, Pi));
  }

  PetscCall(MatCreateSubMatrix(G, ctx->R, nullptr, MAT_INITIAL_MATRIX, &RX));
  PetscCall(MatMatMult(RX, P, MAT_INITIAL_MATRIX, 1.0, &H));
  PetscCall(MatDestroy(&RX));

  PetscCall(MatConvert(H, MATHYPRE, MAT_INITIAL_MATRIX, &ctx->G));
  PetscCall(MatDestroy(&H));

  if (ctx->order > 0) {
    for (i = 0; i < 3; ++i) {
      PetscCall(MatCreateSubMatrix(Pi[i], ctx->R, nullptr, MAT_INITIAL_MATRIX, &RX));
      PetscCall(MatMatMult(RX, P, MAT_INITIAL_MATRIX, 1.0, &H));
      PetscCall(MatDestroy(&RX));
      PetscCall(MatConvert(H, MATHYPRE, MAT_INITIAL_MATRIX, &ctx->ND_Pi[i]));
      PetscCall(MatDestroy(&H));
    }
  }

  PetscCall(MatDestroy(&G));
  PetscCall(MatDestroy(&P));
  PetscCall(MatDestroy(&Pi[0]));
  PetscCall(MatDestroy(&Pi[1]));
  PetscCall(MatDestroy(&Pi[2]));
  PetscCall(MatDestroy(&RX));

  dh_nodal.clear();

  PetscFunctionReturn(0);
}

/**
 * @brief Applies the inner preconditioner.
 *
 * @param pc The preconditioner to apply.
 * @param b The input vector to apply the preconditioner to.
 * @param x The output vector to store the result of the preconditioner.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode pc_apply_a(PC pc, Vec b, Vec x) {
  EMContext *ctx;
  Vec br, bi, xr, xi, bb, xx;

  PetscFunctionBegin;

  PetscCall(PCShellGetContext(pc, (void **)&ctx));
  PetscCall(PetscViewerASCIIPushTab(ctx->log));

  PetscCall(VecNestGetSubVec(b, 0, &br));
  PetscCall(VecNestGetSubVec(b, 1, &bi));
  PetscCall(VecNestGetSubVec(x, 0, &xr));
  PetscCall(VecNestGetSubVec(x, 1, &xi));

  if (ctx->use_ams) {
    PetscCall(VecGetSubVector(br, ctx->R, &bb));
    PetscCall(VecGetSubVector(xr, ctx->R, &xx));
    PetscCall(KSPSolve(ctx->B_ksp, bb, xx));
    PetscCall(VecRestoreSubVector(xr, ctx->R, &xx));
    PetscCall(VecRestoreSubVector(br, ctx->R, &bb));

    PetscCall(VecGetSubVector(bi, ctx->R, &bb));
    PetscCall(VecGetSubVector(xi, ctx->R, &xx));
    PetscCall(KSPSolve(ctx->B_ksp, bb, xx));
    PetscCall(VecRestoreSubVector(xi, ctx->R, &xx));
    PetscCall(VecRestoreSubVector(bi, ctx->R, &bb));
  } else {
    PetscCall(KSPSolve(ctx->B_ksp, br, xr));
    PetscCall(KSPSolve(ctx->B_ksp, bi, xi));
  }

  PetscCall(PetscViewerASCIIPopTab(ctx->log));

  PetscFunctionReturn(0);
}

/**
 * @brief Creates the block-diagonal preconditioner for the global lienar system.
 *
 * @param ctx Pointer to the EMContext object.
 * @return PetscErrorCode indicating success or failure of the function.
 */
PetscErrorCode create_preconditioner(EMContext *ctx) {
  PC A_pc, B_pc;
  Mat A_blocks[4], B;
  PetscViewerAndFormat *vf;

  PetscFunctionBegin;

  LogEventHelper leh(ctx->CreatePC);

  PetscCall(MatDuplicate(ctx->K->block(0, 0), MAT_COPY_VALUES, &B));
  PetscCall(MatAXPY(B, -1.0, ctx->K->block(0, 1), DIFFERENT_NONZERO_PATTERN));

  if (ctx->use_ams) {
    PetscCall(setup_ams(ctx));
    PetscCall(MatCreateSubMatrix(B, ctx->R, ctx->R, MAT_INITIAL_MATRIX, &ctx->B));
  } else {
    ctx->B = B;
    PetscCall(PetscObjectReference((PetscObject)ctx->B));
  }
  PetscCall(MatDestroy(&B));
  PetscCall(PetscObjectSetName((PetscObject)(ctx->B), "B"));
  PetscCall(MatSetOptionsPrefix(ctx->B, "B_"));
  PetscCall(MatSetOption(ctx->B, MAT_SPD, PETSC_TRUE));
  PetscCall(MatSetFromOptions(ctx->B));

  PetscCall(KSPCreate(ctx->world_comm, &ctx->B_ksp));
  PetscCall(KSPSetOptionsPrefix(ctx->B_ksp, "B_"));
  PetscCall(KSPSetOperators(ctx->B_ksp, ctx->B, ctx->B));

  PetscCall(KSPGetPC(ctx->B_ksp, &B_pc));
  if (ctx->use_ams) {
    PetscCall(KSPSetType(ctx->B_ksp, KSPCG));
    PetscCall(KSPSetNormType(ctx->B_ksp, KSP_NORM_UNPRECONDITIONED));
    PetscCall(KSPSetTolerances(ctx->B_ksp, 1E-2, PETSC_DEFAULT, PETSC_DEFAULT, 100));
    PetscCall(PCSetType(B_pc, PCHYPRE));
    PetscCall(PCHYPRESetType(B_pc, "ams"));
    PetscCall(PCHYPRESetDiscreteGradient(B_pc, ctx->G));
    if (ctx->order == 0) {
      PetscCall(PCSetCoordinates(B_pc, 3, (PetscInt)ctx->coords.size() / 3, &ctx->coords[0]));
    } else {
      PetscCall(PCHYPRESetInterpolations(B_pc, 3, nullptr, nullptr, nullptr, ctx->ND_Pi));
    }
  } else {
#ifndef PETSC_HAVE_SUPERLU_DIST
    if (ctx->direct_solver_type == SUPERLU_DIST) {
      ctx->direct_solver_type = MUMPS;
    }
#endif

#ifndef PETSC_HAVE_MKL_CPARDISO
    if (ctx->direct_solver_type == CPARDISO) {
      ctx->direct_solver_type = MUMPS;
    }
#endif

    PetscCall(KSPSetType(ctx->B_ksp, KSPPREONLY));
    PetscCall(PCSetType(B_pc, PCCHOLESKY));
    if (ctx->direct_solver_type == SUPERLU_DIST) {
      PetscCall(PCFactorSetMatSolverType(B_pc, MATSOLVERSUPERLU_DIST));
    } else if (ctx->direct_solver_type == CPARDISO) {
      PetscCall(PCFactorSetMatSolverType(B_pc, MATSOLVERMKL_CPARDISO));
    } else {
      PetscCall(PCFactorSetMatSolverType(B_pc, MATSOLVERMUMPS));
    }
  }
  PetscCall(KSPSetFromOptions(ctx->B_ksp));
  PetscCall(KSPSetUp(ctx->B_ksp));

  if (ctx->use_ams) {
    PetscCall(PetscViewerAndFormatCreate(ctx->log, PETSC_VIEWER_DEFAULT, &vf));
#if PETSC_VERSION_LT(3, 15, 0)
    PetscCall(KSPMonitorSet(ctx->B_ksp, (PetscErrorCode(*)(KSP, PetscInt, PetscReal, void *))KSPMonitorTrueResidualNorm, vf,
                            (PetscErrorCode(*)(void **))PetscViewerAndFormatDestroy));
#else
    PetscCall(KSPMonitorSet(ctx->B_ksp, (PetscErrorCode(*)(KSP, PetscInt, PetscReal, void *))KSPMonitorTrueResidual, vf,
                            (PetscErrorCode(*)(void **))PetscViewerAndFormatDestroy));
#endif
  }

  A_blocks[0] = ctx->K->block(0, 0);
  A_blocks[1] = ctx->K->block(0, 1);
  A_blocks[2] = ctx->K->block(1, 0);
  A_blocks[3] = ctx->K->block(1, 1);
  PetscCall(MatCreateNest(ctx->world_comm, 2, nullptr, 2, nullptr, &A_blocks[0], &ctx->A));
  PetscCall(MatNestSetVecType(ctx->A, VECNEST));
  PetscCall(PetscObjectSetName((PetscObject)(ctx->A), "A"));
  PetscCall(MatSetOptionsPrefix(ctx->A, "A_"));
  PetscCall(MatSetFromOptions(ctx->A));

  PetscCall(KSPCreate(ctx->world_comm, &ctx->A_ksp));
  PetscCall(KSPSetOptionsPrefix(ctx->A_ksp, "A_"));
  PetscCall(KSPSetOperators(ctx->A_ksp, ctx->A, ctx->A));
  PetscCall(KSPSetType(ctx->A_ksp, KSPFGMRES));
  PetscCall(KSPGetPC(ctx->A_ksp, &A_pc));
  PetscCall(PCSetType(A_pc, PCSHELL));
  PetscCall(PCShellSetContext(A_pc, ctx));
  PetscCall(PCShellSetApply(A_pc, pc_apply_a));
  PetscCall(KSPSetFromOptions(ctx->A_ksp));
  PetscCall(KSPSetUp(ctx->A_ksp));

  PetscCall(PetscViewerAndFormatCreate(ctx->log, PETSC_VIEWER_DEFAULT, &vf));
#if PETSC_VERSION_LT(3, 15, 0)
  PetscCall(KSPMonitorSet(ctx->A_ksp, (PetscErrorCode(*)(KSP, PetscInt, PetscReal, void *))KSPMonitorTrueResidualNorm, vf,
                          (PetscErrorCode(*)(void **))PetscViewerAndFormatDestroy));
#else
  PetscCall(KSPMonitorSet(ctx->A_ksp, (PetscErrorCode(*)(KSP, PetscInt, PetscReal, void *))KSPMonitorTrueResidual, vf,
                          (PetscErrorCode(*)(void **))PetscViewerAndFormatDestroy));
#endif

  PetscFunctionReturn(0);
}

/**
 * @brief Destroys the preconditioner for the electromagnetic problem.
 *
 * This function destroys the preconditioner for the electromagnetic problem by destroying the necessary
 * matrices and vectors, and freeing the memory allocated for the Krylov subspace solver.
 *
 * @param ctx Pointer to the EMContext object.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode destroy_preconditioner(EMContext *ctx) {
  PetscFunctionBegin;

  PetscCall(KSPDestroy(&ctx->A_ksp));
  PetscCall(KSPDestroy(&ctx->B_ksp));
  PetscCall(MatDestroy(&ctx->A));
  PetscCall(MatDestroy(&ctx->B));

  PetscCall(MatDestroy(&ctx->G));
  PetscCall(MatDestroy(&ctx->ND_Pi[0]));
  PetscCall(MatDestroy(&ctx->ND_Pi[1]));
  PetscCall(MatDestroy(&ctx->ND_Pi[2]));
  PetscCall(ISDestroy(&ctx->R));

  PetscFunctionReturn(0);
}

/**
 * @brief Solves the linear system Ax = b using the Krylov subspace method.
 *
 * This function solves the linear system Ax = b using the Krylov subspace method.
 * The solution is stored in the PETScBlockVector e.
 *
 * @param ctx Pointer to the EMContext object.
 * @param s The right-hand side vector.
 * @param e The solution vector.
 * @param max_it The maximum number of iterations.
 * @param rtol The relative tolerance.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode solve_linear_system(EMContext *ctx, const PETScBlockVector &s, PETScBlockVector &e, PetscInt max_it,
                                   PetscReal rtol) {
  Vec xx[2], bb[2], x, b;

  PetscFunctionBegin;

  LogEventHelper leh(ctx->SolveLS);

  e = 0.0;

  xx[0] = e.block(0);
  xx[1] = e.block(1);
  PetscCall(VecCreateNest(ctx->world_comm, 2, nullptr, xx, &x));
  bb[0] = s.block(0);
  bb[1] = s.block(1);
  PetscCall(VecCreateNest(ctx->world_comm, 2, nullptr, bb, &b));

  PetscCall(KSPSetTolerances(ctx->A_ksp, rtol, PETSC_DEFAULT, PETSC_DEFAULT, max_it));
  PetscCall(KSPSolve(ctx->A_ksp, b, x));

  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&b));

  PetscFunctionReturn(0);
}
