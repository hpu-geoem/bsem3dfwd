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

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include "em_ctx.h"
#include "em_io.h"
#include "em_la.h"
#include "em_response.h"
#include "em_utils.h"

/**
 * @brief Divides a line current into multiple dipoles.
 *
 * Given a transmitter `tx` and the number of divisions `n_tx_divisions`, this
 * function divides the line current into multiple dipoles and returns a vector
 * of tuples containing the midpoint of each dipole, its direction, current, and
 * weight.
 *
 * @param tx The transmitter to divide.
 * @param n_tx_divisions The number of divisions to make.
 * @return A vector of tuples containing the midpoint of each dipole, its
 * direction, current, and weight.
 */
std::vector<std::tuple<Point, VectorD, double, double>> divide_line_current(const Transmitter &tx, int n_tx_divisions) {
  int i;
  Point A, B;
  VectorD direction;
  std::vector<std::tuple<Point, VectorD, double, double>> dipoles;

  direction[0] = std::cos(tx.azimuth * DTOR) * std::cos(tx.dip * DTOR);
  direction[1] = std::sin(tx.azimuth * DTOR) * std::cos(tx.dip * DTOR);
  direction[2] = std::sin(tx.dip * DTOR);

  if (std::abs(tx.length) <= 0.1) {
    dipoles.push_back(std::make_tuple(tx.center, direction, tx.current, 1.0));
  } else {
    A = tx.center - direction * tx.length / 2;
    for (i = 0; i < n_tx_divisions; ++i) {
      B = A + direction * tx.length / n_tx_divisions;
      dipoles.push_back(std::make_tuple((A + B) / 2, direction, tx.current, 1.0 / n_tx_divisions));
      A = B;
    }
  }

  return dipoles;
}

/**
 * @brief Returns the locally owned degrees of freedom (DOFs) of a given
 * DoFHandler object.
 *
 * This function returns the locally owned DOFs of a given DoFHandler object as
 * an IndexSet. It then creates a new IndexSet that contains both the original
 * DOFs and the DOFs shifted by the number of total DOFs in the system. This is
 * useful for block matrix assembly, where the same DoFHandler is used for both the
 * rea and imaginary parts.
 *
 * @param dh The DoFHandler object to extract the locally owned DOFs from.
 * @return An IndexSet containing the locally owned DOFs and their shifted
 * counterparts.
 */
dealii::IndexSet get_locally_owned_dofs_block(const DoFHandler &dh) {
  dealii::IndexSet dof_set, block_dof_set;

  dof_set = dh.locally_owned_dofs();

  block_dof_set = dealii::IndexSet(dh.n_dofs() * 2);
  block_dof_set.add_indices(dof_set);
  block_dof_set.add_indices(dof_set, dh.n_dofs());
  block_dof_set.compress();

  return block_dof_set;
}

/**
 * @brief Returns the locally relevant degrees of freedom (DOFs) of a given
 * DoFHandler object.
 *
 * This function returns the locally relevant DOFs of a given DoFHandler object as
 * an IndexSet. It then creates a new IndexSet that contains both the original
 * DOFs and the DOFs shifted by the number of total DOFs in the system. This is
 * useful for block matrix assembly, where the same DoFHandler is used for both the
 * real and imaginary parts.
 *
 * @param dh The DoFHandler object to extract the locally relevant DOFs from.
 * @return An IndexSet containing the locally relevant DOFs and their shifted
 * counterparts.
 */
dealii::IndexSet get_locally_relevant_dofs_block(const DoFHandler &dh) {
  dealii::IndexSet dof_set, block_dof_set;

  dealii::DoFTools::extract_locally_relevant_dofs(dh, dof_set);

  block_dof_set = dealii::IndexSet(dh.n_dofs() * 2);
  block_dof_set.add_indices(dof_set);
  block_dof_set.add_indices(dof_set, dh.n_dofs());

  return block_dof_set;
}

/**
 * @brief Finds the active cell around a given point.
 *
 * This function finds the active cell around a given point by first using the
 * `find_active_cell_around_point` function from the `GridTools` namespace to
 * get an initial guess. It then uses `find_all_active_cells_around_point` to
 * find all active cells around the point within a given tolerance. Finally, it
 * returns the active cell with the smallest diameter.
 *
 * @param cache The cache object to use for the search.
 * @param p The point to search around.
 * @return A pair containing an iterator to the active cell and the point within
 * the cell closest to the input point.
 */
std::pair<Triangulation::active_cell_iterator, Point> find_active_cell_around_point(const dealii::GridTools::Cache<3> &cache,
                                                                                    const Point &p) {
  auto cell_point = dealii::GridTools::find_active_cell_around_point(cache, p);
  auto cells =
      dealii::GridTools::find_all_active_cells_around_point(cache.get_mapping(), cache.get_triangulation(), p, 1E-10, cell_point);

  int min_idx = -1;
  double min_size = 1E10;
  for (unsigned int i = 0; i < cells.size(); ++i) {
    if (cells[i].first->diameter() < min_size) {
      min_idx = i;
      min_size = cells[i].first->diameter();
    }
  }

  return cells[min_idx];
}

/**
 * @brief Refines the mesh around the transmitter and receiver areas.
 *
 * This function refines the mesh around the transmitter and receiver areas by
 * finding the active cell around each point and setting its refine flag if its
 * diameter is greater than the minimum cell size. The mesh is then refined
 * according to the refine flags.
 *
 * @param ctx Pointer to the EMContext object.
 * @return PetscErrorCode indicating success or failure of the function.
 */
PetscErrorCode refine_tx_rx_area(EMContext *ctx) {
  PetscFunctionBegin;

  dealii::GridTools::Cache<3> cache(*ctx->mesh);

  for (int r = 0; r < ctx->n_rx_cell_refinements; ++r) {
    for (auto p : ctx->rx) {
      auto cell = find_active_cell_around_point(cache, p).first;

      if (cell->diameter() > ctx->min_rx_cell_size) {
        cell->set_refine_flag();
      }
    }
    ctx->mesh->execute_coarsening_and_refinement();
  }

  std::vector<Point> tx_centers;
  for (auto tx : ctx->tx) {
    auto dipoles = divide_line_current(tx, ctx->n_tx_divisions);
    for (auto d : dipoles) {
      tx_centers.push_back(std::get<0>(d));
    }
  }

  for (int r = 0; r < ctx->n_tx_cell_refinements; ++r) {
    for (auto p : tx_centers) {
      auto cell = find_active_cell_around_point(cache, p).first;

      if (cell->diameter() > ctx->min_tx_cell_size) {
        cell->set_refine_flag();
      }
    }
    ctx->mesh->execute_coarsening_and_refinement();
  }

  PetscFunctionReturn(0);
}

/**
 * @brief Generates the block sparsity pattern for a given DoFHandler object.
 *
 * This function generates the block sparsity pattern for a given DoFHandler object
 * and its associated AffineConstraints object. It creates a BlockDynamicSparsityPattern
 * object and adds entries to it for each locally owned cell in the DoFHandler. The
 * entries are added in blocks of size 2, where the first block corresponds to the
 * real part of the system and the second block corresponds to the imaginary part.
 *
 * @param dh The DoFHandler object to generate the block sparsity pattern for.
 * @param cm The AffineConstraints object associated with the DoFHandler.
 * @param bdsp The BlockDynamicSparsityPattern object to add entries to.
 * @return PetscErrorCode indicating success or failure of the function.
 */
PetscErrorCode make_block_sparsity_pattern(DoFHandler &dh, AffineConstaints &cm, dealii::BlockDynamicSparsityPattern &bdsp) {
  size_t i;
  DoFHandler::active_cell_iterator cell;
  const size_t dofs_per_cell = dh.get_fe().dofs_per_cell;
  std::vector<dealii::types::global_dof_index> dof_indices, block_dof_indices;

  PetscFunctionBegin;

  dof_indices.resize(dofs_per_cell);
  block_dof_indices.resize(dofs_per_cell * 2);

  for (cell = dh.begin_active(); cell != dh.end(); ++cell) {
    if (cell->is_locally_owned()) {
      cell->get_dof_indices(dof_indices);

      for (i = 0; i < dofs_per_cell; ++i) {
        block_dof_indices[i] = dof_indices[i];
        block_dof_indices[i + dofs_per_cell] = dof_indices[i] + dh.n_dofs();
      }

      cm.add_entries_local_to_global(block_dof_indices, bdsp, false);
    }
  }

  PetscFunctionReturn(0);
}

/**
 * @brief Generates affine constraints for a given DoFHandler object.
 *
 * This function generates affine constraints for a given DoFHandler object and
 * its associated AffineConstraints object. It applies the boundary values of
 * the given real and imaginary functions to the DoFHandler object and creates
 * constraints for hanging nodes. The constraints are then added to the
 * AffineConstraints object.
 *
 * @param dh The DoFHandler object to generate affine constraints for.
 * @param cm The AffineConstraints object to add constraints to.
 * @param f_real The real function to apply boundary values from.
 * @param f_imag The imaginary function to apply boundary values from.
 * @return PetscErrorCode indicating success or failure of the function.
 */
PetscErrorCode make_block_affine_constraints(DoFHandler &dh, AffineConstaints &cm, const dealii::Function<3> &f_real,
                                             const dealii::Function<3> &f_imag) {
  size_t i;
  dealii::IndexSet relevant_dofs;
  AffineConstaints cm_real, cm_imag;
  dealii::types::global_dof_index ndofs;
  const std::vector<std::pair<AffineConstaints::size_type, double>> *entries;
  std::vector<std::pair<AffineConstaints::size_type, double>>::const_iterator it;

  PetscFunctionBegin;

  ndofs = dh.n_dofs();
  relevant_dofs = get_locally_relevant_dofs_block(dh);

  cm_real.reinit(relevant_dofs);
  dealii::VectorTools::project_boundary_values_curl_conforming_l2(dh, 0, f_real, 0, cm_real, dealii::StaticMappingQ1<3>::mapping);
  dealii::DoFTools::make_hanging_node_constraints(dh, cm_real);
  cm_real.close();

  cm_imag.reinit(relevant_dofs);
  dealii::VectorTools::project_boundary_values_curl_conforming_l2(dh, 0, f_imag, 0, cm_imag, dealii::StaticMappingQ1<3>::mapping);
  dealii::DoFTools::make_hanging_node_constraints(dh, cm_imag);
  cm_imag.close();

  cm.reinit(relevant_dofs);
  for (i = 0; i < ndofs; ++i) {
    if (cm_real.can_store_line(i) && cm_real.is_constrained(i)) {
      cm.add_line(i);

      entries = cm_real.get_constraint_entries(i);
      for (it = entries->begin(); it != entries->end(); ++it) {
        cm.add_entry(i, it->first, it->second);
      }

      if (cm_real.is_inhomogeneously_constrained(i)) {
        cm.set_inhomogeneity(i, cm_real.get_inhomogeneity(i));
      }
    }
    if (cm_imag.can_store_line(i) && cm_imag.is_constrained(i)) {
      cm.add_line(i + ndofs);

      entries = cm_imag.get_constraint_entries(i);
      for (it = entries->begin(); it != entries->end(); ++it) {
        cm.add_entry(i + ndofs, it->first + ndofs, it->second);
      }

      if (cm_imag.is_inhomogeneously_constrained(i)) {
        cm.set_inhomogeneity(i + ndofs, cm_imag.get_inhomogeneity(i));
      }
    }
  }
  cm.close();

  PetscFunctionReturn(0);
}

/**
 * @brief Sets up the system for the forward problem.
 *
 * This function sets up the system for the electromagnetic forward problem by
 * initializing the DoFHandler object, distributing degrees of freedom, generating
 * affine constraints, creating the block sparsity pattern, and initializing the
 * matrices and vectors. It also sets the use_ams flag based on the inner_pc_type
 * parameter.
 *
 * @param ctx Pointer to the EMContext object.
 * @return PetscErrorCode indicating success or failure of the function.
 */
PetscErrorCode setup_system(EMContext *ctx) {
  dealii::IndexSet locally_owned_dofs_block, locally_relevant_dofs_block;
  std::vector<dealii::IndexSet> local_dofs_per_block(2), relevant_dofs_per_block(2);

  PetscFunctionBegin;

  LogEventHelper leh(ctx->SetupSystem);

  ctx->dh->reinit(*ctx->mesh);
  ctx->dh->distribute_dofs(*ctx->fe_nedelec);

  locally_owned_dofs_block = get_locally_owned_dofs_block(*ctx->dh);
  locally_relevant_dofs_block = get_locally_relevant_dofs_block(*ctx->dh);

  make_block_affine_constraints(*ctx->dh, *ctx->constraints, dealii::Functions::ZeroFunction<3>(3),
                                dealii::Functions::ZeroFunction<3>(3));

  local_dofs_per_block[0] = locally_owned_dofs_block.get_view(0, ctx->dh->n_dofs());
  local_dofs_per_block[1] = locally_owned_dofs_block.get_view(ctx->dh->n_dofs(), ctx->dh->n_dofs() * 2);

  relevant_dofs_per_block[0] = locally_relevant_dofs_block.get_view(0, ctx->dh->n_dofs());
  relevant_dofs_per_block[1] = locally_relevant_dofs_block.get_view(ctx->dh->n_dofs(), ctx->dh->n_dofs() * 2);

  dealii::BlockDynamicSparsityPattern bdsp(relevant_dofs_per_block);
  make_block_sparsity_pattern(*ctx->dh, *ctx->constraints, bdsp);
  dealii::SparsityTools::distribute_sparsity_pattern(bdsp, locally_owned_dofs_block, ctx->world_comm,
                                                     locally_relevant_dofs_block);
  bdsp.compress();

  ctx->K->reinit(local_dofs_per_block, local_dofs_per_block, bdsp, ctx->world_comm);
  ctx->s->reinit(local_dofs_per_block, ctx->world_comm);
  ctx->e->reinit(local_dofs_per_block, ctx->world_comm);
  ctx->ghost_e->reinit(local_dofs_per_block, relevant_dofs_per_block, ctx->world_comm);

  if ((InnerPCType)ctx->inner_pc_type == Mixed) {
    if ((PetscInt)ctx->dh->n_dofs() > ctx->pc_threshold) {
      ctx->use_ams = PETSC_TRUE;
    } else {
      ctx->use_ams = PETSC_FALSE;
    }
  } else if ((InnerPCType)ctx->inner_pc_type == AMS) {
    ctx->use_ams = PETSC_TRUE;
  } else {
    ctx->use_ams = PETSC_FALSE;
  }

  PetscFunctionReturn(0);
}

/**
 * @brief Assembles the matrix for the forward problem.
 *
 * This function assembles the matrix for the electromagnetic forward problem by
 * initializing the quadrature formula, FEValues object, and FullMatrix object.
 * It then loops over all locally owned cells, calculates the conductivity, and
 * computes the cell matrix. Finally, it distributes the local cell matrix to
 * the global matrix.
 *
 * @param ctx Pointer to the EMContext object.
 * @param fidx The index of the frequency to assemble the matrix for.
 * @return PetscErrorCode indicating success or failure of the function.
 */
PetscErrorCode assemble_matrix(EMContext *ctx, int fidx) {
  dealii::QGauss<3> quadrature_formula(ctx->order + 2);
  dealii::FEValues<3> fe_values(*ctx->fe_nedelec, quadrature_formula,
                                dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points |
                                    dealii::update_JxW_values);

  const size_t n_q_points = quadrature_formula.size();
  const size_t dofs_per_cell = ctx->fe_nedelec->dofs_per_cell;

  const dealii::FEValuesExtractors::Vector E(0);
  dealii::FullMatrix<double> cell_matrix(dofs_per_cell * 2, dofs_per_cell * 2);
  std::vector<dealii::types::global_dof_index> dof_indices(dofs_per_cell), block_dof_indices(dofs_per_cell * 2);

  double omega, sigma;
  size_t i, j, q_point;
  DoFHandler::active_cell_iterator cell;

  PetscFunctionBegin;

  LogEventHelper leh(ctx->AssembleMat);

  omega = 2 * PI * ctx->freqs[fidx];

  *(ctx->K) = 0.0;

  for (cell = ctx->dh->begin_active(); cell != ctx->dh->end(); ++cell) {
    if (!cell->is_locally_owned()) {
      continue;
    }

    sigma = 1.0 / ctx->rho[cell->material_id()];

    fe_values.reinit(cell);

    cell->get_dof_indices(dof_indices);

    for (i = 0; i < dofs_per_cell; ++i) {
      block_dof_indices[i] = dof_indices[i];
      block_dof_indices[i + dofs_per_cell] = dof_indices[i] + ctx->dh->n_dofs();
    }

    cell_matrix = 0.0;
    for (i = 0; i < dofs_per_cell; ++i) {
      for (j = 0; j < dofs_per_cell; ++j) {
        for (q_point = 0; q_point < n_q_points; ++q_point) {
          cell_matrix(i, j) += fe_values[E].curl(i, q_point) * fe_values[E].curl(j, q_point) * fe_values.JxW(q_point);
          cell_matrix(i, j + dofs_per_cell) +=
              -omega * MU * (fe_values[E].value(i, q_point) * sigma * fe_values[E].value(j, q_point)) * fe_values.JxW(q_point);
        }
        cell_matrix(i + dofs_per_cell, j + dofs_per_cell) = -cell_matrix(i, j);
        cell_matrix(i + dofs_per_cell, j) = cell_matrix(i, j + dofs_per_cell);
      }
    }

    ctx->constraints->distribute_local_to_global(cell_matrix, block_dof_indices, *ctx->K);
  }
  ctx->K->compress(dealii::VectorOperation::add);

  PetscFunctionReturn(0);
}

/**
 * @brief Assembles the right-hand side vector for the controlled-source
 * electromagnetic problem.
 *
 * This function assembles the right-hand side vector for the controlled-source
 * electromagnetic problem by initializing the quadrature formula, FEValues
 * object, and Vector object. It then find the locally owned cell that contains
 * the transmitter and distributes the local contributions to the global vector.
 *
 * @param ctx Pointer to the EMContext object.
 * @param fidx The index of the frequency to assemble the vector for.
 * @param tidx The index of the transmitter to assemble the vector for.
 * @return PetscErrorCode indicating success or failure of the function.
 */
PetscErrorCode assemble_rhs(EMContext *ctx, int fidx, int tidx) {
  size_t i, t;
  VectorD direction;
  double omega, current, weights;
  DoFHandler::active_cell_iterator cell;
  dealii::FEValuesExtractors::Vector E(0);
  size_t dofs_per_cell = ctx->fe_nedelec->dofs_per_cell;
  dealii::Vector<double> cell_rhs(dofs_per_cell * 2);
  std::vector<std::tuple<Point, VectorD, double, double>> dipoles;
  std::pair<Triangulation::active_cell_iterator, Point> cell_point;
  std::vector<dealii::types::global_dof_index> dof_indices(dofs_per_cell), block_dof_indices(dofs_per_cell * 2);

  dealii::GridTools::Cache<3> cache(*ctx->mesh);

  PetscFunctionBegin;

  LogEventHelper leh(ctx->AssembleRHS);

  omega = 2 * PI * ctx->freqs[fidx];

  *ctx->s = 0.0;

  dipoles = divide_line_current(ctx->tx[tidx], ctx->n_tx_divisions);

  for (t = 0; t < dipoles.size(); ++t) {
    cell_point = find_active_cell_around_point(cache, std::get<0>(dipoles[t]));

    if (cell_point.first.state() != dealii::IteratorState::valid) {
      SETERRQ(ctx->world_comm, EM_ERR_USER, "%s",
              string_format("Transmitter #%d not found in triangulation.", (int)tidx).c_str());
    }

    cell = DoFHandler::active_cell_iterator(*(cell_point.first), ctx->dh.get());

    if (cell->is_locally_owned()) {
      cell->get_dof_indices(dof_indices);

      for (i = 0; i < dofs_per_cell; ++i) {
        block_dof_indices[i] = dof_indices[i];
        block_dof_indices[i + dofs_per_cell] = dof_indices[i] + ctx->dh->n_dofs();
      }

      cell_rhs = 0.0;

      const dealii::Quadrature<3> quadrature(cell_point.second);
      dealii::FEValues<3> fe_values(*ctx->fe_nedelec, quadrature, dealii::update_values);

      fe_values.reinit(cell);

      direction = std::get<1>(dipoles[t]);
      current = std::get<2>(dipoles[t]);
      weights = std::get<3>(dipoles[t]);

      for (i = 0; i < dofs_per_cell; ++i) {
        cell_rhs[i + dofs_per_cell] = omega * MU * current * fe_values[E].value(i, 0) * direction * weights;
      }
      ctx->constraints->distribute_local_to_global(cell_rhs, block_dof_indices, *ctx->s);
    }
  }
  ctx->s->compress(dealii::VectorOperation::add);

  PetscFunctionReturn(0);
}

/**
 * @brief Solves the forward problem for all frequencies and transmitters.
 *
 * This function solves the forward problem for all frequencies and transmitters
 * by calling the necessary functions to set up the system, assemble the matrix
 * and right-hand side vector, solve the linear system, and calculate the
 * response. It loops over all frequencies and transmitters and saves the
 * response to a file.
 *
 * @param ctx Pointer to the EMContext object.
 * @return PetscErrorCode indicating success or failure of the function.
 */
PetscErrorCode forward(EMContext *ctx) {
  PetscFunctionBegin;

  ctx->mesh->refine_global(ctx->n_global_refinements);

  PetscCall(refine_tx_rx_area(ctx));

  for (unsigned int fidx = 0; fidx < ctx->freqs.size(); ++fidx) {
    for (unsigned int tidx = 0; tidx < ctx->tx.size(); ++tidx) {
      LogStageHelper freq_lsh(string_format("Freq-%d-Tx-%d", fidx, tidx));

      PetscCall(PetscViewerASCIIPrintf(ctx->log, "Freq %g Hz, TX %d:\n", ctx->freqs[fidx], tidx));
      PetscCall(PetscViewerASCIIPushTab(ctx->log));

      PetscCall(setup_system(ctx));

      PetscCall(PetscViewerASCIIPrintf(ctx->log, "Cells %d, DoFs %d:\n", ctx->mesh->n_global_active_cells(), ctx->dh->n_dofs()));

      PetscCall(assemble_matrix(ctx, fidx));

      PetscCall(create_preconditioner(ctx));

      PetscCall(assemble_rhs(ctx, fidx, tidx));
      PetscCall(solve_linear_system(ctx, *ctx->s, *ctx->e, ctx->K_max_it, ctx->e_rtol));
      ctx->constraints->distribute(*ctx->e);
      *ctx->ghost_e = *ctx->e;

      PetscCall(destroy_preconditioner(ctx));

      PetscCall(save_mesh(ctx, string_format("%s-csem-%02d-%02d", ctx->oprefix, fidx, tidx).c_str()));

      PetscCall(calculate_response(ctx, fidx, tidx));

      PetscCall(PetscViewerASCIIPopTab(ctx->log));
    }
  }

  PetscCall(save_rsp(ctx, string_format("%s.rsp", ctx->oprefix).c_str()));

  PetscFunctionReturn(0);
}
