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

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include "em_ctx.h"
#include "em_response.h"
#include "em_utils.h"
#include "em_fwd.h"

/**
 * @brief Computes the basis functions for a given point.
 *
 * This function computes the basis functions for a given point, which is used
 * to interpolate the electromagnetic fields at the receiver locations.
 *
 * @param cell_point A pair containing an active cell iterator and a point.
 * @param dof_indices A vector of global degrees of freedom indices.
 * @param basis_functions An array of vectors containing the basis function of the degrees of freedom.
 * @return PetscErrorCode An error code indicating the status of the function.
 */
PetscErrorCode compute_basis_functions(const std::pair<DoFHandler::active_cell_iterator, Point> &cell_point,
                                       std::vector<dealii::types::global_dof_index> &dof_indices,
                                       std::vector<VectorD> &basis_functions) {
  size_t i, dofs_per_cell;
  const dealii::FEValuesExtractors::Vector E(0);
  const dealii::Quadrature<3> quadrature(cell_point.second);
  dealii::FEValues<3> fe_values(cell_point.first->get_fe(), quadrature, dealii::update_values | dealii::update_gradients);

  PetscFunctionBegin;

  fe_values.reinit(cell_point.first);

  dofs_per_cell = fe_values.get_fe().dofs_per_cell;

  dof_indices.resize(dofs_per_cell);
  basis_functions.resize(dofs_per_cell);

  cell_point.first->get_dof_indices(dof_indices);

  for (i = 0; i < dofs_per_cell; ++i) {
    basis_functions[i] = fe_values[E].value(i, 0);
  }

  PetscFunctionReturn(0);
}

/**
 * @brief Interpolates the electromagnetic fields at a given point.
 *
 * This function interpolates the electromagnetic fields at a given point using
 * the provided coefficients and basis functions. The function computes the
 * electric fields at the point by summing the basis function of each degree of
 * freedom, weighted by the corresponding coefficient.
 *
 * @param dof_indices A vector of global degrees of freedom indices.
 * @param basis_functions An array of vectors containing the basis functions of the degrees of freedom.
 * @param e A PETScBlockVector containing interpolation coefficients.
 * @param es A VectorZ containing the interpolated electric field.
 * @return PetscErrorCode An error code indicating the status of the function.
 */
PetscErrorCode interpolate_fields(const std::vector<dealii::types::global_dof_index> &dof_indices,
                                  const std::vector<VectorD> &basis_functions, const PETScBlockVector &e, VectorZ &es) {
  Complex v;
  size_t i;

  PetscFunctionBegin;

  es = 0.0;
  for (i = 0; i < dof_indices.size(); ++i) {
    v = Complex(e.block(0)[dof_indices[i]], e.block(1)[dof_indices[i]]);
    es += v * basis_functions[i];
  }

  PetscFunctionReturn(0);
}

/**
 * @brief Calculates the Er response for a BSEM survey.
 *
 * This function calculates the electromagnetic response for BSEM survey using
 * the provided electromagnetic context, frequency index, and transmitter index.
 * The function computes the basis functions for each receiver location,
 * and then interpolates the electric fields at each receiver
 * using the provided degrees of freedom indices and coefficients. The
 * function then computes the Er response for each receiver location and stores it
 * in the context's response array.
 *
 * @param ctx Pointer to the EMContext object.
 * @param fidx The frequency index.
 * @param tidx The transmitter index.
 * @return PetscErrorCode An error code indicating the status of the function.
 */
PetscErrorCode calculate_response(EMContext *ctx, int fidx, int tidx) {
  VectorZ e;
  double theta;
  std::vector<VectorD> basis_values;
  std::vector<dealii::types::global_dof_index> dof_indices;

  PetscFunctionBegin;

  LogEventHelper leh(ctx->CalculateRSP);

  dealii::GridTools::Cache<3> cache(*ctx->mesh);

  for (unsigned int ridx = 0; ridx < ctx->rx.size(); ++ridx) {
    const auto cell_point = find_active_cell_around_point(cache, ctx->rx[ridx]);

    if (cell_point.first.state() != dealii::IteratorState::valid) {
      SETERRQ(ctx->world_comm, EM_ERR_USER, "%s", string_format("Receiver #%d not found in triangulation.", (int)ridx).c_str());
    }

    if (!cell_point.first->is_locally_owned()) {
      continue;
    }

    const auto cell_point_dh =
        std::make_pair(DoFHandler::active_cell_iterator(*(cell_point.first), ctx->dh.get()), cell_point.second);
    PetscCall(compute_basis_functions(cell_point_dh, dof_indices, basis_values));

    PetscCall(interpolate_fields(dof_indices, basis_values, *ctx->ghost_e, e));

    theta = std::fmod(std::atan2(ctx->rx[ridx][1] - ctx->tx[tidx].center[1], ctx->rx[ridx][0] - ctx->tx[tidx].center[0]) + 2 * PI,
                      2 * PI);

    unsigned int idx = fidx * ctx->tx.size() * ctx->rx.size() + tidx * ctx->rx.size() + ridx;
    ctx->rsp[idx] = std::cos(theta) * e[0] + std::sin(theta) * e[1];
  }

  PetscFunctionReturn(0);
}
