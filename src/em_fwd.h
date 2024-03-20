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

#ifndef _EM_SOLVER_H
#define _EM_SOLVER_H

#include <petsc.h>

#include "em_defs.h"

#include <deal.II/base/function.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/lac/block_sparsity_pattern.h>

struct EMContext;

dealii::IndexSet get_locally_owned_dofs_block(const DoFHandler &);
dealii::IndexSet get_locally_relevant_dofs_block(const DoFHandler &);

std::pair<Triangulation::active_cell_iterator, Point> find_active_cell_around_point(const dealii::GridTools::Cache<3> &,
                                                                                    const Point &);

PetscErrorCode refine_tx_rx_area(EMContext *);

PetscErrorCode make_block_affine_constraints(DoFHandler &, AffineConstaints &, const dealii::Function<3> &,
                                             const dealii::Function<3> &);
PetscErrorCode make_block_sparsity_pattern(DoFHandler &dh, AffineConstaints &, dealii::BlockDynamicSparsityPattern &);

PetscErrorCode setup_system(EMContext *);
PetscErrorCode assemble_matrix(EMContext *, int);
PetscErrorCode assemble_rhs(EMContext *, int, int);

PetscErrorCode forward(EMContext *);

#endif
