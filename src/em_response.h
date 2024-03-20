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

#ifndef _EM_RESPONSE_H_
#define _EM_RESPONSE_H_

#include "em_defs.h"

struct EMContext;

PetscErrorCode compute_basis_functions(const std::pair<DoFHandler::active_cell_iterator, Point> &,
                                       std::vector<dealii::types::global_dof_index> &, std::vector<VectorD> &);
PetscErrorCode interpolate_fields(const std::vector<dealii::types::global_dof_index> &,
                                  const std::vector<VectorD> &basis_functions, const PETScBlockVector &, VectorZ &);

PetscErrorCode calculate_response(EMContext *, int, int);

#endif
