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

#ifndef _EM_LA_H_
#define _EM_LA_H_

#include <petsc.h>

#include "em_defs.h"

struct EMContext;

PetscErrorCode create_preconditioner(EMContext *);
PetscErrorCode destroy_preconditioner(EMContext *);

PetscErrorCode solve_linear_system(EMContext *, const PETScBlockVector &, PETScBlockVector &, PetscInt, PetscReal);

#endif
