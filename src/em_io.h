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

#ifndef _EM_IO_H_
#define _EM_IO_H_

#include <petsc.h>

#include "em_defs.h"

struct EMContext;

PetscErrorCode read_mdl(EMContext *);
PetscErrorCode read_emd(EMContext *);

PetscErrorCode save_mesh(EMContext *, const char *);
PetscErrorCode save_rsp(EMContext *, const char *);

#endif
