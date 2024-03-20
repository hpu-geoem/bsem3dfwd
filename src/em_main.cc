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

#include "em_ctx.h"
#include "em_fwd.h"
#include "em_io.h"

int main(int argc, char **argv) {
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  EMContext ctx;

  PetscFunctionBegin;

  PetscCall(process_options(&ctx));
  PetscCall(create_context(&ctx));

  PetscCall(read_mdl(&ctx));
  PetscCall(read_emd(&ctx));

  PetscCall(forward(&ctx));

  PetscCall(destroy_context(&ctx));

  PetscFunctionReturn(0);
}
