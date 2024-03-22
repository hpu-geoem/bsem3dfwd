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

#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_dgq.h>

#include <deal.II/lac/block_vector.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <string>

#include "em_ctx.h"
#include "em_utils.h"

/**
 * @brief Reads the model data from the input file and stores it in the
 * EMContext object.
 *
 * This function reads the model data from the input file and stores it in the
 * EMContext object. The input file is expected to be in the format specified in
 * the documentation.
 *
 * @param ctx Pointer to the EMContext object.
 * @return PetscErrorCode 0 if successful, non-zero otherwise.
 */
PetscErrorCode read_mdl(EMContext *ctx) {
  PetscFunctionBegin;

  std::ifstream ifs_tria(std::string(ctx->iprefix) + ".tria");

  if (!ifs_tria.good()) {
    SETERRQ(ctx->world_comm, EM_ERR_USER, "%s", string_format("Unable to open file %s.tria.", ctx->iprefix).c_str());
  }
  boost::archive::binary_iarchive ia(ifs_tria);

  ctx->mesh->clear();
  ctx->mesh->load(ia, 0);

  std::ifstream ifs_rho(std::string(ctx->iprefix) + ".rho");

  if (!ifs_rho.good()) {
    SETERRQ(ctx->world_comm, EM_ERR_USER, "%s", string_format("Unable to open file %s.rho.", ctx->iprefix).c_str());
  }

  std::string l;
  std::stringstream ss;
  while (std::getline(ifs_rho, l)) {
    ss << parse_string(l);
  }

  unsigned int n_attrs;
  ss >> n_attrs;

  ctx->rho.resize(n_attrs);
  for (unsigned int i = 0; i < n_attrs; ++i) {
    ss >> ctx->rho[i];
  }

  PetscFunctionReturn(0);
}

/**
 * @brief Reads the frequencies, transmitters, receivers from the input file
 * and stores it in the EMContext object.
 *
 * The input file is expected to be in the format specified in the documentation.
 *
 * @param ctx Pointer to the EMContext object.
 * @return PetscErrorCode 0 if successful, non-zero otherwise.
 */
PetscErrorCode read_emd(EMContext *ctx) {
  unsigned int n_freqs, n_rxes, n_txes, i;

  PetscFunctionBegin;

  std::ifstream ifs(std::string(ctx->iprefix) + ".emd");
  if (!ifs.good()) {
    SETERRQ(ctx->world_comm, EM_ERR_USER, "%s", string_format("Unable to open file %s.emd.", ctx->iprefix).c_str());
  }

  std::string l;
  std::stringstream ss;
  while (std::getline(ifs, l)) {
    ss << parse_string(l);
  }

  ss >> n_freqs;
  ctx->freqs.resize(n_freqs);
  for (i = 0; i < n_freqs; ++i) {
    ss >> ctx->freqs[i];
  }

  ss >> n_txes;
  ctx->tx.resize(n_txes);
  for (i = 0; i < n_txes; ++i) {
    ss >> ctx->tx[i].center[0] >> ctx->tx[i].center[1] >> ctx->tx[i].center[2];
    ss >> ctx->tx[i].azimuth >> ctx->tx[i].dip >> ctx->tx[i].current >> ctx->tx[i].length;
  }

  ss >> n_rxes;
  ctx->rx.resize(n_rxes);
  for (i = 0; i < n_rxes; ++i) {
    ss >> ctx->rx[i][0] >> ctx->rx[i][1] >> ctx->rx[i][2];
  }

  ctx->rsp.resize(n_freqs * n_txes * n_rxes);

  PetscFunctionReturn(0);
}

/**
 * @brief Saves the mesh data to a file in VTU format.
 *
 * This function saves the mesh data to a file in VTU format. The mesh data
 * includes the cell-wise values of the resistivity and the real and imaginary
 * components of the electric field. The output file is named after the input
 * filename with a ".vtu" extension.
 *
 * @param ctx Pointer to the EMContext object.
 * @param fn Name of the input file.
 * @return PetscErrorCode 0 if successful, non-zero otherwise.
 */
PetscErrorCode save_mesh(EMContext *ctx, const char *fn) {
  dealii::Vector<double> partition, rho;
  dealii::DataOut<3> data_out;
  DoFHandler::active_cell_iterator cell;

  PetscFunctionBegin;

  LogEventHelper leh(ctx->SaveMesh);

  rho.reinit(ctx->mesh->n_active_cells());
  partition.reinit(ctx->mesh->n_active_cells());
  for (cell = ctx->dh->begin_active(); cell != ctx->dh->end(); ++cell) {
    rho[cell->active_cell_index()] = ctx->rho[cell->material_id()];
    partition[cell->active_cell_index()] = ctx->mesh->locally_owned_subdomain();
  }

  dealii::DataPostprocessorVector<3> e_real("e_real", dealii::update_gradients | dealii::update_gradients |
                                                          dealii::update_quadrature_points);

  data_out.attach_dof_handler(*ctx->dh);
  data_out.add_data_vector(rho, "rho", dealii::DataOut<3>::type_cell_data);
  data_out.add_data_vector(*ctx->dh, ctx->ghost_e->block(0), "e_real");
  data_out.add_data_vector(*ctx->dh, ctx->ghost_e->block(1), "e_imag");
  data_out.build_patches();

  data_out.write_vtu_in_parallel((std::string(fn) + ".vtu").c_str(), ctx->world_comm);

  PetscFunctionReturn(0);
}

/**
 * @brief Saves the response data to a file.
 *
 * This function saves the response data to a file in the format specified in
 * the documentation. The response data includes the frequencies, transmitters,
 * receivers, and the real and imaginary components of the response.
 *
 * @param ctx Pointer to the EMContext object.
 * @param fn Name of the output file.
 * @return PetscErrorCode 0 if successful, non-zero otherwise.
 */
PetscErrorCode save_rsp(EMContext *ctx, const char *fn) {
  FILE *fp;

  PetscFunctionBegin;

  LogEventHelper leh(ctx->SaveRSP);

  PetscCall(MPI_Allreduce(MPI_IN_PLACE, &ctx->rsp[0], ctx->rsp.size() * 2, MPI_DOUBLE, MPI_SUM, ctx->world_comm));

  PetscCall(PetscFOpen(ctx->world_comm, fn, "w", &fp));

  for (unsigned int fidx = 0; fidx < ctx->freqs.size(); ++fidx) {
    for (unsigned int tidx = 0; tidx < ctx->tx.size(); ++tidx) {
      for (unsigned int ridx = 0; ridx < ctx->rx.size(); ++ridx) {
        unsigned int idx = fidx * ctx->tx.size() * ctx->rx.size() + tidx * ctx->rx.size() + ridx;
        PetscFPrintf(ctx->world_comm, fp, "%E %3d % E % E % E % E % E\n", ctx->freqs[fidx], tidx, ctx->rx[ridx][0],
                     ctx->rx[ridx][1], ctx->rx[ridx][2], std::real(ctx->rsp[idx]), std::imag(ctx->rsp[idx]));
      }
    }
  }

  PetscCall(PetscFClose(ctx->world_comm, fp));

  PetscFunctionReturn(0);
}
