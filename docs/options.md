# Command line arguments

The program's behavior can be controlled using command line arguments, which are
placed into an options file and passed to the program. The available arguments are
as follows:

- `-iprefix` Specifies the prefix of the input files. This argument is required.
  The input files are
  - `<iprefix>.tria` The mesh file.
  - `<iprefix>.rho` The resistivity file.
  - `<iprefix>.emd` The EM data file.

  The mesh file contains the mesh information. The resistivity file contains the
  resistivity values for each material. For the format of the mesh file and the
  resistivity file, please refer to the [Model File](/docs/model.md).

  The EM data file contains the frequencies, transmitters and receivers. For the format
  of the data file, please refer to the [Data File](/docs/data.md#data-file).

- `-oprefix` Specifies the prefix of the output files. This argument is required.
  The output files are
  - `<oprefix>.vtu` The mesh file that constains the solution.
  - `<oprefix>.rsp` The response file.

  The response file contains the response for each frequency-transmitter-receiver pair.
  For the format of the response file, please refer to the
  [Response File](/docs/data.md#response-file).

- `-order` Specifies the order of the basis function. The default is 0, which
  corresponds to the lowest order Nedelec element. The order can range from 0 to 12.
- `-inner_pc_type` This option controls the preconditioner used for solving
  the small linear system when applying the block-diagonal preconditioner.
  The available options are:
  - `direct` Uses a direct solver, with the solver type controlled by `-direct_solver_type`.
  - `ams` Uses the Auxiliary Space Maxwell (AMS) solver.
  - `mixed` Uses a direct solver when the number of DoFs is smaller than a certain
    threshold, and the AMS solver otherwise.
    The threshold is specified by `-pc_threshold`.
- `-n_rx_cell_refinements` and `-min_rx_cell_size`
  These two options specify the number of times the cells where the receivers
  are located will be refined and the minimum cell size. The default values are 0 and -1,
  respectively, meaning no refinement will be performed.
- `-n_tx_cell_refinements` and `-min_tx_cell_size`
  Similar to the above two options, but for the transmitters.
- `-n_tx_divisions` - Specifies the number of parts the finite length dipole will be divided
  into. The default is 20.
- `-n_global_refinements` - Specifies the number of times the mesh will be globally refined.
  The default is 0.

Here is an example of an options file:

```text
-iprefix half-space
-oprefix output-1/half-space

-order 0

-inner_pc_type ams

-n_rx_cell_refinements 5
-min_rx_cell_size 2

-n_tx_cell_refinements 10
-min_tx_cell_size 5
```
