# Steps to run the example

This example demonstrates how to run the forward modeling process for a half-space
model with different orders of elements.

Run the following command to generate the model and data files, and run the forward
modeling process:

```shell
spack load dealii # This command is not needed if the user is using the Docker image
python make-model.py
python make-emdata.py
mpirun -np <number of processors> ../../bsem3dfwd -options_file order-1/half-space.cfg
mpirun -np <number of processors> ../../bsem3dfwd -options_file order-2/half-space.cfg
mpirun -np <number of processors> ../../bsem3dfwd -options_file order-3/half-space.cfg
mpirun -np <number of processors> ../../bsem3dfwd -options_file refine-1/half-space.cfg
mpirun -np <number of processors> ../../bsem3dfwd -options_file refine-2/half-space.cfg
```

Please replace `<number of processors>` with the number of processors you want to use.
