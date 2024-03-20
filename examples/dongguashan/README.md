# Steps to run the example

This example demonstrates how to run the forward modeling process for the
DongGuaShan model.

Step into each subdirectory and run the following command to generate the model
and data files, and run the forward modeling process:

```shell
spack load dealii # This command is not needed if the user is using the Docker image
python make-model.py
python make-emdata.py
mpirun -np <number of processors> ../../../bsem3dfwd -options_file output-1/dongguashan.cfg
```

Please replace `<number of processors>` with the number of processors you want to use.
