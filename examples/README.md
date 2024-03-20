# Examples

This directory contains several examples of the forward modeling process for
different models. Each example is contained in a subdirectory and has its own
README file. The examples are:

- half-space: A simple half-space model to verify the accuracy of the
  implementation
- borehole: A borehole model to investigate the effect of the borehole
  resistivity
- dongguashan: A model of the Donghuashan porphyry copper deposit in China
- oil-gas: A model of a hydrocarbon reservoir

We ran the examples using a Linux machine with an Intel(R) Xeon(R) Silver 4314
CPU and 256 GB of RAM. Below we show the approximate runtime and memory usage
for each example using 1 process:

- half-space
  - order-1: less than 1 minute, 0.5 GB
  - order-2: 10 minutes, 4 GB
  - order-3: 90 minutes, 25 GB
  - refine-1: 2 minutes, 4 GB
  - refine-2: 20 minutes, 35 GB
- borehole: 3 hours, 40 GB
- dongguashan: 5.5 hours, 50 GB
- oil-gas: 6.5 hours, 80 GB

The examples can also be run in parallel using MPI. However, to achieve optimal
speedup, the number of processes should be chosen carefully. Each MPI process
should be assigned with a sufficiently large number of DoFs, at least 100,000,
to ensure the overhead associated with parallelization remains minimal. Note
that the threshold for DoFs can differ across computational platforms due to
variations in memory and network bandwidths. Please feel free to experiment with
different numbers of processes and see what works best for you.
