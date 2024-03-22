# BSEM3DFWD

## About

BSEM3DFWD is a parallel C++ program designed for the forward modeling of the 3D
borehole-to-surface electromagnetic (BSEM) method using the finite element
method. The main features of BSEM3DFWD including:

- Modular design, making it easy to extend
- Utilization of octree mesh to discretize small features and complex structures
- Implementation of high-order finite elements to improve accuracy
- Parallel computation using MPI
- Utilization of an efficient iterative solver to solve the linear system

## Installing Prerequisites and Building

BSEM3DFWD is developed based on the open-source library
[deal.II](https://www.dealii.org/), which is a C++ program library that
provides building blocks for the finite element method.
We provide three ways to install deal.II and its dependencies: using `Docker`,
using `spack`, and using our pre-built binaries.

Note that BSEM3DFWD has only been tested on Linux and macOS systems. For
Windows users, we recommend using Docker, WSL1/WSL2 or a virtual machine.

Please extract the source code of BSEM3DFWD into a directory
`/path/to/bsem3dfwd-sources` before installing the prerequisites and building.
In the following sections, we will use `/path/to/bsem3dfwd-sources` to refer to
the source directory of BSEM3DFWD.

Among the above three approaches, we recommend using `Docker`, since it is the
most robust and portable way to build BSEM3DFWD.

### Using Docker

If the user is familiar with `Docker`, we provide a Docker image that contains
all the dependencies of BSEM3DFWD. For macOS and Windows users, please install
[Docker Desktop](https://www.docker.com/products/docker-desktop) first. For
Linux users, please install `Docker Engine` following the
[installation instructions](https://docs.docker.com/engine/install).

There is an alternative Docker client called `OrbStack` for macOS users.
It is claimed to be faster than `Docker Desktop` and has a better user
interface. Please refer to [OrbStack](https://orbstack.dev) for more details.

After installing `Docker`, the user can pull the Docker image from
[Docker Hub](https://hub.docker.com/r/adamqc/bsem3dfwd-deps) using the following
command:

```shell
docker pull adamqc/bsem3dfwd-deps:latest
```

Then run the Docker container and mount the source directory of BSEM3DFWD to
`/mnt` using the following command:

```shell
docker run -it --rm -v /path/to/bsem3dfwd-sources:/mnt adamqc/bsem3dfwd-deps:latest
```

The above command will start a shell inside the Docker container. Now the user
can build BSEM3DFWD:

```shell
cd /mnt && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make
```

### Using spack

To build deal.II and its dependencies from source, we recommend using
[spack](https://spack.io), which is a package manager for supercomputers,
Linux, and macOS. Please refer to the [spack documentation](https://spack.readthedocs.io/en/latest/)
for more details.

Before we install `spack`, we need to install a C++ compiler and other
dependencies. For debian based Linux distributions, e.g., Ubuntu, Mint, etc.,
please use the following command to install GCC and other necessary packages:

```shell
sudo apt install -y build-essential gfortran git python3
```

For RedHat based Linux distributions, e.g., RHEL, Centos, etc., please use the
following command:

```shell
sudo yum install -y gcc gcc-c++ gcc-gfortran git python3
```

On macOS, the user can install `Xcode` from the App Store, or install
`Command Line Tools for Xcode` manually.

Note that `deal.II` requires a C++ compiler that supports C++17 standard, e.g.,
GCC >= 9.0, Clang >= 10.0, and AppleClang >= 12.0. The user can check the
version of the compiler using `g++ --version` or `clang++ --version`.
Also note that AppleClang >= 15.0 has a problem compiling `deal.II`, so we
recommend using `Xcode` or `Command Line Tools for Xcode` prior to version 14.3.

The next step is to install `spack`. We need the developing version of spack
since it contains the latest version of these packages. The following commands
can be used to install `spack` and `deal.II`:

```shell
# Clone the spack repository
git clone https://github.com/spack/spack.git
# Initialize spack environment
source ./spack/share/spack/setup-env.sh
# Install deal.II and its dependencies
spack install dealii@master+mpi+petsc+python~p4est~arborx~arpack~slepc~gsl~adol-c~hdf5~gmsh~sundials~oce~cgal~assimp~symengine~examples~ginkgo~threads~muparser~vtk build_type=Release ^python@3.10
```

Note that the `+python` and `+petsc` options are required since BSEM3DFWD uses
Python to generate the mesh and PETSc to solve the linear system. This command
may take a while to finish since it will build all the packages from source,
please be patient.

Once all the dependencies are installed, the following commands can be used to
build BSEM3DFWD:

```shell
# Go to the source directory of BSEM3DFWD
cd /path/to/bsem3dfwd-sources
# Load environment variables of deal.II
spack load dealii
# Create a build directory
mkdir build
# Go to the build directory
cd build
# Configure
cmake -DCMAKE_BUILD_TYPE=Release ..
# Build
make
```

The above commands will generate an executable file `bsem3dfwd` in the source
directory, which is the main program of BSEM3DFWD.

### Using pre-built binaries

We also provide pre-built binaries for Linux and macOS users. The binaries can
be downloaded from [Prebuilt Binaries](https://www.dropbox.com/scl/fo/fkfa8vp912gl9klt8qxvx/h?rlkey=2blw6mzrnt1t0lvzeamtbu27p&dl=0).

We provide three versions: `bsem3dfwd-deps-centos7.tar.xz` for RedHat based
Linux distributions, `bsem3dfwd-deps-ubuntu18.04.tar.xz` for Debian based Linux
distributions, and `bsem3dfwd-deps-ventura.tar.xz` for macOS. Note that all the
binaries are built on x86 architecture, so they may not work if the user is
using different architectures, e.g., ARM, etc.

Please download the corresponding file for your system and extract it into
`/usr/local` using the following command:

```shell
cd /usr/local
tar -xf /path/to/bsem3dfwd-deps-<system>.tar.xz
```

Regular users may not have permission to write to `/usr/local`, so the user
may need to use `sudo` to run the `tar` command.

Then load the environment variables by running the following command:

```shell
source /usr/local/bsem3dfwd-deps/setvars.sh
```

The executable file `bsem3dfwd` is installed in `/usr/local/bsem3dfwd-deps/bin`,
which is already in the `PATH` environment variable. The user can run it directly.

Alternatively, the user can also build BSEM3DFWD from source using the following
commands:

```shell
cd /path/to/bsem3dfwd-sources
mkdir build && cd build
CXX=g++-9.5 cmake -DCMAKE_BUILD_TYPE=Release .. && make # For Linux
CXX=clang++ cmake -DCMAKE_BUILD_TYPE=Release .. && make # For macOS
```

## Usage

To run the forward modeling process, the user needs to provide a configuration file
to specify the parameters. Then use the flowing command to run BSEM3DFWD:

```shell
mpirun -np <number of processors> /path/to/bsem3dfwd -options_file path/to/options.cfg
```

For details about the options file, please refer to the [Options file](/docs/options.md).

## Author

BSEM3DFWD is developed by Ce Qin, Henan Polytechnic University, China.

## License

BSEM3DFWD is distributed under the MIT License. Please refer to the [LICENSE](/LICENSE)
file for more details.

## Contributing

Users are encouraged to open an issue for questions or bugs. Pull requests for
any enhancements are also welcome.
