.. _midway3:

================
UChicago-Midway3
================

Midway3 is the HPC cluster of the University of Chicago, maintained by UChicago's `RCC <https://rcc.uchicago.edu/>`_. Midway3 has both GPU-accelerated nodes and CPU-only nodes.

.. code-block:: bash

   $ ssh <username>@midway3.rcc.uchicago.edu

Building WEST (CPU)
~~~~~~~~~~~~~~~~~~~

WEST executables can be compiled using the following script (tested on March 9, 2026):

.. code-block:: bash

   $ cat build_west.sh
   #!/bin/bash

   module load intel/2022.0
   module load intelmpi/2021.5+intel-2022.0
   module load mkl/2023.1
   module load python/anaconda-2024.10

   export MPIF90=mpiifort
   export F90=ifort
   export CC=icc
   export SCALAPACK_LIBS="-lmkl_scalapack_lp64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -Wl,--end-group"
   export ANACONDA_DIR=/software/python-anaconda-2024.10-el8-x86_64

   ./configure --with-scalapack=intel
   make -j 8 pw

   cd West

   make conf PYT=python3 PYT_LDFLAGS="-L$ANACONDA_DIR/lib/ -lpython3.12 -Wl,-rpath,$ANACONDA_DIR/lib/"
   make -j 8 all

To use the script do:

.. code-block:: bash

   $ bash build_west.sh


Running WEST Jobs (CPU)
~~~~~~~~~~~~~~~~~~~~~~~

The following is an example executable script `run_west.sh` to run the `wstat.x` WEST executable on two nodes of Midway3 with 48 MPI ranks per node. The <project_name> and <account_name> must be replaced with an active project allocation.

.. note::

   It is recommended to run the calculation from the scratch space (`$SCRATCH` instead of `/home`).

.. code-block:: bash

   $ cat run_west.sh
   #!/bin/bash
   #SBATCH --time=00:20:00
   #SBATCH --partition=<partition_name>
   #SBATCH --account=<account_name>
   #SBATCH --nodes=2
   #SBATCH --ntasks-per-node=48
   #SBATCH --cpus-per-task=1

   module load intel/2022.0
   module load intelmpi/2021.5+intel-2022.0
   module load mkl/2023.1
   module load python/anaconda-2024.10

   export OMP_NUM_THREADS=1

   ulimit -s unlimited

   mpirun -np 96 ./wstat.x -i wstat.in > wstat.out

Job submission is done with the following:

.. code-block:: bash

   $ sbatch run_west.sh

Building WEST (GPU)
~~~~~~~~~~~~~~~~~~~

WEST executables can be compiled using the following script (tested on March 9, 2026):

.. note::

   `--with-cuda-cc=80,cc90` below is a trick to generate executables targeting both NVIDIA Ampere (compute capability 8.0) and Hopper (compute capability 9.0) architectures, ensuring compatibility with GPUs from either generation.

.. code-block:: bash

   $ cat build_west.sh
   #!/bin/bash

   module load nvhpc/26.1
   module load cuda/12.9
   module load python/anaconda-2024.10

   export BLAS_LIBS=$NVHPC_ROOT/compilers/lib/libblas.a
   export LAPACK_LIBS=$NVHPC_ROOT/compilers/lib/liblapack.a
   export ANACONDA_DIR=/software/python-anaconda-2024.10-el8-x86_64

   ./configure --with-cuda=$NVHPC_ROOT/cuda/12.9 --with-cuda-cc=80,cc90 --with-cuda-runtime=12.9

   make -j 8 pw

   cd West

   make conf PYT=python3 PYT_LDFLAGS="-L$ANACONDA_DIR/lib/ -lpython3.12 -Wl,-rpath,$ANACONDA_DIR/lib/"
   make -j 8 all

To use the script do:

.. code-block:: bash

   $ bash build_west.sh


Running WEST Jobs (GPU)
~~~~~~~~~~~~~~~~~~~~~~~

The following is an example executable script `run_west.sh` to run the `wstat.x` WEST executable on two nodes of Midway3 with 4 MPI ranks and 4 NVIDIA A100 GPUs per node. The <project_name> and <account_name> must be replaced with an active project allocation. The directives can be adjusted to request GPUs of different types.

.. note::

   It is recommended to run the calculation from the scratch space (`$SCRATCH` instead of `/home`).

.. code-block:: bash

   $ cat run_west.sh
   #!/bin/bash
   #SBATCH --time=00:20:00
   #SBATCH --partition=<partition_name>
   #SBATCH --account=<account_name>
   #SBATCH --constraint=A100
   #SBATCH --gres=gpu:4
   #SBATCH --nodes=2
   #SBATCH --ntasks-per-node=4
   #SBATCH --cpus-per-task=8

   module load nvhpc/26.1
   module load cuda/12.9
   module load python/anaconda-2024.10

   export OMP_NUM_THREADS=1

   mpirun -np 8 ./wstat.x -i wstat.in > wstat.out

Job submission is done with the following:

.. code-block:: bash

   $ sbatch run_west.sh

.. seealso::
   For more information, visit the `RCC user guide <https://docs.rcc.uchicago.edu/>`_.
