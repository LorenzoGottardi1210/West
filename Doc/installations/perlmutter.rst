.. _perlmutter:

================
NERSC-Perlmutter
================

Perlmutter is an HPE Cray EX supercomputer located at National Energy Research Scientific Computing Center (`NERSC <https://www.nersc.gov/>`_). Perlmutter has both GPU-accelerated nodes and CPU-only nodes.

.. code-block:: bash

   $ ssh <username>@perlmutter.nersc.gov  # or ssh <username>@saul.nersc.gov

Building WEST (GPU)
~~~~~~~~~~~~~~~~~~~

WEST executables can be compiled using the following script (tested on August 20, 2026):

.. code-block:: bash

   $ cat build_west.sh
   #!/bin/bash

   module unload darshan
   module load gpu
   module load PrgEnv-nvidia
   module load nvidia/26.5
   module load cudatoolkit/13.2
   module load craype-accel-nvidia80
   module load cray-python/3.12.12

   ./configure --with-cuda=$CUDA_HOME --with-cuda-runtime=13.2 --with-cuda-cc=80 --with-cuda-mpi=yes

   # Edit make.inc:
   sed -i 's/^MPIF90 *=.*/MPIF90 = ftn/' make.inc
   sed -i 's/^F90 *=.*/F90 = ftn/' make.inc
   sed -i 's/^CC *=.*/CC = cc/' make.inc
   sed -i 's/^LD *=.*/LD = ftn/' make.inc
   sed -i 's/^BLAS_LIBS *=.*/BLAS_LIBS =/' make.inc
   sed -i 's/^LAPACK_LIBS *=.*/LAPACK_LIBS =/' make.inc

   make -j 8 pw

   cd West

   make conf PYT=python3 PYT_LDFLAGS="-L$PYTHON_PATH/lib/ -lpython3.12 -Wl,-rpath,$PYTHON_PATH/lib/"
   make -j 8 all

To use the script do:

.. code-block:: bash

   $ bash build_west.sh

Running WEST Jobs (GPU)
~~~~~~~~~~~~~~~~~~~~~~~

The following is an example executable script `run_west.sh` to run the `wstat.x` WEST executable on two GPU nodes of Perlmutter with 4 MPI ranks and 4 GPUs per node. The <project_name> must be replaced with an active project allocation.

.. note::

   It is recommended to run the calculation from the Lustre file system (`$SCRATCH` instead of `/home`).

.. code-block:: bash

   $ cat run_west.sh
   #!/bin/bash

   #SBATCH --job-name=WEST
   #SBATCH --time=00:20:00
   #SBATCH --account=<project_name>
   #SBATCH --constraint=gpu
   #SBATCH --qos=debug
   #SBATCH --nodes=2
   #SBATCH --ntasks-per-node=4
   #SBATCH --gpus-per-node=4
   #SBATCH --cpus-per-task=32

   module unload darshan
   module load gpu
   module load PrgEnv-nvidia
   module load nvidia/26.5
   module load cudatoolkit/13.2
   module load craype-accel-nvidia80
   module load cray-python/3.12.12

   export OMP_NUM_THREADS=1
   export SLURM_CPU_BIND=cores
   export MPICH_GPU_SUPPORT_ENABLED=1

   srun -n 8 ./wstat.x -i wstat.in &> wstat.out

Job submission is done with the following:

.. code-block:: bash

   $ sbatch run_west.sh

Building WEST (CPU)
~~~~~~~~~~~~~~~~~~~

WEST executables can be compiled using the following script (tested on August 20, 2026):

.. code-block:: bash

   $ cat build_west.sh
   #!/bin/bash

   module unload darshan
   module load cpu
   module load cray-fftw/3.3.10.11
   module load cray-python/3.12.12

   export MPIF90=ftn
   export F90=ftn
   export CC=cc

   ./configure --with-scalapack

   # Edit make.inc:

   sed -i 's/^DFLAGS *=.*/DFLAGS = -D__FFTW3 -D__MPI -D__MPI_MODULE -D__SCALAPACK/' make.inc
   sed -i 's/^IFLAGS *=.*/IFLAGS = -I. -I\$(TOPDIR)\/include -I\/opt\/cray\/pe\/fftw\/3.3.10.11\/x86_milan\/include/' make.inc
   sed -i 's/^BLAS_LIBS *=.*/BLAS_LIBS =/' make.inc
   sed -i 's/^LAPACK_LIBS *=.*/LAPACK_LIBS =/' make.inc

   make -j 8 pw

   cd West

   make conf PYT=python3 PYT_LDFLAGS="-L$PYTHON_PATH/lib/ -lpython3.12 -Wl,-rpath,$PYTHON_PATH/lib/"
   make -j 8 all

To use the script do:

.. code-block:: bash

   $ bash build_west.sh

Running WEST Jobs (CPU)
~~~~~~~~~~~~~~~~~~~~~~~

The following is an example executable script `run_west.sh` to run the `wstat.x` WEST executable on two CPU nodes of Perlmutter with 128 MPI ranks per node. The <project_name> must be replaced with an active project allocation.

.. note::

   It is recommended to run the calculation from the Lustre file system (`$SCRATCH` instead of `/home`).

.. code-block:: bash

   $ cat run_west.sh
   #!/bin/bash

   #SBATCH --job-name=WEST
   #SBATCH --time=00:20:00
   #SBATCH --account=<project_name>
   #SBATCH --constraint=cpu
   #SBATCH --qos=debug
   #SBATCH --nodes=2
   #SBATCH --ntasks-per-node=128
   #SBATCH --cpus-per-task=2

   module unload darshan
   module load cpu
   module load cray-fftw/3.3.10.11
   module load cray-python/3.12.12

   export OMP_NUM_THREADS=1
   export SLURM_CPU_BIND=cores

   srun -n 256 ./wstat.x -i wstat.in &> wstat.out

Job submission is done with the following:

.. code-block:: bash

   $ sbatch run_west.sh

.. seealso::
   For more information, visit the `NERSC user guide <https://docs.nersc.gov/systems/perlmutter/architecture/>`_.
