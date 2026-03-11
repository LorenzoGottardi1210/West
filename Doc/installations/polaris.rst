.. _polaris:

============
ALCF-Polaris
============

Polaris is a GPU-accelerated supercomputer located at Argonne National Laboratory, maintained by `ALCF <https://www.alcf.anl.gov/>`_.

.. code-block:: bash

   $ ssh <username>@polaris.alcf.anl.gov

Building WEST
~~~~~~~~~~~~~

WEST executables can be compiled using the following script (tested on February 3, 2026):

.. code-block:: bash

   $ cat build_west.sh
   #!/bin/bash

   module load cuda/12.9
   module load craype-accel-nvidia80
   module load nvidia/25.5
   module load cray-libsci/25.09.0
   module load cray-python/3.11.7

   ./configure --with-cuda=$NVIDIA_PATH/cuda/12.9 --with-cuda-runtime=12.9 --with-cuda-cc=80 --with-cuda-mpi=yes

   # Edit make.inc:
   sed -i 's/^MPIF90 *=.*/MPIF90 = ftn/' make.inc
   sed -i 's/^F90 *=.*/F90 = ftn/' make.inc
   sed -i 's/^CC *=.*/CC = cc/' make.inc
   sed -i 's/^LD *=.*/LD = ftn/' make.inc
   sed -i 's/^BLAS_LIBS *=.*/BLAS_LIBS =/' make.inc
   sed -i 's/^LAPACK_LIBS *=.*/LAPACK_LIBS =/' make.inc

   make -j 8 pw

   cd West

   make conf PYT=python3 PYT_LDFLAGS="-L$PYTHON_PATH/lib/ -lpython3.11 -Wl,-rpath,$PYTHON_PATH/lib/"
   make -j 8 all

To use the script do:

.. code-block:: bash

   $ bash build_west.sh

Running WEST Jobs
~~~~~~~~~~~~~~~~~

The following is an example executable script `run_west.sh` to run the `wstat.x` WEST executable on two nodes of Polaris with 4 MPI ranks and 4 GPUs per node. The <project_name> must be replaced with an active project allocation.

.. note::

   It is recommended to run the calculation from one of the Lustre file systems (`/grand` or `/eagle` instead of `/home`).

.. code-block:: bash

   $ cat run_west.sh
   #!/bin/bash -l
   #PBS -l select=1:system=polaris
   #PBS -l place=scatter
   #PBS -l walltime=0:20:00
   #PBS -l filesystems=home:grand
   #PBS -j oe
   #PBS -q debug
   #PBS -A <project_name>
   #PBS -N job_name

   module load cuda/12.9
   module load craype-accel-nvidia80
   module load nvidia/25.5
   module load cray-libsci/25.09.0
   module load cray-python/3.11.7

   export MPICH_GPU_SUPPORT_ENABLED=1
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PYTHON_PATH/lib

   NNODES=`wc -l < $PBS_NODEFILE`
   NRANKS_PER_NODE=$(nvidia-smi -L | wc -l)
   NDEPTH=8
   NTHREADS=1
   NTOTRANKS=$(( NNODES * NRANKS_PER_NODE ))

   cd ${PBS_O_WORKDIR}

   mpiexec -n ${NTOTRANKS} --ppn ${NRANKS_PER_NODE} --depth=${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS=${NTHREADS} -env OMP_PLACES=threads ./wstat.x -i wstat.in &> wstat.out

Job submission is done with the following:

.. code-block:: bash

   $ qsub run_west.sh

.. seealso::
   For more information, visit the `ALCF user guide <https://docs.alcf.anl.gov/polaris/getting-started/>`_.
