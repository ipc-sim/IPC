#!/bin/bash
# Run the files located in `./input/12`
# Submit job as:
#
#   sbatch run.sh
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --time=4:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=IPC_run

# Check if the system uses modules
if type "module" > /dev/null 2>&1; then
    # Load modules
    module purge

    module load cgal/intel/4.10b1
    module swap boost/intel/1.62.0 boost/intel/1.71.0
    module load suitesparse/intel/4.5.4 swig/gnu/3.0.11 ilmbase/intel/2.2.0 \
    openexr/intel/2.2.0 openmpi/intel/3.1.4 fftw/intel/3.3.6-pl2 \
    glew/intel/2.1.0 lapack/gnu/3.7.0
    module swap gcc gcc/9.1.0
    module load cmake/intel/3.11.4 mosek/8.1.0.64 tbb/intel/2017u3
    module load gurobi/9.0.0
    module load gmp/gnu/6.1.2

    export CC=${GCC_ROOT}/bin/gcc
    export CXX=${GCC_ROOT}/bin/g++
    export CX=$CXX
fi

if [ -z "${SLURM_SUBMIT_DIR}" ]; then
    cd "${SLURM_SUBMIT_DIR}"
fi

INPUT_PATH="input"
PROG_PATH="build/IPC_bin"

for num_threads in 12; do
    export MKL_NUM_THREADS=$num_threads  # for Ubuntu or Mac when CHOLMOD is compiled with MKL LAPACK and BLAS
    export OMP_NUM_THREADS=$num_threads  # for Ubuntu when CHOLMOD is compiled with libopenblas
    export VECLIB_MAXIMUM_THREADS=$num_threads # for Mac when CHOLMOD is compiled with default LAPACK and BLAS
    for file in "$INPUT_PATH/$num_threads"/*; do
        if [ ! -d "$file" ]; then
            $PROG_PATH 100 $(realpath $file) 0.999 666 4 t$num_threads
        fi
    done
done
