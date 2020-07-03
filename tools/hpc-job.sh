#!/bin/bash
# Run the file given
# Submit job as:
#
#   sbatch hpc-job.sh {path/to/IPC_bin} {path/to/input}
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
##SBATCH --time=4:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=IPC

# Check if the system uses modules
if type "module" > /dev/null 2>&1; then
    # Load modules
    module purge
    module load cgal/intel/4.10b1
    module swap boost/intel/1.62.0 boost/intel/1.71.0
    module load suitesparse/intel/4.5.4 swig/gnu/3.0.11 ilmbase/intel/2.2.0 \
    openexr/intel/2.2.0 openmpi/intel/3.1.4 fftw/intel/3.3.6-pl2 \
    glew/intel/2.1.0 lapack/gnu/3.7.0
    module load cmake/intel/3.11.4 mosek/8.1.0.64 tbb/intel/2017u3
    module load gurobi/9.0.0
    module load gmp/gnu/6.1.2
    module swap gcc gcc/9.1.0
fi

printf "Job ID: %d\n" $SLURM_JOB_ID # stdout
printf "Job ID: %d\n" $SLURM_JOB_ID 1>&2 # stderr

IPC_BIN="$1"
INPUT=$(realpath "$2")
OUTPUT_SUFFIX="_$3"

NUM_THREADS=$SLURM_CPUS_PER_TASK

export MKL_NUM_THREADS=$NUM_THREADS # for Ubuntu or Mac when CHOLMOD is compiled with MKL LAPACK and BLAS
export OMP_NUM_THREADS=$NUM_THREADS # for Ubuntu when CHOLMOD is compiled with libopenblas
export VECLIB_MAXIMUM_THREADS=$NUM_THREADS # for Mac when CHOLMOD is compiled with default LAPACK and BLAS

TMP_OUTPUT_LOG="output${OUTPUT_SUFFIX}.txt"
${IPC_BIN} 100 "$INPUT" 0.999 666 4 "${OUTPUT_SUFFIX}" --logLevel 0 | tee ${TMP_OUTPUT_LOG}
OUTPUT_PATH=$(grep "output path" ${TMP_OUTPUT_LOG} | sed 's/.*output path: \(.*\)/\1/')
rm ${TMP_OUTPUT_LOG}

# Tar all files
cd $(dirname ${OUTPUT_PATH})
OUTPUT_PATH=$(basename $OUTPUT_PATH)
tar -czf "${OUTPUT_PATH%/}.tar.gz" ${OUTPUT_PATH}
rm -rf "${OUTPUT_PATH}"

# Delete any core dumps
rm core.*
