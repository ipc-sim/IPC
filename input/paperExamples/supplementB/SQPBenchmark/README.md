# Supplement B: SQP Benchmark

This directory contains the complete collection of scenes used to generate the table at the end of "Comparison Supplement to Incremental Potential Contact: Intersection-and Inversion-free, Large-Deformation Dynamics."

Each file contains the default settings for each scene. This does not include the varying time-steps, constraint offsets, constraint types, or constraint solvers. We used these files in combination with `tools/run-comparison-benchmark-hpc.sh` to generate a jobs on our HPC for each cell in the table.

You can run all of these script scenes with default SQP settings using `tools/run-comparison-benchmark-sqp.sh`.
