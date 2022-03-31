[![Build status](https://github.com/ipc-sim/IPC/workflows/Build/badge.svg?event=push)](https://github.com/ipc-sim/IPC/actions?query=workflow%3ABuild+branch%3Amaster+event%3Apush)
[![License](https://img.shields.io/github/license/ipc-sim/IPC.svg?color=blue)](https://github.com/ipc-sim/IPC/blob/master/LICENSE)

# IPC

This is the opensource reference implementation of the SIGGRAPH 2020 paper [Incremental Potential Contact: Intersection- and Inversion-free Large Deformation Dynamics](https://ipc-sim.github.io/).

## Files

* `src/`: source code
* `cmake/` and `CMakeLists.txt`: CMake files
* `Format/`: automatic code formatting (requires `clang-format`)
* `input/`: input data and scripts needed to rerun all examples in our paper, together with some tutorial examples
* `output/`: output data (will be created)
* `tests/`: unit-tests
* `tools/`: Python and Bash scripts for generating and processing results
* `wiki/`: images used in Github wiki page

* `build.py`: a python script to automatically build IPC with default compiler
* `batch.py`: a python script to automatically run a batch of examples (in `input/n/` by default) with `n` threads for each linear solver in the program
* `batch_tetgen.py`: a python script to automatically tetrahedralize a batch of closed surface meshes using [Tetgen](http://wias-berlin.de/software/tetgen/)

## Build

You can either use `python build.py` or
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
```
gcc 7 or later is recommended.

### Dependencies

* [libigl](https://libigl.github.io/), [OSQP](https://osqp.org/), [TBB](https://software.intel.com/content/www/us/en/develop/tools/threading-building-blocks.html), [spdlog](https://github.com/gabime/spdlog), [MshIO](https://github.com/qnzhou/MshIO.git), [Gulrak's Filesystem](https://github.com/gulrak/filesystem.git), and [CLI11](https://github.com/CLIUtils/CLI11): downloaded and built through CMake
* [CCD-Wrapper](https://github.com/Continuous-Collision-Detection/CCD-Wrapper): downloaded and built through CMake
    * Includes [Etienne Vouga's CTCD](https://github.com/Continuous-Collision-Detection/Floating-Point-Root-Finder) and [Tight-Inclusion](https://github.com/Continuous-Collision-Detection/Tight-Inclusion) CCD methods

#### Optional Dependencies

* [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html): must be
installed at a system level
    * Specifically IPC requires CHOLMOD, AMD, CAMD, CCOLAMD, and COLAMD
    which are all part of SuiteSparse and METIS which is sometimes a separate
    library
    * Enabled by default, and it can be disabled using the CMake argument `-DIPC_WITH_CHOLMOD=OFF`
    * **Use of the CHOLMOD linear system solver is highly recommended as it is significantly faster than the alternative solvers provided (Eigen and AMGCL).**
    * macOS: `brew install suite-sparse`
        * Alternatively, you can build SuiteSparse and then copy all files under
        `SuiteSparse/include` and `SuiteSparse/lib` to `/usr/local/include` and
        `/usr/local/lib` respectively.
    * Ubuntu: `sudo apt-get install libsuitesparse-dev libmetis-dev`
        * For program efficiency, we strongly recommend compiling SuiteSparse using [Intel MKL](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html) LAPACK and BLAS (on an Intel machine): `make library BLAS='-lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -lmkl_blas95_lp64 -liomp5' LAPACK='-lmkl_lapack95_lp64' -j 8`
    * IPC also supports using [Eigen](http://eigen.tuxfamily.org/) or [AMGCL](https://github.com/ddemidov/amgcl) as linear solver, which can be set via `IPC_LINSYSSOLVER` in `CMakeLists.txt`. To use custom linear solvers, you can implement a new interface (subclass) to our `LinSysSolver` class.
* [AMGCL](https://github.com/ddemidov/amgcl): downloaded and built through CMake
    * Enabled by default, and it can be disabled using the CMake argument `-DIPC_WITH_AMGCL=OFF`
* [Catch2](https://github.com/catchorg/Catch2): downloaded and built through CMake
    * Used for unit-tests.
    * Enabled by default iff IPC is being built as the top-level project.
<!-- * [GMP LIB](https://gmplib.org/): `sudo apt install libgmp-dev` on Ubuntu. -->

### Without OpenGL

If your system does not have the required OpenGL libraries, you can disable the viewer using the CMake argument `-DIPC_WITH_OPENGL=OFF`. It is important to then run IPC in offline mode (see below).

### 2D?

IPC is only implemented for 3D simulation.

### HPC

If your system uses Slurm like our cluster does, then you can take advantage of our scripts to build IPC as a batch job.

From the IPC root directory use `sbatch build.sh` to compile IPC via a
batch job. You can check the status of the job using
`squeue -u $USER -n IPC_build`.

It is important to not build on the login node because while it supports AVX2,
the compute nodes may not. Building on the login node will enable AVX2
support which may cause an `Illegal instruction` error when run on a compute
node. Alternatively, this can be fixed by using Intel's compilers, `icc` and `icpc`, with the
flag `-axCORE-AVX2`. This enables executables to run on compute nodes with or
without AVX2 instructions. If the compute nodes support AVX2 instructions,
executables will run AVX2 instructions, otherwise it will use other available
instructions such as AVX or SSE4.2.

### Sub-Projects

IPC contains two optional sub-projects for unit tests, debug, and file processing. To enable them use the CMake
flag `-DIPC_BUILD_<sub-project-name>_PROJECT=On` where `sub-project-name` is the name
of the sub-project (e.g. `DIAGNOSTIC` or `MESH_PROCESSING`). You can also set these interactively using `ccmake`.

## Run

Please see our [quick start guide](https://github.com/ipc-sim/IPC/wiki) for a `hello world` example with various settings. The output of each script will be saved into a separate folder in `output/`.

### Offline Mode

It is possible to run IPC with and without the viewer. By default, the `batch.py` script runs IPC with the viewer. If you provide the argument `--offline` to `batch.py` then it will run IPC in offline mode (i.e. without the viewer).

If you are running the `IPC_bin` executable directly then the first argument controls the mode. Mode `10` runs IPC with the viewer and is the default in the `batch.py`. Mode `100` runs IPC in offline mode (without the viewer). See `IPC_bin --help` for more detail.

### HPC

From the IPC root directory use `sbatch run.sh` to processes all the input
files in `input/12/` using 12 CPU cores.

This will process all input scenes in series. I am currently, working on a
script to process all scenes in parallel via separate batch jobs.

### Sub-Projects

The executables for the sub-projects are placed in their respective directory
in `src/Projects/`.

## Output Files
* `<n>.obj`: surface mesh of the simulated objects and volumetric/surface kinematic collision objects in n-th time step
* `n.seg`: mesh file of the segment kinematic collision objects if any in n-th time step
* `pt<n>.obj`: vertex coordinates of the point cloud kinematic collision objects if any in n-th time step
* `status<n>`: saved information of n-th time step that can be used to restart simulation from
* `ACO<m>_<n>.obj`: mesh file of m-th analytical plane if any in n-th (and the following if unchanged) time steps
* `BC_<n>.txt`: node index of Dirichlet nodes if any in n-th (and the following if unchanged)  time steps
* `0.png`: snapshot of initial time step
* `finalResult.png`: snapshot of last time step
* `finalResult_mesh.msh`: rest shape of simulated object
* `anim.gif`: animation preview of the simulation
* `config.txt`: the config file containing all simulation settings that can be used to rerun the simulation
* `iterStats.txt`: simulation/optimization statistics of each iteration
* `info.txt`: timing information
* `log.txt`: debug information
* `sysE.txt`: total energy (kinematic, gravity, and elasticity) of each simulated object
* `sysL.txt`: linear momentum of each simulated object
* `sysM.txt`: angular momentum of each simulated object


## Script Settings

See our [quick start guide](https://github.com/ipc-sim/IPC/wiki) for detailed examples of minimal settings needed for a simulation.

* `energy {FCR | NH}`
    * `NH`: Neo-Hookean elasticity energy (default)
    * `FCR`: Fixed Co-rotational elasticity energy
* `timeIntegration {BE | NM} <β> <γ>`
    * `BE`: Backward Euler (default)
    * `NM`: Newmark (optional: add β and γ to change from default of 0.25 and 0.5)
* `dampingStiff <dampingStiff>`
    * add lagged Rayleigh damping with stiffness `dampingStiff` or not if `dampingStiff` is zero or not specified
    * `dampingStiff` must be greater than or equal to zero
* `dampingRatio <dampingRatio>`
    * set `dampingStiff` to be `0.75 * dampingRatio * h^3` or `0` if `dampingRatio` is zero or not specified
    * `dampingRatio` must be greater than or equal to zero
* `time <num-seconds> <time-step>`
    * total time in seconds, 5 by default
    * time-step, h, in seconds, 0.025 by default
* `density <density>`
    * density of the mesh in Kg/m³, 1000 by default
* `stiffness <young-modulus> <poisson-ratio>`
    * Young's modulus (E) in Pa, 1e5 by default
    * Poisson's ratio (ν), 0.4 by default
* `turnOffGravity`
    * no gravity if specified
* `shapes input <num-of-shapes>`
    * number of input shapes followed by that number of lines specifying the shape paths and initial configuration in model coordinates:
    * `<mesh-file-path>  <x> <y> <z>  <rot-deg-x> <rot-deg-y> <rot-deg-z>  <scale-x> <scale-y> <scale-z>  <extra-command> ...`
    * `<extra-command>`: assigning different materials, initial velocities, boundary conditions, or scripted motions for each shape individually if any
      * `material <density> <young-modulus> <poisson-ratio>`
      * `initVel <lvx> <lvy> <lvz> <avx> <avy> <avz>`
      * `DBC <left-bottom-back-x> <left-bottom-back-y> <left-bottom-back-z> <right-top-front-x> <right-top-front-y> <right-top-front-z>  <lvx> <lvy> <lvz> <avx> <avy> <avz> [<t_begin>] [<t_end>]`
      * `NBC <left-bottom-back-x> <left-bottom-back-y> <left-bottom-back-z> <right-top-front-x> <right-top-front-y> <right-top-front-z>  <ax> <ay> <az> [<t_begin>] [<t_end>]`
      * `linearVelocity <lvx> <lvy> <lvz>`
      * `angularVelocity <avx> <avy> <avz>`
* `DBCTimeRange <t_begin> <t_end>`
  * only apply all Dirichlet BC set in `shapes` command during `(t_begin, t_end]`, by default `(0, +inf]`
* `NBCTimeRange <t_begin> <t_end>`
  * only apply all Neuman BC set in `shapes` command during `(t_begin, t_end]`, by default `(0, +inf]`
* `ground <friction-coefficient> <height>`
* `halfSpace <x> <y> <z>  <nx> <ny> <nz>   <stiffness>  <friction-coefficient>`
    * create a half-space analytical collision object centered at (x, y, z)
      with normal (nx, ny, nz).
    * `stiffness` is unused
* `selfCollisionOff`
    * self-collision is not handled if specified
* `selfFric <selfFricCoeff>`: specify the friction coefficient of self-contact, by default `0` (no friction)
  * `selfFricCoeff` must be larger or equal to `0`, setting it larger than `1` is allowed
* `dHat <dHat>`: set the distance that IPC sees objects as attaching and exert contact forces, by default it is `1e-3` that of the diagonal length of the bounding box of the scene
  * `dHat` must be larger than `0`
* `epsv <epsv>`: set the relative tangent velocity that IPC sees the attached objects as not sliding and exert static friction forces, by default it is `1e-3` that of the diagonal length of the bounding box of the scene per second
  * `epsv` must be larger than `0`
* `fricIterAmt <fricIterAmt>`: set the maximal friction iterations that IPC updates friction tangent and normal forces, by default it is `1`
  * if set to `0` or negative, IPC will update friction tangent and normal forces until convergence (not guaranteed so not recommended)
* `tol positiveInteger`
    * number of time steps to set for Newton tolerance other than default 1e-2 followed by that number of positiveRealNumber
    * `positiveInteger` is also a relative value on the bounding box diagonal length
* `constraintSolver {interiorPoint | QP | SQP}`
    * collision constraint solver
    * `interiorPoint`: IPC (default)
    * `QP`: solve elasticity as a single QP with nonlinear contact constraints
    * `SQP`: solve elasticity as a sequence of quadratic programs with nonlinear contact constraints
* `restart <status-file-path-str>`
    * restart simulation from the time step in the specified status file (other settings must be the same)
* `view {orthographic | perspective}`
* `zoom <positive-real-number>`
* `cameraTracking`
  * automatically place the focus of the camera at the center of the bounding box of the scene at each frame if specified
* `playBackSpeed <playBackSpeed>`
  * playback speed of `anim.gif`, which also controls how many frame data will be output (smaller value results in more output (at most output every frame))
* `appendStr <str>`
    * append result folder with the string specified
* `disableCout`
    * no command line output if specified
* `section {interiorPoint | QP | penalty}`
    * section of settings only used when the specified constraint solver is used
    * note: the constraint solver must be declared before the section
    * a section ends when the line `section end` is encountered or at the end
      of a file.
    * example:
      ```
      constraintSolver QP
      # ...
      section QP
      constraintType graphics
      section end
      ```

### Legacy Settings

Some of our settings are legacy settings or the ones that are not usually used/recommended. We include them [here](input/paperExamples/README.md) for completeness as they have been used by our paper examples or tests.

### QP/SQP Settings

The following settings only affect the QP and SQP contact constraint solvers we use for comparisons.

* `QPSolver {OSQP | Gurobi}`
    * Quadratic programming solver to use when using SQP constraint solver
    * [OSQP](https://osqp.org/) is downloaded through CMake automatically
    * [Gurobi](https://www.gurobi.com/) must be installed manually and requires a license
* `constraintType {volume | graphics | gapFunction | nonsmoothNewmark | CMR | Verschoor}`
    * which type of collision constraint to use with the SQP solver (default: `volume`)
    * `volume`: measure the volume of a tetrahedron comprised of the colliding points
    * `graphics`: distance constraint from Harmon et al. [2008].
    * `nonsmoothNewmark`: variation of volume constraints proposed by Kane et al. [1999].
    * `gapFunction`: signed distance along the normal direction from the point of intersection.
    * `CMR`: constraint manifold refinement (same as standard graphics approach).
    * `Verschoor`: Variant of standard graphics approach by Verschoor and Jalba [2019].
* `constraintOffset <offset>`
    * surface thickness used for mesh collision objects when using `constraintSolver QP`/`SQP`
* `useActiveSetConvergence`
    * Optimization can only converge if the active set converged to a fixed set
      of constraints or no constraints are active.

## Comparisons

We conducted a lengthy comparison of our method and several other commercial and
open-source software, and "vanilla" implementations of other methods. You can
read the full details in [Supplement B](https://ipc-sim.github.io/file/IPC-supplement-B-comparisons.pdf).

You can find the scenes, meshes, and other files used in these comparisons in
`input/paperExamples/supplementB`.

* `SQPBenchmark`: The collection of scene files used to generate the benchmark
of QP and SQP methods.
    * The all benchmark scenes can be run with default settings using
    `tools/run-comparison-benchmark-sqp.sh`
    * The full benchmark sweep can be run using `tools/run-comparison-benchmark-hpc.sh`
* `COMSOL`: The 2D COMSOL files used in our comparison
* `Houdini`: The Houdini files used to generate the chain comparison
* `SOFA`: The SOFA scene files used to compare against [SOFA](https://www.sofa-framework.org/)
* `Utopia`: The mesh and IPC script used to compare against Utopia [Krause and Zulian 2016]

## Tests

The `tests` directory contains a very limited number of unit-tests. These tests
use [Catch2](https://github.com/catchorg/Catch2) which is downloaded through
CMake unless the option `-DIPC_WITH_TESTS=OFF` is provided. Currently, these
test contain test the collision constraints used by the SQP solver.

## Tools

The `tools` directory contains several tools for generating meshes, running
scenes, and processing results.

* `convert_to_ipc_msh.py`: this script uses [PyMesh](https://pymesh.readthedocs.io/en/latest/)
to convert a tetrahedral mesh to the expected `.msh` format of IPC
    * If the input mesh is a triangle mesh, then Tetgen will be used to tetrahedralize it
* `geo_to_msh.py`: this script uses [PyMesh](https://pymesh.readthedocs.io/en/latest/)
to convert a GEO (format exported by Houdini) format mesh to a standard MSH mesh
* `generate_chain_scene.py`: generate a IPC script for a variable length chain
* `run-comparison-benchmark-sqp.sh`: run the comparison SQP benchmark with default settings
* `run-comparison-benchmark-hpc.sh`: run the comparison full sweep of the SQP benchmark
    * dispatches jobs by calling the `hpc-job.sh` script
* `process_*_results.py`: process the benchmark results, generating a CSV of the table
