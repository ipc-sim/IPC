# IPC Legacy Script Setting Commands

* `warmStart {0-5}`: with intersection-free constraints, warmstart barely accelerate the optimization
    * `0`: initial position (no warm start, default) (xⁿ)
    * `1`: initial position + integrated velocity (xⁿ + Δt vⁿ)
    * `2`: initial position + integrated velocity + doubly integrated gravity (xⁿ + Δt vⁿ + Δt² g)
    * `3`: uniformly accelerated warm start
    * `4`: symplectic Euler warm start
    * `5`: Jacobi warm start
* `CCDMethod {FloatingPointRootFinder | RootParity | BSC | RationalRootParity | TightInclusion | ...}`: By default, IPC uses a floating-point version of CCD together with edge-triangle intersection checks as safeguard, which is much faster than exact CCD and does not have any failure cases so far in practice
    * which method of exact CCD to use (default: `FloatingPointRootFinder`)
    * `FloatingPointRootFinder`: use the floating-point root finder method implemented by Etienne Vouga
    * `TightInclusion`: use our new provably conservative CCD [Wang et al. 2021]
    * `RootParity`: use the root parity method of Brochu et al. [2012], which in fact has failure cases in our experiments (e.g. sphere roller)!
    * `BSC`: use the Bernstein sign classification method of Tang et al. [2014], which in fact has failure cases in our experiments (e.g. sphere roller)!
    * `RationalRootParity`: use our reimplementation with rational numbers of the root parity method in Brochu et al. [2012], which in fact does not handle double collisions in a single time-step (a limitation of the original [Brochu et al. 2012])!
    * For more options see https://github.com/Continuous-Collision-Detection/CCD-Wrapper
* `tuning <positive-integer>`: please use `dHat`, `epsv`, `fricIterAmt`, and `tol` to control accuracies instead
    * number of accuracy parameters to be set followed by that number of realNumber
* `meshCO  <mesh-file-path>  <x> <y> <z>  <scale>  <stiffness>  <friction-coefficient>`
  * a legacy way of specifying static collision object, which lacks many important features, so please use `shapes` with `linearVelocity` or `angularVelocity` extra keyword instead
* `shape input <script-file-path>`
  * a legacy way of specifying a single input shape, where assigning materials, bounary conditions, or initial velocities are not supported, so please use `shapes` instead
* `shapeMatrix input  <num-x> <num-y> <num-z>`
    * a legacy way of specifying a matrix of a shape with `<num-x>`, `<num-y>`, and `<num-z>` number of that shape in each axis, followed by a line in the format of:
    * `<mesh-file-path>  <dx> <dy> <dz>  <rot-deg-x> <rot-deg-y> <rot-deg-z>  <scale-x> <scale-y>` where `<dx>`, `<dy>`, and `<dz>` sets the distance between adjacent objects in each axis.
    * assigning materials, boundary conditions, or initial velocities are not supported in this command, so please use `shapes` instead
* `rotateModel <rot-axis-x> <rot-axis-y> <rot-axis-z> <rot-deg>`
    * rotate all models as a single one in model space with rotDeg degrees along direction `[rot-axis-x rot-axis-y rot-axis-z]`
    * please orient shapes using `shapes` instead
* `size <size>`
    * size of the bounding box in meters for IPC to rescale the scene
    * `-1` (default): do not rescale, which adds less complication to scene settings
* `script {null | scaleF | hang | hang2 | swing | stamp | undstamp | stand | slip | topbottomfix | corner | push | tear | upndown | stretch | squash | stretchnsquash | bend | twist | twistnstretch | twistnsns | twistnsns_old | rubberBandPull | fourLegPull | headTailPull | onepoint | random | fall | dragdown | dragright | toggleTop | ...}`
    * Different animation scripts (see `src/AnimScripter.cpp` for details)
    * By default it is set to `null`, where all settings can be set via script in a general way (see our [quick start guide]())
* `handleRatio <ratio>`
    * a legacy way of specifying the proportion of fixed vertices for setting Dirichlet boundary conditions
    * `ratio` must be greater than zero and less than one
    * please use `shapes` command with `DBC` extra keyword instead
