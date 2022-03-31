//
//  Type.hpp
//  IPC
//
//  Created by Minchen Li on 2/5/18.
//

#ifndef Types_hpp
#define Types_hpp

#include <cstdio>

//#define OSQP_USE_MKL_PARDISO

#define USE_TBB
#ifdef USE_TBB
#include <tbb/parallel_for.h>
#include <tbb/info.h>
#define TBB_NUM_THREADS (tbb::info::default_concurrency())
#endif

#define USE_IQRSVD

// interior point method for collision handling
#define HOMOTOPY_VAR 1 // 0: kappa, 1: dHat
#define BARRIER_FUNC_TYPE 2 // 0: C0 clamped log, 2: C2 clamped log

#define USE_SH_LFSS // enable spatial hashing for largestFeasibleStepSize
#define USE_SH_CCS // enable spatial hashing for computeConstraintSets
#define USE_SH_INTERSECTED // enable spatial hashing for isIntersected
// #define PARALLEL_SH_CONSTRUCT // construct spatial hash in parallel
// #define CHECK_ET_INTERSECTION

#define CFL_FOR_CCD 2
// #define OUTPUT_CCD_FAIL
// #define NO_CCD_FAILSAFE
// #define CHECK_RATIONAL_CCD
// #define CHECK_RATIONAL_CCD_GLOBAL

#define USE_DISCRETE_CMS
#define ADAPTIVE_KAPPA
#define SFCLAMPING_ORDER 1 // 0: C0 clamping, 1: C1, 2: C2
// #define CHECK_FRIC_TANGENT

// Export friction data
// WARNING: Enabling this will disable the always project zero DBC in Mesh::isProjectDBCVertex
// #define EXPORT_FRICTION_DATA

#include "Triplet.hpp"

#endif /* Types_hpp */
