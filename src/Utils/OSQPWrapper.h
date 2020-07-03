//
//  OSQP.h
//  IPC
//
//  Created by Minchen Li on 10/15/18.
//

#ifndef OSQP_h
#define OSQP_h

#include <osqp.h>

namespace IPC {

class OSQP {
    OSQPSettings* settings; // Problem settings
    OSQPData* data; // OSQPData
    OSQPWorkspace* work; // Workspace

public:
    OSQP(bool printOutput = false)
    {
        settings = (OSQPSettings*)c_malloc(sizeof(OSQPSettings));
        // Define Solver settings as default
        osqp_set_default_settings(settings);
        settings->verbose = printOutput;
        settings->eps_abs = 1.0e-20;
        settings->eps_rel = 1.0e-3;
        settings->scaled_termination = true; // does not mean eps_abs can be adaptive
#ifdef OSQP_USE_MKL_PARDISO
        settings->linsys_solver = MKL_PARDISO_SOLVER;
#endif

        data = (OSQPData*)c_malloc(sizeof(OSQPData));

        //TODO: update matrix using the osqp_update API to save time
        data->P = NULL;
        data->A = NULL;

        work = NULL;
    }
    ~OSQP(void)
    {
        if (work) {
            osqp_cleanup(work);
        }
        if (data->P) {
            c_free(data->P);
        }
        if (data->A) {
            c_free(data->A);
        }
        c_free(data);
        c_free(settings);
    }

    void setup(c_float* P_x, c_int P_nnz, c_int* P_i, c_int* P_p,
        c_float* q,
        c_float* A_x, c_int A_nnz, c_int* A_i, c_int* A_p,
        c_float* l, c_float* u,
        c_int n, c_int m)
    {
        if (work) {
            osqp_cleanup(work);
        }
        if (data->P) {
            c_free(data->P);
        }
        if (data->A) {
            c_free(data->A);
        }

        // Populate data
        data->n = n;
        data->m = m;
        data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
        data->q = q;
        data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
        data->l = l;
        data->u = u;

        // Setup workspace
        // c_int err = osqp_setup(&work, data, settings);
        work = osqp_setup(data, settings);
        // assert(!err);
    }

    c_int solve(void)
    {
        // Solve Problem
        osqp_solve(work);
        // solution is in work->solution
        // primal: c_float *work->solution->x
        // dual: c_float *work->solution->y
        return work->info->status_val;
    }

    c_float* getPrimal(void) const
    {
        assert(work);
        return work->solution->x;
    }
    c_float* getDual(void) const
    {
        assert(work);
        return work->solution->y;
    }
};

} // namespace IPC

#endif /* OSQP_h */
