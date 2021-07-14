//
//  CHOLMODSolver.cpp
//  IPC
//
//  Created by Minchen Li on 6/22/18.
//

#ifdef USE_CHOLMOD

#include "CHOLMODSolver.hpp"
#include "getRSS.hpp"

#include <iostream>

namespace IPC {

template <typename vectorTypeI, typename vectorTypeS>
CHOLMODSolver<vectorTypeI, vectorTypeS>::CHOLMODSolver(void)
{
    cholmod_start(&cm);
    A = NULL;
    L = NULL;
    b = NULL;
    x_cd = y_cd = NULL;

    Ai = Ap = Ax = NULL;
    bx = NULL;
    solutionx = x_cdx = y_cdx = NULL;
}

template <typename vectorTypeI, typename vectorTypeS>
CHOLMODSolver<vectorTypeI, vectorTypeS>::~CHOLMODSolver(void)
{
    if (A) {
        A->i = Ai;
        A->p = Ap;
        A->x = Ax;
        cholmod_free_sparse(&A, &cm);
    }

    cholmod_free_factor(&L, &cm);

    if (b) {
        b->x = bx;
        cholmod_free_dense(&b, &cm);
    }

    if (x_cd) {
        x_cd->x = x_cdx;
        cholmod_free_dense(&x_cd, &cm);
    }

    if (y_cd) {
        y_cd->x = y_cdx;
        cholmod_free_dense(&y_cd, &cm);
    }

    cholmod_finish(&cm);
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::set_pattern(const std::vector<std::set<int>>& vNeighbor,
    const std::set<int>& fixedVert)
{
    Base::set_pattern(vNeighbor, fixedVert);

    //TODO: directly save into A
    if (!A) {
        A = cholmod_allocate_sparse(Base::numRows, Base::numRows, Base::ja.size(),
            true, true, -1, CHOLMOD_REAL, &cm);
        Ax = A->x;
        Ap = A->p;
        Ai = A->i;
        // -1: upper right part will be ignored during computation
    }
    Base::ia.array() -= 1;
    Base::ja.array() -= 1; // CHOLMOD's index starts from 0
    A->i = Base::ja.data();
    A->p = Base::ia.data();
    A->x = Base::a.data();
}
template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::set_pattern(const Eigen::SparseMatrix<double>& mtr)
{
    Base::set_pattern(mtr);

    if (!A) {
        A = cholmod_allocate_sparse(Base::numRows, Base::numRows, mtr.nonZeros(),
            true, true, -1, CHOLMOD_REAL, &cm);
        Ax = A->x;
        Ap = A->p;
        Ai = A->i;
        // -1: upper right part will be ignored during computation

        A->i = Base::ja.data();
        A->p = Base::ia.data();
        A->x = Base::a.data();
    }
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::load(const char* filePath, Eigen::VectorXd& rhs)
{
    Base::load(filePath, rhs);

    //TODO: directly save into A
    if (!A) {
        A = cholmod_allocate_sparse(Base::numRows, Base::numRows, Base::ja.size(),
            true, true, -1, CHOLMOD_REAL, &cm);
        Ax = A->x;
        Ap = A->p;
        Ai = A->i;
        // -1: upper right part will be ignored during computation
    }
    Base::ia.array() -= 1;
    Base::ja.array() -= 1; // CHOLMOD's index starts from 0
    A->i = Base::ja.data();
    A->p = Base::ia.data();
    A->x = Base::a.data();
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::analyze_pattern(void)
{
    // std::cout << getCurrentRSS() << std::endl;
    cholmod_free_factor(&L, &cm);
    L = cholmod_analyze(A, &cm);
}

template <typename vectorTypeI, typename vectorTypeS>
bool CHOLMODSolver<vectorTypeI, vectorTypeS>::factorize(void)
{
    cholmod_factorize(A, L, &cm);
    // std::cout << getCurrentRSS() << std::endl;
    // exit(0);
    return cm.status != CHOLMOD_NOT_POSDEF;
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::solve(Eigen::VectorXd& rhs,
    Eigen::VectorXd& result)
{
    //TODO: directly point to rhs?
    if (!b) {
        b = cholmod_allocate_dense(Base::numRows, 1, Base::numRows, CHOLMOD_REAL, &cm);
        bx = b->x;
    }
    b->x = rhs.data();
    cholmod_dense* x;
    x = cholmod_solve(CHOLMOD_A, L, b, &cm);
    result.conservativeResize(rhs.size());
    memcpy(result.data(), x->x, result.size() * sizeof(result[0]));
    cholmod_free_dense(&x, &cm);
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::multiply(const Eigen::VectorXd& x,
    Eigen::VectorXd& Ax)
{
    assert(x.size() == Base::numRows);

    if (!x_cd) {
        x_cd = cholmod_allocate_dense(Base::numRows, 1, Base::numRows,
            CHOLMOD_REAL, &cm);
        x_cdx = x_cd->x;
    }
    x_cd->x = (void*)x.data();

    Ax.conservativeResize(Base::numRows);
    if (!y_cd) {
        y_cd = cholmod_allocate_dense(Base::numRows, 1, Base::numRows,
            CHOLMOD_REAL, &cm);
        y_cdx = y_cd->x;
    }
    y_cd->x = (void*)Ax.data();

    double alpha[2] = { 1.0, 1.0 }, beta[2] = { 0.0, 0.0 };

    cholmod_sdmult(A, 0, alpha, beta, x_cd, y_cd, &cm);
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::outputFactorization(const std::string& filePath)
{
    cholmod_sparse* spm = cholmod_factor_to_sparse(L, &cm);

    FILE* out = fopen(filePath.c_str(), "w");
    assert(out);

    cholmod_write_sparse(out, spm, NULL, "", &cm);

    fclose(out);
}

template class CHOLMODSolver<Eigen::VectorXi, Eigen::VectorXd>;

} // namespace IPC

#endif
