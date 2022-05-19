//
//  CHOLMODSolver.cpp
//  IPC
//
//  Created by Minchen Li on 6/22/18.
//

#include "ipc/LinSysSolver/CHOLMODSolver.hpp"
#include "ipc/Utils/getRSS.hpp"

#include <string>
#include <unordered_map>
#include <iostream>
#include <spdlog/spdlog.h>

extern const std::string outputFolderPath;

bool ValidateCholmodStatus(int status, std::string step) {
    static const std::unordered_map<int, std::string>& cholmod_errors = 
        *new std::unordered_map<int, std::string> {{CHOLMOD_NOT_INSTALLED, "method not installed"},
                                                   {CHOLMOD_OUT_OF_MEMORY, "out of memory"},
                                                   {CHOLMOD_TOO_LARGE, "integer overflow occured"},
                                                   {CHOLMOD_INVALID, "invalid input"},
                                                   {CHOLMOD_GPU_PROBLEM, "GPU fatal error"},
                                                   {CHOLMOD_NOT_POSDEF, "matrix not pos"},
                                                   {CHOLMOD_DSMALL, "LL' has tiny absolute value"}};
    if (status == CHOLMOD_OK)
        return true;
    const std::string error_message = cholmod_errors.at(status);
    spdlog::error("CHOLMOD failure at " + step + ": " + error_message);
    FILE* out = fopen((outputFolderPath + "linear_solver_failure.txt").c_str(), "w");
    assert(out);
    fprintf(out, "cholmod\n");
    fprintf(out, "failed at %s \n", step.c_str());
    fprintf(out, "status id %d \n", status);
    fprintf(out, "status error: %s \n", error_message.c_str());
    fclose(out);
    exit(-1);
    return false;
}

namespace IPC {

template <typename vectorTypeI, typename vectorTypeS>
CHOLMODSolver<vectorTypeI, vectorTypeS>::CHOLMODSolver(void)
{
    cholmod_start(&cm);
    ValidateCholmodStatus(cm.status, "CHOLMODSolver::CHOLMODSolver->cholmod_start");
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
        ValidateCholmodStatus(cm.status, "CHOLMODSolver::~CHOLMODSolver->cholmod_free_sparse");
    }

    cholmod_free_factor(&L, &cm);
    ValidateCholmodStatus(cm.status, "CHOLMODSolver::~CHOLMODSolver->cholmod_free_factor");

    if (b) {
        b->x = bx;
        cholmod_free_dense(&b, &cm);
        ValidateCholmodStatus(cm.status, "CHOLMODSolver::~CHOLMODSolver->cholmod_free_dense");
    }

    if (x_cd) {
        x_cd->x = x_cdx;
        cholmod_free_dense(&x_cd, &cm);
        ValidateCholmodStatus(cm.status, "CHOLMODSolver::~CHOLMODSolver->cholmod_free_dense");
    }

    if (y_cd) {
        y_cd->x = y_cdx;
        cholmod_free_dense(&y_cd, &cm);
        ValidateCholmodStatus(cm.status, "CHOLMODSolver::~CHOLMODSolver->cholmod_free_dense");
    }

    cholmod_finish(&cm);
    ValidateCholmodStatus(cm.status, "CHOLMODSolver::~CHOLMODSolver->cholmod_finish");
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
        ValidateCholmodStatus(cm.status, "CHOLMODSolver::set_pattern->cholmod_allocate_sparse");
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
        ValidateCholmodStatus(cm.status, "CHOLMODSolver::set_pattern->cholmod_allocate_sparse");
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
        ValidateCholmodStatus(cm.status, "CHOLMODSolver::load->cholmod_allocate_sparse");
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
    ValidateCholmodStatus(cm.status, "CHOLMODSolver::analyze_pattern->cholmod_free_factor");
    L = cholmod_analyze(A, &cm);
    ValidateCholmodStatus(cm.status, "CHOLMODSolver::analyze_pattern->cholmod_analyze");
}

template <typename vectorTypeI, typename vectorTypeS>
bool CHOLMODSolver<vectorTypeI, vectorTypeS>::factorize(void)
{
    cholmod_factorize(A, L, &cm);
    ValidateCholmodStatus(cm.status, "CHOLMODSolver::factorize->cholmod_factorize");
    // std::cout << getCurrentRSS() << std::endl;
    // exit(0);
    return cm.status == CHOLMOD_OK;
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::solve(Eigen::VectorXd& rhs,
    Eigen::VectorXd& result)
{
    //TODO: directly point to rhs?
    if (!b) {
        b = cholmod_allocate_dense(Base::numRows, 1, Base::numRows, CHOLMOD_REAL, &cm);
        ValidateCholmodStatus(cm.status, "CHOLMODSolver::solve->cholmod_allocate_dense");
        bx = b->x;
    }
    b->x = rhs.data();
    cholmod_dense* x;
    x = cholmod_solve(CHOLMOD_A, L, b, &cm);
    ValidateCholmodStatus(cm.status, "CHOLMODSolver::solve->cholmod_solve");
    result.conservativeResize(rhs.size());
    memcpy(result.data(), x->x, result.size() * sizeof(result[0]));
    cholmod_free_dense(&x, &cm);
    ValidateCholmodStatus(cm.status, "CHOLMODSolver::solve->cholmod_free_dense");
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::multiply(const Eigen::VectorXd& x,
    Eigen::VectorXd& Ax)
{
    assert(x.size() == Base::numRows);

    if (!x_cd) {
        x_cd = cholmod_allocate_dense(Base::numRows, 1, Base::numRows,
            CHOLMOD_REAL, &cm);
        ValidateCholmodStatus(cm.status, "CHOLMODSolver::multiply->cholmod_allocate_dense");
        x_cdx = x_cd->x;
    }
    x_cd->x = (void*)x.data();

    Ax.conservativeResize(Base::numRows);
    if (!y_cd) {
        y_cd = cholmod_allocate_dense(Base::numRows, 1, Base::numRows,
            CHOLMOD_REAL, &cm);
        ValidateCholmodStatus(cm.status, "CHOLMODSolver::multiply->cholmod_allocate_dense");
        y_cdx = y_cd->x;
    }
    y_cd->x = (void*)Ax.data();

    double alpha[2] = { 1.0, 1.0 }, beta[2] = { 0.0, 0.0 };

    cholmod_sdmult(A, 0, alpha, beta, x_cd, y_cd, &cm);
    ValidateCholmodStatus(cm.status, "CHOLMODSolver::multiply->cholmod_sdmult");
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::outputFactorization(const std::string& filePath)
{
    cholmod_sparse* spm = cholmod_factor_to_sparse(L, &cm);
    ValidateCholmodStatus(cm.status, "CHOLMODSolver::outputFactorization->cholmod_factor_to_sparse");

    FILE* out = fopen(filePath.c_str(), "w");
    assert(out);

    cholmod_write_sparse(out, spm, NULL, "", &cm);
    ValidateCholmodStatus(cm.status, "CHOLMODSolver::outputFactorization->cholmod_write_sparse");

    fclose(out);
}

template class CHOLMODSolver<Eigen::VectorXi, Eigen::VectorXd>;

} // namespace IPC
