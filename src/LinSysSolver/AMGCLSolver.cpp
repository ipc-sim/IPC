//
//  AMGCLSolver.hpp
//  IPC
//
//  Created by Minchen Li on 11/06/19.
//

#ifdef USE_AMGCL

#include "AMGCLSolver.hpp"
#include "getRSS.hpp"

#include <amgcl/io/mm.hpp>

#include <iostream>

namespace IPC {

template <typename vectorTypeI, typename vectorTypeS>
AMGCLSolver<vectorTypeI, vectorTypeS>::AMGCLSolver(void)
{
    solver = NULL;
}

template <typename vectorTypeI, typename vectorTypeS>
AMGCLSolver<vectorTypeI, vectorTypeS>::~AMGCLSolver(void)
{
    if (solver) {
        delete solver;
    }
}

template <typename vectorTypeI, typename vectorTypeS>
void AMGCLSolver<vectorTypeI, vectorTypeS>::set_pattern(const std::vector<std::set<int>>& vNeighbor,
    const std::set<int>& fixedVert)
{
    Base::set_pattern(vNeighbor, fixedVert);
    Base::ia.array() -= 1;
    Base::ja.array() -= 1;
    copyOffDiag_IJ();
}

template <typename vectorTypeI, typename vectorTypeS>
void AMGCLSolver<vectorTypeI, vectorTypeS>::load(const char* filePath, Eigen::VectorXd& rhs)
{
    Base::load(filePath, rhs);

    Base::ia.array() -= 1;
    Base::ja.array() -= 1;

    //_ia.resize(Base::ia.size());
    //std::memcpy(_ia.data(), Base::ia.data(), Base::ia.size() * sizeof(Base::ia[0]));
    //_ja.resize(Base::ja.size());
    //std::memcpy(_ja.data(), Base::ja.data(), Base::ja.size() * sizeof(Base::ja[0]));
    //_a.resize(Base::a.size());
    //std::memcpy(_a.data(), Base::a.data(), Base::a.size() * sizeof(Base::a[0]));
    copyOffDiag_IJ();
}

template <typename vectorTypeI, typename vectorTypeS>
void AMGCLSolver<vectorTypeI, vectorTypeS>::load_AMGCL(const char* filePath, Eigen::VectorXd& rhs)
{
    int rows, cols, n, m;
    std::tie(rows, cols) = amgcl::io::mm_reader(filePath)(_ia, _ja, _a);
    Base::numRows = rows;

    std::vector<double> rhs_stdv;
    std::tie(n, m) = amgcl::io::mm_reader(std::string(filePath) + "_rhs")(rhs_stdv);
    rhs.resize(rhs_stdv.size());
    std::memcpy(rhs.data(), rhs_stdv.data(), sizeof(rhs_stdv[0]) * rhs_stdv.size());
}

template <typename vectorTypeI, typename vectorTypeS>
void AMGCLSolver<vectorTypeI, vectorTypeS>::write_AMGCL(const char* filePath, const Eigen::VectorXd& rhs) const
{
    amgcl::io::mm_write(filePath, std::tie(Base::numRows, _ia, _ja, _a));
    amgcl::io::mm_write(std::string(filePath) + "_rhs", rhs.data(), rhs.size());
}

template <typename vectorTypeI, typename vectorTypeS>
void AMGCLSolver<vectorTypeI, vectorTypeS>::copyOffDiag_IJ(void)
{
    std::vector<std::vector<int>> ja_new(Base::numRows);
    for (int rowI = 0; rowI < Base::numRows; ++rowI) {
        int jaStart = Base::ia[rowI];
        while (Base::ja[jaStart] < rowI) {
            ++jaStart;
        }
        int jaEnd = Base::ia[rowI + 1];
        assert(jaEnd >= jaStart);

        int curSize = ja_new[rowI].size();
        ja_new[rowI].resize(curSize + jaEnd - jaStart);
        std::memcpy(ja_new[rowI].data() + curSize, Base::ja.data() + jaStart,
            sizeof(Base::ja[0]) * (jaEnd - jaStart));

        for (int jaI = jaStart + 1; jaI < jaEnd; ++jaI) {
            int colI = Base::ja[jaI];
            ja_new[colI].emplace_back(rowI);
        }
    }

    _ia.resize(Base::numRows + 1);
    _ia[0] = 0;
    _ja.resize(0);
    _ja.reserve(Base::ja.size() * 2);
    for (int rowI = 0; rowI < Base::numRows; ++rowI) {
        _ia[rowI + 1] = _ia[rowI] + ja_new[rowI].size();
        _ja.insert(_ja.end(), ja_new[rowI].begin(), ja_new[rowI].end());
    }
}

template <typename vectorTypeI, typename vectorTypeS>
void AMGCLSolver<vectorTypeI, vectorTypeS>::copyOffDiag_a(void)
{
    std::vector<std::vector<double>> a_new(Base::numRows);
    for (int rowI = 0; rowI < Base::numRows; ++rowI) {
        int jaStart = Base::ia[rowI];
        while (Base::ja[jaStart] < rowI) {
            ++jaStart;
        }
        int jaEnd = Base::ia[rowI + 1];
        assert(jaEnd >= jaStart);

        int curSize = a_new[rowI].size();
        a_new[rowI].resize(curSize + jaEnd - jaStart);
        std::memcpy(a_new[rowI].data() + curSize, Base::a.data() + jaStart,
            sizeof(Base::a[0]) * (jaEnd - jaStart));

        for (int jaI = jaStart + 1; jaI < jaEnd; ++jaI) {
            int colI = Base::ja[jaI];
            a_new[colI].emplace_back(Base::a[jaI]);
        }
    }

    _a.resize(0);
    _a.reserve(Base::a.size() * 2);
    for (int rowI = 0; rowI < Base::numRows; ++rowI) {
        _a.insert(_a.end(), a_new[rowI].begin(), a_new[rowI].end());
    }
}

template <typename vectorTypeI, typename vectorTypeS>
void AMGCLSolver<vectorTypeI, vectorTypeS>::analyze_pattern(void)
{
}

template <typename vectorTypeI, typename vectorTypeS>
bool AMGCLSolver<vectorTypeI, vectorTypeS>::factorize(void)
{
    // std::cout << getCurrentRSS() << std::endl;
    copyOffDiag_a();
    if (solver) {
        delete solver;
    }
    Solver::params prm;
    prm.solver.tol = 1e-5; // relative
    prm.solver.maxiter = 1000;
    prm.precond.coarsening.aggr.eps_strong = 0.0;
    prm.solver.M = 100;
#ifdef USE_BW_BACKEND
    auto A = amgcl::adapter::block_matrix<value_type>(std::tie(Base::numRows, _ia, _ja, _a));
    std::cout << A.rows() << " " << A.cols() << std::endl; //NOTE: [a wierd bug] must keep this line as it is
    solver = new Solver(A, prm);
#else
#if defined(USE_BW_AGGREGATION)
    prm.precond.coarsening.aggr.block_size = DIM;
#endif
    solver = new Solver(std::tie(Base::numRows, _ia, _ja, _a), prm);
#endif
    std::cout << solver->precond() << std::endl;
    // std::cout << getCurrentRSS() << std::endl;
    // exit(0);
    return true;
}

template <typename vectorTypeI, typename vectorTypeS>
void AMGCLSolver<vectorTypeI, vectorTypeS>::solve(Eigen::VectorXd& rhs,
    Eigen::VectorXd& result)
{
    assert(rhs.size() == Base::numRows);

#ifdef USE_BW_BACKEND ///

    result.setZero(Base::numRows);
    size_t iters;
    double resid;
    rhs_type const* fptr = reinterpret_cast<rhs_type const*>(rhs.data());
    rhs_type* xptr = reinterpret_cast<rhs_type*>(result.data());
    amgcl::backend::numa_vector<rhs_type> F(fptr, fptr + rhs.size() / DIM);
    amgcl::backend::numa_vector<rhs_type> X(xptr, xptr + result.size() / DIM);
#ifdef USE_AMG_SOLVER
    solver->precond().cycle(F, X);
#else
    std::tie(iters, resid) = (*solver)(F, X);
    std::cout << "Iterations: " << iters << std::endl
              << "Error:      " << resid << std::endl;
    if (resid != resid) {
        std::cout << "AMGCL failed!" << std::endl;
        Base::write("Mnan", rhs);
        write_AMGCL("M", rhs);
        exit(-1);
    }
#endif
    std::copy(X.data(), X.data() + X.size(), xptr);

#else ///

    std::vector<double> x(Base::numRows, 0.0), _rhs(Base::numRows);
    std::memcpy(_rhs.data(), rhs.data(), sizeof(rhs[0]) * rhs.size());
#ifdef USE_AMG_SOLVER
    solver->precond().cycle(_rhs, x);
#else
    int iters;
    double error;
    std::tie(iters, error) = (*solver)(_rhs, x);
    std::cout << "Iterations: " << iters << std::endl
              << "Error:      " << error << std::endl;
    if (error != error) {
        std::cout << "AMGCL failed!" << std::endl;
        Base::write("Mnan", rhs);
        write_AMGCL("M", rhs);
        exit(-1);
    }
#endif
    result.resize(x.size());
    std::memcpy(result.data(), x.data(), sizeof(x[0]) * x.size());

#endif ///
}

template class AMGCLSolver<Eigen::VectorXi, Eigen::VectorXd>;

} // namespace IPC

#endif
