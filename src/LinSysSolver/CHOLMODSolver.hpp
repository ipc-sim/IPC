//
//  CHOLMODSolver.hpp
//  IPC
//
//  Created by Minchen Li on 6/22/18.
//

#ifdef USE_CHOLMOD

#ifndef CHOLMODSolver_hpp
#define CHOLMODSolver_hpp

#include "LinSysSolver.hpp"

#include "cholmod.h"

#include <Eigen/Eigen>

#include <vector>
#include <set>

namespace IPC {

template <typename vectorTypeI, typename vectorTypeS>
class CHOLMODSolver : public LinSysSolver<vectorTypeI, vectorTypeS> {
    typedef LinSysSolver<vectorTypeI, vectorTypeS> Base;

protected:
    cholmod_common cm;
    cholmod_sparse* A;
    cholmod_factor* L;
    cholmod_dense *b, *solution;
    cholmod_dense *x_cd, *y_cd; // for multiply

    void *Ai, *Ap, *Ax, *bx, *solutionx, *x_cdx, *y_cdx;

public:
    CHOLMODSolver(void);
    ~CHOLMODSolver(void);

    void set_pattern(const std::vector<std::set<int>>& vNeighbor,
        const std::set<int>& fixedVert);
    void set_pattern(const Eigen::SparseMatrix<double>& mtr); //NOTE: mtr must be SPD
    void load(const char* filePath, Eigen::VectorXd& rhs);

    void analyze_pattern(void);

    bool factorize(void);

    void solve(Eigen::VectorXd& rhs,
        Eigen::VectorXd& result);

    virtual void multiply(const Eigen::VectorXd& x,
        Eigen::VectorXd& Ax);

    virtual void outputFactorization(const std::string& filePath);
};

} // namespace IPC

#endif /* CHOLMODSolver_hpp */

#endif /* USE_CHOLMOD */
