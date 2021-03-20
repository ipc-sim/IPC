//
//  EigenLibSolver.hpp
//  IPC
//
//  Created by Minchen Li on 6/30/18.
//

#ifndef EigenLibSolver_hpp
#define EigenLibSolver_hpp

#include "LinSysSolver.hpp"

#include <Eigen/Eigen>

#include <vector>
#include <set>

namespace IPC {

template <typename vectorTypeI, typename vectorTypeS>
class EigenLibSolver : public LinSysSolver<vectorTypeI, vectorTypeS> {
    typedef LinSysSolver<vectorTypeI, vectorTypeS> Base;

protected:
    Eigen::SparseMatrix<double> coefMtr;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> simplicialLDLT;

public:
    void set_pattern(const std::vector<std::set<int>>& vNeighbor,
        const std::set<int>& fixedVert);
    void set_pattern(const Eigen::SparseMatrix<double>& mtr); //NOTE: mtr must be SPD

    void analyze_pattern(void);

    bool factorize(void);

    void solve(Eigen::VectorXd& rhs,
        Eigen::VectorXd& result);

    double coeffMtr(int rowI, int colI) const;

    void setZero(void);

    virtual void setCoeff(int rowI, int colI, double val);

    virtual void addCoeff(int rowI, int colI, double val);

    virtual void setUnit_row(int rowI);

    virtual void setUnit_col(int colI, const std::set<int>& rowVIs);
};

} // namespace IPC

#endif /* EigenLibSolver_hpp */
