//
//  EigenLibSolver.hpp
//  IPC
//
//  Created by Minchen Li on 6/30/18.
//
#pragma once

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
    LinSysSolverType type() const override { return LinSysSolverType::EIGEN; }

    void set_pattern(const std::vector<std::set<int>>& vNeighbor, const std::set<int>& fixedVert) override;
    void set_pattern(const Eigen::SparseMatrix<double>& mtr) override; //NOTE: mtr must be SPD

    void analyze_pattern(void) override;

    bool factorize(void) override;

    void solve(Eigen::VectorXd& rhs, Eigen::VectorXd& result) override;

    double coeffMtr(int rowI, int colI) const override;

    void setZero(void) override;

    virtual void setCoeff(int rowI, int colI, double val) override;

    virtual void addCoeff(int rowI, int colI, double val) override;

    virtual void setUnit_row(int rowI) override;

    virtual void setUnit_col(int colI, const std::set<int>& rowVIs) override;
};

} // namespace IPC
