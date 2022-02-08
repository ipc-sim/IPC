//
//  LinSysSolver.cpp
//  IPC
//

#include "ipc/LinSysSolver/LinSysSolver.hpp"
#include "ipc/LinSysSolver/AMGCLSolver.hpp"
#include "ipc/LinSysSolver/CHOLMODSolver.hpp"
#include "ipc/LinSysSolver/EigenLibSolver.hpp"

#include <spdlog/spdlog.h>

namespace IPC {

template <typename vectorTypeI, typename vectorTypeS>
LinSysSolver<vectorTypeI, vectorTypeS>* LinSysSolver<vectorTypeI, vectorTypeS>::create(const LinSysSolverType type)
{
    switch (type) {
    case LinSysSolverType::CHOLMOD:
        return new CHOLMODSolver<vectorTypeI, vectorTypeS>();
    case LinSysSolverType::AMGCL:
        return new AMGCLSolver<vectorTypeI, vectorTypeS>();
    case LinSysSolverType::EIGEN:
        return new EigenLibSolver<vectorTypeI, vectorTypeS>();
    default:
        spdlog::error("Uknown linear system solver type: {}", type);
        throw "Uknown linear system solver type: " + type;
    }
}

template class LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>;

}
