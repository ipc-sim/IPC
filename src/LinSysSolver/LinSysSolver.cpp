#include "LinSysSolver.hpp"
#include "AMGCLSolver.hpp"
#include "CHOLMODSolver.hpp"
#include "EigenLibSolver.hpp"

#include <spdlog/spdlog.h>

namespace IPC {

template <typename vectorTypeI, typename vectorTypeS>
LinSysSolver<vectorTypeI, vectorTypeS>* LinSysSolver<vectorTypeI, vectorTypeS>::create(const LinSysSolverType type)
{
    switch (type) {
#ifdef IPC_WITH_CHOLMOD
    case LinSysSolverType::CHOLMOD:
        return new CHOLMODSolver<vectorTypeI, vectorTypeS>();
#endif
#ifdef IPC_WITH_AMGCL
    case LinSysSolverType::AMGCL:
        return new AMGCLSolver<vectorTypeI, vectorTypeS>();
#endif
    case LinSysSolverType::EIGEN:
        return new EigenLibSolver<vectorTypeI, vectorTypeS>();
    default:
        spdlog::error("Uknown linear system solver type: {}", type);
        throw fmt::format("Uknown linear system solver type: {}", type);
    }
}

template class LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>;

}
