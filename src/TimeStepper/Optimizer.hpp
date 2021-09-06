//
//  Optimizer.hpp
//  IPC
//
//  Created by Minchen Li on 8/31/17.
//

#ifndef Optimizer_hpp
#define Optimizer_hpp

#include "Types.hpp"
#include "Energy.hpp"
#include "AnimScripter.hpp"
#include "Config.hpp"
#include "LinSysSolver.hpp"
#include "SpatialHash.hpp"

#include "OSQPWrapper.h"
#ifdef USE_GUROBI
#include <Gurobi.h>
#endif

#include <fstream>

namespace IPC {

// a class for solving an optimization problem
template <int dim>
class Optimizer {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    friend class Mesh<dim>;

protected: // referenced data
    const Mesh<dim>& data0; // initial guess
    const std::vector<Energy<dim>*>& energyTerms; // E_0, E_1, E_2, ...
    const std::vector<double>& energyParams; // a_0, a_1, a_2, ...
    // E = \sum_i a_i E_i
    Config animConfig;

protected: // owned data
    bool mute;
    int globalIterNum;
    double relGL2Tol;
    Mesh<dim> result; // intermediate results of each iteration
    bool useGD = false;

    // SPD solver for solving the linear system for search directions
    LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>*linSysSolver, *dampingMtr;

    bool solveWithQP; ///< @brief Use QP with nonlinear contact constraints.
    bool solveWithSQP; ///< @brief Use Sequential QP with nonlinear contact constraints.
    bool solveQPSuccess; ///< @brief Was the previous QP solve successful?
    OSQP OSQPSolver;
#ifdef USE_GUROBI
    Eigen::GurobiSparse gurobiQPSolver;
#endif
    Eigen::SparseMatrix<double> P_QP;
    std::vector<double*> elemPtr_P_QP;
    Eigen::VectorXd dual_QP;
    std::vector<int> constraintStartInds;
    std::vector<std::vector<int>> activeSet, activeSet_next;
    std::vector<std::vector<int>> activeSet_lastH;
    std::vector<Eigen::VectorXd> lambda_lastH;

    bool solveIP, solveFric;
    double kappa, dHatEps, dHat;
    double fbNormTol, bboxDiagSize2, dTolRel, dTol, dHatTarget;
    double fricDHat, fricDHatThres, fricDHat0, fricDHatTarget;
    std::vector<std::vector<MMCVID>> MMActiveSet, MMActiveSet_next;
    std::unordered_map<MMCVID, double, MMCVIDHash> mesh_mmcvid_to_toi;
    std::vector<std::vector<MMCVID>> paraEEMMCVIDSet;
    std::vector<std::vector<std::pair<int, int>>> paraEEeIeJSet;
    std::vector<std::vector<std::pair<int, int>>> MMActiveSet_CCD;
    std::vector<std::vector<MMCVID>> MMActiveSet_lastH;
    std::vector<Eigen::VectorXd> MMLambda_lastH;
    std::vector<std::vector<Eigen::Vector2d>> MMDistCoord;
    std::vector<std::vector<Eigen::Matrix<double, 3, 2>>> MMTanBasis;
    Eigen::VectorXi n_collPairs_sum, n_collPairs_max, n_convCollPairs_sum, n_convCollPairs_max;
    int n_collPairs_total_max, n_convCollPairs_total_max;
    SpatialHash<dim> sh;
    std::vector<std::pair<int, int>> closeConstraintID;
    std::vector<std::pair<int, MMCVID>> closeMConstraintID;
    std::vector<double> closeConstraintVal, closeMConstraintVal;
    double CN_MBC, rho_DBC;
    std::vector<std::set<int>> vNeighbor_IP;
    bool m_projectDBC;
    std::set<Triplet> initSF;

    Eigen::VectorXd gradient; // energy gradient computed in each iteration
    Eigen::VectorXd searchDir; // search direction comptued in each iteration
    double lastEnergyVal; // for output and line search
    double targetGRes;
    std::vector<Eigen::VectorXd> gradient_ET;
    std::vector<double> energyVal_ET;

    std::ofstream file_iterStats;
    std::ofstream file_sysE, file_sysM, file_sysL;

    int numOfLineSearch;

#ifdef EXPORT_FRICTION_DATA
    Eigen::MatrixXd savedFrictionLaggedVertices;
#endif

protected: // dynamic information
    Eigen::VectorXd velocity;
    Eigen::MatrixXd xTilta, dx_Elastic, acceleration;
    double dt, dtSq;
    Eigen::Matrix<double, dim, 1> gravity, gravityDtSq;
    int frameAmt;
    AnimScripter<dim> animScripter;
    int innerIterAmt;
    std::vector<AutoFlipSVD<Eigen::Matrix<double, dim, dim>>> svd;
    std::vector<Eigen::Matrix<double, dim, dim>> F;
    double beta_NM; // \in [1/4,1/2], default 1/4
    double gamma_NM; // \in [0,1], default 1/2

public: // constructor and destructor
    Optimizer(const Mesh<dim>& p_data0,
        const std::vector<Energy<dim>*>& p_energyTerms,
        const std::vector<double>& p_energyParams,
        bool p_mute = false,
        const Eigen::MatrixXd& UV_bnds = Eigen::MatrixXd(),
        const Eigen::MatrixXi& E = Eigen::MatrixXi(),
        const Eigen::VectorXi& bnd = Eigen::VectorXi(),
        const Config& animConfig = Config());
    virtual ~Optimizer(void);

public: // API
    virtual void setTime(double duration, double dt);
    virtual void updateTargetGRes(void);
    // precompute preconditioning matrix and factorize for fast solve, prepare initial guess
    virtual void precompute(void);

    // solve the optimization problem that minimizes E using a hill-climbing method,
    // the final result will be in result
    virtual int solve(int maxIter = 100);

    virtual void updatePrecondMtrAndFactorize(void);

    virtual void getFaceFieldForVis(Eigen::VectorXd& field);
    virtual Mesh<dim>& getResult(void);
    virtual int getIterNum(void) const;
    virtual int getFrameAmt(void) const;
    virtual int getInnerIterAmt(void) const;
    virtual void setRelGL2Tol(double p_relTol = 1.0e-2);
    virtual double getDt(void) const;
    virtual const AnimScripter<dim>& getAnimScripter() const;
    virtual void setAnimScriptType(AnimScriptType animScriptType,
        const std::string& meshSeqFolderPath);

    virtual void saveStatus(const std::string& appendStr = "");
    virtual void outputCollStats(std::ostream& out);
#ifdef EXPORT_FRICTION_DATA
    virtual void saveFrictionData(const Mesh<dim>& mesh) const;
#endif

    virtual void checkGradient(void);
    virtual void checkHessian(void);

protected: // helper functions
    /// @breif Compute the intial values for P_QP
    virtual void precomputeQPObjective(
        const LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* linSys,
        Eigen::SparseMatrix<double>& P,
        std::vector<double*>& elemPtr_P) const;
    /// @brief Update the elements of P_QP
    virtual void updateQPObjective(
        const LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* linSys,
        const std::vector<double*>& elemPtr_P) const;
    /// @brief Clear the next active set.
    virtual void clearActiveSet();
    /// @brief Update the active set used by the SQP and QP solvers.
    virtual bool updateActiveSet_QP();
    /// @brief Count the number of constraints in the active set.
    virtual int countConstraints();
    /// @brief Compute the QP inequality constraints \f$Ax \geq b\f$
    virtual void computeQPInequalityConstraint(
        const Mesh<dim>& mesh,
        const std::vector<std::shared_ptr<CollisionObject<dim>>>& collisionObjects,
        const std::vector<std::vector<int>>& activeSet,
        const int num_vars,
        std::vector<int>& constraintStartInds,
        Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b) const;

    /// @brief Solve the QP of the objective and collision constraints.
    virtual bool solveQP(
        const Mesh<dim>& mesh,
        const std::vector<std::shared_ptr<CollisionObject<dim>>>& collisionObjects,
        const std::vector<std::vector<int>>& activeSet,
        const LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* linSys,
        Eigen::SparseMatrix<double>& P,
        const std::vector<double*>& elemPtr_P,
        Eigen::VectorXd& gradient,
        OSQP& OSQPSolver,
#ifdef USE_GUROBI
        Eigen::GurobiSparse& gurobiQPSolver,
#endif
        std::vector<int>& constraintStartInds,
        Eigen::VectorXd& searchDir,
        Eigen::VectorXd& dual,
        QPSolverType qpSolverType) const;
    /// @brief Solve the QP using OSQP.
    virtual bool solveQP_OSQP(
        Eigen::SparseMatrix<double>& Q, Eigen::VectorXd& c,
        Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b,
        OSQP& QPSolver, Eigen::VectorXd& x, Eigen::VectorXd& dual) const;
#ifdef USE_GUROBI
    /// @brief Solve the QP using Gurobi.
    virtual bool solveQP_Gurobi(
        const Eigen::SparseMatrix<double>& Q, const Eigen::VectorXd& c,
        const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b,
        Eigen::GurobiSparse& QPSolver,
        Eigen::VectorXd& x, Eigen::VectorXd& dual) const;
#endif

    virtual void computeQPResidual(const Mesh<dim>& mesh,
        const std::vector<std::shared_ptr<CollisionObject<dim>>>& collisionObjects,
        const std::vector<std::vector<int>>& activeSet,
        const std::vector<int>& constraintStartInds,
        const Eigen::VectorXd& gradient,
        const Eigen::VectorXd& dual,
        Eigen::VectorXd& grad_KKT,
        Eigen::VectorXd& constraintVal,
        Eigen::VectorXd& fb);

    virtual void initX(int option, std::vector<std::vector<int>>& p_activeSet_next);
    virtual void computeXTilta(void);

    virtual bool fullyImplicit(void);
    virtual bool fullyImplicit_IP(void);

    virtual bool solveSub_IP(double kappa, std::vector<std::vector<int>>& AHat,
        std::vector<std::vector<MMCVID>>& MMAHat);
    virtual void initKappa(double& kappa);
    virtual void upperBoundKappa(double& kappa);
    virtual void suggestKappa(double& kappa);

    virtual void initSubProb_IP(void);
    virtual void computeSearchDir(int k, bool projectDBC = true);
    virtual void postLineSearch(double alpha);

    // solve for new configuration in the next iteration
    //NOTE: must compute current gradient first
    virtual bool solve_oneStep(void);

    virtual void computeConstraintSets(const Mesh<dim>& data, bool rehash = true);
    virtual void buildConstraintStartInds(const std::vector<std::vector<int>>& activeSet,
        std::vector<int>& constraintStartInds);
    virtual void buildConstraintStartIndsWithMM(const std::vector<std::vector<int>>& activeSet,
        const std::vector<std::vector<MMCVID>>& MMActiveSet,
        std::vector<int>& constraintStartInds);

    virtual bool lineSearch(double& stepSize, double armijoParam = 0.0, double lowerBound = 0.0);

    virtual void stepForward(const Eigen::MatrixXd& dataV0,
        Mesh<dim>& data, double stepSize) const;

    virtual void computeEnergyVal(const Mesh<dim>& data,
        int redoSVD, double& energyVal);
    virtual void computeGradient(const Mesh<dim>& data,
        bool redoSVD, Eigen::VectorXd& gradient, bool projectDBC = true);
    virtual void computePrecondMtr(const Mesh<dim>& data,
        bool redoSVD,
        LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* p_linSysSolver,
        bool updateDamping = false, bool projectDBC = true);

    virtual void computeDampingMtr(const Mesh<dim>& data,
        bool redoSVD,
        LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* p_dampingMtr,
        bool projectDBC = true);
    virtual void setupDampingMtr(const Mesh<dim>& data,
        bool redoSVD,
        LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>* p_dampingMtr);

    virtual void computeSystemEnergy(std::vector<double>& sysE,
        std::vector<Eigen::Matrix<double, 1, dim>>& sysM,
        std::vector<Eigen::Matrix<double, 1, dim>>& sysL);

    virtual void initStepSize(const Mesh<dim>& data, double& stepSize);
};

} // namespace IPC

#endif /* Optimizer_hpp */
