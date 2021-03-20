#ifndef get_feasible_steps_hpp
#define get_feasible_steps_hpp

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

double getSmallestPositiveRealQuadRoot(double a, double b, double c,
    double tol);

void computeInjectiveStepSize_2d(const MatrixXi& F,
    const MatrixXd& x,
    const VectorXd& p,
    double tol,
    double slackness,
    double* output);

double getSmallestPositiveRealCubicRoot(double a, double b, double c, double d,
    double tol);

void computeInjectiveStepSize_3d(const MatrixXi& F,
    const MatrixXd& x,
    const VectorXd& p,
    double tol,
    double slackness,
    double* output);

void computeInjectiveStepSize_3d(const RowVector3d& xx1,
    const RowVector3d& xx2,
    const RowVector3d& xx3,
    const RowVector3d& xx4,
    const Vector3d& pp1,
    const Vector3d& pp2,
    const Vector3d& pp3,
    const Vector3d& pp4,
    double tol,
    double& output);

#endif
