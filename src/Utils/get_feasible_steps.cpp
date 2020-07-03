#include "get_feasible_steps.hpp"

#include <Eigen/SparseCore>
#include <complex>

using namespace Eigen;
using namespace std;

double getSmallestPositiveRealQuadRoot(double a, double b, double c,
    double tol)
{
    // return negative value if no positive real root is found
    double t;

    if (abs(a) <= tol)
        t = -c / b;
    else {
        double desc = b * b - 4 * a * c;
        if (desc > 0) {
            t = (-b - sqrt(desc)) / (2 * a);
            if (t < 0)
                t = (-b + sqrt(desc)) / (2 * a);
        }
        else // desv<0 ==> imag
            t = -1;
    }
    return t;
}

// MatrixXi F
void computeInjectiveStepSize_2d(const MatrixXi& F,
    const MatrixXd& x,
    const VectorXd& p,
    double tol,
    double slackness,
    double* output)
{
    int n_tri = F.rows();
    double x1, x2, x3, y1, y2, y3;
    double p1, p2, p3, q1, q2, q3;
    double a, b, c, t;

    for (int ii = 0; ii < n_tri; ii++) {
        x1 = x(F(ii, 0), 0);
        x2 = x(F(ii, 1), 0);
        x3 = x(F(ii, 2), 0);

        y1 = x(F(ii, 0), 1);
        y2 = x(F(ii, 1), 1);
        y3 = x(F(ii, 2), 1);

        p1 = p[F(ii, 0) * 2];
        p2 = p[F(ii, 1) * 2];
        p3 = p[F(ii, 2) * 2];

        q1 = p[F(ii, 0) * 2 + 1];
        q2 = p[F(ii, 1) * 2 + 1];
        q3 = p[F(ii, 2) * 2 + 1];

        a = p1 * q2 - p2 * q1 - p1 * q3 + p3 * q1 + p2 * q3 - p3 * q2;
        b = p1 * y2 - p2 * y1 - q1 * x2 + q2 * x1 - p1 * y3 + p3 * y1 + q1 * x3 - q3 * x1 + p2 * y3 - p3 * y2 - q2 * x3 + q3 * x2;
        c = (1.0 - slackness) * (x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2);

        t = getSmallestPositiveRealQuadRoot(a, b, c, tol);

        if (t >= 0) {
            output[ii] = t;
        }
        else {
            output[ii] = 1e20;
        }
    }
}

double getSmallestPositiveRealCubicRoot(double a, double b, double c, double d,
    double tol)
{
    // return negative value if no positive real root is found
    double t = -1;

    if (abs(a) <= tol)
        t = getSmallestPositiveRealQuadRoot(b, c, d, tol);
    else {
        complex<double> i(0, 1);
        complex<double> delta0(b * b - 3 * a * c, 0);
        complex<double> delta1(2 * b * b * b - 9 * a * b * c + 27 * a * a * d, 0);
        complex<double> C = pow((delta1 + sqrt(delta1 * delta1 - 4.0 * delta0 * delta0 * delta0)) / 2.0, 1.0 / 3.0);
        if (std::abs(C) == 0.0) {
            // a corner case listed by wikipedia found by our collaborate from another project
            C = pow((delta1 - sqrt(delta1 * delta1 - 4.0 * delta0 * delta0 * delta0)) / 2.0, 1.0 / 3.0);
        }

        complex<double> u2 = (-1.0 + sqrt(3.0) * i) / 2.0;
        complex<double> u3 = (-1.0 - sqrt(3.0) * i) / 2.0;

        complex<double> t1 = (b + C + delta0 / C) / (-3.0 * a);
        complex<double> t2 = (b + u2 * C + delta0 / (u2 * C)) / (-3.0 * a);
        complex<double> t3 = (b + u3 * C + delta0 / (u3 * C)) / (-3.0 * a);

        if ((abs(imag(t1)) < tol) && (real(t1) > 0))
            t = real(t1);
        if ((abs(imag(t2)) < tol) && (real(t2) > 0) && ((real(t2) < t) || (t < 0)))
            t = real(t2);
        if ((abs(imag(t3)) < tol) && (real(t3) > 0) && ((real(t3) < t) || (t < 0)))
            t = real(t3);
    }
    return t;
}

void computeInjectiveStepSize_3d(const MatrixXi& F,
    const MatrixXd& x,
    const VectorXd& p,
    double tol,
    double slackness,
    double* output)
{
    int n_tri = F.rows();
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
    double p1, p2, p3, p4, q1, q2, q3, q4, r1, r2, r3, r4;
    double a, b, c, d, t;

    for (int ii = 0; ii < n_tri; ii++) {
        x1 = x(F(ii, 0), 0);
        x2 = x(F(ii, 1), 0);
        x3 = x(F(ii, 2), 0);
        x4 = x(F(ii, 3), 0);

        y1 = x(F(ii, 0), 1);
        y2 = x(F(ii, 1), 1);
        y3 = x(F(ii, 2), 1);
        y4 = x(F(ii, 3), 1);

        z1 = x(F(ii, 0), 2);
        z2 = x(F(ii, 1), 2);
        z3 = x(F(ii, 2), 2);
        z4 = x(F(ii, 3), 2);

        int _3Fii0 = F(ii, 0) * 3;
        int _3Fii1 = F(ii, 1) * 3;
        int _3Fii2 = F(ii, 2) * 3;
        int _3Fii3 = F(ii, 3) * 3;

        p1 = p[_3Fii0];
        p2 = p[_3Fii1];
        p3 = p[_3Fii2];
        p4 = p[_3Fii3];

        q1 = p[_3Fii0 + 1];
        q2 = p[_3Fii1 + 1];
        q3 = p[_3Fii2 + 1];
        q4 = p[_3Fii3 + 1];

        r1 = p[_3Fii0 + 2];
        r2 = p[_3Fii1 + 2];
        r3 = p[_3Fii2 + 2];
        r4 = p[_3Fii3 + 2];

        a = -p1 * q2 * r3 + p1 * r2 * q3 + q1 * p2 * r3 - q1 * r2 * p3 - r1 * p2 * q3 + r1 * q2 * p3 + p1 * q2 * r4 - p1 * r2 * q4 - q1 * p2 * r4 + q1 * r2 * p4 + r1 * p2 * q4 - r1 * q2 * p4 - p1 * q3 * r4 + p1 * r3 * q4 + q1 * p3 * r4 - q1 * r3 * p4 - r1 * p3 * q4 + r1 * q3 * p4 + p2 * q3 * r4 - p2 * r3 * q4 - q2 * p3 * r4 + q2 * r3 * p4 + r2 * p3 * q4 - r2 * q3 * p4;
        b = -x1 * q2 * r3 + x1 * r2 * q3 + y1 * p2 * r3 - y1 * r2 * p3 - z1 * p2 * q3 + z1 * q2 * p3 + x2 * q1 * r3 - x2 * r1 * q3 - y2 * p1 * r3 + y2 * r1 * p3 + z2 * p1 * q3 - z2 * q1 * p3 - x3 * q1 * r2 + x3 * r1 * q2 + y3 * p1 * r2 - y3 * r1 * p2 - z3 * p1 * q2 + z3 * q1 * p2 + x1 * q2 * r4 - x1 * r2 * q4 - y1 * p2 * r4 + y1 * r2 * p4 + z1 * p2 * q4 - z1 * q2 * p4 - x2 * q1 * r4 + x2 * r1 * q4 + y2 * p1 * r4 - y2 * r1 * p4 - z2 * p1 * q4 + z2 * q1 * p4 + x4 * q1 * r2 - x4 * r1 * q2 - y4 * p1 * r2 + y4 * r1 * p2 + z4 * p1 * q2 - z4 * q1 * p2 - x1 * q3 * r4 + x1 * r3 * q4 + y1 * p3 * r4 - y1 * r3 * p4 - z1 * p3 * q4 + z1 * q3 * p4 + x3 * q1 * r4 - x3 * r1 * q4 - y3 * p1 * r4 + y3 * r1 * p4 + z3 * p1 * q4 - z3 * q1 * p4 - x4 * q1 * r3 + x4 * r1 * q3 + y4 * p1 * r3 - y4 * r1 * p3 - z4 * p1 * q3 + z4 * q1 * p3 + x2 * q3 * r4 - x2 * r3 * q4 - y2 * p3 * r4 + y2 * r3 * p4 + z2 * p3 * q4 - z2 * q3 * p4 - x3 * q2 * r4 + x3 * r2 * q4 + y3 * p2 * r4 - y3 * r2 * p4 - z3 * p2 * q4 + z3 * q2 * p4 + x4 * q2 * r3 - x4 * r2 * q3 - y4 * p2 * r3 + y4 * r2 * p3 + z4 * p2 * q3 - z4 * q2 * p3;
        c = -x1 * y2 * r3 + x1 * z2 * q3 + x1 * y3 * r2 - x1 * z3 * q2 + y1 * x2 * r3 - y1 * z2 * p3 - y1 * x3 * r2 + y1 * z3 * p2 - z1 * x2 * q3 + z1 * y2 * p3 + z1 * x3 * q2 - z1 * y3 * p2 - x2 * y3 * r1 + x2 * z3 * q1 + y2 * x3 * r1 - y2 * z3 * p1 - z2 * x3 * q1 + z2 * y3 * p1 + x1 * y2 * r4 - x1 * z2 * q4 - x1 * y4 * r2 + x1 * z4 * q2 - y1 * x2 * r4 + y1 * z2 * p4 + y1 * x4 * r2 - y1 * z4 * p2 + z1 * x2 * q4 - z1 * y2 * p4 - z1 * x4 * q2 + z1 * y4 * p2 + x2 * y4 * r1 - x2 * z4 * q1 - y2 * x4 * r1 + y2 * z4 * p1 + z2 * x4 * q1 - z2 * y4 * p1 - x1 * y3 * r4 + x1 * z3 * q4 + x1 * y4 * r3 - x1 * z4 * q3 + y1 * x3 * r4 - y1 * z3 * p4 - y1 * x4 * r3 + y1 * z4 * p3 - z1 * x3 * q4 + z1 * y3 * p4 + z1 * x4 * q3 - z1 * y4 * p3 - x3 * y4 * r1 + x3 * z4 * q1 + y3 * x4 * r1 - y3 * z4 * p1 - z3 * x4 * q1 + z3 * y4 * p1 + x2 * y3 * r4 - x2 * z3 * q4 - x2 * y4 * r3 + x2 * z4 * q3 - y2 * x3 * r4 + y2 * z3 * p4 + y2 * x4 * r3 - y2 * z4 * p3 + z2 * x3 * q4 - z2 * y3 * p4 - z2 * x4 * q3 + z2 * y4 * p3 + x3 * y4 * r2 - x3 * z4 * q2 - y3 * x4 * r2 + y3 * z4 * p2 + z3 * x4 * q2 - z3 * y4 * p2;
        d = (1.0 - slackness) * (x1 * z2 * y3 - x1 * y2 * z3 + y1 * x2 * z3 - y1 * z2 * x3 - z1 * x2 * y3 + z1 * y2 * x3 + x1 * y2 * z4 - x1 * z2 * y4 - y1 * x2 * z4 + y1 * z2 * x4 + z1 * x2 * y4 - z1 * y2 * x4 - x1 * y3 * z4 + x1 * z3 * y4 + y1 * x3 * z4 - y1 * z3 * x4 - z1 * x3 * y4 + z1 * y3 * x4 + x2 * y3 * z4 - x2 * z3 * y4 - y2 * x3 * z4 + y2 * z3 * x4 + z2 * x3 * y4 - z2 * y3 * x4);

        t = getSmallestPositiveRealCubicRoot(a, b, c, d, tol);

        if (t >= 0) {
            output[ii] = t;
        }
        else {
            output[ii] = 1e20;
        }
    }
}

void computeInjectiveStepSize_3d(const RowVector3d& xx1,
    const RowVector3d& xx2,
    const RowVector3d& xx3,
    const RowVector3d& xx4,
    const Vector3d& pp1,
    const Vector3d& pp2,
    const Vector3d& pp3,
    const Vector3d& pp4,
    double tol,
    double& output)
{
    Eigen::Vector3d v01 = xx2 - xx1;
    Eigen::Vector3d v02 = xx3 - xx1;
    Eigen::Vector3d v03 = xx4 - xx1;
    Eigen::Vector3d p01 = pp2 - pp1;
    Eigen::Vector3d p02 = pp3 - pp1;
    Eigen::Vector3d p03 = pp4 - pp1;

    Eigen::Vector3d p01xp02 = p01.cross(p02);
    Eigen::Vector3d v01xp02pp01xv02 = v01.cross(p02) + p01.cross(v02);
    Eigen::Vector3d v01xv02 = v01.cross(v02);

    double a = p03.dot(p01xp02);
    double b = v03.dot(p01xp02) + p03.dot(v01xp02pp01xv02);
    double c = p03.dot(v01xv02) + v03.dot(v01xp02pp01xv02);
    double d = v03.dot(v01xv02);

    double t = getSmallestPositiveRealCubicRoot(a, b, c, d, tol);

    if (t >= 0.0) {
        output = t;
    }
    else {
        output = 1e20;
    }
}
