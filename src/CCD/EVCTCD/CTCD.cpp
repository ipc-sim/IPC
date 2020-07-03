#include "CTCD.h"
#include <vector>
#include "rpoly.h"
#include <Eigen/Geometry>
#include <iostream>

using namespace Eigen;
using namespace std;

bool TimeInterval::overlap(const TimeInterval& t1, const TimeInterval& t2)
{
    return !(t1.l > t2.u || t2.l > t1.u);
}

bool TimeInterval::overlap(const std::vector<TimeInterval>& intervals)
{
    for (std::vector<TimeInterval>::const_iterator it1 = intervals.begin(); it1 != intervals.end(); ++it1) {
        std::vector<TimeInterval>::const_iterator it2 = it1;
        for (++it2; it2 != intervals.end(); ++it2)
            if (!overlap(*it1, *it2))
                return false;
    }
    return true;
}

TimeInterval TimeInterval::intersect(const std::vector<TimeInterval>& intervals)
{
    TimeInterval isect(0.0, 1.0);
    for (std::vector<TimeInterval>::const_iterator it = intervals.begin(); it != intervals.end(); ++it) {
        isect.l = max(it->l, isect.l);
        isect.u = min(it->u, isect.u);
    }
    return isect;
}

int CTCD::getQuadRoots(double a, double b, double c, double& t0, double& t1)
{
    int roots = 0;
    int sign = 1;

    if (b < 0)
        sign = -1;

    double D = b * b - 4 * a * c;
    if (D >= 0) {
        roots = 2;
        double q = -0.5 * (b + sign * sqrt(D));
        t0 = q / a;
        t1 = c / q;
        if (t0 > t1)
            std::swap(t1, t0);
    }
    return roots;
}

void CTCD::checkInterval(double t1, double t2, double* op, int degree, vector<TimeInterval>& intervals, bool pos)
{
    // clamp values
    t1 = max(0.0, t1);
    t2 = max(0.0, t2);
    t1 = min(1.0, t1);
    t2 = min(1.0, t2);

    double tmid = (t2 + t1) / 2;
    double f = op[0];
    for (int i = 1; i <= degree; i++) {
        f *= tmid;
        f += op[i];
    }

    if (pos && f >= 0)
        intervals.push_back(TimeInterval(t1, t2));
    else if (!pos && f <= 0)
        intervals.push_back(TimeInterval(t1, t2));
}

bool CTCD::couldHaveRoots(double* op, int degree, bool pos)
{
    double result = 0;
    if ((pos && op[0] > 0) || (!pos && op[0] < 0))
        result = op[0];
    for (int i = 1; i < degree; i++) {
        result *= 1.0;
        if ((pos && op[i] > 0) || (!pos && op[i] < 0))
            result += op[i];
    }
    result *= 1.0;
    result += op[degree];
    return !((pos && result < 0) || (!pos && result > 0));
}

void CTCD::findIntervals(double* op, int n, vector<TimeInterval>& intervals, bool pos)
{
    int roots = 0;
    int reducedDegree = n;

    if (n > 6) {
        assert(!"Polynomials of degree > 6 not supported");
        return;
    }

    double time[6];
    // We don't care one bit about these imaginary roots
    double zeroi[6];

    // normalize
    double maxval = 0;
    for (int i = 0; i <= n; i++)
        maxval = std::max(maxval, fabs(op[i]));
    if (maxval != 0)
        for (int i = 0; i <= n; i++)
            op[i] /= maxval;

    for (int i = 0; i < n; i++) {
        if (op[i] == 0)
            reducedDegree--;
        else
            break;
    }

    if (reducedDegree < n) {
        for (int i = 0; i <= reducedDegree; i++)
            op[i] = op[i + n - reducedDegree];
    }

    if (reducedDegree > 2) {
        if (!couldHaveRoots(op, reducedDegree, pos))
            return;

        RootFinder rf;

        roots = rf.rpoly(op, reducedDegree, time, zeroi);
    }
    else if (reducedDegree == 2) {
        roots = getQuadRoots(op[0], op[1], op[2], time[0], time[1]);
    }
    else if (reducedDegree == 1) {
        time[0] = -op[1] / op[0];
        roots = 1;
    }
    else {
        // both points stationary -- check if colliding at t=0
        if ((!pos && op[0] <= 0) || (pos && op[0] >= 0))
            intervals.push_back(TimeInterval(0, 1.0));
        return;
    }

    // check intervals
    if (roots > 0) {
        std::sort(time, time + roots);
        if (time[0] >= 0)
            checkInterval(0, time[0], op, reducedDegree, intervals, pos);
        for (int i = 0; i < roots - 1; i++) {
            if (!((time[i] < 0 && time[i + 1] < 0) || (time[i] > 1.0 && time[i + 1] > 1.0)))
                checkInterval(time[i], time[i + 1], op, reducedDegree, intervals, pos);
        }
        if (time[roots - 1] <= 1.0)
            checkInterval(time[roots - 1], 1.0, op, reducedDegree, intervals, pos);
    }
    else {
        checkInterval(0.0, 1.0, op, reducedDegree, intervals, pos);
    }
}

void CTCD::barycentricPoly3D(const Vector3d& x10,
    const Vector3d& x20,
    const Vector3d& x30,
    const Vector3d& v10,
    const Vector3d& v20,
    const Vector3d& v30,
    vector<TimeInterval>& result)
{
    // alpha > 0
    double A = x10.dot(x10);
    double B = 2 * x10.dot(v10);
    double C = (v10).dot(v10);
    // e0.e1
    double D = (x20).dot(x10);
    double E = (x20).dot(v10) + (v20).dot(x10);
    double F = (v20).dot(v10);
    //(q0-q1).e0
    double G = (x30).dot(x20);
    double H = (x30).dot(v20) + (v30).dot(x20);
    double I = (v30).dot(v20);
    //(q0-q1).e1
    double J = (x30).dot(x10);
    double K = (x30).dot(v10) + (v30).dot(x10);
    double L = (v30).dot(v10);

    double op[5];

    op[0] = F * L - C * I;
    op[1] = F * K + E * L - C * H - B * I;
    op[2] = F * J + D * L + E * K - C * G - A * I - B * H;
    op[3] = D * K + E * J - A * H - B * G;
    op[4] = D * J - A * G;

    findIntervals(op, 4, result, true);
}

void CTCD::planePoly3D(const Vector3d& x10,
    const Vector3d& x20,
    const Vector3d& x30,
    const Vector3d& v10,
    const Vector3d& v20,
    const Vector3d& v30,
    vector<TimeInterval>& result)
{
    double op[4];
    op[0] = v10.dot(v20.cross(v30));
    op[1] = x10.dot(v20.cross(v30)) + v10.dot(x20.cross(v30)) + v10.dot(v20.cross(x30));
    op[2] = x10.dot(x20.cross(v30)) + x10.dot(v20.cross(x30)) + v10.dot(x20.cross(x30));
    op[3] = x10.dot(x20.cross(x30));
    findIntervals(op, 3, result, true);
}

void CTCD::distancePoly3D(const Vector3d& x10,
    const Vector3d& x20,
    const Vector3d& x30,
    const Vector3d& v10,
    const Vector3d& v20,
    const Vector3d& v30,
    double minDSquared,
    vector<TimeInterval>& result)
{
    double A = v10.dot(v20.cross(v30));
    double B = x10.dot(v20.cross(v30)) + v10.dot(x20.cross(v30)) + v10.dot(v20.cross(x30));
    double C = x10.dot(x20.cross(v30)) + x10.dot(v20.cross(x30)) + v10.dot(x20.cross(x30));
    double D = x10.dot(x20.cross(x30));
    Vector3d E = x20.cross(x30);
    Vector3d F = x20.cross(v30) + v20.cross(x30);
    Vector3d G = v20.cross(v30);

    double op[7];
    op[0] = A * A;
    op[1] = 2 * A * B;
    op[2] = B * B + 2 * A * C - G.dot(G) * minDSquared;
    op[3] = 2 * A * D + 2 * B * C - 2 * G.dot(F) * minDSquared;
    op[4] = 2 * B * D + C * C - (2 * G.dot(E) + F.dot(F)) * minDSquared;
    op[5] = 2 * C * D - 2 * F.dot(E) * minDSquared;
    op[6] = D * D - E.dot(E) * minDSquared;
    findIntervals(op, 6, result, false);
}

bool CTCD::edgeEdgeCTCD(const Eigen::Vector3d& q0start,
    const Eigen::Vector3d& p0start,
    const Eigen::Vector3d& q1start,
    const Eigen::Vector3d& p1start,
    const Eigen::Vector3d& q0end,
    const Eigen::Vector3d& p0end,
    const Eigen::Vector3d& q1end,
    const Eigen::Vector3d& p1end,
    double eta,
    double& t)
{
    double minD = eta * eta;

    // time intervals during which v is colinear with the edge, on the side of e1 towards e2, and on the side of e2 towards e1
    std::vector<TimeInterval> rawcoplane, a0, a1, b0, b1;

    Vector3d x10 = p0start - p1start;
    Vector3d x20 = p0start - q0start;
    Vector3d x30 = p1start - q1start;

    Vector3d vp0 = p0end - p0start;
    Vector3d vp1 = p1end - p1start;
    Vector3d vq0 = q0end - q0start;
    Vector3d vq1 = q1end - q1start;

    Vector3d v10 = vp0 - vp1;
    Vector3d v20 = vp0 - vq0;
    Vector3d v30 = vp1 - vq1;

    distancePoly3D(x10, x20, x30, v10, v20, v30, minD, rawcoplane);

    // check for parallel edges
    std::vector<TimeInterval> coplane;
    std::vector<TimeInterval> parallel;

    for (size_t i = 0; i < rawcoplane.size(); i++) {
        double midt = (rawcoplane[i].u + rawcoplane[i].l) / 2;
        x10 = (q0start - p0start) + midt * (vq0 - vp0);
        x20 = (q1start - p1start) + midt * (vq1 - vp1);

        if (x10.cross(x20).norm() < 1e-8) {
            parallel.push_back(rawcoplane[i]);
        }
        else {
            coplane.push_back(rawcoplane[i]);
        }
    }

    if (coplane.empty())
        return false;

    x10 = p1start - q1start;
    v10 = vp1 - vq1;
    x20 = p0start - q0start;
    v20 = vp0 - vq0;
    x30 = q0start - q1start;
    v30 = vq0 - vq1;
    barycentricPoly3D(x10, x20, x30, v10, v20, v30, a0);
    if (a0.empty())
        return false;

    x20 = q0start - p0start;
    v20 = vq0 - vp0;
    x30 = p0start - q1start;
    v30 = vp0 - vq1;
    barycentricPoly3D(x10, x20, x30, v10, v20, v30, a1);
    if (a1.empty())
        return false;

    x10 = p0start - q0start;
    v10 = vp0 - vq0;
    x20 = p1start - q1start;
    v20 = vp1 - vq1;
    x30 = q1start - q0start;
    v30 = vq1 - vq0;
    barycentricPoly3D(x10, x20, x30, v10, v20, v30, b0);
    if (b0.empty())
        return false;

    //x10 = p0 - q0;
    //v10 = vp0 - vq0;
    x20 = q1start - p1start;
    v20 = vq1 - vp1;
    x30 = p1start - q0start;
    v30 = vp1 - vq0;

    barycentricPoly3D(x10, x20, x30, v10, v20, v30, b1);
    if (b1.empty())
        return false;

    // check intervals for overlap
    bool col = false;
    double mint = 1.0;
    for (int i = 0; i < (int)coplane.size(); i++) {
        for (int j = 0; j < (int)a0.size(); j++) {
            for (int k = 0; k < (int)a1.size(); k++) {
                for (int l = 0; l < (int)b0.size(); l++) {
                    for (int m = 0; m < (int)b1.size(); m++) {
                        vector<TimeInterval> intervals;
                        intervals.push_back(coplane[i]);
                        intervals.push_back(a0[j]);
                        intervals.push_back(a1[k]);
                        intervals.push_back(b0[l]);
                        intervals.push_back(b1[m]);
                        if (TimeInterval::overlap(intervals)) {
                            TimeInterval isect = TimeInterval::intersect(intervals);
                            bool skip = false;
                            for (int p = 0; p < (int)parallel.size(); p++) {
                                vector<TimeInterval> pcheck;
                                pcheck.push_back(isect);
                                pcheck.push_back(parallel[p]);
                                if (TimeInterval::overlap(pcheck)) {
                                    skip = true;
                                    break;
                                }
                            }
                            if (!skip) {
                                mint = min(mint, isect.l);
                                col = true;
                            }
                        }
                    }
                }
            }
        }
    }

    // handle parallel edges
    for (int i = 0; i < (int)parallel.size(); i++) {
    }

    if (col) {
        t = mint;
        return true;
    }
    return false;
}

bool CTCD::vertexFaceCTCD(const Vector3d& q0start,
    const Vector3d& q1start,
    const Vector3d& q2start,
    const Vector3d& q3start,
    const Vector3d& q0end,
    const Vector3d& q1end,
    const Vector3d& q2end,
    const Vector3d& q3end,
    double eta, double& t)
{
    double minD = eta * eta;
    Vector3d v0 = q0end - q0start;
    Vector3d v1 = q1end - q1start;
    Vector3d v2 = q2end - q2start;
    Vector3d v3 = q3end - q3start;

    // time intervals during which v is colinear with the edge, on the side of e1 towards e2, and on the side of e2 towards e1
    vector<TimeInterval> coplane, e1, e2, e3;

    // check p.((axb)xb)
    Vector3d x10 = q0start - q1start;
    Vector3d v10 = v0 - v1;
    Vector3d x20 = (q3start - q1start).cross(q2start - q1start);
    Vector3d v20 = (v3 - v1).cross(v2 - v1);
    Vector3d x30 = q3start - q1start;
    Vector3d v30 = v3 - v1;
    planePoly3D(x10, x20, x30, v10, v20, v30, e1);

    if (e1.empty())
        return false;

    x10 = q0start - q2start;
    v10 = v0 - v2;
    x20 = (q1start - q2start).cross(q3start - q2start);
    v20 = (v1 - v2).cross(v3 - v2);
    x30 = q1start - q2start;
    v30 = v1 - v2;
    planePoly3D(x10, x20, x30, v10, v20, v30, e2);

    if (e2.empty())
        return false;

    x10 = q0start - q3start;
    v10 = v0 - v3;
    x20 = (q2start - q3start).cross(q1start - q3start);
    v20 = (v2 - v3).cross(v1 - v3);
    x30 = q2start - q3start;
    v30 = v2 - v3;
    planePoly3D(x10, x20, x30, v10, v20, v30, e3);

    if (e3.empty())
        return false;

    x10 = q0start - q1start;
    x20 = q2start - q1start;
    x30 = q3start - q1start;
    v10 = v0 - v1;
    v20 = v2 - v1;
    v30 = v3 - v1;
    distancePoly3D(x10, x20, x30, v10, v20, v30, minD, coplane);

    if (coplane.empty())
        return false;

    bool col = false;
    double mint = 1.0;
    for (int i = 0; i < (int)coplane.size(); i++) {
        for (int j = 0; j < (int)e1.size(); j++) {
            for (int k = 0; k < (int)e2.size(); k++) {
                for (int l = 0; l < (int)e3.size(); l++) {
                    vector<TimeInterval> intervals;
                    intervals.push_back(coplane[i]);
                    intervals.push_back(e1[j]);
                    intervals.push_back(e2[k]);
                    intervals.push_back(e3[l]);
                    if (TimeInterval::overlap(intervals)) {
                        mint = std::min(TimeInterval::intersect(intervals).l, mint);
                        col = true;
                    }
                }
            }
        }
    }

    if (col) {
        t = mint;
        return true;
    }
    return false;
}

bool CTCD::vertexEdgeCTCD(const Vector3d& q0start,
    const Vector3d& q1start,
    const Vector3d& q2start,
    const Vector3d& q0end,
    const Vector3d& q1end,
    const Vector3d& q2end,
    double eta,
    double& t)
{
    double op[5];
    double minD = eta * eta;
    Vector3d v0 = q0end - q0start;
    Vector3d v1 = q1end - q1start;
    Vector3d v2 = q2end - q2start;

    // time intervals during which v is colinear with the edge, on the side of e1 towards e2, and on the side of e2 towards e1
    vector<TimeInterval> colin, e1, e2;

    Vector3d ab = q2start - q1start;
    Vector3d ac = q0start - q1start;
    Vector3d cb = q2start - q0start;
    Vector3d vab = v2 - v1;
    Vector3d vac = v0 - v1;
    Vector3d vcb = v2 - v0;

    double c = ab.dot(ac);
    double b = ac.dot(vab) + ab.dot(vac);
    double a = vab.dot(vac);
    op[0] = a;
    op[1] = b;
    op[2] = c;
    findIntervals(op, 2, e1, true);
    if (e1.empty())
        return false;

    c = ab.dot(cb);
    b = cb.dot(vab) + ab.dot(vcb);
    a = vab.dot(vcb);

    op[0] = a;
    op[1] = b;
    op[2] = c;
    findIntervals(op, 2, e2, true);
    if (e2.empty())
        return false;

    double A = ab.dot(ab);
    double B = 2 * ab.dot(vab);
    double C = vab.dot(vab);
    double D = ac.dot(ac);
    double E = 2 * ac.dot(vac);
    double F = vac.dot(vac);
    double G = ac.dot(ab);
    double H = vab.dot(ac) + vac.dot(ab);
    double I = vab.dot(vac);
    op[4] = A * D - G * G - minD * A;
    op[3] = B * D + A * E - 2 * G * H - minD * B;
    op[2] = B * E + A * F + C * D - H * H - 2 * G * I - minD * C;
    op[1] = B * F + C * E - 2 * H * I;
    op[0] = C * F - I * I;
    findIntervals(op, 4, colin, false);
    if (colin.empty())
        return false;

    double mint = 1.0;
    bool col = false;
    for (int i = 0; i < (int)colin.size(); i++) {
        for (int j = 0; j < (int)e1.size(); j++) {
            for (int k = 0; k < (int)e2.size(); k++) {
                vector<TimeInterval> intervals;
                intervals.push_back(colin[i]);
                intervals.push_back(e1[j]);
                intervals.push_back(e2[k]);
                if (TimeInterval::overlap(intervals)) {
                    mint = std::min(TimeInterval::intersect(intervals).l, mint);
                    col = true;
                }
            }
        }
    }

    if (col) {
        t = mint;
        return true;
    }
    return false;
}

bool CTCD::vertexVertexCTCD(const Vector3d& q1start,
    const Vector3d& q2start,
    const Vector3d& q1end,
    const Vector3d& q2end,
    double eta, double& t)
{
    int roots = 0;
    double min_d = eta * eta;
    double t1 = 0, t2 = 0;
    Vector3d v1 = q1end - q1start;
    Vector3d v2 = q2end - q2start;

    // t^2 term
    double a = v1.dot(v1) + v2.dot(v2) - 2 * v1.dot(v2);
    // t term
    double b = 2 * (v1.dot(q1start) - v2.dot(q1start) - v1.dot(q2start) + v2.dot(q2start));
    // current distance - min_d
    double c = q1start.dot(q1start) + q2start.dot(q2start) - 2 * q1start.dot(q2start) - min_d;
    if (a != 0) {
        roots = getQuadRoots(a, b, c, t1, t2);
    }
    else if (b != 0) {
        t1 = -c / b;
        roots = 1;
    }
    else {
        if (c <= 0) {
            t = 0;
            return true;
        }
        return false;
    }

    double op[3];
    op[0] = a;
    op[1] = b;
    op[2] = c;
    vector<TimeInterval> interval;

    if (roots == 2) {
        checkInterval(0, t1, op, 2, interval, false);
        if (!interval.empty()) {
            t = 0;
            return true;
        }
        checkInterval(t1, t2, op, 2, interval, false);
        if (!interval.empty()) {
            t = t1;
            return true;
        }
        checkInterval(t2, 1.0, op, 2, interval, false);
        if (!interval.empty()) {
            t = t2;
            return true;
        }
        return false;
    }
    else if (roots == 1) {
        checkInterval(0, t1, op, 2, interval, false);
        if (!interval.empty()) {
            t = 0;
            return true;
        }
        checkInterval(t1, 1.0, op, 2, interval, false);
        if (!interval.empty()) {
            t = t1;
            return true;
        }
        return false;
    }
    checkInterval(0, 1.0, op, 2, interval, false);
    if (!interval.empty()) {
        t = 0;
        return true;
    }
    return false;
}
