#ifndef CTCD_H
#define CTCD_H

#include <Eigen/Core>
#include <vector>

struct TimeInterval {
    TimeInterval(double tl, double tu)
        : l(tl), u(tu)
    {
        if (l > u) std::swap(l, u);
        l = std::max(l, 0.0);
        u = std::min(u, 1.0);
    }

    TimeInterval()
        : l(0), u(0) {}

    // Returns whether or not the intersection of the intervals is nonempty
    static bool overlap(const TimeInterval& t1, const TimeInterval& t2);
    static bool overlap(const std::vector<TimeInterval>& intervals);

    // Returns the intersection of the intervals **asuming the intersection is nonempty**
    static TimeInterval intersect(const std::vector<TimeInterval>& intervals);

    double l, u;
};

class CTCD {
public:
    // Looks for collisions between edges (q0start, p0start) and (q1start, p1start) as they move towards
    // (q0end, p0end) and (q1end, p1end). Returns true if the edges ever come closer than a distance eta to each
    // other, and stores the earliest time (in the interval [0,1]) at which they do so in t.
    // WARNING: Does not work correctly if two edges are parallel at the time of intersection -- vertexEdgeCTCD
    // should catch this case, though.
    static bool edgeEdgeCTCD(const Eigen::Vector3d& q0start,
        const Eigen::Vector3d& p0start,
        const Eigen::Vector3d& q1start,
        const Eigen::Vector3d& p1start,
        const Eigen::Vector3d& q0end,
        const Eigen::Vector3d& p0end,
        const Eigen::Vector3d& q1end,
        const Eigen::Vector3d& p1end, double eta,
        double& t);

    // Looks for collisions between the vertex q0start and the face (q1start, q2start, q3start) as they move
    // towards q0end and (q1end, q2end, q3end). Returns true if the vertex and face ever come closer than a distance
    // eta to each, and stores the earliest time (in the interval [0,1]) at which they do so in t.
    static bool vertexFaceCTCD(const Eigen::Vector3d& q0start,
        const Eigen::Vector3d& q1start,
        const Eigen::Vector3d& q2start,
        const Eigen::Vector3d& q3start,
        const Eigen::Vector3d& q0end,
        const Eigen::Vector3d& q1end,
        const Eigen::Vector3d& q2end,
        const Eigen::Vector3d& q3end,
        double eta,
        double& t);

    // Looks for the degenerate case of collisions between the vertex q0start and the edge (q1start, s2start) as they
    // move towards q0end and (q1end, q2end). Returns true if the vertex and edge ever come closer than a distance
    // eta to each other, and stores the earliest time (in the interval [0,1]) at which they do so in t.
    static bool vertexEdgeCTCD(const Eigen::Vector3d& q0start,
        const Eigen::Vector3d& q1start,
        const Eigen::Vector3d& q2start,
        const Eigen::Vector3d& q0end,
        const Eigen::Vector3d& q1end,
        const Eigen::Vector3d& q2end,
        double eta,
        double& t);

    // Looks for the degenerate case of collisions between the vertices q1start and q2start, as they move towards
    // q1end and q2end. Returns true if the vertices ever come closer than a distane of eta to each other, and stores
    // the earliest time (in the interval [0,1]) at which they do so in t.
    static bool vertexVertexCTCD(const Eigen::Vector3d& q1start,
        const Eigen::Vector3d& q2start,
        const Eigen::Vector3d& q1end,
        const Eigen::Vector3d& q2end,
        double eta, double& t);

private:
    // Solves the quadratic equation ax^2 + bx + c = 0, and puts the roots in t0, t1, in ascending order.
    // Returns the number of real roots found.
    static int getQuadRoots(double a, double b, double c, double& t0, double& t1);

    // Conservatively checks if the polynomial of degree degree, with coefficients op, could be positive (if pos is true)
    // or negative (if pos is negative) on the interval [0,1].
    static bool couldHaveRoots(double* op, int degree, bool pos);

    // Looks at the interval [t1, t2], on which a polynomial of degree degree and coefficients op is assumed to have
    // constant sign, and determines if the polynomial is all positive or all negative on that interval.
    // If positive, and pos is true, or if negative, and pos is false, clamps the interval to [0,1] and adds it to
    // intervals.
    static void checkInterval(double t1, double t2, double* op, int degree, std::vector<TimeInterval>& intervals, bool pos);

    // Computes the intervals of x in [0,1] where the polynomial of degree n, with coefficients in op
    // (given in "natural," descending order of power of x) is positive (when pos = true) or negative (if pos = false).
    static void findIntervals(double* op, int n, std::vector<TimeInterval>& intervals, bool pos);

    static void distancePoly3D(const Eigen::Vector3d& x10,
        const Eigen::Vector3d& x20,
        const Eigen::Vector3d& x30,
        const Eigen::Vector3d& v10,
        const Eigen::Vector3d& v20,
        const Eigen::Vector3d& v30,
        double minDSquared,
        std::vector<TimeInterval>& result);

    static void barycentricPoly3D(const Eigen::Vector3d& x10,
        const Eigen::Vector3d& x20,
        const Eigen::Vector3d& x30,
        const Eigen::Vector3d& v10,
        const Eigen::Vector3d& v20,
        const Eigen::Vector3d& v30,
        std::vector<TimeInterval>& result);

    static void planePoly3D(const Eigen::Vector3d& x10,
        const Eigen::Vector3d& x20,
        const Eigen::Vector3d& x30,
        const Eigen::Vector3d& v10,
        const Eigen::Vector3d& v20,
        const Eigen::Vector3d& v30,
        std::vector<TimeInterval>& result);
};

#endif
