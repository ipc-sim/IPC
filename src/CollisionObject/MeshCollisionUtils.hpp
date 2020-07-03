//
//  MeshCollisionUtils.hpp
//  IPC
//
//  Created by Minchen Li on 6/03/19.
//

#ifndef MeshCollisionUtils_hpp
#define MeshCollisionUtils_hpp

#include "Types.hpp"

#include "Mesh.hpp"

#include <iostream>
#include <array>

#include <boost/functional/hash.hpp>

#include <Eigen/Eigen>

namespace IPC {

class MMCVID {
public:
    std::array<int, 4> data;

public:
    MMCVID(int a, int b, int c, int d)
    {
        data[0] = a;
        data[1] = b;
        data[2] = c;
        data[3] = d;
    }
    MMCVID(int a = -1)
    {
        data[0] = a;
        data[1] = data[2] = data[3] = -1;
    }

    void set(int a, int b, int c, int d)
    {
        data[0] = a;
        data[1] = b;
        data[2] = c;
        data[3] = d;
    }

    const int& operator[](int i) const
    {
        return data[i];
    }
    int& operator[](int i)
    {
        return data[i];
    }
    friend bool operator==(const MMCVID& lhs, const MMCVID& rhs)
    {
        if ((lhs[0] == rhs[0]) && (lhs[1] == rhs[1]) && (lhs[2] == rhs[2]) && (lhs[3] == rhs[3])) {
            return true;
        }
        else {
            return false;
        }
    }
    friend bool operator<(const MMCVID& lhs, const MMCVID& rhs)
    {
        if (lhs[0] < rhs[0]) {
            return true;
        }
        else if (lhs[0] == rhs[0]) {
            if (lhs[1] < rhs[1]) {
                return true;
            }
            else if (lhs[1] == rhs[1]) {
                if (lhs[2] < rhs[2]) {
                    return true;
                }
                else if (lhs[2] == rhs[2]) {
                    if (lhs[3] < rhs[3]) {
                        return true;
                    }
                    else {
                        return false;
                    }
                }
                else {
                    return false;
                }
            }
            else {
                return false;
            }
        }
        else {
            return false;
        }
    }

    std::string str() const
    {
        std::stringstream s;
        s << "(";
        for (int i = 0; i < data.size(); i++) {
            s << this->data[i] << (i != data.size() - 1 ? "," : ")");
        }
        return s.str();
    }
};

class MMCVIDHash {
public:
    size_t operator()(const MMCVID& mmcvid) const
    {
        return boost::hash_range(mmcvid.data.begin(), mmcvid.data.end());
    }
};

template <int dim>
inline void compute_c(
    const Eigen::Matrix<double, 1, dim>& v0, // point
    const Eigen::Matrix<double, 1, dim>& v1, // triangle point 0
    const Eigen::Matrix<double, 1, dim>& v2, // triangle point 2
    const Eigen::Matrix<double, 1, dim>& v3, // triangle point 1
    double& c, double coef)
{
    if constexpr (dim == 3) {
        c = coef * (v3 - v0).dot((v1 - v0).cross(v2 - v0));
    }
    else {
    }
}
inline bool pointInTriangle(double x, double y)
{
    return (x >= 0.0) && (y >= 0.0) && (1.0 - x - y >= 0.0);
}
template <int dim>
inline bool pointInTriangleSweep(const Eigen::Matrix<double, 1, dim>& pt,
    const Eigen::Matrix<double, 1, dim>& v0,
    const Eigen::Matrix<double, 1, dim>& v1,
    const Eigen::Matrix<double, 1, dim>& v2)
{
    //TODO: consider directions!
    Eigen::Matrix<double, 2, dim> basis;
    basis.row(0) = v1 - v0;
    basis.row(1) = v2 - v0;

    Eigen::Matrix<double, 1, dim> vect = pt - v0;

    Eigen::Matrix<double, 2, 1> pt2D01 = (basis * basis.transpose()).ldlt().solve(basis * vect.transpose());

    return pointInTriangle(pt2D01[0], pt2D01[1]);
}

inline void d_PP(const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    double& d)
{
    d = (v0 - v1).squaredNorm();
}

inline void g_PP(const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    Eigen::Matrix<double, 6, 1>& g)
{
    g.template segment<3>(0) = 2.0 * (v0 - v1);
    g.template segment<3>(3) = -g.template segment<3>(0);
}

inline void H_PP(Eigen::Matrix<double, 6, 6>& H)
{
    H.setZero();
    H(0, 0) = H(1, 1) = H(2, 2) = H(3, 3) = H(4, 4) = H(5, 5) = 2.0;
    H(0, 3) = H(1, 4) = H(2, 5) = H(3, 0) = H(4, 1) = H(5, 2) = -2.0;
}

inline void derivTest_PP(void)
{
    Eigen::RowVector3d v[2];
    v[0] << 0.0, 0.0, 0.0;
    v[1] << 1.0, 0.1, 0.0;

    double d0;
    d_PP(v[0], v[1], d0);

    double eps = 1.0e-6;
    Eigen::Matrix<double, 6, 1> fd, s;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
            Eigen::RowVector3d vv[2];
            vv[0] = v[0];
            vv[1] = v[1];
            vv[i][j] += eps;
            double d;
            d_PP(vv[0], vv[1], d);
            fd[i * 3 + j] = (d - d0) / eps;
        }
    }

    g_PP(v[0], v[1], s);

    std::cout << s.transpose() << std::endl;
    std::cout << fd.transpose() << std::endl;
    std::cout << "relError = " << (s - fd).norm() / fd.norm() << std::endl;

    Eigen::Matrix<double, 6, 6> Hfd, Hs;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
            Eigen::RowVector3d vv[2];
            vv[0] = v[0];
            vv[1] = v[1];
            vv[i][j] += eps;
            Eigen::Matrix<double, 6, 1> g;
            g_PP(vv[0], vv[1], g);
            Hfd.col(i * 3 + j) = (g - s) / eps;
        }
    }

    H_PP(Hs);

    std::cout << Hs.transpose() << std::endl;
    std::cout << Hfd.transpose() << std::endl;
    std::cout << "relError = " << (Hs - Hfd).norm() / Hfd.norm() << std::endl;
}

inline void d_PE(const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    double& d)
{
    d = (v1 - v0).cross(v2 - v0).squaredNorm() / (v2 - v1).squaredNorm();
}

inline void g_PE(double v01, double v02, double v03, double v11, double v12, double v13,
    double v21, double v22, double v23, double g[9])
{
    double t17;
    double t18;
    double t19;
    double t20;
    double t21;
    double t22;
    double t23;
    double t24;
    double t25;
    double t42;
    double t44;
    double t45;
    double t46;
    double t43;
    double t50;
    double t51;
    double t52;

    /* G_PE */
    /*     G = G_PE(V01,V02,V03,V11,V12,V13,V21,V22,V23) */
    /*     This function was generated by the Symbolic Math Toolbox version 8.3. */
    /*     10-Jun-2019 18:02:37 */
    t17 = -v11 + v01;
    t18 = -v12 + v02;
    t19 = -v13 + v03;
    t20 = -v21 + v01;
    t21 = -v22 + v02;
    t22 = -v23 + v03;
    t23 = -v21 + v11;
    t24 = -v22 + v12;
    t25 = -v23 + v13;
    t42 = 1.0 / ((t23 * t23 + t24 * t24) + t25 * t25);
    t44 = t17 * t21 + -(t18 * t20);
    t45 = t17 * t22 + -(t19 * t20);
    t46 = t18 * t22 + -(t19 * t21);
    t43 = t42 * t42;
    t50 = (t44 * t44 + t45 * t45) + t46 * t46;
    t51 = (v11 * 2.0 + -(v21 * 2.0)) * t43 * t50;
    t52 = (v12 * 2.0 + -(v22 * 2.0)) * t43 * t50;
    t43 = (v13 * 2.0 + -(v23 * 2.0)) * t43 * t50;
    g[0] = t42 * (t24 * t44 * 2.0 + t25 * t45 * 2.0);
    g[1] = -t42 * (t23 * t44 * 2.0 - t25 * t46 * 2.0);
    g[2] = -t42 * (t23 * t45 * 2.0 + t24 * t46 * 2.0);
    g[3] = -t51 - t42 * (t21 * t44 * 2.0 + t22 * t45 * 2.0);
    g[4] = -t52 + t42 * (t20 * t44 * 2.0 - t22 * t46 * 2.0);
    g[5] = -t43 + t42 * (t20 * t45 * 2.0 + t21 * t46 * 2.0);
    g[6] = t51 + t42 * (t18 * t44 * 2.0 + t19 * t45 * 2.0);
    g[7] = t52 - t42 * (t17 * t44 * 2.0 - t19 * t46 * 2.0);
    g[8] = t43 - t42 * (t17 * t45 * 2.0 + t18 * t46 * 2.0);
}

inline void g_PE(const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    Eigen::Matrix<double, 9, 1>& g)
{
    g_PE(v0[0], v0[1], v0[2],
        v1[0], v1[1], v1[2],
        v2[0], v2[1], v2[2],
        g.data());
}

inline void H_PE(double v01, double v02, double v03, double v11, double v12, double v13,
    double v21, double v22, double v23, double H[81])
{
    double t17;
    double t18;
    double t19;
    double t20;
    double t21;
    double t22;
    double t23;
    double t24;
    double t25;
    double t26;
    double t27;
    double t28;
    double t35;
    double t36;
    double t37;
    double t50;
    double t51;
    double t52;
    double t53;
    double t54;
    double t55;
    double t56;
    double t62;
    double t70;
    double t71;
    double t75;
    double t79;
    double t80;
    double t84;
    double t88;
    double t38;
    double t39;
    double t40;
    double t41;
    double t42;
    double t43;
    double t44;
    double t46;
    double t48;
    double t57;
    double t58;
    double t60;
    double t63;
    double t65;
    double t67;
    double t102;
    double t103;
    double t104;
    double t162;
    double t163;
    double t164;
    double t213;
    double t214;
    double t215;
    double t216;
    double t217;
    double t218;
    double t225;
    double t226;
    double t227;
    double t229;
    double t230;
    double t311;
    double t231;
    double t232;
    double t233;
    double t234;
    double t235;
    double t236;
    double t237;
    double t238;
    double t239;
    double t240;
    double t245;
    double t279;
    double t281;
    double t282;
    double t283;
    double t287;
    double t289;
    double t247;
    double t248;
    double t249;
    double t250;
    double t251;
    double t252;
    double t253;
    double t293;
    double t295;
    double t299;
    double t300;
    double t303;
    double t304;
    double t294;
    double t297;
    double t301;
    double t302;

    /* H_PE */
    /*     H = H_PE(V01,V02,V03,V11,V12,V13,V21,V22,V23) */
    /*     This function was generated by the Symbolic Math Toolbox version 8.3. */
    /*     10-Jun-2019 18:02:39 */
    t17 = -v11 + v01;
    t18 = -v12 + v02;
    t19 = -v13 + v03;
    t20 = -v21 + v01;
    t21 = -v22 + v02;
    t22 = -v23 + v03;
    t23 = -v21 + v11;
    t24 = -v22 + v12;
    t25 = -v23 + v13;
    t26 = v11 * 2.0 + -(v21 * 2.0);
    t27 = v12 * 2.0 + -(v22 * 2.0);
    t28 = v13 * 2.0 + -(v23 * 2.0);
    t35 = t23 * t23;
    t36 = t24 * t24;
    t37 = t25 * t25;
    t50 = t17 * t21;
    t51 = t18 * t20;
    t52 = t17 * t22;
    t53 = t19 * t20;
    t54 = t18 * t22;
    t55 = t19 * t21;
    t56 = t17 * t20 * 2.0;
    t62 = t18 * t21 * 2.0;
    t70 = t19 * t22 * 2.0;
    t71 = t17 * t23 * 2.0;
    t75 = t18 * t24 * 2.0;
    t79 = t19 * t25 * 2.0;
    t80 = t20 * t23 * 2.0;
    t84 = t21 * t24 * 2.0;
    t88 = t22 * t25 * 2.0;
    t38 = t17 * t17 * 2.0;
    t39 = t18 * t18 * 2.0;
    t40 = t19 * t19 * 2.0;
    t41 = t20 * t20 * 2.0;
    t42 = t21 * t21 * 2.0;
    t43 = t22 * t22 * 2.0;
    t44 = t35 * 2.0;
    t46 = t36 * 2.0;
    t48 = t37 * 2.0;
    t57 = t50 * 2.0;
    t58 = t51 * 2.0;
    t60 = t52 * 2.0;
    t63 = t53 * 2.0;
    t65 = t54 * 2.0;
    t67 = t55 * 2.0;
    t102 = 1.0 / ((t35 + t36) + t37);
    t36 = t50 + -t51;
    t35 = t52 + -t53;
    t37 = t54 + -t55;
    t103 = t102 * t102;
    t104 = pow(t102, 3.0);
    t162 = -(t23 * t24 * t102 * 2.0);
    t163 = -(t23 * t25 * t102 * 2.0);
    t164 = -(t24 * t25 * t102 * 2.0);
    t213 = t18 * t36 * 2.0 + t19 * t35 * 2.0;
    t214 = t17 * t35 * 2.0 + t18 * t37 * 2.0;
    t215 = t21 * t36 * 2.0 + t22 * t35 * 2.0;
    t216 = t20 * t35 * 2.0 + t21 * t37 * 2.0;
    t217 = t24 * t36 * 2.0 + t25 * t35 * 2.0;
    t218 = t23 * t35 * 2.0 + t24 * t37 * 2.0;
    t35 = (t36 * t36 + t35 * t35) + t37 * t37;
    t225 = t17 * t36 * 2.0 + -(t19 * t37 * 2.0);
    t226 = t20 * t36 * 2.0 + -(t22 * t37 * 2.0);
    t227 = t23 * t36 * 2.0 + -(t25 * t37 * 2.0);
    t36 = t26 * t103;
    t229 = t36 * t213;
    t37 = t27 * t103;
    t230 = t37 * t213;
    t311 = t28 * t103;
    t231 = t311 * t213;
    t232 = t36 * t214;
    t233 = t37 * t214;
    t234 = t311 * t214;
    t235 = t36 * t215;
    t236 = t37 * t215;
    t237 = t311 * t215;
    t238 = t36 * t216;
    t239 = t37 * t216;
    t240 = t311 * t216;
    t214 = t36 * t217;
    t215 = t37 * t217;
    t216 = t311 * t217;
    t217 = t36 * t218;
    t245 = t37 * t218;
    t213 = t311 * t218;
    t279 = t103 * t35 * 2.0;
    t281 = t26 * t26 * t104 * t35 * 2.0;
    t282 = t27 * t27 * t104 * t35 * 2.0;
    t283 = t28 * t28 * t104 * t35 * 2.0;
    t287 = t26 * t27 * t104 * t35 * 2.0;
    t218 = t26 * t28 * t104 * t35 * 2.0;
    t289 = t27 * t28 * t104 * t35 * 2.0;
    t247 = t36 * t225;
    t248 = t37 * t225;
    t249 = t311 * t225;
    t250 = t36 * t226;
    t251 = t37 * t226;
    t252 = t311 * t226;
    t253 = t36 * t227;
    t35 = t37 * t227;
    t36 = t311 * t227;
    t293 = t102 * (t75 + t79) + t214;
    t295 = -(t102 * (t80 + t84)) + t213;
    t299 = t102 * ((t63 + t22 * t23 * 2.0) + -t60) + t217;
    t300 = t102 * ((t67 + t22 * t24 * 2.0) + -t65) + t245;
    t303 = -(t102 * ((t57 + t17 * t24 * 2.0) + -t58)) + t215;
    t304 = -(t102 * ((t60 + t17 * t25 * 2.0) + -t63)) + t216;
    t294 = t102 * (t71 + t75) + -t213;
    t297 = -(t102 * (t80 + t88)) + t35;
    t88 = -(t102 * (t84 + t88)) + -t214;
    t301 = t102 * ((t58 + t21 * t23 * 2.0) + -t57) + t253;
    t302 = t102 * ((t65 + t21 * t25 * 2.0) + -t67) + t36;
    t84 = t102 * ((t57 + t20 * t24 * 2.0) + -t58) + -t215;
    t80 = t102 * ((t60 + t20 * t25 * 2.0) + -t63) + -t216;
    t75 = -(t102 * ((t63 + t19 * t23 * 2.0) + -t60)) + -t217;
    t227 = -(t102 * ((t67 + t19 * t24 * 2.0) + -t65)) + -t245;
    t311 = ((-(t17 * t19 * t102 * 2.0) + t231) + -t232) + t218;
    t245 = ((-(t20 * t22 * t102 * 2.0) + t237) + -t238) + t218;
    t226 = ((-t102 * (t67 - t54 * 4.0) + t233) + t252) + -t289;
    t28 = ((-t102 * (t63 - t52 * 4.0) + t232) + -t237) + -t218;
    t27 = ((-t102 * (t58 - t50 * 4.0) + t247) + -t236) + -t287;
    t225 = ((-(t102 * (t65 + -(t55 * 4.0))) + t239) + t249) + -t289;
    t26 = ((-(t102 * (t60 + -(t53 * 4.0))) + t238) + -t231) + -t218;
    t103 = ((-(t102 * (t57 + -(t51 * 4.0))) + t250) + -t230) + -t287;
    t104 = (((-(t102 * (t56 + t62)) + t234) + t240) + t279) + -t283;
    t218 = (((-(t102 * (t56 + t70)) + t248) + t251) + t279) + -t282;
    t217 = (((-(t102 * (t62 + t70)) + -t229) + -t235) + t279) + -t281;
    t216 = t102 * (t71 + t79) + -t35;
    t215 = -(t102 * ((t58 + t18 * t23 * 2.0) + -t57)) + -t253;
    t214 = -(t102 * ((t65 + t18 * t25 * 2.0) + -t67)) + -t36;
    t213 = ((-(t17 * t18 * t102 * 2.0) + t230) + -t247) + t287;
    t37 = ((-(t20 * t21 * t102 * 2.0) + t236) + -t250) + t287;
    t36 = ((-(t18 * t19 * t102 * 2.0) + -t233) + -t249) + t289;
    t35 = ((-(t21 * t22 * t102 * 2.0) + -t239) + -t252) + t289;
    H[0] = t102 * (t46 + t48);
    H[1] = t162;
    H[2] = t163;
    H[3] = t88;
    H[4] = t84;
    H[5] = t80;
    H[6] = t293;
    H[7] = t303;
    H[8] = t304;
    H[9] = t162;
    H[10] = t102 * (t44 + t48);
    H[11] = t164;
    H[12] = t301;
    H[13] = t297;
    H[14] = t302;
    H[15] = t215;
    H[16] = t216;
    H[17] = t214;
    H[18] = t163;
    H[19] = t164;
    H[20] = t102 * (t44 + t46);
    H[21] = t299;
    H[22] = t300;
    H[23] = t295;
    H[24] = t75;
    H[25] = t227;
    H[26] = t294;
    H[27] = t88;
    H[28] = t301;
    H[29] = t299;
    H[30] = ((t235 * 2.0 + -t279) + t281) + t102 * (t42 + t43);
    H[31] = t37;
    H[32] = t245;
    H[33] = t217;
    H[34] = t27;
    H[35] = t28;
    H[36] = t84;
    H[37] = t297;
    H[38] = t300;
    H[39] = t37;
    H[40] = ((t251 * -2.0 + -t279) + t282) + t102 * (t41 + t43);
    H[41] = t35;
    H[42] = t103;
    H[43] = t218;
    H[44] = t226;
    H[45] = t80;
    H[46] = t302;
    H[47] = t295;
    H[48] = t245;
    H[49] = t35;
    H[50] = ((t240 * -2.0 + -t279) + t283) + t102 * (t41 + t42);
    H[51] = t26;
    H[52] = t225;
    H[53] = t104;
    H[54] = t293;
    H[55] = t215;
    H[56] = t75;
    H[57] = t217;
    H[58] = t103;
    H[59] = t26;
    H[60] = ((t229 * 2.0 + -t279) + t281) + t102 * (t39 + t40);
    H[61] = t213;
    H[62] = t311;
    H[63] = t303;
    H[64] = t216;
    H[65] = t227;
    H[66] = t27;
    H[67] = t218;
    H[68] = t225;
    H[69] = t213;
    H[70] = ((t248 * -2.0 + -t279) + t282) + t102 * (t38 + t40);
    H[71] = t36;
    H[72] = t304;
    H[73] = t214;
    H[74] = t294;
    H[75] = t28;
    H[76] = t226;
    H[77] = t104;
    H[78] = t311;
    H[79] = t36;
    H[80] = ((t234 * -2.0 + -t279) + t283) + t102 * (t38 + t39);
}

inline void H_PE(const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    Eigen::Matrix<double, 9, 9>& H)
{
    H_PE(v0[0], v0[1], v0[2],
        v1[0], v1[1], v1[2],
        v2[0], v2[1], v2[2],
        H.data());
}

inline void derivTest_PE(void)
{
    Eigen::RowVector3d v[3];
    v[0] << 0.0, 0.0, 0.0;
    v[1] << 1.0, 0.1, 0.0;
    v[2] << 0.0, 1.1, -0.1;

    double d0;
    d_PE(v[0], v[1], v[2], d0);

    double eps = 1.0e-6;
    Eigen::Matrix<double, 9, 1> fd, s;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Eigen::RowVector3d vv[3];
            vv[0] = v[0];
            vv[1] = v[1];
            vv[2] = v[2];
            vv[i][j] += eps;
            double d;
            d_PE(vv[0], vv[1], vv[2], d);
            fd[i * 3 + j] = (d - d0) / eps;
        }
    }

    g_PE(v[0], v[1], v[2], s);

    std::cout << s.transpose() << std::endl;
    std::cout << fd.transpose() << std::endl;
    std::cout << "relError = " << (s - fd).norm() / fd.norm() << std::endl;

    Eigen::Matrix<double, 9, 9> Hfd, Hs;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Eigen::RowVector3d vv[3];
            vv[0] = v[0];
            vv[1] = v[1];
            vv[2] = v[2];
            vv[i][j] += eps;
            Eigen::Matrix<double, 9, 1> g;
            g_PE(vv[0], vv[1], vv[2], g);
            Hfd.col(i * 3 + j) = (g - s) / eps;
        }
    }

    H_PE(v[0], v[1], v[2], Hs);

    std::cout << Hs.transpose() << std::endl;
    std::cout << Hfd.transpose() << std::endl;
    std::cout << "relError = " << (Hs - Hfd).norm() / Hfd.norm() << std::endl;
}

inline void d_PT(const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    const Eigen::RowVector3d& v3,
    double& d)
{
    Eigen::RowVector3d b = (v2 - v1).cross(v3 - v1);
    double aTb = (v0 - v1).dot(b);
    d = aTb * aTb / b.squaredNorm();
}

inline void g_PT(double v01, double v02, double v03, double v11, double v12, double v13,
    double v21, double v22, double v23, double v31, double v32, double v33,
    double g[12])
{
    double t11;
    double t12;
    double t13;
    double t14;
    double t15;
    double t16;
    double t17;
    double t18;
    double t19;
    double t20;
    double t21;
    double t22;
    double t32;
    double t33;
    double t34;
    double t43;
    double t45;
    double t44;
    double t46;

    /* G_PT */
    /*     G = G_PT(V01,V02,V03,V11,V12,V13,V21,V22,V23,V31,V32,V33) */
    /*     This function was generated by the Symbolic Math Toolbox version 8.3. */
    /*     10-Jun-2019 17:42:16 */
    t11 = -v11 + v01;
    t12 = -v12 + v02;
    t13 = -v13 + v03;
    t14 = -v21 + v11;
    t15 = -v22 + v12;
    t16 = -v23 + v13;
    t17 = -v31 + v11;
    t18 = -v32 + v12;
    t19 = -v33 + v13;
    t20 = -v31 + v21;
    t21 = -v32 + v22;
    t22 = -v33 + v23;
    t32 = t14 * t18 + -(t15 * t17);
    t33 = t14 * t19 + -(t16 * t17);
    t34 = t15 * t19 + -(t16 * t18);
    t43 = 1.0 / ((t32 * t32 + t33 * t33) + t34 * t34);
    t45 = (t13 * t32 + t11 * t34) + -(t12 * t33);
    t44 = t43 * t43;
    t46 = t45 * t45;
    g[0] = t34 * t43 * t45 * 2.0;
    g[1] = t33 * t43 * t45 * -2.0;
    g[2] = t32 * t43 * t45 * 2.0;
    t45 *= t43;
    g[3] = -t44 * t46 * (t21 * t32 * 2.0 + t22 * t33 * 2.0) - t45 * ((t34 + t12 * t22) - t13 * t21) * 2.0;
    t43 = t44 * t46;
    g[4] = t43 * (t20 * t32 * 2.0 - t22 * t34 * 2.0) + t45 * ((t33 + t11 * t22) - t13 * t20) * 2.0;
    g[5] = t43 * (t20 * t33 * 2.0 + t21 * t34 * 2.0) - t45 * ((t32 + t11 * t21) - t12 * t20) * 2.0;
    g[6] = t45 * (t12 * t19 - t13 * t18) * 2.0 + t43 * (t18 * t32 * 2.0 + t19 * t33 * 2.0);
    g[7] = t45 * (t11 * t19 - t13 * t17) * -2.0 - t43 * (t17 * t32 * 2.0 - t19 * t34 * 2.0);
    g[8] = t45 * (t11 * t18 - t12 * t17) * 2.0 - t43 * (t17 * t33 * 2.0 + t18 * t34 * 2.0);
    g[9] = t45 * (t12 * t16 - t13 * t15) * -2.0 - t43 * (t15 * t32 * 2.0 + t16 * t33 * 2.0);
    g[10] = t45 * (t11 * t16 - t13 * t14) * 2.0 + t43 * (t14 * t32 * 2.0 - t16 * t34 * 2.0);
    g[11] = t45 * (t11 * t15 - t12 * t14) * -2.0 + t43 * (t14 * t33 * 2.0 + t15 * t34 * 2.0);
}

inline void g_PT(const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    const Eigen::RowVector3d& v3,
    Eigen::Matrix<double, 12, 1>& g)
{
    g_PT(v0[0], v0[1], v0[2],
        v1[0], v1[1], v1[2],
        v2[0], v2[1], v2[2],
        v3[0], v3[1], v3[2],
        g.data());
}

inline void H_PT(double v01, double v02, double v03, double v11, double v12, double v13,
    double v21, double v22, double v23, double v31, double v32, double v33,
    double H[144])
{
    double t11;
    double t12;
    double t13;
    double t18;
    double t20;
    double t22;
    double t23;
    double t24;
    double t25;
    double t26;
    double t27;
    double t28;
    double t65;
    double t66;
    double t67;
    double t68;
    double t69;
    double t70;
    double t71;
    double t77;
    double t85;
    double t86;
    double t90;
    double t94;
    double t95;
    double t99;
    double t103;
    double t38;
    double t39;
    double t40;
    double t41;
    double t42;
    double t43;
    double t44;
    double t45;
    double t46;
    double t72;
    double t73;
    double t75;
    double t78;
    double t80;
    double t82;
    double t125;
    double t126;
    double t127;
    double t128;
    double t129;
    double t130;
    double t131;
    double t133;
    double t135;
    double t149;
    double t150;
    double t151;
    double t189;
    double t190;
    double t191;
    double t192;
    double t193;
    double t194;
    double t195;
    double t196;
    double t197;
    double t198;
    double t199;
    double t200;
    double t202;
    double t205;
    double t203;
    double t204;
    double t206;
    double t241;
    double t309;
    double t310;
    double t312;
    double t313;
    double t314;
    double t315;
    double t316;
    double t317;
    double t318;
    double t319;
    double t321;
    double t322;
    double t323;
    double t324;
    double t325;
    double t261;
    double t262;
    double t599;
    double t600;
    double t602;
    double t605;
    double t609;
    double t610;
    double t611;
    double t613;
    double t615;
    double t616;
    double t621;
    double t622;
    double t623;
    double t625;
    double t645;
    double t646_tmp;
    double t646;
    double t601;
    double t603;
    double t604;
    double t606;
    double t607;
    double t608;
    double t612;
    double t614;
    double t617;
    double t618;
    double t619;
    double t620;
    double t624;
    double t626;
    double t627;
    double t628;
    double t629;
    double t630;
    double t631;
    double t632;
    double t633;
    double t634;
    double t635;
    double t636;
    double t637;
    double t638;

    /* H_PT */
    /*     H = H_PT(V01,V02,V03,V11,V12,V13,V21,V22,V23,V31,V32,V33) */
    /*     This function was generated by the Symbolic Math Toolbox version 8.3. */
    /*     10-Jun-2019 17:42:25 */
    t11 = -v11 + v01;
    t12 = -v12 + v02;
    t13 = -v13 + v03;
    t18 = -v21 + v11;
    t20 = -v22 + v12;
    t22 = -v23 + v13;
    t23 = -v31 + v11;
    t24 = -v32 + v12;
    t25 = -v33 + v13;
    t26 = -v31 + v21;
    t27 = -v32 + v22;
    t28 = -v33 + v23;
    t65 = t18 * t24;
    t66 = t20 * t23;
    t67 = t18 * t25;
    t68 = t22 * t23;
    t69 = t20 * t25;
    t70 = t22 * t24;
    t71 = t18 * t23 * 2.0;
    t77 = t20 * t24 * 2.0;
    t85 = t22 * t25 * 2.0;
    t86 = t18 * t26 * 2.0;
    t90 = t20 * t27 * 2.0;
    t94 = t22 * t28 * 2.0;
    t95 = t23 * t26 * 2.0;
    t99 = t24 * t27 * 2.0;
    t103 = t25 * t28 * 2.0;
    t38 = t18 * t18 * 2.0;
    t39 = t20 * t20 * 2.0;
    t40 = t22 * t22 * 2.0;
    t41 = t23 * t23 * 2.0;
    t42 = t24 * t24 * 2.0;
    t43 = t25 * t25 * 2.0;
    t44 = t26 * t26 * 2.0;
    t45 = t27 * t27 * 2.0;
    t46 = t28 * t28 * 2.0;
    t72 = t65 * 2.0;
    t73 = t66 * 2.0;
    t75 = t67 * 2.0;
    t78 = t68 * 2.0;
    t80 = t69 * 2.0;
    t82 = t70 * 2.0;
    t125 = t11 * t20 + -(t12 * t18);
    t126 = t11 * t22 + -(t13 * t18);
    t127 = t12 * t22 + -(t13 * t20);
    t128 = t11 * t24 + -(t12 * t23);
    t129 = t11 * t25 + -(t13 * t23);
    t130 = t12 * t25 + -(t13 * t24);
    t131 = t65 + -t66;
    t133 = t67 + -t68;
    t135 = t69 + -t70;
    t149 = t131 * t131;
    t150 = t133 * t133;
    t151 = t135 * t135;
    t189 = (t11 * t27 + -(t12 * t26)) + t131;
    t190 = (t11 * t28 + -(t13 * t26)) + t133;
    t191 = (t12 * t28 + -(t13 * t27)) + t135;
    t192 = t20 * t131 * 2.0 + t22 * t133 * 2.0;
    t193 = t18 * t133 * 2.0 + t20 * t135 * 2.0;
    t194 = t24 * t131 * 2.0 + t25 * t133 * 2.0;
    t195 = t23 * t133 * 2.0 + t24 * t135 * 2.0;
    t196 = t27 * t131 * 2.0 + t28 * t133 * 2.0;
    t197 = t26 * t133 * 2.0 + t27 * t135 * 2.0;
    t198 = t18 * t131 * 2.0 + -(t22 * t135 * 2.0);
    t199 = t23 * t131 * 2.0 + -(t25 * t135 * 2.0);
    t200 = t26 * t131 * 2.0 + -(t28 * t135 * 2.0);
    t202 = 1.0 / ((t149 + t150) + t151);
    t205 = (t13 * t131 + t11 * t135) + -(t12 * t133);
    t203 = t202 * t202;
    t204 = pow(t202, 3.0);
    t206 = t205 * t205;
    t241 = t131 * t135 * t202 * 2.0;
    t309 = t11 * t202 * t205 * 2.0;
    t310 = t12 * t202 * t205 * 2.0;
    t13 = t13 * t202 * t205 * 2.0;
    t312 = (-v21 + v01) * t202 * t205 * 2.0;
    t313 = (-v22 + v02) * t202 * t205 * 2.0;
    t314 = (-v23 + v03) * t202 * t205 * 2.0;
    t315 = (-v31 + v01) * t202 * t205 * 2.0;
    t316 = t18 * t202 * t205 * 2.0;
    t317 = (-v32 + v02) * t202 * t205 * 2.0;
    t318 = t20 * t202 * t205 * 2.0;
    t319 = (-v33 + v03) * t202 * t205 * 2.0;
    t11 = t22 * t202 * t205 * 2.0;
    t321 = t23 * t202 * t205 * 2.0;
    t322 = t24 * t202 * t205 * 2.0;
    t323 = t25 * t202 * t205 * 2.0;
    t324 = t26 * t202 * t205 * 2.0;
    t325 = t27 * t202 * t205 * 2.0;
    t12 = t28 * t202 * t205 * 2.0;
    t261 = -(t131 * t133 * t202 * 2.0);
    t262 = -(t133 * t135 * t202 * 2.0);
    t599 = t130 * t135 * t202 * 2.0 + t135 * t194 * t203 * t205 * 2.0;
    t600 = -(t125 * t131 * t202 * 2.0) + t131 * t193 * t203 * t205 * 2.0;
    t602 = t129 * t133 * t202 * 2.0 + t133 * t199 * t203 * t205 * 2.0;
    t605 = -(t131 * t189 * t202 * 2.0) + t131 * t197 * t203 * t205 * 2.0;
    t609 = (t127 * t133 * t202 * 2.0 + -t11) + t133 * t192 * t203 * t205 * 2.0;
    t610 = (t126 * t135 * t202 * 2.0 + t11) + t135 * t198 * t203 * t205 * 2.0;
    t611 = (t130 * t131 * t202 * 2.0 + -t322) + t131 * t194 * t203 * t205 * 2.0;
    t613 = (t126 * t131 * t202 * 2.0 + -t316) + t131 * t198 * t203 * t205 * 2.0;
    t615 = (-(t125 * t135 * t202 * 2.0) + -t318) + t135 * t193 * t203 * t205 * 2.0;
    t616 = (-(t128 * t133 * t202 * 2.0) + -t321) + t133 * t195 * t203 * t205 * 2.0;
    t621 = (t133 * t191 * t202 * 2.0 + -t12) + t133 * t196 * t203 * t205 * 2.0;
    t622 = (t135 * t190 * t202 * 2.0 + t12) + t135 * t200 * t203 * t205 * 2.0;
    t623 = (t131 * t190 * t202 * 2.0 + -t324) + t131 * t200 * t203 * t205 * 2.0;
    t625 = (-(t135 * t189 * t202 * 2.0) + -t325) + t135 * t197 * t203 * t205 * 2.0;
    t645 = ((((t127 * t129 * t202 * 2.0 + -t13) + (t72 + -(t66 * 4.0)) * t203 * t206) + t129 * t192 * t203 * t205 * 2.0) + t127 * t199 * t203 * t205 * 2.0) + t192 * t199 * t204 * t206 * 2.0;
    t646_tmp = t203 * t206;
    t646 = ((((t126 * t130 * t202 * 2.0 + t13) + t646_tmp * (t73 - t65 * 4.0)) + t126 * t194 * t203 * t205 * 2.0) + t130 * t198 * t203 * t205 * 2.0) + t194 * t198 * t204 * t206 * 2.0;
    t601 = t128 * t131 * t202 * 2.0 + -(t131 * t195 * t203 * t205 * 2.0);
    t603 = -(t127 * t135 * t202 * 2.0) + -(t135 * t192 * t203 * t205 * 2.0);
    t604 = -(t126 * t133 * t202 * 2.0) + -(t133 * t198 * t203 * t205 * 2.0);
    t606 = -(t135 * t191 * t202 * 2.0) + -(t135 * t196 * t203 * t205 * 2.0);
    t607 = -(t133 * t190 * t202 * 2.0) + -(t133 * t200 * t203 * t205 * 2.0);
    t608 = (t125 * t133 * t202 * 2.0 + t316) + -(t133 * t193 * t203 * t205 * 2.0);
    t612 = (t128 * t135 * t202 * 2.0 + t322) + -(t135 * t195 * t203 * t205 * 2.0);
    t614 = (-(t127 * t131 * t202 * 2.0) + t318) + -(t131 * t192 * t203 * t205 * 2.0);
    t617 = (-(t130 * t133 * t202 * 2.0) + t323) + -(t133 * t194 * t203 * t205 * 2.0);
    t618 = (-(t129 * t131 * t202 * 2.0) + t321) + -(t131 * t199 * t203 * t205 * 2.0);
    t619 = (-(t129 * t135 * t202 * 2.0) + -t323) + -(t135 * t199 * t203 * t205 * 2.0);
    t620 = (t133 * t189 * t202 * 2.0 + t324) + -(t133 * t197 * t203 * t205 * 2.0);
    t624 = (-(t131 * t191 * t202 * 2.0) + t325) + -(t131 * t196 * t203 * t205 * 2.0);
    t626 = (((t125 * t127 * t202 * 2.0 + t18 * t22 * t203 * t206 * 2.0) + t125 * t192 * t203 * t205 * 2.0) + -(t127 * t193 * t203 * t205 * 2.0)) + -(t192 * t193 * t204 * t206 * 2.0);
    t627 = (((t128 * t130 * t202 * 2.0 + t23 * t25 * t203 * t206 * 2.0) + t128 * t194 * t203 * t205 * 2.0) + -(t130 * t195 * t203 * t205 * 2.0)) + -(t194 * t195 * t204 * t206 * 2.0);
    t628 = (((-(t125 * t126 * t202 * 2.0) + t20 * t22 * t203 * t206 * 2.0) + t126 * t193 * t203 * t205 * 2.0) + -(t125 * t198 * t203 * t205 * 2.0)) + t193 * t198 * t204 * t206 * 2.0;
    t629 = (((-(t128 * t129 * t202 * 2.0) + t24 * t25 * t203 * t206 * 2.0) + t129 * t195 * t203 * t205 * 2.0) + -(t128 * t199 * t203 * t205 * 2.0)) + t195 * t199 * t204 * t206 * 2.0;
    t630 = (((-(t126 * t127 * t202 * 2.0) + t18 * t20 * t203 * t206 * 2.0) + -(t126 * t192 * t203 * t205 * 2.0)) + -(t127 * t198 * t203 * t205 * 2.0)) + -(t192 * t198 * t204 * t206 * 2.0);
    t631 = (((-(t129 * t130 * t202 * 2.0) + t23 * t24 * t203 * t206 * 2.0) + -(t129 * t194 * t203 * t205 * 2.0)) + -(t130 * t199 * t203 * t205 * 2.0)) + -(t194 * t199 * t204 * t206 * 2.0);
    t632 = (((-(t125 * t128 * t202 * 2.0) + (t71 + t77) * t203 * t206) + t128 * t193 * t203 * t205 * 2.0) + t125 * t195 * t203 * t205 * 2.0) + -(t193 * t195 * t204 * t206 * 2.0);
    t633 = (((-(t127 * t130 * t202 * 2.0) + (t77 + t85) * t203 * t206) + -(t130 * t192 * t203 * t205 * 2.0)) + -(t127 * t194 * t203 * t205 * 2.0)) + -(t192 * t194 * t204 * t206 * 2.0);
    t634 = (((-(t126 * t129 * t202 * 2.0) + (t71 + t85) * t203 * t206) + -(t129 * t198 * t203 * t205 * 2.0)) + -(t126 * t199 * t203 * t205 * 2.0)) + -(t198 * t199 * t204 * t206 * 2.0);
    t635 = (((t127 * t191 * t202 * 2.0 + -((t90 + t94) * t203 * t206)) + t127 * t196 * t203 * t205 * 2.0) + t191 * t192 * t203 * t205 * 2.0) + t192 * t196 * t204 * t206 * 2.0;
    t636 = (((-(t128 * t189 * t202 * 2.0) + (t95 + t99) * t203 * t206) + t128 * t197 * t203 * t205 * 2.0) + t189 * t195 * t203 * t205 * 2.0) + -(t195 * t197 * t204 * t206 * 2.0);
    t637 = (((t125 * t189 * t202 * 2.0 + -((t86 + t90) * t203 * t206)) + -(t125 * t197 * t203 * t205 * 2.0)) + -(t189 * t193 * t203 * t205 * 2.0)) + t193 * t197 * t204 * t206 * 2.0;
    t638 = (((-(t130 * t191 * t202 * 2.0) + (t99 + t103) * t203 * t206) + -(t130 * t196 * t203 * t205 * 2.0)) + -(t191 * t194 * t203 * t205 * 2.0)) + -(t194 * t196 * t204 * t206 * 2.0);
    t86 = (((t126 * t190 * t202 * 2.0 + -((t86 + t94) * t203 * t206)) + t126 * t200 * t203 * t205 * 2.0) + t190 * t198 * t203 * t205 * 2.0) + t198 * t200 * t204 * t206 * 2.0;
    t71 = (((-(t129 * t190 * t202 * 2.0) + (t95 + t103) * t203 * t206) + -(t129 * t200 * t203 * t205 * 2.0)) + -(t190 * t199 * t203 * t205 * 2.0)) + -(t199 * t200 * t204 * t206 * 2.0);
    t85 = (((t189 * t191 * t202 * 2.0 + t26 * t28 * t203 * t206 * 2.0) + t189 * t196 * t203 * t205 * 2.0) + -(t191 * t197 * t203 * t205 * 2.0)) + -(t196 * t197 * t204 * t206 * 2.0);
    t90 = (((-(t189 * t190 * t202 * 2.0) + t27 * t28 * t203 * t206 * 2.0) + t190 * t197 * t203 * t205 * 2.0) + -(t189 * t200 * t203 * t205 * 2.0)) + t197 * t200 * t204 * t206 * 2.0;
    t99 = (((-(t190 * t191 * t202 * 2.0) + t26 * t27 * t203 * t206 * 2.0) + -(t190 * t196 * t203 * t205 * 2.0)) + -(t191 * t200 * t203 * t205 * 2.0)) + -(t196 * t200 * t204 * t206 * 2.0);
    t77 = ((((-(t127 * t128 * t202 * 2.0) + t310) + (t75 + -(t68 * 4.0)) * t203 * t206) + t127 * t195 * t203 * t205 * 2.0) + -(t128 * t192 * t203 * t205 * 2.0)) + t192 * t195 * t204 * t206 * 2.0;
    t131 = ((((t126 * t128 * t202 * 2.0 + -t309) + (t80 + -(t70 * 4.0)) * t203 * t206) + t128 * t198 * t203 * t205 * 2.0) + -(t126 * t195 * t203 * t205 * 2.0)) + -(t195 * t198 * t204 * t206 * 2.0);
    t133 = ((((-(t125 * t130 * t202 * 2.0) + -t310) + t646_tmp * (t78 - t67 * 4.0))
                + t130 * t193 * t203 * t205 * 2.0)
               + -(t125 * t194 * t203 * t205 * 2.0))
        + t193 * t194 * t204 * t206 * 2.0;
    t325 = ((((t125 * t129 * t202 * 2.0 + t309) + t646_tmp * (t82 - t69 * 4.0)) + t125 * t199 * t203 * t205 * 2.0) + -(t129 * t193 * t203 * t205 * 2.0))
        + -(t193 * t199 * t204 * t206 * 2.0);
    t135 = ((((t125 * t191 * t202 * 2.0 + t313) + ((t75 + t18 * t28 * 2.0) + -t78) * t203 * t206) + t125 * t196 * t203 * t205 * 2.0) + -(t191 * t193 * t203 * t205 * 2.0)) + -(t193 * t196 * t204 * t206 * 2.0);
    t324 = ((((t127 * t189 * t202 * 2.0 + -t313) + ((t78 + t22 * t26 * 2.0) + -t75) * t203 * t206) + -(t127 * t197 * t203 * t205 * 2.0)) + t189 * t192 * t203 * t205 * 2.0) + -(t192 * t197 * t204 * t206 * 2.0);
    t318 = ((((-(t126 * t189 * t202 * 2.0) + t312) + ((t82 + t22 * t27 * 2.0) + -t80) * t203 * t206) + t126 * t197 * t203 * t205 * 2.0) + -(t189 * t198 * t203 * t205 * 2.0)) + t197 * t198 * t204 * t206 * 2.0;
    t321 = ((((-(t130 * t189 * t202 * 2.0) + t317) + -(((t78 + t25 * t26 * 2.0) + -t75) * t203 * t206)) + t130 * t197 * t203 * t205 * 2.0) + -(t189 * t194 * t203 * t205 * 2.0)) + t194 * t197 * t204 * t206 * 2.0;
    t323 = ((((t129 * t191 * t202 * 2.0 + t319) + -(((t72 + t23 * t27 * 2.0) + -t73) * t203 * t206)) + t129 * t196 * t203 * t205 * 2.0) + t191 * t199 * t203 * t205 * 2.0) + t196 * t199 * t204 * t206 * 2.0;
    t322 = ((((-(t125 * t190 * t202 * 2.0) + -t312) + ((t80 + t20 * t28 * 2.0) + -t82) * t203 * t206) + -(t125 * t200 * t203 * t205 * 2.0)) + t190 * t193 * t203 * t205 * 2.0) + t193 * t200 * t204 * t206 * 2.0;
    t316 = ((((t130 * t190 * t202 * 2.0 + -t319) + -(((t73 + t24 * t26 * 2.0) + -t72) * t203 * t206)) + t130 * t200 * t203 * t205 * 2.0) + t190 * t194 * t203 * t205 * 2.0) + t194 * t200 * t204 * t206 * 2.0;
    t65 = ((((-(t128 * t191 * t202 * 2.0) + -t317) + -(((t75 + t23 * t28 * 2.0) + -t78) * t203 * t206)) + -(t128 * t196 * t203 * t205 * 2.0)) + t191 * t195 * t203 * t205 * 2.0) + t195 * t196 * t204 * t206 * 2.0;
    t66 = ((((-(t127 * t190 * t202 * 2.0) + t314) + ((t73 + t20 * t26 * 2.0) + -t72) * t203 * t206) + -(t127 * t200 * t203 * t205 * 2.0)) + -(t190 * t192 * t203 * t205 * 2.0)) + -(t192 * t200 * t204 * t206 * 2.0);
    t13 = ((((t128 * t190 * t202 * 2.0 + t315) + -(((t80 + t24 * t28 * 2.0) + -t82) * t203 * t206)) + t128 * t200 * t203 * t205 * 2.0) + -(t190 * t195 * t203 * t205 * 2.0)) + -(t195 * t200 * t204 * t206 * 2.0);
    t12 = ((((-(t126 * t191 * t202 * 2.0) + -t314) + ((t72 + t18 * t27 * 2.0) + -t73) * t203 * t206) + -(t126 * t196 * t203 * t205 * 2.0)) + -(t191 * t198 * t203 * t205 * 2.0)) + -(t196 * t198 * t204 * t206 * 2.0);
    t11 = ((((t129 * t189 * t202 * 2.0 + -t315) + -(((t82 + t25 * t27 * 2.0) + -t80) * t203 * t206)) + -(t129 * t197 * t203 * t205 * 2.0)) + t189 * t199 * t203 * t205 * 2.0) + -(t197 * t199 * t204 * t206 * 2.0);
    H[0] = t151 * t202 * 2.0;
    H[1] = t262;
    H[2] = t241;
    H[3] = t606;
    H[4] = t622;
    H[5] = t625;
    H[6] = t599;
    H[7] = t619;
    H[8] = t612;
    H[9] = t603;
    H[10] = t610;
    H[11] = t615;
    H[12] = t262;
    H[13] = t150 * t202 * 2.0;
    H[14] = t261;
    H[15] = t621;
    H[16] = t607;
    H[17] = t620;
    H[18] = t617;
    H[19] = t602;
    H[20] = t616;
    H[21] = t609;
    H[22] = t604;
    H[23] = t608;
    H[24] = t241;
    H[25] = t261;
    H[26] = t149 * t202 * 2.0;
    H[27] = t624;
    H[28] = t623;
    H[29] = t605;
    H[30] = t611;
    H[31] = t618;
    H[32] = t601;
    H[33] = t614;
    H[34] = t613;
    H[35] = t600;
    H[36] = t606;
    H[37] = t621;
    H[38] = t624;
    H[39] = ((t191 * t191 * t202 * 2.0 + t196 * t196 * t204 * t206 * 2.0) - t646_tmp * (t45 + t46)) + t191 * t196 * t203 * t205 * 4.0;
    H[40] = t99;
    H[41] = t85;
    H[42] = t638;
    H[43] = t323;
    H[44] = t65;
    H[45] = t635;
    H[46] = t12;
    H[47] = t135;
    H[48] = t622;
    H[49] = t607;
    H[50] = t623;
    H[51] = t99;
    H[52] = ((t190 * t190 * t202 * 2.0 + t200 * t200 * t204 * t206 * 2.0) - t646_tmp * (t44 + t46)) + t190 * t200 * t203 * t205 * 4.0;
    H[53] = t90;
    H[54] = t316;
    H[55] = t71;
    H[56] = t13;
    H[57] = t66;
    H[58] = t86;
    H[59] = t322;
    H[60] = t625;
    H[61] = t620;
    H[62] = t605;
    H[63] = t85;
    H[64] = t90;
    H[65] = ((t189 * t189 * t202 * 2.0 + t197 * t197 * t204 * t206 * 2.0) - t646_tmp * (t44 + t45)) - t189 * t197 * t203 * t205 * 4.0;
    H[66] = t321;
    H[67] = t11;
    H[68] = t636;
    H[69] = t324;
    H[70] = t318;
    H[71] = t637;
    H[72] = t599;
    H[73] = t617;
    H[74] = t611;
    H[75] = t638;
    H[76] = t316;
    H[77] = t321;
    H[78] = ((t130 * t130 * t202 * 2.0 + t194 * t194 * t204 * t206 * 2.0) - t646_tmp * (t42 + t43)) + t130 * t194 * t203 * t205 * 4.0;
    H[79] = t631;
    H[80] = t627;
    H[81] = t633;
    H[82] = t646;
    H[83] = t133;
    H[84] = t619;
    H[85] = t602;
    H[86] = t618;
    H[87] = t323;
    H[88] = t71;
    H[89] = t11;
    H[90] = t631;
    H[91] = ((t129 * t129 * t202 * 2.0 + t199 * t199 * t204 * t206 * 2.0) - t646_tmp * (t41 + t43)) + t129 * t199 * t203 * t205 * 4.0;
    H[92] = t629;
    H[93] = t645;
    H[94] = t634;
    H[95] = t325;
    H[96] = t612;
    H[97] = t616;
    H[98] = t601;
    H[99] = t65;
    H[100] = t13;
    H[101] = t636;
    H[102] = t627;
    H[103] = t629;
    H[104] = ((t128 * t128 * t202 * 2.0 + t195 * t195 * t204 * t206 * 2.0) - t646_tmp * (t41 + t42)) - t128 * t195 * t203 * t205 * 4.0;
    H[105] = t77;
    H[106] = t131;
    H[107] = t632;
    H[108] = t603;
    H[109] = t609;
    H[110] = t614;
    H[111] = t635;
    H[112] = t66;
    H[113] = t324;
    H[114] = t633;
    H[115] = t645;
    H[116] = t77;
    H[117] = ((t127 * t127 * t202 * 2.0 + t192 * t192 * t204 * t206 * 2.0) - t646_tmp * (t39 + t40)) + t127 * t192 * t203 * t205 * 4.0;
    H[118] = t630;
    H[119] = t626;
    H[120] = t610;
    H[121] = t604;
    H[122] = t613;
    H[123] = t12;
    H[124] = t86;
    H[125] = t318;
    H[126] = t646;
    H[127] = t634;
    H[128] = t131;
    H[129] = t630;
    H[130] = ((t126 * t126 * t202 * 2.0 + t198 * t198 * t204 * t206 * 2.0) - t646_tmp * (t38 + t40)) + t126 * t198 * t203 * t205 * 4.0;
    H[131] = t628;
    H[132] = t615;
    H[133] = t608;
    H[134] = t600;
    H[135] = t135;
    H[136] = t322;
    H[137] = t637;
    H[138] = t133;
    H[139] = t325;
    H[140] = t632;
    H[141] = t626;
    H[142] = t628;
    H[143] = ((t125 * t125 * t202 * 2.0 + t193 * t193 * t204 * t206 * 2.0) - t646_tmp * (t38 + t39)) - t125 * t193 * t203 * t205 * 4.0;
}

inline void H_PT(const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    const Eigen::RowVector3d& v3,
    Eigen::Matrix<double, 12, 12>& H)
{
    H_PT(v0[0], v0[1], v0[2],
        v1[0], v1[1], v1[2],
        v2[0], v2[1], v2[2],
        v3[0], v3[1], v3[2],
        H.data());
}

inline void derivTest_PT(void)
{
    Eigen::RowVector3d v[4];
    v[0] << 0.0, 0.0, 0.0;
    v[1] << 1.0, 0.1, 0.0;
    v[2] << 0.0, 1.1, -0.1;
    v[3] << 0.0, 0.1, -1.1;

    double d0;
    d_PT(v[0], v[1], v[2], v[3], d0);

    double eps = 1.0e-6;
    Eigen::Matrix<double, 12, 1> fd, s;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            Eigen::RowVector3d vv[4];
            vv[0] = v[0];
            vv[1] = v[1];
            vv[2] = v[2];
            vv[3] = v[3];
            vv[i][j] += eps;
            double d;
            d_PT(vv[0], vv[1], vv[2], vv[3], d);
            fd[i * 3 + j] = (d - d0) / eps;
        }
    }

    g_PT(v[0], v[1], v[2], v[3], s);

    std::cout << s.transpose() << std::endl;
    std::cout << fd.transpose() << std::endl;
    std::cout << "relError = " << (s - fd).norm() / fd.norm() << std::endl;

    Eigen::Matrix<double, 12, 12> Hfd, Hs;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            Eigen::RowVector3d vv[4];
            vv[0] = v[0];
            vv[1] = v[1];
            vv[2] = v[2];
            vv[3] = v[3];
            vv[i][j] += eps;
            Eigen::Matrix<double, 12, 1> g;
            g_PT(vv[0], vv[1], vv[2], vv[3], g);
            Hfd.col(i * 3 + j) = (g - s) / eps;
        }
    }

    H_PT(v[0], v[1], v[2], v[3], Hs);

    std::cout << Hs.transpose() << std::endl;
    std::cout << Hfd.transpose() << std::endl;
    std::cout << "relError = " << (Hs - Hfd).norm() / Hfd.norm() << std::endl;
}

inline void d_EE(const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    const Eigen::RowVector3d& v3,
    double& d)
{
    Eigen::RowVector3d b = (v1 - v0).cross(v3 - v2);
    double aTb = (v2 - v0).dot(b);
    d = aTb * aTb / b.squaredNorm();
}

inline void g_EE(double v01, double v02, double v03, double v11, double v12, double v13,
    double v21, double v22, double v23, double v31, double v32, double v33,
    double g[12])
{
    double t11;
    double t12;
    double t13;
    double t14;
    double t15;
    double t16;
    double t17;
    double t18;
    double t19;
    double t32;
    double t33;
    double t34;
    double t35;
    double t36;
    double t37;
    double t44;
    double t45;
    double t46;
    double t75;
    double t77;
    double t76;
    double t78;
    double t79;
    double t80;
    double t81;
    double t83;

    /* G_EE */
    /*     G = G_EE(V01,V02,V03,V11,V12,V13,V21,V22,V23,V31,V32,V33) */
    /*     This function was generated by the Symbolic Math Toolbox version 8.3. */
    /*     14-Jun-2019 13:58:25 */
    t11 = -v11 + v01;
    t12 = -v12 + v02;
    t13 = -v13 + v03;
    t14 = -v21 + v01;
    t15 = -v22 + v02;
    t16 = -v23 + v03;
    t17 = -v31 + v21;
    t18 = -v32 + v22;
    t19 = -v33 + v23;
    t32 = t14 * t18;
    t33 = t15 * t17;
    t34 = t14 * t19;
    t35 = t16 * t17;
    t36 = t15 * t19;
    t37 = t16 * t18;
    t44 = t11 * t18 + -(t12 * t17);
    t45 = t11 * t19 + -(t13 * t17);
    t46 = t12 * t19 + -(t13 * t18);
    t75 = 1.0 / ((t44 * t44 + t45 * t45) + t46 * t46);
    t77 = (t16 * t44 + t14 * t46) + -(t15 * t45);
    t76 = t75 * t75;
    t78 = t77 * t77;
    t79 = (t12 * t44 * 2.0 + t13 * t45 * 2.0) * t76 * t78;
    t80 = (t11 * t45 * 2.0 + t12 * t46 * 2.0) * t76 * t78;
    t81 = (t18 * t44 * 2.0 + t19 * t45 * 2.0) * t76 * t78;
    t18 = (t17 * t45 * 2.0 + t18 * t46 * 2.0) * t76 * t78;
    t83 = (t11 * t44 * 2.0 + -(t13 * t46 * 2.0)) * t76 * t78;
    t19 = (t17 * t44 * 2.0 + -(t19 * t46 * 2.0)) * t76 * t78;
    t76 = t75 * t77;
    g[0] = -t81 + t76 * ((-t36 + t37) + t46) * 2.0;
    g[1] = t19 - t76 * ((-t34 + t35) + t45) * 2.0;
    g[2] = t18 + t76 * ((-t32 + t33) + t44) * 2.0;
    g[3] = t81 + t76 * (t36 - t37) * 2.0;
    g[4] = -t19 - t76 * (t34 - t35) * 2.0;
    g[5] = -t18 + t76 * (t32 - t33) * 2.0;
    t17 = t12 * t16 + -(t13 * t15);
    g[6] = t79 - t76 * (t17 + t46) * 2.0;
    t18 = t11 * t16 + -(t13 * t14);
    g[7] = -t83 + t76 * (t18 + t45) * 2.0;
    t19 = t11 * t15 + -(t12 * t14);
    g[8] = -t80 - t76 * (t19 + t44) * 2.0;
    g[9] = -t79 + t76 * t17 * 2.0;
    g[10] = t83 - t76 * t18 * 2.0;
    g[11] = t80 + t76 * t19 * 2.0;
}

inline void g_EE(const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    const Eigen::RowVector3d& v3,
    Eigen::Matrix<double, 12, 1>& g)
{
    g_EE(v0[0], v0[1], v0[2],
        v1[0], v1[1], v1[2],
        v2[0], v2[1], v2[2],
        v3[0], v3[1], v3[2],
        g.data());
}

inline void H_EE(double v01, double v02, double v03, double v11, double v12, double v13,
    double v21, double v22, double v23, double v31, double v32, double v33,
    double H[144])
{
    double t11;
    double t12;
    double t13;
    double t14;
    double t15;
    double t16;
    double t26;
    double t27;
    double t28;
    double t47;
    double t48;
    double t49;
    double t50;
    double t51;
    double t52;
    double t53;
    double t54;
    double t55;
    double t56;
    double t57;
    double t58;
    double t59;
    double t65;
    double t73;
    double t35;
    double t36;
    double t37;
    double t38;
    double t39;
    double t40;
    double t98;
    double t99;
    double t100;
    double t101;
    double t103;
    double t105;
    double t107;
    double t108;
    double t109;
    double t137;
    double t138;
    double t139;
    double t140;
    double t141;
    double t142;
    double t143;
    double t144;
    double t145;
    double t146;
    double t147;
    double t148;
    double t156;
    double t159;
    double t157;
    double t262;
    double t263;
    double t264;
    double t265;
    double t266;
    double t267;
    double t268;
    double t269;
    double t270;
    double t271;
    double t272;
    double t273;
    double t274;
    double t275;
    double t276;
    double t277;
    double t278;
    double t279;
    double t298;
    double t299;
    double t300;
    double t301;
    double t302;
    double t303;
    double t310;
    double t311;
    double t312;
    double t313;
    double t314;
    double t315;
    double t322;
    double t323;
    double t325;
    double t326;
    double t327;
    double t328;
    double t329;
    double t330;
    double t335;
    double t337;
    double t339;
    double t340;
    double t341;
    double t342;
    double t343;
    double t345;
    double t348;
    double t353;
    double t356;
    double t358;
    double t359;
    double t360;
    double t362;
    double t367;
    double t368;
    double t369;
    double t371;
    double t374;
    double t377;
    double t382;
    double t386;
    double t387;
    double t398;
    double t399;
    double t403;
    double t408;
    double t423;
    double t424;
    double t427;
    double t428;
    double t431;
    double t432;
    double t433;
    double t434;
    double t437;
    double t438;
    double t441;
    double t442;
    double t446;
    double t451;
    double t455;
    double t456;
    double t467;
    double t468;
    double t472;
    double t477;
    double t491;
    double t492;
    double t495;
    double t497;
    double t499;
    double t500;
    double t503;
    double t504;
    double t506;
    double t508;
    double t550;
    double t568;
    double t519_tmp;
    double b_t519_tmp;
    double t519;
    double t520_tmp;
    double b_t520_tmp;
    double t520;
    double t521_tmp;
    double b_t521_tmp;
    double t521;
    double t522_tmp;
    double b_t522_tmp;
    double t522;
    double t523_tmp;
    double b_t523_tmp;
    double t523;
    double t524_tmp;
    double b_t524_tmp;
    double t524;
    double t525;
    double t526;
    double t527;
    double t528;
    double t529;
    double t530;
    double t531;
    double t532;
    double t533;
    double t534;
    double t535;
    double t536;
    double t537;
    double t538;
    double t539;
    double t540;
    double t542;
    double t543;
    double t544;

    /* H_EE */
    /*     H = H_EE(V01,V02,V03,V11,V12,V13,V21,V22,V23,V31,V32,V33) */
    /*     This function was generated by the Symbolic Math Toolbox version 8.3. */
    /*     14-Jun-2019 13:58:38 */
    t11 = -v11 + v01;
    t12 = -v12 + v02;
    t13 = -v13 + v03;
    t14 = -v21 + v01;
    t15 = -v22 + v02;
    t16 = -v23 + v03;
    t26 = -v31 + v21;
    t27 = -v32 + v22;
    t28 = -v33 + v23;
    t47 = t11 * t27;
    t48 = t12 * t26;
    t49 = t11 * t28;
    t50 = t13 * t26;
    t51 = t12 * t28;
    t52 = t13 * t27;
    t53 = t14 * t27;
    t54 = t15 * t26;
    t55 = t14 * t28;
    t56 = t16 * t26;
    t57 = t15 * t28;
    t58 = t16 * t27;
    t59 = t11 * t26 * 2.0;
    t65 = t12 * t27 * 2.0;
    t73 = t13 * t28 * 2.0;
    t35 = t11 * t11 * 2.0;
    t36 = t12 * t12 * 2.0;
    t37 = t13 * t13 * 2.0;
    t38 = t26 * t26 * 2.0;
    t39 = t27 * t27 * 2.0;
    t40 = t28 * t28 * 2.0;
    t98 = t11 * t15 + -(t12 * t14);
    t99 = t11 * t16 + -(t13 * t14);
    t100 = t12 * t16 + -(t13 * t15);
    t101 = t47 + -t48;
    t103 = t49 + -t50;
    t105 = t51 + -t52;
    t107 = t53 + -t54;
    t108 = t55 + -t56;
    t109 = t57 + -t58;
    t137 = t98 + t101;
    t138 = t99 + t103;
    t139 = t100 + t105;
    t140 = (t54 + -t53) + t101;
    t141 = (t56 + -t55) + t103;
    t142 = (t58 + -t57) + t105;
    t143 = t12 * t101 * 2.0 + t13 * t103 * 2.0;
    t144 = t11 * t103 * 2.0 + t12 * t105 * 2.0;
    t145 = t27 * t101 * 2.0 + t28 * t103 * 2.0;
    t146 = t26 * t103 * 2.0 + t27 * t105 * 2.0;
    t147 = t11 * t101 * 2.0 + -(t13 * t105 * 2.0);
    t148 = t26 * t101 * 2.0 + -(t28 * t105 * 2.0);
    t156 = 1.0 / ((t101 * t101 + t103 * t103) + t105 * t105);
    t159 = (t16 * t101 + t14 * t105) + -(t15 * t103);
    t157 = t156 * t156;
    t57 = pow(t156, 3.0);
    t58 = t159 * t159;
    t262 = t11 * t156 * t159 * 2.0;
    t263 = t12 * t156 * t159 * 2.0;
    t264 = t13 * t156 * t159 * 2.0;
    t265 = t14 * t156 * t159 * 2.0;
    t266 = t15 * t156 * t159 * 2.0;
    t267 = t16 * t156 * t159 * 2.0;
    t268 = (-v31 + v01) * t156 * t159 * 2.0;
    t269 = (-v21 + v11) * t156 * t159 * 2.0;
    t270 = (-v32 + v02) * t156 * t159 * 2.0;
    t271 = (-v22 + v12) * t156 * t159 * 2.0;
    t272 = (-v33 + v03) * t156 * t159 * 2.0;
    t273 = (-v23 + v13) * t156 * t159 * 2.0;
    t274 = (-v31 + v11) * t156 * t159 * 2.0;
    t275 = (-v32 + v12) * t156 * t159 * 2.0;
    t276 = (-v33 + v13) * t156 * t159 * 2.0;
    t277 = t26 * t156 * t159 * 2.0;
    t278 = t27 * t156 * t159 * 2.0;
    t279 = t28 * t156 * t159 * 2.0;
    t298 = t11 * t12 * t157 * t58 * 2.0;
    t299 = t11 * t13 * t157 * t58 * 2.0;
    t300 = t12 * t13 * t157 * t58 * 2.0;
    t301 = t26 * t27 * t157 * t58 * 2.0;
    t302 = t26 * t28 * t157 * t58 * 2.0;
    t303 = t27 * t28 * t157 * t58 * 2.0;
    t310 = (t35 + t36) * t157 * t58;
    t311 = (t35 + t37) * t157 * t58;
    t312 = (t36 + t37) * t157 * t58;
    t313 = (t38 + t39) * t157 * t58;
    t314 = (t38 + t40) * t157 * t58;
    t315 = (t39 + t40) * t157 * t58;
    t322 = (t59 + t65) * t157 * t58;
    t323 = (t59 + t73) * t157 * t58;
    t59 = (t65 + t73) * t157 * t58;
    t325 = (t47 * 2.0 + -(t48 * 4.0)) * t157 * t58;
    t53 = -t157 * t58;
    t56 = t48 * 2.0 - t47 * 4.0;
    t326 = t53 * t56;
    t327 = (t49 * 2.0 + -(t50 * 4.0)) * t157 * t58;
    t55 = t50 * 2.0 - t49 * 4.0;
    t328 = t53 * t55;
    t329 = (t51 * 2.0 + -(t52 * 4.0)) * t157 * t58;
    t54 = t52 * 2.0 - t51 * 4.0;
    t330 = t53 * t54;
    t53 = t157 * t58;
    t335 = t53 * t56;
    t337 = t53 * t55;
    t339 = t53 * t54;
    t340 = t143 * t143 * t57 * t58 * 2.0;
    t341 = t144 * t144 * t57 * t58 * 2.0;
    t342 = t145 * t145 * t57 * t58 * 2.0;
    t343 = t146 * t146 * t57 * t58 * 2.0;
    t345 = t147 * t147 * t57 * t58 * 2.0;
    t348 = t148 * t148 * t57 * t58 * 2.0;
    t36 = t98 * t143 * t157 * t159 * 2.0;
    t353 = t99 * t143 * t157 * t159 * 2.0;
    t356 = t99 * t144 * t157 * t159 * 2.0;
    t65 = t100 * t144 * t157 * t159 * 2.0;
    t358 = t107 * t143 * t157 * t159 * 2.0;
    t359 = t98 * t145 * t157 * t159 * 2.0;
    t360 = t108 * t143 * t157 * t159 * 2.0;
    t54 = t107 * t144 * t157 * t159 * 2.0;
    t362 = t99 * t145 * t157 * t159 * 2.0;
    t53 = t98 * t146 * t157 * t159 * 2.0;
    t56 = t109 * t143 * t157 * t159 * 2.0;
    t27 = t108 * t144 * t157 * t159 * 2.0;
    t55 = t100 * t145 * t157 * t159 * 2.0;
    t367 = t99 * t146 * t157 * t159 * 2.0;
    t368 = t109 * t144 * t157 * t159 * 2.0;
    t369 = t100 * t146 * t157 * t159 * 2.0;
    t38 = t107 * t145 * t157 * t159 * 2.0;
    t371 = t108 * t145 * t157 * t159 * 2.0;
    t374 = t108 * t146 * t157 * t159 * 2.0;
    t28 = t109 * t146 * t157 * t159 * 2.0;
    t377 = t98 * t147 * t157 * t159 * 2.0;
    t382 = t100 * t147 * t157 * t159 * 2.0;
    t386 = t107 * t147 * t157 * t159 * 2.0;
    t387 = t98 * t148 * t157 * t159 * 2.0;
    t103 = t108 * t147 * t157 * t159 * 2.0;
    t101 = t99 * t148 * t157 * t159 * 2.0;
    t398 = t109 * t147 * t157 * t159 * 2.0;
    t399 = t100 * t148 * t157 * t159 * 2.0;
    t403 = t107 * t148 * t157 * t159 * 2.0;
    t408 = t109 * t148 * t157 * t159 * 2.0;
    t73 = t137 * t143 * t157 * t159 * 2.0;
    t423 = t138 * t143 * t157 * t159 * 2.0;
    t424 = t138 * t144 * t157 * t159 * 2.0;
    t37 = t139 * t144 * t157 * t159 * 2.0;
    t427 = t140 * t143 * t157 * t159 * 2.0;
    t428 = t137 * t145 * t157 * t159 * 2.0;
    t16 = t140 * t144 * t157 * t159 * 2.0;
    t11 = t137 * t146 * t157 * t159 * 2.0;
    t431 = t141 * t143 * t157 * t159 * 2.0;
    t432 = t138 * t145 * t157 * t159 * 2.0;
    t433 = t141 * t144 * t157 * t159 * 2.0;
    t434 = t138 * t146 * t157 * t159 * 2.0;
    t105 = t142 * t143 * t157 * t159 * 2.0;
    t14 = t139 * t145 * t157 * t159 * 2.0;
    t437 = t142 * t144 * t157 * t159 * 2.0;
    t438 = t139 * t146 * t157 * t159 * 2.0;
    t35 = t140 * t145 * t157 * t159 * 2.0;
    t441 = t141 * t145 * t157 * t159 * 2.0;
    t442 = t141 * t146 * t157 * t159 * 2.0;
    t39 = t142 * t146 * t157 * t159 * 2.0;
    t446 = t137 * t147 * t157 * t159 * 2.0;
    t451 = t139 * t147 * t157 * t159 * 2.0;
    t455 = t140 * t147 * t157 * t159 * 2.0;
    t456 = t137 * t148 * t157 * t159 * 2.0;
    t13 = t141 * t147 * t157 * t159 * 2.0;
    t26 = t138 * t148 * t157 * t159 * 2.0;
    t467 = t142 * t147 * t157 * t159 * 2.0;
    t468 = t139 * t148 * t157 * t159 * 2.0;
    t472 = t140 * t148 * t157 * t159 * 2.0;
    t477 = t142 * t148 * t157 * t159 * 2.0;
    t47 = t143 * t144 * t57 * t58 * 2.0;
    t15 = t143 * t145 * t57 * t58 * 2.0;
    t491 = t143 * t146 * t57 * t58 * 2.0;
    t492 = t144 * t145 * t57 * t58 * 2.0;
    t12 = t144 * t146 * t57 * t58 * 2.0;
    t40 = t145 * t146 * t57 * t58 * 2.0;
    t495 = t143 * t147 * t57 * t58 * 2.0;
    t497 = t144 * t147 * t57 * t58 * 2.0;
    t499 = t143 * t148 * t57 * t58 * 2.0;
    t500 = t145 * t147 * t57 * t58 * 2.0;
    t503 = t146 * t147 * t57 * t58 * 2.0;
    t504 = t144 * t148 * t57 * t58 * 2.0;
    t506 = t145 * t148 * t57 * t58 * 2.0;
    t508 = t146 * t148 * t57 * t58 * 2.0;
    t57 = t147 * t148 * t57 * t58 * 2.0;
    t550 = ((((t98 * t109 * t156 * 2.0 + -t266) + t337) + t359) + t368) + t492;
    t568 = ((((t108 * t137 * t156 * 2.0 + -t268) + t330) + t27) + t456) + t504;
    t519_tmp = t139 * t143 * t157 * t159;
    b_t519_tmp = t100 * t143 * t157 * t159;
    t519 = (((-(t100 * t139 * t156 * 2.0) + t312) + -t340) + b_t519_tmp * 2.0) + t519_tmp * 2.0;
    t520_tmp = t140 * t146 * t157 * t159;
    b_t520_tmp = t107 * t146 * t157 * t159;
    t520 = (((t107 * t140 * t156 * 2.0 + t313) + -t343) + b_t520_tmp * 2.0) + -(t520_tmp * 2.0);
    t521_tmp = t142 * t145 * t157 * t159;
    b_t521_tmp = t109 * t145 * t157 * t159;
    t521 = (((t109 * t142 * t156 * 2.0 + t315) + -t342) + -(b_t521_tmp * 2.0)) + t521_tmp * 2.0;
    t522_tmp = t137 * t144 * t157 * t159;
    b_t522_tmp = t98 * t144 * t157 * t159;
    t522 = (((-(t98 * t137 * t156 * 2.0) + t310) + -t341) + -(b_t522_tmp * 2.0)) + -(t522_tmp * 2.0);
    t523_tmp = t138 * t147 * t157 * t159;
    b_t523_tmp = t99 * t147 * t157 * t159;
    t523 = (((-(t99 * t138 * t156 * 2.0) + t311) + -t345) + b_t523_tmp * 2.0) + t523_tmp * 2.0;
    t524_tmp = t141 * t148 * t157 * t159;
    b_t524_tmp = t108 * t148 * t157 * t159;
    t524 = (((t108 * t141 * t156 * 2.0 + t314) + -t348) + -(b_t524_tmp * 2.0)) + t524_tmp * 2.0;
    t525 = (((t98 * t100 * t156 * 2.0 + t299) + t65) + -t36) + -t47;
    t526 = (((t107 * t109 * t156 * 2.0 + t302) + t38) + -t28) + -t40;
    t527 = (((-(t98 * t99 * t156 * 2.0) + t300) + t377) + -t356) + t497;
    t528 = (((-(t99 * t100 * t156 * 2.0) + t298) + t353) + t382) + -t495;
    t529 = (((-(t107 * t108 * t156 * 2.0) + t303) + t374) + -t403) + t508;
    t530 = (((-(t108 * t109 * t156 * 2.0) + t301) + -t371) + -t408) + -t506;
    t531 = (((t98 * t107 * t156 * 2.0 + t322) + t54) + -t53) + -t12;
    t532 = (((t100 * t109 * t156 * 2.0 + t59) + t55) + -t56) + -t15;
    t533 = (((t99 * t108 * t156 * 2.0 + t323) + t101) + -t103) + -t57;
    t534 = (((t98 * t140 * t156 * 2.0 + -t322) + t53) + t16) + t12;
    t535 = (((-(t107 * t137 * t156 * 2.0) + -t322) + -t54) + t11) + t12;
    t536 = (((t100 * t142 * t156 * 2.0 + -t59) + -t55) + -t105) + t15;
    t537 = (((-(t109 * t139 * t156 * 2.0) + -t59) + t56) + -t14) + t15;
    t538 = (((t99 * t141 * t156 * 2.0 + -t323) + -t101) + -t13) + t57;
    t539 = (((-(t108 * t138 * t156 * 2.0) + -t323) + t103) + -t26) + t57;
    t540 = (((t137 * t139 * t156 * 2.0 + t299) + t37) + -t73) + -t47;
    t148 = (((t140 * t142 * t156 * 2.0 + t302) + t39) + -t35) + -t40;
    t542 = (((-(t137 * t138 * t156 * 2.0) + t300) + t446) + -t424) + t497;
    t543 = (((-(t138 * t139 * t156 * 2.0) + t298) + t423) + t451) + -t495;
    t544 = (((-(t140 * t141 * t156 * 2.0) + t303) + t472) + -t442) + t508;
    t53 = (((-(t141 * t142 * t156 * 2.0) + t301) + t441) + t477) + -t506;
    t157 = (((-(t139 * t142 * t156 * 2.0) + t59) + t105) + t14) + -t15;
    t159 = (((-(t137 * t140 * t156 * 2.0) + t322) + -t16) + -t11) + -t12;
    t147 = (((-(t138 * t141 * t156 * 2.0) + t323) + t13) + t26) + -t57;
    t146 = ((((t100 * t107 * t156 * 2.0 + t266) + t327) + -t358) + -t369) + t491;
    t145 = ((((-(t99 * t107 * t156 * 2.0) + -t265) + t329) + t367) + t386) + -t503;
    t144 = ((((-(t100 * t108 * t156 * 2.0) + -t267) + t325) + t360) + -t399) + t499;
    t143 = ((((-(t99 * t109 * t156 * 2.0) + t267) + t335) + -t362) + t398) + t500;
    t52 = ((((-(t98 * t108 * t156 * 2.0) + t265) + t339) + -t27) + -t387) + -t504;
    t51 = ((((t109 * t140 * t156 * 2.0 + -t278) + -t302) + t28) + t35) + t40;
    t50 = ((((-(t98 * t139 * t156 * 2.0) + t263) + -t299) + t36) + -t37) + t47;
    t49 = ((((t107 * t142 * t156 * 2.0 + t278) + -t302) + -t38) + -t39) + t40;
    t48 = ((((-(t100 * t137 * t156 * 2.0) + -t263) + -t299) + -t65) + t73) + t47;
    t47 = ((((t99 * t137 * t156 * 2.0 + t262) + -t300) + t356) + -t446) + -t497;
    t73 = ((((t100 * t138 * t156 * 2.0 + t264) + -t298) + -t382) + -t423) + t495;
    t65 = ((((-(t109 * t141 * t156 * 2.0) + t279) + -t301) + t408) + -t441) + t506;
    t59 = ((((t98 * t138 * t156 * 2.0 + -t262) + -t300) + -t377) + t424) + -t497;
    t40 = ((((t99 * t139 * t156 * 2.0 + -t264) + -t298) + -t353) + -t451) + t495;
    t39 = ((((-(t107 * t141 * t156 * 2.0) + -t277) + -t303) + t403) + t442) + -t508;
    t38 = ((((-(t108 * t142 * t156 * 2.0) + -t279) + -t301) + t371) + -t477) + t506;
    t37 = ((((-(t108 * t140 * t156 * 2.0) + t277) + -t303) + -t374) + -t472) + -t508;
    t36 = ((((t98 * t142 * t156 * 2.0 + t271) + t328) + -t359) + t437) + -t492;
    t35 = ((((-(t109 * t137 * t156 * 2.0) + t270) + t328) + -t368) + -t428) + -t492;
    t28 = ((((t100 * t140 * t156 * 2.0 + -t271) + -t327) + t369) + -t427) + -t491;
    t27 = ((((-(t98 * t141 * t156 * 2.0) + -t269) + t330) + t387) + -t433) + t504;
    t26 = ((((t109 * t138 * t156 * 2.0 + -t272) + t326) + -t398) + t432) + -t500;
    t13 = ((((-(t107 * t139 * t156 * 2.0) + -t270) + -t327) + t358) + t438) + -t491;
    t12 = ((((-(t99 * t142 * t156 * 2.0) + -t273) + t326) + t362) + t467) + -t500;
    t11 = ((((-(t99 * t140 * t156 * 2.0) + t269) + -t329) + -t367) + t455) + t503;
    t16 = ((((t107 * t138 * t156 * 2.0 + t268) + -t329) + -t386) + -t434) + t503;
    t15 = ((((-(t100 * t141 * t156 * 2.0) + t273) + -t325) + t399) + t431) + -t499;
    t14 = ((((t108 * t139 * t156 * 2.0 + t272) + -t325) + -t360) + t468) + -t499;
    t105 = ((((-(t139 * t140 * t156 * 2.0) + t275) + t327) + t427) + -t438) + t491;
    t103 = ((((t138 * t140 * t156 * 2.0 + -t274) + t329) + t434) + -t455) + -t503;
    t101 = ((((-(t137 * t142 * t156 * 2.0) + -t275) + t337) + t428) + -t437) + t492;
    t58 = ((((t139 * t141 * t156 * 2.0 + -t276) + t325) + -t431) + -t468) + t499;
    t57 = ((((t137 * t141 * t156 * 2.0 + t274) + t339) + t433) + -t456) + -t504;
    t56 = ((((t138 * t142 * t156 * 2.0 + t276) + t335) + -t432) + -t467) + t500;
    t55 = -t315 + t342;
    H[0] = (t55 + t142 * t142 * t156 * 2.0) - t521_tmp * 4.0;
    H[1] = t53;
    H[2] = t148;
    H[3] = t521;
    H[4] = t38;
    H[5] = t49;
    H[6] = t157;
    H[7] = t56;
    H[8] = t101;
    H[9] = t536;
    H[10] = t12;
    H[11] = t36;
    H[12] = t53;
    t54 = -t314 + t348;
    H[13] = (t54 + t141 * t141 * t156 * 2.0) - t524_tmp * 4.0;
    H[14] = t544;
    H[15] = t65;
    H[16] = t524;
    H[17] = t39;
    H[18] = t58;
    H[19] = t147;
    H[20] = t57;
    H[21] = t15;
    H[22] = t538;
    H[23] = t27;
    H[24] = t148;
    H[25] = t544;
    t53 = -t313 + t343;
    H[26] = (t53 + t140 * t140 * t156 * 2.0) + t520_tmp * 4.0;
    H[27] = t51;
    H[28] = t37;
    H[29] = t520;
    H[30] = t105;
    H[31] = t103;
    H[32] = t159;
    H[33] = t28;
    H[34] = t11;
    H[35] = t534;
    H[36] = t521;
    H[37] = t65;
    H[38] = t51;
    H[39] = (t55 + t109 * t109 * t156 * 2.0) + b_t521_tmp * 4.0;
    H[40] = t530;
    H[41] = t526;
    H[42] = t537;
    H[43] = t26;
    H[44] = t35;
    H[45] = t532;
    H[46] = t143;
    H[47] = t550;
    H[48] = t38;
    H[49] = t524;
    H[50] = t37;
    H[51] = t530;
    H[52] = (t54 + t108 * t108 * t156 * 2.0) + b_t524_tmp * 4.0;
    H[53] = t529;
    H[54] = t14;
    H[55] = t539;
    H[56] = t568;
    H[57] = t144;
    H[58] = t533;
    H[59] = t52;
    H[60] = t49;
    H[61] = t39;
    H[62] = t520;
    H[63] = t526;
    H[64] = t529;
    H[65] = (t53 + t107 * t107 * t156 * 2.0) - b_t520_tmp * 4.0;
    H[66] = t13;
    H[67] = t16;
    H[68] = t535;
    H[69] = t146;
    H[70] = t145;
    H[71] = t531;
    H[72] = t157;
    H[73] = t58;
    H[74] = t105;
    H[75] = t537;
    H[76] = t14;
    H[77] = t13;
    t55 = -t312 + t340;
    H[78] = (t55 + t139 * t139 * t156 * 2.0) - t519_tmp * 4.0;
    H[79] = t543;
    H[80] = t540;
    H[81] = t519;
    H[82] = t40;
    H[83] = t50;
    H[84] = t56;
    H[85] = t147;
    H[86] = t103;
    H[87] = t26;
    H[88] = t539;
    H[89] = t16;
    H[90] = t543;
    t54 = -t311 + t345;
    H[91] = (t54 + t138 * t138 * t156 * 2.0) - t523_tmp * 4.0;
    H[92] = t542;
    H[93] = t73;
    H[94] = t523;
    H[95] = t59;
    H[96] = t101;
    H[97] = t57;
    H[98] = t159;
    H[99] = t35;
    H[100] = t568;
    H[101] = t535;
    H[102] = t540;
    H[103] = t542;
    t53 = -t310 + t341;
    H[104] = (t53 + t137 * t137 * t156 * 2.0) + t522_tmp * 4.0;
    H[105] = t48;
    H[106] = t47;
    H[107] = t522;
    H[108] = t536;
    H[109] = t15;
    H[110] = t28;
    H[111] = t532;
    H[112] = t144;
    H[113] = t146;
    H[114] = t519;
    H[115] = t73;
    H[116] = t48;
    H[117] = (t55 + t100 * t100 * t156 * 2.0) - b_t519_tmp * 4.0;
    H[118] = t528;
    H[119] = t525;
    H[120] = t12;
    H[121] = t538;
    H[122] = t11;
    H[123] = t143;
    H[124] = t533;
    H[125] = t145;
    H[126] = t40;
    H[127] = t523;
    H[128] = t47;
    H[129] = t528;
    H[130] = (t54 + t99 * t99 * t156 * 2.0) - b_t523_tmp * 4.0;
    H[131] = t527;
    H[132] = t36;
    H[133] = t27;
    H[134] = t534;
    H[135] = t550;
    H[136] = t52;
    H[137] = t531;
    H[138] = t50;
    H[139] = t59;
    H[140] = t522;
    H[141] = t525;
    H[142] = t527;
    H[143] = (t53 + t98 * t98 * t156 * 2.0) + b_t522_tmp * 4.0;
}

inline void H_EE(const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    const Eigen::RowVector3d& v3,
    Eigen::Matrix<double, 12, 12>& H)
{
    H_EE(v0[0], v0[1], v0[2],
        v1[0], v1[1], v1[2],
        v2[0], v2[1], v2[2],
        v3[0], v3[1], v3[2],
        H.data());
}

inline void derivTest_EE(void)
{
    Eigen::RowVector3d v[4];
    v[0] << 0.0, 0.0, 0.0;
    v[1] << 1.0, 0.1, 0.0;
    v[2] << 0.0, 1.1, -0.1;
    v[3] << 0.0, 0.1, -1.1;

    double d0;
    d_EE(v[0], v[1], v[2], v[3], d0);

    double eps = 1.0e-6;
    Eigen::Matrix<double, 12, 1> fd, s;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            Eigen::RowVector3d vv[4];
            vv[0] = v[0];
            vv[1] = v[1];
            vv[2] = v[2];
            vv[3] = v[3];
            vv[i][j] += eps;
            double d;
            d_EE(vv[0], vv[1], vv[2], vv[3], d);
            fd[i * 3 + j] = (d - d0) / eps;
        }
    }

    g_EE(v[0], v[1], v[2], v[3], s);

    std::cout << s.transpose() << std::endl;
    std::cout << fd.transpose() << std::endl;
    std::cout << "relError = " << (s - fd).norm() / fd.norm() << std::endl;

    Eigen::Matrix<double, 12, 12> Hfd, Hs;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            Eigen::RowVector3d vv[4];
            vv[0] = v[0];
            vv[1] = v[1];
            vv[2] = v[2];
            vv[3] = v[3];
            vv[i][j] += eps;
            Eigen::Matrix<double, 12, 1> g;
            g_EE(vv[0], vv[1], vv[2], vv[3], g);
            Hfd.col(i * 3 + j) = (g - s) / eps;
        }
    }

    H_EE(v[0], v[1], v[2], v[3], Hs);

    std::cout << Hs.transpose() << std::endl;
    std::cout << Hfd.transpose() << std::endl;
    std::cout << "relError = " << (Hs - Hfd).norm() / Hfd.norm() << std::endl;
}

// http://geomalgorithms.com/a07-_distance.html
inline int dType_EE(const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    const Eigen::RowVector3d& v3)
{
    Eigen::RowVector3d u = v1 - v0;
    Eigen::RowVector3d v = v3 - v2;
    Eigen::RowVector3d w = v0 - v2;
    double a = u.squaredNorm(); // always >= 0
    double b = u.dot(v);
    double c = v.squaredNorm(); // always >= 0
    double d = u.dot(w);
    double e = v.dot(w);
    double D = a * c - b * b; // always >= 0
    double tD = D; // tc = tN / tD, default tD = D >= 0
    double sN, tN;

    int defaultCase = 8;

    // compute the line parameters of the two closest points
    // const double SMALL_NUM = 0.0;
    // if (D <= SMALL_NUM) { // the lines are almost parallel
    //     tN = e;
    //     tD = c;
    //     defaultCase = 2;
    // }
    // else { // get the closest points on the infinite lines
    sN = (b * e - c * d);
    if (sN <= 0.0) { // sc < 0 => the s=0 edge is visible
        tN = e;
        tD = c;
        defaultCase = 2;
    }
    else if (sN >= D) { // sc > 1  => the s=1 edge is visible
        tN = e + b;
        tD = c;
        defaultCase = 5;
    }
    else {
        tN = (a * e - b * d);
        if (tN > 0.0 && tN < tD && (u.cross(v).dot(w) == 0.0 || u.cross(v).squaredNorm() < 1.0e-20 * a * c)) {
            // if (tN > 0.0 && tN < tD && (u.cross(v).dot(w) == 0.0 || u.cross(v).squaredNorm() == 0.0)) {
            // std::cout << u.cross(v).squaredNorm() / (a * c) << ": " << sN << " " << D << ", " << tN << " " << tD << std::endl;
            // avoid coplanar or nearly parallel EE
            if (sN < D / 2) {
                tN = e;
                tD = c;
                defaultCase = 2;
            }
            else {
                tN = e + b;
                tD = c;
                defaultCase = 5;
            }
        }
        // else defaultCase stays as 8
    }
    // }

    if (tN <= 0.0) { // tc < 0 => the t=0 edge is visible
        // recompute sc for this edge
        if (-d <= 0.0) {
            return 0;
        }
        else if (-d >= a) {
            return 3;
        }
        else {
            return 6;
        }
    }
    else if (tN >= tD) { // tc > 1  => the t=1 edge is visible
        // recompute sc for this edge
        if ((-d + b) <= 0.0) {
            return 1;
        }
        else if ((-d + b) >= a) {
            return 4;
        }
        else {
            return 7;
        }
    }

    return defaultCase;
}

inline int dType_PT(const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    const Eigen::RowVector3d& v3)
{
    Eigen::Matrix<double, 2, 3> basis;
    basis.row(0) = v2 - v1;
    basis.row(1) = v3 - v1;

    const Eigen::RowVector3d nVec = basis.row(0).cross(basis.row(1));

    Eigen::Matrix<double, 2, 3> param;

    basis.row(1) = basis.row(0).cross(nVec);
    param.col(0) = (basis * basis.transpose()).ldlt().solve(basis * (v0 - v1).transpose());
    if (param(0, 0) > 0.0 && param(0, 0) < 1.0 && param(1, 0) >= 0.0) {
        return 3; // PE v1v2
    }
    else {
        basis.row(0) = v3 - v2;

        basis.row(1) = basis.row(0).cross(nVec);
        param.col(1) = (basis * basis.transpose()).ldlt().solve(basis * (v0 - v2).transpose());
        if (param(0, 1) > 0.0 && param(0, 1) < 1.0 && param(1, 1) >= 0.0) {
            return 4; // PE v2v3
        }
        else {
            basis.row(0) = v1 - v3;

            basis.row(1) = basis.row(0).cross(nVec);
            param.col(2) = (basis * basis.transpose()).ldlt().solve(basis * (v0 - v3).transpose());
            if (param(0, 2) > 0.0 && param(0, 2) < 1.0 && param(1, 2) >= 0.0) {
                return 5; // PE v3v1
            }
            else {
                if (param(0, 0) <= 0.0 && param(0, 2) >= 1.0) {
                    return 0; // PP v1
                }
                else if (param(0, 1) <= 0.0 && param(0, 0) >= 1.0) {
                    return 1; // PP v2
                }
                else if (param(0, 2) <= 0.0 && param(0, 1) >= 1.0) {
                    return 2; // PP v3
                }
                else {
                    return 6; // PT
                }
            }
        }
    }
}

inline void checkDType(void)
{
    std::vector<std::vector<Eigen::RowVector3d>> testData;

    testData.resize(testData.size() + 1);
    testData.back().resize(4);
    testData.back()[0] << 0.0, 0.0, 1.0;
    testData.back()[1] << 0.0, 0.0, 0.0;
    testData.back()[2] << 1.0, 0.0, 0.0;
    testData.back()[3] << 0.0, 1.0, 0.0;

    testData.resize(testData.size() + 1);
    testData.back().resize(4);
    testData.back()[0] << 0.0, 0.0, 0.0;
    testData.back()[1] << 0.0, 0.0, 0.0;
    testData.back()[2] << 1.0, 0.0, 0.0;
    testData.back()[3] << 0.0, 1.0, 0.0;

    srand((unsigned int)time(0));
    for (int i = 0; i < 6; ++i) {
        testData.resize(testData.size() + 1);
        testData.back().resize(4);
        testData.back()[0].setRandom();
        testData.back()[1].setRandom();
        testData.back()[2].setRandom();
        testData.back()[3].setRandom();
    }

    int i = 0;
    for (const auto& testI : testData) {
        std::cout << ">>>>>> Test " << i << " >>>>>>" << std::endl;
        std::cout << testI[0] << std::endl;
        std::cout << testI[1] << std::endl;
        std::cout << testI[2] << std::endl;
        std::cout << testI[3] << std::endl;
        std::cout << "PT: " << dType_PT(testI[0], testI[1], testI[2], testI[3]) << std::endl;
        std::cout << "EE: " << dType_EE(testI[0], testI[1], testI[2], testI[3]) << std::endl;

        ++i;
    }
}

// http://geomalgorithms.com/a02-_lines.html
inline void computePointEdgeD(const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    double& d)
{
    Eigen::RowVector3d v = v2 - v1;
    Eigen::RowVector3d w = v0 - v1;

    double c1 = w.dot(v);
    if (c1 <= 0.0) {
        d_PP(v0, v1, d);
    }
    else {
        double c2 = v.squaredNorm();
        if (c2 <= c1) {
            d_PP(v0, v2, d);
        }
        else {
            double b = c1 / c2;
            d_PP(v0, v1 + b * v, d);
        }
    }
}

inline void computePointTriD(const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    const Eigen::RowVector3d& v3,
    double& d)
{
    switch (dType_PT(v0, v1, v2, v3)) {
    case 0: {
        d_PP(v0, v1, d);
        break;
    }

    case 1: {
        d_PP(v0, v2, d);
        break;
    }

    case 2: {
        d_PP(v0, v3, d);
        break;
    }

    case 3: {
        d_PE(v0, v1, v2, d);
        break;
    }

    case 4: {
        d_PE(v0, v2, v3, d);
        break;
    }

    case 5: {
        d_PE(v0, v3, v1, d);
        break;
    }

    case 6: {
        d_PT(v0, v1, v2, v3, d);
        break;
    }

    default:
        d = -1.0;
        break;
    }
}

inline void computeEdgeEdgeD(const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    const Eigen::RowVector3d& v3,
    double& d)
{
    switch (dType_EE(v0, v1, v2, v3)) {
    case 0: {
        d_PP(v0, v2, d);
        break;
    }

    case 1: {
        d_PP(v0, v3, d);
        break;
    }

    case 2: {
        d_PE(v0, v2, v3, d);
        break;
    }

    case 3: {
        d_PP(v1, v2, d);
        break;
    }

    case 4: {
        d_PP(v1, v3, d);
        break;
    }

    case 5: {
        d_PE(v1, v2, v3, d);
        break;
    }

    case 6: {
        d_PE(v2, v0, v1, d);
        break;
    }

    case 7: {
        d_PE(v3, v0, v1, d);
        break;
    }

    case 8: {
        d_EE(v0, v1, v2, v3, d);
        break;
    }

    default:
        d = -1.0;
        break;
    }
}

inline void checkEdgeEdgeD(void)
{
    Eigen::RowVector3d v0, v1, v20, v3t, v2t;
    v0 << 0.0, 0.0, 0.0;
    v1 << 1.0, 0.0, 0.0;
    v20 << -0.1, 1.0, 1.0;
    v3t << 1.0, 0.0, -1.0;
    v2t << 0.0, 0.0, -1.0;

    Eigen::RowVector3d v30 = v20 + 0.01 * (v3t - v20);
    double step = 0.01;
    for (double t = 0.0; t < 1.0; t += step) {
        double d;
        computeEdgeEdgeD(v0, v1, v20, v30 + t * (v3t - v30), d);
        std::cout << t << " " << d << std::endl;
    }
    for (double t = 1.0; t <= 2.0; t += step) {
        double d;
        computeEdgeEdgeD(v0, v1, v20 + (t - 1.0) * (v2t - v20), v3t, d);
        std::cout << t << " " << d << std::endl;
    }
}

// for fixing tiny step size issue
inline void computeEECrossSqNorm(
    const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    const Eigen::RowVector3d& v3,
    double& result)
{
    result = (v1 - v0).cross(v3 - v2).squaredNorm();
}

inline void computeEECrossSqNormGradient(double v01, double v02, double v03, double v11,
    double v12, double v13, double v21, double v22, double v23, double v31, double v32, double v33, double g[12])
{
    double t8;
    double t9;
    double t10;
    double t11;
    double t12;
    double t13;
    double t23;
    double t24;
    double t25;
    double t26;
    double t27;
    double t28;
    double t29;
    double t30;
    double t31;
    double t32;
    double t33;

    /* COMPUTEEECROSSSQNORMGRADIENT */
    /*     G = COMPUTEEECROSSSQNORMGRADIENT(V01,V02,V03,V11,V12,V13,V21,V22,V23,V31,V32,V33) */
    /*     This function was generated by the Symbolic Math Toolbox version 8.3. */
    /*     01-Nov-2019 16:54:23 */
    t8 = -v11 + v01;
    t9 = -v12 + v02;
    t10 = -v13 + v03;
    t11 = -v31 + v21;
    t12 = -v32 + v22;
    t13 = -v33 + v23;
    t23 = t8 * t12 + -(t9 * t11);
    t24 = t8 * t13 + -(t10 * t11);
    t25 = t9 * t13 + -(t10 * t12);
    t26 = t8 * t23 * 2.0;
    t27 = t9 * t23 * 2.0;
    t28 = t8 * t24 * 2.0;
    t29 = t10 * t24 * 2.0;
    t30 = t9 * t25 * 2.0;
    t31 = t10 * t25 * 2.0;
    t32 = t11 * t23 * 2.0;
    t33 = t12 * t23 * 2.0;
    t23 = t11 * t24 * 2.0;
    t10 = t13 * t24 * 2.0;
    t9 = t12 * t25 * 2.0;
    t8 = t13 * t25 * 2.0;
    g[0] = t33 + t10;
    g[1] = -t32 + t8;
    g[2] = -t23 - t9;
    g[3] = -t33 - t10;
    g[4] = t32 - t8;
    g[5] = t23 + t9;
    g[6] = -t27 - t29;
    g[7] = t26 - t31;
    g[8] = t28 + t30;
    g[9] = t27 + t29;
    g[10] = -t26 + t31;
    g[11] = -t28 - t30;
}

inline void computeEECrossSqNormGradient(
    const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    const Eigen::RowVector3d& v3,
    Eigen::Matrix<double, 12, 1>& grad)
{
    computeEECrossSqNormGradient(
        v0[0], v0[1], v0[2],
        v1[0], v1[1], v1[2],
        v2[0], v2[1], v2[2],
        v3[0], v3[1], v3[2],
        grad.data());
}

inline void computeEECrossSqNormHessian(double v01, double v02, double v03, double v11,
    double v12, double v13, double v21, double v22, double v23, double v31, double v32, double v33, double H[144])
{
    double t8;
    double t9;
    double t10;
    double t11;
    double t12;
    double t13;
    double t32;
    double t33;
    double t34;
    double t35;
    double t48;
    double t36;
    double t49;
    double t37;
    double t38;
    double t39;
    double t40;
    double t41;
    double t42;
    double t43;
    double t44;
    double t45;
    double t46;
    double t47;
    double t50;
    double t51;
    double t52;
    double t20;
    double t23;
    double t24;
    double t25;
    double t86;
    double t87;
    double t88;
    double t74;
    double t75;
    double t76;
    double t77;
    double t78;
    double t79;
    double t89;
    double t90;
    double t91;
    double t92;
    double t93;
    double t94;
    double t95;

    /* COMPUTEEECROSSSQNORMHESSIAN */
    /*     H = COMPUTEEECROSSSQNORMHESSIAN(V01,V02,V03,V11,V12,V13,V21,V22,V23,V31,V32,V33) */
    /*     This function was generated by the Symbolic Math Toolbox version 8.3. */
    /*     01-Nov-2019 16:54:23 */
    t8 = -v11 + v01;
    t9 = -v12 + v02;
    t10 = -v13 + v03;
    t11 = -v31 + v21;
    t12 = -v32 + v22;
    t13 = -v33 + v23;
    t32 = t8 * t9 * 2.0;
    t33 = t8 * t10 * 2.0;
    t34 = t9 * t10 * 2.0;
    t35 = t8 * t11 * 2.0;
    t48 = t8 * t12;
    t36 = t48 * 2.0;
    t49 = t9 * t11;
    t37 = t49 * 2.0;
    t38 = t48 * 4.0;
    t48 = t8 * t13;
    t39 = t48 * 2.0;
    t40 = t49 * 4.0;
    t41 = t9 * t12 * 2.0;
    t49 = t10 * t11;
    t42 = t49 * 2.0;
    t43 = t48 * 4.0;
    t48 = t9 * t13;
    t44 = t48 * 2.0;
    t45 = t49 * 4.0;
    t49 = t10 * t12;
    t46 = t49 * 2.0;
    t47 = t48 * 4.0;
    t48 = t49 * 4.0;
    t49 = t10 * t13 * 2.0;
    t50 = t11 * t12 * 2.0;
    t51 = t11 * t13 * 2.0;
    t52 = t12 * t13 * 2.0;
    t20 = t8 * t8 * 2.0;
    t9 = t9 * t9 * 2.0;
    t8 = t10 * t10 * 2.0;
    t23 = t11 * t11 * 2.0;
    t24 = t12 * t12 * 2.0;
    t25 = t13 * t13 * 2.0;
    t86 = t35 + t41;
    t87 = t35 + t49;
    t88 = t41 + t49;
    t74 = t20 + t9;
    t75 = t20 + t8;
    t76 = t9 + t8;
    t77 = t23 + t24;
    t78 = t23 + t25;
    t79 = t24 + t25;
    t89 = t40 + -t36;
    t90 = t36 + -t40;
    t91 = t37 + -t38;
    t92 = t38 + -t37;
    t93 = t45 + -t39;
    t94 = t39 + -t45;
    t95 = t42 + -t43;
    t37 = t43 + -t42;
    t39 = t48 + -t44;
    t45 = t44 + -t48;
    t38 = t46 + -t47;
    t40 = t47 + -t46;
    t36 = -t35 + -t41;
    t13 = -t35 + -t49;
    t11 = -t41 + -t49;
    t12 = -t20 + -t9;
    t10 = -t20 + -t8;
    t8 = -t9 + -t8;
    t9 = -t23 + -t24;
    t49 = -t23 + -t25;
    t48 = -t24 + -t25;
    H[0] = t79;
    H[1] = -t50;
    H[2] = -t51;
    H[3] = t48;
    H[4] = t50;
    H[5] = t51;
    H[6] = t11;
    H[7] = t92;
    H[8] = t37;
    H[9] = t88;
    H[10] = t91;
    H[11] = t95;
    H[12] = -t50;
    H[13] = t78;
    H[14] = -t52;
    H[15] = t50;
    H[16] = t49;
    H[17] = t52;
    H[18] = t89;
    H[19] = t13;
    H[20] = t40;
    H[21] = t90;
    H[22] = t87;
    H[23] = t38;
    H[24] = -t51;
    H[25] = -t52;
    H[26] = t77;
    H[27] = t51;
    H[28] = t52;
    H[29] = t9;
    H[30] = t93;
    H[31] = t39;
    H[32] = t36;
    H[33] = t94;
    H[34] = t45;
    H[35] = t86;
    H[36] = t48;
    H[37] = t50;
    H[38] = t51;
    H[39] = t79;
    H[40] = -t50;
    H[41] = -t51;
    H[42] = t88;
    H[43] = t91;
    H[44] = t95;
    H[45] = t11;
    H[46] = t92;
    H[47] = t37;
    H[48] = t50;
    H[49] = t49;
    H[50] = t52;
    H[51] = -t50;
    H[52] = t78;
    H[53] = -t52;
    H[54] = t90;
    H[55] = t87;
    H[56] = t38;
    H[57] = t89;
    H[58] = t13;
    H[59] = t40;
    H[60] = t51;
    H[61] = t52;
    H[62] = t9;
    H[63] = -t51;
    H[64] = -t52;
    H[65] = t77;
    H[66] = t94;
    H[67] = t45;
    H[68] = t86;
    H[69] = t93;
    H[70] = t39;
    H[71] = t36;
    H[72] = t11;
    H[73] = t89;
    H[74] = t93;
    H[75] = t88;
    H[76] = t90;
    H[77] = t94;
    H[78] = t76;
    H[79] = -t32;
    H[80] = -t33;
    H[81] = t8;
    H[82] = t32;
    H[83] = t33;
    H[84] = t92;
    H[85] = t13;
    H[86] = t39;
    H[87] = t91;
    H[88] = t87;
    H[89] = t45;
    H[90] = -t32;
    H[91] = t75;
    H[92] = -t34;
    H[93] = t32;
    H[94] = t10;
    H[95] = t34;
    H[96] = t37;
    H[97] = t40;
    H[98] = t36;
    H[99] = t95;
    H[100] = t38;
    H[101] = t86;
    H[102] = -t33;
    H[103] = -t34;
    H[104] = t74;
    H[105] = t33;
    H[106] = t34;
    H[107] = t12;
    H[108] = t88;
    H[109] = t90;
    H[110] = t94;
    H[111] = t11;
    H[112] = t89;
    H[113] = t93;
    H[114] = t8;
    H[115] = t32;
    H[116] = t33;
    H[117] = t76;
    H[118] = -t32;
    H[119] = -t33;
    H[120] = t91;
    H[121] = t87;
    H[122] = t45;
    H[123] = t92;
    H[124] = t13;
    H[125] = t39;
    H[126] = t32;
    H[127] = t10;
    H[128] = t34;
    H[129] = -t32;
    H[130] = t75;
    H[131] = -t34;
    H[132] = t95;
    H[133] = t38;
    H[134] = t86;
    H[135] = t37;
    H[136] = t40;
    H[137] = t36;
    H[138] = t33;
    H[139] = t34;
    H[140] = t12;
    H[141] = -t33;
    H[142] = -t34;
    H[143] = t74;
}

inline void computeEECrossSqNormHessian(
    const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    const Eigen::RowVector3d& v3,
    Eigen::Matrix<double, 12, 12>& Hessian)
{
    computeEECrossSqNormHessian(
        v0[0], v0[1], v0[2],
        v1[0], v1[1], v1[2],
        v2[0], v2[1], v2[2],
        v3[0], v3[1], v3[2],
        Hessian.data());
}

inline void derivTest_EECross(void)
{
    Eigen::RowVector3d v[4];
    v[0] << 0.0, 0.0, 0.0;
    v[1] << 1.0, 0.1, 0.0;
    v[2] << 0.0, 1.1, -0.1;
    v[3] << 0.0, 0.1, -1.1;

    double f0;
    computeEECrossSqNorm(v[0], v[1], v[2], v[3], f0);

    double eps = 1.0e-6;
    Eigen::Matrix<double, 12, 1> fd, s;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            Eigen::RowVector3d vv[4];
            vv[0] = v[0];
            vv[1] = v[1];
            vv[2] = v[2];
            vv[3] = v[3];
            vv[i][j] += eps;
            double f;
            computeEECrossSqNorm(vv[0], vv[1], vv[2], vv[3], f);
            fd[i * 3 + j] = (f - f0) / eps;
        }
    }

    computeEECrossSqNormGradient(v[0], v[1], v[2], v[3], s);

    std::cout << s.transpose() << std::endl;
    std::cout << fd.transpose() << std::endl;
    std::cout << "relError = " << (s - fd).norm() / fd.norm() << std::endl;

    Eigen::Matrix<double, 12, 12> Hfd, Hs;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            Eigen::RowVector3d vv[4];
            vv[0] = v[0];
            vv[1] = v[1];
            vv[2] = v[2];
            vv[3] = v[3];
            vv[i][j] += eps;
            Eigen::Matrix<double, 12, 1> g;
            computeEECrossSqNormGradient(vv[0], vv[1], vv[2], vv[3], g);
            Hfd.col(i * 3 + j) = (g - s) / eps;
        }
    }

    computeEECrossSqNormHessian(v[0], v[1], v[2], v[3], Hs);

    std::cout << Hs.transpose() << std::endl;
    std::cout << Hfd.transpose() << std::endl;
    std::cout << "relError = " << (Hs - Hfd).norm() / Hfd.norm() << std::endl;
}

inline void compute_q(double input, double eps_x, double& e)
{
    double input_div_eps_x = input / eps_x;
    e = (-input_div_eps_x + 2.0) * input_div_eps_x;
}

inline void compute_q_g(double input, double eps_x, double& g)
{
    double one_div_eps_x = 1.0 / eps_x;
    g = 2.0 * one_div_eps_x * (-one_div_eps_x * input + 1.0);
}

inline void compute_q_H(double input, double eps_x, double& H)
{
    H = -2.0 / (eps_x * eps_x);
}

inline void compute_e(
    const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    const Eigen::RowVector3d& v3,
    double eps_x, double& e)
{
    double EECrossSqNorm;
    computeEECrossSqNorm(v0, v1, v2, v3, EECrossSqNorm);
    if (EECrossSqNorm < eps_x) {
        compute_q(EECrossSqNorm, eps_x, e);
    }
    else {
        e = 1.0;
    }
}

inline void compute_e_g(
    const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    const Eigen::RowVector3d& v3,
    double eps_x, Eigen::Matrix<double, 12, 1>& g)
{
    double EECrossSqNorm;
    computeEECrossSqNorm(v0, v1, v2, v3, EECrossSqNorm);
    if (EECrossSqNorm < eps_x) {
        double q_g;
        compute_q_g(EECrossSqNorm, eps_x, q_g);
        computeEECrossSqNormGradient(v0, v1, v2, v3, g);
        g *= q_g;
    }
    else {
        g.setZero();
    }
}

inline void compute_e_H(
    const Eigen::RowVector3d& v0,
    const Eigen::RowVector3d& v1,
    const Eigen::RowVector3d& v2,
    const Eigen::RowVector3d& v3,
    double eps_x, Eigen::Matrix<double, 12, 12>& H)
{
    double EECrossSqNorm;
    computeEECrossSqNorm(v0, v1, v2, v3, EECrossSqNorm);
    if (EECrossSqNorm < eps_x) {
        double q_g, q_H;
        compute_q_g(EECrossSqNorm, eps_x, q_g);
        compute_q_H(EECrossSqNorm, eps_x, q_H);

        Eigen::Matrix<double, 12, 1> g;
        computeEECrossSqNormGradient(v0, v1, v2, v3, g);

        computeEECrossSqNormHessian(v0, v1, v2, v3, H);
        H *= q_g;
        H += (q_H * g) * g.transpose();
    }
    else {
        H.setZero();
    }
}

inline void derivTest_e(double eps_x = 10.0)
{
    Eigen::RowVector3d v[4];
    v[0] << 0.0, 0.0, 0.0;
    v[1] << 1.0, 0.1, 0.0;
    v[2] << 0.0, 1.1, -0.1;
    v[3] << 0.0, 0.1, -1.1;

    double f0;
    compute_e(v[0], v[1], v[2], v[3], eps_x, f0);

    double eps = 1.0e-6;
    Eigen::Matrix<double, 12, 1> fd, s;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            Eigen::RowVector3d vv[4];
            vv[0] = v[0];
            vv[1] = v[1];
            vv[2] = v[2];
            vv[3] = v[3];
            vv[i][j] += eps;
            double f;
            compute_e(vv[0], vv[1], vv[2], vv[3], eps_x, f);
            fd[i * 3 + j] = (f - f0) / eps;
        }
    }

    compute_e_g(v[0], v[1], v[2], v[3], eps_x, s);

    std::cout << s.transpose() << std::endl;
    std::cout << fd.transpose() << std::endl;
    std::cout << "relError = " << (s - fd).norm() / fd.norm() << std::endl;

    Eigen::Matrix<double, 12, 12> Hfd, Hs;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            Eigen::RowVector3d vv[4];
            vv[0] = v[0];
            vv[1] = v[1];
            vv[2] = v[2];
            vv[3] = v[3];
            vv[i][j] += eps;
            Eigen::Matrix<double, 12, 1> g;
            compute_e_g(vv[0], vv[1], vv[2], vv[3], eps_x, g);
            Hfd.col(i * 3 + j) = (g - s) / eps;
        }
    }

    compute_e_H(v[0], v[1], v[2], v[3], eps_x, Hs);

    std::cout << Hs.transpose() << std::endl;
    std::cout << Hfd.transpose() << std::endl;
    std::cout << "relError = " << (Hs - Hfd).norm() / Hfd.norm() << std::endl;
}

template <int dim>
inline void compute_eps_x(const Mesh<dim>& mesh,
    int eI0, int eI1, int eJ0, int eJ1, double& eps_x)
{
    eps_x = 1.0e-3 * (mesh.V_rest.row(eI0) - mesh.V_rest.row(eI1)).squaredNorm() * (mesh.V_rest.row(eJ0) - mesh.V_rest.row(eJ1)).squaredNorm();
}

template <int dim>
inline void compute_eps_x(const Mesh<dim>& mesh, const Eigen::MatrixXd& V_MCO,
    int eI0, int eI1, int eJ0, int eJ1, double& eps_x)
{
    eps_x = 1.0e-3 * (mesh.V_rest.row(eI0) - mesh.V_rest.row(eI1)).squaredNorm() * (V_MCO.row(eJ0) - V_MCO.row(eJ1)).squaredNorm();
}

// broad phase collision detection

// inline bool broadPhaseCollision(const Eigen::Matrix<double, Eigen::Dynamic, 3>& obj0,
//     const Eigen::Matrix<double, Eigen::Dynamic, 3>& obj1)
// {
//     Eigen::Matrix<double, 1, 3> min0 = obj0.colwise().minCoeff();
//     Eigen::Matrix<double, 1, 3> max0 = obj0.colwise().maxCoeff();
//     Eigen::Matrix<double, 1, 3> min1 = obj1.colwise().minCoeff();
//     Eigen::Matrix<double, 1, 3> max1 = obj1.colwise().maxCoeff();
//     return ((min0 - max1).array() < 0).all() && ((min1 - max0).array() < 0).all();
// }

} // namespace IPC

#endif /* MeshCollisionUtils_hpp */
