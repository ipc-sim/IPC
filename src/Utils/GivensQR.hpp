//
//  GivensQR.hpp
//  IPC
//
//  Created by Minchen Li on 10/19/19.
//

#include "Types.hpp"

#include <Eigen/Eigen>

namespace IPC {

template <int dim>
class GivensQR {
public:
    static void compute(const Eigen::Matrix<double, dim, dim>& A,
        Eigen::Matrix<double, dim, dim>& Q,
        Eigen::Matrix<double, dim, dim>& R)
    {
        R = A;
        Q.setIdentity();
        for (int j = 0; j < dim; ++j) {
            for (int i = dim - 1; i > j; --i) {
                double a = R(i - 1, j), b = R(i, j);
                double d = a * a + b * b;
                double c = 1, s = 0;
                double sqrtd = std::sqrt(d);
                if (sqrtd) {
                    double t = 1.0 / sqrtd;
                    c = a * t;
                    s = -b * t;
                }

                for (int k = 0; k < dim; ++k) {
                    double tau1 = R(i - 1, k);
                    double tau2 = R(i, k);
                    R(i - 1, k) = c * tau1 - s * tau2;
                    R(i, k) = s * tau1 + c * tau2;
                }

                for (int k = 0; k < dim; ++k) {
                    double tau1 = Q(i - 1, k);
                    double tau2 = Q(i, k);
                    Q(i - 1, k) = c * tau1 - s * tau2;
                    Q(i, k) = s * tau1 + c * tau2;
                }
            }
        }
        Q.transposeInPlace();
    }
};

} // namespace IPC