//
//  ClosedFormSVDSVD2d.hpp
//  IPC
//
//  Created by Minchen Li on 8/31/18.
//  based on https://www.researchgate.net/publication/263580188_Closed_Form_SVD_Solutions_for_2_x_2_Matrices_-_Rev_2
//

#ifndef ClosedFormSVD2d_hpp
#define ClosedFormSVD2d_hpp

#include <Eigen/Eigen>

#include <iostream>

namespace IPC {

template <typename MatrixType>
class AutoFlipSVD : Eigen::JacobiSVD<MatrixType> {
    typedef Eigen::JacobiSVD<MatrixType> Base;

public:
    AutoFlipSVD(void) {}
    AutoFlipSVD(const MatrixType& mtr, unsigned int computationOptions = 0)
        : Base(2, 2, computationOptions)
    {
        if (MatrixType::RowsAtCompileTime == Eigen::Dynamic) {
            assert(mtr.rows() == 2);
        }
        else {
            assert(MatrixType::RowsAtCompileTime == 2);
        }

        if (MatrixType::ColsAtCompileTime == Eigen::Dynamic) {
            assert(mtr.cols() == 2);
        }
        else {
            assert(MatrixType::ColsAtCompileTime == 2);
        }

        Base::m_isInitialized = true;

        compute(mtr);
    }

public:
    AutoFlipSVD& compute(const MatrixType& A)
    {
        bool computeU = (Base::m_computeFullU || Base::m_computeThinU);
        bool computeV = (Base::m_computeFullV || Base::m_computeThinV);

        const double a = A(0, 0);
        const double b = A(1, 0);
        const double c = A(0, 1);
        const double d = A(1, 1);

        const double ad = a * d;
        const double bc = b * c;
        const double _2admbc = 2.0 * (ad - bc);
        const double sqn = A.squaredNorm();

        const double sum = sqn + _2admbc;
        const double dif = sqn - _2admbc;
        if (dif <= 0.0) {
            // avoid dividing by 0 in general formula
            const double aa = (a + d) / 2.0;
            const double bb = (b - c) / 2.0;
            const double lambda = std::sqrt(aa * aa + bb * bb);
            Base::m_singularValues.setConstant(lambda);

            if (computeU) {
                if (lambda == 0.0) {
                    // avoid dividing by 0
                    Base::m_matrixU.setIdentity();
                }
                else {
                    const double cosl = aa / lambda;
                    const double sinl = bb / lambda;
                    Base::m_matrixU << cosl, -sinl, sinl, cosl;
                }
            }

            if (computeV) {
                Base::m_matrixV.setIdentity();
            }
        }
        else if (sum <= 0.0) {
            // avoid dividing by 0 in general formula
            // symmetric matrix with a=-d
            const double aa = (a - d) / 2.0;
            const double bb = (b + c) / 2.0;
            const double lambda = std::sqrt(aa * aa + bb * bb);
            Base::m_singularValues << lambda, -lambda;

            if (computeU || computeV) {
                if (bb == 0.0) {
                    // avoid dividing by 0 and sqrt(<0)
                    if (computeU) {
                        Base::m_matrixU.setIdentity();
                    }

                    if (computeV) {
                        Base::m_matrixV.setIdentity();
                    }
                }
                else {
                    const double a_div_lambda_half = aa / lambda / 2.0;
                    bool neg_b = (bb < 0.0);
                    const double cos2 = 0.5 + a_div_lambda_half;
                    const double cos = ((cos2 <= 0.0) ? 0.0 : std::sqrt(cos2));
                    const double sin2 = 0.5 - a_div_lambda_half;
                    const double sin = ((sin2 <= 0.0) ? 0.0 : (neg_b ? -std::sqrt(sin2) : std::sqrt(sin2)));
                    if (computeU) {
                        Base::m_matrixU << cos, -sin, sin, cos;
                    }

                    if (computeV) {
                        Base::m_matrixV << cos, -sin, sin, cos;
                    }
                }
            }
        }
        else {
            const double sqrt_sum = std::sqrt(sum); // safe
            const double sqrt_dif = std::sqrt(dif); // safe

            Base::m_singularValues[0] = (sqrt_sum + sqrt_dif) / 2.0;
            Base::m_singularValues[1] = ((_2admbc < 0.0) ? (-std::abs(sqrt_sum - sqrt_dif) / 2.0) : (std::abs(sqrt_sum - sqrt_dif) / 2.0));

            if (computeU || computeV) {
                const double a2 = a * a;
                const double b2 = b * b;
                const double c2 = c * c;
                const double d2 = d * d;

                const double denom = sqrt_sum * sqrt_dif * 2.0;

                const double a2md2 = a2 - d2;
                const double b2mc2 = b2 - c2;

                const double ab = a * b;
                const double cd = c * d;
                const bool neg_ab_p_cd = ((ab + cd) < 0);

                if (computeU) {
                    const double a2md2_m_b2mc2_div_ = (a2md2 - b2mc2) / denom; // safe

                    // avoid sqrt(<0)
                    const double cosl2 = 0.5 + a2md2_m_b2mc2_div_;
                    const double cosl = ((cosl2 <= 0.0) ? 0.0 : std::sqrt(cosl2));
                    const double sinl2 = 0.5 - a2md2_m_b2mc2_div_;
                    const double sinl = ((sinl2 <= 0.0) ? 0.0 : (neg_ab_p_cd ? -std::sqrt(sinl2) : std::sqrt(sinl2)));

                    Base::m_matrixU << cosl, -sinl, sinl, cosl;
                }

                if (computeV) {
                    const double ac = a * c;
                    const double bd = b * d;
                    const bool neg_ac_p_bd = ((ac + bd) < 0);
                    const double a2md2_p_b2mc2_div_ = (a2md2 + b2mc2) / denom; // safe

                    // avoid sqrt(<0)
                    const double cosr2 = 0.5 + a2md2_p_b2mc2_div_;
                    const double cosr = ((cosr2 <= 0.0) ? 0.0 : std::sqrt(cosr2));
                    const double sinr2 = 0.5 - a2md2_p_b2mc2_div_;
                    const double sinr = ((sinr2 <= 0.0) ? 0.0 : (neg_ac_p_bd ? -std::sqrt(sinr2) : std::sqrt(sinr2)));

                    const bool s = neg_ab_p_cd ^ neg_ac_p_bd;
                    const bool neg_apsd = ((a + (s ? -d : d)) < 0.0);
                    if (neg_apsd) {
                        Base::m_matrixV << -cosr, sinr, -sinr, -cosr;
                    }
                    else {
                        Base::m_matrixV << cosr, -sinr, sinr, cosr;
                    }
                }
            }
        }

        return *this;
    }
    AutoFlipSVD& compute(const MatrixType& mtr, unsigned int computationOptions)
    {
        if (MatrixType::RowsAtCompileTime == Eigen::Dynamic) {
            assert(mtr.rows() == 2);
        }
        else {
            assert(MatrixType::RowsAtCompileTime == 2);
        }

        if (MatrixType::ColsAtCompileTime == Eigen::Dynamic) {
            assert(mtr.cols() == 2);
        }
        else {
            assert(MatrixType::ColsAtCompileTime == 2);
        }

        allocate(computationOptions);
        Base::m_isInitialized = true;

        compute(mtr);

        return *this;
    }

protected:
    void allocate(unsigned int computationOptions)
    {
        if (Base::m_isAllocated && 2 == Base::m_rows && 2 == Base::m_cols && computationOptions == Base::m_computationOptions) {
            return;
        }

        Base::m_rows = 2;
        Base::m_cols = 2;
        Base::m_isInitialized = false;
        Base::m_isAllocated = true;
        Base::m_computationOptions = computationOptions;
        Base::m_computeFullU = (computationOptions & Eigen::ComputeFullU) != 0;
        Base::m_computeThinU = (computationOptions & Eigen::ComputeThinU) != 0;
        Base::m_computeFullV = (computationOptions & Eigen::ComputeFullV) != 0;
        Base::m_computeThinV = (computationOptions & Eigen::ComputeThinV) != 0;
        eigen_assert(!(Base::m_computeFullU && Base::m_computeThinU) && "JacobiSVD: you can't ask for both full and thin U");
        eigen_assert(!(Base::m_computeFullV && Base::m_computeThinV) && "JacobiSVD: you can't ask for both full and thin V");
        eigen_assert(EIGEN_IMPLIES(Base::m_computeThinU || Base::m_computeThinV,
                         MatrixType::ColsAtCompileTime == Eigen::Dynamic)
            && "JacobiSVD: thin U and V are only available when your matrix has a dynamic number of columns.");

        Base::m_diagSize = 2;
        Base::m_singularValues.resize(2);
        Base::m_matrixU.resize(2, 2);
        Base::m_matrixV.resize(2, 2);
        Base::m_workMatrix.resize(2, 2);
    }

public:
    const typename Eigen::JacobiSVD<MatrixType>::SingularValuesType& singularValues(void) const
    {
        return Base::m_singularValues;
    }
    const MatrixType& matrixU(void) const
    {
        return Base::m_matrixU;
    }
    const MatrixType& matrixV(void) const
    {
        return Base::m_matrixV;
    }
};

} // namespace IPC

#endif /* ClosedFormSVD2d_hpp */
