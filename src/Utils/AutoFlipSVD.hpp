//
//  AutoFlipSVD.hpp
//  IPC
//
//  Created by Minchen Li on 6/21/18.
//

#ifndef AutoFlipSVD_hpp
#define AutoFlipSVD_hpp

#include "ImplicitQRSVD.h"

#include <Eigen/Eigen>

#include <iostream>

namespace IPC {

template <typename MatrixType>
class AutoFlipSVD : Eigen::JacobiSVD<MatrixType> {
protected:
    bool flipped_U, flipped_V, flipped_sigma;

    typename Eigen::JacobiSVD<MatrixType>::SingularValuesType singularValues_flipped;
    MatrixType matrixU_flipped, matrixV_flipped;

public:
    AutoFlipSVD(void) {}
    AutoFlipSVD(const MatrixType& mtr, unsigned int computationOptions = 0)
    {
        compute(mtr, computationOptions);
    }

public:
    template <int dim = MatrixType::RowsAtCompileTime>
    typename std::enable_if<dim == 3, AutoFlipSVD<MatrixType>>::type&
    compute(const MatrixType& mtr, unsigned int computationOptions)
    {
        flipped_U = flipped_V = flipped_sigma = true;
#ifdef USE_IQRSVD
        JIXIE::singularValueDecomposition(mtr,
            matrixU_flipped,
            singularValues_flipped,
            matrixV_flipped);
#else
        if ((computationOptions & Eigen::ComputeFullU) || (computationOptions & Eigen::ComputeFullV)) {
            fastSVD3d(mtr, matrixU_flipped, singularValues_flipped, matrixV_flipped);
        }
        else {
            fastComputeSingularValues3d(mtr, singularValues_flipped);
        }
#endif
        return *this;
    }
    template <int dim = MatrixType::RowsAtCompileTime>
    typename std::enable_if<dim == 2, AutoFlipSVD<MatrixType>>::type&
    compute(const MatrixType& mtr, unsigned int computationOptions)
    {
#ifdef USE_IQRSVD
        flipped_U = flipped_V = flipped_sigma = true;
        JIXIE::singularValueDecomposition(mtr,
            matrixU_flipped,
            singularValues_flipped,
            matrixV_flipped);
#else
        flipped_U = flipped_V = flipped_sigma = false;
        Eigen::JacobiSVD<MatrixType>::compute(mtr, computationOptions);
        flip2d(mtr, computationOptions);
#endif
        return *this;
    }

    void set(const Eigen::Matrix3d& U, const Eigen::Vector3d& Sigma, const Eigen::Matrix3d& V)
    {
        flipped_U = true, flipped_V = true, flipped_sigma = true;
        matrixU_flipped = U;
        singularValues_flipped = Sigma;
        matrixV_flipped = V;
    }

protected:
    void flip2d(const MatrixType& mtr, unsigned int computationOptions)
    {
        //!!! this flip algorithm is only valid in 2D
        bool fullUComputed = (computationOptions & Eigen::ComputeFullU);
        bool fullVComputed = (computationOptions & Eigen::ComputeFullV);
        if (fullUComputed && fullVComputed) {
            if (Eigen::JacobiSVD<MatrixType>::m_matrixU.determinant() < 0.0) {
                matrixU_flipped = Eigen::JacobiSVD<MatrixType>::m_matrixU;
                matrixU_flipped.col(1) *= -1.0;
                flipped_U = true;

                if (!flipped_sigma) {
                    singularValues_flipped = Eigen::JacobiSVD<MatrixType>::m_singularValues;
                }
                singularValues_flipped[1] *= -1.0;
                flipped_sigma = true;
            }
            if (Eigen::JacobiSVD<MatrixType>::m_matrixV.determinant() < 0.0) {
                matrixV_flipped = Eigen::JacobiSVD<MatrixType>::m_matrixV;
                matrixV_flipped.col(1) *= -1.0;
                flipped_V = true;

                if (!flipped_sigma) {
                    singularValues_flipped = Eigen::JacobiSVD<MatrixType>::m_singularValues;
                }
                singularValues_flipped[1] *= -1.0;
                flipped_sigma = true;
            }
        }
        else if (mtr.determinant() < 0.0) {
            singularValues_flipped = Eigen::JacobiSVD<MatrixType>::m_singularValues;
            singularValues_flipped[1] *= -1.0;
            flipped_sigma = true;
        }

        if (std::isnan(singularValues()[0]) || std::isnan(singularValues()[1])) {
            // degenerated case
            singularValues_flipped.setZero();
            flipped_sigma = true;
            if (fullUComputed && fullVComputed) {
                matrixU_flipped.setIdentity();
                matrixV_flipped.setIdentity();
                flipped_U = flipped_V = true;
            }
        }
    }

    //TODO: merge with IglUtils::computeCofactorMtr
    template <int dim>
    void computeCofactorMtr(const Eigen::Matrix<double, dim, dim>& F,
        Eigen::Matrix<double, dim, dim>& A)
    {
        switch (dim) {
        case 2:
            A(0, 0) = F(1, 1);
            A(0, 1) = -F(1, 0);
            A(1, 0) = -F(0, 1);
            A(1, 1) = F(0, 0);
            break;

        case 3:
            A(0, 0) = F(1, 1) * F(2, 2) - F(1, 2) * F(2, 1);
            A(0, 1) = F(1, 2) * F(2, 0) - F(1, 0) * F(2, 2);
            A(0, 2) = F(1, 0) * F(2, 1) - F(1, 1) * F(2, 0);
            A(1, 0) = F(0, 2) * F(2, 1) - F(0, 1) * F(2, 2);
            A(1, 1) = F(0, 0) * F(2, 2) - F(0, 2) * F(2, 0);
            A(1, 2) = F(0, 1) * F(2, 0) - F(0, 0) * F(2, 1);
            A(2, 0) = F(0, 1) * F(1, 2) - F(0, 2) * F(1, 1);
            A(2, 1) = F(0, 2) * F(1, 0) - F(0, 0) * F(1, 2);
            A(2, 2) = F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0);
            break;

        default:
            assert(0 && "dim not 2 or 3");
            break;
        }
    }
    void fastEigenvalues(const Eigen::Matrix3d& A_Sym,
        Eigen::Vector3d& lambda)
    // 24 mults, 20 adds, 1 atan2, 1 sincos, 2 sqrts
    {
        using T = double;
        using std::max;
        using std::swap;
        T m = ((T)1 / 3) * (A_Sym(0, 0) + A_Sym(1, 1) + A_Sym(2, 2));
        T a00 = A_Sym(0, 0) - m;
        T a11 = A_Sym(1, 1) - m;
        T a22 = A_Sym(2, 2) - m;
        T a12_sqr = A_Sym(0, 1) * A_Sym(0, 1);
        T a13_sqr = A_Sym(0, 2) * A_Sym(0, 2);
        T a23_sqr = A_Sym(1, 2) * A_Sym(1, 2);
        T p = ((T)1 / 6) * (a00 * a00 + a11 * a11 + a22 * a22 + 2 * (a12_sqr + a13_sqr + a23_sqr));
        T q = (T).5 * (a00 * (a11 * a22 - a23_sqr) - a11 * a13_sqr - a22 * a12_sqr) + A_Sym(0, 1) * A_Sym(0, 2) * A_Sym(1, 2);
        T sqrt_p = sqrt(p);
        T disc = p * p * p - q * q;
        T phi = ((T)1 / 3) * atan2(sqrt(max((T)0, disc)), q);
        T c = cos(phi), s = sin(phi);
        T sqrt_p_cos = sqrt_p * c;
        T root_three_sqrt_p_sin = sqrt((T)3) * sqrt_p * s;
        lambda(0) = m + 2 * sqrt_p_cos;
        lambda(1) = m - sqrt_p_cos - root_three_sqrt_p_sin;
        lambda(2) = m - sqrt_p_cos + root_three_sqrt_p_sin;
        if (lambda(0) < lambda(1))
            swap(lambda(0), lambda(1));
        if (lambda(1) < lambda(2))
            swap(lambda(1), lambda(2));
        if (lambda(0) < lambda(1))
            swap(lambda(0), lambda(1));
    }
    void fastEigenvectors(const Eigen::Matrix3d& A_Sym,
        const Eigen::Vector3d& lambda,
        Eigen::Matrix3d& V)
    // 71 mults, 44 adds, 3 divs, 3 sqrts
    {
        // flip if necessary so that first eigenvalue is the most different
        using T = double;
        using std::sqrt;
        using std::swap;
        bool flipped = false;
        Eigen::Vector3d lambda_flip(lambda);
        if (lambda(0) - lambda(1) < lambda(1) - lambda(2)) { // 2a
            swap(lambda_flip(0), lambda_flip(2));
            flipped = true;
        }

        // get first eigenvector
        Eigen::Matrix3d C1;
        computeCofactorMtr<3>(A_Sym - lambda_flip(0) * Eigen::Matrix3d::Identity(), C1);
        Eigen::Matrix3d::Index i;
        T norm2 = C1.colwise().squaredNorm().maxCoeff(&i); // 3a + 12m+6a + 9m+6a+1d+1s = 21m+15a+1d+1s
        Eigen::Vector3d v1;
        if (norm2 != 0) {
            T one_over_sqrt = (T)1 / sqrt(norm2);
            v1 = C1.col(i) * one_over_sqrt;
        }
        else
            v1 << 1, 0, 0;

        // form basis for orthogonal complement to v1, and reduce A to this space
        Eigen::Vector3d v1_orthogonal = v1.unitOrthogonal(); // 6m+2a+1d+1s (tweak: 5m+1a+1d+1s)
        Eigen::Matrix<T, 3, 2> other_v;
        other_v.col(0) = v1_orthogonal;
        other_v.col(1) = v1.cross(v1_orthogonal); // 6m+3a (tweak: 4m+1a)
        Eigen::Matrix2d A_reduced = other_v.transpose() * A_Sym * other_v; // 21m+12a (tweak: 18m+9a)

        // find third eigenvector from A_reduced, and fill in second via cross product
        Eigen::Matrix2d C3;
        computeCofactorMtr<2>(A_reduced - lambda_flip(2) * Eigen::Matrix2d::Identity(), C3);
        Eigen::Matrix2d::Index j;
        norm2 = C3.colwise().squaredNorm().maxCoeff(&j); // 3a + 12m+6a + 9m+6a+1d+1s = 21m+15a+1d+1s
        Eigen::Vector3d v3;
        if (norm2 != 0) {
            T one_over_sqrt = (T)1 / sqrt(norm2);
            v3 = other_v * C3.col(j) * one_over_sqrt;
        }
        else
            v3 = other_v.col(0);

        Eigen::Vector3d v2 = v3.cross(v1); // 6m+3a

        // finish
        if (flipped) {
            V.col(0) = v3;
            V.col(1) = v2;
            V.col(2) = -v1;
        }
        else {
            V.col(0) = v1;
            V.col(1) = v2;
            V.col(2) = v3;
        }
    }
    void fastSolveEigenproblem(const Eigen::Matrix3d& A_Sym,
        Eigen::Vector3d& lambda,
        Eigen::Matrix3d& V)
    // 71 mults, 44 adds, 3 divs, 3 sqrts
    {
        fastEigenvalues(A_Sym, lambda);
        fastEigenvectors(A_Sym, lambda, V);
    }

    void fastSVD3d(const Eigen::Matrix3d& A,
        Eigen::Matrix3d& U,
        Eigen::Vector3d& singular_values,
        Eigen::Matrix3d& V)
    // 182 mults, 112 adds, 6 divs, 11 sqrts, 1 atan2, 1 sincos
    {
        using T = double;
        // decompose normal equations
        Eigen::Vector3d lambda;
        fastSolveEigenproblem(A.transpose() * A, lambda, V);

        // compute singular values
        if (lambda(2) < 0)
            lambda = (lambda.array() >= (T)0).select(lambda, (T)0);
        singular_values = lambda.array().sqrt();
        if (A.determinant() < 0)
            singular_values(2) = -singular_values(2);

        // compute singular vectors
        U.col(0) = A * V.col(0);
        T norm = U.col(0).norm();
        if (norm != 0) {
            T one_over_norm = (T)1 / norm;
            U.col(0) = U.col(0) * one_over_norm;
        }
        else
            U.col(0) << 1, 0, 0;
        Eigen::Vector3d v1_orthogonal = U.col(0).unitOrthogonal();
        Eigen::Matrix<T, 3, 2> other_v;
        other_v.col(0) = v1_orthogonal;
        other_v.col(1) = U.col(0).cross(v1_orthogonal);
        Eigen::Vector2d w = other_v.transpose() * A * V.col(1);
        norm = w.norm();
        if (norm != 0) {
            T one_over_norm = (T)1 / norm;
            w = w * one_over_norm;
        }
        else
            w << 1, 0;
        U.col(1) = other_v * w;
        U.col(2) = U.col(0).cross(U.col(1));
    }

    void fastComputeSingularValues3d(const Eigen::Matrix3d& A,
        Eigen::Vector3d& singular_values)
    {
        using T = double;
        // decompose normal equations
        Eigen::Vector3d lambda;
        fastEigenvalues(A.transpose() * A, lambda);

        // compute singular values
        if (lambda(2) < 0)
            lambda = (lambda.array() >= (T)0).select(lambda, (T)0);
        singular_values = lambda.array().sqrt();
        if (A.determinant() < 0)
            singular_values(2) = -singular_values(2);
    }

public:
    const typename Eigen::JacobiSVD<MatrixType>::SingularValuesType& singularValues(void) const
    {
        if (flipped_sigma) {
            return singularValues_flipped;
        }
        else {
            return Eigen::JacobiSVD<MatrixType>::singularValues();
        }
    }
    const MatrixType& matrixU(void) const
    {
        if (flipped_U) {
            return matrixU_flipped;
        }
        else {
            return Eigen::JacobiSVD<MatrixType>::matrixU();
        }
    }
    const MatrixType& matrixV(void) const
    {
        if (flipped_V) {
            return matrixV_flipped;
        }
        else {
            return Eigen::JacobiSVD<MatrixType>::matrixV();
        }
    }

    void setIdentity(void)
    {
        flipped_sigma = true;
        flipped_V = true;
        flipped_U = true;

        matrixU_flipped.setIdentity();
        matrixV_flipped.setIdentity();
        singularValues_flipped.setOnes();
    }
};

} // namespace IPC

#endif /* AutoFlipSVD_hpp */
