//
//  LinSysSolver.hpp
//  IPC
//
//  Created by Minchen Li on 6/30/18.
//

#ifndef LinSysSolver_hpp
#define LinSysSolver_hpp

#include "Types.hpp"

#include <Eigen/Eigen>
#include <Eigen/Sparse>

#include <set>
#include <map>
#include <unordered_map>
#include <iostream>

#ifdef USE_TBB
#include <tbb/tbb.h>
#endif

namespace IPC {

template <typename vectorTypeI, typename vectorTypeS>
class LinSysSolver {
protected:
    int numRows;
    Eigen::VectorXi ia, ja;
    std::vector<std::map<int, int>> IJ2aI;
    Eigen::VectorXd a;

public:
    virtual ~LinSysSolver(void){};

public:
    virtual void set_pattern(const std::vector<std::set<int>>& vNeighbor,
        const std::set<int>& fixedVert)
    {
        numRows = static_cast<int>(vNeighbor.size()) * DIM;
        ia.resize(numRows + 1);
        ia[0] = 1; // 1 + nnz above row i
        ja.resize(0); // colI of each element
        IJ2aI.resize(0); // map from matrix index to ja index
        IJ2aI.resize(numRows);

        std::vector<Eigen::VectorXi> ja_v(vNeighbor.size());
        std::vector<int> rowNNZ(numRows);
#ifdef USE_TBB
        tbb::parallel_for(0, (int)vNeighbor.size(), 1, [&](int vI)
#else
        for (int vI = 0; vI < vNeighbor.size(); ++vI)
#endif
            {
                ja_v[vI].resize((vNeighbor[vI].size() + 1) * DIM);

                ja_v[vI][0] = vI * DIM;
                ja_v[vI][1] = ja_v[vI][0] + 1;
                IJ2aI[ja_v[vI][0]][ja_v[vI][0]] = 0;
                IJ2aI[ja_v[vI][0]][ja_v[vI][1]] = 1;
                if constexpr (DIM == 3) {
                    ja_v[vI][2] = ja_v[vI][0] + 2;
                    IJ2aI[ja_v[vI][0]][ja_v[vI][2]] = 2;
                }

                int nnz = DIM;
                for (const auto& nbVI : vNeighbor[vI]) {
                    if (nbVI > vI) {
                        ja_v[vI][nnz] = nbVI * DIM;
                        ja_v[vI][nnz + 1] = ja_v[vI][nnz] + 1;
                        IJ2aI[ja_v[vI][0]][ja_v[vI][nnz]] = nnz;
                        IJ2aI[ja_v[vI][0]][ja_v[vI][nnz + 1]] = nnz + 1;
                        if constexpr (DIM == 3) {
                            ja_v[vI][nnz + 2] = ja_v[vI][nnz] + 2;
                            IJ2aI[ja_v[vI][0]][ja_v[vI][nnz + 2]] = nnz + 2;
                        }
                        nnz += DIM;
                    }
                }

                rowNNZ[ja_v[vI][0]] = nnz;
                if constexpr (DIM == 2) {
                    ja_v[vI].conservativeResize(nnz * DIM - 1);
                    ja_v[vI].tail(nnz - 1) = ja_v[vI].segment(1, nnz - 1);

                    IJ2aI[ja_v[vI][0] + 1] = IJ2aI[ja_v[vI][0]];
                    IJ2aI[ja_v[vI][0] + 1].erase(ja_v[vI][0]);

                    rowNNZ[ja_v[vI][0] + 1] = nnz - 1;
                }
                else {
                    ja_v[vI].conservativeResize(nnz * DIM - 3);
                    ja_v[vI].segment(nnz, nnz - 1) = ja_v[vI].segment(1, nnz - 1);
                    ja_v[vI].tail(nnz - 2) = ja_v[vI].segment(2, nnz - 2);

                    IJ2aI[ja_v[vI][0] + 1] = IJ2aI[ja_v[vI][0]];
                    IJ2aI[ja_v[vI][0] + 1].erase(ja_v[vI][0]);
                    IJ2aI[ja_v[vI][0] + 2] = IJ2aI[ja_v[vI][0] + 1];
                    IJ2aI[ja_v[vI][0] + 2].erase(ja_v[vI][0] + 1);

                    rowNNZ[ja_v[vI][0] + 1] = nnz - 1;
                    rowNNZ[ja_v[vI][0] + 2] = nnz - 2;
                }
            }
#ifdef USE_TBB
        );
#endif

        for (int rowI = 0; rowI < numRows; ++rowI) {
            ia[rowI + 1] = ia[rowI] + rowNNZ[rowI];
        }

        ja.resize(ia[numRows] - 1);
#ifdef USE_TBB
        tbb::parallel_for(0, (int)vNeighbor.size(), 1, [&](int vI)
#else
        for (int vI = 0; vI < vNeighbor.size(); ++vI)
#endif
            {
                int rowIStart = vI * DIM;

                ja.segment(ia[rowIStart] - 1, ja_v[vI].size()) = ja_v[vI];

                for (auto& indexI : IJ2aI[rowIStart]) {
                    indexI.second += ia[rowIStart] - 1;
                }
                for (auto& indexI : IJ2aI[rowIStart + 1]) {
                    indexI.second += ia[rowIStart + 1] - 2;
                }
                if constexpr (DIM == 3) {
                    for (auto& indexI : IJ2aI[rowIStart + 2]) {
                        indexI.second += ia[rowIStart + 2] - 3;
                    }
                }
            }
#ifdef USE_TBB
        );
#endif
        ja.array() += 1;
        a.resize(ja.size());

        //NOTE: fixed verts nnz entries are not eliminated
    }

    virtual void load(const char* filePath, Eigen::VectorXd& rhs)
    {
        FILE* in = fopen(filePath, "rb");
        assert(in);

        size_t vecSize;
        fread(&vecSize, sizeof(size_t), 1, in);
        std::cout << "ia size " << vecSize << std::endl;
        ia.resize(vecSize);
        fread(ia.data(), sizeof(ia[0]), vecSize, in);

        fread(&vecSize, sizeof(size_t), 1, in);
        std::cout << "ja size " << vecSize << std::endl;
        ja.resize(vecSize);
        fread(ja.data(), sizeof(ja[0]), vecSize, in);

        if (ia[0] == 0) {
            ia.array() += 1;
            ja.array() += 1;
        }

        fread(&vecSize, sizeof(size_t), 1, in);
        std::cout << "a size " << vecSize << std::endl;
        a.resize(vecSize);
        fread(a.data(), sizeof(a[0]), vecSize, in);

        fread(&vecSize, sizeof(size_t), 1, in);
        std::cout << "rhs size " << vecSize << std::endl;
        rhs.resize(vecSize);
        fread(rhs.data(), sizeof(rhs[0]), vecSize, in);

        numRows = vecSize;

        fclose(in);
        std::cout << "load done" << std::endl;
    }
    virtual void write(const char* filePath, const Eigen::VectorXd& rhs)
    {
        FILE* out = fopen(filePath, "wb");

        size_t vecSize = ia.size();
        fwrite(&vecSize, sizeof(vecSize), 1, out);
        fwrite(ia.data(), sizeof(ia[0]), ia.size(), out);

        vecSize = ja.size();
        fwrite(&vecSize, sizeof(vecSize), 1, out);
        fwrite(ja.data(), sizeof(ja[0]), ja.size(), out);

        vecSize = a.size();
        fwrite(&vecSize, sizeof(vecSize), 1, out);
        fwrite(a.data(), sizeof(a[0]), a.size(), out);

        vecSize = rhs.size();
        fwrite(&vecSize, sizeof(vecSize), 1, out);
        fwrite(rhs.data(), sizeof(rhs[0]), rhs.size(), out);

        fclose(out);
    }

    virtual void set_pattern(const Eigen::SparseMatrix<double>& mtr)
    {
        //NOTE: mtr must be SPD

        numRows = static_cast<int>(mtr.rows());

        ja.conservativeResize(mtr.nonZeros());
        memcpy(ja.data(), mtr.innerIndexPtr(),
            mtr.nonZeros() * sizeof(mtr.innerIndexPtr()[0]));

        ia.conservativeResize(numRows + 1);
        memcpy(ia.data(), mtr.outerIndexPtr(),
            (numRows + 1) * sizeof(mtr.outerIndexPtr()[0]));

        a.conservativeResize(mtr.nonZeros());
        memcpy(a.data(), mtr.valuePtr(),
            mtr.nonZeros() * sizeof(mtr.valuePtr()[0]));
    }

    virtual void analyze_pattern(void) = 0;

    virtual bool factorize(void) = 0;

    virtual void solve(Eigen::VectorXd& rhs,
        Eigen::VectorXd& result)
        = 0;

    virtual void multiply(const Eigen::VectorXd& x,
        Eigen::VectorXd& Ax)
    {
        assert(x.size() == numRows);
        assert(IJ2aI.size() == numRows);

        Ax.setZero(numRows);
        for (int rowI = 0; rowI < numRows; ++rowI) {
            for (const auto& colI : IJ2aI[rowI]) {
                Ax[rowI] += a[colI.second] * x[colI.first];
                if (rowI != colI.first) {
                    Ax[colI.first] += a[colI.second] * x[rowI];
                }
            }
        }
    }

public:
    virtual void outputFactorization(const std::string& filePath)
    {
        assert(0 && "please implement!");
    }
    virtual double coeffMtr(int rowI, int colI) const
    {
        if (rowI > colI) {
            // return only upper right part for symmetric matrix
            int temp = rowI;
            rowI = colI;
            colI = temp;
        }
        assert(rowI < IJ2aI.size());
        const auto finder = IJ2aI[rowI].find(colI);
        if (finder != IJ2aI[rowI].end()) {
            return a[finder->second];
        }
        else {
            return 0.0;
        }
    }
    virtual void getCoeffMtr(Eigen::SparseMatrix<double>& mtr) const
    {
        mtr.resize(numRows, numRows);
        mtr.setZero();
        mtr.reserve(a.size() * 2 - numRows);
        for (int rowI = 0; rowI < numRows; rowI++) {
            for (const auto& colIter : IJ2aI[rowI]) {
                mtr.insert(rowI, colIter.first) = a[colIter.second];
                if (rowI != colIter.first) {
                    mtr.insert(colIter.first, rowI) = a[colIter.second];
                }
            }
        }
    }
    virtual void getCoeffMtr_lower(Eigen::SparseMatrix<double>& mtr) const
    {
        assert(numRows > 0);

        mtr.conservativeResize(numRows, numRows);
        mtr.reserve(a.size());

        memcpy(mtr.innerIndexPtr(), ja.data(), ja.size() * sizeof(ja[0]));
        memcpy(mtr.outerIndexPtr(), ia.data(), ia.size() * sizeof(ia[0]));
        memcpy(mtr.valuePtr(), a.data(), a.size() * sizeof(a[0]));
    }
    virtual void getTriplets(const Eigen::VectorXi& nodeList,
        std::vector<Eigen::Triplet<double>>& triplet) const
    {
        std::map<int, int> rowIMapper;
        for (int i = 0; i < nodeList.size(); ++i) {
            int startI = i * DIM;
            int startRowI = nodeList[i] * DIM;

            rowIMapper[startRowI] = startI;
            rowIMapper[startRowI + 1] = startI + 1;
            if constexpr (DIM == 3) {
                rowIMapper[startRowI + 2] = startI + 2;
            }
        }

        triplet.resize(0);
        for (int rowI = 0; rowI < numRows; rowI++) {
            auto rowIFinder = rowIMapper.find(rowI);
            for (const auto& colIter : IJ2aI[rowI]) {
                auto colIFinder = rowIMapper.find(colIter.first);
                if (rowIFinder != rowIMapper.end() && colIFinder != rowIMapper.end()) {
                    triplet.emplace_back(rowIFinder->second, colIFinder->second, a[colIter.second]);
                    if (rowIFinder->second != colIFinder->second) {
                        triplet.emplace_back(colIFinder->second, rowIFinder->second, a[colIter.second]);
                    }
                }
            }
        }
    }
    virtual void setCoeff(int rowI, int colI, double val)
    {
        if (rowI <= colI) {
            assert(rowI < IJ2aI.size());
            const auto finder = IJ2aI[rowI].find(colI);
            assert(finder != IJ2aI[rowI].end());
            a[finder->second] = val;
        }
    }
    virtual void setCoeff(const LinSysSolver<vectorTypeI, vectorTypeS>* other,
        double multiplier)
    {
        assert(numRows == other->numRows);
        assert(ja.size() == other->a.size());

        a = multiplier * other->a;
    }
    virtual void setZero(void)
    {
        a.setZero();
    }
    virtual void setUnit_row(int rowI)
    {
        assert(numRows == IJ2aI.size());
        assert(rowI < numRows);
        for (const auto& colIter : IJ2aI[rowI]) {
            a[colIter.second] = (colIter.first == rowI);
        }
    }
    virtual void setUnit_row(int rowI, std::unordered_map<int, double>& rowVec)
    {
        assert(numRows == IJ2aI.size());
        assert(rowI < numRows);
        rowVec.clear();
        for (const auto& colIter : IJ2aI[rowI]) {
            rowVec[colIter.first] = a[colIter.second];
            a[colIter.second] = (colIter.first == rowI);
        }
    }
    virtual void setUnit_col(int colI, const std::set<int>& rowVIs)
    {
        assert(numRows == IJ2aI.size());
        assert(colI < numRows);
        for (const auto& rowVI : rowVIs) {
            for (int dimI = 0; dimI < DIM; ++dimI) {
                int rowI = rowVI * DIM + dimI;
                assert(rowI < numRows);
                if (rowI <= colI) {
                    const auto finder = IJ2aI[rowI].find(colI);
                    if (finder != IJ2aI[rowI].end()) {
                        a[finder->second] = (rowI == colI);
                    }
                }
            }
        }
    }
    virtual void setUnit_col_dim1(int colI, const std::set<int>& rowVIs)
    {
        assert(numRows == IJ2aI.size());
        assert(colI < numRows);
        for (const auto& rowI : rowVIs) {
            assert(rowI < numRows);
            if (rowI <= colI) {
                const auto finder = IJ2aI[rowI].find(colI);
                if (finder != IJ2aI[rowI].end()) {
                    a[finder->second] = (rowI == colI);
                }
            }
        }
    }

    virtual void addCoeff(int rowI, int colI, double val)
    {
        if (rowI <= colI) {
            assert(rowI < IJ2aI.size());
            const auto finder = IJ2aI[rowI].find(colI);
            assert(finder != IJ2aI[rowI].end());
            a[finder->second] += val;
        }
    }
    virtual void precondition_diag(const Eigen::VectorXd& input, Eigen::VectorXd& output)
    {
        assert(numRows == input.size());
        output.resize(numRows);
        for (int rowI = 0; rowI < numRows; ++rowI) {
            const auto finder = IJ2aI[rowI].find(rowI);
            assert(finder != IJ2aI[rowI].end());
            output[rowI] = input[rowI] / a[finder->second];
        }
    }
    virtual void getMaxDiag(double& maxDiag)
    {
        maxDiag = -std::numeric_limits<double>::infinity();
        for (int rowI = 0; rowI < numRows; ++rowI) {
            const auto finder = IJ2aI[rowI].find(rowI);
            assert(finder != IJ2aI[rowI].end());
            if (maxDiag < a[finder->second]) {
                maxDiag = a[finder->second];
            }
        }
    }
    virtual void addCoeff(const LinSysSolver<vectorTypeI, vectorTypeS>* other,
        double multiplier)
    {
        assert(numRows == other->numRows);
        if (a.size() == other->a.size()) {
            a += multiplier * other->a;
        }
        else {
            for (int rowI = 0; rowI < numRows; ++rowI) {
                for (const auto& colIter : other->IJ2aI[rowI]) {
                    const auto finder = IJ2aI[rowI].find(colIter.first);
                    if (finder != IJ2aI[rowI].end()) {
                        a[finder->second] += multiplier * other->a[colIter.second];
                    }
                }
            }
        }
    }

    virtual int getNumRows(void) const
    {
        return numRows;
    }
    virtual int getNumNonzeros(void) const
    {
        return a.size();
    }
    virtual const std::vector<std::map<int, int>>& getIJ2aI(void) const
    {
        return IJ2aI;
    }
    virtual Eigen::VectorXi& get_ia(void) { return ia; }
    virtual Eigen::VectorXi& get_ja(void) { return ja; }
    virtual Eigen::VectorXd& get_a(void) { return a; }
    virtual const Eigen::VectorXd& get_a(void) const { return a; }
};

} // namespace IPC

#endif /* LinSysSolver_hpp */
