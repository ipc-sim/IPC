//
//  SpatialHash.hpp
//  IPC
//
//  Created by Minchen Li on 6/26/19.
//

#ifndef SpatialHash_hpp
#define SpatialHash_hpp

#include "Timer.hpp"
#include "Mesh.hpp"

#include <unordered_map>
#include <unordered_set>

extern Timer timer_mt;

namespace IPC {

template <int dim>
class SpatialHash {
public: // data
    Eigen::Matrix<double, 1, dim> leftBottomCorner, rightTopCorner;
    double one_div_voxelSize;
    Eigen::Array<int, 1, dim> voxelCount;
    int voxelCount0x1;

    int surfEdgeStartInd, surfTriStartInd;

    std::unordered_map<int, std::vector<int>> voxel;
    std::vector<std::vector<int>> pointAndEdgeOccupancy;

public: // constructor
    SpatialHash(void) {}
    SpatialHash(const Mesh<dim>& mesh, double voxelSize, bool use_V_prev = false)
    {
        build(mesh, voxelSize, use_V_prev);
    }
    SpatialHash(const Mesh<dim>& mesh, const Eigen::VectorXd& searchDir, double curMaxStepSize, double voxelSize, bool use_V_prev = false)
    {
        build(mesh, searchDir, curMaxStepSize, voxelSize, use_V_prev);
    }

public: // API
    void build(const Mesh<dim>& mesh, double voxelSize, bool use_V_prev = false)
    {
        const Eigen::MatrixXd& V = use_V_prev ? mesh.V_prev : mesh.V;
        leftBottomCorner = V.colwise().minCoeff();
        rightTopCorner = V.colwise().maxCoeff();
        one_div_voxelSize = 1.0 / voxelSize;
        Eigen::Array<double, 1, dim> range = rightTopCorner - leftBottomCorner;
        voxelCount = (range * one_div_voxelSize).ceil().template cast<int>();
        if (voxelCount.minCoeff() <= 0) {
            // cast overflow due to huge search direction
            one_div_voxelSize = 1.0 / (range.maxCoeff() * 1.01);
            voxelCount.setOnes();
        }
        voxelCount0x1 = voxelCount[0] * voxelCount[1];

        surfEdgeStartInd = mesh.SVI.size();
        surfTriStartInd = surfEdgeStartInd + mesh.SFEdges.size();

        timer_mt.start(15);
        // precompute svVAI
        std::vector<Eigen::Array<int, 1, dim>> svVoxelAxisIndex(mesh.SVI.size());
        std::vector<int> vI2SVI(V.rows());
#ifdef USE_TBB
        tbb::parallel_for(0, (int)mesh.SVI.size(), 1, [&](int svI)
#else
        for (int svI = 0; svI < mesh.SVI.size(); ++svI)
#endif
            {
                int vI = mesh.SVI[svI];
                locateVoxelAxisIndex(V.row(vI), svVoxelAxisIndex[svI]);
                vI2SVI[vI] = svI;
            }
#ifdef USE_TBB
        );
#endif

        voxel.clear();

#ifdef PARALLEL_SH_CONSTRUCT
        std::vector<std::pair<int, int>> voxel_tmp;

        for (int svI = 0; svI < mesh.SVI.size(); ++svI) {
            voxel_tmp.emplace_back(locateVoxelIndex(V.row(mesh.SVI[svI])), svI);
        }
#else
        for (int svI = 0; svI < mesh.SVI.size(); ++svI) {
            voxel[locateVoxelIndex(V.row(mesh.SVI[svI]))].emplace_back(svI);
        }
#endif

        timer_mt.start(16);
        std::vector<std::vector<int>> voxelLoc_e(mesh.SFEdges.size());
#ifdef USE_TBB
        tbb::parallel_for(0, (int)mesh.SFEdges.size(), 1, [&](int seCount)
#else
        for (int seCount = 0; seCount < mesh.SFEdges.size(); ++seCount)
#endif
            {
                const auto& seI = mesh.SFEdges[seCount];

                const Eigen::Array<int, 1, dim>& voxelAxisIndex_first = svVoxelAxisIndex[vI2SVI[seI.first]];
                const Eigen::Array<int, 1, dim>& voxelAxisIndex_second = svVoxelAxisIndex[vI2SVI[seI.second]];
                Eigen::Array<int, 1, dim> mins = voxelAxisIndex_first.min(voxelAxisIndex_second);
                Eigen::Array<int, 1, dim> maxs = voxelAxisIndex_first.max(voxelAxisIndex_second);
                for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                    int zOffset = iz * voxelCount0x1;
                    for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                        int yzOffset = iy * voxelCount[0] + zOffset;
                        for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                            voxelLoc_e[seCount].emplace_back(ix + yzOffset);
                        }
                    }
                }
            }
#ifdef USE_TBB
        );
#endif

        timer_mt.start(17);
        std::vector<std::vector<int>> voxelLoc_sf(mesh.SF.rows());
#ifdef USE_TBB
        tbb::parallel_for(0, (int)mesh.SF.rows(), 1, [&](int sfI)
#else
        for (int sfI = 0; sfI < mesh.SF.rows(); ++sfI)
#endif
            {
                const Eigen::Array<int, 1, dim>& voxelAxisIndex0 = svVoxelAxisIndex[vI2SVI[mesh.SF(sfI, 0)]];
                const Eigen::Array<int, 1, dim>& voxelAxisIndex1 = svVoxelAxisIndex[vI2SVI[mesh.SF(sfI, 1)]];
                const Eigen::Array<int, 1, dim>& voxelAxisIndex2 = svVoxelAxisIndex[vI2SVI[mesh.SF(sfI, 2)]];
                Eigen::Array<int, 1, dim> mins = voxelAxisIndex0.min(voxelAxisIndex1).min(voxelAxisIndex2);
                Eigen::Array<int, 1, dim> maxs = voxelAxisIndex0.max(voxelAxisIndex1).max(voxelAxisIndex2);
                for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                    int zOffset = iz * voxelCount0x1;
                    for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                        int yzOffset = iy * voxelCount[0] + zOffset;
                        for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                            voxelLoc_sf[sfI].emplace_back(ix + yzOffset);
                        }
                    }
                }
            }
#ifdef USE_TBB
        );
#endif

        timer_mt.start(18);

#ifdef PARALLEL_SH_CONSTRUCT
        for (int seCount = 0; seCount < voxelLoc_e.size(); ++seCount) {
            for (const auto& voxelI : voxelLoc_e[seCount]) {
                voxel_tmp.emplace_back(voxelI, seCount + surfEdgeStartInd);
            }
        }

        for (int sfI = 0; sfI < voxelLoc_sf.size(); ++sfI) {
            for (const auto& voxelI : voxelLoc_sf[sfI]) {
                voxel_tmp.emplace_back(voxelI, sfI + surfTriStartInd);
            }
        }

#ifdef USE_TBB
        tbb::parallel_sort(voxel_tmp.begin(), voxel_tmp.end(), [](const std::pair<int, int>& f, const std::pair<int, int>& s) { return f.first < s.first; });
#else
        std::sort(voxel_tmp.begin(), voxel_tmp.end(), [](const std::pair<int, int>& f, const std::pair<int, int>& s) { return f.first < s.first; });
#endif
        std::vector<std::pair<int, std::vector<int>>> voxel_tmp_merged;
        voxel_tmp_merged.reserve(voxel_tmp.size());
        int current_voxel = -1;
        for (const auto& v : voxel_tmp) {
            if (current_voxel != v.first) {
                assert(current_voxel < v.first);
                voxel_tmp_merged.emplace_back();
                voxel_tmp_merged.back().first = v.first;
                current_voxel = v.first;
            }

            voxel_tmp_merged.back().second.push_back(v.second);
        }
        assert(voxel_tmp_merged.size() <= voxel_tmp.size());

        voxel.insert(voxel_tmp_merged.begin(), voxel_tmp_merged.end());
#else
        for (int seCount = 0; seCount < voxelLoc_e.size(); ++seCount) {
            for (const auto& voxelI : voxelLoc_e[seCount]) {
                voxel[voxelI].emplace_back(seCount + surfEdgeStartInd);
            }
        }
        for (int sfI = 0; sfI < voxelLoc_sf.size(); ++sfI) {
            for (const auto& voxelI : voxelLoc_sf[sfI]) {
                voxel[voxelI].emplace_back(sfI + surfTriStartInd);
            }
        }
#endif

        timer_mt.stop();
    }

    void queryPointForTriangles(const Eigen::Matrix<double, 1, dim>& pos,
        double radius, std::unordered_set<int>& triInds) const
    {
        Eigen::Array<int, 1, dim> mins, maxs;
        locateVoxelAxisIndex(pos.array() - radius, mins);
        locateVoxelAxisIndex(pos.array() + radius, maxs);
        mins = mins.max(Eigen::Array<int, 1, dim>::Zero());
        maxs = maxs.min(voxelCount - 1);

        triInds.clear();
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    const auto voxelI = voxel.find(ix + yzOffset);
                    if (voxelI != voxel.end()) {
                        for (const auto& indI : voxelI->second) {
                            if (indI >= surfTriStartInd) {
                                triInds.insert(indI - surfTriStartInd);
                            }
                        }
                    }
                }
            }
        }
    }
    void queryPointForTriangles(const Eigen::Matrix<double, 1, dim>& pos,
        const Eigen::Matrix<double, 1, dim>& dir,
        std::unordered_set<int>& triInds, double radius = 0) const
    {
        Eigen::Array<int, 1, dim> mins, maxs;
        locateVoxelAxisIndex(pos.array().min((pos + dir).array()) - radius, mins);
        locateVoxelAxisIndex(pos.array().max((pos + dir).array()) + radius, maxs);
        mins = mins.max(Eigen::Array<int, 1, dim>::Zero());
        maxs = maxs.min(voxelCount - 1);

        triInds.clear();
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    const auto voxelI = voxel.find(ix + yzOffset);
                    if (voxelI != voxel.end()) {
                        for (const auto& indI : voxelI->second) {
                            if (indI >= surfTriStartInd) {
                                triInds.insert(indI - surfTriStartInd);
                            }
                        }
                    }
                }
            }
        }
    }
    void queryPointForPrimitives(const Eigen::Matrix<double, 1, dim>& pos,
        const Eigen::Matrix<double, 1, dim>& dir, std::unordered_set<int>& sVInds,
        std::unordered_set<int>& sEdgeInds, std::unordered_set<int>& sTriInds) const
    {
        Eigen::Array<int, 1, dim> mins, maxs;
        locateVoxelAxisIndex(pos.array().min((pos + dir).array()), mins);
        locateVoxelAxisIndex(pos.array().max((pos + dir).array()), maxs);
        mins = mins.max(Eigen::Array<int, 1, dim>::Zero());
        maxs = maxs.min(voxelCount - 1);

        sVInds.clear();
        sEdgeInds.clear();
        sTriInds.clear();
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    const auto voxelI = voxel.find(ix + yzOffset);
                    if (voxelI != voxel.end()) {
                        for (const auto& indI : voxelI->second) {
                            if (indI < surfEdgeStartInd) {
                                sVInds.insert(indI);
                            }
                            else if (indI < surfTriStartInd) {
                                sEdgeInds.insert(indI - surfEdgeStartInd);
                            }
                            else {
                                sTriInds.insert(indI - surfTriStartInd);
                            }
                        }
                    }
                }
            }
        }
    }

    void queryEdgeForPE(const Eigen::Matrix<double, 1, dim>& vBegin,
        const Eigen::Matrix<double, 1, dim>& vEnd,
        std::vector<int>& svInds, std::vector<int>& edgeInds) const
    {
        // timer_mt.start(19);
        Eigen::Matrix<double, 1, dim> leftBottom = vBegin.array().min(vEnd.array());
        Eigen::Matrix<double, 1, dim> rightTop = vBegin.array().max(vEnd.array());
        Eigen::Array<int, 1, dim> mins, maxs;
        locateVoxelAxisIndex(leftBottom, mins);
        locateVoxelAxisIndex(rightTop, maxs);
        mins = mins.max(Eigen::Array<int, 1, dim>::Zero());
        maxs = maxs.min(voxelCount - 1);

        // timer_mt.start(20);
        svInds.resize(0);
        edgeInds.resize(0);
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    // timer_mt.start(21);
                    const auto voxelI = voxel.find(ix + yzOffset);
                    if (voxelI != voxel.end()) {
                        for (const auto& indI : voxelI->second) {
                            if (indI < surfEdgeStartInd) {
                                svInds.emplace_back(indI);
                            }
                            else if (indI < surfTriStartInd) {
                                edgeInds.emplace_back(indI - surfEdgeStartInd);
                            }
                        }
                    }
                    // timer_mt.start(20);
                }
            }
        }
        std::sort(edgeInds.begin(), edgeInds.end());
        edgeInds.erase(std::unique(edgeInds.begin(), edgeInds.end()), edgeInds.end());
        std::sort(svInds.begin(), svInds.end());
        svInds.erase(std::unique(svInds.begin(), svInds.end()), svInds.end());
        // timer_mt.stop();
    }
    void queryEdgeForEdges(const Eigen::Matrix<double, 1, dim>& vBegin,
        const Eigen::Matrix<double, 1, dim>& vEnd,
        double radius, std::vector<int>& edgeInds, int eIq = -1) const
    {
        // timer_mt.start(19);
        Eigen::Matrix<double, 1, dim> leftBottom = vBegin.array().min(vEnd.array()) - radius;
        Eigen::Matrix<double, 1, dim> rightTop = vBegin.array().max(vEnd.array()) + radius;
        Eigen::Array<int, 1, dim> mins, maxs;
        locateVoxelAxisIndex(leftBottom, mins);
        locateVoxelAxisIndex(rightTop, maxs);
        mins = mins.max(Eigen::Array<int, 1, dim>::Zero());
        maxs = maxs.min(voxelCount - 1);

        // timer_mt.start(20);
        edgeInds.resize(0);
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    // timer_mt.start(21);
                    const auto voxelI = voxel.find(ix + yzOffset);
                    if (voxelI != voxel.end()) {
                        for (const auto& indI : voxelI->second) {
                            if (indI >= surfEdgeStartInd && indI < surfTriStartInd && indI - surfEdgeStartInd > eIq) {
                                edgeInds.emplace_back(indI - surfEdgeStartInd);
                            }
                        }
                    }
                    // timer_mt.start(20);
                }
            }
        }
        std::sort(edgeInds.begin(), edgeInds.end());
        edgeInds.erase(std::unique(edgeInds.begin(), edgeInds.end()), edgeInds.end());
        // timer_mt.stop();
    }
    void queryEdgeForEdgesWithBBoxCheck(
        const Mesh<dim>& mesh,
        const Eigen::Matrix<double, 1, dim>& vBegin,
        const Eigen::Matrix<double, 1, dim>& vEnd,
        double radius, std::vector<int>& edgeInds,
        int eIq = -1) const
    {
        // timer_mt.start(19);
        Eigen::Matrix<double, 1, dim> leftBottom = vBegin.array().min(vEnd.array()) - radius;
        Eigen::Matrix<double, 1, dim> rightTop = vBegin.array().max(vEnd.array()) + radius;
        Eigen::Array<int, 1, dim> mins, maxs;
        locateVoxelAxisIndex(leftBottom, mins);
        locateVoxelAxisIndex(rightTop, maxs);
        mins = mins.max(Eigen::Array<int, 1, dim>::Zero());
        maxs = maxs.min(voxelCount - 1);

        // timer_mt.start(20);
        edgeInds.resize(0);
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    // timer_mt.start(21);
                    const auto voxelI = voxel.find(ix + yzOffset);
                    if (voxelI != voxel.end()) {
                        for (const auto& indI : voxelI->second) {
                            if (indI >= surfEdgeStartInd && indI < surfTriStartInd && indI - surfEdgeStartInd > eIq) {
                                int seJ = indI - surfEdgeStartInd;
                                const Eigen::Matrix<double, 1, dim>& eJ_v0 = mesh.V.row(mesh.SFEdges[seJ].first);
                                const Eigen::Matrix<double, 1, dim>& eJ_v1 = mesh.V.row(mesh.SFEdges[seJ].second);
                                Eigen::Array<double, 1, dim> bboxEJTopRight = eJ_v0.array().max(eJ_v1.array());
                                Eigen::Array<double, 1, dim> bboxEJBottomLeft = eJ_v0.array().min(eJ_v1.array());
                                if (!((bboxEJBottomLeft - rightTop.array() > 0.0).any() || (leftBottom.array() - bboxEJTopRight > 0.0).any())) {
                                    edgeInds.emplace_back(indI - surfEdgeStartInd);
                                }
                            }
                        }
                    }
                    // timer_mt.start(20);
                }
            }
        }
        std::sort(edgeInds.begin(), edgeInds.end());
        edgeInds.erase(std::unique(edgeInds.begin(), edgeInds.end()), edgeInds.end());
        // timer_mt.stop();
    }
    void queryEdgeForEdges(const Eigen::Matrix<double, 1, dim>& vBegin,
        const Eigen::Matrix<double, 1, dim>& vEnd,
        const Eigen::Matrix<double, 1, dim>& pBegin,
        const Eigen::Matrix<double, 1, dim>& pEnd,
        std::vector<int>& edgeInds, double radius = 0, int eIq = -1) const
    {
        Eigen::Matrix<double, 1, dim> vtBegin = vBegin + pBegin;
        Eigen::Matrix<double, 1, dim> vtEnd = vEnd + pEnd;
        Eigen::Matrix<double, 1, dim> leftBottom = vBegin.array().min(vEnd.array()).min(vtBegin.array()).min(vtEnd.array());
        Eigen::Matrix<double, 1, dim> rightTop = vBegin.array().max(vEnd.array()).max(vtBegin.array()).max(vtEnd.array());
        Eigen::Array<int, 1, dim> mins, maxs;
        locateVoxelAxisIndex(leftBottom.array() - radius, mins);
        locateVoxelAxisIndex(rightTop.array() + radius, maxs);
        mins = mins.max(Eigen::Array<int, 1, dim>::Zero());
        maxs = maxs.min(voxelCount - 1);

        edgeInds.resize(0);
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    const auto voxelI = voxel.find(ix + yzOffset);
                    if (voxelI != voxel.end()) {
                        for (const auto& indI : voxelI->second) {
                            if (indI >= surfEdgeStartInd && indI < surfTriStartInd && indI - surfEdgeStartInd > eIq) {
                                edgeInds.emplace_back(indI - surfEdgeStartInd);
                            }
                        }
                    }
                }
            }
        }
        std::sort(edgeInds.begin(), edgeInds.end());
        edgeInds.erase(std::unique(edgeInds.begin(), edgeInds.end()), edgeInds.end());
    }

    void queryTriangleForPoints(const Eigen::Matrix<double, 1, dim>& v0,
        const Eigen::Matrix<double, 1, dim>& v1,
        const Eigen::Matrix<double, 1, dim>& v2,
        double radius, std::unordered_set<int>& pointInds) const
    {
        Eigen::Matrix<double, 1, dim> leftBottom = v0.array().min(v1.array()).min(v2.array());
        Eigen::Matrix<double, 1, dim> rightTop = v0.array().max(v1.array()).max(v2.array());
        Eigen::Array<int, 1, dim> mins, maxs;
        locateVoxelAxisIndex(leftBottom.array() - radius, mins);
        locateVoxelAxisIndex(rightTop.array() + radius, maxs);
        mins = mins.max(Eigen::Array<int, 1, dim>::Zero());
        maxs = maxs.min(voxelCount - 1);

        pointInds.clear();
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    const auto voxelI = voxel.find(ix + yzOffset);
                    if (voxelI != voxel.end()) {
                        for (const auto& indI : voxelI->second) {
                            if (indI < surfEdgeStartInd) {
                                pointInds.insert(indI);
                            }
                        }
                    }
                }
            }
        }
    }
    void queryTriangleForPoints(const Eigen::Matrix<double, 1, dim>& v0,
        const Eigen::Matrix<double, 1, dim>& v1,
        const Eigen::Matrix<double, 1, dim>& v2,
        const Eigen::Matrix<double, 1, dim>& p0,
        const Eigen::Matrix<double, 1, dim>& p1,
        const Eigen::Matrix<double, 1, dim>& p2,
        std::unordered_set<int>& pointInds) const
    {
        Eigen::Matrix<double, 1, dim> v0t = v0 + p0, v1t = v1 + p1, v2t = v2 + p2;
        Eigen::Matrix<double, 1, dim> leftBottom = v0.array().min(v1.array()).min(v2.array()).min(v0t.array()).min(v1t.array()).min(v2t.array());
        Eigen::Matrix<double, 1, dim> rightTop = v0.array().max(v1.array()).max(v2.array()).max(v0t.array()).max(v1t.array()).max(v2t.array());
        Eigen::Array<int, 1, dim> mins, maxs;
        locateVoxelAxisIndex(leftBottom.array(), mins);
        locateVoxelAxisIndex(rightTop.array(), maxs);
        mins = mins.max(Eigen::Array<int, 1, dim>::Zero());
        maxs = maxs.min(voxelCount - 1);

        pointInds.clear();
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    const auto voxelI = voxel.find(ix + yzOffset);
                    if (voxelI != voxel.end()) {
                        for (const auto& indI : voxelI->second) {
                            if (indI < surfEdgeStartInd) {
                                pointInds.insert(indI);
                            }
                        }
                    }
                }
            }
        }
    }

    void queryTriangleForEdges(const Eigen::Matrix<double, 1, dim>& v0,
        const Eigen::Matrix<double, 1, dim>& v1,
        const Eigen::Matrix<double, 1, dim>& v2,
        double radius, std::unordered_set<int>& edgeInds) const
    {
        Eigen::Matrix<double, 1, dim> leftBottom = v0.array().min(v1.array()).min(v2.array());
        Eigen::Matrix<double, 1, dim> rightTop = v0.array().max(v1.array()).max(v2.array());
        Eigen::Array<int, 1, dim> mins, maxs;
        locateVoxelAxisIndex(leftBottom.array() - radius, mins);
        locateVoxelAxisIndex(rightTop.array() + radius, maxs);
        mins = mins.max(Eigen::Array<int, 1, dim>::Zero());
        maxs = maxs.min(voxelCount - 1);

        edgeInds.clear();
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    const auto voxelI = voxel.find(ix + yzOffset);
                    if (voxelI != voxel.end()) {
                        for (const auto& indI : voxelI->second) {
                            if (indI >= surfEdgeStartInd && indI < surfTriStartInd) {
                                edgeInds.insert(indI - surfEdgeStartInd);
                            }
                        }
                    }
                }
            }
        }
    }

    void queryEdgeForTriangles(const Eigen::Matrix<double, 1, dim>& vBegin,
        const Eigen::Matrix<double, 1, dim>& vEnd,
        double radius, std::unordered_set<int>& triInds) const
    {
        Eigen::Matrix<double, 1, dim> leftBottom = vBegin.array().min(vEnd.array());
        Eigen::Matrix<double, 1, dim> rightTop = vBegin.array().max(vEnd.array());
        Eigen::Array<int, 1, dim> mins, maxs;
        locateVoxelAxisIndex(leftBottom.array() - radius, mins);
        locateVoxelAxisIndex(rightTop.array() + radius, maxs);
        mins = mins.max(Eigen::Array<int, 1, dim>::Zero());
        maxs = maxs.min(voxelCount - 1);

        triInds.clear();
        for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
            int zOffset = iz * voxelCount0x1;
            for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                int yzOffset = iy * voxelCount[0] + zOffset;
                for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                    const auto voxelI = voxel.find(ix + yzOffset);
                    if (voxelI != voxel.end()) {
                        for (const auto& indI : voxelI->second) {
                            if (indI >= surfTriStartInd) {
                                triInds.insert(indI - surfTriStartInd);
                            }
                        }
                    }
                }
            }
        }
    }

    void build(const Mesh<dim>& mesh, const Eigen::VectorXd& searchDir, double& curMaxStepSize, double voxelSize, bool use_V_prev = false)
    {
        // Can't do it as we sometimes still need to use the Legacy meshCO
        // Eigen::VectorXd searchDir = p_searchDir;
        // Eigen::Matrix<double, dim, 1> avgSD;
        // avgSD.setZero();
        // for (int svI = 0; svI < mesh.SVI.size(); ++svI) {
        //     avgSD += searchDir.template segment<dim>(mesh.SVI[svI] * dim);
        // }
        // avgSD /= mesh.SVI.size();
        // for (int svI = 0; svI < mesh.SVI.size(); ++svI) {
        //     searchDir.template segment<dim>(mesh.SVI[svI] * dim) -= avgSD;
        // }

        double pSize = 0;
        for (int svI = 0; svI < mesh.SVI.size(); ++svI) {
            int vI = mesh.SVI[svI];
            pSize += std::abs(searchDir[vI * dim]);
            pSize += std::abs(searchDir[vI * dim + 1]);
            if constexpr (dim == 3) {
                pSize += std::abs(searchDir[vI * dim + 2]);
            }
        }
        pSize /= mesh.SVI.size() * dim;

        const double spanSize = curMaxStepSize * pSize / voxelSize;
        if (spanSize > 1) {
            curMaxStepSize /= spanSize;
            // curMaxStepSize reduced for CCD spatial hash efficiency
        }

        const Eigen::MatrixXd& V = use_V_prev ? mesh.V_prev : mesh.V;
        Eigen::MatrixXd SVt(mesh.SVI.size(), dim);
        std::unordered_map<int, int> vI2SVI;
        for (int svI = 0; svI < mesh.SVI.size(); ++svI) {
            int vI = mesh.SVI[svI];
            SVt.row(svI) = V.row(vI) + curMaxStepSize * searchDir.template segment<dim>(vI * dim).transpose();
            vI2SVI[vI] = svI;
        }

        leftBottomCorner = V.colwise().minCoeff().array().min(SVt.colwise().minCoeff().array());
        rightTopCorner = V.colwise().maxCoeff().array().max(SVt.colwise().maxCoeff().array());
        one_div_voxelSize = 1.0 / voxelSize;
        Eigen::Array<double, 1, dim> range = rightTopCorner - leftBottomCorner;
        voxelCount = (range * one_div_voxelSize).ceil().template cast<int>();
        if (voxelCount.minCoeff() <= 0) {
            // cast overflow due to huge search direction
            one_div_voxelSize = 1.0 / (range.maxCoeff() * 1.01);
            voxelCount.setOnes();
        }
        voxelCount0x1 = voxelCount[0] * voxelCount[1];

        surfEdgeStartInd = mesh.SVI.size();
        surfTriStartInd = surfEdgeStartInd + mesh.SFEdges.size();

        // precompute svVAI
        std::vector<Eigen::Array<int, 1, dim>> svMinVAI(mesh.SVI.size());
        std::vector<Eigen::Array<int, 1, dim>> svMaxVAI(mesh.SVI.size());
#ifdef USE_TBB
        tbb::parallel_for(0, (int)mesh.SVI.size(), 1, [&](int svI)
#else
        for (int svI = 0; svI < mesh.SVI.size(); ++svI)
#endif
            {
                int vI = mesh.SVI[svI];
                Eigen::Array<int, 1, dim> v0VAI, vtVAI;
                locateVoxelAxisIndex(V.row(vI), v0VAI);
                locateVoxelAxisIndex(SVt.row(svI), vtVAI);
                svMinVAI[svI] = v0VAI.min(vtVAI);
                svMaxVAI[svI] = v0VAI.max(vtVAI);
            }
#ifdef USE_TBB
        );
#endif

        voxel.clear(); //TODO: parallel insert
        pointAndEdgeOccupancy.resize(0);
        pointAndEdgeOccupancy.resize(surfTriStartInd);

#ifdef USE_TBB
        tbb::parallel_for(0, (int)mesh.SVI.size(), 1, [&](int svI)
#else
        for (int svI = 0; svI < mesh.SVI.size(); ++svI)
#endif
            {
                const Eigen::Array<int, 1, dim>& mins = svMinVAI[svI];
                const Eigen::Array<int, 1, dim>& maxs = svMaxVAI[svI];
                pointAndEdgeOccupancy[svI].reserve((maxs - mins + 1).prod());
                for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                    int zOffset = iz * voxelCount0x1;
                    for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                        int yzOffset = iy * voxelCount[0] + zOffset;
                        for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                            pointAndEdgeOccupancy[svI].emplace_back(ix + yzOffset);
                        }
                    }
                }
            }
#ifdef USE_TBB
        );
#endif

#ifdef USE_TBB
        tbb::parallel_for(0, (int)mesh.SFEdges.size(), 1, [&](int seCount)
#else
        for (int seCount = 0; seCount < mesh.SFEdges.size(); ++seCount)
#endif
            {
                int seIInd = seCount + surfEdgeStartInd;
                const auto& seI = mesh.SFEdges[seCount];

                Eigen::Array<int, 1, dim> mins = svMinVAI[vI2SVI[seI.first]].min(svMinVAI[vI2SVI[seI.second]]);
                Eigen::Array<int, 1, dim> maxs = svMaxVAI[vI2SVI[seI.first]].max(svMaxVAI[vI2SVI[seI.second]]);
                pointAndEdgeOccupancy[seIInd].reserve((maxs - mins + 1).prod());
                for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                    int zOffset = iz * voxelCount0x1;
                    for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                        int yzOffset = iy * voxelCount[0] + zOffset;
                        for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                            pointAndEdgeOccupancy[seIInd].emplace_back(ix + yzOffset);
                        }
                    }
                }
            }
#ifdef USE_TBB
        );
#endif

        std::vector<std::vector<int>> voxelLoc_sf(mesh.SF.rows());
#ifdef USE_TBB
        tbb::parallel_for(0, (int)mesh.SF.rows(), 1, [&](int sfI)
#else
        for (int sfI = 0; sfI < mesh.SF.rows(); ++sfI)
#endif
            {
                Eigen::Array<int, 1, dim> mins = svMinVAI[vI2SVI[mesh.SF(sfI, 0)]].min(svMinVAI[vI2SVI[mesh.SF(sfI, 1)]]).min(svMinVAI[vI2SVI[mesh.SF(sfI, 2)]]);
                Eigen::Array<int, 1, dim> maxs = svMaxVAI[vI2SVI[mesh.SF(sfI, 0)]].max(svMaxVAI[vI2SVI[mesh.SF(sfI, 1)]]).max(svMaxVAI[vI2SVI[mesh.SF(sfI, 2)]]);
                for (int iz = mins[2]; iz <= maxs[2]; ++iz) {
                    int zOffset = iz * voxelCount0x1;
                    for (int iy = mins[1]; iy <= maxs[1]; ++iy) {
                        int yzOffset = iy * voxelCount[0] + zOffset;
                        for (int ix = mins[0]; ix <= maxs[0]; ++ix) {
                            voxelLoc_sf[sfI].emplace_back(ix + yzOffset);
                        }
                    }
                }
            }
#ifdef USE_TBB
        );
#endif

        for (int i = 0; i < pointAndEdgeOccupancy.size(); ++i) {
            for (const auto& voxelI : pointAndEdgeOccupancy[i]) {
                voxel[voxelI].emplace_back(i);
            }
        }
        for (int sfI = 0; sfI < voxelLoc_sf.size(); ++sfI) {
            for (const auto& voxelI : voxelLoc_sf[sfI]) {
                voxel[voxelI].emplace_back(sfI + surfTriStartInd);
            }
        }
    }

    void queryPointForPrimitives(int svI, std::unordered_set<int>& sVInds,
        std::unordered_set<int>& sEdgeInds, std::unordered_set<int>& sTriInds) const
    {
        sVInds.clear();
        sEdgeInds.clear();
        sTriInds.clear();
        for (const auto& voxelInd : pointAndEdgeOccupancy[svI]) {
            const auto& voxelI = voxel.find(voxelInd);
            assert(voxelI != voxel.end());
            for (const auto& indI : voxelI->second) {
                if (indI < surfEdgeStartInd) {
                    sVInds.insert(indI);
                }
                else if (indI < surfTriStartInd) {
                    sEdgeInds.insert(indI - surfEdgeStartInd);
                }
                else {
                    sTriInds.insert(indI - surfTriStartInd);
                }
            }
        }
    }
    void queryPointForTriangles(int svI, std::unordered_set<int>& sTriInds) const
    {
        sTriInds.clear();
        for (const auto& voxelInd : pointAndEdgeOccupancy[svI]) {
            const auto& voxelI = voxel.find(voxelInd);
            assert(voxelI != voxel.end());
            for (const auto& indI : voxelI->second) {
                if (indI >= surfTriStartInd) {
                    sTriInds.insert(indI - surfTriStartInd);
                }
            }
        }
    }

    // will only put edges with larger than seI index into sEdgeInds
    void queryEdgeForEdges(int seI, std::unordered_set<int>& sEdgeInds) const
    {
        sEdgeInds.clear();
        for (const auto& voxelInd : pointAndEdgeOccupancy[seI + surfEdgeStartInd]) {
            const auto& voxelI = voxel.find(voxelInd);
            assert(voxelI != voxel.end());
            for (const auto& indI : voxelI->second) {
                if (indI >= surfEdgeStartInd && indI < surfTriStartInd && indI - surfEdgeStartInd > seI) {
                    sEdgeInds.insert(indI - surfEdgeStartInd);
                }
            }
        }
    }

    void queryEdgeForEdgesWithBBoxCheck(const Mesh<dim>& mesh,
        const Eigen::VectorXd& searchDir, double curMaxStepSize,
        int seI, std::unordered_set<int>& sEdgeInds) const
    {
        const Eigen::Matrix<double, 1, dim>& eI_v0 = mesh.V.row(mesh.SFEdges[seI].first);
        const Eigen::Matrix<double, 1, dim>& eI_v1 = mesh.V.row(mesh.SFEdges[seI].second);
        Eigen::Matrix<double, 1, dim> eI_v0t = eI_v0 + curMaxStepSize * searchDir.template segment<dim>(mesh.SFEdges[seI].first * dim).transpose();
        Eigen::Matrix<double, 1, dim> eI_v1t = eI_v1 + curMaxStepSize * searchDir.template segment<dim>(mesh.SFEdges[seI].second * dim).transpose();
        Eigen::Array<double, 1, dim> bboxEITopRight = eI_v0.array().max(eI_v0t.array()).max(eI_v1.array()).max(eI_v1t.array());
        Eigen::Array<double, 1, dim> bboxEIBottomLeft = eI_v0.array().min(eI_v0t.array()).min(eI_v1.array()).min(eI_v1t.array());
        sEdgeInds.clear();
        for (const auto& voxelInd : pointAndEdgeOccupancy[seI + surfEdgeStartInd]) {
            const auto& voxelI = voxel.find(voxelInd);
            assert(voxelI != voxel.end());
            for (const auto& indI : voxelI->second) {
                if (indI >= surfEdgeStartInd && indI < surfTriStartInd && indI - surfEdgeStartInd > seI) {
                    int seJ = indI - surfEdgeStartInd;
                    const Eigen::Matrix<double, 1, dim>& eJ_v0 = mesh.V.row(mesh.SFEdges[seJ].first);
                    const Eigen::Matrix<double, 1, dim>& eJ_v1 = mesh.V.row(mesh.SFEdges[seJ].second);
                    Eigen::Matrix<double, 1, dim> eJ_v0t = eJ_v0 + curMaxStepSize * searchDir.template segment<dim>(mesh.SFEdges[seJ].first * dim).transpose();
                    Eigen::Matrix<double, 1, dim> eJ_v1t = eJ_v1 + curMaxStepSize * searchDir.template segment<dim>(mesh.SFEdges[seJ].second * dim).transpose();
                    Eigen::Array<double, 1, dim> bboxEJTopRight = eJ_v0.array().max(eJ_v0t.array()).max(eJ_v1.array()).max(eJ_v1t.array());
                    Eigen::Array<double, 1, dim> bboxEJBottomLeft = eJ_v0.array().min(eJ_v0t.array()).min(eJ_v1.array()).min(eJ_v1t.array());
                    if (!((bboxEJBottomLeft - bboxEITopRight > 0.0).any() || (bboxEIBottomLeft - bboxEJTopRight > 0.0).any())) {
                        sEdgeInds.insert(indI - surfEdgeStartInd);
                    }
                }
            }
        }
    }

public: // helper functions
    int locateVoxelIndex(const Eigen::Matrix<double, 1, dim>& pos) const
    {
        Eigen::Array<int, 1, dim> voxelAxisIndex;
        locateVoxelAxisIndex(pos, voxelAxisIndex);
        return voxelAxisIndex2VoxelIndex(voxelAxisIndex.data());
    }
    void locateVoxelAxisIndex(const Eigen::Matrix<double, 1, dim>& pos,
        Eigen::Array<int, 1, dim>& voxelAxisIndex) const
    {
        voxelAxisIndex = ((pos - leftBottomCorner) * one_div_voxelSize).array().floor().template cast<int>();
    }
    int voxelAxisIndex2VoxelIndex(const int voxelAxisIndex[3]) const
    {
        return voxelAxisIndex2VoxelIndex(voxelAxisIndex[0], voxelAxisIndex[1], voxelAxisIndex[2]);
    }
    int voxelAxisIndex2VoxelIndex(int ix, int iy, int iz) const
    {
        assert(ix >= 0 && iy >= 0 && iz >= 0 && ix < voxelCount[0] && iy < voxelCount[1] && iz < voxelCount[2]);
        return ix + iy * voxelCount[0] + iz * voxelCount0x1;
    }
}; // class SpatialHash

} // namespace IPC

#endif // SpatialHash_hpp
