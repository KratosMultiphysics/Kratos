//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License: BSD License (kratos/license.txt)
//

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"

namespace Kratos {

/**
 * @brief Samples the physical support of a captured two-fluid Q1 phase.
 * @details The input is a snapshot, not a reference to the live model. Gap
 * polygons retain their Q1 parent, so interpolation may extend beyond that
 * parent. At a seam each distinct parent contributes exactly one trace.
 */
class TwoFluidPhaseVelocitySampler
{
public:
    using PointType = std::array<double, 2>;
    using VelocityType = std::array<double, 3>;
    using NodeType = std::tuple<std::size_t, double, double, VelocityType>;
    using ParentType = std::array<NodeType, 4>;
    using CellType = std::tuple<std::vector<PointType>, ParentType, std::string>;
    using CellsType = std::vector<CellType>;

    TwoFluidPhaseVelocitySampler(CellsType Cells, std::string PhaseName, double Tolerance)
        : mCells(std::move(Cells)), mPhaseName(std::move(PhaseName)), mTolerance(Tolerance)
    {
        KRATOS_ERROR_IF_NOT(std::isfinite(mTolerance) && mTolerance > 0.0)
            << "Phase velocity sampling requires a positive finite tolerance." << std::endl;
        for (const auto& r_cell : mCells) {
            const auto& r_parent = std::get<1>(r_cell);
            for (std::size_t i = 0; i < 4; ++i) {
                const auto& r_a = r_parent[i];
                const auto& r_b = r_parent[(i + 1) % 4];
                mWidth = std::max(mWidth, std::hypot(std::get<1>(r_b) - std::get<1>(r_a),
                                                     std::get<2>(r_b) - std::get<2>(r_a)));
            }
        }
        if (mWidth == 0.0) {
            return;
        }
        for (std::size_t i = 0; i < mCells.size(); ++i) {
            const auto& r_polygon = std::get<0>(mCells[i]);
            KRATOS_ERROR_IF(r_polygon.empty()) << "Empty historical phase polygon." << std::endl;
            double xmin = r_polygon[0][0], xmax = xmin;
            double ymin = r_polygon[0][1], ymax = ymin;
            for (const auto& r_point : r_polygon) {
                xmin = std::min(xmin, r_point[0]);
                xmax = std::max(xmax, r_point[0]);
                ymin = std::min(ymin, r_point[1]);
                ymax = std::max(ymax, r_point[1]);
            }
            for (long ix = Bin(xmin); ix <= Bin(xmax); ++ix) {
                for (long iy = Bin(ymin); iy <= Bin(ymax); ++iy) {
                    mBins[{ix, iy}].push_back(i);
                }
            }
        }
    }

    VelocityType Sample(const double X, const double Y) const
    {
        KRATOS_TRY
        KRATOS_ERROR_IF_NOT(std::isfinite(X) && std::isfinite(Y))
            << "Phase velocity sampling requires a finite point." << std::endl;
        std::set<std::size_t> candidates;
        if (mWidth == 0.0) {
            for (std::size_t i = 0; i < mCells.size(); ++i) {
                candidates.insert(i);
            }
        } else {
            // Match the Python history index, including its one-parent-width
            // search radius and original cell ordering.
            const double radius = mWidth + mTolerance;
            for (long ix = Bin(X - radius); ix <= Bin(X + radius); ++ix) {
                for (long iy = Bin(Y - radius); iy <= Bin(Y + radius); ++iy) {
                    const auto it = mBins.find({ix, iy});
                    if (it != mBins.end()) {
                        candidates.insert(it->second.begin(), it->second.end());
                    }
                }
            }
        }

        using ParentKeyType = std::array<std::size_t, 4>;
        std::map<ParentKeyType, std::size_t> parent_positions;
        std::vector<const ParentType*> parents;
        for (const auto index : candidates) {
            const auto& r_cell = mCells[index];
            if (!Contains(std::get<0>(r_cell), X, Y)) {
                continue;
            }
            const auto& r_parent = std::get<1>(r_cell);
            ParentKeyType key;
            for (std::size_t i = 0; i < 4; ++i) {
                key[i] = std::get<0>(r_parent[i]);
            }
            const auto insertion = parent_positions.emplace(key, parents.size());
            if (insertion.second) {
                parents.push_back(&r_parent);
            } else {
                // Python dictionaries retain the first position and replace
                // the value when multiple polygons have the same parent key.
                parents[insertion.first->second] = &r_parent;
            }
        }
        KRATOS_ERROR_IF(parents.empty()) << "No physical " << mPhaseName
            << " velocity support at (" << X << ", " << Y << ")" << std::endl;
        VelocityType result{0.0, 0.0, 0.0};
        for (const auto* p_parent : parents) {
            const auto value = Interpolate(*p_parent, X, Y);
            for (std::size_t d = 0; d < 3; ++d) {
                result[d] += value[d];
            }
        }
        for (double& r_value : result) {
            r_value /= parents.size();
        }
        return result;
        KRATOS_CATCH("")
    }

    /** Sample all level-set nodes without modifying the model part. */
    static std::vector<VelocityType> SampleNodes(
        const ModelPart& rModelPart,
        const TwoFluidPhaseVelocitySampler& rWater,
        const TwoFluidPhaseVelocitySampler& rAir,
        const double Tolerance)
    {
        std::vector<VelocityType> samples;
        samples.reserve(rModelPart.NumberOfNodes());
        for (const auto& r_node : rModelPart.Nodes()) {
            const double distance = r_node.FastGetSolutionStepValue(DISTANCE);
            VelocityType value;
            if (std::abs(distance) <= Tolerance) {
                const auto water = rWater.Sample(r_node.X(), r_node.Y());
                const auto air = rAir.Sample(r_node.X(), r_node.Y());
                ValidateVelocity(water, r_node.Id());
                ValidateVelocity(air, r_node.Id());
                value = {(water[0] + air[0]) / 2.0, (water[1] + air[1]) / 2.0, 0.0};
            } else {
                value = (distance < 0.0 ? rWater : rAir).Sample(r_node.X(), r_node.Y());
                ValidateVelocity(value, r_node.Id());
                value[2] = 0.0;
            }
            samples.push_back(value);
        }
        return samples;
    }

private:
    CellsType mCells;
    std::string mPhaseName;
    double mTolerance;
    double mWidth = 0.0;
    std::map<std::pair<long, long>, std::vector<std::size_t>> mBins;

    long Bin(const double Coordinate) const
    {
        return static_cast<long>(std::floor(Coordinate / mWidth));
    }

    bool Contains(const std::vector<PointType>& rPolygon, const double X, const double Y) const
    {
        bool positive = false, negative = false;
        double distance = std::numeric_limits<double>::infinity();
        const double sign_tolerance = mTolerance * mTolerance;
        for (std::size_t i = 0; i < rPolygon.size(); ++i) {
            const auto& r_a = rPolygon[i];
            const auto& r_b = rPolygon[(i + 1) % rPolygon.size()];
            const double dx = r_b[0] - r_a[0], dy = r_b[1] - r_a[1];
            const double length2 = dx * dx + dy * dy;
            KRATOS_ERROR_IF(length2 <= 0.0) << "Degenerate historical phase polygon" << std::endl;
            const double cross = dx * (Y - r_a[1]) - dy * (X - r_a[0]);
            positive = positive || cross > sign_tolerance;
            negative = negative || cross < -sign_tolerance;
            const double fraction = std::max(0.0, std::min(1.0,
                ((X - r_a[0]) * dx + (Y - r_a[1]) * dy) / length2));
            distance = std::min(distance, std::hypot(X - r_a[0] - fraction * dx,
                                                     Y - r_a[1] - fraction * dy));
        }
        return !positive || !negative || distance <= mTolerance;
    }

    static void ValidateVelocity(const VelocityType& rVelocity, const std::size_t NodeId)
    {
        if (!std::isfinite(rVelocity[0]) || !std::isfinite(rVelocity[1])) {
            throw std::invalid_argument("Invalid phase velocity at node " + std::to_string(NodeId));
        }
    }

    static VelocityType Interpolate(const ParentType& rParent, const double X, const double Y)
    {
        const auto& a = rParent[0];
        const auto& b = rParent[1];
        const auto& d = rParent[3];
        const double ex = std::get<1>(b) - std::get<1>(a), ey = std::get<2>(b) - std::get<2>(a);
        const double fx = std::get<1>(d) - std::get<1>(a), fy = std::get<2>(d) - std::get<2>(a);
        const double determinant = ex * fy - ey * fx;
        KRATOS_ERROR_IF(std::abs(determinant) <= 1e-30) << "Degenerate phase-history parent" << std::endl;
        const double s = ((X - std::get<1>(a)) * fy - (Y - std::get<2>(a)) * fx) / determinant;
        const double t = (ex * (Y - std::get<2>(a)) - ey * (X - std::get<1>(a))) / determinant;
        const std::array<double, 4> shape{(1-s)*(1-t), s*(1-t), s*t, (1-s)*t};
        VelocityType value{0.0, 0.0, 0.0};
        for (std::size_t i = 0; i < 4; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                value[j] += shape[i] * std::get<3>(rParent[i])[j];
            }
        }
        return value;
    }
};
} // namespace Kratos
