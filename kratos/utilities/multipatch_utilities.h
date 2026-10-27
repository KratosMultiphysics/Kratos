//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License: BSD License (kratos/license.txt)

#pragma once

// System includes
#include <array>
#include <map>
#include <utility>
#include <tuple>
#include <vector>

// Project includes
#include "includes/model_part.h"

namespace Kratos
{

/**
 * @brief Geometry-only activation and interface extraction for 2D FEM patches.
 * @details The background can contain P1 triangles and Q1 quadrilaterals.
 * Edge adjacency is cached independently of ACTIVE flags and box contours.
 * Call RebuildTopology after adding/removing elements or changing connectivity.
 * Boxes must be disjoint, strictly internal and aligned with element faces;
 * placement validation belongs to the patch manager. No elements are removed.
 */
class KRATOS_API(KRATOS_CORE) FemMultipatchUtilities
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(FemMultipatchUtilities);

    using IndexType = ModelPart::IndexType;
    using BoxType = std::array<double, 4>; // xmin, ymin, xmax, ymax
    using NodePairType = std::array<IndexType, 2>;
    using OwnerMapType = std::map<IndexType, IndexType>;
    using SurrogateEdgesType = std::map<IndexType, std::vector<NodePairType>>;
    using GaussPointsType = std::vector<std::pair<IndexType, std::array<double, 3>>>;

    explicit FemMultipatchUtilities(ModelPart& rBackground);

    /// Cache edge neighbours and element pointers for the current mesh.
    void RebuildTopology();

    /// Set ACTIVE and return inactive element ID -> one-based internal patch index.
    OwnerMapType UpdateBackgroundActivation(
        const std::vector<BoxType>& rBoxes,
        const double Tolerance);

    /// Return oriented background edge node IDs grouped by internal patch index.
    SurrogateEdgesType ExtractSurrogateEdges(const OwnerMapType& rInactiveOwners) const;

    /// Read each element's actual quadrature rule and map its points to global space.
    static GaussPointsType GetGaussPoints(const ModelPart& rModelPart, const bool ActiveOnly = true);

    /// Element ID, physical coordinates and physical quadrature weight.
    static std::vector<std::tuple<IndexType, std::array<double, 3>, double>>
    GetGaussPointsWithWeights(const ModelPart& rModelPart, bool ActiveOnly = true);

private:
    struct EdgeData
    {
        std::array<Element::Pointer, 2> Elements{};
        std::array<NodePairType, 2> Nodes{};
    };

    ModelPart& mrBackground;
    std::vector<Element::Pointer> mElements;
    std::vector<EdgeData> mEdges;
};

} // namespace Kratos
