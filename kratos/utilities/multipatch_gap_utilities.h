//     |  /           |
//     ' /   __| _` | __|  _ \   __|
//     . \  |   (   | |   (   |\__ `
//    _|\_\_|  \__,_|\__|\___/ ____/
//                    Multi-Physics
//
//   License:         BSD License
//                    Kratos default license: kratos/license.txt
//
//   Main author:     Nicolò Antonelli
//

#pragma once

#include <array>
#include <vector>
#include "utilities/multipatch_utilities.h"

namespace Kratos {
/** @brief Conservative, positive-gap geometry for Cartesian Q1 backgrounds.
 * Contours are recovered from independent fitted meshes, never from boxes.
 * The initial implementation requires convex contours (the square use case).
 * Polygon vertices used for quadrature/visualization carry no solution DOFs.
 */
class KRATOS_API(KRATOS_CORE) FemMultipatchGapUtilities
{
public:
    using Coordinates = array_1d<double, 3>;
    struct Segment {
        Element::Pointer pFirst;
        Element::Pointer pSecond;
        Coordinates Start;
        Coordinates End;
    };
    struct GapCell {
        Element::Pointer pParent;
        std::array<Coordinates, 3> Vertices;
        Geometry<Node>::Pointer pQuadrature;
    };
    struct GapData {
        std::vector<GapCell> Cells;
        std::vector<Segment> Surrogate;
        std::vector<Segment> Physical;
        std::vector<Segment> Connectors;
    };

    static FemMultipatchUtilities::OwnerMapType UpdateActivation(
        ModelPart& rBackground, const std::vector<ModelPart*>& rInnerPatches,
        double Tolerance);
    static GapData Build(ModelPart& rBackground, ModelPart& rInner, double Tolerance);
    /// Materialize independent views: TrueBoundary, SurrogateBoundary,
    /// GapVolumes, GapTriangles, PhysicalInterfaces, GapInterfaces.
    static void CreateGeometry(ModelPart& rBackground, ModelPart& rInner,
                               ModelPart& rInterface, double Tolerance);
};
} // namespace Kratos
