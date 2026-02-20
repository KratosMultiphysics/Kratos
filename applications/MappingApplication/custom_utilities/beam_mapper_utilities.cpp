//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
// Main authors:    E. G. Loera Villeda
// Contributor:     Juan I. Camarotti
//

// System includes

// External includes

// Project includes
#include "beam_mapper_utilities.h"

namespace Kratos
{
namespace BeamMapperUtilities
{
    using GeometryType = Geometry<Node>;
    
    void HermitianShapeFunctionsValues(
        Vector& rHermitianShapeFunctions, 
        Vector& rHermitianShapeFunctionsDer, 
        const array_1d<double, 3>& rCoordinates) 
    {
        if(rHermitianShapeFunctions.size() != 4) {
            rHermitianShapeFunctions.resize(4, false);
        }

        rHermitianShapeFunctions[0] =  0.25 * ( 1.0 - rCoordinates[0]) * ( 1.0 - rCoordinates[0]) * ( 2.0 + rCoordinates[0]);
        rHermitianShapeFunctions[1] =  0.125 * ( 1.0 - rCoordinates[0]) * ( 1.0 - rCoordinates[0]) * ( 1.0 + rCoordinates[0]);
        rHermitianShapeFunctions[2] =  0.25 * ( 1.0 + rCoordinates[0]) * ( 1.0 + rCoordinates[0]) * ( 2.0 - rCoordinates[0]);
        rHermitianShapeFunctions[3] =  -0.125  * ( 1.0 + rCoordinates[0]) * ( 1.0 + rCoordinates[0]) * ( 1.0 - rCoordinates[0]);

        if(rHermitianShapeFunctionsDer.size() != 4) {
            rHermitianShapeFunctionsDer.resize(4, false);
        }

        rHermitianShapeFunctionsDer[0] = -(3.0/2.0) * ( 1.0 - rCoordinates[0] ) * ( 1.0 + rCoordinates[0] ); 
        rHermitianShapeFunctionsDer[1] = -0.25 * ( 1.0 - rCoordinates[0] ) * ( 1.0 + 3 * rCoordinates[0]) ;
        rHermitianShapeFunctionsDer[2] = (3.0/2.0) * ( 1.0 + rCoordinates[0] ) * ( 1.0 - rCoordinates[0] );
        rHermitianShapeFunctionsDer[3] = -0.25 * ( 1.0 + rCoordinates[0] ) * ( 1.0 - 3 * rCoordinates[0]) ;    
    }

    ProjectionUtilities::PairingIndex ProjectOnLineHermitian(
        const GeometryType& rGeometry,
        const Point& rPointToProject,
        const double LocalCoordTol,
        Vector& rHermitianShapeFunctionValues,
        Vector& rHermitianShapeFunctionValuesDer,
        double& rProjectionDistance,
        Point& rProjectionOfPoint)
    {
        rProjectionDistance = std::abs(GeometricalProjectionUtilities::FastProjectOnLine(rGeometry, rPointToProject, rProjectionOfPoint));
        array_1d<double, 3> local_coords;
        ProjectionUtilities::PairingIndex pairing_index = ProjectionUtilities::PairingIndex::Unspecified;

        if (rGeometry.IsInside(rProjectionOfPoint, local_coords, 1e-14)) {
            pairing_index = ProjectionUtilities::PairingIndex::Line_Inside;
            HermitianShapeFunctionsValues(rHermitianShapeFunctionValues, rHermitianShapeFunctionValuesDer, local_coords);
        } 

        return pairing_index;
    }
} // namespace BeamMapperUtilities

} // namespace Kratos.
