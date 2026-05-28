//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Andreas Apostolatos
//                   Pooyan Dadvand
//                   Philipp Bucher
//                   Nicolò Antonelli
//                   Andrea Gorgi
//

#if !defined(KRATOS_BREP_FACE_3D_H_INCLUDED )
#define  KRATOS_BREP_FACE_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/brep_curve_on_surface.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"

// trimming integration
#include "utilities/geometry_utilities/brep_trimming_utilities.h"
#include "utilities/geometry_utilities/brep_sbm_utilities.h"

namespace Kratos
{
class THBBrepSurface : public BrepSurface
{
public:

    using BaseType = BrepSurface;

    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        const IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) override;

private:

    THBSurfaceGeometry::Pointer mpTHBSurfaceGeometry;
};
}// namespace Kratos.

#endif // KRATOS_BREP_FACE_3D_H_INCLUDED  defined
