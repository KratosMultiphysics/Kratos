//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


#if !defined(KRATOS_TESSELLATION_GEOMETRY_INTERFACE_H_INCLUDED )
#define  KRATOS_TESSELLATION_GEOMETRY_INTERFACE_H_INCLUDED

#include "tessellation.h"
#include "curve_tessellation.h"
#include "geometries/geometry.h"

namespace Kratos {

template <class TContainerPointType>
class TessellationGeometryInterface
{
public:

    typedef Geometry<typename TContainerPointType::value_type> GeometryType;
    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    ///@name Constructor
    ///@{

    /// Conctructor for tessellation of a nurbs curve
    TessellationGeometryInterface(
        const GeometryType& rGeometry,
        double Tolerance,
        const int NumberOfGuessesPerInterval = 1)
    {
        if (rGeometry.Dimension() == 1) {
            mTessellation = CurveTessellation<TContainerPointType>();
            mTessellation.Tessellate(rGeometry, Tolerance, NumberOfGuessesPerInterval);
        }
        else {
            KRATOS_ERROR << "Other dimensions than 1 are not implemented yet." << std::endl;
        }
    }

    ///@}
    ///@name Access Functions
    ///@{

    /// This method returns closest point within the tesselation.
    void Tessellate(
        const GeometryType& rGeometry,
        const double Tolerance,
        const int NumberOfGuessesPerInterval = 1)
    {
        mTessellation.Tessellate(
            rGeometry, Tolerance, NumberOfGuessesPerInterval);
    }

    /// This method returns closest point within the tesselation.
    void GetClosestPoint(
        const CoordinatesArrayType& rPointGlobalCoordinates,
        CoordinatesArrayType& rClosestPointGlobalCoordinates,
        CoordinatesArrayType& rClosestPointLocalCoordinates)
    {
        mTessellation.GetClosestPoint(
            rPointGlobalCoordinates,
            rClosestPointGlobalCoordinates,
            rClosestPointLocalCoordinates);
    }

    ///@}
private:
    ///@name Private Member Variables
    ///@{

    CurveTessellation<TContainerPointType> mTessellation;

    ///@}
};

} // namespace CurveTessellation

#endif // KRATOS_TESSELLATION_GEOMETRY_INTERFACE_H_INCLUDED defined
