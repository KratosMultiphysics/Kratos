//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

#if !defined(KRATOS_NURBS_VOLUME_REFINEMENT_UTILITIES_H_INCLUDED )
#define  KRATOS_NURBS_VOLUME_REFINEMENT_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/nurbs_volume_geometry.h"
#include "nurbs_utilities.h"
#include "includes/node.h"

namespace Kratos {

class NurbsVolumeRefinementUtilities
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef Node NodeType;

    typedef NurbsVolumeGeometry<PointerVector<NodeType>> NurbsVolumeGeometryType;
    typedef typename NurbsVolumeGeometryType::Pointer NurbsVolumeGeometryPointerType;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Refines the u-knotvector of a NurbsVolumeGeometry.
     * @details This function adopts the surface knot-refinement algorithm from Piegl1995 (p.164 Algorithm A5.5).
     *          The algorithm is modified to suit trivariant B-Spline volumes.
     * @param rGeometry surface to be refined.
     * @param rKnotsUToInsert Knots to be inserted.
     * @param rPointsRefined the nodes for the refined geometry.
     * @param rKnotsURefined the new knot vector.
     * @note This function does not consider weights, thus only B-Spline-Volumes can be refined.
     **/
    static void KnotRefinementU(NurbsVolumeGeometryType& rGeometry, std::vector<double>& rKnotsUToInsert,
                            PointerVector<NodeType>& rPointsRefined, Vector& rKnotsURefined );

    static void KnotRefinementV(NurbsVolumeGeometryType& rGeometry, std::vector<double>& rKnotsVToInsert,
                                PointerVector<NodeType>& rPointsRefined, Vector& rKnotsVRefined );

    static void KnotRefinementW(NurbsVolumeGeometryType& rGeometry, std::vector<double>& rKnotsWToInsert,
                                PointerVector<NodeType>& rPointsRefined, Vector& rKnotsWRefined );

    static void DegreeElevationU(NurbsVolumeGeometryType& rGeometry,SizeType& rDegreeUToElevate,
                                PointerVector<NodeType>& rPointsRefined, Vector& rKnotsURefined, Vector& rWeightsRefined);
    
    static void DegreeElevationV(NurbsVolumeGeometryType& rGeometry,SizeType& rDegreeVToElevate,
                                PointerVector<NodeType>& rPointsRefined, Vector& rKnotsVRefined, Vector& rWeightsRefined);
    
    static void DegreeElevationW(NurbsVolumeGeometryType& rGeometry,SizeType& rDegreeWToElevate,
                                PointerVector<NodeType>& rPointsRefined, Vector& rKnotsWRefined, Vector& rWeightsRefined);

    ///@}
};

} // End namespace kratos

#endif // KRATOS_NURBS_VOLUME_REFINEMENT_UTILITIES_H_INCLUDED