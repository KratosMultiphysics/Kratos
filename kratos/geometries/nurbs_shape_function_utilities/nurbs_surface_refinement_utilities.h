//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_NURBS_SURFACE_REFINEMENT_UTILITIES_H_INCLUDED )
#define  KRATOS_NURBS_SURFACE_REFINEMENT_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/nurbs_surface_geometry.h"
#include "nurbs_utilities.h"
#include "includes/node.h"

namespace Kratos {

class KRATOS_API(KRATOS_CORE) NurbsSurfaceRefinementUtilities
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef Node NodeType;

    typedef NurbsSurfaceGeometry<3, PointerVector<NodeType>> NurbsSurfaceGeometryType;
    typedef typename NurbsSurfaceGeometryType::Pointer NurbsSurfaceGeometryPointerType;

    ///@}
    ///@name Operations
    ///@{

    /*
     * @brief Refines the u-knot vector of a NurbsSurfaceGeometry.
     * @details This function adopts the surface knot-refinement algorithm from
     *          Piegl 1995, Algorithm A5.5.
     *
     * Portet from carat++ (https://www.bgu.tum.de/st/software/forschung/carat/)
     * and ANurbs (https://github.com/oberbichler/ANurbs)
     *
     * @param rGeometry surface to be refined.
     * @param rKnotsUToInsert Knots to be inserted.
     * @param rPointsRefined the nodes for the refined geometry.
     * @param rKnotsURefined the new knot vector.
     * @param rWeightsRefined the new weight vector.
     */
    static void KnotRefinementU(
        NurbsSurfaceGeometryType& rGeometry,
        std::vector<double>& rKnotsUToInsert,
        PointerVector<NodeType>& rPointsRefined,
        Vector& rKnotsURefined,
        Vector& rWeightsRefined);

    /*
     * @brief Refines the v-knot vector of a NurbsSurfaceGeometry.
     * @details This function adopts the surface knot-refinement algorithm from
     *          Piegl 1995, Algorithm A5.5.
     *
     * Portet from carat++ (https://www.bgu.tum.de/st/software/forschung/carat/)
     * and ANurbs (https://github.com/oberbichler/ANurbs)
     *
     * @param rGeometry surface to be refined.
     * @param rKnotsVToInsert Knots to be inserted.
     * @param rPointsRefined the nodes for the refined geometry.
     * @param rKnotsVRefined the new knot vector.
     * @param rWeightsRefined the new weight vector.
     */
    static void KnotRefinementV(
        NurbsSurfaceGeometryType& rGeometry,
        std::vector<double>& rKnotsVToInsert,
        PointerVector<NodeType>& rPointsRefined,
        Vector& rKnotsVRefined,
        Vector& rWeightsRefined);

    /*
     * @brief Elevates the degree-u of a NurbsSurfaceGeometry.
     * @details This function adopts the surface degree elevation algorithm from
     *          Piegl 1995, Algorithm A5.10.
     *
     * Portet from carat++ (https://www.bgu.tum.de/st/software/forschung/carat/)
     *
     * @param rGeometry surface to be refined.
     * @param rDegreeUToElevate Degree to be elevated.
     * @param rPointsRefined the nodes for the refined geometry.
     * @param rKnotsURefined the new knot vector.
     * @param rWeightsRefined the new weight vector.
     */
    static void DegreeElevationU(
        NurbsSurfaceGeometryType& rGeometry,
        SizeType& rDegreeUToElevate,
        PointerVector<NodeType>& rPointsRefined,
        Vector& rKnotsURefined,
        Vector& rWeightsRefined);

    /*
     * @brief Elevates the degree-v of a NurbsSurfaceGeometry.
     * @details This function adopts the surface degree elevation algorithm from
     *          Piegl 1995, Algorithm A5.10.
     *
     * Portet from carat++ (https://www.bgu.tum.de/st/software/forschung/carat/)
     *
     * @param rGeometry surface to be refined.
     * @param rDegreeVToElevate Degree to be elevated.
     * @param rPointsRefined the nodes for the refined geometry.
     * @param rKnotsVRefined the new knot vector.
     * @param rWeightsRefined the new weight vector.
     */
    static void DegreeElevationV(
        NurbsSurfaceGeometryType& rGeometry,
        SizeType& rDegreeVToElevate,
        PointerVector<NodeType>& rPointsRefined,
        Vector& rKnotsVRefined,
        Vector& rWeightsRefined);

    ///@}
    ///@name Utilities
    ///@{

    /*
     * @brief Sorts the knot vector, which is required for the refinement.
     *        Filters the knots to insert which are within the former knot span
     *        and accordingly coordinates of the former geometry.
     *
     * @param rKnotsToInsert the knot vector to be sorted and filtered.
     * @param rKnotsOld the old knot vector, defining the limits.
     */
    static void SortAndFilter(
        std::vector<double>&rKnotsToInsert,
        const Vector& rKnotsOld);

};

} // End namespace kratos

#endif // KRATOS_NURBS_SURFACE_REFINEMENT_UTILITIES_H_INCLUDED