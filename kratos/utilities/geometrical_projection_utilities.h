//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_GEOMETRICAL_PROJECTION_UTILITIES)
#define KRATOS_GEOMETRICAL_PROJECTION_UTILITIES

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/checks.h"
#include "includes/node.h"
#include "includes/variables.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class GeometricalProjectionUtilities
 * @ingroup KratosCore
 * @brief This is a class that provides auxiliar utilities for projections
 * @details This is a class that provides auxiliar utilities for the projections. Check the documentation for more details
 * @author Vicente Mataix Ferrandiz
 * Contact: vmataix@cimne.upc.edu
 */
class GeometricalProjectionUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometricalProjectionUtilities
    KRATOS_CLASS_POINTER_DEFINITION( GeometricalProjectionUtilities );

    // Some geometrical definitions
    typedef Node<3>                                              NodeType;
    typedef Point                                               PointType;
    typedef PointType::CoordinatesArrayType          CoordinatesArrayType;

    /// Definition of geometries
    typedef Geometry<NodeType>                               GeometryType;
    typedef Geometry<PointType>                         GeometryPointType;

    /// Index type definition
    typedef std::size_t                                         IndexType;

    /// Size type definition
    typedef std::size_t                                          SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Project a point over a line/plane following an arbitrary direction
     * @param rGeom The geometry where to be projected
     * @param rPointDestiny The point to be projected
     * @param rPointProjected The point pojected over the plane
     * @param rNormal The normal of the geometry
     * @param rVector The direction to project
     * @return Distance The distance between surfaces
     */

    static inline double FastProjectDirection(
        const GeometryType& rGeom,
        const PointType& rPointDestiny,
        PointType& rPointProjected,
        const array_1d<double,3>& rNormal,
        const array_1d<double,3>& rVector
        )
    {
        // Zero tolerance
        const double zero_tolerance = std::numeric_limits<double>::epsilon();

        // We define the distance
        double distance = 0.0;

        const array_1d<double,3> vector_points = rGeom[0].Coordinates() - rPointDestiny.Coordinates();

        if( norm_2( rVector ) < zero_tolerance && norm_2( rNormal ) > zero_tolerance ) {
            distance = inner_prod(vector_points, rNormal)/norm_2(rNormal);
            noalias(rPointProjected.Coordinates()) = rPointDestiny.Coordinates() + rVector * distance;
            KRATOS_WARNING("Warning: Zero projection vector.") << " Projection using the condition vector instead." << std::endl;
        } else if (std::abs(inner_prod(rVector, rNormal) ) > zero_tolerance) {
            distance = inner_prod(vector_points, rNormal)/inner_prod(rVector, rNormal);
            noalias(rPointProjected.Coordinates()) = rPointDestiny.Coordinates() + rVector * distance;
        } else {
            noalias(rPointProjected.Coordinates()) = rPointDestiny.Coordinates();
            KRATOS_WARNING("Warning: The line and the plane are coplanar.")  << " Something wrong happened " << std::endl;
        }

        return distance;
    }

    /**
     * @brief Project a point over a plane (avoiding some steps)
     * @param rPointOrigin A point in the plane
     * @param rPointDestiny The point to be projected
     * @param rNormal The normal of the plane
     * @param rDistance The distance to the projection
     * @return PointProjected The point pojected over the plane
     */
    static inline PointType FastProject(
        const PointType& rPointOrigin,
        const PointType& rPointDestiny,
        const array_1d<double,3>& rNormal,
        double& rDistance
        )
    {
        const array_1d<double,3> vector_points = rPointDestiny.Coordinates() - rPointOrigin.Coordinates();

        rDistance = inner_prod(vector_points, rNormal);

        PointType point_projected;
    #ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
        point_projected.Coordinates() = rPointDestiny.Coordinates() - rNormal * rDistance;
    #else
        noalias(point_projected.Coordinates()) = rPointDestiny.Coordinates() - rNormal * rDistance;
    #endif // ifdef KRATOS_USE_AMATRIX

        return point_projected;
    }

    /**
     * @brief Projects iteratively to get the coordinate
     * @param rGeomOrigin The origin geometry
     * @param rPointDestiny The destination point
     * @param rResultingPoint The distance between the point and the plane
     * @return Inside True is inside, false not
     */
    static inline bool ProjectIterativeLine2D(
        GeometryType& rGeomOrigin,
        const GeometryType::CoordinatesArrayType& rPointDestiny,
        GeometryType::CoordinatesArrayType& rResultingPoint,
        const array_1d<double, 3>& rNormal,
        const double Tolerance = 1.0e-8,
        double DeltaXi = 0.5
        )
    {
//         rResultingPoint.clear();

        double old_delta_xi = 0.0;

        array_1d<double, 3> current_global_coords;

        KRATOS_DEBUG_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL, rGeomOrigin[0])
        KRATOS_DEBUG_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL, rGeomOrigin[1])

        array_1d<array_1d<double, 3>, 2> normals;
        normals[0] = rGeomOrigin[0].FastGetSolutionStepValue(NORMAL);
        normals[1] = rGeomOrigin[1].FastGetSolutionStepValue(NORMAL);

        BoundedMatrix<double,2,2> X;
        BoundedMatrix<double,2,1> DN;
        for(unsigned int i=0; i<2;++i) {
            X(0,i) = rGeomOrigin[i].X();
            X(1,i) = rGeomOrigin[i].Y();
        }

        Matrix J = ZeroMatrix( 1, 1 );

        //Newton iteration:

        const unsigned int max_iter = 20;

        for ( unsigned int k = 0; k < max_iter; ++k ) {
            array_1d<double, 2> N_origin;
            N_origin[0] = 0.5 * ( 1.0 - rResultingPoint[0]);
            N_origin[1] = 0.5 * ( 1.0 + rResultingPoint[0]);

            array_1d<double,3> normal_xi = N_origin[0] * normals[0] + N_origin[1] * normals[1];
            normal_xi = normal_xi/norm_2(normal_xi);

            current_global_coords = ZeroVector(3);
            for( IndexType i_node = 0; i_node < 2; ++i_node )
                current_global_coords += N_origin[i_node] * rGeomOrigin[i_node].Coordinates();

            const array_1d<double,3> VectorPoints = rGeomOrigin.Center() - rPointDestiny;
            const double distance = inner_prod(VectorPoints, rNormal)/inner_prod(-normal_xi, rNormal);
            const array_1d<double, 3> current_destiny_global_coords = rPointDestiny - normal_xi * distance;

            // Derivatives of shape functions
            Matrix ShapeFunctionsGradients;
            ShapeFunctionsGradients = rGeomOrigin.ShapeFunctionsLocalGradients(ShapeFunctionsGradients, rResultingPoint );

        #ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
            DN = prod(X,ShapeFunctionsGradients);

            J = prod(trans(DN),DN); // TODO: Add the non linearity concerning the normal
        #else
            noalias(DN) = prod(X,ShapeFunctionsGradients);

            noalias(J) = prod(trans(DN),DN); // TODO: Add the non linearity concerning the normal
        #endif // ifdef KRATOS_USE_AMATRIX

            const Vector RHS = prod(trans(DN),subrange(current_destiny_global_coords - current_global_coords,0,2));

            old_delta_xi = DeltaXi;
            DeltaXi = RHS[0]/J(0, 0);

            rResultingPoint[0] += DeltaXi;

            if (rResultingPoint[0] <= -1.0)
                rResultingPoint[0] = -1.0;
            else if (rResultingPoint[0] >= 1.0)
                rResultingPoint[0] = 1.0;

            if ( std::abs(DeltaXi - old_delta_xi) < Tolerance )
                return true;
        }

        return false;
    }

private:
};// class GeometricalProjectionUtilities

}
#endif /* KRATOS_GEOMETRICAL_PROJECTION_UTILITIES defined */
