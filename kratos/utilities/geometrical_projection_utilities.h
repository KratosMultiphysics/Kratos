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
#include "utilities/math_utils.h"
#include "includes/model_part.h"

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
     * @param Geom The geometry where to be projected
     * @param PointDestiny The point to be projected
     * @param PointProjected The point pojected over the plane
     * @param Normal The normal of the geometry
     * @param Vector The direction to project
     * @return Distance The distance between surfaces
     */

    static inline double FastProjectDirection(
        const GeometryType& Geom,
        const PointType& PointDestiny,
        PointType& PointProjected,
        const array_1d<double,3>& Normal,
        const array_1d<double,3>& Vector
        )
    {
        // Zero tolerance
        const double zero_tolerance = std::numeric_limits<double>::epsilon();

        // We define the distance
        double distance = 0.0;
        
        const array_1d<double,3> vector_points = Geom[0].Coordinates() - PointDestiny.Coordinates();

        if( norm_2( Vector ) < zero_tolerance && norm_2( Normal ) > zero_tolerance ) {
            distance = inner_prod(vector_points, Normal)/norm_2(Normal);

            PointProjected.Coordinates() = PointDestiny.Coordinates() + Vector * distance;
            KRATOS_WARNING("Warning: Zero projection vector.") << " Projection using the condition vector instead." << std::endl;
        } else if (std::abs(inner_prod(Vector, Normal) ) > zero_tolerance) {
            distance = inner_prod(vector_points, Normal)/inner_prod(Vector, Normal); 

            PointProjected.Coordinates() = PointDestiny.Coordinates() + Vector * distance;
        } else {
            PointProjected.Coordinates() = PointDestiny.Coordinates();
            KRATOS_WARNING("Warning: The line and the plane are coplanar.")  << " Something wrong happened " << std::endl;
        }
        
        return distance;
    }
    
    /**
     * @brief Project a point over a plane (avoiding some steps)
     * @param PointOrigin A point in the plane
     * @param PointDestiny The point to be projected
     * @param Normal The normal of the plane
     * @param Distance The distance to the projection
     * @return PointProjected The point pojected over the plane
     */
    
    static inline PointType FastProject(
        const PointType& PointOrigin,
        const PointType& PointDestiny,
        const array_1d<double,3>& Normal,
        double& Distance
        )
    {
        const array_1d<double,3> vector_points = PointDestiny.Coordinates() - PointOrigin.Coordinates();

        Distance = inner_prod(vector_points, Normal); 
        
        PointType point_projected;
        noalias(point_projected.Coordinates()) = PointDestiny.Coordinates() - Normal * Distance;
        
        return point_projected;
    }
    
    /**
     * @brief Projects iteratively to get the coordinate
     * @param GeomOrigin The origin geometry
     * @param PointDestiny The destination point
     * @param ResultingPoint The distance between the point and the plane
     * @return Inside True is inside, false not
     */
    
    static inline bool ProjectIterativeLine2D(
        GeometryType& GeomOrigin,
        const GeometryType::CoordinatesArrayType& PointDestiny,
        GeometryType::CoordinatesArrayType& ResultingPoint,
        const array_1d<double, 3>& Normal,
        const double Tolerance = 1.0e-8,
        double DeltaXi = 0.5
        )
    {
//         ResultingPoint.clear();
        
        double old_delta_xi = 0.0;

        array_1d<double, 3> current_global_coords;

        array_1d<array_1d<double, 3>, 2> normals;
        normals[0] = GeomOrigin[0].FastGetSolutionStepValue(NORMAL);
        normals[1] = GeomOrigin[1].FastGetSolutionStepValue(NORMAL);
        
        BoundedMatrix<double,2,2> X;
        BoundedMatrix<double,2,1> DN;
        for(unsigned int i=0; i<2;++i) {
            X(0,i) = GeomOrigin[i].X();
            X(1,i) = GeomOrigin[i].Y();
        }

        Matrix J = ZeroMatrix( 1, 1 );
        
        //Newton iteration:

        const unsigned int max_iter = 20;

        for ( unsigned int k = 0; k < max_iter; ++k ) {
            array_1d<double, 2> N_origin;
            N_origin[0] = 0.5 * ( 1.0 - ResultingPoint[0]);
            N_origin[1] = 0.5 * ( 1.0 + ResultingPoint[0]);
            
            array_1d<double,3> normal_xi = N_origin[0] * normals[0] + N_origin[1] * normals[1];
            normal_xi = normal_xi/norm_2(normal_xi); 
            
            current_global_coords = ZeroVector(3);
            for( IndexType i_node = 0; i_node < 2; ++i_node )
                current_global_coords += N_origin[i_node] * GeomOrigin[i_node].Coordinates(); 
            
            const array_1d<double,3> VectorPoints = GeomOrigin.Center() - PointDestiny;
            const double distance = inner_prod(VectorPoints, Normal)/inner_prod(-normal_xi, Normal); 
            const array_1d<double, 3> current_destiny_global_coords = PointDestiny - normal_xi * distance;
            
            // Derivatives of shape functions
            Matrix ShapeFunctionsGradients;
            ShapeFunctionsGradients = GeomOrigin.ShapeFunctionsLocalGradients(ShapeFunctionsGradients, ResultingPoint );
            noalias(DN) = prod(X,ShapeFunctionsGradients);

            noalias(J) = prod(trans(DN),DN); // TODO: Add the non linearity concerning the normal
            Vector RHS = prod(trans(DN),subrange(current_destiny_global_coords - current_global_coords,0,2));
            
            old_delta_xi = DeltaXi;
            DeltaXi = RHS[0]/J(0, 0);
            
            ResultingPoint[0] += DeltaXi;
            
            if (ResultingPoint[0] <= -1.0)
                ResultingPoint[0] = -1.0;
            else if (ResultingPoint[0] >= 1.0)
                ResultingPoint[0] = 1.0;
            
            if ( std::abs(DeltaXi - old_delta_xi) < Tolerance )
                return true;
        }
        
        return false;
    }
    
private:
};// class GeometricalProjectionUtilities

}
#endif /* KRATOS_GEOMETRICAL_PROJECTION_UTILITIES defined */
