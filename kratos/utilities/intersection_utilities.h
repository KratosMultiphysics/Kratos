//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_INTERSECTION_UTILITIES )
#define  KRATOS_INTERSECTION_UTILITIES

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/node.h"
#include "utilities/math_utils.h"

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
 * @class IntersectionUtilities
 * @ingroup KratosCore
 * @brief Utilities to compute intersections between different geometries
 * @details This class provides static methods to check if there is 
 * intersections between different entities, and if there is, give back
 * the intersection points.
 * @author Ruben Zorrilla
 */
class IntersectionUtilities
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of IntersectionUtilities
    KRATOS_CLASS_POINTER_DEFINITION( IntersectionUtilities );
    
    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Default constructor
     */
    IntersectionUtilities(){}

    /// Destructor
    virtual ~IntersectionUtilities(){}

    ///@}
    ///@name Operations
    ///@{

	/**
     * Find the 3D intersection of a ray with a triangle
     * @param rTriangleGeometry Is the triangle to intersect
     * @param rLinePoint1 Coordinates of the first point of the intersecting line
     * @param rLinePoint2 Coordinates of the second point of the intersecting line
     * @return rIntersectionPoint The intersection point coordinates
     * @return The intersection type index: 
     * -1 (the triangle is degenerate)
     * 0 (disjoint - no intersection) 
     * 1 (intersect in a unique point)
     * 2 (are in the same plane)
     */
	template <class TGeometryType>
	static int ComputeTriangleLineIntersection(
        const TGeometryType& rTriangleGeometry, 
        const array_1d<double,3>& rLinePoint1, 
        const array_1d<double,3>& rLinePoint2, 
        array_1d<double,3>& rIntersectionPoint,
		const double epsilon = 1e-12) {

		// This is the adaption of the implemnetation provided in:
        // http://www.softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm#intersect_RayTriangle()

		// Get triangle edge vectors and plane normal
		const array_1d<double,3> u = rTriangleGeometry[1] - rTriangleGeometry[0];
		const array_1d<double,3> v = rTriangleGeometry[2] - rTriangleGeometry[0];
		array_1d<double,3> n;
		MathUtils<double>::CrossProduct<array_1d<double,3>,array_1d<double,3>,array_1d<double,3>>(n,u,v);

		// Check if the triangle is degenerate (do not deal with this case)
		if (MathUtils<double>::Norm3(n) < epsilon){
			return -1;
		}

		const array_1d<double,3> dir = rLinePoint2 - rLinePoint1; // Edge direction vector
		const array_1d<double,3> w_0 = rLinePoint1 - rTriangleGeometry[0];
		const double a = -inner_prod(n,w_0);
		const double b = inner_prod(n,dir);

		// Check if the ray is parallel to the triangle plane
		if (std::abs(b) < epsilon){
			if (a == 0.0){	
				return 2;	// Edge lies in the triangle plane
			} else {
				return 0;	// Edge does not lie in the triangle plane
			}
		}

		// If the edge is not parallel, compute the intersection point
		const double r = a / b;
		if (r < 0.0){
			return 0;	// Edge goes away from triangle
		} else if (r > 1.0) {
			return 0;	// Edge goes away from triangle
		}

		rIntersectionPoint = rLinePoint1 + r*dir;

		// Check if the intersection point is inside the triangle
		const double uu = inner_prod(u,u);
		const double uv = inner_prod(u,v);
		const double vv = inner_prod(v,v);
		const array_1d<double,3> w = rIntersectionPoint - rTriangleGeometry[0];
		const double wu = inner_prod(w,u);
		const double wv = inner_prod(w,v);
		const double D = uv * uv - uu * vv;

		// Get and test parametric coords
		const double s = (uv * wv - vv * wu) / D;
		if (s < 0.0 - epsilon || s > 1.0 + epsilon){
			return 0;	// Intersection point is outside the triangle
		}			
		const double t = (uv * wu - uu * wv) / D;
		if (t < 0.0 - epsilon || (s + t) > 1.0 + epsilon){
			return 0;	// Intersection point is outside the triangle
		}  	

		return 1;	// Intersection point is inside the triangle
	}

	/**
     * Find the 2D intersection of two lines
     * @param rLineGeometry Is the line to intersect
     * @param rLinePoint1 Coordinates of the first point of the intersecting line
     * @param rLinePoint2 Coordinates of the second point of the intersecting line
     * @return rIntersectionPoint The intersection point coordinates
     * @return The intersection type index: 
     * 0 (disjoint - no intersection) 
     * 1 (intersect in a unique point)
     * 2 (overlap)
     */
	template <class TGeometryType>
	static int ComputeLineLineIntersection(
        const TGeometryType& rLineGeometry, 
        const array_1d<double,3>& rLinePoint0, 
        const array_1d<double,3>& rLinePoint1, 
        array_1d<double,3>& rIntersectionPoint,
		const double epsilon = 1e-12) {

		const array_1d<double,3> r = rLineGeometry[1] - rLineGeometry[0];
		const array_1d<double,3> s = rLinePoint1 - rLinePoint0;
		const array_1d<double,3> p_q = rLineGeometry[0] - rLinePoint0;		// p - q 
		const array_1d<double,3> q_p = rLinePoint0 - rLineGeometry[0];		// q - p

		const double aux_1 = CrossProd2D(r,s);
		const double aux_2 = CrossProd2D(q_p,r);
		const double aux_3 = CrossProd2D(q_p,s);
		
		// Check intersection cases
		// Check that the lines are not collinear
		if (std::abs(aux_1) < epsilon && std::abs(aux_2) < epsilon){
			const double aux_4 = inner_prod(r,r);
			const double aux_5 = inner_prod(s,r);
			const double t_0 = inner_prod(q_p,r)/aux_4;
			const double t_1 = t_0 + aux_5/aux_4;
			if (aux_5 < 0.0){
				if (t_1 >= 0.0 && t_0 <= 1.0){
					return 2;	// The lines are collinear and overlapping
				}
			} else {
				if (t_0 >= 0.0 && t_1 <= 1.0){
					return 2;	// The lines are collinear and overlapping
				}
			}
		// Check that the lines are parallel and non-intersecting
		} else if (std::abs(aux_1) < epsilon && std::abs(aux_2) > epsilon){
			return 0;
		// Check that the lines intersect
		} else if (std::abs(aux_1) > epsilon){
			const double u = aux_2/aux_1;
			const double t = aux_3/aux_1;
			if (((u >= 0.0) && (u <= 1.0)) && ((t >= 0.0) && (t <= 1.0))){
				rIntersectionPoint = rLinePoint0 + u*s;
				return 1;
			}
		}
		// Otherwise, the lines are non-parallel but do not intersect
		return 0;														
	}

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
private:

    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This inline function computes the 2D cross product between two arrays
     * @param a First vector
     * @param b Second vector
     * @return The 2D cross product value
     */
	static inline double CrossProd2D(const array_1d<double,3> &a, const array_1d<double,3> &b){
		return (a(0)*b(1) - a(1)*b(0));
	}

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Private LifeCycle
    ///@{

    ///@}
    ///@name Unaccessible methods
    ///@{

}; /* Class IntersectionUtilities */

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/

#endif /* KRATOS_INTERSECTION_UTILITIES  defined */

