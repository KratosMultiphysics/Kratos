// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef AABB_PRIMITIVE_INCLUDE_H
#define AABB_PRIMITIVE_INCLUDE_H

//// External includes
#include "aabb_tree/AABB_base.h"
//// Project includes
#include "embedding/aabb_primitive_base.h"

namespace queso {

// Alias to external base class.
typedef aabb_base::AABB_base AABB_lohedges;

///@name QuESo Classes
///@{

/**
 * @class  AABB primitive
 * @author Manuel Messmer
 * @brief  AABB primitive to be used in AABB_tree. Provides functions to check for aabb-triangle
 *         and aabb-aabb intersections.
*/
class AABB_primitive : public AABB_primitive_base, public AABB_lohedges
{
public:
    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Life cycle
    ///@{
    AABB_primitive(const Vector3d& rLowerBound, const Vector3d& rUpperBound) :
        AABB_lohedges(rLowerBound, rUpperBound)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    ///@brief Returns true if AABB intersects with aabb.
    ///@param aabb const AABB_primitive &aabb.
    ///@return bool.
    bool intersect(const AABB_primitive &aabb) const override;

    ///@brief Returns true if AABB intersect with triangle.
    ///@brief This function uses the seperating axis theorem. https://dyn4j.org/2010/01/sat/
    ///@param v0 Vertex 1 of triangle.
    ///@param v1 Vertex 2 of triangle.
    ///@param v2 Vertex 3 of triangle.
    ///@param tolerance. Reduces extent of aabb. If tolerance > 0 -> touching triangles are not detected as intersected.
    bool intersect(const Vector3d &v0, const Vector3d &v1, const Vector3d &v2,
                   double tolerance) const;

private:

    ///@}
    ///@name Private Operations
    ///@{

    ///@brief Returns true if axes intersect.
    ///@details This function uses the seperating axis theorem.
    ///@param u0 Normal vector of aabb (1,0,0).
    ///@param u1 Normal vector of aabb (0,1,0).
    ///@param u2 Normal vector of aabb (0,0,1).
    ///@param v0 Vertex triangle 1. Triangle must be shifted to origin.
    ///@param v0 Vertex triangle 2. Triangle must be shifted to origin.
    ///@param v0 Vertex triangle 3. Triangle must be shifted to origin.
    ///@param extent Extent of aabb. Half of edge length in each direction.
    ///@param test_axis Axis to be tested.
    ///@return bool
    bool check_axis( const Vector3d &u0, const Vector3d &u1, const Vector3d &u2,
                     const Vector3d &v0, const Vector3d &v1, const Vector3d &v2,
                     const Vector3d &extent, const Vector3d& test_axis ) const;

    ///@}
}; // End AABB_primitive class
///@} // End QuESo classes

} // End namespace queso

#endif // AABB_PRIMITIVE_INCLUDE_H