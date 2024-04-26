// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef TRIANGLE_MESH_KRATOS_INCLUDE_H
#define TRIANGLE_MESH_KRATOS_INCLUDE_H


//// Project includes
#include "containers/triangle_mesh.hpp"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  TriangleMeshKratos
 * @author Manuel Messmer
 * @brief  ....
*/
template<class ContainerType>
class TriangleMeshKratos : public TriangleMesh
{
public:
    ///@name Type Definitions
    ///@{
    void SetContainer(ContainerType& rContainer) {
        mContainer = &rContainer;
    }

    ///@brief Get triangle vertex 1
    ///@param TriangleId
    ///@return const Vector3d&
    const Vector3d& P1(IndexType TriangleId) const {
        return mContainer.Elements()[TriangleId].GetGeometry()[0];
    }

    ///@brief Get triangle vertex 2
    ///@param TriangleId
    ///@return const Vector3d&
    const Vector3d& P2(IndexType TriangleId) const {
        return mContainer.Elements()[TriangleId].GetGeometry()[1];
    }

    ///@brief Get triangle vertex 3
    ///@param TriangleId
    ///@return const Vector3d&
    const Vector3d& P3(IndexType TriangleId) const {
        return mContainer.Elements()[TriangleId].GetGeometry()[2];
    }

private:
    ContainerType* mContainer;
    ///@}
}; // End AABB_primitive class
///@} // End QuESo classes

} // End namespace queso

#endif // AABB_PRIMITIVE_INCLUDE_H