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
    TriangleMeshKratos(const ContainerType& rContainer) : mContainer(rContainer)
    {
        auto& r_normals = this->GetNormals();
        r_normals.reserve(rContainer.NumberOfElements());
        for( auto& r_el : rContainer.Elements() ){
            r_normals.push_back( TriangleMesh::Normal( r_el.GetGeometry()[0].Coordinates().data(),
                                                       r_el.GetGeometry()[1].Coordinates().data(),
                                                       r_el.GetGeometry()[2].Coordinates().data() ));
        }
    }

    ///@brief Get triangle vertex 1
    ///@param TriangleId
    ///@return const Vector3d&
    const Vector3d P1(IndexType TriangleId) const override {
        return (mContainer.ElementsBegin() + TriangleId)->GetGeometry()[0].Coordinates().data();
    }

    ///@brief Get triangle vertex 2
    ///@param TriangleId
    ///@return const Vector3d&
    const Vector3d P2(IndexType TriangleId) const override{
        return (mContainer.ElementsBegin() + TriangleId)->GetGeometry()[1].Coordinates().data();
    }

    ///@brief Get triangle vertex 3
    ///@param TriangleId
    ///@return const Vector3d&
    const Vector3d P3(IndexType TriangleId) const override{
        return (mContainer.ElementsBegin() + TriangleId)->GetGeometry()[2].Coordinates().data();
    }

    ///@brief Get number of triangles in mesh.
    IndexType NumOfTriangles() const override{
        return mContainer.NumberOfElements();
    }

    ///@brief Get number of vertices in mesh.
    IndexType NumOfVertices() const override{
        return mContainer.NumberOfNodes();
    }

private:

    const ContainerType& mContainer;
    ///@}
}; // End AABB_primitive class
///@} // End QuESo classes

} // End namespace queso

#endif // AABB_PRIMITIVE_INCLUDE_H