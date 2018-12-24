//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:           BSD License
//                          Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//
//

// System includes
#include <unordered_map>
#include <utility>

// External includes

// Project includes
#include "processes/tetrahedral_mesh_orientation_check.h"
#include "utilities/math_utils.h"
#include "includes/key_hash.h"

namespace Kratos
{
KRATOS_CREATE_LOCAL_FLAG(TetrahedralMeshOrientationCheck,ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS, 0);
KRATOS_CREATE_LOCAL_FLAG(TetrahedralMeshOrientationCheck,COMPUTE_NODAL_NORMALS, 1);
KRATOS_CREATE_LOCAL_FLAG(TetrahedralMeshOrientationCheck,COMPUTE_CONDITION_NORMALS, 2);

/***********************************************************************************/
/***********************************************************************************/

void TetrahedralMeshOrientationCheck::Execute()
{
    KRATOS_TRY;

    if(mrOptions.Is(COMPUTE_NODAL_NORMALS)) {
        KRATOS_ERROR_IF_NOT(mrModelPart.NodesBegin()->SolutionStepsDataHas(NORMAL)) << "Missing NORMAL variable on solution step data" << std::endl;
        for(ModelPart::NodesContainerType::iterator it_node = mrModelPart.NodesBegin(); it_node != mrModelPart.NodesEnd(); it_node++) {
            noalias(it_node->FastGetSolutionStepValue(NORMAL)) = ZeroVector(3);
        }
    }

    //********************************************************
    //begin by orienting all of the elements in the volume
    unsigned int ElemSwitchCount = 0;

    for (ModelPart::ElementIterator it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); it_elem++) {
        ElementType::GeometryType& rGeom = it_elem->GetGeometry();
        GeometryData::KratosGeometryType GeoType = rGeom.GetGeometryType();

        if (GeoType == GeometryData::Kratos_Tetrahedra3D4  || GeoType == GeometryData::Kratos_Triangle2D3) {
            bool Switched = this->Orient(rGeom);
            if (Switched)
                ElemSwitchCount++;
        }
    }

    // Generate output message, throw error if necessary
    std::stringstream OutMsg;
    if (ElemSwitchCount > 0) {
        OutMsg << "Mesh orientation check found " << ElemSwitchCount << " inverted elements." << std::endl;
    } else {
        OutMsg << "No inverted elements found" << std::endl;
    }


    //********************************************************
    //reset the flag BOUNDARY on all of the nodes
    for(ModelPart::NodesContainerType::iterator it_node = mrModelPart.NodesBegin(); it_node != mrModelPart.NodesEnd(); it_node++)
    {
        it_node->Set(BOUNDARY, false);
    }

    //********************************************************
    //next check that the conditions are oriented accordingly

    //to do so begin by putting all of the conditions in a map

    typedef std::unordered_map<DenseVector<int>, Condition::Pointer, KeyHasherRange<DenseVector<int>>, KeyComparorRange<DenseVector<int>> > hashmap;
    hashmap faces_map;

    for (ModelPart::ConditionIterator itCond = mrModelPart.ConditionsBegin(); itCond != mrModelPart.ConditionsEnd(); itCond++)
    {
        itCond->Set(VISITED,false); //mark

        Geometry< Node<3> >& geom = itCond->GetGeometry();
        DenseVector<int> ids(geom.size());

        for(unsigned int i=0; i<ids.size(); i++)
        {
            geom[i].Set(BOUNDARY,true);
            ids[i] = geom[i].Id();
        }

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids.begin(), ids.end());

        //insert a pointer to the condition identified by the hash value ids
        //faces_map.insert( std::make_pair<DenseVector<int>, Condition::Pointer >(ids, *itCond.base()) );
        faces_map.insert( hashmap::value_type(ids, *itCond.base()) );
        //faces_map[ids] = *itCond.base();
    }

    //now loop for all the elements and for each face of the element check if it is in the "faces_map"
    //if it happens to be there check the orientation
    unsigned int CondSwitchCount = 0;
    for (ModelPart::ElementIterator it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); it_elem++)
    {
        ElementType::GeometryType& rGeom = it_elem->GetGeometry();
        GeometryData::KratosGeometryType GeoType = rGeom.GetGeometryType();

        if (GeoType == GeometryData::Kratos_Tetrahedra3D4  || GeoType == GeometryData::Kratos_Triangle2D3)
        {
            //allocate a work array long enough to contain the Ids of a face
            DenseVector<int> aux( rGeom.size() - 1);

            //loop over the faces
            for(unsigned int outer_node_index=0; outer_node_index< rGeom.size(); outer_node_index++)
            {
                unsigned int localindex_node_on_face = -1;
                //we put in "aux" the indices of all of the nodes which do not
                //coincide with the face_index we are currently considering
                //telling in other words:
                //face_index will contain the local_index of the node which is NOT on the face
                //localindex_node_on_face the local_index of one of the nodes on the face
                unsigned int counter = 0;
                for(unsigned int i=0; i<rGeom.size(); i++)
                {
                    if(i != outer_node_index)
                    {
                        aux[counter++] = rGeom[i].Id();
                        localindex_node_on_face = i;
                    }
                }

                //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
                std::sort(aux.begin(), aux.end());

                hashmap::iterator it_face = faces_map.find(aux);
                if(it_face != faces_map.end() ) //it was actually found!!
                {
                    //mark the condition as visited. This will be useful for a check at the endif
                    (it_face->second)->Set(VISITED,true);

                    if(mrOptions.Is(ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS))
                    {
                        WeakPointerVector< Element > VectorOfNeighbours;
                        VectorOfNeighbours.resize(1);
                        VectorOfNeighbours(0) = Element::WeakPointer( *it_elem.base() );
                        (it_face->second)->SetValue(NEIGHBOUR_ELEMENTS, VectorOfNeighbours);
                    }

                    //compute the normal of the face
                    array_1d<double,3> FaceNormal(3,0.0);
                    Geometry<Node<3> >& rFaceGeom = (it_face->second)->GetGeometry();

                    if ( rFaceGeom.GetGeometryType() == GeometryData::Kratos_Triangle3D3 )
                        FaceNormal3D(FaceNormal,rFaceGeom);
                    else if ( rFaceGeom.GetGeometryType()  == GeometryData::Kratos_Line2D2 )
                        FaceNormal2D(FaceNormal,rFaceGeom);



                    //do a dotproduct with the DenseVector that goes from
                    //"outer_node_index" to any of the nodes in aux;
                    array_1d<double,3> lvec = rGeom[outer_node_index]-rGeom[localindex_node_on_face];

                    double dotprod = inner_prod(lvec, FaceNormal);

                    //if dotprod > 0 then the normal to the face goes in the same half space as
                    //an edge that goes from the space to the node not on the face
                    //hence the face need to be swapped
                    if(dotprod > 0) {
                        rFaceGeom(0).swap(rFaceGeom(1));
                        FaceNormal = -FaceNormal;

                        CondSwitchCount++;
                    }

                    if(mrOptions.Is(COMPUTE_NODAL_NORMALS)) {
                        double factor = 1.0/static_cast<double>(rFaceGeom.size());
                        for(unsigned int i=0; i<rFaceGeom.size(); i++)
                            rFaceGeom[i].FastGetSolutionStepValue(NORMAL) += factor*FaceNormal;
                    }
                    if(mrOptions.Is(COMPUTE_CONDITION_NORMALS)) {
                        (it_face->second)->SetValue(NORMAL, FaceNormal );
                    }

                }

            }
        }
    }

    //check that all of the conditions belong to at least an element. Throw an error otherwise (this is particularly useful in mpi)
    for (ModelPart::ConditionIterator itCond = mrModelPart.ConditionsBegin(); itCond != mrModelPart.ConditionsEnd(); itCond++) {
        KRATOS_ERROR_IF(itCond->IsNot(VISITED)) << "Found a condition without any corresponding element. ID of condition = " << itCond->Id() << std::endl;
    }


    if (CondSwitchCount > 0) {
        OutMsg << "Mesh orientation check found " << CondSwitchCount << " inverted conditions." << std::endl;
    } else {
        OutMsg << "No inverted conditions found" << std::endl;
    }


    if (mThrowErrors && (ElemSwitchCount+CondSwitchCount) > 0) {
        KRATOS_ERROR << OutMsg.str() << std::endl;
    } else {
        KRATOS_INFO("TetrahedralMeshOrientationCheck") << OutMsg.str();
    }


    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void TetrahedralMeshOrientationCheck::SwapAll()
{
    for (ModelPart::ConditionIterator itCond = mrModelPart.ConditionsBegin(); itCond != mrModelPart.ConditionsEnd(); itCond++)
    {
        ConditionType::GeometryType& rGeom = itCond->GetGeometry();
        GeometryData::KratosGeometryType GeoType = rGeom.GetGeometryType();

        if ( GeoType == GeometryData::Kratos_Triangle3D3 )
            rGeom(0).swap(rGeom(1));
    }
}

/***********************************************************************************/
/***********************************************************************************/

void TetrahedralMeshOrientationCheck::SwapNegativeElements()
{
    for (ModelPart::ElementIterator it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); it_elem++)
    {
        if(it_elem->GetGeometry().Volume() < 0.0) {
            it_elem->GetGeometry()(0).swap(it_elem->GetGeometry()(1));
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool TetrahedralMeshOrientationCheck::Orient(Geometry< Node<3> >& rGeom)
{
    const unsigned int PointIndex = 0;
    const GeometryData::IntegrationMethod Method = GeometryData::GI_GAUSS_1;

    // Re-orient the element if needed
    double DetJ = rGeom.DeterminantOfJacobian(PointIndex,Method);
    if (DetJ < 0.0) {
        // swap two nodes to change orientation
        rGeom(0).swap(rGeom(1));
        return true;
    } else
        return false;
}

/***********************************************************************************/
/***********************************************************************************/

void TetrahedralMeshOrientationCheck::FaceNormal2D(
    array_1d<double,3>& An,
    Geometry<Node<3> >& rGeometry
    )
{
    An[0] =   rGeometry[1].Y() - rGeometry[0].Y();
    An[1] = - (rGeometry[1].X() - rGeometry[0].X());
    An[2] =    0.00;

}

/***********************************************************************************/
/***********************************************************************************/

void TetrahedralMeshOrientationCheck::FaceNormal3D(
    array_1d<double,3>& An,
    Geometry<Node<3> >& rGeometry
    )
{

    array_1d<double,3> v1,v2;
    v1[0] = rGeometry[1].X() - rGeometry[0].X();
    v1[1] = rGeometry[1].Y() - rGeometry[0].Y();
    v1[2] = rGeometry[1].Z() - rGeometry[0].Z();

    v2[0] = rGeometry[2].X() - rGeometry[0].X();
    v2[1] = rGeometry[2].Y() - rGeometry[0].Y();
    v2[2] = rGeometry[2].Z() - rGeometry[0].Z();

    MathUtils<double>::CrossProduct(An,v1,v2);
}

} // namespace Kratos
