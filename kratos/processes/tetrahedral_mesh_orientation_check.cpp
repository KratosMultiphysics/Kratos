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
#include "utilities/variable_utils.h"
#include "includes/key_hash.h"

namespace Kratos
{
KRATOS_CREATE_LOCAL_FLAG(TetrahedralMeshOrientationCheck,ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS, 0);
KRATOS_CREATE_LOCAL_FLAG(TetrahedralMeshOrientationCheck,COMPUTE_NODAL_NORMALS, 1);
KRATOS_CREATE_LOCAL_FLAG(TetrahedralMeshOrientationCheck,COMPUTE_CONDITION_NORMALS, 2);
KRATOS_CREATE_LOCAL_FLAG(TetrahedralMeshOrientationCheck,MAKE_VOLUMES_POSITIVE, 3);
KRATOS_CREATE_LOCAL_FLAG(TetrahedralMeshOrientationCheck,ALLOW_CONDITIONS_WITH_SAME_GEOMETRY, 4);

/***********************************************************************************/
/***********************************************************************************/

void TetrahedralMeshOrientationCheck::Execute()
{
    KRATOS_TRY;

    if(mrOptions.Is(COMPUTE_NODAL_NORMALS)) {
        KRATOS_ERROR_IF_NOT(mrModelPart.NodesBegin()->SolutionStepsDataHas(NORMAL)) << "Missing NORMAL variable on solution step data" << std::endl;
        VariableUtils().SetVectorVar(NORMAL, ZeroVector(3), mrModelPart.Nodes());
    }

    // Begin by orienting all of the elements in the volume
    unsigned int ElemSwitchCount = 0;

    for (ModelPart::ElementIterator it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); it_elem++) {
        GeometryType& rGeom = it_elem->GetGeometry();
        GeometryData::KratosGeometryType GeoType = rGeom.GetGeometryType();

        if (GeoType == GeometryData::Kratos_Tetrahedra3D4  || GeoType == GeometryData::Kratos_Triangle2D3) {
            const bool Switched = this->Orient(rGeom);
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

    // Reset the flag BOUNDARY on all of the nodes
    VariableUtils().SetFlag(BOUNDARY, false, mrModelPart.Nodes());

    // Next check that the conditions are oriented accordingly
    // to do so begin by putting all of the conditions in a map
    typedef std::unordered_map<DenseVector<int>, std::vector<Condition::Pointer>, KeyHasherRange<DenseVector<int>>, KeyComparorRange<DenseVector<int>> > hashmap;
    hashmap faces_map;

    for (ModelPart::ConditionIterator itCond = mrModelPart.ConditionsBegin(); itCond != mrModelPart.ConditionsEnd(); itCond++)
    {
        itCond->Set(VISITED, false); //mark

        GeometryType& geom = itCond->GetGeometry();
        DenseVector<int> ids(geom.size());

        for(unsigned int i=0; i<ids.size(); i++)
        {
            geom[i].Set(BOUNDARY,true);
            ids[i] = geom[i].Id();
        }

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids.begin(), ids.end());

        // Insert a pointer to the condition identified by the hash value ids
        hashmap::iterator it_face = faces_map.find(ids);
        if(it_face != faces_map.end() ) { // Already defined geometry
            KRATOS_ERROR_IF_NOT(mrOptions.Is(ALLOW_CONDITIONS_WITH_SAME_GEOMETRY)) << "The condition of ID:\t" << itCond->Id() << " shares the same geometry as the condition ID:\t" << it_face->second[0]->Id() << " this is not allowed. Please, check your mesh" << std::endl;
            it_face->second.push_back(*itCond.base());
        } else {
            faces_map.insert( hashmap::value_type(ids, std::vector<Condition::Pointer>({*itCond.base()})) );
        }
    }

    // Now loop for all the elements and for each face of the element check if it is in the "faces_map"
    // if it happens to be there check the orientation
    unsigned int CondSwitchCount = 0;
    for (ModelPart::ElementIterator it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); it_elem++)
    {
        GeometryType& rGeom = it_elem->GetGeometry();
        GeometryData::KratosGeometryType GeoType = rGeom.GetGeometryType();

        if (GeoType == GeometryData::Kratos_Tetrahedra3D4  || GeoType == GeometryData::Kratos_Triangle2D3)
        {
            // Allocate a work array long enough to contain the Ids of a face
            DenseVector<int> aux( rGeom.size() - 1);

            // Loop over the faces
            for(unsigned int outer_node_index=0; outer_node_index< rGeom.size(); outer_node_index++)
            {
                unsigned int localindex_node_on_face = -1;
                // We put in "aux" the indices of all of the nodes which do not
                // coincide with the face_index we are currently considering telling in other words:
                // face_index will contain the local_index of the node which is NOT on the face
                // localindex_node_on_face the local_index of one of the nodes on the face
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
                if(it_face != faces_map.end() ) // It was actually found!!
                {
                    // Mark the condition as visited. This will be useful for a check at the endif
                    std::vector<Condition::Pointer>& list_conditions = it_face->second;
                    for (Condition::Pointer p_cond : list_conditions) {
                        p_cond->Set(VISITED,true);
                    }

                    if(mrOptions.Is(ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS))
                    {
                        WeakPointerVector< Element > VectorOfNeighbours;
                        VectorOfNeighbours.resize(1);
                        VectorOfNeighbours(0) = Element::WeakPointer( *it_elem.base() );
                        for (Condition::Pointer p_cond : list_conditions) {
                            p_cond->SetValue(NEIGHBOUR_ELEMENTS, VectorOfNeighbours);
                        }
                    }

                    // Compute the normal of the face
                    array_1d<double,3> face_normal = ZeroVector(3);
                    GeometryType& r_face_geom = (list_conditions[0])->GetGeometry(); // The geometry is shared, we just take the first one

                    Point::CoordinatesArrayType local_coords;
                    local_coords.clear();
                    noalias(face_normal) = r_face_geom.Normal(local_coords);

                    // Do a dotproduct with the DenseVector that goes from
                    // "outer_node_index" to any of the nodes in aux;
                    array_1d<double,3> lvec = rGeom[outer_node_index]-rGeom[localindex_node_on_face];

                    const double dotprod = inner_prod(lvec, face_normal);

                    // If dotprod > 0 then the normal to the face goes in the same half space as
                    // an edge that goes from the space to the node not on the face hence the face need to be swapped
                    if(dotprod > 0) {
                        r_face_geom(0).swap(r_face_geom(1));
                        face_normal = -face_normal;

                        CondSwitchCount++;
                    }

                    if(mrOptions.Is(COMPUTE_NODAL_NORMALS)) {
                        double factor = 1.0/static_cast<double>(r_face_geom.size());
                        for(unsigned int i=0; i<r_face_geom.size(); i++)
                            r_face_geom[i].FastGetSolutionStepValue(NORMAL) += factor*face_normal;
                    }
                    if(mrOptions.Is(COMPUTE_CONDITION_NORMALS)) {
                        for (Condition::Pointer p_cond : list_conditions) {
                            p_cond->SetValue(NORMAL, face_normal );
                        }
                    }

                }

            }
        }
    }

    //check that all of the conditions belong to at least an element. Throw an error otherwise (this is particularly useful in mpi)
    for (auto& r_cond : mrModelPart.Conditions()) {
        KRATOS_ERROR_IF(r_cond.IsNot(VISITED)) << "Found a condition without any corresponding element. ID of condition = " << r_cond.Id() << std::endl;
    }


    if (CondSwitchCount > 0) {
        OutMsg << "Mesh orientation check found " << CondSwitchCount << " inverted conditions." << std::endl;
    } else {
        OutMsg << "No inverted conditions found" << std::endl;
    }


    if (mThrowErrors && (ElemSwitchCount+CondSwitchCount) > 0) {
        mrModelPart.GetProcessInfo().SetValue(FLAG_VARIABLE, 0.0); //Set flag variable as check, this is not supposed to reach here anyway
        KRATOS_ERROR << OutMsg.str() << std::endl;
    } else {
        KRATOS_INFO("TetrahedralMeshOrientationCheck") << OutMsg.str();
        mrModelPart.GetProcessInfo().SetValue(FLAG_VARIABLE, 1.0); //Set flag variable as check
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

bool TetrahedralMeshOrientationCheck::Orient(GeometryType& rGeom)
{
    const unsigned int PointIndex = 0;
    const GeometryData::IntegrationMethod Method = GeometryData::GI_GAUSS_1;

    // Re-orient the element if needed
    double DetJ = rGeom.DeterminantOfJacobian(PointIndex,Method);
    if (DetJ < 0.0) {
        // Swap two nodes to change orientation
        rGeom(0).swap(rGeom(1));
        return true;
    } else
        return false;
}

} // namespace Kratos
