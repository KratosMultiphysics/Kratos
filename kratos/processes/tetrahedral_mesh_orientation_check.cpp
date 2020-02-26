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
KRATOS_CREATE_LOCAL_FLAG(TetrahedralMeshOrientationCheck,ALLOW_REPEATED_CONDITIONS, 4);

/***********************************************************************************/
/***********************************************************************************/

void TetrahedralMeshOrientationCheck::Execute()
{
    KRATOS_TRY;

    if(mrOptions.Is(COMPUTE_NODAL_NORMALS)) {
        KRATOS_ERROR_IF_NOT(mrModelPart.NodesBegin()->SolutionStepsDataHas(NORMAL)) << "Missing NORMAL variable on solution step data" << std::endl;
        VariableUtils().SetVariable(NORMAL, ZeroVector(3), mrModelPart.Nodes());
    }

    // Begin by orienting all of the elements in the volume
    SizeType elem_switch_count = 0;

    for (auto it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); it_elem++) {
        GeometryType& r_geometry = it_elem->GetGeometry();
        const GeometryData::KratosGeometryType geometry_type = r_geometry.GetGeometryType();

        if (geometry_type == GeometryData::Kratos_Tetrahedra3D4  || geometry_type == GeometryData::Kratos_Triangle2D3) {
            const bool switched = this->Orient(r_geometry);
            if (switched)
                elem_switch_count++;
        }
    }

    // Generate output message, throw error if necessary
    std::stringstream out_message;
    if (elem_switch_count > 0) {
        out_message << "Mesh orientation check found " << elem_switch_count << " inverted elements." << std::endl;
    } else {
        out_message << "No inverted elements found" << std::endl;
    }

    // Reset the flag BOUNDARY on all of the nodes
    VariableUtils().SetFlag(BOUNDARY, false, mrModelPart.Nodes());

    // Next check that the conditions are oriented accordingly
    // to do so begin by putting all of the conditions in a map
    typedef std::unordered_map<DenseVector<int>, std::vector<Condition::Pointer>, KeyHasherRange<DenseVector<int>>, KeyComparorRange<DenseVector<int>> > hashmap;
    hashmap faces_map;

    for (auto it_cond = mrModelPart.ConditionsBegin(); it_cond != mrModelPart.ConditionsEnd(); it_cond++) {
        it_cond->Set(VISITED, false); //mark

        GeometryType& geom = it_cond->GetGeometry();
        DenseVector<int> ids(geom.size());

        for(IndexType i=0; i<ids.size(); i++) {
            geom[i].Set(BOUNDARY,true);
            ids[i] = geom[i].Id();
        }

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids.begin(), ids.end());

        // Insert a pointer to the condition identified by the hash value ids
        hashmap::iterator it_face = faces_map.find(ids);
        if(it_face != faces_map.end() ) { // Already defined geometry
            KRATOS_ERROR_IF_NOT(mrOptions.Is(ALLOW_REPEATED_CONDITIONS)) << "The condition of ID:\t" << it_cond->Id() << " shares the same geometry as the condition ID:\t" << it_face->second[0]->Id() << " this is not allowed. Please, check your mesh" << std::endl;
            it_face->second.push_back(*it_cond.base());
        } else {
            faces_map.insert( hashmap::value_type(ids, std::vector<Condition::Pointer>({*it_cond.base()})) );
        }
    }

    // Now loop for all the elements and for each face of the element check if it is in the "faces_map"
    // if it happens to be there check the orientation
    SizeType cond_switch_count = 0;
    for (auto it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); it_elem++) {
        GeometryType& r_geometry = it_elem->GetGeometry();
        const GeometryData::KratosGeometryType geometry_type = r_geometry.GetGeometryType();

        if (geometry_type == GeometryData::Kratos_Tetrahedra3D4  || geometry_type == GeometryData::Kratos_Triangle2D3) {
            // Allocate a work array long enough to contain the Ids of a face
            DenseVector<int> aux( r_geometry.size() - 1);

            // Loop over the faces
            for(IndexType outer_node_index=0; outer_node_index< r_geometry.size(); outer_node_index++) {
                IndexType localindex_node_on_face = -1;
                // We put in "aux" the indices of all of the nodes which do not
                // coincide with the face_index we are currently considering telling in other words:
                // face_index will contain the local_index of the node which is NOT on the face
                // localindex_node_on_face the local_index of one of the nodes on the face
                IndexType counter = 0;
                for(IndexType i=0; i<r_geometry.size(); i++) {
                    if(i != outer_node_index) {
                        aux[counter++] = r_geometry[i].Id();
                        localindex_node_on_face = i;
                    }
                }

                //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
                std::sort(aux.begin(), aux.end());

                hashmap::iterator it_face = faces_map.find(aux);
                if(it_face != faces_map.end()) { // It was actually found!!
                    // Mark the condition as visited. This will be useful for a check at the endif
                    std::vector<Condition::Pointer>& list_conditions = it_face->second;
                    for (Condition::Pointer p_cond : list_conditions) {
                        p_cond->Set(VISITED,true);
                    }

                    if(mrOptions.Is(ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS)) {
                        GlobalPointersVector< Element > vector_of_neighbours;
                        vector_of_neighbours.resize(1);
                        vector_of_neighbours(0) = Element::WeakPointer( *it_elem.base() );
                        for (Condition::Pointer p_cond : list_conditions) {
                            p_cond->SetValue(NEIGHBOUR_ELEMENTS, vector_of_neighbours);
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
                    array_1d<double,3> lvec = r_geometry[outer_node_index]-r_geometry[localindex_node_on_face];

                    const double dotprod = inner_prod(lvec, face_normal);

                    // If dotprod > 0 then the normal to the face goes in the same half space as
                    // an edge that goes from the space to the node not on the face hence the face need to be swapped
                    if(dotprod > 0.0) {
                        r_face_geom(0).swap(r_face_geom(1));
                        face_normal = -face_normal;

                        cond_switch_count++;
                    }

                    if(mrOptions.Is(COMPUTE_NODAL_NORMALS)) {
                        for(IndexType i=0; i<r_face_geom.size(); i++) {
                            r_face_geom.PointLocalCoordinates(local_coords, r_face_geom[i].Coordinates());
                            r_face_geom[i].FastGetSolutionStepValue(NORMAL) += r_face_geom.Normal(local_coords);
                        }
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

    if (cond_switch_count > 0) {
        out_message << "Mesh orientation check found " << cond_switch_count << " inverted conditions." << std::endl;
    } else {
        out_message << "No inverted conditions found" << std::endl;
    }

    if (mThrowErrors && (elem_switch_count+cond_switch_count) > 0) {
        mrModelPart.GetProcessInfo().SetValue(FLAG_VARIABLE, 0.0); //Set flag variable as check, this is not supposed to reach here anyway
        KRATOS_ERROR << out_message.str() << std::endl;
    } else {
        KRATOS_INFO("TetrahedralMeshOrientationCheck") << out_message.str();
        mrModelPart.GetProcessInfo().SetValue(FLAG_VARIABLE, 1.0); //Set flag variable as check
    }


    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void TetrahedralMeshOrientationCheck::SwapAll()
{
    for (auto& r_cond : mrModelPart.Conditions()) {
        GeometryType& r_geometry = r_cond.GetGeometry();
        const GeometryData::KratosGeometryType geometry_type = r_geometry.GetGeometryType();

        if ( geometry_type == GeometryData::Kratos_Triangle3D3 )
            r_geometry(0).swap(r_geometry(1));
    }
}

/***********************************************************************************/
/***********************************************************************************/

void TetrahedralMeshOrientationCheck::SwapNegativeElements()
{
    for (auto& r_elem : mrModelPart.Elements()) {
        if(r_elem.GetGeometry().Volume() < 0.0) {
            r_elem.GetGeometry()(0).swap(r_elem.GetGeometry()(1));
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool TetrahedralMeshOrientationCheck::Orient(GeometryType& rGeometry)
{
    const IndexType point_index = 0;
    const GeometryData::IntegrationMethod integration_method = GeometryData::GI_GAUSS_1;

    // Re-orient the element if needed
    const double DetJ = rGeometry.DeterminantOfJacobian(point_index,integration_method);
    if (DetJ < 0.0) {
        // Swap two nodes to change orientation
        rGeometry(0).swap(rGeometry(1));
        return true;
    } else {
        return false;
    }
}

} // namespace Kratos
