//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//
//

// System includes
#include <unordered_map>
#include <utility>
#include <sstream>

// External includes

// Project includes
#include "processes/tetrahedral_mesh_orientation_check.h"
#include "utilities/math_utils.h"
#include "utilities/variable_utils.h"
#include "includes/key_hash.h"

namespace Kratos
{
KRATOS_CREATE_LOCAL_FLAG(TetrahedralMeshOrientationCheck, THROW_ERRORS,                            0);
KRATOS_CREATE_LOCAL_FLAG(TetrahedralMeshOrientationCheck, ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS, 1);
KRATOS_CREATE_LOCAL_FLAG(TetrahedralMeshOrientationCheck, COMPUTE_NODAL_NORMALS,                   2);
KRATOS_CREATE_LOCAL_FLAG(TetrahedralMeshOrientationCheck, COMPUTE_CONDITION_NORMALS,               3);
KRATOS_CREATE_LOCAL_FLAG(TetrahedralMeshOrientationCheck, MAKE_VOLUMES_POSITIVE,                   4);
KRATOS_CREATE_LOCAL_FLAG(TetrahedralMeshOrientationCheck, ALLOW_REPEATED_CONDITIONS,               5);
KRATOS_CREATE_LOCAL_FLAG(TetrahedralMeshOrientationCheck, DETAILED_INVERTED_ENTITIES_MESSAGE,      6);

/***********************************************************************************/
/***********************************************************************************/

TetrahedralMeshOrientationCheck::TetrahedralMeshOrientationCheck(
    ModelPart& rModelPart,
    const bool ThrowErrors,
    const Flags Options
    ):  Process(),
        mrModelPart(rModelPart),
        mOptions(Options)

{
    // Set the THROW_ERRORS flag based on the input parameter
    mOptions.Set(THROW_ERRORS, ThrowErrors);
}

/***********************************************************************************/
/***********************************************************************************/

TetrahedralMeshOrientationCheck::TetrahedralMeshOrientationCheck(
    ModelPart& rModelPart,
    const Flags Options
    ):  Process(),
        mrModelPart(rModelPart),
        mOptions(Options)
{
    // By default, do not throw errors unless specified
    mOptions.Set(THROW_ERRORS, false);
}

/***********************************************************************************/
/***********************************************************************************/

TetrahedralMeshOrientationCheck::TetrahedralMeshOrientationCheck(
    Model& rModel,
    Parameters ThisParameters
    ): Process(),
        mrModelPart(rModel.GetModelPart(ThisParameters["model_part_name"].GetString()))
{
    // Validate and assign default parameters
    ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());

    // Store the number of inverted entities to print
    mNumberOfInvertedEntitiesToPrint = ThisParameters["number_of_inverted_entities_to_print"].GetInt();

    // Create and set options based on parameters
    Flags options;

    // Check for throw_errors flag
    if (ThisParameters.Has("throw_errors")) {
        if (ThisParameters["throw_errors"].GetBool()) {
            options |= THROW_ERRORS;
        }
    }

    // Check for compute_nodal_normals flag
    if (ThisParameters.Has("compute_nodal_normals") ) {
        if (ThisParameters["compute_nodal_normals"].GetBool()) {
            options |= COMPUTE_NODAL_NORMALS;
        }
    }

    // Check for compute_condition_normals flag
    if (ThisParameters.Has("compute_condition_normals")) {
        if (ThisParameters["compute_condition_normals"].GetBool()) {
            options |= COMPUTE_CONDITION_NORMALS;
        }
    }

    // Check for assign_neighbour_elements_to_conditions flag
    if (ThisParameters.Has("assign_neighbour_elements_to_conditions") ) {
        if (ThisParameters["assign_neighbour_elements_to_conditions"].GetBool()) {
            options |= ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS;
        }
    }

    // Check for allow_repeated_conditions flag
    if (ThisParameters.Has("allow_repeated_conditions")) {
        if (ThisParameters["allow_repeated_conditions"].GetBool()) {
            options |= ALLOW_REPEATED_CONDITIONS;
        }
    }

    // Check for make_volumes_positive flag
    if (ThisParameters.Has("make_volumes_positive")) {
        if (ThisParameters["make_volumes_positive"].GetBool()) {
            options |= MAKE_VOLUMES_POSITIVE;
        }
    }

    // Check for detailed_inverted_entities_message flag
    if (ThisParameters.Has("detailed_inverted_entities_message")) {
        if (ThisParameters["detailed_inverted_entities_message"].GetBool()) {
            options |= DETAILED_INVERTED_ENTITIES_MESSAGE;
        }
    }

    // Set the options
    mOptions = options;
}

/***********************************************************************************/
/***********************************************************************************/

void TetrahedralMeshOrientationCheck::Execute()
{
    KRATOS_TRY;

    if(mOptions.Is(COMPUTE_NODAL_NORMALS)) {
        KRATOS_ERROR_IF_NOT(mrModelPart.NodesBegin()->SolutionStepsDataHas(NORMAL)) << "Missing NORMAL variable on solution step data" << std::endl;
        VariableUtils().SetVariable(NORMAL, ZeroVector(3), mrModelPart.Nodes());
    }

    // Create a stringstream object for the error message
    std::stringstream ss_elements, ss_conditions;

    // Begin by orienting all of the elements in the volume
    SizeType elem_switch_count = 0;

    for (auto it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); it_elem++) {
        GeometryType& r_geometry = it_elem->GetGeometry();
        if (SupportedElement(r_geometry)) { // Acceptable quadratic geoms
            const bool switched = this->Orient(r_geometry);
            if (switched) {
                ++elem_switch_count;
                if (mNumberOfInvertedEntitiesToPrint != -1 && static_cast<int>(elem_switch_count) <= mNumberOfInvertedEntitiesToPrint) {
                    // Append the element ID to the stringstream
                    if (mOptions.Is(DETAILED_INVERTED_ENTITIES_MESSAGE)) {
                        ss_elements << "Element number: " << it_elem->Id() << " is inverted. Volume: " << r_geometry.DomainSize() << "\n";
                    } else {
                        ss_elements << "\t" << it_elem->Id();
                        if (static_cast<int>(elem_switch_count) == mNumberOfInvertedEntitiesToPrint) {
                            ss_elements << " ...";
                        } else {
                            ss_elements << ",";
                        }
                    }
                }
            }
        }
    }

    // Generate output message, throw error if necessary
    std::stringstream out_message;
    if (elem_switch_count > 0) {
        out_message << "Mesh orientation check found " << elem_switch_count << " inverted elements.\n";
        if (mOptions.Is(DETAILED_INVERTED_ENTITIES_MESSAGE)) {
            out_message << "The following elements were inverted:\n"
                        << ss_elements.str();
        } else {
            out_message << "The IDs of the inverted elements are:\n"
                        << ss_elements.str() << "\n";
        }
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

        GeometryType& r_geometry = it_cond->GetGeometry();

        if (SupportedCondition(r_geometry)) { 
            DenseVector<int> ids(r_geometry.size());

            for(IndexType i=0; i<ids.size(); i++) {
                r_geometry[i].Set(BOUNDARY,true);
                ids[i] = r_geometry[i].Id();
            }

            //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
            std::sort(ids.begin(), ids.end());

            // Insert a pointer to the condition identified by the hash value ids
            hashmap::iterator it_face = faces_map.find(ids);
            if (it_face != faces_map.end() ) { // Already defined geometry
                KRATOS_ERROR_IF_NOT(mOptions.Is(ALLOW_REPEATED_CONDITIONS)) << "The condition of ID:\t" << it_cond->Id() << " shares the same geometry as the condition ID:\t" << it_face->second[0]->Id() << " this is not allowed. Please, check your mesh" << std::endl;
                it_face->second.push_back(*it_cond.base());
            } else {
                faces_map.insert( hashmap::value_type(ids, std::vector<Condition::Pointer>({*it_cond.base()})) );
            }
        }
    }

    // Now loop for all the elements and for each face of the element check if it is in the "faces_map"
    // if it happens to be there check the orientation
    SizeType cond_switch_count = 0;
    // Allocate a work array long enough to contain the Ids of a face
    DenseVector<int> aux;
    DenseMatrix<int> boundaries_nodes;
    for (auto it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); it_elem++) {
        GeometryType& r_geometry = it_elem->GetGeometry();
        if (SupportedElement(r_geometry)){ //acceptable quadratic geoms
            IndexType n_boundaries = BoundariesEntitiesNumber(r_geometry);
            aux.resize(NumberOfNodesInEachBoundary(r_geometry));
            NodesOfBoundaries(r_geometry, boundaries_nodes); //getting the nodes of each of the n_boundaries of the element
            const array_1d<double,3> elem_center = r_geometry.Center();
            const SizeType first_node_id_of_elem = r_geometry[0].Id();
            // Loop over the faces(edges in 2d geoms)
            for (IndexType i_face=0; i_face< n_boundaries; i_face++) {
                aux = row(boundaries_nodes, i_face);

                //finding local index of node that is part of the face
                IndexType localindex_node_on_face = 1; //default option. but if first node is contained, then we assign it as 0:
                for (IndexType face_node_id : aux ) {
                    if(face_node_id==first_node_id_of_elem){
                        localindex_node_on_face=0;
                    }
                }

                //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
                std::sort(aux.begin(), aux.end());

                hashmap::iterator it_face = faces_map.find(aux);
                if (it_face != faces_map.end()) { // It was actually found!!
                    // Mark the condition as visited. This will be useful for a check at the endif
                    std::vector<Condition::Pointer>& list_conditions = it_face->second;
                    for (const Condition::Pointer& p_cond : list_conditions) {
                        p_cond->Set(VISITED,true);
                    }

                    if (mOptions.Is(ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS)) {
                        GlobalPointersVector< Element > vector_of_neighbours;
                        vector_of_neighbours.resize(1);
                        vector_of_neighbours(0) = Element::WeakPointer( *it_elem.base() );
                        for (const Condition::Pointer& p_cond : list_conditions) {
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
                    array_1d<double,3> lvec = elem_center-r_geometry[localindex_node_on_face];

                    const double dotprod = inner_prod(lvec, face_normal);

                    // If dotprod > 0 then the normal to the face goes in the same half space as
                    // an edge that goes from the space to the node not on the face hence the face need to be swapped
                    if (dotprod > 0.0) {
                        r_face_geom(0).swap(r_face_geom(1));
                        face_normal = -face_normal;
                        ++cond_switch_count;
                        if (mNumberOfInvertedEntitiesToPrint != -1 && static_cast<int>(cond_switch_count) <= mNumberOfInvertedEntitiesToPrint) {
                            if (mOptions.Is(DETAILED_INVERTED_ENTITIES_MESSAGE)) {
                                ss_conditions << "Condition number: " << (list_conditions[0])->Id() << " is inverted. Area: " << r_face_geom.DomainSize() << "\n";
                            } else {
                                ss_conditions << "\t" << (list_conditions[0])->Id();
                                if (static_cast<int>(cond_switch_count) == mNumberOfInvertedEntitiesToPrint) {
                                    ss_conditions << " ...";
                                } else {
                                    ss_conditions << ",";
                                }
                            }
                        }
                    }

                    if (mOptions.Is(COMPUTE_NODAL_NORMALS)) {
                        for(IndexType i=0; i<r_face_geom.size(); i++) {
                            r_face_geom.PointLocalCoordinates(local_coords, r_face_geom[i].Coordinates());
                            r_face_geom[i].FastGetSolutionStepValue(NORMAL) += r_face_geom.Normal(local_coords);
                        }
                    }
                    if (mOptions.Is(COMPUTE_CONDITION_NORMALS)) {
                        for (const Condition::Pointer& p_cond : list_conditions) {
                            p_cond->SetValue(NORMAL, face_normal);
                        }
                    }

                }

            }
        }
    }

    // Check that all of the conditions belong to at least an element. Throw an error otherwise (this is particularly useful in mpi)
    for (auto& r_cond : mrModelPart.Conditions()) {
        const GeometryType& r_geometry = r_cond.GetGeometry();
        if(SupportedCondition(r_geometry)){ //acceptable quadratic geoms
            KRATOS_ERROR_IF(r_cond.IsNot(VISITED)) << "Found a condition without any corresponding element. ID of condition = " << r_cond.Id() << std::endl;
        }
    }

    if (cond_switch_count > 0) {
        out_message << "Mesh orientation check found " << cond_switch_count << " inverted conditions.\n";
        if (mOptions.Is(DETAILED_INVERTED_ENTITIES_MESSAGE)) {
            out_message << "The following conditions were inverted:\n"
                        << ss_conditions.str();
        } else {
            out_message << "The IDs of the inverted conditions are:\n"
                        << ss_conditions.str() << "\n";
        }
    } else {
        out_message << "No inverted conditions found" << std::endl;
    }

    if (mOptions.Is(THROW_ERRORS) && (elem_switch_count+cond_switch_count) > 0) {
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

        if (geometry_type == GeometryData::KratosGeometryType::Kratos_Triangle3D3) {
            r_geometry(0).swap(r_geometry(1));
        }
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

const Parameters TetrahedralMeshOrientationCheck::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"({
        "model_part_name"                         : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "throw_errors"                            : false,
        "compute_nodal_normals"                   : false,
        "compute_condition_normals"               : false,
        "assign_neighbour_elements_to_conditions" : false,
        "allow_repeated_conditions"               : false,
        "make_volumes_positive"                   : false,
        "detailed_inverted_entities_message"      : false,
        "number_of_inverted_entities_to_print"    : 20 // -1 to unlimited
    })");

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

bool TetrahedralMeshOrientationCheck::Orient(GeometryType& rGeometry)
{
    const IndexType point_index = 0;
    const GeometryData::IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_1;

    // Re-orient the element if needed
    const double det_J = rGeometry.DeterminantOfJacobian(point_index,integration_method);
    if (det_J < 0.0) {
        // Swap two nodes to change orientation
        rGeometry(0).swap(rGeometry(1));
        return true;
    } else {
        return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool TetrahedralMeshOrientationCheck::LinearElement(const GeometryType& rGeometry)
{
    const GeometryData::KratosGeometryType geometry_type = rGeometry.GetGeometryType();
    return (geometry_type == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4 || geometry_type == GeometryData::KratosGeometryType::Kratos_Triangle2D3);
}

/***********************************************************************************/
/***********************************************************************************/

bool TetrahedralMeshOrientationCheck::SupportedElement(const GeometryType& rGeometry)
{
    const GeometryData::KratosGeometryType geometry_type = rGeometry.GetGeometryType();
    return (geometry_type == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4  || geometry_type == GeometryData::KratosGeometryType::Kratos_Triangle2D3 || //acceptable linear geoms
            geometry_type == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10 || geometry_type == GeometryData::KratosGeometryType::Kratos_Triangle2D6); //acceptable quadratic geoms
}

/***********************************************************************************/
/***********************************************************************************/

bool TetrahedralMeshOrientationCheck::SupportedCondition(const GeometryType& rGeometry)
{
    const GeometryData::KratosGeometryType geometry_type = rGeometry.GetGeometryType();
    return (geometry_type == GeometryData::KratosGeometryType::Kratos_Triangle3D3 || geometry_type == GeometryData::KratosGeometryType::Kratos_Line2D2 || //acceptable linear geoms
            geometry_type == GeometryData::KratosGeometryType::Kratos_Triangle3D6 || geometry_type == GeometryData::KratosGeometryType::Kratos_Line2D3); //acceptable quadratic geoms
}

/***********************************************************************************/
/***********************************************************************************/

TetrahedralMeshOrientationCheck::SizeType TetrahedralMeshOrientationCheck::BoundariesEntitiesNumber(const GeometryType& rGeometry) 
{
    const SizeType dimension = rGeometry.LocalSpaceDimension();
    if (dimension == 3) {
        return rGeometry.FacesNumber();
    } else if (dimension == 2) {
        return rGeometry.EdgesNumber();
    }
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

TetrahedralMeshOrientationCheck::SizeType TetrahedralMeshOrientationCheck::NumberOfNodesInEachBoundary(const GeometryType& rGeometry) 
{
    const GeometryData::KratosGeometryType geometry_type = rGeometry.GetGeometryType();
    if  (geometry_type == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4){
        return 3;
    } else if (geometry_type == GeometryData::KratosGeometryType::Kratos_Triangle2D3){
        return 2;
    } else if (geometry_type == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10){
        return 6;
    } else if (geometry_type == GeometryData::KratosGeometryType::Kratos_Triangle2D6){
        return 3;
    }
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

void TetrahedralMeshOrientationCheck::NodesOfBoundaries(const GeometryType& rGeometry, DenseMatrix<int>& rNodesIds){
    IndexType n_boundaries = BoundariesEntitiesNumber(rGeometry);
    IndexType nodes_in_boundary = NumberOfNodesInEachBoundary(rGeometry);
    rNodesIds.resize(n_boundaries,nodes_in_boundary);
    if(LinearElement(rGeometry)){
        for(IndexType outer_node_index=0; outer_node_index< BoundariesEntitiesNumber(rGeometry); outer_node_index++) {
            // We put in rNodesIds the indices of all of the nodes which do not
            // coincide with the face_index we are currently considering telling in other words:
            // face_index will contain the local_index of the node which is NOT on the face
            IndexType counter = 0;
            for(IndexType i=0; i<rGeometry.size(); i++) {
                if(i != outer_node_index) {
                    rNodesIds(outer_node_index,counter++) = rGeometry[i].Id();
                }
            }
        }
    } else { //quadratic element
        const auto& faces = rGeometry.GenerateFaces();
        for(IndexType outer_node_index=0; outer_node_index< faces.size(); outer_node_index++) {
            for(IndexType i=0; i<faces[outer_node_index].size(); i++) {
                rNodesIds(outer_node_index,i) = faces[outer_node_index][i].Id();  
            }
        }
    }
}

} // namespace Kratos