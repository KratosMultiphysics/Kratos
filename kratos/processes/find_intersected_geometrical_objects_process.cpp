//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Davand
//  Collaborators:   Ruben Zorrilla Martinez
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "geometries/line_2d_2.h"
#include "processes/find_intersected_geometrical_objects_process.h"
#include "utilities/intersection_utilities.h"

namespace Kratos
{
/// Local Flags
KRATOS_CREATE_LOCAL_FLAG(FindIntersectedGeometricalObjectsProcess, INTERSECTING_CONDITIONS, 0);
KRATOS_CREATE_LOCAL_FLAG(FindIntersectedGeometricalObjectsProcess, INTERSECTING_ELEMENTS,   1);
KRATOS_CREATE_LOCAL_FLAG(FindIntersectedGeometricalObjectsProcess, INTERSECTED_CONDITIONS,  2);
KRATOS_CREATE_LOCAL_FLAG(FindIntersectedGeometricalObjectsProcess, INTERSECTED_ELEMENTS,    3);

/***********************************************************************************/
/***********************************************************************************/

FindIntersectedGeometricalObjectsProcess::FindIntersectedGeometricalObjectsProcess(
    ModelPart& rModelPartIntersected,
    ModelPart& rModelPartIntersecting,
    const Flags Options
    ) : mrModelPartIntersected(rModelPartIntersected),
        mrModelPartIntersecting(rModelPartIntersecting),
        mOptions(Options),
        mpOctree(new OctreeType())
{
}

/***********************************************************************************/
/***********************************************************************************/

FindIntersectedGeometricalObjectsProcess::FindIntersectedGeometricalObjectsProcess(
    Model& rModel,
    Parameters ThisParameters
    ) : mrModelPartIntersected(rModel.GetModelPart(ThisParameters["intersected_model_part_name"].GetString())),
        mrModelPartIntersecting(rModel.GetModelPart(ThisParameters["intersecting_model_part_name"].GetString())),
        mpOctree(new OctreeType())
{
    const Parameters default_parameters = GetDefaultParameters();
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Checking that the names of the model parts are not empty (this is supposed to be already declared)
    const std::string& r_intersected_model_part_name = ThisParameters["intersected_model_part_name"].GetString();
    const std::string& r_intersecting_model_part_name = ThisParameters["intersecting_model_part_name"].GetString();

    KRATOS_ERROR_IF(r_intersected_model_part_name == "") << "intersected_model_part_name must be defined on parameters" << std::endl;
    KRATOS_ERROR_IF(r_intersecting_model_part_name == "") << "intersecting_model_part_name must be defined on parameters" << std::endl;

    // Setting flags
    const bool intersecting_conditions = ThisParameters["intersecting_conditions"].GetBool();
    const bool intersecting_elements = ThisParameters["intersecting_elements"].GetBool();
    const bool intersected_conditions = ThisParameters["intersected_conditions"].GetBool();
    const bool intersected_elements = ThisParameters["intersected_elements"].GetBool();

    mOptions.Set(FindIntersectedGeometricalObjectsProcess::INTERSECTING_CONDITIONS, intersecting_conditions);
    mOptions.Set(FindIntersectedGeometricalObjectsProcess::INTERSECTING_ELEMENTS, intersecting_elements);
    mOptions.Set(FindIntersectedGeometricalObjectsProcess::INTERSECTED_CONDITIONS, intersected_conditions);
    mOptions.Set(FindIntersectedGeometricalObjectsProcess::INTERSECTED_ELEMENTS, intersected_elements);
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsProcess::Initialize()
{
    GenerateOctree();
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsProcess::FindIntersectedSkinObjects(std::vector<PointerVector<GeometricalObject>>& rResults)
{
    auto& r_elements_array = mrModelPartIntersected.ElementsArray();
    auto& r_conditions_array = mrModelPartIntersected.ConditionsArray();
    SizeType number_of_entities = 0;
    if (mOptions.Is(FindIntersectedGeometricalObjectsProcess::INTERSECTED_ELEMENTS)) {
        number_of_entities += r_elements_array.size();
    }
    if (mOptions.Is(FindIntersectedGeometricalObjectsProcess::INTERSECTED_CONDITIONS)) {
        number_of_entities += r_conditions_array.size();
    }
    OtreeCellVectorType leaves;

    IndexType counter = 0;
    rResults.resize(number_of_entities);
    if (mOptions.Is(FindIntersectedGeometricalObjectsProcess::INTERSECTED_ELEMENTS)) {
        for (auto& r_element : r_elements_array) {
            leaves.clear();
            mpOctree->GetIntersectedLeaves(r_element, leaves);
            FindIntersectedSkinObjects(*r_element, leaves, rResults[counter]);
            ++counter;
        }
    }
    if (mOptions.Is(FindIntersectedGeometricalObjectsProcess::INTERSECTED_CONDITIONS)) {
        for (auto& r_condition : r_conditions_array) {
            leaves.clear();
            mpOctree->GetIntersectedLeaves(r_condition, leaves);
            FindIntersectedSkinObjects(*r_condition, leaves, rResults[counter]);
            ++counter;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsProcess::FindIntersections()
{
    this->FindIntersectedSkinObjects(mIntersectedObjects);
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<PointerVector<GeometricalObject>>& FindIntersectedGeometricalObjectsProcess::GetIntersections()
{
    return mIntersectedObjects;
}

/***********************************************************************************/
/***********************************************************************************/

ModelPart& FindIntersectedGeometricalObjectsProcess::GetModelPart1()
{
    return mrModelPartIntersected;
}

/***********************************************************************************/
/***********************************************************************************/

ModelPart& FindIntersectedGeometricalObjectsProcess::GetModelPart2()
{
    return mrModelPartIntersecting;
}

/***********************************************************************************/
/***********************************************************************************/

FindIntersectedGeometricalObjectsProcess::OctreePointerType& FindIntersectedGeometricalObjectsProcess::GetOctreePointer()
{
    return mpOctree;
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsProcess::Clear()
{
    mIntersectedObjects.clear();
    OctreePointerType aux_ptr = Kratos::make_unique<OctreeType>();
    mpOctree.swap(aux_ptr);
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsProcess::Execute()
{
    // Calling initialize first (initialize Octree)
    ExecuteInitialize();

    OtreeCellVectorType leaves;

    // Iterate over elements
    if (mOptions.Is(FindIntersectedGeometricalObjectsProcess::INTERSECTED_ELEMENTS)) {
        auto& r_elements_array = mrModelPartIntersected.Elements();
        const SizeType number_of_elements = r_elements_array.size();

        const auto it_elements_begin = r_elements_array.begin();

        #pragma omp parallel for private(leaves)
        for (int i = 0; i < static_cast<int>(number_of_elements); i++) {
            auto it_elem = it_elements_begin + i;
            leaves.clear();
            IdentifyNearEntitiesAndCheckEntityForIntersection(*(it_elem.base()), leaves);
        }
    }

    // Iterate over conditions
    if (mOptions.Is(FindIntersectedGeometricalObjectsProcess::INTERSECTED_CONDITIONS)) {
        auto& r_conditions_array = mrModelPartIntersected.Conditions();
        const SizeType number_of_conditions = r_conditions_array.size();

        const auto it_conditions_begin = r_conditions_array.begin();

        #pragma omp parallel for private(leaves)
        for (int i = 0; i < static_cast<int>(number_of_conditions); i++) {
            auto it_cond = it_conditions_begin + i;
            leaves.clear();
            IdentifyNearEntitiesAndCheckEntityForIntersection(*(it_cond.base()), leaves);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsProcess::ExecuteInitialize()
{
    GenerateOctree();
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t FindIntersectedGeometricalObjectsProcess::WorkingSpaceDimension()
{
    const auto& r_elements_array = mrModelPartIntersected.Elements();
    if (r_elements_array.size() > 0) {
        const auto it_elements_begin = r_elements_array.begin();
        const auto& r_geometry = it_elements_begin->GetGeometry();
        return r_geometry.WorkingSpaceDimension();
    } else {
        auto& r_conditions_array = mrModelPartIntersected.Conditions();
        const auto it_conditions_begin = r_conditions_array.begin();
        const auto& r_geometry = it_conditions_begin->GetGeometry();
        return r_geometry.WorkingSpaceDimension();
    }
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsProcess::GenerateOctree()
{
    this->SetOctreeBoundingBox();

    // Adding mrModelPart2 to the octree
    for (auto it_node = mrModelPartIntersecting.NodesBegin(); it_node != mrModelPartIntersecting.NodesEnd(); it_node++) {
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
        mpOctree->Insert(it_node->Coordinates().data());

#else
        mpOctree->Insert(it_node->Coordinates().data().data());
#endif // ifdef KRATOS_USE_AMATRIX
    }

    // Add entities
    if (mOptions.Is(FindIntersectedGeometricalObjectsProcess::INTERSECTING_ELEMENTS)) {
        auto& r_elements_array = mrModelPartIntersecting.Elements();
        const auto it_elements_begin = r_elements_array.begin();
        const SizeType number_of_elements = r_elements_array.size();

        // Iterate over the elements
        for (int i = 0; i < static_cast<int>(number_of_elements); i++) {
            auto it_elements = it_elements_begin + i;
            mpOctree->Insert(*(it_elements).base());
        }
    }

    if (mOptions.Is(FindIntersectedGeometricalObjectsProcess::INTERSECTING_CONDITIONS)) {
        auto& r_conditions_array = mrModelPartIntersecting.Conditions();
        const auto it_conditions_begin = r_conditions_array.begin();
        const SizeType number_of_conditions = r_conditions_array.size();

        // Iterate over the conditions
        for (int i = 0; i < static_cast<int>(number_of_conditions); i++) {
            auto it_conditions = it_conditions_begin + i;
            mpOctree->Insert(*(it_conditions).base());
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void  FindIntersectedGeometricalObjectsProcess::SetOctreeBoundingBox()
{
    PointType low(mrModelPartIntersected.NodesBegin()->Coordinates());
    PointType high(mrModelPartIntersected.NodesBegin()->Coordinates());

    // Loop over all nodes in first modelpart
    for (auto it_node = mrModelPartIntersected.NodesBegin(); it_node != mrModelPartIntersected.NodesEnd(); it_node++) {
        const array_1d<double,3>& r_coordinates = it_node->Coordinates();
        for (IndexType i = 0; i < 3; i++) {
            low[i] = r_coordinates[i] < low[i] ? r_coordinates[i] : low[i];
            high[i] = r_coordinates[i] > high[i] ? r_coordinates[i] : high[i];
        }
    }

    // Loop over all skin nodes
    for (auto it_node = mrModelPartIntersecting.NodesBegin(); it_node != mrModelPartIntersecting.NodesEnd(); it_node++) {
        const array_1d<double,3>& r_coordinates = it_node->Coordinates();
        for (IndexType i = 0; i < 3; i++) {
            low[i] = r_coordinates[i] < low[i] ? r_coordinates[i] : low[i];
            high[i] = r_coordinates[i] > high[i] ? r_coordinates[i] : high[i];
        }
    }

    // Slightly increase the bounding box size to avoid problems with geometries in the borders
    // Note that std::numeric_limits<double>::double() is added for the 2D cases. Otherwise, the
    // third component will be 0, breaking the octree behaviour.
    for(IndexType i = 0 ; i < 3; i++) {
        low[i] -= std::abs(high[i] - low[i])*1e-3 + std::numeric_limits<double>::epsilon();
        high[i] += std::abs(high[i] - low[i])*1e-3 + std::numeric_limits<double>::epsilon();
    }

    // TODO: Octree needs refactoring to work with BoundingBox. Pooyan.
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
    mpOctree->SetBoundingBox(low.data(), high.data());
#else
    mpOctree->SetBoundingBox(low.data().data(), high.data().data());
#endif // ifdef KRATOS_USE_AMATRIX
}

/***********************************************************************************/
/***********************************************************************************/

void  FindIntersectedGeometricalObjectsProcess::IdentifyNearEntitiesAndCheckEntityForIntersection(
    GeometricalObject::Pointer pGeometricalObject,
    OtreeCellVectorType& rLeaves
    )
{
    mpOctree->GetIntersectedLeaves(pGeometricalObject, rLeaves);
    MarkIfIntersected(*pGeometricalObject, rLeaves);
}

/***********************************************************************************/
/***********************************************************************************/

void  FindIntersectedGeometricalObjectsProcess::MarkIfIntersected(
    GeometricalObject& rIntersectedElement,
    OtreeCellVectorType& rLeaves
    )
{
    for (auto p_leaf : rLeaves) {
        auto& r_leaf = *(p_leaf->pGetObjects());
        for (auto p_intersecting_entity : r_leaf) {
            if (HasIntersection(rIntersectedElement.GetGeometry(),p_intersecting_entity->GetGeometry())) {
                rIntersectedElement.Set(SELECTED);
                return;
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool FindIntersectedGeometricalObjectsProcess::HasIntersection(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    const IndexType working_space_dimension = rFirstGeometry.WorkingSpaceDimension(); // TODO: DOMAIN_SIZE should be considered for consistency with other implementations
    const IndexType local_space_dimension = rFirstGeometry.LocalSpaceDimension();
    if (working_space_dimension == 2) {
        if (local_space_dimension == 2) {
            return this->HasIntersection2D(rFirstGeometry, rSecondGeometry);
        } else {
            return this->HasDirectIntersection2D(rFirstGeometry, rSecondGeometry);
        }
    } else {
        if (local_space_dimension == 3) {
            return this->HasIntersection3D(rFirstGeometry, rSecondGeometry);
        } else {
            return this->HasDirectIntersection3D(rFirstGeometry, rSecondGeometry);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool FindIntersectedGeometricalObjectsProcess::HasIntersection2D(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    // Check the intersection of each edge against the intersecting object
    const array_1d<double, 3>& r_coordinates_second_geometry_1 = rSecondGeometry[0].Coordinates();
    const array_1d<double, 3>& r_coordinates_second_geometry_2 = rSecondGeometry[1].Coordinates();
    const auto edges = rFirstGeometry.GenerateEdges();
    PointType int_pt(0.0,0.0,0.0);
    for (auto& edge : edges) {
        const int int_id = IntersectionUtilities::ComputeLineLineIntersection<Line2D2<NodeType>>(
            Line2D2<NodeType>{edge},
            r_coordinates_second_geometry_1,
            r_coordinates_second_geometry_2,
            int_pt.Coordinates());

        if (int_id != 0){
            return true;
        }
    }

    // Let check second geometry is inside the first one.
    // Considering that there are no intersection, if one point is inside all of it is inside.
    array_1d<double, 3> local_point;
    if (rFirstGeometry.IsInside(rSecondGeometry.GetPoint(0), local_point)){
        return true;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool FindIntersectedGeometricalObjectsProcess::HasDirectIntersection2D(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    // Check the intersection of each edge against the intersecting object
    const array_1d<double, 3>& r_coordinates_second_geometry_1 = rSecondGeometry[0].Coordinates();
    const array_1d<double, 3>& r_coordinates_second_geometry_2 = rSecondGeometry[1].Coordinates();
    PointType int_pt(0.0,0.0,0.0);
    const int int_id = IntersectionUtilities::ComputeLineLineIntersection<Line2D2<NodeType>>(
        Line2D2<NodeType>{rFirstGeometry},
        r_coordinates_second_geometry_1,
        r_coordinates_second_geometry_2,
        int_pt.Coordinates());

    if (int_id != 0){
        return true;
    }

    // Let check second geometry is inside the first one.
    // Considering that there are no intersection, if one point is inside all of it is inside.
    array_1d<double, 3> local_point;
    if (rFirstGeometry.IsInside(rSecondGeometry.GetPoint(0), local_point)){
        return true;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool FindIntersectedGeometricalObjectsProcess::HasIntersection3D(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    // Check the intersection of each face against the intersecting object
    const auto faces = rFirstGeometry.GenerateFaces();
    for (auto& face : faces) {
        if (face.HasIntersection(rSecondGeometry)){
            return true;
        }
    }

    // Let check second geometry is inside the first one.
    // Considering that there are no intersection, if one point is inside all of it is inside.
    array_1d<double, 3> local_point;
    if (rFirstGeometry.IsInside(rSecondGeometry.GetPoint(0), local_point)){
        return true;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool FindIntersectedGeometricalObjectsProcess::HasDirectIntersection3D(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    // Check the intersection of each face against the intersecting object
    if (rFirstGeometry.HasIntersection(rSecondGeometry)){
        return true;
    }

    // Let check second geometry is inside the first one.
    // Considering that there are no intersection, if one point is inside all of it is inside.
    array_1d<double, 3> local_point;
    if (rFirstGeometry.IsInside(rSecondGeometry.GetPoint(0), local_point)){
        return true;
    }

    return false;
}


/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsProcess::FindIntersectedSkinObjects(
    GeometricalObject& rIntersectedEntity,
    OtreeCellVectorType& rLeaves,
    PointerVector<GeometricalObject>& rResults
    )
{
    for (auto p_leaf : rLeaves) {
        for (auto p_intersecting_entity : *(p_leaf->pGetObjects())) {
            if (HasIntersection(rIntersectedEntity.GetGeometry(), p_intersecting_entity->GetGeometry())) {
                rIntersectedEntity.Set(SELECTED);
                if(std::find(rResults.ptr_begin(), rResults.ptr_end(), p_intersecting_entity) == rResults.ptr_end())
                    rResults.push_back(p_intersecting_entity);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

Parameters FindIntersectedGeometricalObjectsProcess::GetDefaultParameters()
{
    Parameters default_parameters = Parameters(R"(
    {
        "intersected_model_part_name"  : "",
        "intersecting_model_part_name" : "",
        "intersecting_conditions"      : true,
        "intersecting_elements"        : true,
        "intersected_conditions"       : true,
        "intersected_elements"         : true
    })" );

    return default_parameters;
}

}  // namespace Kratos.
