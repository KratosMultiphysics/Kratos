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
#include "utilities/parallel_utilities.h"

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
    KRATOS_WARNING("Using FindIntersectedGeometricalObjectsProcess.Initialize() method is deprecated. Please use the ExecuteInitialize() method instead.");
    this->ExecuteInitialize();
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

    // Iterate over elements
    if (mOptions.Is(FindIntersectedGeometricalObjectsProcess::INTERSECTED_ELEMENTS)) {
        auto& r_elements_array = mrModelPartIntersected.Elements();
        const SizeType number_of_elements = r_elements_array.size();

        const auto it_elements_begin = r_elements_array.begin();

        IndexPartition<std::size_t>(number_of_elements).for_each(OtreeCellVectorType() ,[&](std::size_t index, OtreeCellVectorType& rLocalLeaves){
            auto it_elem = it_elements_begin + index;
            rLocalLeaves.clear();
            IdentifyNearEntitiesAndCheckEntityForIntersection(*(it_elem.base()), rLocalLeaves);
        });
    }

    // Iterate over conditions
    if (mOptions.Is(FindIntersectedGeometricalObjectsProcess::INTERSECTED_CONDITIONS)) {
        auto& r_conditions_array = mrModelPartIntersected.Conditions();
        const SizeType number_of_conditions = r_conditions_array.size();

        const auto it_conditions_begin = r_conditions_array.begin();

        IndexPartition<std::size_t>(number_of_conditions).for_each(OtreeCellVectorType(), [&](std::size_t index, OtreeCellVectorType& rLocalLeaves){
            auto it_cond = it_conditions_begin + index;
            rLocalLeaves.clear();
            IdentifyNearEntitiesAndCheckEntityForIntersection(*(it_cond.base()), rLocalLeaves);
        });
    }
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsProcess::ExecuteInitialize()
{
    if (mrModelPartIntersected.NumberOfNodes() > 0) {
        GenerateOctree();
    }
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
        mpOctree->Insert(it_node->Coordinates().data().data());
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
    mpOctree->SetBoundingBox(low.data().data(), high.data().data());
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
        for (const auto& p_intersecting_entity : r_leaf) {
            if (HasIntersection(rIntersectedElement.GetGeometry(), p_intersecting_entity->GetGeometry())) {
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
    return rFirstGeometry.HasIntersection(rSecondGeometry);
}

void FindIntersectedGeometricalObjectsProcess::FindIntersectedSkinObjects(
    GeometricalObject& rIntersectedEntity,
    OtreeCellVectorType& rLeaves,
    PointerVector<GeometricalObject>& rResults
    )
{
    for (auto p_leaf : rLeaves) {
        for (const auto& p_intersecting_entity : *(p_leaf->pGetObjects())) {
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

const Parameters FindIntersectedGeometricalObjectsProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
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
