//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


// System includes


// External includes


// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "find_intersected_objects_utility.h"

namespace Kratos
{

/// Local Flags
KRATOS_CREATE_LOCAL_FLAG(FindIntersectedObjectsUtility, INTERSECT_ELEMENTS,   0);
KRATOS_CREATE_LOCAL_FLAG(FindIntersectedObjectsUtility, INTERSECT_CONDITIONS, 1);

FindIntersectedObjectsUtility::FindIntersectedObjectsUtility(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
) : mrModelPart(rThisModelPart),
    mpOctree(new OctreeType())
{
    Parameters default_parameters(R"({
        "intersect_elements"   : true,
        "intersect_conditions" : false
    })");
    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    // Set the flags
    const bool intersect_elements = ThisParameters["intersect_elements"].GetBool();
    const bool intersect_conditions = ThisParameters["intersect_conditions"].GetBool();

    mOptions.Set(INTERSECT_ELEMENTS, intersect_elements);
    mOptions.Set(INTERSECT_CONDITIONS, intersect_conditions);

    if (mrModelPart.NumberOfNodes() > 0) {
        GenerateOctree();
    }
}

void FindIntersectedObjectsUtility::UpdateSearchStructure()
{
    if (mrModelPart.NumberOfNodes() > 0) {
        GenerateOctree();
    }
}

void FindIntersectedObjectsUtility::FindIntersectedObjects(
    GeometryType::Pointer pGeometry,
    PointerVector<GeometricalObject>& rResults)
{
    OctreeCellVectorType leaves;
    leaves.clear();
    auto p_geometrical_object = Kratos::make_intrusive<GeometricalObject>(0, pGeometry);
    mpOctree->GetIntersectedLeaves(p_geometrical_object, leaves);
    FindIntersectedObjects(*pGeometry, leaves, rResults);
}

void FindIntersectedObjectsUtility::FindIntersectedObjects(
    const GeometryType& rIntersectionGeometry,
    OctreeCellVectorType& rLeaves,
    PointerVector<GeometricalObject>& rResults)
{
    for (auto p_leaf : rLeaves) {
        for (auto p_volume_object : *(p_leaf->pGetObjects())) {
            if (p_volume_object->GetGeometry().HasIntersection(rIntersectionGeometry)) {
                rResults.push_back(p_volume_object);
            }
        }
    }
}

void FindIntersectedObjectsUtility::GenerateOctree()
{
    SetOctreeBoundingBox();

    // Adding the volume model part to the octree
    for (auto it_node = mrModelPart.NodesBegin(); it_node != mrModelPart.NodesEnd(); it_node++) {
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
        mpOctree->Insert(it_node->Coordinates().data());
#else
        mpOctree->Insert(it_node->Coordinates().data().data());
#endif // ifdef KRATOS_USE_AMATRIX
    }

    // Add elements
    if (mOptions.Is(INTERSECT_ELEMENTS)) {
        auto& r_elements_array = mrModelPart.Elements();
        const auto it_elements_begin = r_elements_array.begin();
        const auto number_of_elements = r_elements_array.size();

        // Iterate over the elements
        for (int i = 0; i < static_cast<int>(number_of_elements); i++) {
            auto it_elements = it_elements_begin + i;
            mpOctree->Insert(*(it_elements).base());
        }
    }

    // Add conditions
    if (mOptions.Is(INTERSECT_CONDITIONS)) {
        auto& r_conditions_array = mrModelPart.Conditions();
        const auto it_conditions_begin = r_conditions_array.begin();
        const auto number_of_conditions = r_conditions_array.size();

        // Iterate over the conditions
        for (int i = 0; i < static_cast<int>(number_of_conditions); i++) {
            auto it_conditions = it_conditions_begin + i;
            mpOctree->Insert(*(it_conditions).base());
        }
    }
}

void FindIntersectedObjectsUtility::SetOctreeBoundingBox()
{
    PointType low(mrModelPart.NodesBegin()->Coordinates());
    PointType high(mrModelPart.NodesBegin()->Coordinates());

    // Loop over all nodes in the volume model part
    for (auto it_node = mrModelPart.NodesBegin(); it_node != mrModelPart.NodesEnd(); it_node++) {
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

}  // namespace Kratos.
