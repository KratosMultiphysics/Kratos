// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:		 BSD License
//					 license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_processes/find_intersected_geometrical_objects_with_obb_for_contact_search_process.h"

namespace Kratos
{
FindIntersectedGeometricalObjectsWithOBBContactSearchProcess::FindIntersectedGeometricalObjectsWithOBBContactSearchProcess(
    ModelPart& rPart1,
    ModelPart& rPart2,
    const double BoundingBoxFactor,
    const Flags Options
    ) : BaseType(rPart1, rPart2, BoundingBoxFactor, Options)
{
}

/***********************************************************************************/
/***********************************************************************************/

FindIntersectedGeometricalObjectsWithOBBContactSearchProcess::FindIntersectedGeometricalObjectsWithOBBContactSearchProcess(
    Model& rModel,
    Parameters ThisParameters
    ) : BaseType(rModel.GetModelPart(ThisParameters["intersected_model_part_name"].GetString()),
        rModel.GetModelPart(ThisParameters["intersecting_model_part_name"].GetString()))
{
    const Parameters default_parameters = GetDefaultParameters();
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    BaseType::mThisParameters = ThisParameters;

    // Checking that the names of the model parts are not empty (this is supposed to be already declared)
    const std::string& r_intersected_model_part_name = BaseType::mThisParameters["intersected_model_part_name"].GetString();
    const std::string& r_intersecting_model_part_name = BaseType::mThisParameters["intersecting_model_part_name"].GetString();

    KRATOS_ERROR_IF(r_intersected_model_part_name == "") << "intersected_model_part_name must be defined on parameters" << std::endl;
    KRATOS_ERROR_IF(r_intersecting_model_part_name == "") << "intersecting_model_part_name must be defined on parameters" << std::endl;

    // Setting flags
    const bool intersecting_conditions = ThisParameters["intersecting_conditions"].GetBool();
    const bool intersecting_elements = ThisParameters["intersecting_elements"].GetBool();
    const bool intersected_conditions = ThisParameters["intersected_conditions"].GetBool();
    const bool intersected_elements = ThisParameters["intersected_elements"].GetBool();

    BaseProcessType::mOptions.Set(BaseProcessType::INTERSECTING_CONDITIONS, intersecting_conditions);
    BaseProcessType::mOptions.Set(BaseProcessType::INTERSECTING_ELEMENTS, intersecting_elements);
    BaseProcessType::mOptions.Set(BaseProcessType::INTERSECTED_CONDITIONS, intersected_conditions);
    BaseProcessType::mOptions.Set(BaseProcessType::INTERSECTED_ELEMENTS, intersected_elements);

    // Setting the bounding box factor
    mBoundingBoxFactor = BaseType::mThisParameters["bounding_box_factor"].GetDouble();

    // In case we consider the bounding box we set the BOUNDARY flag
    if (mBoundingBoxFactor > 0.0)
        this->Set(BOUNDARY, true);
    else
        this->Set(BOUNDARY, false);

    // If we debug OBB
    BaseProcessType::mOptions.Set(BaseType::DEBUG_OBB, mThisParameters["debug_obb"].GetBool());

    // If we build the OBB from the geometry BB
    BaseType::mOptions.Set(BaseType::BUILD_OBB_FROM_BB, mThisParameters["build_from_bounding_box"].GetBool());

    // The intersection type
    ConvertIntersection(mThisParameters["OBB_intersection_type"].GetString());

    // We create new properties for debugging
#ifdef  KRATOS_DEBUG
    if (BaseProcessType::mOptions.Is(FindIntersectedGeometricalObjectsWithOBBProcess::DEBUG_OBB)) {
        this->GetModelPart1().CreateNewProperties(1001);
        this->GetModelPart1().CreateSubModelPart(this->GetModelPart1().Name() + "_AUXILIAR_DEBUG_OBB");
        this->GetModelPart2().CreateNewProperties(1002);
        this->GetModelPart2().CreateSubModelPart(this->GetModelPart2().Name() + "_AUXILIAR_DEBUG_OBB");
    }
#endif

    mLowerBBCoefficient = ThisParameters["lower_bounding_box_coefficient"].GetDouble();
    mHigherBBCoefficient = ThisParameters["higher_bounding_box_coefficient"].GetDouble();
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsWithOBBContactSearchProcess::SetOctreeBoundingBox()
{
    // Getting first iterators
    const auto it_node_begin_1 = BaseType::mrModelPartIntersected.NodesBegin();
    const auto it_node_begin_2 = BaseType::mrModelPartIntersecting.NodesBegin();

    // Setting initial guess
    PointType low(it_node_begin_1->Coordinates());
    PointType high(it_node_begin_1->Coordinates());

    // Loop over all nodes in first modelpart
    for (IndexType i_node = 0 ; i_node < BaseType::mrModelPartIntersected.NumberOfNodes(); ++i_node) {
        auto it_node = it_node_begin_1 + i_node;
        const array_1d<double,3>& r_coordinates = it_node->Coordinates();
        for (IndexType i = 0; i < 3; i++) {
            low[i] = r_coordinates[i] < low[i] ? r_coordinates[i] : low[i];
            high[i] = r_coordinates[i] > high[i] ? r_coordinates[i] : high[i];
        }
    }

    // Loop over all skin nodes
    for (IndexType i_node = 0 ; i_node < BaseType::mrModelPartIntersecting.NumberOfNodes(); ++i_node) {
        auto it_node = it_node_begin_2 + i_node;
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

    // In case we consider the bounding box we extend box
    if (mBoundingBoxFactor > 0.0) {
        const std::size_t dimension = BaseType::mrModelPartIntersected.GetProcessInfo()[DOMAIN_SIZE];
        for(IndexType i = 0 ; i < dimension; i++) {
            if (std::abs(high[i] - low[i]) < mBoundingBoxFactor) {
                low[i] -= mLowerBBCoefficient * mBoundingBoxFactor;
                high[i] += mHigherBBCoefficient * mBoundingBoxFactor;
            }
        }
    }

    // TODO: Octree needs refactoring to work with BoundingBox. Pooyan.
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
    GetOctreePointer()->SetBoundingBox(low.data(), high.data());
#else
    GetOctreePointer()->SetBoundingBox(low.data().data(), high.data().data());
#endif // ifdef KRATOS_USE_AMATRIX
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsWithOBBContactSearchProcess::MarkIfIntersected(
    GeometricalObject& rIntersectedGeometricalObject,
    OtreeCellVectorType& rLeaves
    )
{
    if (rIntersectedGeometricalObject.Is(SLAVE)) {
        // We clear previously assign flag
        rIntersectedGeometricalObject.Reset(SELECTED);

        // The difference with the normal one is that we assign the flag to all "paired" conditions
        for (auto p_leaf : rLeaves) {
            auto& r_leaf = *(p_leaf->pGetObjects());
            // We clear previously assign flags
            for (auto p_condition_2 : r_leaf) {
                if (p_condition_2->IsNot(VISITED)) {
                    p_condition_2->Reset(SELECTED);
                }
            }

            // We iterate again and check intersection
            for (auto p_condition_2 : r_leaf) {
                if (p_condition_2->IsNot(VISITED)) {
                    if (HasIntersection(rIntersectedGeometricalObject.GetGeometry(), p_condition_2->GetGeometry())) {
                        rIntersectedGeometricalObject.Set(SELECTED);
                        p_condition_2->Set(SELECTED);
                    }
                    p_condition_2->Set(VISITED);
                }
            }
        }
        // Reset VISITED flag
        for (auto p_leaf : rLeaves) {
            auto& r_leaf = *(p_leaf->pGetObjects());
            for (auto p_condition_2 : r_leaf) {
                p_condition_2->Reset(VISITED);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters FindIntersectedGeometricalObjectsWithOBBContactSearchProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "intersected_model_part_name"     : "",
        "intersecting_model_part_name"    : "",
        "bounding_box_factor"             : -1.0,
        "debug_obb"                       : false,
        "OBB_intersection_type"           : "SeparatingAxisTheorem",
        "build_from_bounding_box"         : true,
        "lower_bounding_box_coefficient"  : 0.0,
        "higher_bounding_box_coefficient" : 1.0,
        "intersecting_conditions"         : true,
        "intersecting_elements"           : false,
        "intersected_conditions"          : true,
        "intersected_elements"            : false
    })" );

    return default_parameters;
}

}  // namespace Kratos.
