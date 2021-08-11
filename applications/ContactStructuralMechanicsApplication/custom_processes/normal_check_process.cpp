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
#include "utilities/variable_utils.h"
#include "utilities/mortar_utilities.h"
#include "utilities/normal_calculation_utils.h"
#include "utilities/parallel_utilities.h"
#include "custom_processes/normal_check_process.h"
#include "contact_structural_mechanics_application_variables.h"

namespace Kratos
{
void NormalCheckProcess::Execute()
{
    KRATOS_TRY

    // The proportion of length considered in the normal check
    const double length_proportion = mParameters["length_proportion"].GetDouble();

    // The check inside threshold
    const double check_threshold = mParameters["check_threshold"].GetDouble();

    // First we compute the normals
    NormalCalculationUtils().CalculateUnitNormals<Condition>(mrModelPart, true);

    // Getting the nodes array
    auto& r_nodes_array = mrModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();
    const int number_of_nodes = static_cast<int>(r_nodes_array.size());

    // Getting elements array
    auto& r_elements_array = mrModelPart.Elements();
    const auto it_elem_begin = r_elements_array.begin();
    const int number_of_elements = static_cast<int>(r_elements_array.size());

    // Getting conditions array
    auto& r_conditions_array = mrModelPart.Conditions();
    const auto it_cond_begin = r_conditions_array.begin();
    const int number_of_conditions = static_cast<int>(r_conditions_array.size());

    // Backup flags
    std::vector<IndexType> nodes_marker_backup, elements_marker_backup, conditions_marker_backup;
    std::vector<IndexType> nodes_not_marker_backup, elements_not_marker_backup, conditions_not_marker_backup;

    #pragma omp parallel
    {
        // The buffer
        std::vector<IndexType> nodes_marker_backup_buffer, elements_marker_backup_buffer, conditions_marker_backup_buffer;
        std::vector<IndexType> nodes_not_marker_backup_buffer, elements_not_marker_backup_buffer, conditions_not_marker_backup_buffer;

        #pragma omp for schedule(guided, 512) nowait
        for (int i = 0; i < number_of_nodes; ++i) {
            auto it_node = it_node_begin + i;

            if (it_node->IsDefined(MARKER)) {
                if (it_node->Is(MARKER)) {
                    nodes_marker_backup_buffer.push_back(it_node->Id());
                } else {
                    nodes_not_marker_backup_buffer.push_back(it_node->Id());
                }
            }
        }

        #pragma omp for schedule(guided, 512) nowait
        for (int i = 0; i < number_of_elements; ++i) {
            auto it_elem = it_elem_begin + i;

            if (it_elem->IsDefined(MARKER)) {
                if (it_elem->Is(MARKER)) {
                    elements_marker_backup_buffer.push_back(it_elem->Id());
                } else {
                    elements_not_marker_backup_buffer.push_back(it_elem->Id());
                }
            }
        }

        #pragma omp for schedule(guided, 512) nowait
        for (int i = 0; i < number_of_conditions; ++i) {
            auto it_cond = it_cond_begin + i;

            if (it_cond->IsDefined(MARKER)) {
                if (it_cond->Is(MARKER)) {
                    conditions_marker_backup_buffer.push_back(it_cond->Id());
                } else {
                    conditions_not_marker_backup_buffer.push_back(it_cond->Id());
                }
            }
        }

        // Combine buffers together
        #pragma omp critical
        {
            std::move(nodes_marker_backup_buffer.begin(),nodes_marker_backup_buffer.end(),back_inserter(nodes_marker_backup));
            std::move(nodes_not_marker_backup_buffer.begin(),nodes_not_marker_backup_buffer.end(),back_inserter(nodes_not_marker_backup));
            std::move(elements_marker_backup_buffer.begin(),elements_marker_backup_buffer.end(),back_inserter(elements_marker_backup));
            std::move(elements_not_marker_backup_buffer.begin(),elements_not_marker_backup_buffer.end(),back_inserter(elements_not_marker_backup));
            std::move(conditions_marker_backup_buffer.begin(),conditions_marker_backup_buffer.end(),back_inserter(conditions_marker_backup));
            std::move(conditions_not_marker_backup_buffer.begin(),conditions_not_marker_backup_buffer.end(),back_inserter(conditions_not_marker_backup));
        }
    }

    // Reset flags to avoid potential conflict
    VariableUtils().ResetFlag(MARKER, r_nodes_array);
    VariableUtils().ResetFlag(MARKER, r_elements_array);
    VariableUtils().ResetFlag(MARKER, r_conditions_array);

    // Declare auxiliar coordinates
    CoordinatesArrayType aux_coords;

    // Iterate over elements
    block_for_each(r_elements_array, aux_coords, [&length_proportion, &check_threshold](Element& rElem, CoordinatesArrayType& aux_coords) {
        const auto& r_geometry = rElem.GetGeometry();
        bool check_element = true;
        if (r_geometry.WorkingSpaceDimension() > r_geometry.LocalSpaceDimension()) {
            KRATOS_INFO("NormalCheckProcess") << "The element: " <<  rElem.Id() << " is a slender element (beam, shell, membrane...). It will be assumed that the normal is properly oriented" << std::endl;
            check_element = false;
        }

        // Only solid elements
        if (check_element) {
            const GeometryType::GeometriesArrayType& r_boundary = r_geometry.GenerateBoundariesEntities();
            for (GeometryType& r_face : r_boundary) {
                for (auto& r_node : r_face) {
                    if (r_node.IsNot(INTERFACE)) {
                        continue;
                    }
                }
                // First check if the elements are inverted
                r_face.PointLocalCoordinates(aux_coords, r_face.Center());
                const array_1d<double, 3> normal = r_face.UnitNormal(aux_coords);
                CoordinatesArrayType aux_perturbed_coords = r_face.Center() + length_proportion * r_face.Length() * normal;
                if (r_geometry.IsInside(aux_perturbed_coords, aux_coords, check_threshold)) {
                    rElem.Set(MARKER);
                    for (auto& r_node : r_face) {
                        r_node.Set(MARKER);
                    }
                    KRATOS_INFO("NormalCheckProcess") << "Normal inverted in element: " << rElem.Id() << " the corresponding element will be inverted" << std::endl;
                    continue; // Element inverted just once
                }
            }
        }
    });

    // Check conditions
    block_for_each(r_conditions_array, [&](Condition& rCond) {
        rCond.Set(MARKER);
        const auto& r_geometry = rCond.GetGeometry();

        for (auto& r_node : r_geometry) {
            if (r_node.IsNot(MARKER)) {
                rCond.Set(MARKER, false);
                break;
            }
        }
    });

    // Invert elements
    MortarUtilities::InvertNormalForFlag<PointerVectorSet<Element, IndexedObject>>(r_elements_array, MARKER);

    // Invert conditions
    MortarUtilities::InvertNormalForFlag<PointerVectorSet<Condition, IndexedObject>>(r_conditions_array, MARKER);

    // We re-compute the normals
    NormalCalculationUtils().CalculateUnitNormals<Condition>(mrModelPart, true);

    // Reset flags
    VariableUtils().ResetFlag(MARKER, r_nodes_array);
    VariableUtils().ResetFlag(MARKER, r_elements_array);
    VariableUtils().ResetFlag(MARKER, r_conditions_array);

    // After recompute the normals we check also the conditions
    block_for_each(r_elements_array, aux_coords, [&length_proportion, &check_threshold](Element& rElem, CoordinatesArrayType& aux_coords) {
        const auto& r_geometry = rElem.GetGeometry();
        bool check_element = true;
        if (r_geometry.WorkingSpaceDimension() > r_geometry.LocalSpaceDimension()) {
            KRATOS_INFO("NormalCheckProcess") << "The element: " <<  rElem.Id() << " is a slender element (beam, shell, membrane...). It will be assumed that the normal is properly oriented" << std::endl;
            check_element = false;
        }

        // Only solid elements
        if (check_element) {
            // Check nodes resulting normal, so the conditions are inverted
            for (auto& r_node : r_geometry) {
                if (r_node.Is(INTERFACE)) {
                    const array_1d<double, 3>& r_normal = r_node.FastGetSolutionStepValue(NORMAL);
                    CoordinatesArrayType aux_perturbed_coords = r_node.Coordinates() + length_proportion * r_geometry.Length() * r_normal;

                    if (r_geometry.IsInside(aux_perturbed_coords, aux_coords, check_threshold)) {
                        r_node.SetLock();
                        r_node.Set(MARKER);
                        r_node.UnSetLock();
                        KRATOS_INFO("NormalCheckProcess") << "Normal inverted in node: " << r_node.Id() << " the corresponding condition will be inverted" << std::endl;
                    }
                }
            }
        }
    });

    // Check conditions
    block_for_each(r_conditions_array, [&](Condition& rCond) {
        rCond.Set(MARKER);
        const auto& r_geometry = rCond.GetGeometry();

        for (auto& r_node : r_geometry) {
            if (r_node.IsNot(MARKER)) {
                rCond.Set(MARKER, false);
                break;
            }
        }
    });

    // Invert conditions
    MortarUtilities::InvertNormalForFlag<PointerVectorSet<Condition, IndexedObject>>(r_conditions_array, MARKER);

    // Reset flags
    VariableUtils().ResetFlag(MARKER, r_nodes_array);
    VariableUtils().ResetFlag(MARKER, r_conditions_array);

    // We re-compute the normals
    NormalCalculationUtils().CalculateUnitNormals<Condition>(mrModelPart, true);

    // Reassign flags
    block_for_each(nodes_marker_backup, [&](IndexType& rId) {
        mrModelPart.pGetNode(rId)->Set(MARKER, false);
    });
    block_for_each(nodes_not_marker_backup, [&](IndexType& rId) {
        mrModelPart.pGetNode(rId)->Set(MARKER, false);
    });
    block_for_each(elements_marker_backup, [&](IndexType& rId) {
        mrModelPart.pGetElement(rId)->Set(MARKER, false);
    });
    block_for_each(elements_not_marker_backup, [&](IndexType& rId) {
        mrModelPart.pGetElement(rId)->Set(MARKER, false);
    });
    block_for_each(conditions_marker_backup, [&](IndexType& rId) {
        mrModelPart.pGetCondition(rId)->Set(MARKER, false);
    });
    block_for_each(conditions_not_marker_backup, [&](IndexType& rId) {
        mrModelPart.pGetCondition(rId)->Set(MARKER, false);
    });

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters NormalCheckProcess::GetDefaultParameters() const
{
    KRATOS_TRY

    const Parameters default_parameters = Parameters(R"(
    {
        "length_proportion" : 0.1,
        "check_threshold"   : 5.0e-7
    })" );

    return default_parameters;

    KRATOS_CATCH("")
}

}  // namespace Kratos.
