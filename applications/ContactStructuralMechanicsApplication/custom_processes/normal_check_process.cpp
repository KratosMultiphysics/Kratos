// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/variable_utils.h"
#include "utilities/mortar_utilities.h"
#include "custom_processes/normal_check_process.h"
#include "contact_structural_mechanics_application_variables.h"

namespace Kratos
{
void NormalCheckProcess::Execute()
{
    KRATOS_TRY

    // First we compute the normals
    MortarUtilities::ComputeNodesMeanNormalModelPart(mrModelPart);

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
    CoordinatesArrayType aux_coords, aux_perturbed_coords;

    // Auxiliar boolean
    bool boundary_face = true;

    #pragma omp parallel for firstprivate(boundary_face,aux_coords,aux_perturbed_coords)
    for(int i = 0; i < number_of_elements; ++i) {
        auto it_elem = it_elem_begin + i;
        const auto& r_geometry = it_elem->GetGeometry();
        if (r_geometry.WorkingSpaceDimension() > r_geometry.LocalSpaceDimension()) {
            KRATOS_INFO("NormalCheckProcess") << "The element: " <<  it_elem->Id() << " is a slender element (beam, shell, membrane...). It will be assumed that the normal is properly oriented" << std::endl;
            continue;
        }
        const GeometryType::GeometriesArrayType& r_boundary = r_geometry.GenerateBoundariesEntities();
        for (GeometryType& r_face : r_boundary) {
            for (auto& r_node : r_face) {
                if (r_node.IsNot(INTERFACE)) {
                    boundary_face = false;
                    break;
                }
            }
            if (boundary_face) {
                // First check if the elements are inverted
                r_face.PointLocalCoordinates(aux_coords, r_face.Center());
                const array_1d<double, 3> normal = r_face.UnitNormal(aux_coords);
                aux_perturbed_coords = r_face.Center() + 1.0e-1 * r_face.Length() * normal;
                if (r_geometry.IsInside(aux_perturbed_coords, aux_coords, 5.0e-7)) {
                    it_elem->Set(MARKER);
                    KRATOS_INFO("NormalCheckProcess") << "Normal inverted in element: " << it_elem->Id() << " the corresponding element will be inverted" << std::endl;
                }
            }
        }

        // Check nodes resulting normal, so the conditions are inverted
        for (auto& r_node : r_geometry) {
            if (r_node.Is(INTERFACE)) {
                const array_1d<double, 3>& r_normal = r_node.FastGetSolutionStepValue(NORMAL);
                aux_perturbed_coords = r_node.Coordinates() + 1.0e-1 * r_geometry.Length() * r_normal;

                if (r_geometry.IsInside(aux_perturbed_coords, aux_coords, 5.0e-7)) {
                    r_node.SetLock();
                    r_node.Set(MARKER);
                    r_node.UnSetLock();
                    KRATOS_INFO("NormalCheckProcess") << "Normal inverted in node: " << r_node.Id() << " the corresponding condition will be inverted" << std::endl;
                }
            }
        }

        boundary_face = true;
    }

    // Invert elements
    MortarUtilities::InvertNormalForFlag<PointerVectorSet<Element, IndexedObject>>(r_elements_array, MARKER);

    // Auxiliar boolean
    bool inverted_condition = true;

    #pragma omp parallel for firstprivate(inverted_condition)
    for(int i = 0; i < number_of_conditions; ++i) {
        auto it_cond = it_cond_begin + i;
        const auto& r_geometry = it_cond->GetGeometry();

        for (auto& r_node : r_geometry) {
            if (r_node.IsNot(MARKER)) {
                inverted_condition = false;
                break;
            }
        }
        if (inverted_condition) {
            it_cond->Set(MARKER);
        }

        inverted_condition = true;
    }

    // Invert conditions
    MortarUtilities::InvertNormalForFlag<PointerVectorSet<Condition, IndexedObject>>(r_conditions_array, MARKER);

    // Reset flags
    VariableUtils().ResetFlag(MARKER, r_nodes_array);
    VariableUtils().ResetFlag(MARKER, r_elements_array);
    VariableUtils().ResetFlag(MARKER, r_conditions_array);

    // Reassign flags
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(nodes_marker_backup.size()); ++i) {
        mrModelPart.pGetNode(nodes_marker_backup[i])->Set(MARKER, true);
    }
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(nodes_not_marker_backup.size()); ++i) {
        mrModelPart.pGetNode(nodes_not_marker_backup[i])->Set(MARKER, false);
    }
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(elements_marker_backup.size()); ++i) {
        mrModelPart.pGetElement(elements_marker_backup[i])->Set(MARKER, true);
    }
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(elements_not_marker_backup.size()); ++i) {
        mrModelPart.pGetElement(elements_not_marker_backup[i])->Set(MARKER, false);
    }
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(conditions_marker_backup.size()); ++i) {
        mrModelPart.pGetCondition(conditions_marker_backup[i])->Set(MARKER, true);
    }
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(conditions_not_marker_backup.size()); ++i) {
        mrModelPart.pGetCondition(conditions_not_marker_backup[i])->Set(MARKER, false);
    }

    // We re-compute the normals
    MortarUtilities::ComputeNodesMeanNormalModelPart(mrModelPart);

    KRATOS_CATCH("")
}

}  // namespace Kratos.
