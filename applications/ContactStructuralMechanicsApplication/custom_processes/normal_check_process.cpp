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

    // Getting elements array
    auto& r_elements_array = mrModelPart.Elements();
    const auto it_elem_begin = r_elements_array.begin();

    // Declare auxiliar coordinates
    CoordinatesArrayType aux_coords, aux_perturbed_coords;

    // Auxiliar boolean
    bool boundary_face = true;

    #pragma omp parallel for firstprivate(boundary_face,aux_coords,aux_perturbed_coords)
    for(int i = 0; i < static_cast<int>(r_elements_array.size()); ++i) {
        auto it_elem = it_elem_begin + i;
        const auto& r_geometry = it_elem->GetGeometry();
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
                if (r_geometry.IsInside(aux_perturbed_coords, aux_coords, 1.0e-12)) {
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

                if (r_geometry.IsInside(aux_perturbed_coords, aux_coords, 1.0e-12)) {
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

    // Getting conditions array
    auto& r_conditions_array = mrModelPart.Conditions();
    const auto it_cond_begin = r_conditions_array.begin();

    // Auxiliar boolean
    bool inverted_condition = true;

    #pragma omp parallel for firstprivate(inverted_condition)
    for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
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
    VariableUtils().ResetFlag(MARKER, mrModelPart.Nodes());
    VariableUtils().ResetFlag(MARKER, r_elements_array);
    VariableUtils().ResetFlag(MARKER, r_conditions_array);

    // We re-compute the normals
    MortarUtilities::ComputeNodesMeanNormalModelPart(mrModelPart);

    KRATOS_CATCH("")
}

}  // namespace Kratos.
