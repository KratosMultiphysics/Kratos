//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    DamApplication developers
//

#if !defined(KRATOS_DAM_NODAL_CAUCHY_STRESS_EXTRAPOLATION_PROCESS_H_INCLUDED)
#define KRATOS_DAM_NODAL_CAUCHY_STRESS_EXTRAPOLATION_PROCESS_H_INCLUDED

// Project includes
#include "processes/process.h"
#include "includes/kratos_flags.h"
#include "includes/model_part.h"
#include "utilities/atomic_utilities.h"
#include "utilities/math_utils.h"
#include "utilities/parallel_utilities.h"
#include "custom_utilities/poro_element_utilities.hpp"

// Application includes
#include "dam_application_variables.h"

namespace Kratos
{

/**
 * @class DamNodalCauchyStressExtrapolationProcess
 * @ingroup DamApplication
 * @brief Extrapolates the Gauss-point Cauchy stress tensor to the nodes,
 * reproducing exactly the element-level accumulation currently embedded in
 * SmallDisplacementThermoMechanicElement::FinalizeSolutionStep.
 *
 * The process works with any element that provides CAUCHY_STRESS_TENSOR at its
 * integration points through the standard
 * Element::CalculateOnIntegrationPoints interface (e.g. the
 * StructuralMechanicsApplication small-displacement elements). It reproduces the
 * legacy element-level extrapolation operator:
 *   Triangle2D3 (1 GP) and Tetrahedra3D4 (1 GP): the single Gauss-point stress
 *   is copied to every node;
 *   Quadrilateral2D4 (4 GP): E_Q4 * sigma  (PoroElementUtilities::Calculate2DExtrapolationMatrix);
 *   Hexahedra3D8 (8 GP):  E_H8 * sigma  (PoroElementUtilities::Calculate3DExtrapolationMatrix);
 * following exactly the legacy arithmetic path: the Gauss-point stress tensors
 * are stored in a Voigt stress container with the legacy component ordering, the
 * element extrapolation is computed as
 *   prod(ExtrapolationMatrix, StressContainer)
 * (the same ublas operation used by the legacy element), and every extrapolated
 * nodal row is converted back to a tensor with
 *   MathUtils<double>::StressVectorToTensor.
 *
 * The raw area-weighted accumulation is then performed with the element measure
 * (Geometry::Area(), which returns the volume for 3D solids):
 *   NODAL_CAUCHY_STRESS_TENSOR += element_measure * extrapolated_stress
 *   NODAL_AREA                 += element_measure
 * using historical nodal storage (FastGetSolutionStepValue).
 *
 * Responsibility boundaries:
 * - The process only EXTRAPOLATES AND ACCUMULATES (see ExtrapolateAndAccumulate).
 *   It does NOT reset NODAL_AREA / NODAL_CAUCHY_STRESS_TENSOR and it does NOT
 *   normalize them. In the legacy workflow the reset and the normalization are
 *   owned by IncrementalUpdateStaticSmoothingScheme::FinalizeSolutionStep, and
 *   the caller is responsible for initializing the accumulators before invoking
 *   ExtrapolateAndAccumulate.
 * - ExecuteFinalizeSolutionStep is provided only as a Process API entry point
 *   for standalone/manual use (and for the characterization tests). A standard
 *   AnalysisStage calls process ExecuteFinalizeSolutionStep AFTER
 *   Solver::FinalizeSolutionStep, whereas the legacy extrapolation runs inside
 *   the solver/scheme finalization sequence between the nodal reset and the
 *   nodal normalization. It is therefore NOT yet the definitive production
 *   execution point; the phase-3D.3 integration is expected to invoke the
 *   extrapolation from the Dam smoothing workflow at the same logical position
 *   previously occupied by the element-based accumulation.
 */
class DamNodalCauchyStressExtrapolationProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(DamNodalCauchyStressExtrapolationProcess);

    /// Constructor
    explicit DamNodalCauchyStressExtrapolationProcess(ModelPart& rModelPart)
        : Process(Flags()), mrModelPart(rModelPart) {}

    /// Destructor
    ~DamNodalCauchyStressExtrapolationProcess() override = default;

    /**
     * @brief Process API entry point (standalone/manual use only; not the final
     * production placement).
     * @note Assumes the nodal accumulators have been initialized by the workflow
     * owner before this is called.
     */
    void ExecuteFinalizeSolutionStep() override
    {
        ExtrapolateAndAccumulate();
    }

    /**
     * @brief The reusable production operation: extrapolates the integration-point
     * Cauchy stress to the nodes and performs the raw area-weighted accumulation.
     * @note The caller (e.g. the Dam smoothing workflow in phase 3D.3) must have
     * initialized NODAL_CAUCHY_STRESS_TENSOR and NODAL_AREA before calling this.
     */
    void ExtrapolateAndAccumulate()
    {
        KRATOS_TRY

        const std::size_t dimension = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

        block_for_each(mrModelPart.Elements(), [this, &dimension](Element& rElement) {
            if (!rElement.IsActive()) return;

            auto& r_geometry = rElement.GetGeometry();
            const std::size_t number_of_nodes = r_geometry.size();

            // Gauss-point Cauchy stress tensor (2x2 in 2D, 3x3 in 3D).
            std::vector<Matrix> gauss_stress;
            rElement.CalculateOnIntegrationPoints(
                CAUCHY_STRESS_TENSOR, gauss_stress, mrModelPart.GetProcessInfo());
            const std::size_t number_of_gps = gauss_stress.size();

            // Skip elements that do not provide a Gauss-point Cauchy stress
            // (e.g. interface elements). This mirrors the legacy workflow, in
            // which only the small-displacement thermo-mechanical elements
            // accumulate the nodal stress and every other element does not
            // participate in the nodal smoothing.
            if (number_of_gps == 0) {
                return;
            }

            // Element-level extrapolation operator (matching the legacy element).
            bool single_gp = false;
            Matrix extrapolation_matrix;
            if (dimension == 2) {
                if (number_of_nodes == 3) {
                    single_gp = true;
                } else if (number_of_nodes == 4) {
                    BoundedMatrix<double, 4, 4> q4;
                    PoroElementUtilities::Calculate2DExtrapolationMatrix(q4);
                    extrapolation_matrix = Matrix(q4);
                } else {
                    KRATOS_ERROR << "DamNodalCauchyStressExtrapolationProcess: unsupported 2D geometry with "
                                 << number_of_nodes << " nodes (supported: Triangle2D3, Quadrilateral2D4)"
                                 << std::endl;
                }
            } else if (dimension == 3) {
                if (number_of_nodes == 4) {
                    single_gp = true;
                } else if (number_of_nodes == 8) {
                    BoundedMatrix<double, 8, 8> h8;
                    PoroElementUtilities::Calculate3DExtrapolationMatrix(h8);
                    extrapolation_matrix = Matrix(h8);
                } else {
                    KRATOS_ERROR << "DamNodalCauchyStressExtrapolationProcess: unsupported 3D geometry with "
                                 << number_of_nodes << " nodes (supported: Tetrahedra3D4, Hexahedra3D8)"
                                 << std::endl;
                }
            } else {
                KRATOS_ERROR << "DamNodalCauchyStressExtrapolationProcess: unsupported dimension "
                             << dimension << std::endl;
            }

            // Verify that the element integration scheme matches the legacy one.
            const std::size_t expected_number_of_gps =
                single_gp ? 1 : extrapolation_matrix.size2();
            KRATOS_ERROR_IF_NOT(number_of_gps == expected_number_of_gps)
                << "DamNodalCauchyStressExtrapolationProcess: element " << rElement.Id()
                << " provides " << number_of_gps
                << " integration points but the legacy extrapolation expects "
                << expected_number_of_gps
                << " for this geometry (incompatible integration scheme)." << std::endl;

            // Build the Gauss-point Voigt stress container with the legacy
            // component ordering, exactly as SmallDisplacementThermoMechanicElement
            // does (2D: [sxx, syy, sxy]; 3D: [sxx, syy, szz, sxy, syz, sxz]).
            const std::size_t voigt_size = (dimension == 2) ? 3 : 6;
            Matrix stress_container(number_of_gps, voigt_size);
            noalias(stress_container) = ZeroMatrix(number_of_gps, voigt_size);
            for (std::size_t gp = 0; gp < number_of_gps; ++gp) {
                stress_container(gp, 0) = gauss_stress[gp](0, 0);
                stress_container(gp, 1) = gauss_stress[gp](1, 1);
                if (dimension == 2) {
                    stress_container(gp, 2) = gauss_stress[gp](0, 1);
                } else {
                    stress_container(gp, 2) = gauss_stress[gp](2, 2);
                    stress_container(gp, 3) = gauss_stress[gp](0, 1);
                    stress_container(gp, 4) = gauss_stress[gp](1, 2);
                    stress_container(gp, 5) = gauss_stress[gp](0, 2);
                }
            }

            // Element-level extrapolation using the same ublas operation as the
            // legacy element (for single-GP geometries every node receives the
            // single Gauss-point stress).
            Matrix aux_nodal_stress;
            if (single_gp) {
                aux_nodal_stress = ZeroMatrix(number_of_nodes, voigt_size);
                for (std::size_t i = 0; i < number_of_nodes; ++i) {
                    for (std::size_t comp = 0; comp < voigt_size; ++comp) {
                        aux_nodal_stress(i, comp) = stress_container(0, comp);
                    }
                }
            } else {
                aux_nodal_stress = prod(extrapolation_matrix, stress_container);
            }

            // Element measure used by the legacy element (Geometry::Area(); the
            // volume for 3D solids).
            const double element_measure = r_geometry.Area();

            // Convert every extrapolated nodal row back to a tensor and accumulate
            // with the element measure.
            for (std::size_t i = 0; i < number_of_nodes; ++i) {
                Vector nodal_voigt(voigt_size);
                for (std::size_t comp = 0; comp < voigt_size; ++comp) {
                    nodal_voigt(comp) = aux_nodal_stress(i, comp);
                }
                const Matrix nodal_tensor = MathUtils<double>::StressVectorToTensor(nodal_voigt);

                Matrix& r_nodal_stress = r_geometry[i].FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
                for (std::size_t a = 0; a < dimension; ++a) {
                    for (std::size_t b = 0; b < dimension; ++b) {
                        AtomicAdd(r_nodal_stress(a, b), element_measure * nodal_tensor(a, b));
                    }
                }
                AtomicAdd(r_geometry[i].FastGetSolutionStepValue(NODAL_AREA), element_measure);
            }
        });

        KRATOS_CATCH("")
    }

private:
    ModelPart& mrModelPart;
};

} // namespace Kratos

#endif // KRATOS_DAM_NODAL_CAUCHY_STRESS_EXTRAPOLATION_PROCESS_H_INCLUDED
