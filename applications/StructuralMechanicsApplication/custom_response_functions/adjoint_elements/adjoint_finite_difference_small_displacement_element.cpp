// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//


#include "adjoint_finite_difference_small_displacement_element.h"
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/response_utilities/stress_response_definitions.h"
#include "includes/checks.h"


namespace Kratos
{

AdjointFiniteDifferenceSmallDisplacementElement::AdjointFiniteDifferenceSmallDisplacementElement(Element::Pointer pPrimalElement)
    : AdjointFiniteDifferencingBaseElement(pPrimalElement)
{
}

AdjointFiniteDifferenceSmallDisplacementElement::~AdjointFiniteDifferenceSmallDisplacementElement()
{
}

void AdjointFiniteDifferenceSmallDisplacementElement::CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                                                                            Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // First check applicability of this function
    KRATOS_ERROR_IF(rStressVariable != STRESS_ON_GP)
        << "AdjointFiniteDifferenceSmallDisplacementElement::CalculateStressDisplacementDerivative: Invalid stress variable! Stress variable not supported for this element!" << std::endl;

    TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));
    KRATOS_ERROR_IF(traced_stress_type != TracedStressType::VON_MISES_STRESS)
        << "AdjointFiniteDifferenceSmallDisplacementElement::CalculateStressDisplacementDerivative: Invalid stress type! Stress type not supported for this element!" << std::endl;

    KRATOS_ERROR_IF(rCurrentProcessInfo.Has(NL_ITERATION_NUMBER))
        << "This stress displacement derivative computation is only usable for linear cases!" << std::endl;

    // Some working variables
    const SizeType num_nodes = mpPrimalElement->GetGeometry().PointsNumber();
    const SizeType dimension = mpPrimalElement->GetGeometry().WorkingSpaceDimension();
    const SizeType num_dofs_per_node = (mHasRotationDofs) ?  2 * dimension : dimension;
    const SizeType num_dofs = num_nodes * num_dofs_per_node;

    // Build vector of variables containing the DOF-variables of the primal problem
    std::vector<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>> primal_solution_variable_list;
    primal_solution_variable_list.reserve(num_dofs_per_node);
    primal_solution_variable_list.push_back(DISPLACEMENT_X);
    primal_solution_variable_list.push_back(DISPLACEMENT_Y);
    primal_solution_variable_list.push_back(DISPLACEMENT_Z);

    std::vector<Matrix> stress_tensor;
    mpPrimalElement->CalculateOnIntegrationPoints(PK2_STRESS_TENSOR, stress_tensor, rCurrentProcessInfo);

    const unsigned int number_integration_points = stress_tensor.size();

    // Computation of prefactors
    std::vector<double> sensitivity_prefactors(number_integration_points);

    for(IndexType k = 0; k < number_integration_points; ++k)
    {
        Matrix & PK2_k = stress_tensor[k];
        double radicant = 0.0;
        radicant += PK2_k(0,0)*PK2_k(0,0) + PK2_k(1,1)*PK2_k(1,1) + PK2_k(2,2)*PK2_k(2,2);
        radicant -= PK2_k(0,0)*PK2_k(1,1) + PK2_k(0,0)*PK2_k(2,2) + PK2_k(1,1)*PK2_k(2,2);
        radicant += 3*PK2_k(0,1)*PK2_k(0,1) + 3*PK2_k(0,2)*PK2_k(0,2) + 3*PK2_k(1,2)*PK2_k(1,2);
        sensitivity_prefactors[k] = 0.5/sqrt(radicant);
    }

    // Store primal results and initialize deformation
    Vector initial_state_variables;
    initial_state_variables.resize(num_dofs, false);

    for (IndexType i = 0; i < num_nodes; ++i)
    {
        const IndexType index = i * num_dofs_per_node;
        for(IndexType j = 0; j < primal_solution_variable_list.size(); ++j)
        {
            initial_state_variables[index + j] = mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]);
            mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = 0.0;
        }
    }

    // Compute gradient using unit deformation states
    rOutput.resize(num_dofs, number_integration_points, false);
    rOutput.clear();

    std::vector<Matrix> partial_stress_derivatives;

    for (IndexType i = 0; i < num_nodes; ++i)
    {
        const IndexType index = i * num_dofs_per_node;

        for(IndexType j = 0; j < primal_solution_variable_list.size(); ++j)
        {
            mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = 1.0;
            mpPrimalElement->CalculateOnIntegrationPoints(PK2_STRESS_TENSOR, partial_stress_derivatives, rCurrentProcessInfo);

            for(IndexType k = 0; k < number_integration_points; ++k)
            {
                double sensitivity_entry_k = 0.0;
                Matrix & PK2_k = stress_tensor[k];

                sensitivity_entry_k += 2*PK2_k(0,0)*partial_stress_derivatives[k](0,0);
                sensitivity_entry_k += 2*PK2_k(1,1)*partial_stress_derivatives[k](1,1);
                sensitivity_entry_k += 2*PK2_k(2,2)*partial_stress_derivatives[k](2,2);

                sensitivity_entry_k -= PK2_k(1,1)*partial_stress_derivatives[k](0,0);
                sensitivity_entry_k -= PK2_k(0,0)*partial_stress_derivatives[k](1,1);

                sensitivity_entry_k -= PK2_k(0,0)*partial_stress_derivatives[k](2,2);
                sensitivity_entry_k -= PK2_k(2,2)*partial_stress_derivatives[k](0,0);

                sensitivity_entry_k -= PK2_k(1,1)*partial_stress_derivatives[k](2,2);
                sensitivity_entry_k -= PK2_k(2,2)*partial_stress_derivatives[k](1,1);

                sensitivity_entry_k += 6*PK2_k(0,1)*partial_stress_derivatives[k](0,1);
                sensitivity_entry_k += 6*PK2_k(0,2)*partial_stress_derivatives[k](0,2);
                sensitivity_entry_k += 6*PK2_k(1,2)*partial_stress_derivatives[k](1,2);

                sensitivity_entry_k *= sensitivity_prefactors[k];

                rOutput(index+j, k) = sensitivity_entry_k;
            }

            mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = 0.0;
        }
    }

    // Recall primal solution
    for (IndexType i = 0; i < num_nodes; ++i)
    {
        const IndexType index = i * num_dofs_per_node;
        for(IndexType j = 0; j < primal_solution_variable_list.size(); ++j)
            mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(primal_solution_variable_list[j]) = initial_state_variables[index + j];
    }

    KRATOS_CATCH("")
}

int AdjointFiniteDifferenceSmallDisplacementElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    Element::Check(rCurrentProcessInfo);

    const SizeType number_of_nodes = this->GetGeometry().size();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    // Verify that the variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(DENSITY)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( IndexType i = 0; i < number_of_nodes; i++ ) {
        const NodeType &rnode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)

        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z, rnode)
    }

    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(this->GetProperties().Has( CONSTITUTIVE_LAW )) << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;

    // Verify that the constitutive law has the correct dimension
    const SizeType strain_size = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
    if ( dimension == 2 ) {
        KRATOS_ERROR_IF( strain_size < 3 || strain_size > 4) << "Wrong constitutive law used. This is a 2D element! expected strain size is 3 or 4 (el id = ) " << this->Id() << std::endl;
    } else {
        KRATOS_ERROR_IF_NOT(strain_size == 6) << "Wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) "<<  this->Id() << std::endl;
    }

    return 0;

    KRATOS_CATCH( "" );
}

void AdjointFiniteDifferenceSmallDisplacementElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, AdjointFiniteDifferencingBaseElement);
}

void AdjointFiniteDifferenceSmallDisplacementElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, AdjointFiniteDifferencingBaseElement);
}

} // namespace Kratos.


