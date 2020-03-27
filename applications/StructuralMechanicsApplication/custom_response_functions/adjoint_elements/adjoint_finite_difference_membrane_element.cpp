// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//

// System includes

// External includes

// Project includes
#include "adjoint_finite_difference_membrane_element.h"
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/response_utilities/stress_response_definitions.h"
#include "includes/checks.h"
#include "custom_elements/membrane_element.hpp"


namespace Kratos
{

template <class TPrimalElement>
void AdjointFiniteDifferencingMembraneElement<TPrimalElement>::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if((rDesignVariable == PRESTRESS_VECTOR_X_SENSITIVITY || rDesignVariable == PRESTRESS_VECTOR_Y_SENSITIVITY) && this->GetProperties().Has(PRESTRESS_VECTOR)) {
        // define working variables
        ProcessInfo process_info = rCurrentProcessInfo;
        Vector RHS;
        Vector RHS_perturbed;
        int component_index = 0;
        if (rDesignVariable == PRESTRESS_VECTOR_Y_SENSITIVITY) {
            component_index = 1;
        }

        // compute unperturbed RHS
        this->pGetPrimalElement()->CalculateRightHandSide(RHS, process_info);
  
        // Save property pointer
        Properties::Pointer p_global_properties = this->pGetPrimalElement()->pGetProperties();

        // Create new property and assign it to the primal element
        Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
        this->pGetPrimalElement()->SetProperties(p_local_property);

        auto pre_stress = this->pGetPrimalElement()->GetProperties()[PRESTRESS_VECTOR];

        // Get perturbation size
        double delta = rCurrentProcessInfo[PERTURBATION_SIZE];
        if (rCurrentProcessInfo[ADAPT_PERTURBATION_SIZE]) {
            delta *= pre_stress[component_index];
        }
        KRATOS_DEBUG_ERROR_IF_NOT(delta > 0) << "The perturbation size is not > 0!";

        // perturb the design variable
        pre_stress[component_index] += delta;
        p_local_property->SetValue(PRESTRESS_VECTOR, pre_stress);

        // Compute RHS after perturbation
        this->pGetPrimalElement()->CalculateRightHandSide(RHS_perturbed, process_info);

        if( (rOutput.size1() != 1) || (rOutput.size2() != RHS.size() ) ) {
            rOutput.resize(1, RHS.size(), false);
        }

        // Compute derivative of RHS w.r.t. design variable with finite differences
        row(rOutput, 0) = (RHS_perturbed - RHS) / delta;
            
        // Give element original properties back
        this->pGetPrimalElement()->SetProperties(p_global_properties);
    }
    else {
        BaseType::CalculateSensitivityMatrix(rDesignVariable, rOutput, rCurrentProcessInfo);
    }
    
    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferencingMembraneElement<TPrimalElement>::CalculateStressDesignVariableDerivative(const Variable<double>& rDesignVariable,
                                                const Variable<Vector>& rStressVariable, Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if((rDesignVariable == PRESTRESS_VECTOR_X_SENSITIVITY || rDesignVariable == PRESTRESS_VECTOR_Y_SENSITIVITY) && this->GetProperties().Has(PRESTRESS_VECTOR)) {
        // Define working variables
        Vector stress_vector_undist;
        Vector stress_vector_dist;
        int component_index = 0;
        if (rDesignVariable == PRESTRESS_VECTOR_Y_SENSITIVITY) {
            component_index = 1;
        }
 
        // Compute stress on GP before perturbation
        TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));
        if (rStressVariable == STRESS_ON_GP) {
            StressCalculation::CalculateStressOnGP(*(this->pGetPrimalElement()), traced_stress_type, stress_vector_undist, rCurrentProcessInfo);
        }
        else {
            StressCalculation::CalculateStressOnNode(*(this->pGetPrimalElement()), traced_stress_type, stress_vector_undist, rCurrentProcessInfo);
        }

        // Save property pointer
        Properties::Pointer p_global_properties = this->pGetPrimalElement()->pGetProperties();

        // Create new property and assign it to the primal element
        Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
        this->pGetPrimalElement()->SetProperties(p_local_property);

        auto pre_stress = this->pGetPrimalElement()->GetProperties()[PRESTRESS_VECTOR];

        // Get perturbation size
        double delta = rCurrentProcessInfo[PERTURBATION_SIZE];
        if (rCurrentProcessInfo[ADAPT_PERTURBATION_SIZE]) {
            delta *= pre_stress[component_index];
        }
        KRATOS_DEBUG_ERROR_IF_NOT(delta > 0) << "The perturbation size is not > 0!";

        // perturb the design variable
        pre_stress[component_index] += delta;
        p_local_property->SetValue(PRESTRESS_VECTOR, pre_stress);

        // Compute stress on GP after perturbation
        if (rStressVariable == STRESS_ON_GP) {
            StressCalculation::CalculateStressOnGP(*(this->pGetPrimalElement()), traced_stress_type, stress_vector_dist, rCurrentProcessInfo);
        }
        else {
            StressCalculation::CalculateStressOnNode(*(this->pGetPrimalElement()), traced_stress_type, stress_vector_dist, rCurrentProcessInfo);
        }

        if( (rOutput.size1() != 1) || (rOutput.size2() != stress_vector_undist.size() ) ) {
            rOutput.resize(1, stress_vector_undist.size(), false);
        }
        // Compute derivative of stress w.r.t. design variable with finite differences
        row(rOutput, 0) = (stress_vector_dist - stress_vector_undist) / delta;

        // Give element original properties back
        this->pGetPrimalElement()->SetProperties(p_global_properties);
    }
    else {
        BaseType::CalculateStressDesignVariableDerivative(rDesignVariable, rStressVariable, rOutput, rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
}

template <class TPrimalElement>
int AdjointFiniteDifferencingMembraneElement<TPrimalElement>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // TODO MFusseder: which checks are already done by the basetype and can be removed??

    int return_value = BaseType::Check(rCurrentProcessInfo);

    const SizeType number_of_nodes = this->GetGeometry().size();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo[DOMAIN_SIZE]==3) << "DOMAIN_SIZE in element " << this->Id() << " is not 3" << std::endl;
    KRATOS_ERROR_IF_NOT(dimension==3) << "dimension in element " << this->Id() << " is not 3" << std::endl;

    // Verify that the variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(DENSITY)
    KRATOS_CHECK_VARIABLE_KEY(VOLUME_ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(THICKNESS)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( SizeType i = 0; i < number_of_nodes; i++ ) {
        const Node<3> &r_node = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT,r_node)

        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y, r_node)
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z, r_node)
    }

    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(this->GetProperties().Has( CONSTITUTIVE_LAW ))
        << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;

    // Verify that the constitutive law has the correct dimension
    const SizeType strain_size = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
    KRATOS_ERROR_IF( strain_size != 3) << "Wrong constitutive law used. This is a membrane element! "
        << "Expected strain size is 3 (el id = " << this->Id() << ")" << std::endl;

    return return_value;

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferencingMembraneElement<TPrimalElement>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TPrimalElement>
void AdjointFiniteDifferencingMembraneElement<TPrimalElement>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType);

}

template class AdjointFiniteDifferencingMembraneElement<MembraneElement>;

} // namespace Kratos

