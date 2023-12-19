// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//


#include "adjoint_finite_difference_truss_element_linear_3D2N.h"
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/response_utilities/stress_response_definitions.h"
#include "custom_elements/truss_element_linear_3D2N.hpp"


namespace Kratos
{

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

        auto& r_geometry = this->pGetPrimalElement()->GetGeometry();
        const SizeType number_of_nodes = r_geometry.PointsNumber();
        const SizeType dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
        const SizeType local_size = number_of_nodes * dimension;
        
    if ( this->pGetPrimalElement()->GetProperties().Has(rDesignVariable) )
    {
        KRATOS_WATCH(rDesignVariable.Name())

        if ( (rOutput.size1() != 1) || (rOutput.size2() != local_size ) )
              rOutput.resize(1, local_size, false);

        Vector DispalacementVector = ZeroVector(local_size);
        Vector VelocityVector = ZeroVector(local_size);
        Vector AccelerationVector = ZeroVector(local_size);
        Matrix DampingMatrix;
        Matrix MassMatrix;
        Vector RHS;
        Vector DampingForce;
        Vector InertiaForce;
    
        for(IndexType n_i = 0; n_i < number_of_nodes; ++n_i) {
            for(IndexType dir_i = 0; dir_i < dimension; ++dir_i) {
                DispalacementVector[dir_i + n_i * dimension] = r_geometry[n_i].FastGetSolutionStepValue(DISPLACEMENT)[dir_i];
                VelocityVector[dir_i + n_i * dimension] = r_geometry[n_i].FastGetSolutionStepValue(VELOCITY)[dir_i];
                AccelerationVector[dir_i + n_i * dimension] = r_geometry[n_i].FastGetSolutionStepValue(ACCELERATION)[dir_i];
            }
        }
    
        this->pGetPrimalElement()->CalculateDampingMatrix(DampingMatrix, rCurrentProcessInfo);
        this->pGetPrimalElement()->CalculateMassMatrix(MassMatrix, rCurrentProcessInfo);

        this->pGetPrimalElement()->CalculateRightHandSide(RHS, rCurrentProcessInfo);
        DampingForce = -prod(DampingMatrix, VelocityVector);
        InertiaForce = -prod(MassMatrix, AccelerationVector);

        KRATOS_WATCH(DispalacementVector)
        KRATOS_WATCH(VelocityVector)
        KRATOS_WATCH(AccelerationVector)
        KRATOS_WATCH(MassMatrix)
        KRATOS_WATCH(DampingMatrix)
        Matrix StiffnessMatrix;
        this->pGetPrimalElement()->CalculateMassMatrix(StiffnessMatrix, rCurrentProcessInfo);
        KRATOS_WATCH(StiffnessMatrix)
        KRATOS_WATCH(RHS)
        KRATOS_WATCH(InertiaForce)

        // Save property pointer
        Properties::Pointer p_global_properties = this->pGetPrimalElement()->pGetProperties();

        // Create new  shared property pointer and assign it to the element
        Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
        this->pGetPrimalElement()->SetProperties(p_local_property);
        
        // Store current property value, calculate perturbation size and perturb the design variable
        const double current_property_value = this->pGetPrimalElement()->GetProperties()[rDesignVariable];
        double delta = rCurrentProcessInfo[PERTURBATION_SIZE] * current_property_value;
        KRATOS_DEBUG_ERROR_IF_NOT(delta > 0) << "The perturbation size is not > 0!";
        p_local_property->SetValue(rDesignVariable, (current_property_value + delta));

        KRATOS_WATCH(current_property_value)
        KRATOS_WATCH(delta)

        // Calculate RHS derivative
        Vector RHS_perturbed;
        this->pGetPrimalElement()->CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
        for(IndexType i = 0; i < RHS_perturbed.size(); ++i)
            rOutput(0, i) = (RHS_perturbed[i] - RHS[i]) / delta;

        KRATOS_WATCH(rOutput)

        // Calculate damping force derivative
        Matrix DampingMatrix_perturbed;
        Vector DampingForce_perturbed;
        this->pGetPrimalElement()->CalculateDampingMatrix(DampingMatrix_perturbed, rCurrentProcessInfo);
        DampingForce_perturbed = -prod(DampingMatrix_perturbed, VelocityVector);
        for(IndexType i = 0; i < DampingForce_perturbed.size(); ++i)
            rOutput(0, i) += (DampingForce_perturbed[i] - DampingForce[i]) / delta;

        KRATOS_WATCH(rOutput)

        // Calculate inertia force derivative
        Matrix MassMatrix_perturbed;
        Vector InertiaForce_perturbed;
        this->pGetPrimalElement()->CalculateMassMatrix(MassMatrix_perturbed, rCurrentProcessInfo);
        InertiaForce_perturbed = -prod(MassMatrix_perturbed, AccelerationVector);
        for(IndexType i = 0; i < InertiaForce_perturbed.size(); ++i)
            rOutput(0, i) += (InertiaForce_perturbed[i] - InertiaForce[i]) / delta;

        KRATOS_WATCH(rOutput)

        // Give element original properties back
        this->pGetPrimalElement()->SetProperties(p_global_properties);
    }
    else
    {
        rOutput = ZeroMatrix(0, local_size);
    }

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
					      std::vector< array_1d<double, 3 > >& rOutput,
					      const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == ADJOINT_STRAIN) {
        std::vector<Vector> strain_vector;
        this->CalculateAdjointFieldOnIntegrationPoints(STRAIN, strain_vector, rCurrentProcessInfo);
        if (rOutput.size() != strain_vector.size()) {
            rOutput.resize(strain_vector.size());
        }

        KRATOS_ERROR_IF(strain_vector[0].size() != 3) << "Dimension of strain vector not as expected!" << std::endl;

        for(IndexType i = 0; i < strain_vector.size(); ++i) {
            for (IndexType j = 0; j < 3 ; ++j) {
                rOutput[i][j] = strain_vector[i][j];
            }
        }
    } else {
        this->CalculateAdjointFieldOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
}


template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));

    if(traced_stress_type  == TracedStressType::FX)
    {
        // ensure that adjoint load is determined without influence of pre-stress
        // pre-stress does not cancel out when computing this derivative with unit-displacements!
        Properties::Pointer p_global_properties = this->mpPrimalElement->pGetProperties();

        Properties::Pointer p_local_property(Kratos::make_shared<Properties>(*p_global_properties));
        this->mpPrimalElement->SetProperties(p_local_property);

        p_local_property->SetValue(TRUSS_PRESTRESS_PK2, 0.0);

        AdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateStressDisplacementDerivative(rStressVariable,
                                           rOutput, rCurrentProcessInfo);

        this->mpPrimalElement->SetProperties(p_global_properties);
    }
    else
        AdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateStressDisplacementDerivative(rStressVariable,
                                   rOutput, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

template class AdjointFiniteDifferenceTrussElementLinear<TrussElementLinear3D2N>;

} // namespace Kratos.


