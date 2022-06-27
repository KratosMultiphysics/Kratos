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
#include "custom_utilities/structural_mechanics_element_utilities.h"


namespace Kratos
{

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType num_GP = (this->mpPrimalElement->GetGeometry().IntegrationPoints()).size();
    if (rOutput.size() != num_GP) {
        rOutput.resize(num_GP);
    }

    if(rVariable == CROSS_AREA_VARIATIONAL_SENSITIVITY) {
        std::vector< array_1d<double, 3 > > pseudo_force;
        this->CalculateOnIntegrationPoints(CROSS_AREA_PSEUDO_FORCE, pseudo_force, rCurrentProcessInfo);
        std::vector< array_1d<double, 3 > > adjoint_strain;
        this->CalculateOnIntegrationPoints(ADJOINT_STRAIN, adjoint_strain, rCurrentProcessInfo);
        const double l_0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
        rOutput[0] = -1 * pseudo_force[0][0] * adjoint_strain[0][0] * l_0;
    }

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
					      std::vector< array_1d<double, 3 > >& rOutput,
					      const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == CROSS_AREA_PSEUDO_FORCE) {
        std::vector<Vector> strain_vector;
        this->mpPrimalElement->CalculateOnIntegrationPoints(STRAIN, strain_vector, rCurrentProcessInfo);
        if (rOutput.size() != strain_vector.size()) {
            rOutput.resize(strain_vector.size());
        }
        KRATOS_ERROR_IF(strain_vector[0].size() != 3) << "Dimension of strain vector not as expected!" << std::endl;

        double prestress = 0.00;
        if (this->mpPrimalElement->GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
            prestress = this->mpPrimalElement->GetProperties()[TRUSS_PRESTRESS_PK2];
        }
        const double youngs_modulus = this->mpPrimalElement->GetProperties()[YOUNG_MODULUS];

        double pseudo_force = youngs_modulus * strain_vector[0][0] + prestress;

        rOutput[0][0] = pseudo_force;

    } else if (rVariable == ADJOINT_STRAIN) {
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


