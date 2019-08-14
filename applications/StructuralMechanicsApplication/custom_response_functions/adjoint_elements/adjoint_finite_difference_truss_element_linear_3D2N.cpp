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
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::Calculate(const Variable<Vector >& rVariable, Vector& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // The particular soltuion of the influence function is rotated in global direction
    if (rVariable == ADJOINT_PARTICULAR_DISPLACEMENT)
    {
        static constexpr int number_of_nodes = 2;
        static constexpr int dimension = 3;
        static constexpr unsigned int element_size = number_of_nodes * dimension;

        KRATOS_ERROR_IF(rOutput.size() != element_size) << "Size of particular solution does not fit!" << std::endl;

        TrussElement3D2N::Pointer p_primal_beam_element = dynamic_pointer_cast<TrussElement3D2N>(this->pGetPrimalElement());
        BoundedMatrix<double, element_size, element_size> transformation_matrix = ZeroMatrix(element_size, element_size);
        p_primal_beam_element->CreateTransformationMatrix(transformation_matrix);

        rOutput = prod(transformation_matrix, rOutput);
    }

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
					      std::vector< array_1d<double, 3 > >& rOutput,
					      const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == ADJOINT_STRAIN)
    {
        std::vector<Vector> strain_vector;
        AdjointFiniteDifferencingBaseElement<TPrimalElement>::CalculateAdjointFieldOnIntegrationPoints(STRAIN, strain_vector, rCurrentProcessInfo);
        if (rOutput.size() != strain_vector.size())
            rOutput.resize(strain_vector.size());

        KRATOS_ERROR_IF(strain_vector[0].size() != 3) << "Dimension of strain vector not as expected!" << std::endl;

        double response_value = 1.0;
        if (rCurrentProcessInfo.Has(RESPONSE_VALUE)) {response_value = rCurrentProcessInfo.GetValue(RESPONSE_VALUE);}

        for(IndexType i = 0; i < strain_vector.size(); ++i)
            for (IndexType j = 0; j < 3 ; ++j)
                rOutput[i][j] = strain_vector[i][j] / response_value;
    }
    else if (rVariable == PSEUDO_FORCE)
    {
        if(this->Has(INFLUENCE_FUNCTIONS_EXTENSIONS))
        {
            GeneralizedInfluenceFunctionsExtension my_extension = *(this->GetValue(INFLUENCE_FUNCTIONS_EXTENSIONS));
            my_extension.CalculatePseudoQuantityOnIntegrationPoints(*this->mpPrimalElement, rVariable, rOutput, rCurrentProcessInfo);
        }
        else
            KRATOS_ERROR << "'GeneralizedInfluenceFunctionsExtension' is necessary to compute "<< rVariable.Name() << "!" << std::endl;
        /*{
            const SizeType  write_points_number = GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
            if (rOutput.size() != write_points_number)
                rOutput.resize(write_points_number);
            for(IndexType i = 0; i < write_points_number; ++i)
                rOutput[i].clear();
        }*/
    }
    else
        this->CalculateAdjointFieldOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);

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


