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
#include "utilities/indirect_scalar.h"


namespace Kratos
{

template <class TPrimalElement>
AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::ThisExtensions::ThisExtensions(Element* pElement)
    : mpElement{pElement}
{
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::ThisExtensions::GetFirstDerivativesVector(
    std::size_t NodeId, std::vector<IndirectScalar<double>>& rVector, std::size_t Step)
{
    auto& r_node = mpElement->GetGeometry()[NodeId];
    rVector.resize(mpElement->GetGeometry().WorkingSpaceDimension());
    std::size_t index = 0;
    rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_VECTOR_2_X, Step);
    rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_VECTOR_2_Y, Step);
    if (mpElement->GetGeometry().WorkingSpaceDimension() == 3)
    {
        rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_VECTOR_2_Z, Step);
    }

}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::ThisExtensions::GetSecondDerivativesVector(
    std::size_t NodeId, std::vector<IndirectScalar<double>>& rVector, std::size_t Step)
{
    auto& r_node = mpElement->GetGeometry()[NodeId];
    rVector.resize(mpElement->GetGeometry().WorkingSpaceDimension());
    std::size_t index = 0;
    rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_VECTOR_3_X, Step);
    rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_VECTOR_3_Y, Step);
    if (mpElement->GetGeometry().WorkingSpaceDimension() == 3)
    {
        rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_VECTOR_3_Z, Step);
    }
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::ThisExtensions::GetAuxiliaryVector(
    std::size_t NodeId, std::vector<IndirectScalar<double>>& rVector, std::size_t Step)
{
    auto& r_node = mpElement->GetGeometry()[NodeId];
    rVector.resize(mpElement->GetGeometry().WorkingSpaceDimension());
    std::size_t index = 0;
    rVector[index++] = MakeIndirectScalar(r_node, AUX_ADJOINT_VECTOR_1_X, Step);
    rVector[index++] = MakeIndirectScalar(r_node, AUX_ADJOINT_VECTOR_1_Y, Step);
    if (mpElement->GetGeometry().WorkingSpaceDimension() == 3)
    {
        rVector[index++] = MakeIndirectScalar(r_node, AUX_ADJOINT_VECTOR_1_Z, Step);
    }
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::ThisExtensions::GetFirstDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    if (rVariables.size() != 1)
    {
        rVariables.resize(1);
    }
    rVariables[0] = &ADJOINT_VECTOR_2;
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::ThisExtensions::GetSecondDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    if (rVariables.size() != 1)
    {
        rVariables.resize(1);
    }
    rVariables[0] = &ADJOINT_VECTOR_3;
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::ThisExtensions::GetAuxiliaryVariables(
    std::vector<VariableData const*>& rVariables) const
{
    if (rVariables.size() != 1)
    {
        rVariables.resize(1);
    }
    rVariables[0] = &AUX_ADJOINT_VECTOR_1;
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::EquationIdVector(EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{  
    KRATOS_TRY
    const GeometryType& geom = this->GetGeometry();

    const SizeType number_of_nodes = geom.PointsNumber();
    const SizeType dimension = geom.WorkingSpaceDimension();
    const SizeType num_dofs = number_of_nodes * dimension;

    if(rResult.size() != num_dofs)
        rResult.resize(num_dofs, false);

    for(IndexType i = 0; i < geom.size(); ++i)
    {
        const IndexType index = i * dimension;
        const NodeType& iNode = geom[i];

        rResult[index    ] = iNode.GetDof(ADJOINT_DISPLACEMENT_X).EquationId();
        rResult[index + 1] = iNode.GetDof(ADJOINT_DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = iNode.GetDof(ADJOINT_DISPLACEMENT_Z).EquationId();
    }
    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::GetDofList(DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType & geom = this->GetGeometry();

    const SizeType number_of_nodes = geom.PointsNumber();
    const SizeType dimension =  geom.WorkingSpaceDimension();
    const SizeType num_dofs = number_of_nodes * dimension;

    if (rElementalDofList.size() != num_dofs)
        rElementalDofList.resize(num_dofs);

    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        const IndexType index = i * dimension;
        const NodeType& iNode = geom[i];

        rElementalDofList[index    ] = iNode.pGetDof(ADJOINT_DISPLACEMENT_X);
        rElementalDofList[index + 1] = iNode.pGetDof(ADJOINT_DISPLACEMENT_Y);
        rElementalDofList[index + 2] = iNode.pGetDof(ADJOINT_DISPLACEMENT_Z);
    }

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::GetValuesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY

    const GeometryType & geom = this->GetGeometry();

    const SizeType number_of_nodes = geom.PointsNumber();
    const SizeType dimension =  geom.WorkingSpaceDimension();
    const SizeType num_dofs = number_of_nodes * dimension;

    if(rValues.size() != num_dofs)
        rValues.resize(num_dofs, false);

    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        const IndexType index = i * dimension;
        const NodeType& iNode = geom[i];

        rValues[index    ] = iNode.FastGetSolutionStepValue(ADJOINT_DISPLACEMENT_X, Step);
        rValues[index + 1] = iNode.FastGetSolutionStepValue(ADJOINT_DISPLACEMENT_Y, Step);
        rValues[index + 2] = iNode.FastGetSolutionStepValue(ADJOINT_DISPLACEMENT_Z, Step);
    }
    
    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::Initialize(rCurrentProcessInfo);
    this->SetValue(ADJOINT_EXTENSIONS, Kratos::make_shared<ThisExtensions>(this));
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                                ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    BaseType::CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    noalias(rLeftHandSideMatrix) = -rLeftHandSideMatrix;
    KRATOS_CATCH("");
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY

    const GeometryType & geom = this->GetGeometry();

    const SizeType number_of_nodes = geom.PointsNumber();
    const SizeType dimension =  geom.WorkingSpaceDimension();
    const SizeType num_dofs = number_of_nodes * dimension;

    if(rValues.size() != num_dofs)
        rValues.resize(num_dofs, false);

    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        const IndexType index = i * dimension;
        const NodeType & iNode = geom[i];

        rValues[index    ] = iNode.FastGetSolutionStepValue(ADJOINT_VECTOR_2_X, Step);
        rValues[index + 1] = iNode.FastGetSolutionStepValue(ADJOINT_VECTOR_2_Y, Step);
        rValues[index + 2] = iNode.FastGetSolutionStepValue(ADJOINT_VECTOR_2_Z, Step);
    }

    KRATOS_CATCH("")  
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY

    const GeometryType & geom = this->GetGeometry();

    const SizeType number_of_nodes = geom.PointsNumber();
    const SizeType dimension =  geom.WorkingSpaceDimension();
    const SizeType num_dofs = number_of_nodes * dimension;

    if(rValues.size() != num_dofs)
        rValues.resize(num_dofs, false);

    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        const IndexType index = i * dimension;
        const NodeType & iNode = geom[i];

        rValues[index    ] = iNode.FastGetSolutionStepValue(ADJOINT_VECTOR_3_X, Step);
        rValues[index + 1] = iNode.FastGetSolutionStepValue(ADJOINT_VECTOR_3_Y, Step);
        rValues[index + 2] = iNode.FastGetSolutionStepValue(ADJOINT_VECTOR_3_Z, Step);
    }

    KRATOS_CATCH("")  
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::CalculateFirstDerivativesLHS(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    BaseType::CalculateDampingMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
    noalias(rLeftHandSideMatrix) = -rLeftHandSideMatrix;
    KRATOS_CATCH("");
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElementLinear<TPrimalElement>::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    BaseType::CalculateMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
    noalias(rLeftHandSideMatrix) = -rLeftHandSideMatrix;
    KRATOS_CATCH("");
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


