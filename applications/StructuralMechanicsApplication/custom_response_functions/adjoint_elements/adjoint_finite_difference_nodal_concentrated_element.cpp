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


#include "adjoint_finite_difference_nodal_concentrated_element.h"
#include "custom_elements/nodal_concentrated_element.hpp"
#include "structural_mechanics_application_variables.h"
#include "includes/checks.h"
#include "utilities/indirect_scalar.h"

namespace Kratos
{
template <class TPrimalElement>
AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::ThisExtensions::ThisExtensions(Element* pElement)
    : mpElement{pElement}
{
}

template <class TPrimalElement>
void AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::ThisExtensions::GetFirstDerivativesVector(
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
void AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::ThisExtensions::GetSecondDerivativesVector(
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
void AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::ThisExtensions::GetAuxiliaryVector(
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
void AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::ThisExtensions::GetFirstDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    if (rVariables.size() != 1)
    {
        rVariables.resize(1);
    }
    rVariables[0] = &ADJOINT_VECTOR_2;
}

template <class TPrimalElement>
void AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::ThisExtensions::GetSecondDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    if (rVariables.size() != 1)
    {
        rVariables.resize(1);
    }
    rVariables[0] = &ADJOINT_VECTOR_3;
}

template <class TPrimalElement>
void AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::ThisExtensions::GetAuxiliaryVariables(
    std::vector<VariableData const*>& rVariables) const
{
    if (rVariables.size() != 1)
    {
        rVariables.resize(1);
    }
    rVariables[0] = &AUX_ADJOINT_VECTOR_1;
}

template <class TPrimalElement>
void AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::EquationIdVector(EquationIdVectorType& rResult,
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

        rResult[index    ] = iNode.GetDof(ADJOINT_VECTOR_1_X).EquationId();
        rResult[index + 1] = iNode.GetDof(ADJOINT_VECTOR_1_Y).EquationId();
        rResult[index + 2] = iNode.GetDof(ADJOINT_VECTOR_1_Z).EquationId();
    }
    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::GetDofList(DofsVectorType& rElementalDofList,
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

        rElementalDofList[index    ] = iNode.pGetDof(ADJOINT_VECTOR_1_X);
        rElementalDofList[index + 1] = iNode.pGetDof(ADJOINT_VECTOR_1_Y);
        rElementalDofList[index + 2] = iNode.pGetDof(ADJOINT_VECTOR_1_Z);
    }

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::GetValuesVector(Vector& rValues, int Step) const
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

        rValues[index    ] = iNode.FastGetSolutionStepValue(ADJOINT_VECTOR_1_X, Step);
        rValues[index + 1] = iNode.FastGetSolutionStepValue(ADJOINT_VECTOR_1_Y, Step);
        rValues[index + 2] = iNode.FastGetSolutionStepValue(ADJOINT_VECTOR_1_Z, Step);
    }
    
    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::Initialize(rCurrentProcessInfo);
    this->SetValue(ADJOINT_EXTENSIONS, Kratos::make_shared<ThisExtensions>(this));
}

template <class TPrimalElement>
void AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::GetFirstDerivativesVector(Vector& rValues, int Step) const
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
void AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::GetSecondDerivativesVector(Vector& rValues, int Step) const
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
void AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                                ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    BaseType::CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    noalias(rLeftHandSideMatrix) = -rLeftHandSideMatrix;
    KRATOS_CATCH("");
}

template <class TPrimalElement>
void AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::CalculateFirstDerivativesLHS(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    BaseType::CalculateDampingMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
    noalias(rLeftHandSideMatrix) = -rLeftHandSideMatrix;
    KRATOS_CATCH("");
}

template <class TPrimalElement>
void AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    BaseType::CalculateMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
    noalias(rLeftHandSideMatrix) = -rLeftHandSideMatrix;
    KRATOS_CATCH("");
}

template <class TPrimalElement>
void AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::CalculateSensitivityMatrix(
                                            const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    if (rDesignVariable == SCALAR_SENSITIVITY)
    {
        const SizeType number_of_nodes = this->GetGeometry().PointsNumber();
        const SizeType local_size = number_of_nodes * 1.0;
        rOutput = ZeroMatrix(1, local_size);
        const double& x = this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
        rOutput(0, 2) = -x;
    }
    else
    {
        KRATOS_ERROR << "Invalid variable: " << rDesignVariable << std::endl;
    }
    KRATOS_CATCH("");
}

template <class TPrimalElement>
void AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::CalculateSensitivityMatrix(
                                            const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType number_of_nodes = this->GetGeometry().PointsNumber();
    const SizeType dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
    const SizeType num_dofs_per_node = dimension;
    const SizeType local_size = number_of_nodes * num_dofs_per_node;

    if(rDesignVariable == SHAPE_SENSITIVITY) {
        if ((rOutput.size1() != dimension*number_of_nodes) || (rOutput.size2() != local_size)) {
                rOutput.resize(dimension*number_of_nodes, local_size, false);
        }
        noalias(rOutput) = ZeroMatrix(dimension*number_of_nodes, local_size);
    }
    else if (this->Has(rDesignVariable) && rDesignVariable == NODAL_DISPLACEMENT_STIFFNESS) {
        if ((rOutput.size1() != dimension) || (rOutput.size2() != local_size)) {
                rOutput.resize(dimension, local_size, false);
        }

        // save original stiffness parameters
        const auto variable_value = this->pGetPrimalElement()->GetValue(rDesignVariable);

        // reset original stiffness parameters before computing the derivatives
        this->pGetPrimalElement()->SetValue(rDesignVariable, rDesignVariable.Zero());

        ProcessInfo process_info = rCurrentProcessInfo;
        Vector RHS;
        for(IndexType dir_i = 0; dir_i < dimension; ++dir_i) {
            // The following approach assumes a linear dependency between RHS and spring stiffness
            array_1d<double, 3> perturbed_nodal_stiffness = ZeroVector(3);
            perturbed_nodal_stiffness[dir_i] = 1.0;
            this->pGetPrimalElement()->SetValue(rDesignVariable, perturbed_nodal_stiffness);
            this->pGetPrimalElement()->CalculateRightHandSide(RHS, process_info);

            KRATOS_ERROR_IF_NOT(RHS.size() == local_size) << "Size of the pseudo-load does not fit!" << std::endl;
            for(IndexType i = 0; i < RHS.size(); ++i) {
                rOutput(dir_i, i) = RHS[i];
            }
        }

        // give original stiffness parameters back
        this->pGetPrimalElement()->SetValue(rDesignVariable, variable_value);
    }
    else {
        if ((rOutput.size1() != 0) || (rOutput.size2() != local_size)) {
            rOutput.resize(0, local_size, false);
        }
        noalias(rOutput) = ZeroMatrix(0, local_size);
    }

    KRATOS_CATCH("")
}

template <class TPrimalElement>
int AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(this->pGetPrimalElement()) << "Primal element pointer is nullptr!" << std::endl;

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_MASS)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_DISPLACEMENT_STIFFNESS)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_DAMPING_RATIO)
    KRATOS_CHECK_VARIABLE_KEY(VOLUME_ACCELERATION)
    /*
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_DISPLACEMENT)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( std::size_t i = 0; i < this->GetGeometry().size(); ++i ) {
        const Node<3>& rnode = this->GetGeometry()[i];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT,rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VOLUME_ACCELERATION,rnode)

        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X,rnode)
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y,rnode)
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z,rnode)
        
    }*/

    return 0;

    KRATOS_CATCH( "Problem in the Check in the AdjointFiniteDifferenceNodalConcentratedElement!" )
}

template <class TPrimalElement>
void AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TPrimalElement>
void AdjointFiniteDifferenceNodalConcentratedElement<TPrimalElement>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

template class AdjointFiniteDifferenceNodalConcentratedElement<NodalConcentratedElement>;

} // namespace Kratos.


