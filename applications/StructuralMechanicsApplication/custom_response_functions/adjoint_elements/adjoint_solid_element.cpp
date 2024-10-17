// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:
//

// System includes
#include <vector>

// External include

// Project includes
#include "includes/checks.h"
#include "containers/variable_data.h"
#include "utilities/indirect_scalar.h"

// Application includes
#include "custom_elements/total_lagrangian.h"
#include "custom_response_functions/adjoint_elements/adjoint_solid_element.h"

namespace Kratos
{
template <class TPrimalElement>
AdjointSolidElement<TPrimalElement>::ThisExtensions::ThisExtensions(Element* pElement)
    : mpElement{pElement}
{
}

template <class TPrimalElement>
void AdjointSolidElement<TPrimalElement>::ThisExtensions::GetFirstDerivativesVector(
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
void AdjointSolidElement<TPrimalElement>::ThisExtensions::GetSecondDerivativesVector(
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
void AdjointSolidElement<TPrimalElement>::ThisExtensions::GetAuxiliaryVector(
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
void AdjointSolidElement<TPrimalElement>::ThisExtensions::GetFirstDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    if (rVariables.size() != 1)
    {
        rVariables.resize(1);
    }
    rVariables[0] = &ADJOINT_VECTOR_2;
}

template <class TPrimalElement>
void AdjointSolidElement<TPrimalElement>::ThisExtensions::GetSecondDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    if (rVariables.size() != 1)
    {
        rVariables.resize(1);
    }
    rVariables[0] = &ADJOINT_VECTOR_3;
}

template <class TPrimalElement>
void AdjointSolidElement<TPrimalElement>::ThisExtensions::GetAuxiliaryVariables(
    std::vector<VariableData const*>& rVariables) const
{
    if (rVariables.size() != 1)
    {
        rVariables.resize(1);
    }
    rVariables[0] = &AUX_ADJOINT_VECTOR_1;
}

template <class TPrimalElement>
AdjointSolidElement<TPrimalElement>::AdjointSolidElement(IndexType NewId)
    : Element(NewId), mPrimalElement(NewId, pGetGeometry())
{
}

template <class TPrimalElement>
AdjointSolidElement<TPrimalElement>::AdjointSolidElement(IndexType NewId,
                                                         GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry), mPrimalElement(NewId, pGeometry)
{
}

template <class TPrimalElement>
AdjointSolidElement<TPrimalElement>::AdjointSolidElement(IndexType NewId,
                                                         GeometryType::Pointer pGeometry,
                                                         PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties), mPrimalElement(NewId, pGeometry, pProperties)
{
}

template <class TPrimalElement>
Element::Pointer AdjointSolidElement<TPrimalElement>::Create(IndexType NewId,
                                                             NodesArrayType const& ThisNodes,
                                                             PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<AdjointSolidElement<TPrimalElement>>(
        NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template <class TPrimalElement>
Element::Pointer AdjointSolidElement<TPrimalElement>::Create(IndexType NewId,
                                                             GeometryType::Pointer pGeom,
                                                             PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<AdjointSolidElement<TPrimalElement>>(NewId, pGeom, pProperties);
}

template <class TPrimalElement>
void AdjointSolidElement<TPrimalElement>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    mPrimalElement.Initialize(rCurrentProcessInfo);
    this->SetValue(ADJOINT_EXTENSIONS, Kratos::make_shared<ThisExtensions>(this));
    KRATOS_CATCH("");
}

template <class TPrimalElement>
void AdjointSolidElement<TPrimalElement>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    mPrimalElement.InitializeSolutionStep(rCurrentProcessInfo);
    KRATOS_CATCH("");
}

template <class TPrimalElement>
void AdjointSolidElement<TPrimalElement>::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    mPrimalElement.InitializeNonLinearIteration(rCurrentProcessInfo);
    KRATOS_CATCH("");
}

template <class TPrimalElement>
void AdjointSolidElement<TPrimalElement>::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    mPrimalElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
    KRATOS_CATCH("");
}

template <class TPrimalElement>
void AdjointSolidElement<TPrimalElement>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    mPrimalElement.FinalizeSolutionStep(rCurrentProcessInfo);
    KRATOS_CATCH("");
}

template <class TPrimalElement>
void AdjointSolidElement<TPrimalElement>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                                const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    mPrimalElement.CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    noalias(rLeftHandSideMatrix) = -rLeftHandSideMatrix;
    KRATOS_CATCH("");
}

template <class TPrimalElement>
void AdjointSolidElement<TPrimalElement>::CalculateFirstDerivativesLHS(
    MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    mPrimalElement.CalculateDampingMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
    noalias(rLeftHandSideMatrix) = -rLeftHandSideMatrix;
    KRATOS_CATCH("");
}

template <class TPrimalElement>
void AdjointSolidElement<TPrimalElement>::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    mPrimalElement.CalculateMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
    noalias(rLeftHandSideMatrix) = -rLeftHandSideMatrix;
    KRATOS_CATCH("");
}

template <class TPrimalElement>
void AdjointSolidElement<TPrimalElement>::GetValuesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY;
    const auto& r_geom = mPrimalElement.GetGeometry();
    const unsigned dimension = r_geom.WorkingSpaceDimension();
    const unsigned mat_size = r_geom.size() * dimension;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);
    for (unsigned int i = 0; i < r_geom.size(); ++i)
    {
        const array_1d<double, 3>& adjoint_displacement =
            r_geom[i].FastGetSolutionStepValue(ADJOINT_DISPLACEMENT, Step);
        const unsigned index = i * dimension;
        for (unsigned k = 0; k < dimension; ++k)
            rValues[index + k] = adjoint_displacement[k];
    }
    KRATOS_CATCH("");
}

template <class TPrimalElement>
void AdjointSolidElement<TPrimalElement>::EquationIdVector(EquationIdVectorType& rResult,
                                                           const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    auto& r_geom = mPrimalElement.GetGeometry();
    const unsigned number_of_nodes = r_geom.size();
    const unsigned dimension = r_geom.WorkingSpaceDimension();

    if (rResult.size() != dimension * number_of_nodes)
        rResult.resize(dimension * number_of_nodes, false);

    const unsigned pos = r_geom[0].GetDofPosition(ADJOINT_DISPLACEMENT_X);

    if (dimension == 2)
    {
        for (unsigned i = 0; i < number_of_nodes; ++i)
        {
            const unsigned index = i * 2;
            rResult[index] = r_geom[i].GetDof(ADJOINT_DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] =
                r_geom[i].GetDof(ADJOINT_DISPLACEMENT_Y, pos + 1).EquationId();
        }
    }
    else
    {
        for (unsigned i = 0; i < number_of_nodes; ++i)
        {
            const unsigned index = i * 3;
            rResult[index] = r_geom[i].GetDof(ADJOINT_DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] =
                r_geom[i].GetDof(ADJOINT_DISPLACEMENT_Y, pos + 1).EquationId();
            rResult[index + 2] =
                r_geom[i].GetDof(ADJOINT_DISPLACEMENT_Z, pos + 2).EquationId();
        }
    }
    KRATOS_CATCH("");
}

template <class TPrimalElement>
void AdjointSolidElement<TPrimalElement>::GetDofList(DofsVectorType& rElementalDofList,
                                                     const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    auto& r_geom = mPrimalElement.GetGeometry();
    const unsigned number_of_nodes = r_geom.size();
    const unsigned dimension = r_geom.WorkingSpaceDimension();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dimension * number_of_nodes);

    if (dimension == 2)
    {
        for (unsigned i = 0; i < number_of_nodes; ++i)
        {
            rElementalDofList.push_back(r_geom[i].pGetDof(ADJOINT_DISPLACEMENT_X));
            rElementalDofList.push_back(r_geom[i].pGetDof(ADJOINT_DISPLACEMENT_Y));
        }
    }
    else
    {
        for (unsigned int i = 0; i < number_of_nodes; ++i)
        {
            rElementalDofList.push_back(r_geom[i].pGetDof(ADJOINT_DISPLACEMENT_X));
            rElementalDofList.push_back(r_geom[i].pGetDof(ADJOINT_DISPLACEMENT_Y));
            rElementalDofList.push_back(r_geom[i].pGetDof(ADJOINT_DISPLACEMENT_Z));
        }
    }
    KRATOS_CATCH("");
}

template <class TPrimalElement>
int AdjointSolidElement<TPrimalElement>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    for (const auto& r_node : GetGeometry()) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z, r_node);
    }
    return 0;
    KRATOS_CATCH("");
}

template <class TPrimalElement>
void AdjointSolidElement<TPrimalElement>::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rDesignVariable, Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    mPrimalElement.CalculateSensitivityMatrix(rDesignVariable, rOutput, rCurrentProcessInfo);
    KRATOS_CATCH("");
}

template class AdjointSolidElement<TotalLagrangian>;

} // namespace Kratos.
