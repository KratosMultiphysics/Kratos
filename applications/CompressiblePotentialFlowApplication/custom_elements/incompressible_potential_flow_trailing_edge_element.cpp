//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez and Riccardo Rossi
//

#include "incompressible_potential_flow_trailing_edge_element.h"

namespace Kratos
{
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <int Dim, int NumNodes>
Element::Pointer IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Element::Pointer(Kratos::make_shared<IncompressiblePotentialFlowTrailingEdgeElement>(
        NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Element::Pointer(Kratos::make_shared<IncompressiblePotentialFlowTrailingEdgeElement>(
        NewId, pGeom, pProperties));
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::Clone(
    IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Element::Pointer(Kratos::make_shared<IncompressiblePotentialFlowTrailingEdgeElement>(
        NewId, GetGeometry().Create(ThisNodes), pGetProperties()));
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    // Note that the lhs and rhs have double the size
    if (rLeftHandSideMatrix.size1() != 2 * NumNodes ||
        rLeftHandSideMatrix.size2() != 2 * NumNodes)
        rLeftHandSideMatrix.resize(2 * NumNodes, 2 * NumNodes, false);
    if (rRightHandSideVector.size() != 2 * NumNodes)
        rRightHandSideVector.resize(2 * NumNodes, false);
    rLeftHandSideMatrix.clear();
    rRightHandSideVector.clear();

    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    GetWakeDistances(data.distances);

    Matrix lhs_total = ZeroMatrix(NumNodes, NumNodes);

    ComputeLHSGaussPointContribution(data.vol, lhs_total, data);

    Matrix lhs_positive = ZeroMatrix(NumNodes, NumNodes);
    Matrix lhs_negative = ZeroMatrix(NumNodes, NumNodes);

    CalculateLocalSystemSubdividedElement(lhs_positive, lhs_negative);
    AssignLocalSystemSubdividedElement(rLeftHandSideMatrix, lhs_positive,
                                       lhs_negative, lhs_total, data);

    Vector split_element_values(2*NumNodes);
    GetPotential(split_element_values,data.distances);
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, split_element_values);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <int Dim, int NumNodes>
int IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Generic geometry check
    int out = Element::Check(rCurrentProcessInfo);
    if (out != 0)
    {
        return out;
    }

    KRATOS_ERROR_IF(GetGeometry().Area() <= 0.0)
        << this->Id() << "Area cannot be less than or equal to 0" << std::endl;

    for (unsigned int i = 0; i < this->GetGeometry().size(); i++)
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY_POTENTIAL,this->GetGeometry()[i]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(AUXILIARY_VELOCITY_POTENTIAL,this->GetGeometry()[i]);
    }

    return out;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::GetValueOnIntegrationPoints(
    const Variable<bool>& rVariable, std::vector<bool>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);
    if (rVariable == TRAILING_EDGE)
        rValues[0] = this->GetValue(TRAILING_EDGE);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <int Dim, int NumNodes>
std::string IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "IncompressiblePotentialFlowTrailingEdgeElement #" << Id();
    return buffer.str();
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "IncompressiblePotentialFlowTrailingEdgeElement #" << Id();
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::PrintData(std::ostream& rOStream) const
{
    pGetGeometry()->PrintData(rOStream);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::CalculateLocalSystemSubdividedElement(
    Matrix& lhs_positive, Matrix& lhs_negative)
{
    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    GetWakeDistances(data.distances);

    // Subdivide the element
    constexpr unsigned int nvolumes = 3 * (Dim - 1);
    BoundedMatrix<double, NumNodes, Dim> Points;
    array_1d<double, nvolumes> PartitionsSign;
    BoundedMatrix<double, nvolumes, NumNodes> GPShapeFunctionValues;
    array_1d<double, nvolumes> Volumes;
    std::vector<Matrix> GradientsValue(nvolumes);
    BoundedMatrix<double, nvolumes, 2> NEnriched;
    for (unsigned int i = 0; i < GradientsValue.size(); ++i)
        GradientsValue[i].resize(2, Dim, false);
    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        const array_1d<double, 3>& coords = GetGeometry()[i].Coordinates();
        for (unsigned int k = 0; k < Dim; ++k)
        {
            Points(i, k) = coords[k];
        }
    }

    const unsigned int nsubdivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(
        Points, data.DN_DX, data.distances, Volumes, GPShapeFunctionValues,
        PartitionsSign, GradientsValue, NEnriched);

    // Compute the lhs and rhs that would correspond to it being divided
    for (unsigned int i = 0; i < nsubdivisions; ++i)
    {
        if (PartitionsSign[i] > 0)
            ComputeLHSGaussPointContribution(Volumes[i], lhs_positive, data);
        else
            ComputeLHSGaussPointContribution(Volumes[i], lhs_negative, data);
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::AssignLocalSystemSubdividedElement(
    MatrixType& rLeftHandSideMatrix,
    Matrix& lhs_positive,
    Matrix& lhs_negative,
    Matrix& lhs_total,
    const ElementalData<NumNodes, Dim>& data) const
{
    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        // The TE node takes the contribution of the subdivided element and
        // we do not apply the wake condition on the TE node
        if (GetGeometry()[i].GetValue(TRAILING_EDGE))
        {
            for (unsigned int j = 0; j < NumNodes; ++j)
            {
                rLeftHandSideMatrix(i, j) = lhs_positive(i, j);
                rLeftHandSideMatrix(i + NumNodes, j + NumNodes) = lhs_negative(i, j);
            }
        }
        else
            AssignLocalSystemWakeNode(rLeftHandSideMatrix, lhs_total, data, i);
    }
}

// serializer

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class IncompressiblePotentialFlowTrailingEdgeElement<2, 3>;

} // namespace Kratos
