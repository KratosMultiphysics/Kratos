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
    std::cout << "TE ELEMENT = " << this->Id() << std::endl;
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

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    // TODO: improve speed
    Matrix tmp;
    CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::EquationIdVector(
    EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    if (rResult.size() != 2 * NumNodes)
        rResult.resize(2 * NumNodes, false);

    array_1d<double, NumNodes> distances;
    GetWakeDistances(distances);

    // Positive part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] > 0.0)
            rResult[i] = GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
        else
            rResult[i] = GetGeometry()[i].GetDof(AUXILIARY_VELOCITY_POTENTIAL).EquationId();
    }

    // Negative part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] < 0.0)
            rResult[NumNodes + i] =
                GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
        else
            rResult[NumNodes + i] =
                GetGeometry()[i].GetDof(AUXILIARY_VELOCITY_POTENTIAL).EquationId();
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::GetDofList(
    DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    if (rElementalDofList.size() != 2 * NumNodes)
        rElementalDofList.resize(2 * NumNodes);

    array_1d<double, NumNodes> distances;
    GetWakeDistances(distances);

    // Positive part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] > 0)
            rElementalDofList[i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
        else
            rElementalDofList[i] = GetGeometry()[i].pGetDof(AUXILIARY_VELOCITY_POTENTIAL);
    }

    // Negative part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] < 0)
            rElementalDofList[NumNodes + i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
        else
            rElementalDofList[NumNodes + i] =
                GetGeometry()[i].pGetDof(AUXILIARY_VELOCITY_POTENTIAL);
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    CheckWakeCondition();
    ComputePotentialJump(rCurrentProcessInfo);
    ComputeElementInternalEnergy();
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
    const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == PRESSURE)
    {
        double p = ComputeUpperPressure(rCurrentProcessInfo);
        rValues[0] = p;
    }
    else if (rVariable == PRESSURE_LOWER)
    {
        double p = ComputeLowerPressure(rCurrentProcessInfo);
        rValues[0] = p;
    }
    else if (rVariable == WAKE)
    {
        const IncompressiblePotentialFlowTrailingEdgeElement& r_this = *this;
        rValues[0] = r_this.GetValue(WAKE);
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::GetValueOnIntegrationPoints(
    const Variable<int>& rVariable, std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == TRAILING_EDGE)
        rValues[0] = this->GetValue(TRAILING_EDGE);
    else if (rVariable == WAKE)
        rValues[0] = this->GetValue(WAKE);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);
    if (rVariable == VELOCITY)
    {
        array_1d<double, 3> v(3, 0.0);
        array_1d<double, Dim> vaux;
        ComputeUpperVelocity(vaux);
        for (unsigned int k = 0; k < Dim; k++)
            v[k] = vaux[k];
        rValues[0] = v;
    }
    else if (rVariable == VELOCITY_LOWER)
    {
        array_1d<double, 3> v(3, 0.0);
        array_1d<double, Dim> vaux;
        ComputeLowerVelocity(vaux);
        for (unsigned int k = 0; k < Dim; k++)
            v[k] = vaux[k];
        rValues[0] = v;
    }
}

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
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::GetWakeDistances(array_1d<double, NumNodes>& distances) const
{
    noalias(distances) = GetValue(ELEMENTAL_DISTANCES);
}

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
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::ComputeLHSGaussPointContribution(
    const double weight, Matrix& lhs, const ElementalData<NumNodes, Dim>& data) const
{
    noalias(lhs) += weight * prod(data.DN_DX, trans(data.DN_DX));
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

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::AssignLocalSystemWakeElement(
    MatrixType& rLeftHandSideMatrix, Matrix& lhs_total, const ElementalData<NumNodes, Dim>& data) const
{
    for (unsigned int row = 0; row < NumNodes; ++row)
        AssignLocalSystemWakeNode(rLeftHandSideMatrix, lhs_total, data, row);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::AssignLocalSystemWakeNode(
    MatrixType& rLeftHandSideMatrix,
    Matrix& lhs_total,
    const ElementalData<NumNodes, Dim>& data,
    unsigned int& row) const
{
    // Filling the diagonal blocks (i.e. decoupling upper and lower dofs)
    for (unsigned int column = 0; column < NumNodes; ++column)
    {
        rLeftHandSideMatrix(row, column) = lhs_total(row, column);
        rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = lhs_total(row, column);
    }

    // Applying wake condition on the AUXILIARY_VELOCITY_POTENTIAL dofs
    if (data.distances[row] < 0.0)
        for (unsigned int column = 0; column < NumNodes; ++column)
            rLeftHandSideMatrix(row, column + NumNodes) = -lhs_total(row, column); // Side 1
    else if (data.distances[row] > 0.0)
        for (unsigned int column = 0; column < NumNodes; ++column)
            rLeftHandSideMatrix(row + NumNodes, column) = -lhs_total(row, column); // Side 2
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::CheckWakeCondition() const
{
    array_1d<double, Dim> upper_wake_velocity;
    ComputeUpperVelocity(upper_wake_velocity);
    const double vupnorm = inner_prod(upper_wake_velocity, upper_wake_velocity);

    array_1d<double, Dim> lower_wake_velocity;
    ComputeLowerVelocity(lower_wake_velocity);
    const double vlownorm = inner_prod(lower_wake_velocity, lower_wake_velocity);

    if (std::abs(vupnorm - vlownorm) > 0.1)
        std::cout << "WAKE CONDITION NOT FULFILLED IN ELEMENT # " << this->Id()
                  << std::endl;
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::ComputePotentialJump(ProcessInfo& rCurrentProcessInfo)
{
    const array_1d<double, 3> vinfinity = rCurrentProcessInfo[VELOCITY_INFINITY];
    const double vinfinity_norm = sqrt(inner_prod(vinfinity, vinfinity));

    array_1d<double, NumNodes> distances;
    GetWakeDistances(distances);

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        double aux_potential = GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL);
        double potential = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
        double potential_jump = aux_potential - potential;

        if (distances[i] > 0)
        {
            GetGeometry()[i].SetValue(POTENTIAL_JUMP, -2.0 / vinfinity_norm * (potential_jump));
        }
        else
        {
            GetGeometry()[i].SetValue(POTENTIAL_JUMP, 2.0 / vinfinity_norm * (potential_jump));
        }
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::ComputeElementInternalEnergy()
{
    double internal_energy = 0.0;
    array_1d<double, Dim> velocity;

    ComputeUpperVelocity(velocity);

    internal_energy = 0.5 * inner_prod(velocity, velocity);
    this->SetValue(INTERNAL_ENERGY, std::abs(internal_energy));
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::GetPotential(
    Vector& split_element_values, const array_1d<double, NumNodes>& distances) const
{
    array_1d<double, NumNodes> upper_phis;
    GetPotentialOnUpperWakeElement(upper_phis, distances);

    array_1d<double, NumNodes> lower_phis;
    GetPotentialOnLowerWakeElement(lower_phis, distances);

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        split_element_values[i] = upper_phis[i];
        split_element_values[NumNodes + i] = lower_phis[i];
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::GetPotentialOnUpperWakeElement(
    array_1d<double, NumNodes>& upper_phis, const array_1d<double, NumNodes>& distances) const
{
    for (unsigned int i = 0; i < NumNodes; i++)
        if (distances[i] > 0)
            upper_phis[i] = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
        else
            upper_phis[i] = GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::GetPotentialOnLowerWakeElement(
    array_1d<double, NumNodes>& lower_phis, const array_1d<double, NumNodes>& distances) const
{
    for (unsigned int i = 0; i < NumNodes; i++)
        if (distances[i] < 0)
            lower_phis[i] = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
        else
            lower_phis[i] = GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::ComputeUpperVelocity(array_1d<double, Dim>& velocity) const
{
    velocity.clear();

    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    array_1d<double, NumNodes> distances;
    GetWakeDistances(distances);

    GetPotentialOnUpperWakeElement(data.phis, distances);

    noalias(velocity) = prod(trans(data.DN_DX), data.phis);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::ComputeLowerVelocity(array_1d<double, Dim>& velocity) const
{
    velocity.clear();

    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    array_1d<double, NumNodes> distances;
    GetWakeDistances(distances);

    GetPotentialOnLowerWakeElement(data.phis, distances);

    noalias(velocity) = prod(trans(data.DN_DX), data.phis);
}

template <int Dim, int NumNodes>
double IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::ComputeUpperPressure(const ProcessInfo& rCurrentProcessInfo) const
{
    const array_1d<double, 3> vinfinity = rCurrentProcessInfo[VELOCITY_INFINITY];
    const double vinfinity_norm2 = inner_prod(vinfinity, vinfinity);

    KRATOS_ERROR_IF(vinfinity_norm2 < std::numeric_limits<double>::epsilon())
        << "Error on element -> " << this->Id() << "\n"
        << "vinfinity_norm2 must be larger than zero." << std::endl;

    array_1d<double, Dim> v;
    ComputeUpperVelocity(v);

    double pressure_coefficient = (vinfinity_norm2 - inner_prod(v, v)) /
               vinfinity_norm2; // 0.5*(norm_2(vinfinity) - norm_2(v));

    return pressure_coefficient;
}

template <int Dim, int NumNodes>
double IncompressiblePotentialFlowTrailingEdgeElement<Dim, NumNodes>::ComputeLowerPressure(const ProcessInfo& rCurrentProcessInfo) const
{
    const array_1d<double, 3> vinfinity = rCurrentProcessInfo[VELOCITY_INFINITY];
    const double vinfinity_norm2 = inner_prod(vinfinity, vinfinity);

    KRATOS_ERROR_IF(vinfinity_norm2 < std::numeric_limits<double>::epsilon())
        << "Error on element -> " << this->Id() << "\n"
        << "vinfinity_norm2 must be larger than zero." << std::endl;

    array_1d<double, Dim> v;
    ComputeLowerVelocity(v);

    double pressure_coefficient = (vinfinity_norm2 - inner_prod(v, v)) /
               vinfinity_norm2; // 0.5*(norm_2(vinfinity) - norm_2(v));

    return pressure_coefficient;
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
