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

#include "custom_elements/incompressible_potential_flow_functions.h"

namespace Kratos
{
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <int Dim, int NumNodes>
Element::Pointer IncompressiblePotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_shared<IncompressiblePotentialFlowElement>(
        NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer IncompressiblePotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_shared<IncompressiblePotentialFlowElement>(
        NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer IncompressiblePotentialFlowElement<Dim, NumNodes>::Clone(
    IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_shared<IncompressiblePotentialFlowElement>(
        NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != NumNodes || rLeftHandSideMatrix.size2() != NumNodes)
        rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);
    if (rRightHandSideVector.size() != NumNodes)
        rRightHandSideVector.resize(NumNodes, false);
    rLeftHandSideMatrix.clear();

    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    ComputeLHSGaussPointContribution<Dim,NumNodes>(data.vol, rLeftHandSideMatrix, data);

    GetPotential<2,3>(this,data.phis);
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, data.phis);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    // TODO: improve speed
    Matrix tmp;
    CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::EquationIdVector(
    EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    if (rResult.size() != NumNodes)
        rResult.resize(NumNodes, false);

    for (unsigned int i = 0; i < NumNodes; i++)
        rResult[i] = GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                                                   ProcessInfo& CurrentProcessInfo)
{
    if (rElementalDofList.size() != NumNodes)
        rElementalDofList.resize(NumNodes);

    for (unsigned int i = 0; i < NumNodes; i++)
        rElementalDofList[i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    ComputeElementInternalEnergy();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <int Dim, int NumNodes>
int IncompressiblePotentialFlowElement<Dim, NumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
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
    }

    return out;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::GetValueOnIntegrationPoints(
    const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == PRESSURE)
    {
        double p = ComputePressure(rCurrentProcessInfo);
        rValues[0] = p;
    }
    else if (rVariable == WAKE)
    {
        const IncompressiblePotentialFlowElement& r_this = *this;
        rValues[0] = r_this.GetValue(WAKE);
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::GetValueOnIntegrationPoints(
    const Variable<int>& rVariable, std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);
    if (rVariable == TRAILING_EDGE)
        rValues[0] = this->GetValue(TRAILING_EDGE);
    else if (rVariable == KUTTA)
        rValues[0] = this->GetValue(KUTTA);
    else if (rVariable == WAKE)
        rValues[0] = this->GetValue(WAKE);
    else if (rVariable == ZERO_VELOCITY_CONDITION)
        rValues[0] = this->GetValue(ZERO_VELOCITY_CONDITION);
    else if (rVariable == DECOUPLED_TRAILING_EDGE_ELEMENT)
        rValues[0] = this->GetValue(DECOUPLED_TRAILING_EDGE_ELEMENT);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::GetValueOnIntegrationPoints(
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
        ComputeVelocity<Dim,NumNodes>(this,vaux);
        for (unsigned int k = 0; k < Dim; k++)
            v[k] = vaux[k];
        rValues[0] = v;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <int Dim, int NumNodes>
std::string IncompressiblePotentialFlowElement<Dim, NumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "IncompressiblePotentialFlowElement #" << Id();
    return buffer.str();
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "IncompressiblePotentialFlowElement #" << Id();
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::PrintData(std::ostream& rOStream) const
{
    pGetGeometry()->PrintData(rOStream);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// template <int Dim, int NumNodes>
// void IncompressiblePotentialFlowElement<Dim, NumNodes>::ComputeLHSGaussPointContribution(
//     const double weight, Matrix& lhs, const ElementalData<NumNodes, Dim>& data) const
// {
//     noalias(lhs) += weight * prod(data.DN_DX, trans(data.DN_DX));
// }

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::ComputeElementInternalEnergy()
{
    double internal_energy = 0.0;
    array_1d<double, Dim> velocity;

    ComputeVelocity<Dim,NumNodes>(this,velocity);

    internal_energy = 0.5 * inner_prod(velocity, velocity);
    this->SetValue(INTERNAL_ENERGY, std::abs(internal_energy));
}

// template <int Dim, int NumNodes>
// void IncompressiblePotentialFlowElement<Dim, NumNodes>::GetPotential(
//     array_1d<double, NumNodes>& phis) const
// {
//     for (unsigned int i = 0; i < NumNodes; i++)
//         phis[i] = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
// }

// template <int Dim, int NumNodes>
// void IncompressiblePotentialFlowElement<Dim, NumNodes>::ComputeVelocity(
//     array_1d<double, Dim>& velocity) const
// {
//     ElementalData<NumNodes, Dim> data;

//     // Calculate shape functions
//     GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

//     GetPotential<2,3>(this,data.phis);

//     noalias(velocity) = prod(trans(data.DN_DX), data.phis);
// }

template <int Dim, int NumNodes>
double IncompressiblePotentialFlowElement<Dim, NumNodes>::ComputePressure(const ProcessInfo& rCurrentProcessInfo) const
{
    const array_1d<double, 3> vinfinity = rCurrentProcessInfo[VELOCITY_INFINITY];
    const double vinfinity_norm2 = inner_prod(vinfinity, vinfinity);

    KRATOS_ERROR_IF(vinfinity_norm2 < std::numeric_limits<double>::epsilon())
        << "Error on element -> " << this->Id() << "\n"
        << "vinfinity_norm2 must be larger than zero." << std::endl;

    array_1d<double, Dim> v;
    ComputeVelocity<Dim,NumNodes>(this,v);

    double pressure_coefficient = (vinfinity_norm2 - inner_prod(v, v)) /
               vinfinity_norm2; // 0.5*(norm_2(vinfinity) - norm_2(v));
    return pressure_coefficient;
}

// serializer

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class IncompressiblePotentialFlowElement<2, 3>;
template struct ElementalData<3, 2>;

} // namespace Kratos
