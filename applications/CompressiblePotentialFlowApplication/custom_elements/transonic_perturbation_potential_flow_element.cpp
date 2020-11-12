//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:     Eloisa Baez Jones, Inigo Lopez and Riccardo Rossi
//

#include "transonic_perturbation_potential_flow_element.h"
#include "compressible_potential_flow_application_variables.h"
#include "includes/cfd_variables.h"
#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/potential_flow_utilities.h"
#include "utilities/geometry_utilities.h"
#include "utilities/enrichment_utilities.h"

namespace Kratos
{
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <int TDim, int TNumNodes>
Element::Pointer TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<TransonicPerturbationPotentialFlowElement>(
        NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

template <int TDim, int TNumNodes>
Element::Pointer TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<TransonicPerturbationPotentialFlowElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

template <int TDim, int TNumNodes>
Element::Pointer TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::Clone(
    IndexType NewId,
    NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<TransonicPerturbationPotentialFlowElement>(
        NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    FindUpwindElement(rCurrentProcessInfo);
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);
    CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const TransonicPerturbationPotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) { // Normal element (non-wake) - eventually an embedded
        CalculateRightHandSideNormalElement(rRightHandSideVector, rCurrentProcessInfo);
    }
    else { // Wake element
        CalculateRightHandSideWakeElement(rRightHandSideVector, rCurrentProcessInfo);
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    const TransonicPerturbationPotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) { // Normal element (non-wake) - eventually an embedded
        if (r_this.IsNot(INLET)) {
            if (rLeftHandSideMatrix.size1() != TNumNodes + 1 ||
                rLeftHandSideMatrix.size2() != TNumNodes + 1) {
                    rLeftHandSideMatrix.resize(TNumNodes + 1, TNumNodes + 1, false);
                }

            rLeftHandSideMatrix.clear();
            CalculateLeftHandSideNormalElement(rLeftHandSideMatrix, rCurrentProcessInfo);
        } else {
            if (rLeftHandSideMatrix.size1() != TNumNodes ||
                rLeftHandSideMatrix.size2() != TNumNodes) {
                    rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
                }

            rLeftHandSideMatrix.clear();
            CalculateLeftHandSideSubsonicElement(rLeftHandSideMatrix, rCurrentProcessInfo);
        }

    }
    else { // Wake element
        CalculateLeftHandSideWakeElement(rLeftHandSideMatrix, rCurrentProcessInfo);
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& CurrentProcessInfo) const
{
    const TransonicPerturbationPotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) // Normal element
    {
        if (r_this.IsNot(INLET)) {
            if (rResult.size() != TNumNodes + 1) {
                rResult.resize(TNumNodes + 1, false);
            }
            GetEquationIdVectorExtendedElement(rResult);
        }
        else {
            if (rResult.size() != TNumNodes) {
                rResult.resize(TNumNodes, false);
            }
            GetEquationIdVectorNormalElement(rResult);
        }
    }
    else { // Wake element
        if (rResult.size() != 2 * TNumNodes)
        {
            rResult.resize(2 * TNumNodes, false);
        }

        GetEquationIdVectorWakeElement(rResult);
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                                                 const ProcessInfo& CurrentProcessInfo) const
{
    const TransonicPerturbationPotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) // Normal element
    {
        if (rElementalDofList.size() != TNumNodes)
        {
            rElementalDofList.resize(TNumNodes);
        }

        const int kutta = r_this.GetValue(KUTTA);

        if (kutta == 0)
        {
            GetDofListNormalElement(rElementalDofList);

        }
        else
        {
            GetDofListKuttaElement(rElementalDofList);
        }
    }
    else // wake element
    {
        if (rElementalDofList.size() != 2 * TNumNodes)
        {
            rElementalDofList.resize(2 * TNumNodes);
        }

        GetDofListWakeElement(rElementalDofList);
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    const bool active = Is(ACTIVE);

    const TransonicPerturbationPotentialFlowElement& r_const_this = *this;
    const int wake = r_const_this.GetValue(WAKE);

    if (wake != 0 && active == true)
    {
        ComputePotentialJump(rCurrentProcessInfo);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <int TDim, int TNumNodes>
int TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
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
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY_POTENTIAL, this->GetGeometry()[i]);
    }

    return out;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
    {
        rValues.resize(1);
    }
    if (rVariable == PRESSURE_COEFFICIENT)
    {
        rValues[0] = PotentialFlowUtilities::ComputePerturbationCompressiblePressureCoefficient<TDim, TNumNodes>(*this, rCurrentProcessInfo);
    }
    else if (rVariable == DENSITY)
    {
        const array_1d<double, TDim> velocity = PotentialFlowUtilities::ComputePerturbedVelocity<TDim,TNumNodes>(*this, rCurrentProcessInfo);
        const double local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<TDim, TNumNodes>(velocity, rCurrentProcessInfo);
        rValues[0] = PotentialFlowUtilities::ComputeDensity<TDim, TNumNodes>(local_mach_number_squared, rCurrentProcessInfo);
    }
    else if (rVariable == MACH)
    {
        const array_1d<double, TDim> velocity = PotentialFlowUtilities::ComputePerturbedVelocity<TDim,TNumNodes>(*this, rCurrentProcessInfo);
        rValues[0] = std::sqrt(PotentialFlowUtilities::ComputeLocalMachNumberSquared<TDim, TNumNodes>(velocity, rCurrentProcessInfo));
    }
    else if (rVariable == SOUND_VELOCITY)
    {
        const array_1d<double, TDim> velocity = PotentialFlowUtilities::ComputePerturbedVelocity<TDim,TNumNodes>(*this, rCurrentProcessInfo);
        rValues[0] = std::sqrt(PotentialFlowUtilities::ComputeLocalSpeedofSoundSquared<TDim, TNumNodes>(velocity, rCurrentProcessInfo));
    }
    else if (rVariable == WAKE)
    {
        const TransonicPerturbationPotentialFlowElement& r_this = *this;
        rValues[0] = r_this.GetValue(WAKE);
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<int>& rVariable,
    std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
    {
        rValues.resize(1);
    }
    if (rVariable == TRAILING_EDGE)
    {
        rValues[0] = this->GetValue(TRAILING_EDGE);
    }
    else if (rVariable == KUTTA)
    {
        rValues[0] = this->GetValue(KUTTA);
    }
    else if (rVariable == WAKE)
    {
        rValues[0] = this->GetValue(WAKE);
    }
    else if (rVariable == ZERO_VELOCITY_CONDITION)
    {
        rValues[0] = this->GetValue(ZERO_VELOCITY_CONDITION);
    }
    else if (rVariable == TRAILING_EDGE_ELEMENT)
    {
        rValues[0] = this->GetValue(TRAILING_EDGE_ELEMENT);
    }
    else if (rVariable == DECOUPLED_TRAILING_EDGE_ELEMENT)
    {
        rValues[0] = this->GetValue(DECOUPLED_TRAILING_EDGE_ELEMENT);
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
    {
        rValues.resize(1);
    }
    if (rVariable == VELOCITY)
    {
        array_1d<double, 3> v(3, 0.0);
        array_1d<double, TDim> velocity = PotentialFlowUtilities::ComputePerturbedVelocity<TDim,TNumNodes>(*this, rCurrentProcessInfo);
        for (unsigned int k = 0; k < TDim; k++)
        {
            v[k] = velocity[k];
        }
        rValues[0] = v;
    }
    else if (rVariable == PERTURBATION_VELOCITY)
    {
        array_1d<double, 3> v(3, 0.0);
        array_1d<double, TDim> vaux = PotentialFlowUtilities::ComputeVelocity<TDim,TNumNodes>(*this);
        for (unsigned int k = 0; k < TDim; k++)
        {
            v[k] = vaux[k];
        }
        rValues[0] = v;
    }
    else if (rVariable == VECTOR_TO_UPWIND_ELEMENT)
    {
        rValues[0] = pGetUpwindElement()->GetGeometry().Center() - this->GetGeometry().Center();
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <int TDim, int TNumNodes>
std::string TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "TransonicPerturbationPotentialFlowElement #" << Id();
    return buffer.str();
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "TransonicPerturbationPotentialFlowElement #" << Id();
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::PrintData(std::ostream& rOStream) const
{
    pGetGeometry()->PrintData(rOStream);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template <int TDim, int TNumNodes>
inline GlobalPointer<Element> TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::pGetUpwindElement() const
{
    KRATOS_ERROR_IF(mpUpwindElement.get() == nullptr)
        << "No upwind element found for element #" << this->Id() << std::endl;
    return mpUpwindElement;
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::GetWakeDistances(
    array_1d<double,
    TNumNodes>& distances) const
{
    noalias(distances) = GetValue(WAKE_ELEMENTAL_DISTANCES);
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::GetEquationIdVectorExtendedElement(
    EquationIdVectorType& rResult) const
{
    // Adding normal element contribution
    const TransonicPerturbationPotentialFlowElement& r_this = *this;
    const int kutta = r_this.GetValue(KUTTA);
    if(kutta == 0) {
        GetEquationIdVectorNormalElement(rResult);
    } else {
        GetEquationIdVectorKuttaElement(rResult);
    }

    // Adding the additional node equation Id
    AddUpwindEquationId(rResult);
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::AddUpwindEquationId(
    EquationIdVectorType& rResult) const
{
    const int additional_upwind_node_index = GetAdditionalUpwindNodeIndex();
    const auto& r_upstream_element = *pGetUpwindElement();
    const auto& r_upwind_geomtery = r_upstream_element.GetGeometry();
    const int upstream_kutta = r_upstream_element.GetValue(KUTTA);
    if (upstream_kutta == 0) { // upwind element is not kutta
        // TODO special treatment for upwind wake elements
        rResult[TNumNodes] = r_upwind_geomtery[additional_upwind_node_index].GetDof(VELOCITY_POTENTIAL).EquationId();
    } else { // upwind element is kutta
        if (!r_upwind_geomtery[additional_upwind_node_index].GetValue(TRAILING_EDGE)) {
            // upwind node is not trailing edge
            rResult[TNumNodes] = r_upwind_geomtery[additional_upwind_node_index].GetDof(VELOCITY_POTENTIAL).EquationId();
        } else {
            // upwind node is trailing edge
            rResult[TNumNodes] = r_upwind_geomtery[additional_upwind_node_index].GetDof(AUXILIARY_VELOCITY_POTENTIAL).EquationId();
        }
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::GetEquationIdVectorNormalElement(
    EquationIdVectorType& rResult) const
{
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        rResult[i] = GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::GetEquationIdVectorKuttaElement(
    EquationIdVectorType& rResult) const
{
    const auto& r_geometry = this->GetGeometry();
    // Kutta elements have only negative part
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        if (!r_geometry[i].GetValue(TRAILING_EDGE))
        {
            rResult[i] = r_geometry[i].GetDof(VELOCITY_POTENTIAL).EquationId();
        }
        else
        {
            rResult[i] = r_geometry[i].GetDof(AUXILIARY_VELOCITY_POTENTIAL).EquationId();
        }
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::GetEquationIdVectorWakeElement(
    EquationIdVectorType& rResult) const
{
    array_1d<double, TNumNodes> distances;
    GetWakeDistances(distances);

    // Positive part
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        if (distances[i] > 0.0)
        {
            rResult[i] = GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
        }
        else
        {
            rResult[i] =
                GetGeometry()[i].GetDof(AUXILIARY_VELOCITY_POTENTIAL, 0).EquationId();
        }
    }

    // Negative part - sign is opposite to the previous case
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        if (distances[i] < 0.0)
        {
            rResult[TNumNodes + i] =
                GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
        }
        else
        {
            rResult[TNumNodes + i] =
                GetGeometry()[i].GetDof(AUXILIARY_VELOCITY_POTENTIAL).EquationId();
        }
    }

}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::GetDofListNormalElement(DofsVectorType& rElementalDofList) const
{
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        rElementalDofList[i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::GetDofListKuttaElement(DofsVectorType& rElementalDofList) const
{
    const auto& r_geometry = this->GetGeometry();
    // Kutta elements have only negative part
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        if (!r_geometry[i].GetValue(TRAILING_EDGE))
        {
                rElementalDofList[i] = r_geometry[i].pGetDof(VELOCITY_POTENTIAL);
        }
        else
        {
                rElementalDofList[i] = r_geometry[i].pGetDof(AUXILIARY_VELOCITY_POTENTIAL);
        }
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::GetDofListWakeElement(DofsVectorType& rElementalDofList) const
{
    array_1d<double, TNumNodes> distances;
    GetWakeDistances(distances);

    // Positive part
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        if (distances[i] > 0)
        {
            rElementalDofList[i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
        }
        else
        {
            rElementalDofList[i] = GetGeometry()[i].pGetDof(AUXILIARY_VELOCITY_POTENTIAL);
        }
    }

    // Negative part - sign is opposite to the previous case
    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        if (distances[i] < 0)
        {
            rElementalDofList[TNumNodes + i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
        }
        else
        {
            rElementalDofList[TNumNodes + i] =
                GetGeometry()[i].pGetDof(AUXILIARY_VELOCITY_POTENTIAL);
        }
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateLeftHandSideSubsonicElement(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Calculate shape functions
    ElementalData<TNumNodes, TDim> data;
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    const array_1d<double, TDim> velocity = PotentialFlowUtilities::ComputePerturbedVelocity<TDim,TNumNodes>(*this, rCurrentProcessInfo);

    const double local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<TDim, TNumNodes>(velocity, rCurrentProcessInfo);

    const double density = PotentialFlowUtilities::ComputeDensity<TDim, TNumNodes>(local_mach_number_squared, rCurrentProcessInfo);

    const double DrhoDu2 = PotentialFlowUtilities::ComputeDensityDerivativeWRTVelocitySquared<TDim, TNumNodes>(local_mach_number_squared, rCurrentProcessInfo);

    const BoundedVector<double, TNumNodes> DNV = prod(data.DN_DX, velocity);

    BoundedMatrix<double, TNumNodes, TNumNodes> current_lhs = data.vol * density * prod(data.DN_DX, trans(data.DN_DX));
    current_lhs += data.vol * 2 * DrhoDu2 * outer_prod(DNV, trans(DNV));

    for (int i = 0; i < TNumNodes; i++)
    {
        for (int j = 0; j < TNumNodes; j++)
        {
            rLeftHandSideMatrix(i, j) = current_lhs(i, j);
        }
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateRightHandSideNormalElement(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const TransonicPerturbationPotentialFlowElement& r_this = *this;

    const array_1d<double, TDim> velocity = PotentialFlowUtilities::ComputePerturbedVelocity<TDim,TNumNodes>(r_this, rCurrentProcessInfo);

    double density = 0.0;

    if(r_this.IsNot(INLET)) {
        if (rRightHandSideVector.size() != TNumNodes + 1) {
            rRightHandSideVector.resize(TNumNodes + 1, false);
        }
        rRightHandSideVector.clear();

        const array_1d<double, TDim> upwind_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<TDim,TNumNodes>(*pGetUpwindElement(), rCurrentProcessInfo);
        density = PotentialFlowUtilities::ComputeUpwindedDensity<TDim, TNumNodes>(velocity, upwind_velocity, rCurrentProcessInfo);

    } else {
        if (rRightHandSideVector.size() != TNumNodes) {
            rRightHandSideVector.resize(TNumNodes, false);
        }
        rRightHandSideVector.clear();

        const double local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<TDim, TNumNodes>(velocity, rCurrentProcessInfo);
        density = PotentialFlowUtilities::ComputeDensity<TDim, TNumNodes>(local_mach_number_squared, rCurrentProcessInfo);
    }
    // Calculate shape functions
    ElementalData<TNumNodes, TDim> data;
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    const BoundedVector<double, TNumNodes> current_rhs = - data.vol * density * prod(data.DN_DX, velocity);

    for (int i = 0; i < TNumNodes; i++)
    {
        rRightHandSideVector[i] = current_rhs[i];
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateLeftHandSideNormalElement(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    const TransonicPerturbationPotentialFlowElement& r_this = *this;

    const array_1d<double, TDim> velocity = PotentialFlowUtilities::ComputePerturbedVelocity<TDim,TNumNodes>(r_this, rCurrentProcessInfo);
    const array_1d<double, TDim> upwind_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<TDim,TNumNodes>(*pGetUpwindElement(), rCurrentProcessInfo);

    const double local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<TDim, TNumNodes>(velocity, rCurrentProcessInfo);
    const double upwind_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<TDim, TNumNodes>(upwind_velocity, rCurrentProcessInfo);

    double DrhoDu2 = 0.0;
    double DrhoDu2_up = 0.0;

    const double critical_mach_sq = std::pow(rCurrentProcessInfo[CRITICAL_MACH], 2.0);

    if (local_mach_number_squared < critical_mach_sq) { // subsonic, not inlet
        // gets [TNumNodes + 1, TNumNodes + 1] size matrix
        CalculateLeftHandSideSubsonicElement(rLeftHandSideMatrix, rCurrentProcessInfo);
        return;
    }
    else if (local_mach_number_squared >= upwind_mach_number_squared) { // supersonic, accelerating
        // density derivatives
        DrhoDu2 = PotentialFlowUtilities::ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicAccelerating<TDim, TNumNodes>(
            velocity, local_mach_number_squared, upwind_mach_number_squared, rCurrentProcessInfo);
        DrhoDu2_up = PotentialFlowUtilities::ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicAccelerating<TDim, TNumNodes>(
            local_mach_number_squared, upwind_mach_number_squared, rCurrentProcessInfo);
    }
    else { // supersonic, deaccelerating
        // density derivatives
        DrhoDu2 = PotentialFlowUtilities::ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicDeaccelerating<TDim, TNumNodes>(
            local_mach_number_squared, upwind_mach_number_squared, rCurrentProcessInfo);
        DrhoDu2_up = PotentialFlowUtilities::ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicDeaccelerating<TDim, TNumNodes>(
            upwind_velocity, local_mach_number_squared, upwind_mach_number_squared, rCurrentProcessInfo);
    }

    AssembleSupersonicLeftHandSide(rLeftHandSideMatrix, DrhoDu2, DrhoDu2_up, velocity, upwind_velocity, rCurrentProcessInfo);
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateLeftHandSideWakeElement(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Note that the lhs and rhs have double the size
    if (rLeftHandSideMatrix.size1() != 2 * TNumNodes ||
        rLeftHandSideMatrix.size2() != 2 * TNumNodes)
    {
        rLeftHandSideMatrix.resize(2 * TNumNodes, 2 * TNumNodes, false);
    }
    rLeftHandSideMatrix.clear();

    ElementalData<TNumNodes, TDim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);
    GetWakeDistances(data.distances);

    // computing local velocity
    const array_1d<double, TDim> velocity = PotentialFlowUtilities::ComputePerturbedVelocity<TDim,TNumNodes>(*this, rCurrentProcessInfo);

    const double local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<TDim, TNumNodes>(velocity, rCurrentProcessInfo);

    const double density = PotentialFlowUtilities::ComputeDensity<TDim, TNumNodes>(local_mach_number_squared, rCurrentProcessInfo);

    const double DrhoDu2 = PotentialFlowUtilities::ComputeDensityDerivativeWRTVelocitySquared<TDim, TNumNodes>(local_mach_number_squared, rCurrentProcessInfo);

    const BoundedVector<double, TNumNodes> DNV = prod(data.DN_DX, velocity);

    const BoundedMatrix<double, TNumNodes, TNumNodes> lhs_total =
        data.vol * density * prod(data.DN_DX, trans(data.DN_DX)) +
        data.vol * 2 * DrhoDu2 * outer_prod(DNV, trans(DNV));

    if (this->Is(STRUCTURE))
    {
        Matrix lhs_positive = ZeroMatrix(TNumNodes, TNumNodes);
        Matrix lhs_negative = ZeroMatrix(TNumNodes, TNumNodes);

        CalculateLeftHandSideSubdividedElement(lhs_positive, lhs_negative, rCurrentProcessInfo);
        AssignLeftHandSideSubdividedElement(rLeftHandSideMatrix, lhs_positive,
                                            lhs_negative, lhs_total, data);
    }
    else
    {
        AssignLeftHandSideWakeElement(rLeftHandSideMatrix, lhs_total, data);
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateRightHandSideWakeElement(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Note that the rhs has double the size
    if (rRightHandSideVector.size() != 2 * TNumNodes)
    {
        rRightHandSideVector.resize(2 * TNumNodes, false);
    }
    rRightHandSideVector.clear();

    ElementalData<TNumNodes, TDim> data;

    // Calculate shape functions
    const auto& r_geometry = this->GetGeometry();
    GeometryUtils::CalculateGeometryData(r_geometry, data.DN_DX, data.N, data.vol);
    GetWakeDistances(data.distances);

    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    array_1d<double, TDim> upper_velocity = PotentialFlowUtilities::ComputeVelocityUpperWakeElement<TDim,TNumNodes>(*this);
    array_1d<double, TDim> lower_velocity = PotentialFlowUtilities::ComputeVelocityLowerWakeElement<TDim,TNumNodes>(*this);

    for (unsigned int i = 0; i < TDim; i++)
    {
        upper_velocity[i] += free_stream_velocity[i];
        lower_velocity[i] += free_stream_velocity[i];
    }

    const double local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<TDim, TNumNodes>(upper_velocity, rCurrentProcessInfo);
    const double density = PotentialFlowUtilities::ComputeDensity<TDim, TNumNodes>(local_mach_number_squared, rCurrentProcessInfo);

    const array_1d<double, TDim> diff_velocity = upper_velocity - lower_velocity;

    const BoundedVector<double, TNumNodes> upper_rhs = - data.vol * density * prod(data.DN_DX, upper_velocity);
    const BoundedVector<double, TNumNodes> lower_rhs = - data.vol * density * prod(data.DN_DX, lower_velocity);
    const BoundedVector<double, TNumNodes> wake_rhs = - data.vol * density * prod(data.DN_DX, diff_velocity);

    if (this->Is(STRUCTURE))
    {
        double upper_vol = 0.0;
        double lower_vol = 0.0;

        CalculateVolumesSubdividedElement(upper_vol, lower_vol, rCurrentProcessInfo);
        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            if (r_geometry[i].GetValue(TRAILING_EDGE))
            {
                rRightHandSideVector[i] = upper_rhs(i)*upper_vol/data.vol;
                rRightHandSideVector[i + TNumNodes] = lower_rhs(i)*lower_vol/data.vol;
            }
            else
            {
                AssignRightHandSideWakeNode(rRightHandSideVector, upper_rhs, lower_rhs, wake_rhs, data, i);
            }
        }
    }
    else
    {
        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            AssignRightHandSideWakeNode(rRightHandSideVector, upper_rhs, lower_rhs, wake_rhs, data, i);
        }
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateLeftHandSideSubdividedElement(
    Matrix& lhs_positive,
    Matrix& lhs_negative,
    const ProcessInfo& rCurrentProcessInfo)
{
    ElementalData<TNumNodes, TDim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    GetWakeDistances(data.distances);

    // Subdivide the element
    constexpr unsigned int nvolumes = 3 * (TDim - 1);
    BoundedMatrix<double, TNumNodes, TDim> Points;
    array_1d<double, nvolumes> partitions_sign;
    BoundedMatrix<double, nvolumes, TNumNodes> gp_shape_function_values;
    array_1d<double, nvolumes> Volumes;
    std::vector<Matrix> GradientsValue(nvolumes);
    BoundedMatrix<double, nvolumes, 2> N_enriched;
    for (unsigned int i = 0; i < GradientsValue.size(); ++i)
    {
        GradientsValue[i].resize(2, TDim, false);
    }
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        const auto& r_coords = GetGeometry()[i].Coordinates();
        for (unsigned int k = 0; k < TDim; ++k)
        {
            Points(i, k) = r_coords[k];
        }
    }

    const unsigned int nsubdivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(
        Points, data.DN_DX, data.distances, Volumes, gp_shape_function_values,
        partitions_sign, GradientsValue, N_enriched);

    // Computing local velocity
    const array_1d<double, TDim> velocity = PotentialFlowUtilities::ComputePerturbedVelocity<TDim,TNumNodes>(*this, rCurrentProcessInfo);

    const double local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<TDim, TNumNodes>(velocity, rCurrentProcessInfo);

    const double density = PotentialFlowUtilities::ComputeDensity<TDim, TNumNodes>(local_mach_number_squared, rCurrentProcessInfo);

    const double DrhoDu2 = PotentialFlowUtilities::ComputeDensityDerivativeWRTVelocitySquared<TDim, TNumNodes>(local_mach_number_squared, rCurrentProcessInfo);

    const BoundedVector<double, TNumNodes> DNV = prod(data.DN_DX, velocity);

    // Compute the lhs and rhs that would correspond to it being divided
    for (unsigned int i = 0; i < nsubdivisions; ++i)
    {
        if (partitions_sign[i] > 0)
        {
            noalias(lhs_positive) +=
                Volumes[i] * density * prod(data.DN_DX, trans(data.DN_DX));
            noalias(lhs_positive) +=
                Volumes[i] * 2 * DrhoDu2 * outer_prod(DNV, trans(DNV));
        }
        else
        {
            noalias(lhs_negative) +=
                Volumes[i] * density * prod(data.DN_DX, trans(data.DN_DX));
            noalias(lhs_negative) +=
                Volumes[i] * 2 * DrhoDu2 * outer_prod(DNV, trans(DNV));
        }
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateVolumesSubdividedElement(
    double& rUpper_vol,
    double& rLower_vol,
    const ProcessInfo& rCurrentProcessInfo)
{
    ElementalData<TNumNodes, TDim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    GetWakeDistances(data.distances);

    // Subdivide the element
    constexpr unsigned int nvolumes = 3 * (TDim - 1);
    BoundedMatrix<double, TNumNodes, TDim> Points;
    array_1d<double, nvolumes> PartitionsSign;
    BoundedMatrix<double, nvolumes, TNumNodes> GPShapeFunctionValues;
    array_1d<double, nvolumes> Volumes;
    std::vector<Matrix> GradientsValue(nvolumes);
    BoundedMatrix<double, nvolumes, 2> NEnriched;
    for (unsigned int i = 0; i < GradientsValue.size(); ++i)
    {
        GradientsValue[i].resize(2, TDim, false);
    }
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        const array_1d<double, 3>& coords = GetGeometry()[i].Coordinates();
        for (unsigned int k = 0; k < TDim; ++k)
        {
            Points(i, k) = coords[k];
        }
    }

    const unsigned int nsubdivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(
        Points, data.DN_DX, data.distances, Volumes, GPShapeFunctionValues,
        PartitionsSign, GradientsValue, NEnriched);

    // Compute the volumes that would correspond to it being divided
    for (unsigned int i = 0; i < nsubdivisions; ++i)
    {
        if (PartitionsSign[i] > 0)
        {
            rUpper_vol += Volumes[i];
        }
        else
        {
            rLower_vol += Volumes[i];
        }
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::ComputeLHSGaussPointContribution(
    const double weight,
    Matrix& lhs,
    const ElementalData<TNumNodes, TDim>& data) const
{
    noalias(lhs) += weight * prod(data.DN_DX, trans(data.DN_DX));
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::AssignLeftHandSideSubdividedElement(
    Matrix& rLeftHandSideMatrix,
    Matrix& lhs_positive,
    Matrix& lhs_negative,
    const BoundedMatrix<double, TNumNodes, TNumNodes>& lhs_total,
    const ElementalData<TNumNodes, TDim>& data) const
{
    const auto& r_geometry = this->GetGeometry();
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        // The TE node takes the contribution of the subdivided element and
        // we do not apply the wake condition on the TE node
        if (r_geometry[i].GetValue(TRAILING_EDGE))
        {
            for (unsigned int j = 0; j < TNumNodes; ++j)
            {
                rLeftHandSideMatrix(i, j) = lhs_positive(i, j);
                rLeftHandSideMatrix(i + TNumNodes, j + TNumNodes) = lhs_negative(i, j);
            }
        }
        else
        {
            AssignLeftHandSideWakeNode(rLeftHandSideMatrix, lhs_total, data, i);
        }
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::AssignLeftHandSideWakeElement(
    MatrixType& rLeftHandSideMatrix,
    const BoundedMatrix<double, TNumNodes, TNumNodes>& lhs_total,
    const ElementalData<TNumNodes, TDim>& data) const
{
    for (unsigned int row = 0; row < TNumNodes; ++row)
    {
        AssignLeftHandSideWakeNode(rLeftHandSideMatrix, lhs_total, data, row);
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::AssignLeftHandSideWakeNode(
    MatrixType& rLeftHandSideMatrix,
    const BoundedMatrix<double, TNumNodes, TNumNodes>& lhs_total,
    const ElementalData<TNumNodes, TDim>& data,
    unsigned int& row) const
{
    // Filling the diagonal blocks (i.e. decoupling upper and lower dofs)
    for (unsigned int column = 0; column < TNumNodes; ++column)
    {
        rLeftHandSideMatrix(row, column) = lhs_total(row, column);
        rLeftHandSideMatrix(row + TNumNodes, column + TNumNodes) = lhs_total(row, column);
    }

    // Applying wake condition on the AUXILIARY_VELOCITY_POTENTIAL dofs
    if (data.distances[row] < 0.0)
    {
        for (unsigned int column = 0; column < TNumNodes; ++column)
        {
            rLeftHandSideMatrix(row, column + TNumNodes) = -lhs_total(row, column); // Side 1
        }
    }
    else if (data.distances[row] > 0.0)
    {
        for (unsigned int column = 0; column < TNumNodes; ++column)
        {
            rLeftHandSideMatrix(row + TNumNodes, column) = -lhs_total(row, column); // Side 2
        }
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::AssignRightHandSideWakeNode(
    VectorType& rRightHandSideVector,
    const BoundedVector<double, TNumNodes>& rUpper_rhs,
    const BoundedVector<double, TNumNodes>& rLower_rhs,
    const BoundedVector<double, TNumNodes>& rWake_rhs,
    const ElementalData<TNumNodes, TDim>& rData,
    unsigned int& rRow) const
{
    if (rData.distances[rRow] > 0.0)
    {
        rRightHandSideVector[rRow] = rUpper_rhs(rRow);
        // TODO: Check theory behind setting rhs to 0
        rRightHandSideVector[rRow + TNumNodes] = 0.0; //rWake_rhs(rRow);
    }
    else
    {
        // TODO: Check theory behind setting rhs to 0
        rRightHandSideVector[rRow] = 0.0; //rWake_rhs(rRow);
        rRightHandSideVector[rRow + TNumNodes] = rLower_rhs(rRow);
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::AssembleSupersonicLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const double densityDerivativeWRTVelocity,
    const double densityDerivativeWRTUpwindVelocity,
    const array_1d<double, TDim> velocity,
    const array_1d<double, TDim> upwindVelocity,
    const ProcessInfo& rCurrentProcessInfo)
{
    BoundedVector<double, TNumNodes + 1> DNV_assembly = AssembleDensityDerivativeAndShapeFunctions(
        densityDerivativeWRTVelocity, densityDerivativeWRTUpwindVelocity, velocity, upwindVelocity, rCurrentProcessInfo);

    // Calculate shape functions
    ElementalData<TNumNodes, TDim> data;
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    const double density = PotentialFlowUtilities::ComputeUpwindedDensity<TDim, TNumNodes>(velocity, upwindVelocity, rCurrentProcessInfo);

    const BoundedVector<double, TNumNodes> DNV = prod(data.DN_DX, velocity);
    BoundedVector<double, TNumNodes + 1> DNV_extended;
    DNV_extended.clear();
    for (int i = 0; i < TNumNodes; i++) {
        DNV_extended[i] = DNV[i];
    }

    BoundedMatrix<double, TNumNodes, TNumNodes> linear_term = data.vol * density * prod(data.DN_DX, trans(data.DN_DX));

    rLeftHandSideMatrix = data.vol * 2.0 * outer_prod(DNV_extended, trans(DNV_assembly));

    for (int i = 0; i < TNumNodes; i++)
    {
        for (int j = 0; j < TNumNodes; j++)
        {
            rLeftHandSideMatrix(i, j) += linear_term(i, j);
        }
    }
}

template <int TDim, int TNumNodes>
BoundedVector<double, TNumNodes + 1> TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::AssembleDensityDerivativeAndShapeFunctions(
    const double densityDerivativeWRTVelocitySquared,
    const double densityDerivativeWRTUpwindVelocitySquared,
    const array_1d<double, TDim>& velocity,
    const array_1d<double, TDim>& upwindVelocity,
    const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& r_geom = this->GetGeometry();
    const GeometryType& r_upwind_geom = pGetUpwindElement()->GetGeometry();

    const array_1d<size_t, TNumNodes> upwind_node_key = GetAssemblyKey(r_geom, r_upwind_geom, rCurrentProcessInfo);

    ElementalData<TNumNodes, TDim> currentElementdata;
    ElementalData<TNumNodes, TDim> upwindElementdata;

    GeometryUtils::CalculateGeometryData(r_geom, currentElementdata.DN_DX, currentElementdata.N, currentElementdata.vol);
    GeometryUtils::CalculateGeometryData(r_upwind_geom, upwindElementdata.DN_DX, upwindElementdata.N, upwindElementdata.vol);

    const BoundedVector<double, TNumNodes> current_DNV = densityDerivativeWRTVelocitySquared * prod(currentElementdata.DN_DX, velocity);
    const BoundedVector<double, TNumNodes> upwind_DNV = densityDerivativeWRTUpwindVelocitySquared * prod(upwindElementdata.DN_DX, upwindVelocity);

    BoundedVector<double, TNumNodes + 1> assembly_DNV;
    assembly_DNV.clear();

    for (int i = 0; i < TNumNodes; i++)
    {
        assembly_DNV[i] += current_DNV[i];
        assembly_DNV[upwind_node_key[i]] += upwind_DNV[i];
    }

    return assembly_DNV;
}

template <int TDim, int TNumNodes>
array_1d<size_t, TNumNodes> TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::GetAssemblyKey(
    const GeometryType& rGeom,
    const GeometryType& rUpwindGeom,
    const ProcessInfo& rCurrentProcessInfo)
{
    array_1d<size_t, TNumNodes> upwind_node_key;
    EquationIdVectorType upwind_element_ids, current_element_ids;

    upwind_node_key.clear();

    pGetUpwindElement()->EquationIdVector(upwind_element_ids, rCurrentProcessInfo);
    this->EquationIdVector(current_element_ids, rCurrentProcessInfo);

    for (int i = 0; i < TNumNodes; i++) {
        auto current_id = std::find(current_element_ids.begin(), current_element_ids.end(),
                upwind_element_ids[i]);

        upwind_node_key[i] = std::distance(current_element_ids.begin(), current_id);
    }

    return upwind_node_key;
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::ComputePotentialJump(const ProcessInfo& rCurrentProcessInfo)
{
    const array_1d<double, 3>& vinfinity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    const double v_infinity_norm = sqrt(inner_prod(vinfinity, vinfinity));

    array_1d<double, TNumNodes> distances;
    GetWakeDistances(distances);

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        double aux_potential =
            GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL);
        double potential = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
        double potential_jump = aux_potential - potential;
        if (distances[i] > 0)
        {
            GetGeometry()[i].SetValue(POTENTIAL_JUMP,
                                      -2.0 / v_infinity_norm * (potential_jump));
        }
        else
        {
            GetGeometry()[i].SetValue(POTENTIAL_JUMP, 2.0 / v_infinity_norm * (potential_jump));
        }
    }
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::FindUpwindElement(const ProcessInfo& rCurrentProcessInfo)
{
    GeometryType upwind_element_boundary;
    FindUpwindEdge(upwind_element_boundary, rCurrentProcessInfo);
    std::vector<size_t> upwind_element_nodes;
    PotentialFlowUtilities::GetSortedIds<TDim, TNumNodes>(upwind_element_nodes, upwind_element_boundary);

    GlobalPointersVector<Element> upwind_element_candidates;
    PotentialFlowUtilities::GetNodeNeighborElementCandidates<TDim, TNumNodes>(upwind_element_candidates, upwind_element_boundary);
    SelectUpwindElement(upwind_element_nodes, upwind_element_candidates);
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::FindUpwindEdge(GeometryType& rUpwindEdge,
    const ProcessInfo& rCurrentProcessInfo)
{

    GeometriesArrayType element_boundary_geometry;
    GetElementGeometryBoundary(element_boundary_geometry);

    // free stream values
    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];

    double minimum_edge_flux = 0.0;
    for (SizeType i = 0; i < element_boundary_geometry.size(); i++)
    {
        const auto edge_normal = GetEdgeNormal(element_boundary_geometry[i]);

        const double edge_flux = inner_prod(edge_normal, free_stream_velocity);

        if(edge_flux < minimum_edge_flux)
        {
            minimum_edge_flux = edge_flux;
            rUpwindEdge = element_boundary_geometry[i];
        }
    }
}

template<int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::GetElementGeometryBoundary(GeometriesArrayType& rElementGeometryBoundary)
{
    const TransonicPerturbationPotentialFlowElement& r_this = *this;

    // current element geometry
    const GeometryType& r_geom = r_this.GetGeometry();

    // get element edges or faces depending on dimension of the problem
    if(TDim == 2)
    {
        // current element edges
        rElementGeometryBoundary = r_geom.GenerateEdges();
    }
    else if(TDim == 3)
    {
        // current element faces
        rElementGeometryBoundary = r_geom.GenerateFaces();
    }
}

template <int TDim, int TNumNodes>
array_1d<double, 3> TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::GetEdgeNormal(const GeometryType& rEdge)
{
    // get local coordinates of edge center
    array_1d<double, 3> edge_center_coordinates;
    rEdge.PointLocalCoordinates(edge_center_coordinates, rEdge.Center());

    // outward pointing normals of each edge
    return rEdge.Normal(edge_center_coordinates);
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::SelectUpwindElement(
    std::vector<IndexType>& rUpwindElementNodeIds,
    GlobalPointersVector<Element>& rUpwindElementCandidates)
{
    for (SizeType i = 0; i < rUpwindElementCandidates.size(); i++)
    {
        // get sorted node ids of neighbording elements
        std::vector<size_t> neighbor_element_ids;
        PotentialFlowUtilities::GetSortedIds<TDim, TNumNodes>(neighbor_element_ids, rUpwindElementCandidates[i].GetGeometry());

        // find element which shares the upwind element nodes with current element
        // but is not the current element
        if(std::includes(neighbor_element_ids.begin(), neighbor_element_ids.end(),
            rUpwindElementNodeIds.begin(), rUpwindElementNodeIds.end())
            && rUpwindElementCandidates[i].Id() != this->Id())
        {
            mpUpwindElement = rUpwindElementCandidates(i);
            break;
        }
    }

    // If no upwind element is found, the element is an INLET element and the
    // upwind element pointer points to itself
    if (mpUpwindElement.get() == nullptr)
    {
        mpUpwindElement = this;
        this->SetFlags(INLET);
    }
}

template<int TDim, int TNumNodes>
int TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::GetAdditionalUpwindNodeIndex() const
{
    // current and upwind element geometry
    const GeometryType& r_geom = this->GetGeometry();
    const GeometryType& r_upwind_geom = pGetUpwindElement()->GetGeometry();
    std::vector<size_t> element_nodes_ids;
    PotentialFlowUtilities::GetSortedIds<TDim, TNumNodes>(element_nodes_ids, r_geom);

    // Search for the Id of the upwind element node that
    // is not contained in the current element
    bool upstream_element_id_found = false;
    // loop over upwind element nodes
    for (unsigned int i = 0; i < TNumNodes; i++) {
        if( std::find(element_nodes_ids.begin(), element_nodes_ids.end(),
            r_upwind_geom[i].Id()) == element_nodes_ids.end() )  {
                upstream_element_id_found = true;
                return i;
            }
    }

    KRATOS_ERROR_IF(!upstream_element_id_found) << "No upstream element id found for element #"
            << this->Id() << std::endl;

    return -1;
}

// serializer

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}

template <int TDim, int TNumNodes>
void TransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class TransonicPerturbationPotentialFlowElement<2, 3>;
template class TransonicPerturbationPotentialFlowElement<3, 4>;

} // namespace Kratos
