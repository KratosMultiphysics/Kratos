//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:     Marc Nuñez, Eloisa Baez Jones, Inigo Lopez and Riccardo Rossi
//

#include "embedded_transonic_perturbation_potential_flow_element.h"
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
Element::Pointer EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    typename PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EmbeddedTransonicPerturbationPotentialFlowElement>(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

template <int TDim, int TNumNodes>
Element::Pointer EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::Create(
    IndexType NewId,
    typename GeometryType::Pointer pGeom,
    typename PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EmbeddedTransonicPerturbationPotentialFlowElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

template <int TDim, int TNumNodes>
Element::Pointer EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::Clone(
    IndexType NewId,
    NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EmbeddedTransonicPerturbationPotentialFlowElement>(
        NewId, this->GetGeometry().Create(ThisNodes), this->pGetProperties());
    KRATOS_CATCH("");
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{

    BaseType::Initialize(rCurrentProcessInfo);

}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const EmbeddedTransonicPerturbationPotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    BoundedVector<double,TNumNodes> distances;
    for(unsigned int i_node = 0; i_node<TNumNodes; i_node++){
        distances[i_node] = this->GetGeometry()[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
    }
    const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<TDim,TNumNodes>(distances);

    if (is_embedded && wake == 0) {
        CalculateRightHandSideNormalElement(rRightHandSideVector,rCurrentProcessInfo);
        // if (std::abs(rCurrentProcessInfo[STABILIZATION_FACTOR]) > std::numeric_limits<double>::epsilon()) {
        //     PotentialFlowUtilities::AddPotentialGradientStabilizationTerm<Dim, NumNodes>(*this,rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
        // }
    }
    else {
        if (wake==0) {
            BaseType::CalculateRightHandSideNormalElement(rRightHandSideVector, rCurrentProcessInfo);
            // CalculateKuttaWakeLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
        } else {

            if (this->Is(STRUCTURE)) {
                CalculateRightHandSideKuttaWakeElement(rRightHandSideVector, rCurrentProcessInfo);
            }
            else {
                BaseType::CalculateRightHandSideWakeElement(rRightHandSideVector, rCurrentProcessInfo);
            }

        }
    }

    if (std::abs(rCurrentProcessInfo[PENALTY_COEFFICIENT]) > std::numeric_limits<double>::epsilon()) {
        PotentialFlowUtilities::AddKuttaConditionPenaltyRightHandSideTerm<TDim,TNumNodes>(r_this,rRightHandSideVector,rCurrentProcessInfo);
    }
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    const EmbeddedTransonicPerturbationPotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    BoundedVector<double,TNumNodes> distances;
    for(unsigned int i_node = 0; i_node<TNumNodes; i_node++){
        distances[i_node] = this->GetGeometry()[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
    }
    const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<TDim,TNumNodes>(distances);

    if (wake == 0) { // Normal element (non-wake) - eventually an embedded
        if (r_this.IsNot(INLET)) {
            if (rLeftHandSideMatrix.size1() != TNumNodes + 1 ||
                rLeftHandSideMatrix.size2() != TNumNodes + 1) {
                    rLeftHandSideMatrix.resize(TNumNodes + 1, TNumNodes + 1, false);
                }

            rLeftHandSideMatrix.clear();

            if (is_embedded)
                CalculateLeftHandSideNormalElement(rLeftHandSideMatrix, rCurrentProcessInfo);
            else
                BaseType::CalculateLeftHandSideNormalElement(rLeftHandSideMatrix, rCurrentProcessInfo);
        } else {
            if (rLeftHandSideMatrix.size1() != TNumNodes ||
                rLeftHandSideMatrix.size2() != TNumNodes) {
                    rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
                }

            rLeftHandSideMatrix.clear();
            if (is_embedded)
                CalculateLeftHandSideSubsonicElement(rLeftHandSideMatrix, rCurrentProcessInfo);
            else
                BaseType::CalculateLeftHandSideNormalElement(rLeftHandSideMatrix, rCurrentProcessInfo);

        }

    }
    else { // Wake element
        if (this->Is(STRUCTURE)) {
            CalculateLeftHandSideKuttaWakeElement(rLeftHandSideMatrix, rCurrentProcessInfo);
        }
        else {
            BaseType::CalculateLeftHandSideWakeElement(rLeftHandSideMatrix, rCurrentProcessInfo);
        }
    }

    if (std::abs(rCurrentProcessInfo[PENALTY_COEFFICIENT]) > std::numeric_limits<double>::epsilon()) {
        PotentialFlowUtilities::AddKuttaConditionPenaltyLeftHandSideTerm<TDim,TNumNodes>(r_this,rLeftHandSideMatrix,rCurrentProcessInfo);
    }
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateRightHandSideKuttaWakeElement(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Note that the rhs has double the size
    if (rRightHandSideVector.size() != 2 * TNumNodes)
    {
        rRightHandSideVector.resize(2 * TNumNodes, false);
    }
    rRightHandSideVector.clear();

    ElementalData data;

    // Calculate shape functions
    const auto& r_geometry = this->GetGeometry();
    GeometryUtils::CalculateGeometryData(r_geometry, data.DN_DX, data.N, data.vol);

    const array_1d<double, 3>& free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    array_1d<double, TDim> upper_velocity = PotentialFlowUtilities::ComputeVelocityUpperWakeElement<TDim,TNumNodes>(*this);
    array_1d<double, TDim> lower_velocity = PotentialFlowUtilities::ComputeVelocityLowerWakeElement<TDim,TNumNodes>(*this);

    for (unsigned int i = 0; i < TDim; i++)
    {
        upper_velocity[i] += free_stream_velocity[i];
        lower_velocity[i] += free_stream_velocity[i];
    }

    const double local_mach_number_squared_upper = PotentialFlowUtilities::ComputeLocalMachNumberSquared<TDim, TNumNodes>(upper_velocity, rCurrentProcessInfo);
    const double density_upper = PotentialFlowUtilities::ComputeDensity<TDim, TNumNodes>(local_mach_number_squared_upper, rCurrentProcessInfo);

    const double local_mach_number_squared_lower = PotentialFlowUtilities::ComputeLocalMachNumberSquared<TDim, TNumNodes>(lower_velocity, rCurrentProcessInfo);
    const double density_lower = PotentialFlowUtilities::ComputeDensity<TDim, TNumNodes>(local_mach_number_squared_lower, rCurrentProcessInfo);

    const BoundedVector<double, TNumNodes> upper_rhs = - data.vol * density_upper * prod(data.DN_DX, upper_velocity);
    const BoundedVector<double, TNumNodes> lower_rhs = - data.vol * density_lower * prod(data.DN_DX, lower_velocity);

    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        rRightHandSideVector[i] = upper_rhs(i);
        rRightHandSideVector[i + TNumNodes] = lower_rhs(i);
    }
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateLeftHandSideKuttaWakeElement(
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

    ElementalData data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), data.DN_DX, data.N, data.vol);

    // Compute upper and lower velocities
    const array_1d<double, TDim> upper_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<TDim,TNumNodes>(*this, rCurrentProcessInfo);
    const array_1d<double, TDim> lower_velocity = PotentialFlowUtilities::ComputePerturbedVelocityLowerElement<TDim,TNumNodes>(*this, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TNumNodes> upper_lhs_total = ZeroMatrix(TNumNodes,TNumNodes);
    BoundedMatrix<double, TNumNodes, TNumNodes> lower_lhs_total = ZeroMatrix(TNumNodes,TNumNodes);

    CalculateLeftHandSideContribution(upper_lhs_total, rCurrentProcessInfo, upper_velocity, data);
    CalculateLeftHandSideContribution(lower_lhs_total, rCurrentProcessInfo, lower_velocity, data);

    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        for (unsigned int j = 0; j < TNumNodes; ++j)
        {
            rLeftHandSideMatrix(i, j) = upper_lhs_total(i, j);
            rLeftHandSideMatrix(i + TNumNodes, j + TNumNodes) = lower_lhs_total(i, j);

        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <int TDim, int TNumNodes>
int EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Generic geometry check
    int out = Element::Check(rCurrentProcessInfo);
    if (out != 0)
    {
        return out;
    }

    KRATOS_ERROR_IF(this->GetGeometry().Area() <= 0.0)
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
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<int>& rVariable,
    std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <int TDim, int TNumNodes>
std::string EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "EmbeddedTransonicPerturbationPotentialFlowElement #" << this->Id();
    return buffer.str();
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "EmbeddedTransonicPerturbationPotentialFlowElement #" << this->Id();
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::PrintData(std::ostream& rOStream) const
{
    this->pGetGeometry()->PrintData(rOStream);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateLeftHandSideSubsonicElement(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Calculate shape functions
    ElementalData data;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), data.DN_DX, data.N, data.vol);

    const array_1d<double, TDim> velocity = PotentialFlowUtilities::ComputePerturbedVelocity<TDim,TNumNodes>(*this, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TNumNodes> current_lhs = ZeroMatrix(TNumNodes,TNumNodes);

    CalculateLeftHandSideContribution(current_lhs, rCurrentProcessInfo, velocity, data);

    for (int i = 0; i < TNumNodes; i++)
    {
        for (int j = 0; j < TNumNodes; j++)
        {
            rLeftHandSideMatrix(i, j) = current_lhs(i, j);
        }
    }
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateRightHandSideNormalElement(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const EmbeddedTransonicPerturbationPotentialFlowElement& r_this = *this;

    const array_1d<double, TDim> velocity = PotentialFlowUtilities::ComputePerturbedVelocity<TDim,TNumNodes>(r_this, rCurrentProcessInfo);

    double density = 0.0;

    if(r_this.IsNot(INLET)) {
        if (rRightHandSideVector.size() != TNumNodes + 1) {
            rRightHandSideVector.resize(TNumNodes + 1, false);
        }
        rRightHandSideVector.clear();
        const array_1d<double, TDim> upwind_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<TDim,TNumNodes>(*this->pGetUpwindElement(), rCurrentProcessInfo);
        density = PotentialFlowUtilities::ComputeUpwindedDensity<TDim, TNumNodes>(velocity, upwind_velocity, rCurrentProcessInfo);

    } else {
        if (rRightHandSideVector.size() != TNumNodes) {
            rRightHandSideVector.resize(TNumNodes, false);
        }
        rRightHandSideVector.clear();

        const double local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<TDim, TNumNodes>(velocity, rCurrentProcessInfo);
        density = PotentialFlowUtilities::ComputeDensity<TDim, TNumNodes>(local_mach_number_squared, rCurrentProcessInfo);
    }
    Vector distances(TNumNodes);
    for(unsigned int i_node = 0; i_node<TNumNodes; i_node++)
        distances(i_node) = this->GetGeometry()[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);

    ModifiedShapeFunctions::Pointer pModifiedShFunc = this->pGetModifiedShapeFunctions(distances);
    Matrix positive_side_sh_func;
    ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients;
    Vector positive_side_weights;
    pModifiedShFunc -> ComputePositiveSideShapeFunctionsAndGradientsValues(
        positive_side_sh_func,
        positive_side_sh_func_gradients,
        positive_side_weights,
        GeometryData::GI_GAUSS_2);
    BoundedVector<double, TNumNodes> current_rhs = ZeroVector(TNumNodes);
    for (unsigned int i_gauss=0;i_gauss<positive_side_sh_func_gradients.size();i_gauss++){
        BoundedMatrix<double, TNumNodes, TDim> DN_DX = positive_side_sh_func_gradients(i_gauss);
        current_rhs += - positive_side_weights(i_gauss) * density * prod(DN_DX, velocity);
    }

    for (int i = 0; i < TNumNodes; i++)
    {
        rRightHandSideVector[i] = current_rhs[i];
    }
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateLeftHandSideNormalElement(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    const EmbeddedTransonicPerturbationPotentialFlowElement& r_this = *this;

    const array_1d<double, TDim> velocity = PotentialFlowUtilities::ComputePerturbedVelocity<TDim,TNumNodes>(r_this, rCurrentProcessInfo);
    const array_1d<double, TDim> upwind_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<TDim,TNumNodes>(*this->pGetUpwindElement(), rCurrentProcessInfo);

    const double local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<TDim, TNumNodes>(velocity, rCurrentProcessInfo);
    const double upwind_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<TDim, TNumNodes>(upwind_velocity, rCurrentProcessInfo);

    double DrhoDu2 = 0.0;
    double DrhoDu2_up = 0.0;

    const double critical_mach_sq = std::pow(rCurrentProcessInfo[CRITICAL_MACH], 2.0);

    // To check clamping
    const double max_velocity_squared = PotentialFlowUtilities::ComputeMaximumVelocitySquared<TDim, TNumNodes>(rCurrentProcessInfo);
    const double local_velocity_squared = inner_prod(velocity, velocity);
    const double upwind_velocity_squared = inner_prod(upwind_velocity, upwind_velocity);
    // KRATOS_WATCH(local_mach_number_squared)
    // KRATOS_WATCH(upwind_mach_number_squared)
    if (local_mach_number_squared < critical_mach_sq or this->pGetUpwindElement()->GetValue(WAKE)) { // subsonic, not inlet
        // gets [TNumNodes + 1, TNumNodes + 1] size matrix

        if (this->Id()==80000000) {
                KRATOS_WATCH(this->Id())
                KRATOS_WATCH("SUBSONIC LEFT")
            }
        CalculateLeftHandSideSubsonicElement(rLeftHandSideMatrix, rCurrentProcessInfo);
        return;
    }
    else if (local_mach_number_squared >= upwind_mach_number_squared) { // supersonic, accelerating
        // density derivatives
        if (local_velocity_squared < max_velocity_squared){
            DrhoDu2 = PotentialFlowUtilities::ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicAccelerating<TDim, TNumNodes>(
            velocity, local_mach_number_squared, upwind_mach_number_squared, rCurrentProcessInfo);
        }
        if (upwind_velocity_squared < max_velocity_squared){
            DrhoDu2_up = PotentialFlowUtilities::ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicAccelerating<TDim, TNumNodes>(
            local_mach_number_squared, upwind_mach_number_squared, rCurrentProcessInfo);
        }
    }
    else { // supersonic, deaccelerating
        // density derivatives
        if (local_velocity_squared < max_velocity_squared){
            DrhoDu2 = PotentialFlowUtilities::ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicDeaccelerating<TDim, TNumNodes>(
            local_mach_number_squared, upwind_mach_number_squared, rCurrentProcessInfo);
        }
        if (upwind_velocity_squared < max_velocity_squared){
            DrhoDu2_up = PotentialFlowUtilities::ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicDeaccelerating<TDim, TNumNodes>(
            upwind_velocity, local_mach_number_squared, upwind_mach_number_squared, rCurrentProcessInfo);
        }
    }

    // KRATOS_WATCH("assemble super")
    AssembleSupersonicLeftHandSide(rLeftHandSideMatrix, DrhoDu2, DrhoDu2_up, velocity, upwind_velocity, rCurrentProcessInfo);
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateLeftHandSideWakeElement(
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

    ElementalData data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), data.DN_DX, data.N, data.vol);
    data.distances = PotentialFlowUtilities::GetWakeDistances<TDim, TNumNodes>(*this);

    // Compute upper and lower velocities
    const array_1d<double, TDim> upper_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<TDim,TNumNodes>(*this, rCurrentProcessInfo);
    const array_1d<double, TDim> lower_velocity = PotentialFlowUtilities::ComputePerturbedVelocityLowerElement<TDim,TNumNodes>(*this, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TNumNodes> upper_lhs_total = ZeroMatrix(TNumNodes,TNumNodes);
    BoundedMatrix<double, TNumNodes, TNumNodes> lower_lhs_total = ZeroMatrix(TNumNodes,TNumNodes);

    CalculateLeftHandSideContribution(upper_lhs_total, rCurrentProcessInfo, upper_velocity, data);
    CalculateLeftHandSideContribution(lower_lhs_total, rCurrentProcessInfo, lower_velocity, data);

    // Compute lhs wake condition
    const BoundedMatrix<double, TNumNodes, TNumNodes> lhs_wake_condition =
        CalculateLeftHandSideWakeConditions(data, rCurrentProcessInfo);

    if (this->Is(STRUCTURE))
    {
        Matrix lhs_positive = ZeroMatrix(TNumNodes, TNumNodes);
        Matrix lhs_negative = ZeroMatrix(TNumNodes, TNumNodes);

        CalculateLeftHandSideSubdividedElement(lhs_positive, lhs_negative, rCurrentProcessInfo);
        AssignLeftHandSideSubdividedElement(
            rLeftHandSideMatrix, lhs_positive, lhs_negative, upper_lhs_total,
            lower_lhs_total, lhs_wake_condition, data);
    }
    else
    {
        AssignLeftHandSideWakeElement(rLeftHandSideMatrix, upper_lhs_total, lower_lhs_total, lhs_wake_condition, data);
    }
}

// In 2D
template <>
BoundedMatrix<double, 3, 3> EmbeddedTransonicPerturbationPotentialFlowElement<2, 3>::CalculateLeftHandSideWakeConditions(
    const ElementalData& rData, const ProcessInfo& rCurrentProcessInfo)
{
    const double free_stream_density = rCurrentProcessInfo[FREE_STREAM_DENSITY];
    return rData.vol * free_stream_density * prod(rData.DN_DX, trans(rData.DN_DX));
}

// In 3D
template <>
BoundedMatrix<double, 4, 4> EmbeddedTransonicPerturbationPotentialFlowElement<3, 4>::CalculateLeftHandSideWakeConditions(
    const ElementalData& rData, const ProcessInfo& rCurrentProcessInfo)
{
    // Computing linearized pressure equality condition lhs
    const array_1d<double, 3>& free_stream_velocity_direction =
        rCurrentProcessInfo[FREE_STREAM_VELOCITY_DIRECTION];
    const BoundedVector<double, 4> DNv = prod(rData.DN_DX, free_stream_velocity_direction);
    const BoundedMatrix<double, 4, 4> pressure_equality_lhs =
        outer_prod(DNv, trans(DNv));

    // Computing wake normal condition lhs
    // Attention: this only works for straight trailing edges
    // TODO: Make it work for curved trailing edges, i.e., find the way to store
    // the local normal vector of the skin in the element.
    const array_1d<double, 3>& wake_normal = rCurrentProcessInfo[WAKE_NORMAL];
    const BoundedVector<double, 4> DNn = prod(rData.DN_DX, wake_normal);
    const BoundedMatrix<double, 4, 4> normal_condition_lhs = outer_prod(DNn, trans(DNn));

    // Adding contributions
    return rData.vol * (pressure_equality_lhs + normal_condition_lhs);
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateRightHandSideWakeElement(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
KRATOS_ERROR << "Computing wake RHS from embedded";
}

// In 2D
template <>
BoundedVector<double, 3> EmbeddedTransonicPerturbationPotentialFlowElement<2, 3>::CalculateRightHandSideWakeConditions(
    const ElementalData& rData,
    const ProcessInfo& rCurrentProcessInfo,
    const array_1d<double, 2>& rDiff_velocity)
{
    const double free_stream_density = rCurrentProcessInfo[FREE_STREAM_DENSITY];
    return -rData.vol * free_stream_density * prod(rData.DN_DX, rDiff_velocity);
}

// In 3D
template <>
BoundedVector<double, 4> EmbeddedTransonicPerturbationPotentialFlowElement<3, 4>::CalculateRightHandSideWakeConditions(
    const ElementalData& rData,
    const ProcessInfo& rCurrentProcessInfo,
    const array_1d<double, 3>& rDiff_velocity)
{
    // Computing linearized pressure equality condition
    const array_1d<double, 3>& free_stream_velocity_direction =
        rCurrentProcessInfo[FREE_STREAM_VELOCITY_DIRECTION];
    const double pressure_equality_condition =
        inner_prod(free_stream_velocity_direction, rDiff_velocity);
    const array_1d<double, 3> vp = free_stream_velocity_direction * pressure_equality_condition;

    // Computing wake normal flux condition
    // Attention: this only works for straight trailing edges
    // TODO: Make it work for curved trailing edges, i.e., find the way to store
    // the local normal vector of the skin in the element.
    const array_1d<double, 3>& wake_normal = rCurrentProcessInfo[WAKE_NORMAL];
    const double wake_normal_condition = inner_prod(wake_normal, rDiff_velocity);
    const array_1d<double, 3> nn = wake_normal * wake_normal_condition;

    return -rData.vol * prod(rData.DN_DX, vp + nn);
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateLeftHandSideContribution(
    BoundedMatrix<double, TNumNodes, TNumNodes>& rLhs_total,
    const ProcessInfo& rCurrentProcessInfo,
    const array_1d<double, TDim>& rVelocity,
    const ElementalData& rData)
{
    // Embedded terms
    Vector distances(TNumNodes);
    for(unsigned int i_node = 0; i_node<TNumNodes; i_node++)
        distances(i_node) = this->GetGeometry()[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
    const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<TDim,TNumNodes>(distances);
    if(this->Id() ==80000000)
    {
        KRATOS_WATCH(is_embedded)
    }
    if (!is_embedded) {
        BaseType::CalculateLeftHandSideContribution(rLhs_total,
            rCurrentProcessInfo,
            rVelocity,
            rData);
        return;
    }
    // Compute density
    const double local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<TDim, TNumNodes>(rVelocity, rCurrentProcessInfo);
    const double density = PotentialFlowUtilities::ComputeDensity<TDim, TNumNodes>(local_mach_number_squared, rCurrentProcessInfo);

    // Compute density derivative
    const double DrhoDu2 = PotentialFlowUtilities::ComputeDensityDerivativeWRTVelocitySquared<TDim, TNumNodes>(local_mach_number_squared, rCurrentProcessInfo);

    // Embedded terms
    ModifiedShapeFunctions::Pointer pModifiedShFunc = this->pGetModifiedShapeFunctions(distances);
    Matrix positive_side_sh_func;
    ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients;
    Vector positive_side_weights;
    pModifiedShFunc -> ComputePositiveSideShapeFunctionsAndGradientsValues(
        positive_side_sh_func,
        positive_side_sh_func_gradients,
        positive_side_weights,
        GeometryData::GI_GAUSS_2);

    // Compute lhs
    const double local_velocity_squared = inner_prod(rVelocity, rVelocity);
    const double max_velocity_squared = PotentialFlowUtilities::ComputeMaximumVelocitySquared<TDim, TNumNodes>(rCurrentProcessInfo);
    for (unsigned int i_gauss=0;i_gauss<positive_side_sh_func_gradients.size();i_gauss++){
        BoundedMatrix<double,TNumNodes,TDim> DN_DX = positive_side_sh_func_gradients(i_gauss);
        const BoundedVector<double, TNumNodes> DNV = prod(DN_DX, rVelocity);

        rLhs_total += positive_side_weights(i_gauss) * density * prod(DN_DX, trans(DN_DX));


        if (local_velocity_squared < max_velocity_squared){
            rLhs_total += positive_side_weights(i_gauss) * 2 * DrhoDu2 * outer_prod(DNV, trans(DNV));
        }
    }
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateLeftHandSideSubdividedElement(
    Matrix& lhs_positive,
    Matrix& lhs_negative,
    const ProcessInfo& rCurrentProcessInfo)
{
    ElementalData data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), data.DN_DX, data.N, data.vol);

    data.distances = PotentialFlowUtilities::GetWakeDistances<TDim, TNumNodes>(*this);

    // Subdivide the element
    constexpr unsigned int nvolumes = 3 * (TDim - 1);
    BoundedMatrix<double, TNumNodes, TDim> Points;
    array_1d<double, nvolumes> PartitionsSign;
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
        const auto& r_coords = this->GetGeometry()[i].Coordinates();
        for (unsigned int k = 0; k < TDim; ++k)
        {
            Points(i, k) = r_coords[k];
        }
    }

    const unsigned int nsubdivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(
        Points, data.DN_DX, data.distances, Volumes, gp_shape_function_values,
        PartitionsSign, GradientsValue, N_enriched);

    // Compute upper and lower velocities
    const array_1d<double, TDim> upper_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<TDim,TNumNodes>(*this, rCurrentProcessInfo);
    const array_1d<double, TDim> lower_velocity = PotentialFlowUtilities::ComputePerturbedVelocityLowerElement<TDim,TNumNodes>(*this, rCurrentProcessInfo);

    // Compute upper and lower densities
    const double upper_local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<TDim, TNumNodes>(upper_velocity, rCurrentProcessInfo);
    const double upper_density = PotentialFlowUtilities::ComputeDensity<TDim, TNumNodes>(upper_local_mach_number_squared, rCurrentProcessInfo);

    const double lower_local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<TDim, TNumNodes>(lower_velocity, rCurrentProcessInfo);
    const double lower_density = PotentialFlowUtilities::ComputeDensity<TDim, TNumNodes>(lower_local_mach_number_squared, rCurrentProcessInfo);

    // Compute upper and lower density derivatives
    const double upper_DrhoDu2 = PotentialFlowUtilities::ComputeDensityDerivativeWRTVelocitySquared<TDim, TNumNodes>(upper_local_mach_number_squared, rCurrentProcessInfo);
    const double lower_DrhoDu2 = PotentialFlowUtilities::ComputeDensityDerivativeWRTVelocitySquared<TDim, TNumNodes>(lower_local_mach_number_squared, rCurrentProcessInfo);

    // Compute upper and lower lhs
    const BoundedVector<double, TNumNodes> upper_DNV = prod(data.DN_DX, upper_velocity);
    const BoundedVector<double, TNumNodes> lower_DNV = prod(data.DN_DX, lower_velocity);

    const double upper_local_velocity_squared = inner_prod(upper_velocity, upper_velocity);
    const double lower_local_velocity_squared = inner_prod(lower_velocity, lower_velocity);
    const double max_velocity_squared = PotentialFlowUtilities::ComputeMaximumVelocitySquared<TDim, TNumNodes>(rCurrentProcessInfo);

    // Compute the lhs and rhs that would correspond to it being divided
    for (unsigned int i = 0; i < nsubdivisions; ++i)
    {
        if (PartitionsSign[i] > 0)
        {
            noalias(lhs_positive) += Volumes[i] * upper_density * prod(data.DN_DX, trans(data.DN_DX));
            if(upper_local_velocity_squared < max_velocity_squared) {
                noalias(lhs_positive) += Volumes[i] * 2 * upper_DrhoDu2 * outer_prod(upper_DNV, trans(upper_DNV));
            }

        }
        else
        {
            noalias(lhs_negative) += Volumes[i] * lower_density * prod(data.DN_DX, trans(data.DN_DX));
            if(lower_local_velocity_squared < max_velocity_squared) {
                noalias(lhs_negative) += Volumes[i] * 2 * lower_DrhoDu2 * outer_prod(lower_DNV, trans(lower_DNV));
            }
        }
    }
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateVolumesSubdividedElement(
    double& rUpper_vol,
    double& rLower_vol,
    const ProcessInfo& rCurrentProcessInfo)
{
    ElementalData data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), data.DN_DX, data.N, data.vol);

    data.distances = PotentialFlowUtilities::GetWakeDistances<TDim, TNumNodes>(*this);

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
        const array_1d<double, 3>& coords = this->GetGeometry()[i].Coordinates();
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
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::ComputeLHSGaussPointContribution(
    const double weight,
    Matrix& lhs,
    const ElementalData& data) const
{
    noalias(lhs) += weight * prod(data.DN_DX, trans(data.DN_DX));
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::AssignLeftHandSideSubdividedElement(
    Matrix& rLeftHandSideMatrix,
    Matrix& lhs_positive,
    Matrix& lhs_negative,
    const BoundedMatrix<double, TNumNodes, TNumNodes>& rUpper_lhs_total,
    const BoundedMatrix<double, TNumNodes, TNumNodes>& rLower_lhs_total,
    const BoundedMatrix<double, TNumNodes, TNumNodes>& rLhs_wake_condition,
    const ElementalData& data) const
{
    const auto& r_geometry = this->GetGeometry();
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        // The TE node takes the contribution of the subdivided element and
        // we do not apply the wake condition on the TE node
        if (r_geometry[i].GetValue(TRAILING_EDGE)) {
            for (unsigned int j = 0; j < TNumNodes; ++j) {
                rLeftHandSideMatrix(i, j) = lhs_positive(i, j);
                rLeftHandSideMatrix(i + TNumNodes, j + TNumNodes) = lhs_negative(i, j);
            }
        }
        else{
            AssignLeftHandSideWakeNode(rLeftHandSideMatrix, rUpper_lhs_total, rLower_lhs_total, rLhs_wake_condition, data, i);
        }
    }
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::AssignLeftHandSideWakeElement(
    MatrixType& rLeftHandSideMatrix,
    const BoundedMatrix<double, TNumNodes, TNumNodes>& rUpper_lhs_total,
    const BoundedMatrix<double, TNumNodes, TNumNodes>& rLower_lhs_total,
    const BoundedMatrix<double, TNumNodes, TNumNodes>& rLhs_wake_condition,
    const ElementalData& rData) const
{
    for (unsigned int row = 0; row < TNumNodes; ++row){
        AssignLeftHandSideWakeNode(rLeftHandSideMatrix, rUpper_lhs_total,
                                   rLower_lhs_total, rLhs_wake_condition, rData, row);
    }
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::AssignLeftHandSideWakeNode(
    MatrixType& rLeftHandSideMatrix,
    const BoundedMatrix<double, TNumNodes, TNumNodes>& rUpper_lhs_total,
    const BoundedMatrix<double, TNumNodes, TNumNodes>& rLower_lhs_total,
    const BoundedMatrix<double, TNumNodes, TNumNodes>& rLhs_wake_condition,
    const ElementalData& rData,
    unsigned int row) const
{
    // Applying wake condition on the AUXILIARY_VELOCITY_POTENTIAL dofs
    if (rData.distances[row] < 0.0){
        for (unsigned int column = 0; column < TNumNodes; ++column){
            // Conservation of mass
            rLeftHandSideMatrix(row + TNumNodes, column + TNumNodes) = rLower_lhs_total(row, column);
            // Wake condition
            rLeftHandSideMatrix(row, column) = rLhs_wake_condition(row, column); // Diagonal
            rLeftHandSideMatrix(row, column + TNumNodes) = -rLhs_wake_condition(row, column); // Off diagonal
        }
    }
    else{ // else if (data.distances[row] > 0.0)
        for (unsigned int column = 0; column < TNumNodes; ++column){
            // Conservation of mass
            rLeftHandSideMatrix(row, column) = rUpper_lhs_total(row, column);
            // Wake condition
            rLeftHandSideMatrix(row + TNumNodes, column + TNumNodes) = rLhs_wake_condition(row, column); // Diagonal
            rLeftHandSideMatrix(row + TNumNodes, column) = -rLhs_wake_condition(row, column); // Off diagonal
        }
    }
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::AssignRightHandSideWakeNode(
    VectorType& rRightHandSideVector,
    const BoundedVector<double, TNumNodes>& rUpper_rhs,
    const BoundedVector<double, TNumNodes>& rLower_rhs,
    const BoundedVector<double, TNumNodes>& rWake_rhs,
    const ElementalData& rData,
    unsigned int& rRow) const
{
    if (rData.distances[rRow] > 0.0)
    {
        rRightHandSideVector[rRow] = rUpper_rhs(rRow);
        rRightHandSideVector[rRow + TNumNodes] = -rWake_rhs(rRow);
    }
    else
    {
        rRightHandSideVector[rRow] = rWake_rhs(rRow);
        rRightHandSideVector[rRow + TNumNodes] = rLower_rhs(rRow);
    }
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::AssembleSupersonicLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const double densityDerivativeWRTVelocity,
    const double densityDerivativeWRTUpwindVelocity,
    const array_1d<double, TDim> velocity,
    const array_1d<double, TDim> upwindVelocity,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (this->Id()==80000000) KRATOS_WATCH("ASSEMBLING EMB LEFT HAND SIDE BASE Normal element")

    Vector distances(TNumNodes);
    for(unsigned int i_node = 0; i_node<TNumNodes; i_node++)
        distances(i_node) = this->GetGeometry()[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
    const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<TDim,TNumNodes>(distances);
    if (!is_embedded) {
        BaseType::AssembleSupersonicLeftHandSide(rLeftHandSideMatrix,densityDerivativeWRTVelocity,
        densityDerivativeWRTUpwindVelocity,
        velocity,
        upwindVelocity,
        rCurrentProcessInfo);
        return;
    }

    // KRATOS_WATCH(rLeftHandSideMatrix)
    // KRATOS_WATCH("assemble density")
    BoundedVector<double, TNumNodes + 1> DNV_assembly = this->AssembleDensityDerivativeAndShapeFunctions(
        densityDerivativeWRTVelocity, densityDerivativeWRTUpwindVelocity, velocity, upwindVelocity, rCurrentProcessInfo);

    // Calculate shape functions
    ElementalData data;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), data.DN_DX, data.N, data.vol);
    // KRATOS_WATCH("compute upwdind density")

    const double density = PotentialFlowUtilities::ComputeUpwindedDensity<TDim, TNumNodes>(velocity, upwindVelocity, rCurrentProcessInfo);

    const BoundedVector<double, TNumNodes> DNV = prod(data.DN_DX, velocity);
    BoundedVector<double, TNumNodes + 1> DNV_extended;
    DNV_extended.clear();
    for (int i = 0; i < TNumNodes; i++) {
        DNV_extended[i] = DNV[i];
    }
    // Embedded terms

    ModifiedShapeFunctions::Pointer pModifiedShFunc = this->pGetModifiedShapeFunctions(distances);
    Matrix positive_side_sh_func;
    ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients;
    Vector positive_side_weights;
    pModifiedShFunc -> ComputePositiveSideShapeFunctionsAndGradientsValues(
        positive_side_sh_func,
        positive_side_sh_func_gradients,
        positive_side_weights,
        GeometryData::GI_GAUSS_2);
    BoundedMatrix<double, TNumNodes, TNumNodes> linear_term = ZeroMatrix(TNumNodes, TNumNodes);

    BoundedMatrix<double,TNumNodes,TDim> DN_DX;
    for (unsigned int i_gauss=0;i_gauss<positive_side_sh_func_gradients.size();i_gauss++){
        DN_DX = positive_side_sh_func_gradients(i_gauss);
        linear_term += positive_side_weights(i_gauss) * density * prod(DN_DX, trans(DN_DX));

        rLeftHandSideMatrix += positive_side_weights(i_gauss) * 2.0 * outer_prod(DNV_extended, trans(DNV_assembly));
    }

    for (int i = 0; i < TNumNodes; i++)
    {
        for (int j = 0; j < TNumNodes; j++)
        {
            rLeftHandSideMatrix(i, j) += linear_term(i, j);
        }
    }
}

template <>
ModifiedShapeFunctions::Pointer EmbeddedTransonicPerturbationPotentialFlowElement<2,3>::pGetModifiedShapeFunctions(Vector& rDistances) {
    return Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(this->pGetGeometry(), rDistances);
}

template <>
ModifiedShapeFunctions::Pointer EmbeddedTransonicPerturbationPotentialFlowElement<3,4>::pGetModifiedShapeFunctions(Vector& rDistances) {
    return Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(this->pGetGeometry(), rDistances);
}


// serializer

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class EmbeddedTransonicPerturbationPotentialFlowElement<2, 3>;
template class EmbeddedTransonicPerturbationPotentialFlowElement<3, 4>;

} // namespace Kratos
