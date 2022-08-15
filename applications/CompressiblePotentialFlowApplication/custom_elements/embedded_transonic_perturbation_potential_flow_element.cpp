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
    KRATOS_TRY
    BaseType::Initialize(rCurrentProcessInfo);
    KRATOS_CATCH("");
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
        BaseType::CalculateRightHandSideNormalElement(rRightHandSideVector,rCurrentProcessInfo);
    }
    else {
        if (wake==0) {
            BaseType::CalculateRightHandSideNormalElement(rRightHandSideVector, rCurrentProcessInfo);
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
        PotentialFlowUtilities::AddKuttaConditionPenaltyPerturbationRHS<TDim,TNumNodes>(r_this,rRightHandSideVector,rCurrentProcessInfo);
    }
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    const EmbeddedTransonicPerturbationPotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) { // Normal element (non-wake)
        if (r_this.IsNot(INLET)) {
            if (rLeftHandSideMatrix.size1() != TNumNodes + 1 ||
                rLeftHandSideMatrix.size2() != TNumNodes + 1) {
                    rLeftHandSideMatrix.resize(TNumNodes + 1, TNumNodes + 1, false);
                }

            rLeftHandSideMatrix.clear();
            BaseType::CalculateLeftHandSideNormalElement(rLeftHandSideMatrix, rCurrentProcessInfo);
        } else {
            if (rLeftHandSideMatrix.size1() != TNumNodes ||
                rLeftHandSideMatrix.size2() != TNumNodes) {
                    rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
                }

            rLeftHandSideMatrix.clear();
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
        PotentialFlowUtilities::AddKuttaConditionPenaltyPerturbationLHS<TDim,TNumNodes>(r_this,rLeftHandSideMatrix,rCurrentProcessInfo);
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

    const auto& r_geometry = this->GetGeometry();
    ElementalData data{r_geometry};

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

    BoundedVector<double, TNumNodes> upper_rhs;
    BoundedVector<double, TNumNodes> lower_rhs;
    CalculateRightHandSideContribution(upper_rhs, density_upper, upper_velocity);
    CalculateRightHandSideContribution(lower_rhs, density_lower, lower_velocity);

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

    ElementalData data{this->GetGeometry()};

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

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::CalculateRightHandSideContribution(
    BoundedVector<double, TNumNodes>& rRhs_total,
    const double rDensity,
    const array_1d<double, TDim>& rVelocity)
{
    Vector distances(TNumNodes);
    for(unsigned int i_node = 0; i_node<TNumNodes; i_node++)
        distances(i_node) = this->GetGeometry()[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);

    const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<TDim,TNumNodes>(distances);
    if (!is_embedded) {
        BaseType::CalculateRightHandSideContribution(rRhs_total,
            rDensity,
            rVelocity);
        return;
    }

    rRhs_total.clear();
    ModifiedShapeFunctions::Pointer pModifiedShFunc = this->pGetModifiedShapeFunctions(distances);
    Matrix positive_side_sh_func;
    ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients;
    Vector positive_side_weights;
    pModifiedShFunc -> ComputePositiveSideShapeFunctionsAndGradientsValues(
        positive_side_sh_func,
        positive_side_sh_func_gradients,
        positive_side_weights,
        GeometryData::IntegrationMethod::GI_GAUSS_2);

    for (unsigned int i_gauss=0;i_gauss<positive_side_sh_func_gradients.size();i_gauss++){
        BoundedMatrix<double, TNumNodes, TDim> DN_DX = positive_side_sh_func_gradients(i_gauss);
        rRhs_total += - positive_side_weights(i_gauss) * rDensity * prod(DN_DX, rVelocity);
    }
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
        GeometryData::IntegrationMethod::GI_GAUSS_2);

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
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::AssembleSupersonicLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const double densityDerivativeWRTVelocity,
    const double densityDerivativeWRTUpwindVelocity,
    const array_1d<double, TDim> velocity,
    const array_1d<double, TDim> upwindVelocity,
    const ProcessInfo& rCurrentProcessInfo)
{
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

    BoundedVector<double, TNumNodes + 1> DNV_assembly = this->AssembleDensityDerivativeAndShapeFunctions(
        densityDerivativeWRTVelocity, densityDerivativeWRTUpwindVelocity, velocity, upwindVelocity, rCurrentProcessInfo);

    // Calculate shape functions
    ElementalData data{this->GetGeometry()};

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
        GeometryData::IntegrationMethod::GI_GAUSS_2);
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

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::FindUpwindElement(const ProcessInfo& rCurrentProcessInfo)
{
    GeometryType upwind_element_boundary=this->GetGeometry();
    std::vector<size_t> upwind_element_nodes;
    GlobalPointersVector<Element> upwind_element_candidates;
    PotentialFlowUtilities::GetNodeNeighborElementCandidates<TDim, TNumNodes>(upwind_element_candidates, upwind_element_boundary);

    SelectUpwindElement(upwind_element_nodes, upwind_element_candidates, rCurrentProcessInfo);
}

template <int TDim, int TNumNodes>
void EmbeddedTransonicPerturbationPotentialFlowElement<TDim, TNumNodes>::SelectUpwindElement(
    std::vector<IndexType>& rUpwindElementNodeIds,
    GlobalPointersVector<Element>& rUpwindElementCandidates,
    const ProcessInfo& rCurrentProcessInfo)
{
    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    double minimum_edge_flux = std::numeric_limits<double>::max();
    bool any_flux_negative = false;
    const auto& element_boundary_geometry = this->GetGeometry().GenerateBoundariesEntities();

    for (SizeType i = 0; i < rUpwindElementCandidates.size(); i++)
    {
        // get sorted node ids of neighbording elements
        std::vector<size_t> neighbor_element_ids;
        PotentialFlowUtilities::GetSortedIds<TDim, TNumNodes>(neighbor_element_ids, rUpwindElementCandidates[i].GetGeometry());

        // find element which shares the upwind element nodes with current element
        // but is not the current element
        // extra check if its active

        if (rUpwindElementCandidates[i].Id()!=this->Id()) {
            for (auto& r_geom_boundary : element_boundary_geometry) {
                std::vector<size_t> this_sorted_ids;
                PotentialFlowUtilities::GetSortedIds<TDim, TNumNodes>(this_sorted_ids, r_geom_boundary);


                if(std::includes(neighbor_element_ids.begin(), neighbor_element_ids.end(),
                this_sorted_ids.begin(), this_sorted_ids.end()))
                {
                    const auto edge_normal = BaseType::GetEdgeNormal(r_geom_boundary);

                    const double edge_flux = inner_prod(edge_normal, free_stream_velocity);

                    if (edge_flux < 0) {
                        any_flux_negative = true;
                    }

                    if(edge_flux < minimum_edge_flux && rUpwindElementCandidates[i].Is(ACTIVE)) {
                        minimum_edge_flux = edge_flux;
                        this->pSetUpwindElement(rUpwindElementCandidates(i));
                    }
                }
            }
        }
    }

    // If no upwind element is found, the element is an INLET element and the
    // upwind element pointer points to itself
    if (this->CheckUpwindElement() || !any_flux_negative)
    {
        this->pSetUpwindElement(this);
        this->SetFlags(INLET);
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
