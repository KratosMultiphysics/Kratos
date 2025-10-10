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

#include "compressible_perturbation_potential_flow_element.h"
#include "compressible_potential_flow_application_variables.h"
#include "fluid_dynamics_application_variables.h"
#include "includes/cfd_variables.h"
#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/potential_flow_utilities.h"

namespace Kratos
{
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <int Dim, int NumNodes>
Element::Pointer CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<CompressiblePerturbationPotentialFlowElement>(
        NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<CompressiblePerturbationPotentialFlowElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::Clone(
    IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<CompressiblePerturbationPotentialFlowElement>(
        NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    FindUpwindElement(rCurrentProcessInfo);
    this->GetValue(UPWIND_CONTRIBUTION) = ZeroMatrix(NumNodes, NumNodes);
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);
    CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    const CompressiblePerturbationPotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) // Normal element (non-wake) - eventually an embedded
        CalculateRightHandSideNormalElement(rRightHandSideVector, rCurrentProcessInfo);
    else // Wake element
        CalculateRightHandSideWakeElement(rRightHandSideVector, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    FindDownwindElement(rCurrentProcessInfo);
    const CompressiblePerturbationPotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    KRATOS_INFO("Compressible Element") << "inicio \n";
    if (wake == 0){ // Normal element (non-wake) - eventually an embedded
        KRATOS_INFO("Compressible Element") << "normal \n";
        CalculateLeftHandSideNormalElement(rLeftHandSideMatrix, rCurrentProcessInfo);
    }else{ // Wake element
        KRATOS_INFO("Compressible Element") << "wake \n";
        CalculateLeftHandSideWakeElement(rLeftHandSideMatrix, rCurrentProcessInfo);
    }
    KRATOS_INFO("Compressible Element") << "pre apuntar puntero \n";
    MatrixType lhs;
    // VectorType rhs;
    // pGetDownwindElement()->CalculateLocalSystem(lhs, rhs, rCurrentProcessInfo);
    pGetDownwindElement()->CalculateLeftHandSide(lhs, rCurrentProcessInfo);
    auto upwind_contribution = pGetDownwindElement()->GetValue(UPWIND_CONTRIBUTION);

    KRATOS_WATCH(pGetDownwindElement()->GetValue(UPWIND_CONTRIBUTION));

    for (int i = 0; i < TNumNodes; i++)
    {
        for (int j = 0; j < TNumNodes; j++)
        {
            rLeftHandSideMatrix(i, j) += upwind_contribution(i, j);
        }
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::EquationIdVector(
    EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const
{
    const CompressiblePerturbationPotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) // Normal element
    {
        if (rResult.size() != NumNodes)
            rResult.resize(NumNodes, false);

        const int kutta = r_this.GetValue(KUTTA);

        if (kutta == 0)
            GetEquationIdVectorNormalElement(rResult);
        else
            GetEquationIdVectorKuttaElement(rResult);
    }
    else // Wake element
    {
        if (rResult.size() != 2 * NumNodes)
            rResult.resize(2 * NumNodes, false);

        GetEquationIdVectorWakeElement(rResult);
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                                                 const ProcessInfo& CurrentProcessInfo) const
{
    const CompressiblePerturbationPotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) // Normal element
    {
        if (rElementalDofList.size() != NumNodes)
            rElementalDofList.resize(NumNodes);

        const int kutta = r_this.GetValue(KUTTA);

        if (kutta == 0)
            GetDofListNormalElement(rElementalDofList);
        else
            GetDofListKuttaElement(rElementalDofList);
    }
    else // wake element
    {
        if (rElementalDofList.size() != 2 * NumNodes)
            rElementalDofList.resize(2 * NumNodes);

        GetDofListWakeElement(rElementalDofList);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <int Dim, int NumNodes>
int CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
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

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == PRESSURE_COEFFICIENT)
    {
        rValues[0] = PotentialFlowUtilities::ComputePerturbationCompressiblePressureCoefficient<Dim, NumNodes>(*this, rCurrentProcessInfo);
    }
    else if (rVariable == DENSITY)
    {
        const array_1d<double, Dim>& velocity = PotentialFlowUtilities::ComputePerturbedVelocity<Dim,NumNodes>(*this, rCurrentProcessInfo);

        const double local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<Dim, NumNodes>(velocity, rCurrentProcessInfo);

        rValues[0] = PotentialFlowUtilities::ComputeDensity<Dim, NumNodes>(local_mach_number_squared, rCurrentProcessInfo);
    }
    else if (rVariable == MACH)
    {
        rValues[0] = PotentialFlowUtilities::ComputePerturbationLocalMachNumber<Dim, NumNodes>(*this, rCurrentProcessInfo);
    }
    else if (rVariable == SOUND_VELOCITY)
    {
        rValues[0] = PotentialFlowUtilities::ComputePerturbationLocalSpeedOfSound<Dim, NumNodes>(*this, rCurrentProcessInfo);
    }
    else if (rVariable == WAKE)
    {
        const CompressiblePerturbationPotentialFlowElement& r_this = *this;
        rValues[0] = r_this.GetValue(WAKE);
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateOnIntegrationPoints(
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
    else if (rVariable == TRAILING_EDGE_ELEMENT)
        rValues[0] = this->GetValue(TRAILING_EDGE_ELEMENT);
    else if (rVariable == DECOUPLED_TRAILING_EDGE_ELEMENT)
        rValues[0] = this->GetValue(DECOUPLED_TRAILING_EDGE_ELEMENT);
    else if (rVariable == ID_DOWNWIND_ELEMENT)
        rValues[0] = this->GetValue(ID_DOWNWIND_ELEMENT);
    else if (rVariable == DOWNWIND_ELEMENT)
        rValues[0] = this->GetValue(DOWNWIND_ELEMENT);
    else if (rVariable == UPWIND_ELEMENT)
        rValues[0] = this->GetValue(UPWIND_ELEMENT);
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);
    if (rVariable == VELOCITY){
        const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
        array_1d<double, 3> v(3, 0.0);
        array_1d<double, Dim> vaux = PotentialFlowUtilities::ComputeVelocity<Dim, NumNodes>(*this);
        for (unsigned int k = 0; k < Dim; k++)
            v[k] = vaux[k] + free_stream_velocity[k];
        rValues[0] = v;
    }
    else if (rVariable == PERTURBATION_VELOCITY)
    {
        array_1d<double, 3> v(3, 0.0);
        array_1d<double, Dim> vaux = PotentialFlowUtilities::ComputeVelocity<Dim,NumNodes>(*this);
        for (unsigned int k = 0; k < Dim; k++)
            v[k] = vaux[k];
        rValues[0] = v;
    }
    else if (rVariable == VECTOR_TO_UPWIND_ELEMENT)
    {
        rValues[0] = pGetUpwindElement()->GetGeometry().Center() - this->GetGeometry().Center();
    }
    else if (rVariable == VECTOR_TO_DOWNWIND_ELEMENT)
    {
        rValues[0] = pGetDownwindElement()->GetGeometry().Center() - this->GetGeometry().Center();
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <int Dim, int NumNodes>
std::string CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "CompressiblePerturbationPotentialFlowElement #" << Id();
    return buffer.str();
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "CompressiblePerturbationPotentialFlowElement #" << Id();
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::PrintData(std::ostream& rOStream) const
{
    pGetGeometry()->PrintData(rOStream);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template <int Dim, int NumNodes>
inline GlobalPointer<Element> CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::pGetUpwindElement() const
{
    KRATOS_ERROR_IF(mpUpwindElement.get() == nullptr)
        << "No upwind element found for element #" << this->Id() << std::endl;
    return mpUpwindElement;
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::pSetUpwindElement(GlobalPointer<Element> pUpwindElement)
{
    mpUpwindElement = pUpwindElement;
}

template <int Dim, int NumNodes>
bool CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CheckUpwindElement()
{
    return mpUpwindElement.get() == nullptr;
}

template <int Dim, int NumNodes>
inline GlobalPointer<Element> CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::pGetDownwindElement() const
{
    KRATOS_ERROR_IF(mpDownwindElement.get() == nullptr)
        << "No Downwind element found for element #" << this->Id() << std::endl;
    return mpDownwindElement;
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::pSetDownwindElement(GlobalPointer<Element> pDownwindElement)
{
    mpDownwindElement = pDownwindElement;
}

template <int Dim, int NumNodes>
bool CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CheckDownwindElement()
{
    return mpDownwindElement.get() == nullptr;
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::GetWakeDistances(
    array_1d<double, NumNodes>& distances) const
{
    noalias(distances) = GetValue(WAKE_ELEMENTAL_DISTANCES);
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::GetEquationIdVectorNormalElement(
    EquationIdVectorType& rResult) const
{
    for (unsigned int i = 0; i < NumNodes; i++)
        rResult[i] = GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::GetEquationIdVectorKuttaElement(
    EquationIdVectorType& rResult) const
{
    const auto& r_geometry = this->GetGeometry();
    // Kutta elements have only negative part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (!r_geometry[i].GetValue(TRAILING_EDGE))
            rResult[i] = r_geometry[i].GetDof(VELOCITY_POTENTIAL).EquationId();
        else
            rResult[i] = r_geometry[i].GetDof(AUXILIARY_VELOCITY_POTENTIAL).EquationId();
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::GetEquationIdVectorWakeElement(
    EquationIdVectorType& rResult) const
{
    array_1d<double, NumNodes> distances;
    GetWakeDistances(distances);

    // Positive part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] > 0.0)
            rResult[i] = GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
        else
            rResult[i] =
                GetGeometry()[i].GetDof(AUXILIARY_VELOCITY_POTENTIAL, 0).EquationId();
    }

    // Negative part - sign is opposite to the previous case
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
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::GetDofListNormalElement(DofsVectorType& rElementalDofList) const
{
    for (unsigned int i = 0; i < NumNodes; i++)
        rElementalDofList[i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::GetDofListKuttaElement(DofsVectorType& rElementalDofList) const
{
    const auto& r_geometry = this->GetGeometry();
    // Kutta elements have only negative part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (!r_geometry[i].GetValue(TRAILING_EDGE))
            rElementalDofList[i] = r_geometry[i].pGetDof(VELOCITY_POTENTIAL);
        else
            rElementalDofList[i] = r_geometry[i].pGetDof(AUXILIARY_VELOCITY_POTENTIAL);
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::GetDofListWakeElement(DofsVectorType& rElementalDofList) const
{
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

    // Negative part - sign is opposite to the previous case
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
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateLeftHandSideNormalElement(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    const array_1d<double, Dim> velocity = PotentialFlowUtilities::ComputePerturbedVelocity<Dim,NumNodes>(*this, rCurrentProcessInfo);
    const double local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<Dim, NumNodes>(velocity, rCurrentProcessInfo);
    const double critical_mach_sq = std::pow(rCurrentProcessInfo[CRITICAL_MACH], 2.0);

    if (local_mach_number_squared < critical_mach_sq){ // subsonic
        KRATOS_INFO("Compressible Element") << "subsonic \n";
        CalculateLeftHandSideSubsonicElement(rLeftHandSideMatrix, rCurrentProcessInfo);
    }else {// supersonic
        KRATOS_INFO("Compressible Element") << "supersonic \n";
        CalculateLeftHandSideSupersonicElement(rLeftHandSideMatrix, rCurrentProcessInfo);
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateLeftHandSideSubsonicElement(
    MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != NumNodes || rLeftHandSideMatrix.size2() != NumNodes)
        rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);
    rLeftHandSideMatrix.clear();

    ElementalData data{GetGeometry()};

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    const array_1d<double, Dim>& velocity = PotentialFlowUtilities::ComputePerturbedVelocity<Dim,NumNodes>(*this, rCurrentProcessInfo);

    BoundedMatrix<double, NumNodes, NumNodes> lhs = ZeroMatrix(NumNodes,NumNodes);

    CalculateLeftHandSideContribution(lhs, rCurrentProcessInfo, velocity, data);

    noalias(rLeftHandSideMatrix) = lhs;
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateLeftHandSideSupersonicElement(
    MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    const array_1d<double, Dim> velocity = PotentialFlowUtilities::ComputePerturbedVelocity<Dim,NumNodes>(*this, rCurrentProcessInfo);
    const double local_velocity_squared = inner_prod(velocity, velocity);
    const double local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<Dim, NumNodes>(velocity, rCurrentProcessInfo);

    const array_1d<double, Dim> upwind_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<Dim,NumNodes>(*pGetUpwindElement(), rCurrentProcessInfo);
    const double upwind_velocity_squared = inner_prod(upwind_velocity, upwind_velocity);
    const double upwind_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<Dim, NumNodes>(upwind_velocity, rCurrentProcessInfo);

    const double max_velocity_squared = PotentialFlowUtilities::ComputeMaximumVelocitySquared<Dim, NumNodes>(rCurrentProcessInfo);

    double DrhoDu2 = 0.0;
    double DrhoDu2_up = 0.0;

    if (local_mach_number_squared >= upwind_mach_number_squared) { // supersonic, accelerating
        // density derivatives
        if (local_velocity_squared < max_velocity_squared){
            DrhoDu2 = PotentialFlowUtilities::ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicAccelerating<Dim, NumNodes>(
            velocity, local_mach_number_squared, upwind_mach_number_squared, rCurrentProcessInfo);
        }
        if (upwind_velocity_squared < max_velocity_squared){
            DrhoDu2_up = PotentialFlowUtilities::ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicAccelerating<Dim, NumNodes>(
            local_mach_number_squared, upwind_mach_number_squared, rCurrentProcessInfo);
        }
    }
    else { // supersonic, deaccelerating
        // density derivatives
        if (local_velocity_squared < max_velocity_squared){
            DrhoDu2 = PotentialFlowUtilities::ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicDeaccelerating<Dim, NumNodes>(
            local_mach_number_squared, upwind_mach_number_squared, rCurrentProcessInfo);
        }
        if (upwind_velocity_squared < max_velocity_squared){
            DrhoDu2_up = PotentialFlowUtilities::ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicDeaccelerating<Dim, NumNodes>(
            upwind_velocity, local_mach_number_squared, upwind_mach_number_squared, rCurrentProcessInfo);
        }
    }

    AssembleSupersonicLeftHandSide(rLeftHandSideMatrix, DrhoDu2, DrhoDu2_up, velocity, upwind_velocity, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::AssembleSupersonicLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const double densityDerivativeWRTVelocity,
    const double densityDerivativeWRTUpwindVelocity,
    const array_1d<double, Dim> velocity,
    const array_1d<double, Dim> upwindVelocity,
    const ProcessInfo& rCurrentProcessInfo)
{
    pGetUpwindElement()->GetValue(UPWIND_ELEMENT) = true;
    this->GetValue(DOWNWIND_ELEMENT) = true;

    // Current element data
    ElementalData current_data{GetGeometry()};
    const double density = PotentialFlowUtilities::ComputeUpwindedDensity<Dim, NumNodes>( velocity, upwindVelocity, rCurrentProcessInfo);

    // ====== LOCAL ======
    const BoundedVector<double, NumNodes> DNV_local = prod(current_data.DN_DX, velocity);
    const BoundedMatrix<double, NumNodes, NumNodes> linear_term = current_data.vol * density * prod(current_data.DN_DX, trans(current_data.DN_DX));

    const BoundedVector<double, NumNodes> density_derivative_local = densityDerivativeWRTVelocity * DNV_local;

    rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);
    rLeftHandSideMatrix.clear();

    // local-local
    for (int i = 0; i < NumNodes; ++i)
        for (int j = 0; j < NumNodes; ++j)
            rLeftHandSideMatrix(i, j) = 2.0 * current_data.vol * DNV_local[i] * density_derivative_local[j] + linear_term(i, j);

    // ====== UPWIND ======
    const GeometryType& r_upwind_geom = pGetUpwindElement()->GetGeometry();
    ElementalData upwind_data{r_upwind_geom};

    const BoundedVector<double, NumNodes> DNV_upwind = prod(upwind_data.DN_DX, upwindVelocity);
    const BoundedVector<double, NumNodes> density_derivative_upwind =
        densityDerivativeWRTUpwindVelocity * DNV_upwind;

    // local × upwind
    BoundedMatrix<double, NumNodes, NumNodes> upwind_coupling_matrix = ZeroMatrix(NumNodes, NumNodes);

    for (int i = 0; i < NumNodes; ++i)
        for (int j = 0; j < NumNodes; ++j)
            upwind_coupling_matrix(i, j) = 2.0 * current_data.vol * DNV_local[i] * density_derivative_upwind[j];

    this->SetValue(UPWIND_CONTRIBUTION, upwind_coupling_matrix);
}

template <int Dim, int NumNodes>
BoundedVector<double, NumNodes + 1> CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::AssembleDensityDerivativeAndShapeFunctions(
    const double densityDerivativeWRTVelocitySquared,
    const double densityDerivativeWRTUpwindVelocitySquared,
    const array_1d<double, Dim>& velocity,
    const array_1d<double, Dim>& upwindVelocity,
    const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& r_geom = this->GetGeometry();
    const GeometryType& r_upwind_geom = pGetUpwindElement()->GetGeometry();

    const array_1d<size_t, NumNodes> upwind_node_key = GetAssemblyKey(r_geom, r_upwind_geom, rCurrentProcessInfo);

    ElementalData currentElementdata{r_geom};
    ElementalData upwindElementdata{r_upwind_geom};

    const BoundedVector<double, NumNodes> current_DNV = densityDerivativeWRTVelocitySquared * prod(currentElementdata.DN_DX, velocity);
    const BoundedVector<double, NumNodes> upwind_DNV = densityDerivativeWRTUpwindVelocitySquared * prod(upwindElementdata.DN_DX, upwindVelocity);

    BoundedVector<double, NumNodes + 1> assembly_DNV;
    assembly_DNV.clear();

    for (int i = 0; i < NumNodes; i++)
    {
        assembly_DNV[i] += current_DNV[i];
        assembly_DNV[upwind_node_key[i]] += upwind_DNV[i];
    }

    return assembly_DNV;
}

template <int Dim, int NumNodes>
array_1d<size_t, NumNodes> CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::GetAssemblyKey(
    const GeometryType& rGeom,
    const GeometryType& rUpwindGeom,
    const ProcessInfo& rCurrentProcessInfo)
{
    array_1d<size_t, NumNodes> upwind_node_key;
    EquationIdVectorType upwind_element_ids, current_element_ids;

    upwind_node_key.clear();

    pGetUpwindElement()->EquationIdVector(upwind_element_ids, rCurrentProcessInfo);
    this->EquationIdVector(current_element_ids, rCurrentProcessInfo);

    for (int i = 0; i < NumNodes; i++) {
        auto current_id = std::find(current_element_ids.begin(), current_element_ids.end(),
                upwind_element_ids[i]);

        upwind_node_key[i] = std::distance(current_element_ids.begin(), current_id);
    }

    return upwind_node_key;
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateRightHandSideNormalElement(
    VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != NumNodes)
        rRightHandSideVector.resize(NumNodes, false);
    rRightHandSideVector.clear();

    ElementalData data{GetGeometry()};

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    const CompressiblePerturbationPotentialFlowElement& r_this = *this;

    array_1d<double, Dim> velocity = PotentialFlowUtilities::ComputePerturbedVelocity<Dim,NumNodes>(r_this, rCurrentProcessInfo);

    BoundedVector<double, NumNodes> rhs = ZeroVector(NumNodes);
    CalculateRightHandSideContribution(rhs, rCurrentProcessInfo, velocity, data);

    noalias(rRightHandSideVector) = rhs;
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateLeftHandSideWakeElement(
    MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    // Note that the lhs and rhs have double the size
    if (rLeftHandSideMatrix.size1() != 2 * NumNodes ||
        rLeftHandSideMatrix.size2() != 2 * NumNodes)
        rLeftHandSideMatrix.resize(2 * NumNodes, 2 * NumNodes, false);
    rLeftHandSideMatrix.clear();

    ElementalData data{GetGeometry()};

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);
    GetWakeDistances(data.distances);

    // Compute upper and lower velocities
    const array_1d<double, Dim>& upper_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<Dim,NumNodes>(*this, rCurrentProcessInfo);
    const array_1d<double, Dim>& lower_velocity = PotentialFlowUtilities::ComputePerturbedVelocityLowerElement<Dim,NumNodes>(*this, rCurrentProcessInfo);

    BoundedMatrix<double, NumNodes, NumNodes> upper_lhs_total = ZeroMatrix(NumNodes,NumNodes);
    BoundedMatrix<double, NumNodes, NumNodes> lower_lhs_total = ZeroMatrix(NumNodes,NumNodes);

    CalculateLeftHandSideContribution(upper_lhs_total, rCurrentProcessInfo, upper_velocity, data);
    CalculateLeftHandSideContribution(lower_lhs_total, rCurrentProcessInfo, lower_velocity, data);

    // Compute lhs wake condition
    const BoundedMatrix<double, NumNodes, NumNodes> lhs_wake_condition =
        CalculateLeftHandSideWakeConditions(data, rCurrentProcessInfo);

    if (this->Is(STRUCTURE)){
        Matrix lhs_positive = ZeroMatrix(NumNodes, NumNodes);
        Matrix lhs_negative = ZeroMatrix(NumNodes, NumNodes);

        CalculateLeftHandSideSubdividedElement(lhs_positive, lhs_negative, rCurrentProcessInfo);
        AssignLeftHandSideSubdividedElement(rLeftHandSideMatrix, lhs_positive,
                                            lhs_negative, upper_lhs_total, lower_lhs_total, lhs_wake_condition, data);
    }
    else{
        AssignLeftHandSideWakeElement(rLeftHandSideMatrix, upper_lhs_total, lower_lhs_total, lhs_wake_condition, data);
    }
}

// In 2D
template <>
BoundedMatrix<double, 3, 3> CompressiblePerturbationPotentialFlowElement<2, 3>::CalculateLeftHandSideWakeConditions(
    const ElementalData& rData, const ProcessInfo& rCurrentProcessInfo)
{
    const double free_stream_density = rCurrentProcessInfo[FREE_STREAM_DENSITY];
    return rData.vol * free_stream_density * prod(rData.DN_DX, trans(rData.DN_DX));
}

// In 3D
template <>
BoundedMatrix<double, 4, 4> CompressiblePerturbationPotentialFlowElement<3, 4>::CalculateLeftHandSideWakeConditions(
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

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateRightHandSideWakeElement(
    VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    // Note that the rhs has double the size
    if (rRightHandSideVector.size() != 2 * NumNodes)
        rRightHandSideVector.resize(2 * NumNodes, false);
    rRightHandSideVector.clear();

    ElementalData data{GetGeometry()};

    // Calculate shape functions
    const auto& r_geometry = this->GetGeometry();
    GeometryUtils::CalculateGeometryData(r_geometry, data.DN_DX, data.N, data.vol);
    GetWakeDistances(data.distances);

    // Compute upper and lower velocities
    const array_1d<double, Dim>& upper_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<Dim,NumNodes>(*this, rCurrentProcessInfo);
    const array_1d<double, Dim>& lower_velocity = PotentialFlowUtilities::ComputePerturbedVelocityLowerElement<Dim,NumNodes>(*this, rCurrentProcessInfo);

    // Compute upper and lower rhs
    BoundedVector<double, NumNodes> upper_rhs = ZeroVector(NumNodes);
    BoundedVector<double, NumNodes> lower_rhs = ZeroVector(NumNodes);
    CalculateRightHandSideContribution(upper_rhs, rCurrentProcessInfo, upper_velocity, data);
    CalculateRightHandSideContribution(lower_rhs, rCurrentProcessInfo, lower_velocity, data);

    const array_1d<double, Dim>& diff_velocity = upper_velocity - lower_velocity;

    // Compute wake condition rhs
    const BoundedVector<double, NumNodes> wake_rhs =
        CalculateRightHandSideWakeConditions(data, rCurrentProcessInfo, diff_velocity);

    if (this->Is(STRUCTURE)){
        double upper_vol = 0.0;
        double lower_vol = 0.0;

        CalculateVolumesSubdividedElement(upper_vol, lower_vol, rCurrentProcessInfo);
        for (unsigned int i = 0; i < NumNodes; ++i){
            if (r_geometry[i].GetValue(TRAILING_EDGE)){
                rRightHandSideVector[i] = upper_rhs(i)*upper_vol/data.vol;
                rRightHandSideVector[i + NumNodes] = lower_rhs(i)*lower_vol/data.vol;
            }
            else{
                AssignRightHandSideWakeNode(rRightHandSideVector, upper_rhs, lower_rhs, wake_rhs, data, i);
            }
        }
    }
    else{
        for (unsigned int i = 0; i < NumNodes; ++i){
            AssignRightHandSideWakeNode(rRightHandSideVector, upper_rhs, lower_rhs, wake_rhs, data, i);
        }
    }
}

// In 2D
template <>
BoundedVector<double, 3> CompressiblePerturbationPotentialFlowElement<2, 3>::CalculateRightHandSideWakeConditions(
    const ElementalData& rData,
    const ProcessInfo& rCurrentProcessInfo,
    const array_1d<double, 2>& rDiff_velocity)
{
    const double free_stream_density = rCurrentProcessInfo[FREE_STREAM_DENSITY];
    return -rData.vol * free_stream_density * prod(rData.DN_DX, rDiff_velocity);
}

// In 3D
template <>
BoundedVector<double, 4> CompressiblePerturbationPotentialFlowElement<3, 4>::CalculateRightHandSideWakeConditions(
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

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateLeftHandSideContribution(
    BoundedMatrix<double, NumNodes, NumNodes>& rLhs_total,
    const ProcessInfo& rCurrentProcessInfo,
    const array_1d<double, Dim>& rVelocity,
    const ElementalData& rData)
{
    // Compute density
    const double local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<Dim, NumNodes>(rVelocity, rCurrentProcessInfo);
    const double density = PotentialFlowUtilities::ComputeDensity<Dim, NumNodes>(local_mach_number_squared, rCurrentProcessInfo);

    // Compute density derivative
    const double DrhoDu2 = PotentialFlowUtilities::ComputeDensityDerivativeWRTVelocitySquared<Dim, NumNodes>(local_mach_number_squared, rCurrentProcessInfo);

    // Compute lhs
    const BoundedVector<double, NumNodes> DNV = prod(rData.DN_DX, rVelocity);

    rLhs_total = rData.vol * density * prod(rData.DN_DX, trans(rData.DN_DX));

    const double local_velocity_squared = inner_prod(rVelocity, rVelocity);

    const double max_velocity_squared = PotentialFlowUtilities::ComputeMaximumVelocitySquared<Dim, NumNodes>(rCurrentProcessInfo);
    if (local_velocity_squared < max_velocity_squared){
        rLhs_total += rData.vol * 2 * DrhoDu2 * outer_prod(DNV, trans(DNV));
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateRightHandSideContribution(
    BoundedVector<double, NumNodes>& rRhs_total,
    const ProcessInfo& rCurrentProcessInfo,
    const array_1d<double, Dim>& rVelocity,
    const ElementalData& rData)
{
    // Compute density
    double density = 0.0;
    const double local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<Dim, NumNodes>(rVelocity, rCurrentProcessInfo);
    const double critical_mach_sq = std::pow(rCurrentProcessInfo[CRITICAL_MACH], 2.0);

    if (local_mach_number_squared < critical_mach_sq) {// subsonic
        density = PotentialFlowUtilities::ComputeDensity<Dim, NumNodes>(local_mach_number_squared, rCurrentProcessInfo);
    }else {// supersonic
        const array_1d<double, Dim> upwind_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<Dim,NumNodes>(*pGetUpwindElement(), rCurrentProcessInfo);
        density = PotentialFlowUtilities::ComputeUpwindedDensity<Dim, NumNodes>(rVelocity, upwind_velocity, rCurrentProcessInfo);
    }

    rRhs_total = - rData.vol * density * prod(rData.DN_DX, rVelocity);
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateLeftHandSideSubdividedElement(
    Matrix& lhs_positive,
    Matrix& lhs_negative,
    const ProcessInfo& rCurrentProcessInfo)
{
    ElementalData data{GetGeometry()};

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

    // Compute upper and lower velocities
    const array_1d<double, Dim>& upper_velocity = PotentialFlowUtilities::ComputePerturbedVelocity<Dim,NumNodes>(*this, rCurrentProcessInfo);
    const array_1d<double, Dim>& lower_velocity = PotentialFlowUtilities::ComputePerturbedVelocityLowerElement<Dim,NumNodes>(*this, rCurrentProcessInfo);

    // Compute upper and lower densities
    const double upper_local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<Dim, NumNodes>(upper_velocity, rCurrentProcessInfo);
    const double upper_density = PotentialFlowUtilities::ComputeDensity<Dim, NumNodes>(upper_local_mach_number_squared, rCurrentProcessInfo);

    const double lower_local_mach_number_squared = PotentialFlowUtilities::ComputeLocalMachNumberSquared<Dim, NumNodes>(lower_velocity, rCurrentProcessInfo);
    const double lower_density = PotentialFlowUtilities::ComputeDensity<Dim, NumNodes>(lower_local_mach_number_squared, rCurrentProcessInfo);

    // Compute upper and lower density derivatives
    const double upper_DrhoDu2 = PotentialFlowUtilities::ComputeDensityDerivativeWRTVelocitySquared<Dim, NumNodes>(upper_local_mach_number_squared, rCurrentProcessInfo);
    const double lower_DrhoDu2 = PotentialFlowUtilities::ComputeDensityDerivativeWRTVelocitySquared<Dim, NumNodes>(lower_local_mach_number_squared, rCurrentProcessInfo);

    // Compute upper and lower lhs
    const BoundedVector<double, NumNodes> upper_DNV = prod(data.DN_DX, upper_velocity);
    const BoundedVector<double, NumNodes> lower_DNV = prod(data.DN_DX, lower_velocity);

    const double upper_local_velocity_squared = inner_prod(upper_velocity, upper_velocity);
    const double lower_local_velocity_squared = inner_prod(lower_velocity, lower_velocity);
    const double max_velocity_squared = PotentialFlowUtilities::ComputeMaximumVelocitySquared<Dim, NumNodes>(rCurrentProcessInfo);

    // Compute the lhs and rhs that would correspond to it being divided
    for (unsigned int i = 0; i < nsubdivisions; ++i)
    {
        if (PartitionsSign[i] > 0)
        {
            noalias(lhs_positive) += Volumes[i] * upper_density * prod(data.DN_DX, trans(data.DN_DX));
            if( upper_local_velocity_squared < max_velocity_squared){
                noalias(lhs_positive) += Volumes[i] * 2 * upper_DrhoDu2 * outer_prod(upper_DNV, trans(upper_DNV));
            }

        }
        else
        {
            noalias(lhs_negative) += Volumes[i] * lower_density * prod(data.DN_DX, trans(data.DN_DX));
            if( lower_local_velocity_squared < max_velocity_squared){
                noalias(lhs_negative) += Volumes[i] * 2 * lower_DrhoDu2 * outer_prod(lower_DNV, trans(lower_DNV));
            }
        }
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::CalculateVolumesSubdividedElement(
    double& rUpper_vol,
    double& rLower_vol,
    const ProcessInfo& rCurrentProcessInfo)
{
    ElementalData data{GetGeometry()};

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

    // Compute the volumes that would correspond to it being divided
    for (unsigned int i = 0; i < nsubdivisions; ++i)
    {
        if (PartitionsSign[i] > 0){
            rUpper_vol += Volumes[i];
        }
        else{
            rLower_vol += Volumes[i];
        }
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::ComputeLHSGaussPointContribution(
    const double weight, Matrix& lhs, const ElementalData& data) const
{
    noalias(lhs) += weight * prod(data.DN_DX, trans(data.DN_DX));
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::AssignLeftHandSideSubdividedElement(
    Matrix& rLeftHandSideMatrix,
    Matrix& lhs_positive,
    Matrix& lhs_negative,
    const BoundedMatrix<double, NumNodes, NumNodes>& rUpper_lhs_total,
    const BoundedMatrix<double, NumNodes, NumNodes>& rLower_lhs_total,
    const BoundedMatrix<double, NumNodes, NumNodes>& rLhs_wake_condition,
    const ElementalData& data) const
{
    const auto& r_geometry = this->GetGeometry();
    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        // The TE node takes the contribution of the subdivided element and
        // we do not apply the wake condition on the TE node
        if (r_geometry[i].GetValue(TRAILING_EDGE)){
            for (unsigned int j = 0; j < NumNodes; ++j){
                rLeftHandSideMatrix(i, j) = lhs_positive(i, j);
                rLeftHandSideMatrix(i + NumNodes, j + NumNodes) = lhs_negative(i, j);
            }
        }
        else{
            AssignLeftHandSideWakeNode(rLeftHandSideMatrix, rUpper_lhs_total, rLower_lhs_total, rLhs_wake_condition, data, i);
        }
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::AssignLeftHandSideWakeElement(
    MatrixType& rLeftHandSideMatrix,
    const BoundedMatrix<double, NumNodes, NumNodes>& rUpper_lhs_total,
    const BoundedMatrix<double, NumNodes, NumNodes>& rLower_lhs_total,
    const BoundedMatrix<double, NumNodes, NumNodes>& rLhs_wake_condition,
    const ElementalData& rData) const
{
    for (unsigned int row = 0; row < NumNodes; ++row){
        AssignLeftHandSideWakeNode(rLeftHandSideMatrix, rUpper_lhs_total, rLower_lhs_total, rLhs_wake_condition, rData, row);
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::AssignLeftHandSideWakeNode(
    MatrixType& rLeftHandSideMatrix,
    const BoundedMatrix<double, NumNodes, NumNodes>& rUpper_lhs_total,
    const BoundedMatrix<double, NumNodes, NumNodes>& rLower_lhs_total,
    const BoundedMatrix<double, NumNodes, NumNodes>& rLhs_wake_condition,
    const ElementalData& rData,
    unsigned int row) const
{
    // Applying wake condition on the AUXILIARY_VELOCITY_POTENTIAL dofs
    if (rData.distances[row] < 0.0){
        for (unsigned int column = 0; column < NumNodes; ++column){
            // Conservation of mass
            rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = rLower_lhs_total(row, column);
            // Wake condition
            rLeftHandSideMatrix(row, column) = rLhs_wake_condition(row, column); // Diagonal
            rLeftHandSideMatrix(row, column + NumNodes) = -rLhs_wake_condition(row, column); // Off diagonal
        }
    }
    else{ // else if (data.distances[row] > 0.0)
        for (unsigned int column = 0; column < NumNodes; ++column){
            // Conservation of mass
            rLeftHandSideMatrix(row, column) = rUpper_lhs_total(row, column);
            // Wake condition
            rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = rLhs_wake_condition(row, column); // Diagonal
            rLeftHandSideMatrix(row + NumNodes, column) = -rLhs_wake_condition(row, column); // Off diagonal
        }
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::AssignRightHandSideWakeNode(
    VectorType& rRightHandSideVector,
    const BoundedVector<double, NumNodes>& rUpper_rhs,
    const BoundedVector<double, NumNodes>& rLower_rhs,
    const BoundedVector<double, NumNodes>& rWake_rhs,
    const ElementalData& rData,
    unsigned int& rRow) const
{
    if (rData.distances[rRow] > 0.0){
        rRightHandSideVector[rRow] = rUpper_rhs(rRow);
        rRightHandSideVector[rRow + NumNodes] = -rWake_rhs(rRow);
    }
    else{
        rRightHandSideVector[rRow] = rWake_rhs(rRow);
        rRightHandSideVector[rRow + NumNodes] = rLower_rhs(rRow);
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::FindUpwindElement(const ProcessInfo& rCurrentProcessInfo)
{
    GeometryType upwind_element_boundary;
    FindUpwindEdge(upwind_element_boundary, rCurrentProcessInfo);
    std::vector<size_t> upwind_element_nodes;
    PotentialFlowUtilities::GetSortedIds<Dim, NumNodes>(upwind_element_nodes, upwind_element_boundary);

    GlobalPointersVector<Element> upwind_element_candidates;
    PotentialFlowUtilities::GetNodeNeighborElementCandidates<Dim, NumNodes>(upwind_element_candidates, upwind_element_boundary);
    SelectUpwindElement(upwind_element_nodes, upwind_element_candidates);

    pGetUpwindElement()->GetValue(ID_DOWNWIND_ELEMENT) = this->Id();
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::FindDownwindElement(const ProcessInfo& rCurrentProcessInfo)
{
    GlobalPointersVector<Element> downwind_element_candidates;
    PotentialFlowUtilities::GetNodeNeighborElementCandidates<Dim, NumNodes>(downwind_element_candidates, this->GetGeometry());

    for (SizeType i = 0; i < downwind_element_candidates.size(); i++)
    {
        if(downwind_element_candidates[i].Id() ==  static_cast<std::size_t>(this->GetValue(ID_DOWNWIND_ELEMENT)))
        {
            mpDownwindElement = downwind_element_candidates(i);
            break;
        }
    }

    // If no downwind element is found the downwind element pointer points to itself
    if (this->CheckDownwindElement())
    {
        this->pSetDownwindElement(this);
    }
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::FindUpwindEdge(GeometryType& rUpwindEdge,
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

template<int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::GetElementGeometryBoundary(GeometriesArrayType& rElementGeometryBoundary)
{
    const CompressiblePerturbationPotentialFlowElement& r_this = *this;

    // current element geometry
    const GeometryType& r_geom = r_this.GetGeometry();

    // get element edges or faces depending on dimension of the problem
    if constexpr (Dim == 2)
    {
        // current element edges
        rElementGeometryBoundary = r_geom.GenerateEdges();
    }
    else if constexpr (Dim == 3)
    {
        // current element faces
        rElementGeometryBoundary = r_geom.GenerateFaces();
    }
}

template <int Dim, int NumNodes>
array_1d<double, 3> CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::GetEdgeNormal(const GeometryType& rEdge)
{
    // get local coordinates of edge center
    array_1d<double, 3> edge_center_coordinates;
    rEdge.PointLocalCoordinates(edge_center_coordinates, rEdge.Center());

    // outward pointing normals of each edge
    return rEdge.Normal(edge_center_coordinates);
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::SelectUpwindElement(
    std::vector<IndexType>& rUpwindElementNodeIds,
    GlobalPointersVector<Element>& rUpwindElementCandidates)
{
    for (SizeType i = 0; i < rUpwindElementCandidates.size(); i++)
    {
        // get sorted node ids of neighbording elements
        std::vector<size_t> neighbor_element_ids;
        PotentialFlowUtilities::GetSortedIds<Dim, NumNodes>(neighbor_element_ids, rUpwindElementCandidates[i].GetGeometry());

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
    if (this->CheckUpwindElement())
    {
        this->pSetUpwindElement(this);
        this->SetFlags(INLET);
    }
}

// serializer

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}

template <int Dim, int NumNodes>
void CompressiblePerturbationPotentialFlowElement<Dim, NumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class CompressiblePerturbationPotentialFlowElement<2, 3>;
template class CompressiblePerturbationPotentialFlowElement<3, 4>;

} // namespace Kratos
