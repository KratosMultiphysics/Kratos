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

#include "transonic_perturbation_potential_flow_element.h"
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
Element::Pointer TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<TransonicPerturbationPotentialFlowElement>(
        NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<TransonicPerturbationPotentialFlowElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::Clone(
    IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<TransonicPerturbationPotentialFlowElement>(
        NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);
    CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    const TransonicPerturbationPotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) // Normal element (non-wake) - eventually an embedded
        CalculateRightHandSideNormalElement(rRightHandSideVector, rCurrentProcessInfo);
    else // Wake element
        CalculateRightHandSideWakeElement(rRightHandSideVector, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    const TransonicPerturbationPotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) // Normal element (non-wake) - eventually an embedded
        CalculateLeftHandSideNormalElement(rLeftHandSideMatrix, rCurrentProcessInfo);
    else // Wake element
        CalculateLeftHandSideWakeElement(rLeftHandSideMatrix, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::EquationIdVector(
    EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    const TransonicPerturbationPotentialFlowElement& r_this = *this;
    const Kratos::Element& r_upstream_element = *pGetUpstreamElement();
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) // Normal element
    {
        if(r_this.Is(INLET)){
            if (rResult.size() != NumNodes)
                rResult.resize(NumNodes, false);

            const int kutta = r_this.GetValue(KUTTA);

            if (kutta == 0)
                PotentialFlowUtilities::GetEquationIdVectorNormalElement<Dim,NumNodes>(*this, rResult);
            else
                KRATOS_ERROR << "Found kutta element at the inle with ID #" << r_this.Id() << std::endl;
        }
        else{
            if (rResult.size() != NumNodes + 1)
                rResult.resize(NumNodes + 1, false);
            for(unsigned int i = 0; i < rResult.size(); i++){
                rResult[i] = -1;
            }
            EquationIdVectorType upstream_element_ids(NumNodes, 0);

            const int kutta = r_this.GetValue(KUTTA);
            if (kutta == 0){
                PotentialFlowUtilities::GetEquationIdVectorNormalElement<Dim,NumNodes>(*this, rResult);
            }
            else{
                PotentialFlowUtilities::GetEquationIdVectorKuttaElement<Dim,NumNodes>(*this, rResult);
            }

            const int upstream_kutta = r_upstream_element.GetValue(KUTTA);
            if (upstream_kutta == 0){
                PotentialFlowUtilities::GetEquationIdVectorNormalElement<Dim,NumNodes>(*pGetUpstreamElement(), upstream_element_ids);
            }
            else{
                PotentialFlowUtilities::GetEquationIdVectorKuttaElement<Dim,NumNodes>(*pGetUpstreamElement(), upstream_element_ids);
            }

            bool upstream_element_id_found = false;
            for(unsigned int i = 0; i < upstream_element_ids.size(); i++){
                if(std::find(rResult.begin(), rResult.end(), upstream_element_ids[i]) == rResult.end()){
                    rResult[NumNodes] = upstream_element_ids[i];
                    upstream_element_id_found = true;
                }
            }
            KRATOS_ERROR_IF(!upstream_element_id_found) << "No upstream element id found for element #" << this->Id() << std::endl;
        }

    }
    else // Wake element
    {
        if (rResult.size() != 2 * NumNodes)
            rResult.resize(2 * NumNodes, false);

        GetEquationIdVectorWakeElement(rResult);
    }
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                                                 ProcessInfo& CurrentProcessInfo)
{
    const TransonicPerturbationPotentialFlowElement& r_this = *this;
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
    if(this->Id()==333){
        KRATOS_WATCH(rElementalDofList)
    }
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    FindUpstreamElementSharingFace(rCurrentProcessInfo);
    //FindUpstreamElementSharingNode(rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    FindUpstreamElementSharingFace(rCurrentProcessInfo);
    //FindUpstreamElementSharingNode(rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    bool active = true;
    if ((this)->IsDefined(ACTIVE))
        active = (this)->Is(ACTIVE);

    const TransonicPerturbationPotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake != 0 && active == true)
    {
        ComputePotentialJump(rCurrentProcessInfo);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <int Dim, int NumNodes>
int TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
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
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::GetValueOnIntegrationPoints(
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
        rValues[0] = ComputeDensity(rCurrentProcessInfo);
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
        const TransonicPerturbationPotentialFlowElement& r_this = *this;
        rValues[0] = r_this.GetValue(WAKE);
    }
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::GetValueOnIntegrationPoints(
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
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::GetValueOnIntegrationPoints(
    const Variable<bool>& rVariable, std::vector<bool>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);
    if (rVariable == INLET_ELEMENT)
        rValues[0] = this->Is(INLET);
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);
    if (rVariable == VELOCITY){
        array_1d<double, 3> v(3, 0.0);
        array_1d<double, Dim> velocity = ComputeVelocity(rCurrentProcessInfo);
        for (unsigned int k = 0; k < Dim; k++)
            v[k] = velocity[k];
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
    else if (rVariable == UPSTREAM_ELEMENT_POSITION_VECTOR)
    {
        rValues[0] = pGetUpstreamElement()->GetGeometry().Center()-this->GetGeometry().Center();
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <int Dim, int NumNodes>
std::string TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "TransonicPerturbationPotentialFlowElement #" << Id();
    return buffer.str();
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "TransonicPerturbationPotentialFlowElement #" << Id();
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::PrintData(std::ostream& rOStream) const
{
    pGetGeometry()->PrintData(rOStream);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::GetWakeDistances(
    array_1d<double, NumNodes>& distances) const
{
    noalias(distances) = GetValue(WAKE_ELEMENTAL_DISTANCES);
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::GetEquationIdVectorNormalElement(
    EquationIdVectorType& rResult) const
{
    for (unsigned int i = 0; i < NumNodes; i++)
        rResult[i] = GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::GetEquationIdVectorKuttaElement(
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
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::GetEquationIdVectorWakeElement(
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
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::GetDofListNormalElement(DofsVectorType& rElementalDofList) const
{
    for (unsigned int i = 0; i < NumNodes; i++)
        rElementalDofList[i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::GetDofListKuttaElement(DofsVectorType& rElementalDofList) const
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
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::GetDofListWakeElement(DofsVectorType& rElementalDofList) const
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
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::CalculateLeftHandSideNormalElement(
    MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    const TransonicPerturbationPotentialFlowElement& r_this = *this;
    if(r_this.Is(INLET)){
        if (rLeftHandSideMatrix.size1() != NumNodes|| rLeftHandSideMatrix.size2() != NumNodes)
            rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);
        rLeftHandSideMatrix.clear();

        ElementalData<NumNodes, Dim> data;

        // Calculate shape functions
        GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

        const double density = ComputeDensity(rCurrentProcessInfo);
        const double DrhoDu2 = ComputeDensityDerivative(density, rCurrentProcessInfo);
        array_1d<double, Dim> velocity = ComputeVelocity(rCurrentProcessInfo);
        const BoundedVector<double, NumNodes> DNV = prod(data.DN_DX, velocity);

        noalias(rLeftHandSideMatrix) += data.vol * density * prod(data.DN_DX, trans(data.DN_DX));
        noalias(rLeftHandSideMatrix) += data.vol * 2 * DrhoDu2 * outer_prod(DNV, trans(DNV));
    }
    else{
        if (rLeftHandSideMatrix.size1() != NumNodes + 1 || rLeftHandSideMatrix.size2() != NumNodes + 1)
            rLeftHandSideMatrix.resize(NumNodes + 1, NumNodes + 1, false);
        rLeftHandSideMatrix.clear();

        ElementalData<NumNodes, Dim> data;

        // Calculate shape functions
        GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

        const Kratos::Element& r_upstream_element = *pGetUpstreamElement();
        const double upwind_density = PotentialFlowUtilities::ComputeUpwindDensity<Dim, NumNodes>(*this, r_upstream_element, rCurrentProcessInfo);
        const double density = PotentialFlowUtilities::ComputePerturbationDensity<Dim, NumNodes>(*this, rCurrentProcessInfo);
        const double Drho_Dv2 = ComputeDensityDerivative(density, rCurrentProcessInfo);
        array_1d<double, Dim> velocity = ComputeVelocity(rCurrentProcessInfo);
        const BoundedVector<double, NumNodes> DNV = prod(data.DN_DX, velocity);

        const double upwind_factor = PotentialFlowUtilities::ComputeUpwindFactor<Dim, NumNodes>(*this, rCurrentProcessInfo);
        const double upstream_upwind_factor = PotentialFlowUtilities::ComputeUpwindFactor<Dim, NumNodes>(r_upstream_element, rCurrentProcessInfo);

        const BoundedMatrix<double, NumNodes, NumNodes> term_matrix_laplacian =
            data.vol * upwind_density * prod(data.DN_DX, trans(data.DN_DX));

        // Subsonic flow (local_mach_number < mach_number_limit)
        if(upwind_factor < 0.0){
            const BoundedMatrix<double, NumNodes, NumNodes> term_matrix_nonlinear =
                data.vol * 2 * Drho_Dv2 * outer_prod(DNV, trans(DNV));

            for(unsigned int i = 0; i < NumNodes; i++){
                for(unsigned int j = 0; j < NumNodes; j++){
                    rLeftHandSideMatrix(i,j)  = term_matrix_laplacian(i,j);
                    rLeftHandSideMatrix(i,j) += term_matrix_nonlinear(i,j);
                }
            }
        }
        // Supersonic flow and accelerating (local_mach_number > upstream_mach_number)
        else if( upwind_factor >= upstream_upwind_factor){
            BoundedVector<double, NumNodes> Drho_DPhi_current =
                PotentialFlowUtilities::ComputeDrhoDphiSupersonicAccelerating<Dim, NumNodes>(
                    *this, r_upstream_element, rCurrentProcessInfo);

            // Vector containing this and upstream element contributions
            BoundedVector<double, NumNodes + 1> Drho_DPhi = ZeroVector(NumNodes + 1);
            // Assembling contributions from this element
            for(unsigned int i = 0; i < NumNodes; i++) {
                Drho_DPhi(i) = Drho_DPhi_current(i);
            }

            // Computing upstream contributions
            BoundedVector<double, NumNodes> Drho_DPhi_upstream =
                PotentialFlowUtilities::ComputeDrhoDphiUpSupersonicAccelerating<Dim, NumNodes>(
                    *this, r_upstream_element, rCurrentProcessInfo);

            EquationIdVectorType equation_id_vector(NumNodes, 0);
            EquationIdVectorType upstream_equation_id_vector(NumNodes, 0);

            const int kutta = r_this.GetValue(KUTTA);
            if (kutta == 0){
                PotentialFlowUtilities::GetEquationIdVectorNormalElement<Dim,NumNodes>(*this, equation_id_vector);
            }
            else{
                PotentialFlowUtilities::GetEquationIdVectorKuttaElement<Dim,NumNodes>(*this, equation_id_vector);
            }

            const int upstream_kutta = r_upstream_element.GetValue(KUTTA);
            if (upstream_kutta == 0){
                PotentialFlowUtilities::GetEquationIdVectorNormalElement<Dim,NumNodes>(r_upstream_element, upstream_equation_id_vector);
            }
            else{
                PotentialFlowUtilities::GetEquationIdVectorKuttaElement<Dim,NumNodes>(r_upstream_element, upstream_equation_id_vector);
            }

            // Loop over upstream_equation_id_vector
            for(unsigned int i = 0; i < upstream_equation_id_vector.size(); i++){
                bool position_found = false;
                // Loop over equation_id_vector
                for(unsigned int j = 0; j < equation_id_vector.size(); j++){
                    if (abs(upstream_equation_id_vector[i] - equation_id_vector[j]) < 1e-3) {
                        Drho_DPhi(j) += Drho_DPhi_upstream(i);
                        position_found = true;
                        break;
                    }
                }
                if(!position_found){
                    Drho_DPhi(NumNodes) = Drho_DPhi_upstream(i);
                }
            }

            const BoundedMatrix<double, NumNodes, NumNodes + 1> term_matrix_nonlinear =
                data.vol * outer_prod(DNV, trans(Drho_DPhi));

            for(unsigned int i = 0; i < NumNodes; i++){
                for(unsigned int j = 0; j < NumNodes; j++){
                    rLeftHandSideMatrix(i,j)  = term_matrix_laplacian(i,j);
                }
            }

            for(unsigned int i = 0; i < NumNodes; i++){
                for(unsigned int j = 0; j < NumNodes + 1; j++){
                    rLeftHandSideMatrix(i,j) += term_matrix_nonlinear(i,j);
                }
            }

        }
        // Supersonic flow and decelerating (local_mach_number < upstream_mach_number)
        else{
        }


    }
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::CalculateRightHandSideNormalElement(
    VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    const TransonicPerturbationPotentialFlowElement& r_this = *this;
    if(r_this.Is(INLET)){
        if (rRightHandSideVector.size() != NumNodes)
            rRightHandSideVector.resize(NumNodes, false);
        rRightHandSideVector.clear();

        ElementalData<NumNodes, Dim> data;

        // Calculate shape functions
        GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

        const double density = ComputeDensity(rCurrentProcessInfo);
        array_1d<double, Dim> velocity = ComputeVelocity(rCurrentProcessInfo);

        noalias(rRightHandSideVector) = - data.vol * density * prod(data.DN_DX, velocity);
    }
    else{
        if (rRightHandSideVector.size() != NumNodes + 1)
            rRightHandSideVector.resize(NumNodes + 1, false);
        rRightHandSideVector.clear();

        ElementalData<NumNodes, Dim> data;

        // Calculate shape functions
        GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

        const double upwind_density = PotentialFlowUtilities::ComputeUpwindDensity<Dim, NumNodes>(*this, *pGetUpstreamElement(), rCurrentProcessInfo);
        array_1d<double, Dim> velocity = ComputeVelocity(rCurrentProcessInfo);

        const BoundedVector<double, NumNodes> rhs = - data.vol * upwind_density * prod(data.DN_DX, velocity);

        for(unsigned int i = 0; i < NumNodes; i++){
            rRightHandSideVector(i) = rhs(i);
        }
    }

}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::CalculateLeftHandSideWakeElement(
    MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    // Note that the lhs and rhs have double the size
    if (rLeftHandSideMatrix.size1() != 2 * NumNodes ||
        rLeftHandSideMatrix.size2() != 2 * NumNodes)
        rLeftHandSideMatrix.resize(2 * NumNodes, 2 * NumNodes, false);
    rLeftHandSideMatrix.clear();

    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);
    GetWakeDistances(data.distances);

    const double density = ComputeDensity(rCurrentProcessInfo);
    const double DrhoDu2 = ComputeDensityDerivative(density, rCurrentProcessInfo);
    array_1d<double, Dim> velocity = ComputeVelocity(rCurrentProcessInfo);
    const BoundedVector<double, NumNodes> DNV = prod(data.DN_DX, velocity);

    const BoundedMatrix<double, NumNodes, NumNodes> lhs_total =
        data.vol * density * prod(data.DN_DX, trans(data.DN_DX)) +
        data.vol * 2 * DrhoDu2 * outer_prod(DNV, trans(DNV));

    if (this->Is(STRUCTURE)){
        Matrix lhs_positive = ZeroMatrix(NumNodes, NumNodes);
        Matrix lhs_negative = ZeroMatrix(NumNodes, NumNodes);

        CalculateLeftHandSideSubdividedElement(lhs_positive, lhs_negative, rCurrentProcessInfo);
        AssignLeftHandSideSubdividedElement(rLeftHandSideMatrix, lhs_positive,
                                            lhs_negative, lhs_total, data);
    }
    else{
        AssignLeftHandSideWakeElement(rLeftHandSideMatrix, lhs_total, data);
    }
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::CalculateRightHandSideWakeElement(
    VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    // Note that the rhs has double the size
    if (rRightHandSideVector.size() != 2 * NumNodes)
        rRightHandSideVector.resize(2 * NumNodes, false);
    rRightHandSideVector.clear();

    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    const auto& r_geometry = this->GetGeometry();
    GeometryUtils::CalculateGeometryData(r_geometry, data.DN_DX, data.N, data.vol);
    GetWakeDistances(data.distances);

    const double density = ComputeDensity(rCurrentProcessInfo);

    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    array_1d<double, Dim> upper_velocity = PotentialFlowUtilities::ComputeVelocityUpperWakeElement<Dim,NumNodes>(*this);
    array_1d<double, Dim> lower_velocity = PotentialFlowUtilities::ComputeVelocityLowerWakeElement<Dim,NumNodes>(*this);

    for (unsigned int i = 0; i < Dim; i++){
        upper_velocity[i] += free_stream_velocity[i];
        lower_velocity[i] += free_stream_velocity[i];
    }
    const array_1d<double, Dim> diff_velocity = upper_velocity - lower_velocity;

    const BoundedVector<double, NumNodes> upper_rhs = - data.vol * density * prod(data.DN_DX, upper_velocity);
    const BoundedVector<double, NumNodes> lower_rhs = - data.vol * density * prod(data.DN_DX, lower_velocity);
    const BoundedVector<double, NumNodes> wake_rhs = - data.vol * density * prod(data.DN_DX, diff_velocity);

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

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::CalculateLeftHandSideSubdividedElement(
    Matrix& lhs_positive,
    Matrix& lhs_negative,
    const ProcessInfo& rCurrentProcessInfo)
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

    const double density = ComputeDensity(rCurrentProcessInfo);
    const double DrhoDu2 = ComputeDensityDerivative(density, rCurrentProcessInfo);
    array_1d<double, Dim> velocity = ComputeVelocity(rCurrentProcessInfo);

    const BoundedVector<double, NumNodes> DNV = prod(data.DN_DX, velocity);

    // Compute the lhs and rhs that would correspond to it being divided
    for (unsigned int i = 0; i < nsubdivisions; ++i)
    {
        if (PartitionsSign[i] > 0)
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

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::CalculateVolumesSubdividedElement(
    double& rUpper_vol,
    double& rLower_vol,
    const ProcessInfo& rCurrentProcessInfo)
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
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::ComputeLHSGaussPointContribution(
    const double weight, Matrix& lhs, const ElementalData<NumNodes, Dim>& data) const
{
    noalias(lhs) += weight * prod(data.DN_DX, trans(data.DN_DX));
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::AssignLeftHandSideSubdividedElement(
    Matrix& rLeftHandSideMatrix,
    Matrix& lhs_positive,
    Matrix& lhs_negative,
    const BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
    const ElementalData<NumNodes, Dim>& data) const
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
            AssignLeftHandSideWakeNode(rLeftHandSideMatrix, lhs_total, data, i);
        }
    }
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::AssignLeftHandSideWakeElement(
    MatrixType& rLeftHandSideMatrix,
    const BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
    const ElementalData<NumNodes, Dim>& data) const
{
    for (unsigned int row = 0; row < NumNodes; ++row)
        AssignLeftHandSideWakeNode(rLeftHandSideMatrix, lhs_total, data, row);
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::AssignLeftHandSideWakeNode(
    MatrixType& rLeftHandSideMatrix,
    const BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
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
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::AssignRightHandSideWakeNode(
    VectorType& rRightHandSideVector,
    const BoundedVector<double, NumNodes>& rUpper_rhs,
    const BoundedVector<double, NumNodes>& rLower_rhs,
    const BoundedVector<double, NumNodes>& rWake_rhs,
    const ElementalData<NumNodes, Dim>& rData,
    unsigned int& rRow) const
{
    if (rData.distances[rRow] > 0.0){
        rRightHandSideVector[rRow] = rUpper_rhs(rRow);
        rRightHandSideVector[rRow + NumNodes] = rWake_rhs(rRow);
    }
    else{
        rRightHandSideVector[rRow] = rWake_rhs(rRow);
        rRightHandSideVector[rRow + NumNodes] = rLower_rhs(rRow);
    }
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::ComputePotentialJump(const ProcessInfo& rCurrentProcessInfo)
{
    const array_1d<double, 3>& vinfinity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    const double vinfinity_norm = sqrt(inner_prod(vinfinity, vinfinity));

    array_1d<double, NumNodes> distances;
    GetWakeDistances(distances);

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        double aux_potential =
            GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL);
        double potential = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
        double potential_jump = aux_potential - potential;

        if (distances[i] > 0)
        {
            GetGeometry()[i].SetValue(POTENTIAL_JUMP,
                                      -2.0 / vinfinity_norm * (potential_jump));
        }
        else
        {
            GetGeometry()[i].SetValue(POTENTIAL_JUMP, 2.0 / vinfinity_norm * (potential_jump));
        }
    }
}

template <int Dim, int NumNodes>
double TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::ComputeDensity(const ProcessInfo& rCurrentProcessInfo) const
{
    // Reading free stream conditions
    const double rho_inf = rCurrentProcessInfo[FREE_STREAM_DENSITY];
    const double M_inf = rCurrentProcessInfo[FREE_STREAM_MACH];
    const double heat_capacity_ratio = rCurrentProcessInfo[HEAT_CAPACITY_RATIO];
    const double mach_number_limit = rCurrentProcessInfo[MACH_LIMIT];

    // Computing local mach number
    double local_mach_number = PotentialFlowUtilities::ComputePerturbationLocalMachNumber<Dim, NumNodes>(*this, rCurrentProcessInfo);

    if (local_mach_number > mach_number_limit)
    { // Clamping the mach number to mach_number_limit
        KRATOS_WARNING("ComputeDensity") << "Clamping the local mach number to " << mach_number_limit << std::endl;
        local_mach_number = mach_number_limit;
    }

    // Computing squares
    const double M_inf_2 = M_inf * M_inf;
    const double M_2 = local_mach_number * local_mach_number;

    // Computing density according to Equation 8.9 of Drela, M. (2014) Flight Vehicle
    // Aerodynamics, The MIT Press, London
    const double numerator = 1 + (heat_capacity_ratio - 1) * M_inf_2 / 2;
    const double denominator = 1 + (heat_capacity_ratio - 1) * M_2 / 2;
    const double base = numerator / denominator;

    if (base > 0.0)
    {
        return rho_inf * pow(base, 1 / (heat_capacity_ratio - 1));
    }
    else
    {
        KRATOS_WARNING("ComputeDensity") << "Using density correction" << std::endl;
        return rho_inf * 0.00001;
    }
}

template <int Dim, int NumNodes>
double TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::ComputeDensityDerivative(
    const double rho, const ProcessInfo& rCurrentProcessInfo) const
{
    // Reading free stream conditions
    const double rho_inf = rCurrentProcessInfo[FREE_STREAM_DENSITY];
    const double heat_capacity_ratio = rCurrentProcessInfo[HEAT_CAPACITY_RATIO];
    const double a_inf = rCurrentProcessInfo[SOUND_VELOCITY];

    return -pow(rho_inf, heat_capacity_ratio - 1) *
           pow(rho, 2 - heat_capacity_ratio) / (2 * a_inf * a_inf);
}

template <int Dim, int NumNodes>
array_1d<double, Dim> TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::ComputeVelocity(
    const ProcessInfo& rCurrentProcessInfo) const
{

    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    array_1d<double, Dim> velocity = PotentialFlowUtilities::ComputeVelocity<Dim,NumNodes>(*this);
    for (unsigned int i = 0; i < Dim; i++){
        velocity[i] += free_stream_velocity[i];
    }

    return velocity;
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::FindUpstreamElementSharingNode(const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& rGeom = this->GetGeometry();
    GlobalPointersVector<Element> ElementsSharingNode;
    GetElementsSharingNode(ElementsSharingNode, rGeom);
    FindUpstreamElement(rCurrentProcessInfo, ElementsSharingNode, rGeom);

    if(mpUpstreamElement.get() == nullptr){
        mpUpstreamElement = this;
        this->Set(INLET);
    }
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::FindUpstreamElementSharingFace(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const GeometryType& rGeom = this->GetGeometry();
    GlobalPointersVector<Element> ElementsSharingFace;
    GetElementsSharingFace(ElementsSharingFace, rGeom);
    FindUpstreamElement(rCurrentProcessInfo, ElementsSharingFace, rGeom);

    // If no upstream element is found, the element should be INLET and the
    // upstream element pointer should point to itself.
    if (mpUpstreamElement.get() == nullptr) {
        mpUpstreamElement = this;
        this->Set(INLET);
    }
    KRATOS_CATCH("")
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::FindUpstreamElement(
    const ProcessInfo& rCurrentProcessInfo,
    GlobalPointersVector<Element>& rElementCandidates,
    const GeometryType& rGeom)
{
    KRATOS_TRY
    // Find the upstream element whose center has the minimum cross flow distance
    const array_1d<double, Dim> velocity = ComputeVelocity(rCurrentProcessInfo);
    const double velocity_norm = sqrt(inner_prod(velocity, velocity));
    array_1d<double, 3> unit_velocity_vector = ZeroVector(3);
    for (unsigned int i = 0; i < Dim; i++){
        unit_velocity_vector[i] = velocity[i] / velocity_norm;
    }
    double minimum_cross_flow_distance = std::numeric_limits<double>::max();

    for (SizeType i = 0; i < rElementCandidates.size(); i++){
        const auto distance_vector = rElementCandidates(i)->GetGeometry().Center() - rGeom.Center();
        const double projection_length = inner_prod(distance_vector, unit_velocity_vector);
        // If the projection is negative the element is upstream
        if(projection_length < 0.0){
            const auto projected_vector = unit_velocity_vector * projection_length;
            const auto cross_flow_vector = distance_vector - projected_vector;
            const double cross_flow_distance = sqrt(inner_prod(cross_flow_vector,cross_flow_vector));
            if(cross_flow_distance < minimum_cross_flow_distance){
                minimum_cross_flow_distance = cross_flow_distance;
                mpUpstreamElement = rElementCandidates(i);
            }
        }
    }
    KRATOS_CATCH("")
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::GetElementsSharingFace(
    GlobalPointersVector<Element>& rElementsSharingFace, const GeometryType& rGeom) const
{
    KRATOS_TRY
    GlobalPointersVector<Element> ElementsSharingNode;
    GetElementsSharingNode(ElementsSharingNode, rGeom);

    std::vector<IndexType> ElementNodesIds, CandidateElementNodesIds;
    GetSortedIds(ElementNodesIds, rGeom);

    for (SizeType i = 0; i < ElementsSharingNode.size(); i++){
        GeometryType& rElemGeom = ElementsSharingNode[i].GetGeometry();
        GetSortedIds(CandidateElementNodesIds, rElemGeom);
        std::vector<IndexType> IntersectionIds;
        std::set_intersection(ElementNodesIds.begin(), ElementNodesIds.end(),
                              CandidateElementNodesIds.begin(),
                              CandidateElementNodesIds.end(),
                              std::back_inserter(IntersectionIds));
        if (IntersectionIds.size() > 1) {
            rElementsSharingFace.push_back(ElementsSharingNode(i));
        }
    }
    rElementsSharingFace.Unique();
    KRATOS_CATCH("")
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::GetElementsSharingNode(
    GlobalPointersVector<Element>& rElementsSharingNode, const GeometryType& rGeom) const
{
    for (SizeType i = 0; i < NumNodes; i++){
        const GlobalPointersVector<Element>& rNodeElementCandidates =
            rGeom[i].GetValue(NEIGHBOUR_ELEMENTS);

        for (SizeType j = 0; j < rNodeElementCandidates.size(); j++){
            rElementsSharingNode.push_back(rNodeElementCandidates(j));
        }
    }
    rElementsSharingNode.Unique();
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::GetSortedIds(std::vector<IndexType>& Ids,
                                                                        const GeometryType& rGeom) const
{
    Ids.resize(rGeom.PointsNumber());
    for (SizeType i = 0; i < Ids.size(); i++)
        Ids[i] = rGeom[i].Id();
    std::sort(Ids.begin(), Ids.end());
}

template <int Dim, int NumNodes>
inline GlobalPointer<Element> TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::pGetUpstreamElement() const
{
    KRATOS_TRY
    KRATOS_ERROR_IF(mpUpstreamElement.get() == nullptr)
        << "No upstream element found for element #" << this->Id() << std::endl;
    return mpUpstreamElement;
    KRATOS_CATCH("");
}

// serializer

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}

template <int Dim, int NumNodes>
void TransonicPerturbationPotentialFlowElement<Dim, NumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class TransonicPerturbationPotentialFlowElement<2, 3>;
template class TransonicPerturbationPotentialFlowElement<3, 4>;

} // namespace Kratos
