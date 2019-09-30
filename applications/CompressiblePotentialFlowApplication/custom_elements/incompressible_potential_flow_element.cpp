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

#include "incompressible_potential_flow_element.h"
#include "compressible_potential_flow_application_variables.h"
#include "includes/cfd_variables.h"
#include "custom_utilities/potential_flow_utilities.h"

namespace Kratos
{
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <int Dim, int NumNodes>
Element::Pointer IncompressiblePotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<IncompressiblePotentialFlowElement>(
        NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer IncompressiblePotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<IncompressiblePotentialFlowElement>(
        NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer IncompressiblePotentialFlowElement<Dim, NumNodes>::Clone(
    IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<IncompressiblePotentialFlowElement>(
        NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    const IncompressiblePotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) // Normal element (non-wake) - eventually an embedded
        CalculateLocalSystemNormalElement(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
    else // Wake element
        CalculateLocalSystemWakeElement(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
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
void IncompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    // TODO: improve speed
    VectorType tmp;
    CalculateLocalSystem(rLeftHandSideMatrix, tmp, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::EquationIdVector(
    EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    const IncompressiblePotentialFlowElement& r_this = *this;
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
        if(this->GetValue(WING_TIP)){
            if (rResult.size() != 2 * NumNodes + 2)
                rResult.resize(2 * NumNodes + 2, false);
        }
        else if(this->GetValue(ZERO_VELOCITY_CONDITION)){
            if (rResult.size() != 2 * NumNodes + 1)
                rResult.resize(2 * NumNodes + 1, false);
        }
        else{
            if (rResult.size() != 2 * NumNodes)
                rResult.resize(2 * NumNodes, false);
        }
        // else if(this->GetValue(ZERO_VELOCITY_CONDITION)){
        //     if (rResult.size() != 2 * NumNodes + 1)
        //         rResult.resize(2 * NumNodes + 1, false);
        // }
        // if (rResult.size() != 2 * NumNodes)
        //     rResult.resize(2 * NumNodes, false);

        GetEquationIdVectorWakeElement(rResult);
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                                                   ProcessInfo& CurrentProcessInfo)
{
    const IncompressiblePotentialFlowElement& r_this = *this;
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
        if(this->GetValue(WING_TIP)){
            if (rElementalDofList.size() != 2 * NumNodes + 2)
                rElementalDofList.resize(2 * NumNodes + 2);
        }
        else if(this->GetValue(ZERO_VELOCITY_CONDITION)){
            if (rElementalDofList.size() != 2 * NumNodes + 1)
                rElementalDofList.resize(2 * NumNodes + 1);
        }
        else{
            if (rElementalDofList.size() != 2 * NumNodes)
                rElementalDofList.resize(2 * NumNodes);
        }
        // else if(this->GetValue(ZERO_VELOCITY_CONDITION)){
        //     if (rElementalDofList.size() != 2 * NumNodes + 1)
        //         rElementalDofList.resize(2 * NumNodes + 1);
        // }
        // if (rElementalDofList.size() != 2 * NumNodes)
        //         rElementalDofList.resize(2 * NumNodes);
        GetDofListWakeElement(rElementalDofList);
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    bool active = true;
    if ((this)->IsDefined(ACTIVE))
        active = (this)->Is(ACTIVE);

    const IncompressiblePotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake != 0 && active == true)
    {
        ComputePotentialJump(rCurrentProcessInfo);

        // for (unsigned int i = 0; i < NumNodes; i++){
        //     if(GetGeometry()[i].GetValue(TRAILING_EDGE)){
        //         const int te_element_counter = GetGeometry()[i].GetValue(TE_ELEMENT_COUNTER);
        //         const int number_of_neighbour_elements = GetGeometry()[i].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
        //         KRATOS_WARNING_IF("element", te_element_counter != number_of_neighbour_elements)
        //             << " NODE WITH ID " << GetGeometry()[i].Id()
        //             << " BELONGING TO ELEMENT " << this-> Id()
        //             << " TE ELEMENTS COUNTERS DO NOT MATCH. \n te_element_counter = " << te_element_counter
        //             << " number_of_neighbour_elements = " << number_of_neighbour_elements << std::endl;
        //     }
        // }
    }
    ComputeElementInternalEnergy();

    // if (r_this.GetValue(WING_TIP)){
    //     ElementalData<NumNodes, Dim> data;

    //     // Calculate shape functions
    //     GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    //     const double free_stream_density = rCurrentProcessInfo[FREE_STREAM_DENSITY];
    //     const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    //     array_1d<double, Dim> v = ZeroVector(Dim);
    //     for (unsigned int i = 0; i < Dim; i++){
    //         v[i] = free_stream_velocity[i];
    //     }

    //     BoundedMatrix<double, NumNodes, NumNodes> lhs_total = ZeroMatrix(NumNodes, NumNodes);

    //     ComputeLHSGaussPointContribution(data.vol*free_stream_density, lhs_total, data);
    //     BoundedVector<double, TNumNodes> right_hand_side_2 = ZeroVector(NumNodes);
    //     const auto& r_distances = PotentialFlowUtilities::GetWakeDistances<Dim,NumNodes>(*this);

    //     data.potentials = PotentialFlowUtilities::GetPotentialOnLowerWakeElement<Dim,NumNodes>(*this, r_distances);
    //     right_hand_side_2 = prod(lhs_total, data.potentials);

    //     BoundedVector<double, TNumNodes> right_hand_side = ZeroVector(NumNodes);
    //     for (unsigned int i = 0; i< NumNodes; i++){
    //         if(GetGeometry()[i].GetValue(WING_TIP)){
    //             for (unsigned int j = 0; j < Dim; j++){
    //                 right_hand_side[i] += data.DN_DX(i,j)*v[j];
    //             }
    //         }
    //     }
    //     auto velocity = PotentialFlowUtilities::ComputeVelocityLowerWakeElement<Dim,NumNodes>(*this);

    //     // noalias(rLeftHandSideMatrix) =
    //     //     data.vol * free_stream_density * prod(data.DN_DX, trans(data.DN_DX));
    //     KRATOS_WATCH(this->Id())
    //     KRATOS_WATCH(v)
    //     KRATOS_WATCH(data.DN_DX)
    //     KRATOS_WATCH(right_hand_side)
    //     KRATOS_WATCH(right_hand_side_2)
    //     KRATOS_WATCH(velocity)
    //     KRATOS_WATCH(lhs_total)
    // }
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

    if (rVariable == PRESSURE_COEFFICIENT)
    {
        rValues[0] = PotentialFlowUtilities::ComputeIncompressiblePressureCoefficient<Dim,NumNodes>(*this,rCurrentProcessInfo);
    }
    else if (rVariable == DENSITY)
    {
        rValues[0] = rCurrentProcessInfo[FREE_STREAM_DENSITY];
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
    else if (rVariable == TRAILING_EDGE_ELEMENT)
        rValues[0] = this->GetValue(TRAILING_EDGE_ELEMENT);
    else if (rVariable == DECOUPLED_TRAILING_EDGE_ELEMENT)
        rValues[0] = this->GetValue(DECOUPLED_TRAILING_EDGE_ELEMENT);
    else if (rVariable == WING_TIP)
        rValues[0] = this->GetValue(WING_TIP);
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
        array_1d<double, Dim> vaux = PotentialFlowUtilities::ComputeVelocity<Dim,NumNodes>(*this);
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

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::GetWakeDistances(array_1d<double, NumNodes>& distances) const
{
    noalias(distances) = GetValue(WAKE_ELEMENTAL_DISTANCES);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::GetEquationIdVectorNormalElement(EquationIdVectorType& rResult) const
{
    for (unsigned int i = 0; i < NumNodes; i++)
        rResult[i] = GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::GetEquationIdVectorKuttaElement(EquationIdVectorType& rResult) const
{
    // Kutta elements have only negative part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (!GetGeometry()[i].GetValue(TRAILING_EDGE))
            rResult[i] = GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
        else{
            rResult[i] = GetGeometry()[i].GetDof(AUXILIARY_VELOCITY_POTENTIAL).EquationId();
            //rResult[i] = GetGeometry()[i].GetDof(PSI).EquationId();
        }
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::GetEquationIdVectorWakeElement(EquationIdVectorType& rResult)
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

    unsigned int te_counter = 0;
    // Negative part - sign is opposite to the previous case
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] < 0.0)
            rResult[NumNodes + i] =
                GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
        else{
            //rResult[NumNodes + i] = GetGeometry()[i].GetDof(AUXILIARY_VELOCITY_POTENTIAL).EquationId();
            if(GetGeometry()[i].GetValue(TRAILING_EDGE)){
                rResult[NumNodes + i] = GetGeometry()[i].GetDof(AUXILIARY_VELOCITY_POTENTIAL).EquationId();
                //rResult[2*NumNodes + te_counter] = GetGeometry()[i].GetDof(PSI).EquationId();

                if(this->GetValue(WING_TIP)){
                    rResult[2*NumNodes + te_counter] = GetGeometry()[i].GetDof(PSI).EquationId();
                }
                else{
                    rResult[2*NumNodes + te_counter] = GetGeometry()[i].GetDof(LAGRANGE_MULTIPLIER_0).EquationId();
                    // int te_element_counter;
                    // if(te_counter == 0){
                    //     te_element_counter = mTeElementCounter0;
                    // }
                    // else if(te_counter == 1){
                    //     te_element_counter = mTeElementCounter1;
                    // }

                    // if(te_element_counter == 1){
                    //     rResult[2*NumNodes + te_counter] = GetGeometry()[i].GetDof(LAGRANGE_MULTIPLIER_0).EquationId();
                    // }
                    // else if(te_element_counter == 2){
                    //     rResult[2*NumNodes + te_counter] = GetGeometry()[i].GetDof(LAGRANGE_MULTIPLIER_1).EquationId();
                    // }
                    // else if(te_element_counter == 3){
                    //     rResult[2*NumNodes + te_counter] = GetGeometry()[i].GetDof(LAGRANGE_MULTIPLIER_2).EquationId();
                    // }
                    // else if(te_element_counter == 4){
                    //     rResult[2*NumNodes + te_counter] = GetGeometry()[i].GetDof(LAGRANGE_MULTIPLIER_3).EquationId();
                    // }
                    // else if(te_element_counter == 5){
                    //     rResult[2*NumNodes + te_counter] = GetGeometry()[i].GetDof(LAGRANGE_MULTIPLIER_4).EquationId();
                    // }
                    // KRATOS_ERROR_IF(te_element_counter > 5)
                    //     << this->Id() << "TRAILING EDGE NODE WITH MORE THAN 5 NEIGHBOUR ELEMENTS" << std::endl;
                }

                te_counter +=1;
            }
            else{
                rResult[NumNodes + i] = GetGeometry()[i].GetDof(AUXILIARY_VELOCITY_POTENTIAL).EquationId();
            }
        }
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::GetDofListNormalElement(DofsVectorType& rElementalDofList) const
{
    for (unsigned int i = 0; i < NumNodes; i++)
        rElementalDofList[i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::GetDofListKuttaElement(DofsVectorType& rElementalDofList) const
{
    // Kutta elements have only negative part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (!GetGeometry()[i].GetValue(TRAILING_EDGE))
            rElementalDofList[i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
        else{
            rElementalDofList[i] = GetGeometry()[i].pGetDof(AUXILIARY_VELOCITY_POTENTIAL);
            //rElementalDofList[i] = GetGeometry()[i].pGetDof(PSI);
        }
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::GetDofListWakeElement(DofsVectorType& rElementalDofList)
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

    unsigned int te_counter = 0;
    // Negative part - sign is opposite to the previous case
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] < 0)
            rElementalDofList[NumNodes + i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
        else{
            //rElementalDofList[NumNodes + i] = GetGeometry()[i].pGetDof(AUXILIARY_VELOCITY_POTENTIAL);
            if(GetGeometry()[i].GetValue(TRAILING_EDGE)){
                rElementalDofList[NumNodes + i] = GetGeometry()[i].pGetDof(AUXILIARY_VELOCITY_POTENTIAL);
                //rElementalDofList[2*NumNodes + te_counter] = GetGeometry()[i].pGetDof(PSI);
                if(this->GetValue(WING_TIP)){
                    rElementalDofList[2*NumNodes + te_counter] = GetGeometry()[i].pGetDof(PSI);
                }
                else{
                    rElementalDofList[2*NumNodes + te_counter] = GetGeometry()[i].pGetDof(LAGRANGE_MULTIPLIER_0);
                    // int& te_element_counter = GetGeometry()[i].GetValue(TE_ELEMENT_COUNTER);
                    // #pragma omp critical
                    // {
                    //     te_element_counter += 1;


                    // if(te_element_counter == 1){
                    //     rElementalDofList[2*NumNodes + te_counter] = GetGeometry()[i].pGetDof(LAGRANGE_MULTIPLIER_0);
                    //     if(GetGeometry()[i].Id() == 475){
                    //         KRATOS_WATCH(this->Id())
                    //         KRATOS_WATCH(te_element_counter)
                    //     }
                    // }
                    // else if(te_element_counter == 2){
                    //     rElementalDofList[2*NumNodes + te_counter] = GetGeometry()[i].pGetDof(LAGRANGE_MULTIPLIER_1);
                    // }
                    // else if(te_element_counter == 3){
                    //     rElementalDofList[2*NumNodes + te_counter] = GetGeometry()[i].pGetDof(LAGRANGE_MULTIPLIER_2);
                    // }
                    // else if(te_element_counter == 4){
                    //     rElementalDofList[2*NumNodes + te_counter] = GetGeometry()[i].pGetDof(LAGRANGE_MULTIPLIER_3);
                    //     if(GetGeometry()[i].Id() == 475){
                    //         KRATOS_WATCH(this->Id())
                    //         KRATOS_WATCH(te_element_counter)
                    //     }
                    // }
                    // else if(te_element_counter == 5){
                    //     rElementalDofList[2*NumNodes + te_counter] = GetGeometry()[i].pGetDof(LAGRANGE_MULTIPLIER_4);
                    // }

                    // // if(GetGeometry()[i].Id() == 12236){
                    // //         KRATOS_WATCH(this->Id())
                    // //         KRATOS_WATCH(te_element_counter)
                    // // }

                    // if(this->Id() == 12236){
                    //         KRATOS_WATCH(this->Id())
                    //         KRATOS_WATCH(te_element_counter)
                    // }

                    // KRATOS_ERROR_IF(te_element_counter > 5)
                    //     << " NODE WITH ID " << GetGeometry()[i].Id()
                    //     << " BELONGING TO ELEMENT " << this-> Id()
                    //     << " HAS MORE THAN 5 NEIGHBOUR ELEMENTS, te_element_counter = " << te_element_counter << std::endl;

                    // if(te_counter == 0){
                    //     mTeElementCounter0 = te_element_counter;
                    // }
                    // else if(te_counter == 1){
                    //     mTeElementCounter1 = te_element_counter;
                    // }
                    // }

                }

                te_counter +=1;
            }
            else{
                rElementalDofList[NumNodes + i] = GetGeometry()[i].pGetDof(AUXILIARY_VELOCITY_POTENTIAL);
            }
        }
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLocalSystemNormalElement(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != NumNodes || rLeftHandSideMatrix.size2() != NumNodes)
        rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);
    if (rRightHandSideVector.size() != NumNodes)
        rRightHandSideVector.resize(NumNodes, false);
    rLeftHandSideMatrix.clear();

    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    const double free_stream_density = rCurrentProcessInfo[FREE_STREAM_DENSITY];

    noalias(rLeftHandSideMatrix) =
        data.vol * free_stream_density * prod(data.DN_DX, trans(data.DN_DX));

    // const IncompressiblePotentialFlowElement& r_this = *this;
    // const int kutta = r_this.GetValue(KUTTA);
    // if(kutta){
    //     //KRATOS_WATCH(this->Id())
    //     for (unsigned int i = 0; i < NumNodes; ++i){
    //         // The TE node takes the contribution of the subdivided element and
    //         // we do not apply the wake condition on the TE node
    //         if (GetGeometry()[i].GetValue(TRAILING_EDGE))// && !GetGeometry()[i].GetValue(WING_TIP))
    //         {
    //             //to_be_decoupled = false;
    //             //KRATOS_WATCH(this->Id())
    //             // KRATOS_WATCH(GetGeometry()[i].Id())
    //             // Kutta elements do not contribute to the TE node
    //             for (unsigned int j = 0; j < NumNodes; ++j){
    //                 rLeftHandSideMatrix(i, j) = 0.0;
    //             }
    //         }
    //     }
    // }

    // for (unsigned int i = 0; i < NumNodes; ++i){
    //     // The TE node takes the contribution of the subdivided element and
    //     // we do not apply the wake condition on the TE node
    //     if (GetGeometry()[i].GetValue(TRAILING_EDGE))// && !GetGeometry()[i].GetValue(WING_TIP))
    //     {
    //         // Elements do not contribute to the TE node
    //         for (unsigned int j = 0; j < NumNodes; ++j){
    //             rLeftHandSideMatrix(i, j) = 0.0;
    //         }
    //     }
    // }

    data.potentials = PotentialFlowUtilities::GetPotentialOnNormalElement<Dim,NumNodes>(*this);
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, data.potentials);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLocalSystemWakeElement(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if(this->GetValue(WING_TIP)){
        if (rLeftHandSideMatrix.size1() != 2 * NumNodes + 2 ||
            rLeftHandSideMatrix.size2() != 2 * NumNodes + 2)
            rLeftHandSideMatrix.resize(2 * NumNodes + 2, 2 * NumNodes + 2, false);
        if (rRightHandSideVector.size() != 2 * NumNodes + 2)
            rRightHandSideVector.resize(2 * NumNodes + 2, false);
    }
    else if(this->GetValue(ZERO_VELOCITY_CONDITION)){
        if (rLeftHandSideMatrix.size1() != 2 * NumNodes + 1 ||
            rLeftHandSideMatrix.size2() != 2 * NumNodes + 1)
            rLeftHandSideMatrix.resize(2 * NumNodes + 1, 2 * NumNodes + 1, false);
        if (rRightHandSideVector.size() != 2 * NumNodes + 1)
            rRightHandSideVector.resize(2 * NumNodes + 1, false);
    }
    else{
        if (rLeftHandSideMatrix.size1() != 2 * NumNodes ||
            rLeftHandSideMatrix.size2() != 2 * NumNodes)
            rLeftHandSideMatrix.resize(2 * NumNodes, 2 * NumNodes, false);
        if (rRightHandSideVector.size() != 2 * NumNodes)
            rRightHandSideVector.resize(2 * NumNodes, false);
    }

    // else if(this->GetValue(ZERO_VELOCITY_CONDITION)){
    //     if (rLeftHandSideMatrix.size1() != 2 * NumNodes + 1 ||
    //         rLeftHandSideMatrix.size2() != 2 * NumNodes + 1)
    //         rLeftHandSideMatrix.resize(2 * NumNodes + 1, 2 * NumNodes + 1, false);
    //     if (rRightHandSideVector.size() != 2 * NumNodes + 1)
    //         rRightHandSideVector.resize(2 * NumNodes + 1, false);
    // }

    // // Note that the lhs and rhs have double the size
    // if (rLeftHandSideMatrix.size1() != 2 * NumNodes ||
    //     rLeftHandSideMatrix.size2() != 2 * NumNodes)
    //     rLeftHandSideMatrix.resize(2 * NumNodes, 2 * NumNodes, false);
    // if (rRightHandSideVector.size() != 2 * NumNodes)
    //     rRightHandSideVector.resize(2 * NumNodes, false);
    rLeftHandSideMatrix.clear();
    rRightHandSideVector.clear();

    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    const double free_stream_density = rCurrentProcessInfo[FREE_STREAM_DENSITY];
    //const double penalty = 1.0/data.vol;
    //const double penalty = 1.0;
    const double penalty = rCurrentProcessInfo[WATER_PRESSURE];
    const double penalty_te = rCurrentProcessInfo[PSI];
    //const double penalty = 1.0;
    // if(penalty > psi){
    //     KRATOS_WATCH(this->Id())
    //     KRATOS_WATCH(penalty)
    // }
    // if(this->Id() == 77307){
    //     KRATOS_WATCH(penalty_te)
    //     KRATOS_WATCH(penalty)
    // }

    GetWakeDistances(data.distances);

    BoundedMatrix<double, NumNodes, NumNodes> lhs_total = ZeroMatrix(NumNodes, NumNodes);

    ComputeLHSGaussPointContribution(data.vol*free_stream_density, lhs_total, data);

    BoundedMatrix<double, Dim, Dim> condition_matrix = IdentityMatrix(Dim,Dim);
    // condition_matrix(0,0) = 1.0;
    // condition_matrix(1,1) = 0.0;
    // condition_matrix(2,2) = 1.0;

    BoundedMatrix<double, Dim, Dim> condition_matrix2 = IdentityMatrix(Dim,Dim);
    // condition_matrix2(0,0) = 1.0;
    // condition_matrix2(1,1) = 0.0;
    // condition_matrix2(2,2) = 1.0;

    auto filter = prod(condition_matrix, trans(data.DN_DX));
    auto filter2 = prod(condition_matrix2, trans(data.DN_DX));

    BoundedMatrix<double, NumNodes, NumNodes> lhs_total2 = ZeroMatrix(NumNodes, NumNodes);
    for (unsigned int i = 0; i < NumNodes; i++){
        for(unsigned int j = 0; j < NumNodes; j++){
            for(unsigned int k = 0; k < Dim; k++){
                lhs_total2(i,j) += data.vol*free_stream_density*data.DN_DX(i,k)*filter(k,j);
            }
        }
    }

    BoundedMatrix<double, NumNodes, NumNodes> lhs_total3 = ZeroMatrix(NumNodes, NumNodes);
    for (unsigned int i = 0; i < NumNodes; i++){
        for(unsigned int j = 0; j < NumNodes; j++){
            for(unsigned int k = 0; k < Dim; k++){
                lhs_total3(i,j) += data.vol*free_stream_density*data.DN_DX(i,k)*filter2(k,j);
            }
        }
    }

    // if(this->Id()==4){
    //     KRATOS_WATCH(condition_matrix)
    //     KRATOS_WATCH(condition_matrix2)
    //     KRATOS_WATCH(data.DN_DX)
    //     KRATOS_WATCH(trans(data.DN_DX))
    //     KRATOS_WATCH(filter)
    //     KRATOS_WATCH(filter2)
    //     KRATOS_WATCH(lhs_total2)
    //     KRATOS_WATCH(lhs_total3)
    // }

    // CalculateLocalSystemSubdividedElement(lhs_positive, lhs_negative, rCurrentProcessInfo);
    if (this->GetValue(WING_TIP))
    //if (this->Is(STRUCTURE))
    {
        array_1d<bool, TNumNodes> free_nodes = this->GetValue(FREE_NODES);
        // KRATOS_WATCH(this->Id())
        //KRATOS_WATCH(free_nodes)
        unsigned int number_of_trailing_edge_nodes = 0;
        unsigned int number_of_free_nodes = 0;
        BoundedMatrix<double, NumNodes, NumNodes> lhs_positive = ZeroMatrix(NumNodes, NumNodes);
        BoundedMatrix<double, NumNodes, NumNodes> lhs_negative = ZeroMatrix(NumNodes, NumNodes);

        CalculateLocalSystemSubdividedElement(lhs_positive, lhs_negative, rCurrentProcessInfo);
        unsigned int te_counter = 0;
        for (unsigned int row = 0; row < NumNodes; ++row){
            // The TE node takes the contribution of the subdivided element and
            // we do not apply the wake condition on the TE node
            if (GetGeometry()[row].GetValue(TRAILING_EDGE)){//} && free_nodes[row]){
                number_of_trailing_edge_nodes += 1;
                if(free_nodes[row]){
                    ++number_of_free_nodes;
                }
                for (unsigned int j = 0; j < NumNodes; ++j){
                    rLeftHandSideMatrix(row, j) = lhs_positive(row, j) + penalty_te*lhs_total2(row, j);//lhs_positive(row, j);
                    rLeftHandSideMatrix(row + NumNodes, j + NumNodes) = lhs_negative(row, j) + penalty_te*lhs_total2(row, j);//lhs_negative(row, j);// penalty_te*lhs_total2(row, j);//(penalty+1)*lhs_negative(row, j);//penalty*lhs_negative(row, j);//2*lhs_negative(row, j);//0.0;//lhs_negative(row, j);
                    // Off-diagonal block
                    rLeftHandSideMatrix(row , j + NumNodes) = -penalty_te*lhs_total2(row, j);//0.0;// -(penalty-1)*lhs_negative(row, j); // Side 2
                    rLeftHandSideMatrix(row + NumNodes, j) = -penalty_te*lhs_total2(row, j);//0.0;// -(penalty-1)*lhs_negative(row, j); // Side 2
                    // rLeftHandSideMatrix(2*NumNodes + te_counter, j + NumNodes) = penalty_te*lhs_total2(row, j);//lhs_negative(row, j);
                    // rLeftHandSideMatrix(2*NumNodes + te_counter, j) = -penalty_te*lhs_total2(row, j);
                    // //rLeftHandSideMatrix(2*NumNodes + te_counter, j) = -lhs_negative(row, j);
                    // rLeftHandSideMatrix(j, 2*NumNodes + te_counter) = -lhs_total2(row, j);
                    // rLeftHandSideMatrix(j + NumNodes, 2*NumNodes + te_counter) = lhs_total2(row, j);
                    // if(free_nodes[row] && number_of_free_nodes < 2){
                    //     rLeftHandSideMatrix(2*NumNodes + te_counter, j + NumNodes) = penalty_te*lhs_total2(row, j);//lhs_negative(row, j);
                    //     rLeftHandSideMatrix(2*NumNodes + te_counter, j) = -penalty_te*lhs_total2(row, j);
                    //     //rLeftHandSideMatrix(2*NumNodes + te_counter, j) = -lhs_negative(row, j);
                    //     rLeftHandSideMatrix(j, 2*NumNodes + te_counter) = -lhs_total2(row, j);
                    //     rLeftHandSideMatrix(j + NumNodes, 2*NumNodes + te_counter) = lhs_total2(row, j);
                    //     // if (!GetGeometry()[j].GetValue(TRAILING_EDGE)){
                    //     //     //rLeftHandSideMatrix(2*NumNodes + te_counter, j) = -lhs_negative(row, j);
                    //     //     rLeftHandSideMatrix(j, 2*NumNodes + te_counter) = -lhs_total2(row, j);
                    //     //     rLeftHandSideMatrix(j + NumNodes, 2*NumNodes + te_counter) = lhs_total2(row, j);
                    //     // }
                    // }
                }
                // Diagonal term
                //rLeftHandSideMatrix(2*NumNodes + te_counter, 2*NumNodes + te_counter) = lhs_total(row, row);// lhs_negative(row, row) + lhs_total(row, row);
                // if(te_counter == 0){
                //     rLeftHandSideMatrix(2*NumNodes + te_counter, 2*NumNodes + te_counter) = 2*lhs_total(row, row);
                //     rLeftHandSideMatrix(2*NumNodes + te_counter, 2*NumNodes + te_counter + 1) = -lhs_total(row, row);
                // }
                // else{
                //     rLeftHandSideMatrix(2*NumNodes + te_counter, 2*NumNodes + te_counter) = lhs_total(row, row);
                // }
                // if(this->Id()==10156){
                //     if(te_counter == 0){
                //         rLeftHandSideMatrix(2*NumNodes + te_counter, 2*NumNodes + te_counter) = 2*lhs_total(row, row);
                //         rLeftHandSideMatrix(2*NumNodes + te_counter, 2*NumNodes + te_counter + 1) = -lhs_total(row, row);
                //     }
                //     else{
                //         rLeftHandSideMatrix(2*NumNodes + te_counter, 2*NumNodes + te_counter) = lhs_total(row, row);
                //     }
                // }
                // else{
                //     rLeftHandSideMatrix(2*NumNodes + te_counter, 2*NumNodes + te_counter) = 2*lhs_total(row, row);
                //     if(te_counter == 0){
                //         rLeftHandSideMatrix(2*NumNodes + te_counter, 2*NumNodes + te_counter + 1) = -lhs_total(row, row);
                //     }
                //     else{
                //         rLeftHandSideMatrix(2*NumNodes + te_counter, 2*NumNodes + te_counter - 1) = -lhs_total(row, row);
                //     }
                // }


                // rLeftHandSideMatrix(2*NumNodes + te_counter, row + NumNodes) = -lhs_total(row, row);
                // rLeftHandSideMatrix(row + NumNodes, row + NumNodes) += lhs_total(row, row);
                // rLeftHandSideMatrix(row + NumNodes, 2*NumNodes + te_counter) -= lhs_total(row, row);
                te_counter +=1;
                if(this->Id()==10156){
                    KRATOS_WATCH(row)
                    // if(free_nodes[row]){
                    //     KRATOS_WATCH(row)
                    //     KRATOS_WATCH(GetGeometry()[row].Id())
                    // }
                }
                //KRATOS_WATCH(this->Id())
                //KRATOS_WATCH(GetGeometry()[row].Id())
            }
            // else if (GetGeometry()[row].GetValue(TRAILING_EDGE)){
            //     for (unsigned int j = 0; j < NumNodes; ++j){
            //         rLeftHandSideMatrix(row, j) = lhs_positive(row, j);
            //         rLeftHandSideMatrix(row + NumNodes, j + NumNodes) = lhs_negative(row, j);// + penalty*lhs_total(row, j);//(penalty+1)*lhs_negative(row, j);//penalty*lhs_negative(row, j);//2*lhs_negative(row, j);//0.0;//lhs_negative(row, j);
            //         // Off-diagonal block
            //         rLeftHandSideMatrix(row + NumNodes, j) = 0.0;//-penalty*lhs_total(row, j);//-penalty*lhs_negative(row, j);//0.0;// -(penalty-1)*lhs_negative(row, j); // Side 2
            //     }
            // }
            else{
                // Applying wake condition on the AUXILIARY_VELOCITY_POTENTIAL dofs
                if (data.distances[row] < 0.0){
                    //unsigned int te_counter_c = 0;
                    for (unsigned int column = 0; column < NumNodes; ++column){
                        // Diagonal block
                        rLeftHandSideMatrix(row, column) = penalty*lhs_total2(row, column);
                        rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = lhs_total(row, column);
                        // Off-diagonal block
                        rLeftHandSideMatrix(row, column + NumNodes) = -penalty*lhs_total2(row, column); // Side 1
                        // //rLeftHandSideMatrix(row + NumNodes, column) = -lhs_total(row, column);
                        // if (GetGeometry()[column].GetValue(TRAILING_EDGE)){
                        //     rLeftHandSideMatrix(row + NumNodes, column + NumNodes) += lhs_total(row, column);
                        //     rLeftHandSideMatrix(row + NumNodes, 2*NumNodes + te_counter_c) = -lhs_total(row, column);
                        //     // rLeftHandSideMatrix(row, column + NumNodes) += penalty*lhs_total2(row, column);
                        //     // rLeftHandSideMatrix(row,  2*NumNodes + te_counter_c) = -penalty*lhs_total2(row, column);
                        //     // rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = 0.0;
                        //     // rLeftHandSideMatrix(row + NumNodes, 2*NumNodes + te_counter) = lhs_total(row, column);
                        //     // rLeftHandSideMatrix(row, 2*NumNodes + te_counter) = -penalty*lhs_total2(row, column);
                        //     // rLeftHandSideMatrix(row, column + NumNodes) += penalty*lhs_total2(row, column);
                        //     te_counter_c += 1;
                        // }
                    }
                }
                else if (data.distances[row] > 0.0){
                    //unsigned int te_counter_c = 0;
                    for (unsigned int column = 0; column < NumNodes; ++column){
                        // Diagonal block
                        rLeftHandSideMatrix(row, column) = lhs_total(row, column);
                        rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = penalty*lhs_total2(row, column);
                        // Off-diagonal block
                        rLeftHandSideMatrix(row + NumNodes, column) = -penalty*lhs_total2(row, column); // Side 2
                        //rLeftHandSideMatrix(row, column + NumNodes) = -lhs_total(row, column);
                        // if (GetGeometry()[column].GetValue(TRAILING_EDGE)){
                        //     // Applying equality in phi dof
                        //     // rLeftHandSideMatrix(row, column + NumNodes) = lhs_total(row, column);
                        //     // rLeftHandSideMatrix(row, 2*NumNodes + te_counter_c) = -lhs_total(row, column);
                        //     //Aplying equality in theta dof
                        //     // rLeftHandSideMatrix(row + NumNodes, column + NumNodes) += penalty*lhs_total2(row, column);
                        //     // rLeftHandSideMatrix(row + NumNodes,  2*NumNodes + te_counter_c) = -penalty*lhs_total2(row, column);

                        //     // rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = 0.0;
                        //     // rLeftHandSideMatrix(row + NumNodes, 2*NumNodes + te_counter) = penalty*lhs_total2(row, column);
                        //     te_counter_c += 1;
                        // }
                    }

                }
            }
        }
         KRATOS_WARNING_IF("element", number_of_free_nodes > 1)
        << "NUMBER OF FREE NODES LARGER THAN 1 " << this->Id() << std::endl;
    }
    else if (this->GetValue(ZERO_VELOCITY_CONDITION))
    {
        BoundedMatrix<double, NumNodes, NumNodes> lhs_positive = ZeroMatrix(NumNodes, NumNodes);
        BoundedMatrix<double, NumNodes, NumNodes> lhs_negative = ZeroMatrix(NumNodes, NumNodes);

        CalculateLocalSystemSubdividedElement(lhs_positive, lhs_negative, rCurrentProcessInfo);
        unsigned int number_of_free_nodes = 0;
        unsigned int number_of_trailing_edge_nodes = 0;
        unsigned int te_counter = 0;
        for (unsigned int row = 0; row < NumNodes; ++row){
            // The TE node takes the contribution of the subdivided element and
            // we do not apply the wake condition on the TE node
            if (GetGeometry()[row].GetValue(TRAILING_EDGE)){//} && free_nodes[row]){
                number_of_trailing_edge_nodes += 1;
                //++number_of_free_nodes;
                for (unsigned int j = 0; j < NumNodes; ++j){
                    rLeftHandSideMatrix(row, j) = lhs_positive(row, j) + penalty_te*lhs_total2(row, j);//lhs_positive(row, j);
                    rLeftHandSideMatrix(row + NumNodes, j + NumNodes) = lhs_negative(row, j) + penalty_te*lhs_total2(row, j);//lhs_negative(row, j) + penalty_te*lhs_total2(row, j);//(penalty+1)*lhs_negative(row, j);//penalty*lhs_negative(row, j);//2*lhs_negative(row, j);//0.0;//lhs_negative(row, j);
                    // Off-diagonal block
                    rLeftHandSideMatrix(row , j + NumNodes) = -penalty_te*lhs_total2(row, j);//0.0;// -(penalty-1)*lhs_negative(row, j); // Side 2
                    rLeftHandSideMatrix(row + NumNodes, j) = -penalty_te*lhs_total2(row, j);//0.0;// -(penalty-1)*lhs_negative(row, j); // Side 2
                    // rLeftHandSideMatrix(2*NumNodes + te_counter, j + NumNodes) = lhs_total2(row, j);//lhs_negative(row, j);
                    // rLeftHandSideMatrix(2*NumNodes + te_counter, j) = -lhs_total2(row, j);
                    // //rLeftHandSideMatrix(2*NumNodes + te_counter, j) = -lhs_negative(row, j);
                    // rLeftHandSideMatrix(j, 2*NumNodes + te_counter) = -lhs_total2(row, j);
                    // rLeftHandSideMatrix(j + NumNodes, 2*NumNodes + te_counter) = lhs_total2(row, j);
                }
                // Diagonal term
                //rLeftHandSideMatrix(2*NumNodes + te_counter, 2*NumNodes + te_counter) = lhs_total(row, row);//lhs_negative(row, row) + lhs_total(row, row);
                //rLeftHandSideMatrix(2*NumNodes + te_counter, row + NumNodes) = -lhs_total(row, row);
                // rLeftHandSideMatrix(row + NumNodes, row + NumNodes) += lhs_total(row, row);
                // rLeftHandSideMatrix(row + NumNodes, 2*NumNodes + te_counter) -= lhs_total(row, row);
                //te_counter +=1;
                //KRATOS_WATCH(this->Id())
                //KRATOS_WATCH(GetGeometry()[row].Id())
            }
            // else if (GetGeometry()[row].GetValue(TRAILING_EDGE)){
            //     for (unsigned int j = 0; j < NumNodes; ++j){
            //         rLeftHandSideMatrix(row, j) = lhs_positive(row, j);
            //         rLeftHandSideMatrix(row + NumNodes, j + NumNodes) = lhs_negative(row, j);// + penalty*lhs_total(row, j);//(penalty+1)*lhs_negative(row, j);//penalty*lhs_negative(row, j);//2*lhs_negative(row, j);//0.0;//lhs_negative(row, j);
            //         // Off-diagonal block
            //         rLeftHandSideMatrix(row + NumNodes, j) = 0.0;//-penalty*lhs_total(row, j);//-penalty*lhs_negative(row, j);//0.0;// -(penalty-1)*lhs_negative(row, j); // Side 2
            //     }
            // }
            else{
                // Applying wake condition on the AUXILIARY_VELOCITY_POTENTIAL dofs
                if (data.distances[row] < 0.0){
                    //unsigned int te_counter_c = 0;
                    for (unsigned int column = 0; column < NumNodes; ++column){
                        // Diagonal block
                        rLeftHandSideMatrix(row, column) = penalty*lhs_total2(row, column);
                        rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = lhs_total(row, column);
                        // Off-diagonal block
                        rLeftHandSideMatrix(row, column + NumNodes) = -penalty*lhs_total2(row, column); // Side 1
                        // //rLeftHandSideMatrix(row + NumNodes, column) = -lhs_total(row, column);
                        // if (GetGeometry()[column].GetValue(TRAILING_EDGE)){
                        //     rLeftHandSideMatrix(row + NumNodes, column + NumNodes) += lhs_total(row, column);
                        //     rLeftHandSideMatrix(row + NumNodes, 2*NumNodes + te_counter_c) = -lhs_total(row, column);
                        //     // rLeftHandSideMatrix(row, column + NumNodes) += penalty*lhs_total2(row, column);
                        //     // rLeftHandSideMatrix(row,  2*NumNodes + te_counter_c) = -penalty*lhs_total2(row, column);
                        //     // rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = 0.0;
                        //     // rLeftHandSideMatrix(row + NumNodes, 2*NumNodes + te_counter) = lhs_total(row, column);
                        //     // rLeftHandSideMatrix(row, 2*NumNodes + te_counter) = -penalty*lhs_total2(row, column);
                        //     // rLeftHandSideMatrix(row, column + NumNodes) += penalty*lhs_total2(row, column);
                        //     te_counter_c += 1;
                        // }
                    }
                }
                else if (data.distances[row] > 0.0){
                    //unsigned int te_counter_c = 0;
                    for (unsigned int column = 0; column < NumNodes; ++column){
                        // Diagonal block
                        rLeftHandSideMatrix(row, column) = lhs_total(row, column);
                        rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = penalty*lhs_total2(row, column);
                        // Off-diagonal block
                        rLeftHandSideMatrix(row + NumNodes, column) = -penalty*lhs_total2(row, column); // Side 2
                        //rLeftHandSideMatrix(row, column + NumNodes) = -lhs_total(row, column);
                        // if (GetGeometry()[column].GetValue(TRAILING_EDGE)){
                        //     // Applying equality in phi dof
                        //     // rLeftHandSideMatrix(row, column + NumNodes) = lhs_total(row, column);
                        //     // rLeftHandSideMatrix(row, 2*NumNodes + te_counter_c) = -lhs_total(row, column);
                        //     //Aplying equality in theta dof
                        //     rLeftHandSideMatrix(row + NumNodes, column + NumNodes) += penalty*lhs_total2(row, column);
                        //     rLeftHandSideMatrix(row + NumNodes,  2*NumNodes + te_counter_c) = -penalty*lhs_total2(row, column);

                        //     // rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = 0.0;
                        //     // rLeftHandSideMatrix(row + NumNodes, 2*NumNodes + te_counter) = penalty*lhs_total2(row, column);
                        //     te_counter_c += 1;
                        // }
                    }

                }
            }
        }
        KRATOS_WARNING_IF("element", number_of_free_nodes > 1)
        << "NUMBER OF FREE NODES LARGER THAN 1 " << this->Id() << std::endl;
        // if(number_of_trailing_edge_nodes > 1){
        //     KRATOS_WATCH(this->Id())
        // }
    }
    else{
        for (unsigned int row = 0; row < NumNodes; ++row){

            // if(GetGeometry()[row].GetValue(WING_TIP) && data.distances[row] < 0.0){
            //     for (unsigned int j = 0; j < NumNodes; ++j){
            //         rLeftHandSideMatrix(row, j) = lhs_total(row, j);//lhs_positive(row, j);
            //         rLeftHandSideMatrix(row + NumNodes, j + NumNodes) = lhs_total(row, j);//0.0;//lhs_negative(row, j);
            //         // Off-diagonal block
            //         rLeftHandSideMatrix(row + NumNodes, j) = 0.0;//-lhs_negative(row, j); // Side 2
            //     }
            // }
            // else{
            // Filling the diagonal blocks (i.e. decoupling upper and lower dofs)
            // for (unsigned int column = 0; column < NumNodes; ++column){
            //     rLeftHandSideMatrix(row, column) = lhs_total3(row, column);
            //     rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = lhs_total3(row, column);
            // }

            // Applying wake condition on the AUXILIARY_VELOCITY_POTENTIAL dofs
            if (data.distances[row] < 0.0){
                for (unsigned int column = 0; column < NumNodes; ++column){
                    // Diagonal block
                    rLeftHandSideMatrix(row, column) = penalty*lhs_total2(row, column);
                    rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = lhs_total(row, column);
                    // Off-diagonal block
                    rLeftHandSideMatrix(row, column + NumNodes) = -penalty*lhs_total2(row, column); // Side 1
                }
            }
            else if (data.distances[row] > 0.0){
                for (unsigned int column = 0; column < NumNodes; ++column){
                    // Diagonal block
                    rLeftHandSideMatrix(row, column) = lhs_total(row, column);
                    rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = penalty*lhs_total2(row, column);
                    // Off-diagonal block
                    rLeftHandSideMatrix(row + NumNodes, column) = -penalty*lhs_total2(row, column); // Side 2
                }
            }
            //}
        }
    }

    // if (this->GetValue(WING_TIP)){
    //     const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    //     array_1d<double, Dim> v = ZeroVector(Dim);
    //     for (unsigned int i = 0; i < Dim; i++){
    //         v[i] = free_stream_velocity[i];
    //     }
    //     v[0] = 20.3999;
    //     v[1] = -0.0270314;
    //     KRATOS_WATCH(rLeftHandSideMatrix)
    //     KRATOS_WATCH(rRightHandSideVector)
    //     BoundedVector<double, TNumNodes> right_hand_side = ZeroVector(NumNodes);
    //     for (unsigned int i = 0; i< NumNodes; i++){
    //         if(GetGeometry()[i].GetValue(WING_TIP)){
    //             if (data.distances[i] < 0.0){
    //                 for (unsigned int j = 0; j < Dim; j++){
    //                     rRightHandSideVector[i] += data.vol*free_stream_density*data.DN_DX(i,j)*v[j];
    //                     //right_hand_side[i] =
    //                     //rLeftHandSideMatrix(i, j) += lhs_total(i, j);
    //                 }
    //             }
    //             // else{
    //             //     for (unsigned int j = 0; j < Dim; j++){
    //             //         rRightHandSideVector[i + NumNodes] += data.DN_DX(i,j)*v[j];
    //             //         rLeftHandSideMatrix(i + NumNodes, j + NumNodes) += lhs_total(i, j);
    //             //     }
    //             // }

    //         }
    //     }
    //     KRATOS_WATCH(rRightHandSideVector)
    //     KRATOS_WATCH(rLeftHandSideMatrix)
    //     KRATOS_WATCH(data.vol)
    // }


    //std::fixed;
    std::cout.precision(5);
    std::cout << std::scientific;
    std::cout << std::showpos;
    if(this->Id()==74756 || this->Id()==10156){
        std::cout << std::endl;
        KRATOS_WATCH(this->Id())
        for(unsigned int row = 0; row < rLeftHandSideMatrix.size1(); ++row){
            for(unsigned int column = 0; column < rLeftHandSideMatrix.size2(); column++){
                if(column == 3 || column == 7){
                    std::cout << " " << rLeftHandSideMatrix(row, column) << " |";
                }
                else{
                    std::cout << " " << rLeftHandSideMatrix(row, column) << " ";
                }
            }

            std::cout << std::endl;
            // for(unsigned int j = 0; j < 14*3; j++){
            //     std::cout << " " << std::endl;
            // }
            if(row ==3|| row == 7){
                for(unsigned int j = 0; j < 15*rLeftHandSideMatrix.size1(); j++){
                std::cout << "_" ;
                }
                std::cout << " " << std::endl;
            }
            else{
                for(unsigned int i = 0; i < 3; i++){
                    for(unsigned int j = 0; j < 14*4; j++){
                        std::cout << " " ;
                    }
                    if(i != 2){
                        std::cout << "|" ;
                    }
                }
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    // BoundedVector<double, 2*NumNodes> split_element_values;
    // split_element_values = PotentialFlowUtilities::GetPotentialOnWakeElement<Dim, NumNodes>(*this, data.distances);
    // noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, split_element_values);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLocalSystemSubdividedElement(
    BoundedMatrix<double, NumNodes, NumNodes>& lhs_positive,
    BoundedMatrix<double, NumNodes, NumNodes>& lhs_negative,
    const ProcessInfo& rCurrentProcessInfo)
{
    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    const double free_stream_density = rCurrentProcessInfo[FREE_STREAM_DENSITY];

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

    BoundedMatrix<double, Dim, Dim> condition_matrix = IdentityMatrix(Dim,Dim);
    // condition_matrix(0,0) = 1.0;
    // condition_matrix(1,1) = 0.0;
    // condition_matrix(2,2) = 1.0;

    auto filter = prod(condition_matrix, trans(data.DN_DX));

    // Compute the lhs and rhs that would correspond to it being divided
    for (unsigned int i = 0; i < nsubdivisions; ++i)
    {
        if (PartitionsSign[i] > 0)
            ComputeLHSGaussPointContribution(Volumes[i]*free_stream_density, lhs_positive, data);
        else{
            ComputeLHSGaussPointContribution(Volumes[i]*free_stream_density, lhs_negative, data);
            // for (unsigned int m = 0; m < NumNodes; m++){
            //     for(unsigned int j = 0; j < NumNodes; j++){
            //         for(unsigned int k = 0; k < Dim; k++){
            //             lhs_negative(m,j) += Volumes[i]*free_stream_density*data.DN_DX(m,k)*filter(k,j);
            //         }
            //     }
            // }
        }
            //ComputeLHSGaussPointContribution(Volumes[i]*free_stream_density, lhs_negative, data);
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::ComputeLHSGaussPointContribution(
    const double weight,
    BoundedMatrix<double, NumNodes, NumNodes>& lhs,
    const ElementalData<NumNodes, Dim>& data) const
{
    noalias(lhs) += weight * prod(data.DN_DX, trans(data.DN_DX));
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::AssignLocalSystemSubdividedElement(
    MatrixType& rLeftHandSideMatrix,
    BoundedMatrix<double, NumNodes, NumNodes>& lhs_positive,
    BoundedMatrix<double, NumNodes, NumNodes>& lhs_negative,
    BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
    const ElementalData<NumNodes, Dim>& data,
    const ProcessInfo& rCurrentProcessInfo) const
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
            AssignLocalSystemWakeNode(rLeftHandSideMatrix, lhs_total, data, i, rCurrentProcessInfo);
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::AssignLocalSystemWakeElement(
    MatrixType& rLeftHandSideMatrix,
    BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
    const ElementalData<NumNodes, Dim>& data,
    const ProcessInfo& rCurrentProcessInfo) const
{
    for (unsigned int row = 0; row < NumNodes; ++row)
        AssignLocalSystemWakeNode(rLeftHandSideMatrix, lhs_total, data, row, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::AssignLocalSystemWakeNode(
    MatrixType& rLeftHandSideMatrix,
    BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
    const ElementalData<NumNodes, Dim>& data,
    unsigned int& row,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const double free_stream_density = rCurrentProcessInfo[FREE_STREAM_DENSITY];

    BoundedMatrix<double, 3, 3> condition_matrix = ZeroMatrix(3,3);
    for (unsigned int i = 0; i < 3; i++){
        condition_matrix(i,i) = 1.0;
    }

    auto filter = prod(condition_matrix, trans(data.DN_DX));

    if(this->Id()==4){
        KRATOS_WATCH(condition_matrix)
        KRATOS_WATCH(data.DN_DX)
        KRATOS_WATCH(trans(data.DN_DX))
        KRATOS_WATCH(filter)
    }

    BoundedMatrix<double, NumNodes, NumNodes> lhs_total2 =
        data.vol * free_stream_density * prod(data.DN_DX, trans(data.DN_DX));

    // Filling the diagonal blocks (i.e. decoupling upper and lower dofs)
    for (unsigned int column = 0; column < NumNodes; ++column)
    {
        rLeftHandSideMatrix(row, column) = lhs_total2(row, column);
        rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = lhs_total2(row, column);
    }

    // Applying wake condition on the AUXILIARY_VELOCITY_POTENTIAL dofs
    if (data.distances[row] < 0.0)
        for (unsigned int column = 0; column < NumNodes; ++column)
            rLeftHandSideMatrix(row, column + NumNodes) = -lhs_total2(row, column); // Side 1
    else if (data.distances[row] > 0.0)
        for (unsigned int column = 0; column < NumNodes; ++column)
            rLeftHandSideMatrix(row + NumNodes, column) = -lhs_total2(row, column); // Side 2
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::ComputePotentialJump(const ProcessInfo& rCurrentProcessInfo)
{
    const array_1d<double, 3> free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    const double free_stream_velocity_norm = sqrt(inner_prod(free_stream_velocity, free_stream_velocity));
    const double reference_chord = rCurrentProcessInfo[REFERENCE_CHORD];

    array_1d<double, NumNodes> distances;
    GetWakeDistances(distances);

    auto r_geometry = GetGeometry();
    for (unsigned int i = 0; i < NumNodes; i++){
        double aux_potential = r_geometry[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL);
        double potential = r_geometry[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
        double potential_jump = aux_potential - potential;

        if (distances[i] > 0){
            r_geometry[i].SetLock();
            r_geometry[i].SetValue(POTENTIAL_JUMP, - (2.0 * potential_jump) / (free_stream_velocity_norm * reference_chord));
            r_geometry[i].UnSetLock();
        }
        else{
            r_geometry[i].SetLock();
            r_geometry[i].SetValue(POTENTIAL_JUMP, (2.0 * potential_jump) / (free_stream_velocity_norm * reference_chord));
            r_geometry[i].UnSetLock();
        }
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::ComputeElementInternalEnergy()
{
    double internal_energy = 0.0;
    array_1d<double, Dim> velocity;

    const IncompressiblePotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) // Normal element (non-wake) - eventually an embedded
        velocity = PotentialFlowUtilities::ComputeVelocityNormalElement<Dim,NumNodes>(*this);
    else // Wake element
        velocity = PotentialFlowUtilities::ComputeVelocityUpperWakeElement<Dim,NumNodes>(*this);

    internal_energy = 0.5 * inner_prod(velocity, velocity);
    this->SetValue(INTERNAL_ENERGY, std::abs(internal_energy));
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
template class IncompressiblePotentialFlowElement<3, 4>;

} // namespace Kratos
