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

#include "incompressible_potential_flow_pressure_element.h"
#include "compressible_potential_flow_application_variables.h"
#include "includes/cfd_variables.h"

namespace Kratos
{
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <int Dim, int NumNodes>
Element::Pointer IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<IncompressiblePotentialFlowPressureElement>(
        NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<IncompressiblePotentialFlowPressureElement>(
        NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::Clone(
    IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<IncompressiblePotentialFlowPressureElement>(
        NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    const IncompressiblePotentialFlowPressureElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) // Normal element (non-wake) - eventually an embedded
        CalculateLocalSystemNormalElement(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
    else // Wake element
        CalculateLocalSystemWakeElement(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    // TODO: improve speed
    Matrix tmp;
    CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    // TODO: improve speed
    VectorType tmp;
    CalculateLocalSystem(rLeftHandSideMatrix, tmp, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::EquationIdVector(
    EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    const IncompressiblePotentialFlowPressureElement& r_this = *this;
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
        // if(this->Is(STRUCTURE)){
        //     if (rResult.size() != 2 * NumNodes + 1)
        //         rResult.resize(2 * NumNodes + 1, false);
        // }
        if(this->GetValue(WING_TIP)){
            if (rResult.size() != 2 * NumNodes + 2)
                rResult.resize(2 * NumNodes + 2, false);
        }
        else if(this->GetValue(ZERO_VELOCITY_CONDITION)){
            if (rResult.size() != 2 * NumNodes + 1)
                rResult.resize(2 * NumNodes + 1, false);
        }
        else{
            if (rResult.size() != 2 * NumNodes + 0)
                rResult.resize(2 * NumNodes + 0, false);
        }

        GetEquationIdVectorWakeElement(rResult);
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                                                   ProcessInfo& CurrentProcessInfo)
{
    const IncompressiblePotentialFlowPressureElement& r_this = *this;
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
        // if(this->Is(STRUCTURE)){
        //     if (rElementalDofList.size() != 2 * NumNodes + 1)
        //         rElementalDofList.resize(2 * NumNodes + 1);
        // }
        if(this->GetValue(WING_TIP)){
            if (rElementalDofList.size() != 2 * NumNodes + 2)
                rElementalDofList.resize(2 * NumNodes + 2);
        }
        else if(this->GetValue(ZERO_VELOCITY_CONDITION)){
            if (rElementalDofList.size() != 2 * NumNodes + 1)
                rElementalDofList.resize(2 * NumNodes + 1);
        }
        else{
            if (rElementalDofList.size() != 2 * NumNodes + 0)
                rElementalDofList.resize(2 * NumNodes + 0);
        }

        GetDofListWakeElement(rElementalDofList);
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    bool active = true;
    if ((this)->IsDefined(ACTIVE))
        active = (this)->Is(ACTIVE);

    const IncompressiblePotentialFlowPressureElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake != 0 && active == true)
    {
        // CheckWakeCondition();
        ComputePotentialJump(rCurrentProcessInfo);
    }
    // ComputeElementInternalEnergy();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <int Dim, int NumNodes>
int IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
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
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::GetValueOnIntegrationPoints(
    const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == PRESSURE_COEFFICIENT)
    {
        rValues[0] = PotentialFlowUtilities::ComputeIncompressiblePressureCoefficient<Dim,NumNodes>(*this,rCurrentProcessInfo);
    }
    else if (rVariable == WAKE)
    {
        const IncompressiblePotentialFlowPressureElement& r_this = *this;
        rValues[0] = r_this.GetValue(WAKE);
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::GetValueOnIntegrationPoints(
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
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::GetValueOnIntegrationPoints(
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
std::string IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "IncompressiblePotentialFlowPressureElement #" << Id();
    return buffer.str();
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "IncompressiblePotentialFlowPressureElement #" << Id();
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::PrintData(std::ostream& rOStream) const
{
    pGetGeometry()->PrintData(rOStream);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::GetWakeDistances(array_1d<double, NumNodes>& distances) const
{
    noalias(distances) = GetValue(WAKE_ELEMENTAL_DISTANCES);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::GetEquationIdVectorNormalElement(EquationIdVectorType& rResult) const
{
    for (unsigned int i = 0; i < NumNodes; i++)
        rResult[i] = GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::GetEquationIdVectorKuttaElement(EquationIdVectorType& rResult) const
{
    // Kutta elements have only negative part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (!GetGeometry()[i].GetValue(TRAILING_EDGE))
            rResult[i] = GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
        else
            rResult[i] = GetGeometry()[i].GetDof(AUXILIARY_VELOCITY_POTENTIAL).EquationId();
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::GetEquationIdVectorWakeElement(EquationIdVectorType& rResult) const
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

    unsigned int te_counter = 0;
    if(this->Is(STRUCTURE)){
    //if(this->GetValue(WING_TIP)){
        for (unsigned int i = 0; i < NumNodes; i++){
            if(GetGeometry()[i].GetValue(TRAILING_EDGE)){
                // rResult[2*NumNodes + te_counter] = GetGeometry()[i].GetDof(LAGRANGE_MULTIPLIER_0).EquationId();
                // te_counter += 1;

                #pragma omp critical
                {
                    int te_element_counter;
                    if(te_counter == 0){
                        te_element_counter = mTeElementCounter0;
                    }
                    else if(te_counter == 1){
                        te_element_counter = mTeElementCounter1;
                    }

                    if(te_element_counter == 1){
                        rResult[2*NumNodes + te_counter] = GetGeometry()[i].GetDof(LAGRANGE_MULTIPLIER_0).EquationId();
                    }
                    else if(te_element_counter == 2){
                        rResult[2*NumNodes + te_counter] = GetGeometry()[i].GetDof(LAGRANGE_MULTIPLIER_1).EquationId();
                    }
                    else if(te_element_counter == 3){
                        rResult[2*NumNodes + te_counter] = GetGeometry()[i].GetDof(LAGRANGE_MULTIPLIER_2).EquationId();
                    }
                    else if(te_element_counter == 4){
                        rResult[2*NumNodes + te_counter] = GetGeometry()[i].GetDof(LAGRANGE_MULTIPLIER_3).EquationId();
                    }
                    else if(te_element_counter == 5){
                        rResult[2*NumNodes + te_counter] = GetGeometry()[i].GetDof(LAGRANGE_MULTIPLIER_4).EquationId();
                    }
                    else if(te_element_counter == 6){
                        rResult[2*NumNodes + te_counter] = GetGeometry()[i].GetDof(LAGRANGE_MULTIPLIER_5).EquationId();
                    }
                    else if(te_element_counter == 7){
                        rResult[2*NumNodes + te_counter] = GetGeometry()[i].GetDof(LAGRANGE_MULTIPLIER_6).EquationId();
                    }
                    else if(te_element_counter == 8){
                        rResult[2*NumNodes + te_counter] = GetGeometry()[i].GetDof(LAGRANGE_MULTIPLIER_7).EquationId();
                    }
                }
                te_counter += 1;
            }
        }

    }

}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::GetDofListNormalElement(DofsVectorType& rElementalDofList) const
{
    for (unsigned int i = 0; i < NumNodes; i++)
        rElementalDofList[i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::GetDofListKuttaElement(DofsVectorType& rElementalDofList) const
{
    // Kutta elements have only negative part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (!GetGeometry()[i].GetValue(TRAILING_EDGE))
            rElementalDofList[i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
        else
            rElementalDofList[i] = GetGeometry()[i].pGetDof(AUXILIARY_VELOCITY_POTENTIAL);
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::GetDofListWakeElement(DofsVectorType& rElementalDofList)
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

    unsigned int te_counter = 0;
    if(this->Is(STRUCTURE)){
    //if(this->GetValue(WING_TIP)){
        for (unsigned int i = 0; i < NumNodes; i++){
            if(GetGeometry()[i].GetValue(TRAILING_EDGE)){
                // rElementalDofList[2*NumNodes + te_counter] = GetGeometry()[i].pGetDof(LAGRANGE_MULTIPLIER_0);
                // te_counter += 1;

                #pragma omp critical
                {
                    int& te_element_counter = GetGeometry()[i].GetValue(TE_ELEMENT_COUNTER);
                    te_element_counter += 1;
                    if(te_element_counter == 1){
                        rElementalDofList[2*NumNodes + te_counter] = GetGeometry()[i].pGetDof(LAGRANGE_MULTIPLIER_0);
                    }
                    else if(te_element_counter == 2){
                        rElementalDofList[2*NumNodes + te_counter] = GetGeometry()[i].pGetDof(LAGRANGE_MULTIPLIER_1);
                    }
                    else if(te_element_counter == 3){
                        rElementalDofList[2*NumNodes + te_counter] = GetGeometry()[i].pGetDof(LAGRANGE_MULTIPLIER_2);
                    }
                    else if(te_element_counter == 4){
                        rElementalDofList[2*NumNodes + te_counter] = GetGeometry()[i].pGetDof(LAGRANGE_MULTIPLIER_3);
                    }
                    else if(te_element_counter == 5){
                        rElementalDofList[2*NumNodes + te_counter] = GetGeometry()[i].pGetDof(LAGRANGE_MULTIPLIER_4);
                    }
                    else if(te_element_counter == 6){
                        rElementalDofList[2*NumNodes + te_counter] = GetGeometry()[i].pGetDof(LAGRANGE_MULTIPLIER_5);
                    }
                    else if(te_element_counter == 7){
                        rElementalDofList[2*NumNodes + te_counter] = GetGeometry()[i].pGetDof(LAGRANGE_MULTIPLIER_6);
                    }
                    else if(te_element_counter == 8){
                        rElementalDofList[2*NumNodes + te_counter] = GetGeometry()[i].pGetDof(LAGRANGE_MULTIPLIER_7);
                    }

                    KRATOS_ERROR_IF(te_element_counter > 8)
                        << " NODE WITH ID " << GetGeometry()[i].Id()
                        << " BELONGING TO ELEMENT " << this-> Id()
                        << " HAS MORE THAN 8 NEIGHBOUR ELEMENTS, te_element_counter = " << te_element_counter << std::endl;

                    if(te_counter == 0){
                        mTeElementCounter0 = te_element_counter;
                    }
                    else if(te_counter == 1){
                        mTeElementCounter1 = te_element_counter;
                    }
                }
                te_counter += 1;
            }
        }
    }

}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::CalculateLocalSystemNormalElement(
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

    data.potentials = PotentialFlowUtilities::GetPotentialOnNormalElement<Dim,NumNodes>(*this);
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, data.potentials);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::CalculateLocalSystemWakeElement(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    // if(this->Is(STRUCTURE)){
    //     if (rLeftHandSideMatrix.size1() != 2 * NumNodes + 1 ||
    //         rLeftHandSideMatrix.size2() != 2 * NumNodes + 1)
    //         rLeftHandSideMatrix.resize(2 * NumNodes + 1, 2 * NumNodes + 1, false);
    //     if (rRightHandSideVector.size() != 2 * NumNodes + 1)
    //         rRightHandSideVector.resize(2 * NumNodes + 1, false);
    // }
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
        if (rLeftHandSideMatrix.size1() != 2 * NumNodes + 0 ||
            rLeftHandSideMatrix.size2() != 2 * NumNodes + 0)
            rLeftHandSideMatrix.resize(2 * NumNodes + 0, 2 * NumNodes + 0, false);
        if (rRightHandSideVector.size() != 2 * NumNodes + 0)
            rRightHandSideVector.resize(2 * NumNodes + 0, false);
    }
    rLeftHandSideMatrix.clear();
    rRightHandSideVector.clear();

    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    const double free_stream_density = rCurrentProcessInfo[FREE_STREAM_DENSITY];
    const double penalty_te = rCurrentProcessInfo[PSI];

    GetWakeDistances(data.distances);

    BoundedMatrix<double, NumNodes, NumNodes> lhs_total = ZeroMatrix(NumNodes, NumNodes);

    ComputeLHSGaussPointContribution(data.vol*free_stream_density, lhs_total, data);

    BoundedMatrix<double, Dim, Dim> condition_matrix = IdentityMatrix(Dim,Dim);
    condition_matrix(0,0) = 1.0;
    condition_matrix(1,1) = 0.0;
    condition_matrix(2,2) = 1.0;

    auto filter = prod(condition_matrix, trans(data.DN_DX));

    BoundedMatrix<double, NumNodes, NumNodes> lhs_total2 = ZeroMatrix(NumNodes, NumNodes);

    for (unsigned int i = 0; i < NumNodes; i++){
        for(unsigned int j = 0; j < NumNodes; j++){
            for(unsigned int k = 0; k < Dim; k++){
                lhs_total2(i,j) += data.vol*free_stream_density*data.DN_DX(i,k)*filter(k,j);
            }
        }
    }

    BoundedMatrix<double, Dim, Dim> condition_matrix2 = IdentityMatrix(Dim,Dim);
    condition_matrix(0,0) = 1.0;
    condition_matrix(1,1) = -1.0;
    condition_matrix(2,2) = 1.0;

    auto filter2 = prod(condition_matrix2, trans(data.DN_DX));

    BoundedMatrix<double, NumNodes, NumNodes> lhs_total3 = ZeroMatrix(NumNodes, NumNodes);

    for (unsigned int i = 0; i < NumNodes; i++){
        for(unsigned int j = 0; j < NumNodes; j++){
            for(unsigned int k = 0; k < Dim; k++){
                lhs_total3(i,j) += data.vol*free_stream_density*data.DN_DX(i,k)*filter2(k,j);
            }
        }
    }

    BoundedMatrix<double, NumNodes, NumNodes> lhs_positive = ZeroMatrix(NumNodes, NumNodes);
    BoundedMatrix<double, NumNodes, NumNodes> lhs_negative = ZeroMatrix(NumNodes, NumNodes);

    CalculateLocalSystemSubdividedElement(lhs_positive, lhs_negative, rCurrentProcessInfo);

    for (unsigned int row = 0; row < NumNodes; ++row){
        // The TE node takes the contribution of the subdivided element and
        // we do not apply the wake condition on the TE node
        if (GetGeometry()[row].GetValue(TRAILING_EDGE)){
            for (unsigned int j = 0; j < NumNodes; ++j)
            {
                rLeftHandSideMatrix(row, j) = lhs_positive(row, j);
                rLeftHandSideMatrix(row + NumNodes, j + NumNodes) = lhs_negative(row, j);
            }
        }
        else{
            // Applying wake condition on the AUXILIARY_VELOCITY_POTENTIAL dofs
            if (data.distances[row] < 0.0){
                for (unsigned int column = 0; column < NumNodes; ++column){
                    // Conservation of mass
                    rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = lhs_total(row, column);
                    // Constraint
                    // Diagonal
                    rLeftHandSideMatrix(row, column) = lhs_total2(row, column);
                    // Off diagonal
                    rLeftHandSideMatrix(row, column + NumNodes) = -lhs_total2(row, column); // Side 1
                }
            }
            else if (data.distances[row] > 0.0){
                for (unsigned int column = 0; column < NumNodes; ++column){
                    // Conservation of mass
                    rLeftHandSideMatrix(row, column) = lhs_total(row, column);
                    // Constraint
                    // Diagonal
                    rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = lhs_total2(row, column);
                    // Off diagonal
                    rLeftHandSideMatrix(row + NumNodes, column) = -lhs_total2(row, column); // Side 2
                }
            }
        }
    }

    unsigned int number_of_trailing_edge_nodes = 0;
    Vector lagrange_multipliers = ZeroVector(number_of_trailing_edge_nodes);
    if(this->GetValue(WING_TIP)){
        if(lagrange_multipliers.size() != 2){
            lagrange_multipliers.resize(2);
        }
        number_of_trailing_edge_nodes = 2;
    }
    else if(this->GetValue(ZERO_VELOCITY_CONDITION)){
        if(lagrange_multipliers.size() != 1){
            lagrange_multipliers.resize(1);
        }
        number_of_trailing_edge_nodes = 1;
    }

    double lagrange_multiplier = 5.0;
    BoundedVector<double, 2*NumNodes> split_element_values;

    array_1d<double, Dim> velocity_upper = PotentialFlowUtilities::ComputeVelocityUpperWakeElement<Dim,NumNodes>(*this);
    BoundedVector<double, NumNodes> dU2dPhi_upper = 2*data.vol*prod(data.DN_DX, velocity_upper);
    const double velocity_upper_2 = inner_prod(velocity_upper, velocity_upper);

    array_1d<double, Dim> velocity_lower = PotentialFlowUtilities::ComputeVelocityLowerWakeElement<Dim,NumNodes>(*this);
    BoundedVector<double, NumNodes> dU2dPhi_lower = 2*data.vol*prod(data.DN_DX, velocity_lower);
    const double velocity_lower_2 = inner_prod(velocity_lower, velocity_lower);

    BoundedMatrix<double, NumNodes, NumNodes> lhs_pressure_upper = ZeroMatrix(NumNodes, NumNodes);
    BoundedMatrix<double, NumNodes, NumNodes> lhs_pressure_lower = ZeroMatrix(NumNodes, NumNodes);
    BoundedVector<double, NumNodes> residual_pressure = ZeroVector(NumNodes);

    for(unsigned int row = 0; row < NumNodes; ++row){
        for(unsigned int column = 0; column < NumNodes; ++column){
            lhs_pressure_upper(row, column) = data.N[row] * dU2dPhi_upper[column];
            lhs_pressure_lower(row, column) = data.N[row] * dU2dPhi_lower[column];
        }
        residual_pressure(row) = data.N[row] * (velocity_upper_2 - velocity_lower_2) * data.vol;
    }

    unsigned int te_counter = 0;
    if(this->Is(STRUCTURE)){
    //if(this->GetValue(WING_TIP)){

        // for (unsigned int column = 0; column < NumNodes; ++column){
        //     rLeftHandSideMatrix(2*NumNodes, column) = lhs_pressure_upper(0, column);
        //     rLeftHandSideMatrix(2*NumNodes, column + NumNodes) = - lhs_pressure_lower(0, column);

        //     rLeftHandSideMatrix(column, 2*NumNodes) = lhs_pressure_upper(0, column);
        //     rLeftHandSideMatrix(column + NumNodes, 2*NumNodes) = - lhs_pressure_lower(0, column);
        // }

        // te_counter = 0;
        // for (unsigned int row = 0; row < NumNodes; ++row){
        //     if (GetGeometry()[row].GetValue(TRAILING_EDGE)){
        //         if(!GetGeometry()[row].GetValue(WING_TIP)){
        //             for(unsigned int column = 0; column < NumNodes; ++column){
        //                 rLeftHandSideMatrix(2*NumNodes + te_counter, column + NumNodes) = lhs_total2(row, column);//lhs_pressure_upper(0, column);//lhs_total2(row, column);
        //                 rLeftHandSideMatrix(2*NumNodes + te_counter, column) = -lhs_total2(row, column);//-lhs_pressure_lower(0, column);//-lhs_total2(row, column);
        //                 rLeftHandSideMatrix(column, 2*NumNodes + te_counter) = -lhs_total2(row, column);//-lhs_pressure_lower(0, column);//-lhs_total2(row, column);
        //                 rLeftHandSideMatrix(column + NumNodes, 2*NumNodes + te_counter) = lhs_total2(row, column);//lhs_pressure_upper(0, column);//lhs_total2(row, column);
        //             }
        //         }
        //         KRATOS_ERROR_IF(te_counter > 1)
        //             << " NODE WITH ID " << GetGeometry()[row].Id()
        //             << " BELONGING TO ELEMENT " << this-> Id()
        //             << " HAS MORE THAN 2 TE NODES, te_element_counter = " << te_counter << std::endl;
        //         te_counter +=1;
        //     }
        // }

        //if(this->GetValue(WING_TIP_ELEMENT) || this->Id()==94253){
        // if(this->Id()==87225){
        //     te_counter = 0;
        //     for (unsigned int row = 0; row < NumNodes; ++row){
        //         if (GetGeometry()[row].GetValue(TRAILING_EDGE)){
        //             if(!GetGeometry()[row].GetValue(WING_TIP)){
        //                 for(unsigned int column = 0; column < NumNodes; ++column){
        //                     rLeftHandSideMatrix(2*NumNodes + te_counter, column + NumNodes) = lhs_total2(row, column);//lhs_pressure_upper(0, column);//lhs_total2(row, column);
        //                     rLeftHandSideMatrix(2*NumNodes + te_counter, column) = -lhs_total2(row, column);//-lhs_pressure_lower(0, column);//-lhs_total2(row, column);
        //                     rLeftHandSideMatrix(column, 2*NumNodes + te_counter) = -lhs_total2(row, column);//-lhs_pressure_lower(0, column);//-lhs_total2(row, column);
        //                     rLeftHandSideMatrix(column + NumNodes, 2*NumNodes + te_counter) = lhs_total2(row, column);//lhs_pressure_upper(0, column);//lhs_total2(row, column);
        //                 }
        //             }
        //             te_counter +=1;
        //         }
        //     }
        // }

        // // Apply pressure constraing in only some elements
        // //if(this->Id()==79624){
        // if(this->GetValue(WING_TIP_ELEMENT) || this->Id()==94253){
        //     te_counter = 0;
        //     for (unsigned int row = 0; row < NumNodes; ++row){
        //         if (GetGeometry()[row].GetValue(TRAILING_EDGE)){
        //             if(!GetGeometry()[row].GetValue(WING_TIP)){
        //                 for(unsigned int column = 0; column < NumNodes; ++column){
        //                     rLeftHandSideMatrix(2*NumNodes + te_counter, column + NumNodes) = -lhs_pressure_lower(0, column);
        //                     rLeftHandSideMatrix(2*NumNodes + te_counter, column) = lhs_pressure_upper(0, column);
        //                     rLeftHandSideMatrix(column, 2*NumNodes + te_counter) = lhs_pressure_upper(0, column);
        //                     rLeftHandSideMatrix(column + NumNodes, 2*NumNodes + te_counter) = -lhs_pressure_lower(0, column);
        //                 }
        //                 rRightHandSideVector(2*NumNodes + te_counter) = - residual_pressure[0];
        //             }
        //             te_counter +=1;
        //         }
        //     }
        // }

        // Apply pressure constraint
        //if(this->GetValue(WING_TIP)){
        //if(this->Id()==844246){
        //if(this->GetValue(WING_TIP_ELEMENT)){// || this->Id()==94253){
        //if(this->GetValue(ZERO_VELOCITY_CONDITION)){
            te_counter = 0;
            int constraint_counter = 0;
            for (unsigned int row = 0; row < NumNodes; ++row){
                if (GetGeometry()[row].GetValue(TRAILING_EDGE)){
                    if(!GetGeometry()[row].GetValue(WING_TIP) && constraint_counter < 1){
                        for(unsigned int column = 0; column < NumNodes; ++column){
                            rLeftHandSideMatrix(2*NumNodes + te_counter, column + NumNodes) = -lhs_pressure_lower(0, column);
                            rLeftHandSideMatrix(2*NumNodes + te_counter, column) = lhs_pressure_upper(0, column);
                            rLeftHandSideMatrix(column, 2*NumNodes + te_counter) = lhs_pressure_upper(0, column);
                            rLeftHandSideMatrix(column + NumNodes, 2*NumNodes + te_counter) = -lhs_pressure_lower(0, column);
                        }
                        rRightHandSideVector(2*NumNodes + te_counter) = - residual_pressure[0];
                        constraint_counter +=1;
                    }
                    te_counter +=1;
                }
            }
        //}


        // Wing tip constraint
        te_counter = 0;
        if(this->GetValue(WING_TIP_ELEMENT)){
            for (unsigned int row = 0; row < NumNodes; ++row){
                if (GetGeometry()[row].GetValue(TRAILING_EDGE)){
                    if(GetGeometry()[row].GetValue(WING_TIP)){
                        rLeftHandSideMatrix(2*NumNodes + te_counter, row) = -lhs_total(row, row);
                        rLeftHandSideMatrix(2*NumNodes + te_counter, row + NumNodes) = lhs_total(row, row);

                        rLeftHandSideMatrix(row, 2*NumNodes + te_counter) = -lhs_total(row, row);
                        rLeftHandSideMatrix(row + NumNodes, 2*NumNodes + te_counter) = lhs_total(row, row);
                    }
                    te_counter +=1;
                }
            }
        }

        te_counter = 0;
        for (unsigned int row = 0; row < NumNodes; ++row){
            if (GetGeometry()[row].GetValue(TRAILING_EDGE)){
                // lagrange_multipliers(te_counter) = GetGeometry()[row].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_0);
                // te_counter +=1;

                int te_element_counter;
                if(te_counter == 0){
                    te_element_counter = mTeElementCounter0;
                }
                else if(te_counter == 1){
                    te_element_counter = mTeElementCounter1;
                }

                if(te_element_counter == 1){
                    lagrange_multipliers(te_counter) = GetGeometry()[row].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_0);
                }
                else if(te_element_counter == 2){
                    lagrange_multipliers(te_counter) = GetGeometry()[row].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_1);
                }
                else if(te_element_counter == 3){
                    lagrange_multipliers(te_counter) = GetGeometry()[row].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_2);
                }
                else if(te_element_counter == 4){
                    lagrange_multipliers(te_counter) = GetGeometry()[row].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_3);
                }
                else if(te_element_counter == 5){
                    lagrange_multipliers(te_counter) = GetGeometry()[row].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_4);
                }
                else if(te_element_counter == 6){
                    lagrange_multipliers(te_counter) = GetGeometry()[row].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_5);
                }
                else if(te_element_counter == 7){
                    lagrange_multipliers(te_counter) = GetGeometry()[row].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_6);
                }
                else if(te_element_counter == 8){
                    lagrange_multipliers(te_counter) = GetGeometry()[row].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_7);
                }
                te_counter +=1;
            }
        }

        split_element_values = PotentialFlowUtilities::GetPotentialOnWakeElement<Dim, NumNodes>(*this, data.distances);

        for (unsigned int row = 0; row < 2*NumNodes; ++row){
        //for (unsigned int row = 0; row < rLeftHandSideMatrix.size1(); ++row){
            // if(this->Id()==79624){
            //     KRATOS_WATCH(row)
            // }
            for (unsigned int column = 0; column < 2*NumNodes; ++column){
                rRightHandSideVector(row) -= rLeftHandSideMatrix(row, column) * split_element_values(column);
                // if(this->Id()==79624){
                //     KRATOS_WATCH(column)
                //     KRATOS_WATCH(rRightHandSideVector(row))
                // }
            }
            for (unsigned int te_node = 0; te_node < lagrange_multipliers.size(); ++te_node){
                rRightHandSideVector(row) -= rLeftHandSideMatrix(row, 2*NumNodes + te_node) * lagrange_multipliers(te_node);
                // if(this->Id()==79624){
                //     KRATOS_WATCH(te_node)
                //     KRATOS_WATCH(rRightHandSideVector(row))
                // }
            }
            //rRightHandSideVector(row) -= rLeftHandSideMatrix(row, 2*NumNodes) * lagrange_multiplier;
        }
        // for (unsigned int te_node = 0; te_node < lagrange_multipliers.size(); ++te_node){
        //     rRightHandSideVector(2*NumNodes + te_node) = - residual_pressure[0];
        // }
        //rRightHandSideVector(2*NumNodes) = - residual_pressure[0];
    }
    else{
        split_element_values = PotentialFlowUtilities::GetPotentialOnWakeElement<Dim, NumNodes>(*this, data.distances);
        noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, split_element_values);
    }

    BoundedVector<double, NumNodes> potential_jump = ZeroVector(NumNodes);
    for(unsigned int i = 0; i < NumNodes; i++){
        const double potential = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
        const double aux_potential = GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL);
        const double diff = std::abs(potential-aux_potential);
        potential_jump(i) = diff;
    }

    BoundedVector<double, 3> vel_diff = ZeroVector(3);
    for(unsigned int i = 0; i < 3; i++){
        vel_diff(i) = velocity_upper(i) - velocity_lower(i);
    }

    const double vel_2_diff = velocity_upper_2 - velocity_lower_2;

    std::cout.precision(5);
    std::cout << std::scientific;
    std::cout << std::showpos;
    // if(this->GetValue(WING_TIP_ELEMENT)){
    //     KRATOS_WATCH(this->Id())
    // }
    //if(this->Id()==12750 || this->Id()==26346 || this->Id()==265451 || this->Id()==87225){
    //if(this->Id()==79624){//} || this->Id()==79624){
    // if(this->GetValue(WING_TIP_ELEMENT) || this->Id()==94253 || this->Id()==1){
    //     std::cout << std::endl;
    //     KRATOS_WATCH(this->Id())
    //     // KRATOS_WATCH(lhs_pressure_upper)
    //     // KRATOS_WATCH(lhs_pressure_lower)
    //     KRATOS_WATCH(residual_pressure[0])
    //     // //KRATOS_WATCH(trans(data.DN_DX))
    //     //KRATOS_WATCH(rRightHandSideVector)
    //     KRATOS_WATCH(split_element_values)
    //     KRATOS_WATCH(lagrange_multipliers)
    //     // //KRATOS_WATCH(potential_jump)
    //     KRATOS_WATCH(velocity_upper)
    //     KRATOS_WATCH(velocity_lower)
    //     // KRATOS_WATCH(vel_diff)
    //     //KRATOS_WATCH(data.distances)
    //     //KRATOS_WATCH(velocity_upper_2)
    //     //KRATOS_WATCH(velocity_lower_2)
    //     KRATOS_WATCH(vel_2_diff)
    //     // for(unsigned int row = 0; row < rLeftHandSideMatrix.size1(); ++row){
    //     //     for(unsigned int column = 0; column < rLeftHandSideMatrix.size2(); column++){
    //     //         if(column == 3 || column == 7){
    //     //             std::cout << " " << rLeftHandSideMatrix(row, column) << " |";
    //     //         }
    //     //         else{
    //     //             std::cout << " " << rLeftHandSideMatrix(row, column) << " ";
    //     //         }
    //     //     }

    //     //     std::cout << std::endl;

    //     //     if(row ==3|| row == 7){
    //     //         for(unsigned int j = 0; j < 14*rLeftHandSideMatrix.size1(); j++){
    //     //         std::cout << "_" ;
    //     //         }
    //     //         std::cout << " " << std::endl;
    //     //     }
    //     //     else{
    //     //         for(unsigned int i = 0; i < 3; i++){
    //     //             for(unsigned int j = 0; j < 14*4; j++){
    //     //                 std::cout << " " ;
    //     //             }
    //     //             if(i != 2){
    //     //                 std::cout << "|" ;
    //     //             }
    //     //         }
    //     //     }
    //     //     std::cout << std::endl;
    //     // }
    //     std::cout << std::endl;
    // }

}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::CalculateLocalSystemSubdividedElement(
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

    // Compute the lhs and rhs that would correspond to it being divided
    for (unsigned int i = 0; i < nsubdivisions; ++i)
    {
        if (PartitionsSign[i] > 0)
            ComputeLHSGaussPointContribution(Volumes[i]*free_stream_density, lhs_positive, data);
        else
            ComputeLHSGaussPointContribution(Volumes[i]*free_stream_density, lhs_negative, data);
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::ComputeLHSGaussPointContribution(
    const double weight,
    BoundedMatrix<double, NumNodes, NumNodes>& lhs,
    const ElementalData<NumNodes, Dim>& data) const
{
    noalias(lhs) += weight * prod(data.DN_DX, trans(data.DN_DX));
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::AssignLocalSystemSubdividedElement(
    MatrixType& rLeftHandSideMatrix,
    BoundedMatrix<double, NumNodes, NumNodes>& lhs_positive,
    BoundedMatrix<double, NumNodes, NumNodes>& lhs_negative,
    BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
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
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::AssignLocalSystemWakeElement(
    MatrixType& rLeftHandSideMatrix,
    BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
    const ElementalData<NumNodes, Dim>& data) const
{
    for (unsigned int row = 0; row < NumNodes; ++row)
        AssignLocalSystemWakeNode(rLeftHandSideMatrix, lhs_total, data, row);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::AssignLocalSystemWakeNode(
    MatrixType& rLeftHandSideMatrix,
    BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
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
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::CheckWakeCondition() const
{
    array_1d<double, Dim> upper_wake_velocity;
    upper_wake_velocity = PotentialFlowUtilities::ComputeVelocityUpperWakeElement<Dim,NumNodes>(*this);
    const double vupnorm = inner_prod(upper_wake_velocity, upper_wake_velocity);

    array_1d<double, Dim> lower_wake_velocity;
    lower_wake_velocity = PotentialFlowUtilities::ComputeVelocityLowerWakeElement<Dim,NumNodes>(*this);
    const double vlownorm = inner_prod(lower_wake_velocity, lower_wake_velocity);

    KRATOS_WARNING_IF("IncompressibleElement", std::abs(vupnorm - vlownorm) > 0.1) << "WAKE CONDITION NOT FULFILLED IN ELEMENT # " << this->Id() << std::endl;
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::ComputePotentialJump(const ProcessInfo& rCurrentProcessInfo)
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
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::ComputeElementInternalEnergy()
{
    double internal_energy = 0.0;
    array_1d<double, Dim> velocity;

    const IncompressiblePotentialFlowPressureElement& r_this = *this;
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
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowPressureElement<Dim, NumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class IncompressiblePotentialFlowPressureElement<2, 3>;
template class IncompressiblePotentialFlowPressureElement<3, 4>;

} // namespace Kratos
