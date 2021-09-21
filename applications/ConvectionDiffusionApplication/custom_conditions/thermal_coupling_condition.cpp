// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External libraries

// Project includes
#include "includes/convection_diffusion_settings.h"
#include "utilities/integration_utilities.h"

// Application includes
#include "thermal_coupling_condition.h"
#include "convection_diffusion_application_variables.h"

namespace Kratos
{

// Public Life Cycle //////////////////////////////////////////////////////////

template<std::size_t TDim, std::size_t TNumNodes>
ThermalCouplingCondition<TDim,TNumNodes>::ThermalCouplingCondition(
    IndexType NewId,
    Geometry< Node<3> >::Pointer pGeometry)
    : Condition(NewId,pGeometry)
{
}

template<std::size_t TDim, std::size_t TNumNodes>
ThermalCouplingCondition<TDim,TNumNodes>::ThermalCouplingCondition(
    IndexType NewId,
    Geometry< Node<3> >::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Condition(NewId,pGeometry,pProperties)
{
}

template<std::size_t TDim, std::size_t TNumNodes>
ThermalCouplingCondition<TDim,TNumNodes>::~ThermalCouplingCondition()
{
}

// Public Operations //////////////////////////////////////////////////////////

template<std::size_t TDim, std::size_t TNumNodes>
Condition::Pointer ThermalCouplingCondition<TDim,TNumNodes>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<ThermalCouplingCondition<TDim,TNumNodes>>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template<std::size_t TDim, std::size_t TNumNodes>
Condition::Pointer ThermalCouplingCondition<TDim,TNumNodes>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<ThermalCouplingCondition<TDim,TNumNodes>>(NewId, pGeom, pProperties);
}

template<std::size_t TDim, std::size_t TNumNodes>
void ThermalCouplingCondition<TDim,TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Get unknown variable from convection diffusion settings
    const ProcessInfo& r_process_info = rCurrentProcessInfo; // To ensure that Gets are threadsafe
    ConvectionDiffusionSettings &r_conv_diff_settings = *(r_process_info[CONVECTION_DIFFUSION_SETTINGS]);
    const Variable<double> &r_unknown_var = r_conv_diff_settings.GetUnknownVariable();

    // Resize the equation ids. vector
    if (rResult.size() != TNumNodes) {
        rResult.resize(TNumNodes, false);
    }

    // Fill the equation ids. vector from the condition DOFs
    const auto &r_geometry = this->GetGeometry();
    for (unsigned int i = 0; i < TNumNodes; ++i){
        rResult[i] = r_geometry[i].GetDof(r_unknown_var).EquationId();
    }

    KRATOS_CATCH("")
}

template<std::size_t TDim, std::size_t TNumNodes>
void ThermalCouplingCondition<TDim,TNumNodes>::GetDofList(
    DofsVectorType& rConditionDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Get unknown variable from convection diffusion settings
    const ProcessInfo &r_process_info = rCurrentProcessInfo; // To ensure that Gets are threadsafe
    ConvectionDiffusionSettings &r_conv_diff_settings = *(r_process_info[CONVECTION_DIFFUSION_SETTINGS]);
    const Variable<double> &r_unknown_var = r_conv_diff_settings.GetUnknownVariable();

    // Resize the DOFs vector
    if (rConditionDofList.size() != TNumNodes){
        rConditionDofList.resize(TNumNodes);
    }

    // Fill the DOFs vector from the condition nodes
    const auto &r_geometry = this->GetGeometry();
    for (unsigned int i = 0; i < TNumNodes; ++i){
        rConditionDofList[i] = r_geometry[i].pGetDof(r_unknown_var);
    }

    KRATOS_CATCH("")
}

template<std::size_t TDim, std::size_t TNumNodes>
void ThermalCouplingCondition<TDim,TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

template<std::size_t TDim, std::size_t TNumNodes>
void ThermalCouplingCondition<TDim,TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check (and resize) LHS matrix
    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes) {
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
    }

    // Get the transmission coefficient from the properties container
    const auto &r_prop = this->GetProperties();
    const double transfer_coefficient = r_prop[TRANSFER_COEFFICIENT];

    // Calculate the nodal integration weight
    const double nodal_weight = CalculateNodalIntegrationWeight();

    // Fill the LHS matrix
    FillLeftHandSideMatrix(transfer_coefficient*nodal_weight, rLeftHandSideMatrix);

    KRATOS_CATCH("")
}

template<std::size_t TDim, std::size_t TNumNodes>
void ThermalCouplingCondition<TDim,TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check (and resize) RHS vector
    if (rRightHandSideVector.size() != TNumNodes) {
        rRightHandSideVector.resize(TNumNodes, false);
    }

    // Get the transmission coefficient from the properties container
    const auto &r_prop = this->GetProperties();
    const double transfer_coefficient = r_prop[TRANSFER_COEFFICIENT];

    // Calculate the nodal integration weight
    const double nodal_weight = CalculateNodalIntegrationWeight();

    // Get unknown variable
    ConvectionDiffusionSettings &r_conv_diff_settings = *(rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS]);
    const Variable<double> &r_unknown_var = r_conv_diff_settings.GetUnknownVariable();

    // Fill the LHS matrix
    FillRightHandSideVector(transfer_coefficient*nodal_weight, r_unknown_var, rRightHandSideVector);

    KRATOS_CATCH("")
}

template<std::size_t TDim, std::size_t TNumNodes>
int ThermalCouplingCondition<TDim,TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Note that we intentionally avoid calling the base class to skip the geometry domain size check
    KRATOS_ERROR_IF( this->Id() < 1 ) << "Condition found with Id " << this->Id() << std::endl;

    return 0;

    KRATOS_CATCH("")
}

// Input and Output ///////////////////////////////////////////////////////////

template<std::size_t TDim, std::size_t TNumNodes>
std::string ThermalCouplingCondition<TDim,TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "ThermalCouplingCondition #" << Id();
    return buffer.str();
}

template<std::size_t TDim, std::size_t TNumNodes>
void ThermalCouplingCondition<TDim,TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "ThermalCouplingCondition #" << Id();
}

template<std::size_t TDim, std::size_t TNumNodes>
void ThermalCouplingCondition<TDim,TNumNodes>::PrintData(std::ostream& rOStream) const
{
    rOStream << "ThermalCouplingCondition #" << Id() << std::endl;
    this->GetGeometry().PrintData(rOStream);
}


// Finite element functions ///////////////////////////////////////////////////

template<>
double ThermalCouplingCondition<2,4>::CalculateNodalIntegrationWeight() const
{
    const auto& r_geom = GetGeometry();
    return 0.5*norm_2(r_geom[0].Coordinates()-r_geom[1].Coordinates());
}

template<>
double ThermalCouplingCondition<3,6>::CalculateNodalIntegrationWeight() const
{
    const auto& r_geom = GetGeometry();
    const double a = norm_2(r_geom[0].Coordinates()-r_geom[1].Coordinates());
    const double b = norm_2(r_geom[1].Coordinates()-r_geom[2].Coordinates());
    const double c = norm_2(r_geom[2].Coordinates()-r_geom[0].Coordinates());
    const double s = 0.5*(a+b+c);
    return std::sqrt(s*(s-a)*(s-b)*(s-c))/3.0;
}

template<>
double ThermalCouplingCondition<3,8>::CalculateNodalIntegrationWeight() const
{
    const auto& r_geom = GetGeometry();

    // First subtriangle area
    const double a_1 = norm_2(r_geom[0].Coordinates()-r_geom[1].Coordinates());
    const double b_1 = norm_2(r_geom[1].Coordinates()-r_geom[2].Coordinates());
    const double c_1 = norm_2(r_geom[2].Coordinates()-r_geom[0].Coordinates());
    const double s_1 = 0.5*(a_1+b_1+c_1);
    const double area_1 = std::sqrt(s_1*(s_1-a_1)*(s_1-b_1)*(s_1-c_1));

    // Second subtriangle area
    const double a_2 = norm_2(r_geom[0].Coordinates()-r_geom[2].Coordinates());
    const double b_2 = norm_2(r_geom[2].Coordinates()-r_geom[3].Coordinates());
    const double c_2 = norm_2(r_geom[3].Coordinates()-r_geom[0].Coordinates());
    const double s_2 = 0.5*(a_2+b_2+c_2);
    const double area_2 = std::sqrt(s_2*(s_2-a_2)*(s_2-b_2)*(s_2-c_2));

    return 0.25*(area_1+area_2);
}

template<>
void ThermalCouplingCondition<2,4>::FillLeftHandSideMatrix(
    const double LeftHandSideCoefficient,
    MatrixType& rLeftHandSideMatrix) const
{
    // Fill the diagonal LHS matrix terms
    noalias(rLeftHandSideMatrix) = LeftHandSideCoefficient * IdentityMatrix(4,4);

    // Fill the off-diagonal LHS matrix terms
    rLeftHandSideMatrix(0,2) = -LeftHandSideCoefficient;
    rLeftHandSideMatrix(1,3) = -LeftHandSideCoefficient;
    rLeftHandSideMatrix(2,0) = -LeftHandSideCoefficient;
    rLeftHandSideMatrix(3,1) = -LeftHandSideCoefficient;
}

template<>
void ThermalCouplingCondition<3,6>::FillLeftHandSideMatrix(
    const double LeftHandSideCoefficient,
    MatrixType& rLeftHandSideMatrix) const
{
    // Fill the diagonal LHS matrix terms
    noalias(rLeftHandSideMatrix) = LeftHandSideCoefficient * IdentityMatrix(6,6);

    // Fill the off-diagonal LHS matrix terms
    rLeftHandSideMatrix(0,3) = -LeftHandSideCoefficient;
    rLeftHandSideMatrix(1,4) = -LeftHandSideCoefficient;
    rLeftHandSideMatrix(2,5) = -LeftHandSideCoefficient;
    rLeftHandSideMatrix(3,0) = -LeftHandSideCoefficient;
    rLeftHandSideMatrix(4,1) = -LeftHandSideCoefficient;
    rLeftHandSideMatrix(5,2) = -LeftHandSideCoefficient;
}

template<>
void ThermalCouplingCondition<3,8>::FillLeftHandSideMatrix(
    const double LeftHandSideCoefficient,
    MatrixType& rLeftHandSideMatrix) const
{
    // Fill the diagonal LHS matrix terms
    noalias(rLeftHandSideMatrix) = LeftHandSideCoefficient * IdentityMatrix(6,6);

    // Fill the off-diagonal LHS matrix terms
    rLeftHandSideMatrix(0,4) = -LeftHandSideCoefficient;
    rLeftHandSideMatrix(1,5) = -LeftHandSideCoefficient;
    rLeftHandSideMatrix(2,6) = -LeftHandSideCoefficient;
    rLeftHandSideMatrix(3,7) = -LeftHandSideCoefficient;
    rLeftHandSideMatrix(4,0) = -LeftHandSideCoefficient;
    rLeftHandSideMatrix(5,1) = -LeftHandSideCoefficient;
    rLeftHandSideMatrix(6,2) = -LeftHandSideCoefficient;
    rLeftHandSideMatrix(6,3) = -LeftHandSideCoefficient;
}

template<>
void ThermalCouplingCondition<2,4>::FillRightHandSideVector(
    const double RightHandSideCoefficient,
    const Variable<double>& rUnknownVariable,
    VectorType& rRightHandSideVector) const
{
    // Initialize RHS vector
    noalias(rRightHandSideVector) = ZeroVector(4);

    // Fill the RHS vector terms
    const auto &r_geometry = this->GetGeometry();
    rRightHandSideVector(0) = RightHandSideCoefficient*(r_geometry[2].FastGetSolutionStepValue(rUnknownVariable)-r_geometry[0].FastGetSolutionStepValue(rUnknownVariable));
    rRightHandSideVector(1) = RightHandSideCoefficient*(r_geometry[3].FastGetSolutionStepValue(rUnknownVariable)-r_geometry[1].FastGetSolutionStepValue(rUnknownVariable));
    rRightHandSideVector(2) = RightHandSideCoefficient*(r_geometry[0].FastGetSolutionStepValue(rUnknownVariable)-r_geometry[2].FastGetSolutionStepValue(rUnknownVariable));
    rRightHandSideVector(3) = RightHandSideCoefficient*(r_geometry[1].FastGetSolutionStepValue(rUnknownVariable)-r_geometry[3].FastGetSolutionStepValue(rUnknownVariable));
}

template<>
void ThermalCouplingCondition<3,6>::FillRightHandSideVector(
    const double RightHandSideCoefficient,
    const Variable<double>& rUnknownVariable,
    VectorType& rRightHandSideVector) const
{
    // Initialize RHS vector
    noalias(rRightHandSideVector) = ZeroVector(6);

    // Fill the RHS vector terms
    const auto &r_geometry = this->GetGeometry();
    rRightHandSideVector(0) = RightHandSideCoefficient*(r_geometry[3].FastGetSolutionStepValue(rUnknownVariable)-r_geometry[0].FastGetSolutionStepValue(rUnknownVariable));
    rRightHandSideVector(1) = RightHandSideCoefficient*(r_geometry[4].FastGetSolutionStepValue(rUnknownVariable)-r_geometry[1].FastGetSolutionStepValue(rUnknownVariable));
    rRightHandSideVector(2) = RightHandSideCoefficient*(r_geometry[5].FastGetSolutionStepValue(rUnknownVariable)-r_geometry[2].FastGetSolutionStepValue(rUnknownVariable));
    rRightHandSideVector(3) = RightHandSideCoefficient*(r_geometry[0].FastGetSolutionStepValue(rUnknownVariable)-r_geometry[3].FastGetSolutionStepValue(rUnknownVariable));
    rRightHandSideVector(4) = RightHandSideCoefficient*(r_geometry[1].FastGetSolutionStepValue(rUnknownVariable)-r_geometry[4].FastGetSolutionStepValue(rUnknownVariable));
    rRightHandSideVector(5) = RightHandSideCoefficient*(r_geometry[2].FastGetSolutionStepValue(rUnknownVariable)-r_geometry[5].FastGetSolutionStepValue(rUnknownVariable));
}

template<>
void ThermalCouplingCondition<3,8>::FillRightHandSideVector(
    const double RightHandSideCoefficient,
    const Variable<double>& rUnknownVariable,
    VectorType& rRightHandSideVector) const
{
    // Initialize RHS vector
    noalias(rRightHandSideVector) = ZeroVector(8);

    // Fill the RHS vector terms
    const auto &r_geometry = this->GetGeometry();
    rRightHandSideVector(0) = RightHandSideCoefficient*(r_geometry[4].FastGetSolutionStepValue(rUnknownVariable)-r_geometry[0].FastGetSolutionStepValue(rUnknownVariable));
    rRightHandSideVector(1) = RightHandSideCoefficient*(r_geometry[5].FastGetSolutionStepValue(rUnknownVariable)-r_geometry[1].FastGetSolutionStepValue(rUnknownVariable));
    rRightHandSideVector(2) = RightHandSideCoefficient*(r_geometry[6].FastGetSolutionStepValue(rUnknownVariable)-r_geometry[2].FastGetSolutionStepValue(rUnknownVariable));
    rRightHandSideVector(3) = RightHandSideCoefficient*(r_geometry[7].FastGetSolutionStepValue(rUnknownVariable)-r_geometry[3].FastGetSolutionStepValue(rUnknownVariable));
    rRightHandSideVector(4) = RightHandSideCoefficient*(r_geometry[0].FastGetSolutionStepValue(rUnknownVariable)-r_geometry[4].FastGetSolutionStepValue(rUnknownVariable));
    rRightHandSideVector(5) = RightHandSideCoefficient*(r_geometry[1].FastGetSolutionStepValue(rUnknownVariable)-r_geometry[5].FastGetSolutionStepValue(rUnknownVariable));
    rRightHandSideVector(6) = RightHandSideCoefficient*(r_geometry[2].FastGetSolutionStepValue(rUnknownVariable)-r_geometry[6].FastGetSolutionStepValue(rUnknownVariable));
    rRightHandSideVector(7) = RightHandSideCoefficient*(r_geometry[3].FastGetSolutionStepValue(rUnknownVariable)-r_geometry[7].FastGetSolutionStepValue(rUnknownVariable));
}

// Serialization //////////////////////////////////////////////////////////////

template<std::size_t TDim, std::size_t TNumNodes>
ThermalCouplingCondition<TDim,TNumNodes>::ThermalCouplingCondition()
    : Condition()
{
}

template<std::size_t TDim, std::size_t TNumNodes>
void ThermalCouplingCondition<TDim,TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

template<std::size_t TDim, std::size_t TNumNodes>
void ThermalCouplingCondition<TDim,TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

// Explicit template instantiation
template class ThermalCouplingCondition<2,4>;
template class ThermalCouplingCondition<3,6>;
template class ThermalCouplingCondition<3,8>;

}
