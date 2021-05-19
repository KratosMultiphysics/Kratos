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


#include "laplacian_shifted_boundary_condition.h"
#include "utilities/integration_utilities.h"
#include "includes/convection_diffusion_settings.h"

namespace Kratos
{

// Public Life Cycle //////////////////////////////////////////////////////////

LaplacianShiftedBoundaryCondition::LaplacianShiftedBoundaryCondition(IndexType NewId, Geometry< Node<3> >::Pointer pGeometry):
    Condition(NewId,pGeometry)
{
}

LaplacianShiftedBoundaryCondition::LaplacianShiftedBoundaryCondition(
    IndexType NewId,
    Geometry< Node<3> >::Pointer pGeometry,
    Properties::Pointer pProperties):
    Condition(NewId,pGeometry,pProperties)
{
}

LaplacianShiftedBoundaryCondition::~LaplacianShiftedBoundaryCondition()
{
}

// Public Operations //////////////////////////////////////////////////////////

Condition::Pointer LaplacianShiftedBoundaryCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianShiftedBoundaryCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Condition::Pointer LaplacianShiftedBoundaryCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianShiftedBoundaryCondition>(NewId, pGeom, pProperties);
}

void LaplacianShiftedBoundaryCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // // Get unknown variable from convection diffusion settings
    // const ProcessInfo& r_process_info = rCurrentProcessInfo; // To ensure that Gets are threadsafe
    // ConvectionDiffusionSettings &r_conv_diff_settings = *(r_process_info[CONVECTION_DIFFUSION_SETTINGS]);
    // const Variable<double> &r_unknown_var = r_conv_diff_settings.GetUnknownVariable();

    // // Resize the equation ids. vector
    // auto &r_geometry = this->GetGeometry();
    // const unsigned int n_nodes = r_geometry.PointsNumber();
    // if (rResult.size() != n_nodes) {
    //     rResult.resize(n_nodes, false);
    // }

    // // Fill the equation ids. vector from the condition DOFs
    // for (unsigned int i = 0; i < n_nodes; ++i){
    //     rResult[i] = r_geometry[i].GetDof(r_unknown_var).EquationId();
    // }

    KRATOS_CATCH("")
}

void LaplacianShiftedBoundaryCondition::GetDofList(
    DofsVectorType& rConditionalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // // Get unknown variable from convection diffusion settings
    // const ProcessInfo &r_process_info = rCurrentProcessInfo; // To ensure that Gets are threadsafe
    // ConvectionDiffusionSettings &r_conv_diff_settings = *(r_process_info[CONVECTION_DIFFUSION_SETTINGS]);
    // const Variable<double> &r_unknown_var = r_conv_diff_settings.GetUnknownVariable();

    // // Resize the DOFs vector
    // auto &r_geometry = this->GetGeometry();
    // const unsigned int n_nodes = r_geometry.PointsNumber();
    // if (rConditionalDofList.size() != n_nodes){
    //     rConditionalDofList.resize(n_nodes);
    // }

    // // Fill the DOFs vector from the condition nodes
    // for (unsigned int i = 0; i < n_nodes; ++i){
    //     rConditionalDofList[i] = r_geometry[i].pGetDof(r_unknown_var);
    // }

    KRATOS_CATCH("")
}

void LaplacianShiftedBoundaryCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

void LaplacianShiftedBoundaryCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check (and resize) LHS matrix
    auto &r_geometry = this->GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    if (rLeftHandSideMatrix.size1() != n_nodes || rLeftHandSideMatrix.size2() != n_nodes) {
        rLeftHandSideMatrix.resize(n_nodes, n_nodes, false);
    }

    // Set LHS to zero
    noalias(rLeftHandSideMatrix) = ZeroMatrix(n_nodes,n_nodes);

    KRATOS_CATCH("")
}

void LaplacianShiftedBoundaryCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check (and resize) RHS vector
    auto &r_geometry = this->GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    if (rRightHandSideVector.size() != n_nodes) {
        rRightHandSideVector.resize(n_nodes,false);
    }

    // Initialize RHS vector
    noalias(rRightHandSideVector) = ZeroVector(n_nodes);

    KRATOS_CATCH("")
}

// Input and Output ///////////////////////////////////////////////////////////

std::string LaplacianShiftedBoundaryCondition::Info() const
{
    std::stringstream buffer;
    buffer << "LaplacianShiftedBoundaryCondition #" << Id();
    return buffer.str();
}

void LaplacianShiftedBoundaryCondition::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "LaplacianShiftedBoundaryCondition #" << Id();
}

void LaplacianShiftedBoundaryCondition::PrintData(std::ostream& rOStream) const
{
    rOStream << "LaplacianShiftedBoundaryCondition #" << Id() << std::endl;
    this->GetGeometry().PrintData(rOStream);
}

// Serialization //////////////////////////////////////////////////////////////

LaplacianShiftedBoundaryCondition::LaplacianShiftedBoundaryCondition():
    Condition()
{
}

void LaplacianShiftedBoundaryCondition::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

void LaplacianShiftedBoundaryCondition::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

}
