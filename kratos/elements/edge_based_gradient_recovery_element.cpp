//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes


// External includes


// Project includes
#include "elements/edge_based_gradient_recovery_element.h"


namespace Kratos
{

template<std::size_t TDim>
void EdgeBasedGradientRecoveryElement<TDim>::CalculateLocalSystem(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRigthHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Check sizes
    if (rRigthHandSideVector.size() != LocalSize) {
        rRigthHandSideVector.resize(LocalSize, false);
    }
    if (rLeftHandSideMatrix.size1() != LocalSize || rLeftHandSideMatrix.size2() != LocalSize) {
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
    }

    // Get the current edge element data
    const auto& r_geom = this->GetGeometry();
    const double h_e = r_geom.Length();
    const double kappa = rCurrentProcessInfo[GRADIENT_PENALTY_COEFFICIENT] * h_e;
    const double val_jump = r_geom[0].GetValue(NODAL_MAUX) - r_geom[1].GetValue(NODAL_MAUX);
    const array_1d<double,3> l_e = (r_geom[1].Coordinates() - r_geom[0].Coordinates()) / h_e;

    // Assemble the current edge RHS contribution
    const double aux = 2.0/h_e;
    std::array<double, 2> aux_stab{{1.0, -1.0}};
    noalias(rRigthHandSideVector) = ZeroVector(LocalSize);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    for (std::size_t i = 0; i < NumNodes; ++i) {
        const auto& r_grad_i = r_geom[i].FastGetSolutionStepValue(NODAL_VAUX);
        for (std::size_t d1 = 0; d1 < TDim; ++d1) {
            rRigthHandSideVector(i*BlockSize + d1) -= aux * l_e[d1] * val_jump;
            for (std::size_t j = 0; j < NumNodes; ++j) {
                const auto& r_grad_j = r_geom[j].FastGetSolutionStepValue(NODAL_VAUX);
                rLeftHandSideMatrix(i*BlockSize + d1, j*BlockSize + d1) += kappa * aux_stab[i] * aux_stab[j];
                rRigthHandSideVector(i*BlockSize + d1) -= kappa * (aux_stab[i] * r_grad_i[d1] - aux_stab[j] * r_grad_j[d1]);
                for (std::size_t d2 = 0; d2 < TDim; ++d2) {
                    rLeftHandSideMatrix(i*BlockSize + d1, j*BlockSize + d2) += l_e[d1]*l_e[d2];
                    rRigthHandSideVector(i*BlockSize + d1) -= l_e[d1] * l_e[d2] * r_grad_j[d2];
                }
            }
        }
    }
}

template<std::size_t TDim>
int EdgeBasedGradientRecoveryElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Perform basic element checks
    int error_code = Kratos::Element::Check(rCurrentProcessInfo);
    if(error_code != 0) {
        return error_code;
    }

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for(auto &r_node : this->GetGeometry()) {
        KRATOS_ERROR_IF(!r_node.SolutionStepsDataHas(NODAL_VAUX)) << "Missing NODAL_VAUX variable on solution step data for node " << r_node.Id() << std::endl;;
    }

    return 0;

    KRATOS_CATCH("");
}

template<std::size_t TDim>
void EdgeBasedGradientRecoveryElement<TDim>::CalculateLeftHandSide(
    MatrixType &rLeftHandSideMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Check size
    if (rLeftHandSideMatrix.size1() != LocalSize || rLeftHandSideMatrix.size2() != LocalSize) {
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
    }

    // Get the current edge element data
    const auto& r_geom = this->GetGeometry();
    const double h_e = r_geom.Length();
    const double kappa = rCurrentProcessInfo[GRADIENT_PENALTY_COEFFICIENT] * h_e;
    const array_1d<double,3> l_e = (r_geom[0].Coordinates() - r_geom[1].Coordinates()) / h_e;

    // Assemble the current edge LHS
    std::array<double, 2> aux_stab{{1.0, -1.0}};
    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    for (std::size_t i = 0; i < NumNodes; ++i) {
        for (std::size_t d1 = 0; d1 < TDim; ++d1) {
            for (std::size_t j = 0; j < NumNodes; ++j) {
                rLeftHandSideMatrix(i*BlockSize + d1, j*BlockSize + d1) += kappa * aux_stab[i] * aux_stab[j];
                for (std::size_t d2 = 0; d2 < TDim; ++d2) {
                    rLeftHandSideMatrix(i*BlockSize + d1, j*BlockSize + d2) += l_e[d1]*l_e[d2];
                }
            }
        }
    }
}

template<std::size_t TDim>
void EdgeBasedGradientRecoveryElement<TDim>::CalculateRightHandSide(
    VectorType &rRigthHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Check size
    if (rRigthHandSideVector.size() != LocalSize) {
        rRigthHandSideVector.resize(LocalSize, false);
    }

    // Get the current edge element data
    const auto& r_geom = this->GetGeometry();
    const double h_e = r_geom.Length();
    const double kappa = rCurrentProcessInfo[GRADIENT_PENALTY_COEFFICIENT] * h_e;
    const double val_jump = r_geom[0].GetValue(NODAL_MAUX) - r_geom[1].GetValue(NODAL_MAUX);
    const array_1d<double,3> l_e = (r_geom[0].Coordinates() - r_geom[1].Coordinates()) / h_e;

    // Assemble the current edge RHS contribution
    const double aux = 2.0/h_e;
    std::array<double, 2> aux_stab{{1.0, -1.0}};
    noalias(rRigthHandSideVector) = ZeroVector(LocalSize);
    for (std::size_t i = 0; i < NumNodes; ++i) {
        const auto& r_grad_i = r_geom[i].FastGetSolutionStepValue(NODAL_VAUX);
        for (std::size_t d1 = 0; d1 < TDim; ++d1) {
            rRigthHandSideVector(i*BlockSize + d1) -= aux * l_e[d1] * val_jump;
            for (std::size_t d2 = 0; d2 < TDim; ++d2) {
                rRigthHandSideVector(i*BlockSize + d1) -= l_e[d1] * l_e[d2] * r_grad_i[d2];
            }
            for (std::size_t j = 0; j < NumNodes; ++j) {
                const auto& r_grad_j = r_geom[j].FastGetSolutionStepValue(NODAL_VAUX);
                rRigthHandSideVector(i*BlockSize + d1) -= kappa * (aux_stab[i] * r_grad_i[d1] - aux_stab[j] * r_grad_j[d1]);
            }
        }
    }
}

template<>
void EdgeBasedGradientRecoveryElement<2>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rResult.size() != LocalSize) {
        rResult.resize(LocalSize, false);
    }

    const auto& r_geom = GetGeometry();
    const unsigned int x_pos = r_geom[0].GetDofPosition(NODAL_VAUX_X);
    for (unsigned int i = 0; i < NumNodes; i++) {
        rResult[i * 2] = r_geom[i].GetDof(NODAL_VAUX_X, x_pos).EquationId();
        rResult[i * 2 + 1] = r_geom[i].GetDof(NODAL_VAUX_Y, x_pos + 1).EquationId();
    }
}

template<>
void EdgeBasedGradientRecoveryElement<3>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rResult.size() != LocalSize) {
        rResult.resize(LocalSize, false);
    }

    const auto& r_geom = GetGeometry();
    const unsigned int x_pos = r_geom[0].GetDofPosition(NODAL_VAUX_X);
    for (unsigned int i = 0; i < NumNodes; i++) {
        rResult[i * 3] = r_geom[i].GetDof(NODAL_VAUX_X, x_pos).EquationId();
        rResult[i * 3 + 1] = r_geom[i].GetDof(NODAL_VAUX_Y, x_pos + 1).EquationId();
        rResult[i * 3 + 2] = r_geom[i].GetDof(NODAL_VAUX_Z, x_pos + 2).EquationId();
    }
}

template<>
void EdgeBasedGradientRecoveryElement<2>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rElementalDofList.size() != LocalSize) {
        rElementalDofList.resize(LocalSize);
    }

    const auto& r_geom = GetGeometry();
    for (unsigned int i = 0; i < NumNodes; i++) {
        rElementalDofList[i * 2] = r_geom[i].pGetDof(NODAL_VAUX_X);
        rElementalDofList[i * 2 + 1] = r_geom[i].pGetDof(NODAL_VAUX_Y);
    }
}

template<>
void EdgeBasedGradientRecoveryElement<3>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rElementalDofList.size() != LocalSize) {
        rElementalDofList.resize(LocalSize);
    }

    const auto& r_geom = GetGeometry();
    for (unsigned int i = 0; i < NumNodes; i++) {
        rElementalDofList[i * 3] = r_geom[i].pGetDof(NODAL_VAUX_X);
        rElementalDofList[i * 3 + 1] = r_geom[i].pGetDof(NODAL_VAUX_Y);
        rElementalDofList[i * 3 + 2] = r_geom[i].pGetDof(NODAL_VAUX_Z);
    }
}

template class EdgeBasedGradientRecoveryElement<2>;
template class EdgeBasedGradientRecoveryElement<3>;

} // namespace Kratos.
