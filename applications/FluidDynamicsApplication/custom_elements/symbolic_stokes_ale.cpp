//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Nicolò Antonelli
//

// System includes

// Project includes
#include "includes/checks.h"

// Application includes
#include "custom_elements/symbolic_stokes_ale.h"

namespace Kratos
{

Element::Pointer SymbolicStokesALE2D4N::Create(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicStokesALE2D4N>(NewId, pGeometry, pProperties);
}

Element::Pointer SymbolicStokesALE2D4N::Create(
    IndexType NewId,
    const NodesArrayType& rNodes,
    Properties::Pointer pProperties) const
{
    return Create(NewId, GetGeometry().Create(rNodes), pProperties);
}

int SymbolicStokesALE2D4N::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const int result = BaseType::Check(rCurrentProcessInfo);
    const auto& r_geometry = GetGeometry();
    KRATOS_ERROR_IF(r_geometry[0].GetBufferSize() < 3)
        << "ALE Stokes requires a three-step buffer." << std::endl;
    for (const auto& r_node : r_geometry) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_VELOCITY, r_node);
        KRATOS_ERROR_IF(std::abs(r_node.FastGetSolutionStepValue(MESH_VELOCITY)[2]) > 1e-12)
            << "ALE Stokes requires planar mesh velocity." << std::endl;
        KRATOS_ERROR_IF(norm_2(r_node.FastGetSolutionStepValue(MESH_VELOCITY) -
            r_geometry[0].FastGetSolutionStepValue(MESH_VELOCITY)) > 1e-12)
            << "SymbolicStokesALE2D4N currently supports rigid translation only." << std::endl;
    }
    return result;

    KRATOS_CATCH("")
}

void SymbolicStokesALE2D4N::AddTimeIntegratedSystem(
    DataType& rData,
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector)
{
    BaseType::AddTimeIntegratedSystem(rData, rLeftHandSideMatrix, rRightHandSideVector);
    AddMeshTransport(rData, &rLeftHandSideMatrix, &rRightHandSideVector);
}

void SymbolicStokesALE2D4N::AddTimeIntegratedLHS(
    DataType& rData,
    MatrixType& rLeftHandSideMatrix)
{
    BaseType::AddTimeIntegratedLHS(rData, rLeftHandSideMatrix);
    AddMeshTransport(rData, &rLeftHandSideMatrix, nullptr);
}

void SymbolicStokesALE2D4N::AddTimeIntegratedRHS(
    DataType& rData,
    VectorType& rRightHandSideVector)
{
    BaseType::AddTimeIntegratedRHS(rData, rRightHandSideVector);
    AddMeshTransport(rData, nullptr, &rRightHandSideVector);
}

void SymbolicStokesALE2D4N::AddMeshTransport(
    const DataType& rData,
    MatrixType* pLeftHandSideMatrix,
    VectorType* pRightHandSideVector) const
{
    array_1d<double, 3> mesh_velocity = ZeroVector(3);
    for (std::size_t i = 0; i < 4; ++i) {
        mesh_velocity += rData.N[i] * GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY);
    }
    const double tau = 1.0 / (rData.DynamicTau * rData.Density / rData.DeltaTime
        + 4.0 * rData.EffectiveViscosity / (rData.ElementSize * rData.ElementSize));

    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
            const double transport = -rData.Density * (mesh_velocity[0] * rData.DN_DX(j, 0)
                + mesh_velocity[1] * rData.DN_DX(j, 1)) * rData.Weight;
            for (std::size_t d = 0; d < 2; ++d) {
                const double galerkin = rData.N[i] * transport;
                const double pspg = tau * rData.DN_DX(i, d) * transport;
                if (pLeftHandSideMatrix) {
                    (*pLeftHandSideMatrix)(3*i+d, 3*j+d) += galerkin;
                    (*pLeftHandSideMatrix)(3*i+2, 3*j+d) += pspg;
                }
                if (pRightHandSideVector) {
                    (*pRightHandSideVector)[3*i+d] -= galerkin * rData.Velocity(j, d);
                    (*pRightHandSideVector)[3*i+2] -= pspg * rData.Velocity(j, d);
                }
            }
        }
    }
}

void SymbolicStokesALE2D4N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

void SymbolicStokesALE2D4N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

} // namespace Kratos
