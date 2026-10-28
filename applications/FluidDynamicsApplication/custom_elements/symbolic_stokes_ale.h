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

#pragma once

// Application includes
#include "custom_elements/symbolic_stokes.h"
#include "custom_elements/data_containers/stokes/stokes_data.h"

namespace Kratos
{

///@addtogroup FluidDynamicsApplication
///@{

/**
 * @brief Transient Stokes element on a rigidly translating Q1 ALE mesh.
 * @details Retains the BDF mass, constitutive operator and ASGS scales of
 * SymbolicStokes. Adds the mesh transport -rho (w.grad)u to the Galerkin
 * momentum equations and the PSPG residual. Physical convection, SUPG and
 * pressure penalty terms are not added. Background and gap elements use
 * the ordinary SymbolicStokes element.
 */
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) SymbolicStokesALE2D4N
    : public SymbolicStokes<SymbolicStokesData<2, 4>>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SymbolicStokesALE2D4N);

    using BaseType = SymbolicStokes<SymbolicStokesData<2, 4>>;
    using DataType = SymbolicStokesData<2, 4>;

    ///@}
    ///@name Life Cycle
    ///@{

    using BaseType::BaseType;

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        Properties::Pointer pProperties) const override;

    Element::Pointer Create(
        IndexType NewId,
        const NodesArrayType& rNodes,
        Properties::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    void AddTimeIntegratedSystem(
        DataType& rData,
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector) override;

    void AddTimeIntegratedLHS(
        DataType& rData,
        MatrixType& rLeftHandSideMatrix) override;

    void AddTimeIntegratedRHS(
        DataType& rData,
        VectorType& rRightHandSideVector) override;

    ///@}

private:
    ///@name Private Operations
    ///@{

    /** @brief Add the ALE mesh-transport and corresponding PSPG terms. */
    void AddMeshTransport(
        const DataType& rData,
        MatrixType* pLeftHandSideMatrix,
        VectorType* pRightHandSideVector) const;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

    ///@}
};

///@}

} // namespace Kratos
