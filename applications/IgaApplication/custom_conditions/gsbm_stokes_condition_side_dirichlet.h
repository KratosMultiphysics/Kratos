//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mario Caballero
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/condition.h"
#include "includes/constitutive_law.h"
#include "includes/define.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(IGA_APPLICATION) GsbmStokesConditionSideDirichlet
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GsbmStokesConditionSideDirichlet);

    using SizeType = std::size_t;
    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Id and geometry.
    GsbmStokesConditionSideDirichlet(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {
    }

    /// Constructor with Id, geometry and properties.
    GsbmStokesConditionSideDirichlet(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {
    }

    /// Default constructor.
    GsbmStokesConditionSideDirichlet()
        : Condition()
    {
    }

    /// Destructor.
    ~GsbmStokesConditionSideDirichlet() override = default;

    /// Create with Id, pointer to geometry and pointer to properties
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GsbmStokesConditionSideDirichlet>(
            NewId, pGeom, pProperties
        );
    }

    /// Create with Id, nodes and pointer to properties
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GsbmStokesConditionSideDirichlet>(
            NewId, GetGeometry().Create(ThisNodes), pProperties
        );
    }

    ///@}
    ///@name Operations
    ///@{

    /// Clones the constitutive law of the properties
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /// Dirichlet condition on a side of the patches by Nitsche's method, the
    /// GSBM one on a fitted side (alpha = 1, no shift), in residual form: the
    /// consistency term -(v_h, \sigma_h n), the viscous adjoint term with
    /// \theta = ADJOINT_FACTOR (-1 skew-symmetric by default), the q term
    /// unless ADJOINT_ON_PRESSURE is false, and the penalty if PENALTY_FACTOR
    /// is positive, with n the outward normal, not weighted by the volume
    /// fraction. The stress comes from the constitutive law and the imposed
    /// velocity from the VELOCITY of the condition
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /// Left hand side of CalculateLocalSystem
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /// Right hand side of CalculateLocalSystem
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /// Equation Ids of u_x, u_y and p
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /// Dofs of u_x, u_y and p
    void GetDofList(
        DofsVectorType& rConditionDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Check
    ///@{

    /// Checks the geometry parts, the constitutive law and, with the penalty,
    /// DYNAMIC_VISCOSITY
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "GsbmStokesConditionSideDirichlet #" + std::to_string(Id());
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    ConstitutiveLaw::Pointer mpConstitutiveLaw = nullptr;

    ///@}
    ///@name Private Operations
    ///@{

    /// Clones the constitutive law of the properties, at the point of u_x's
    /// space
    void InitializeMaterial();

    /// Constitutive matrix D of the Voigt strain [xx, yy, xy], at the current
    /// strain of the side point
    Matrix CalculateConstitutiveMatrix(const ProcessInfo& rCurrentProcessInfo) const;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
        rSerializer.save("ConstitutiveLaw", mpConstitutiveLaw);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
        rSerializer.load("ConstitutiveLaw", mpConstitutiveLaw);
    }

    ///@}
};

///@}

} // namespace Kratos
