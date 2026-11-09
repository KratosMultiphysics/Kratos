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

class KRATOS_API(IGA_APPLICATION) GsbmStokesConditionDirichlet
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GsbmStokesConditionDirichlet);

    using SizeType = std::size_t;
    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Id and geometry.
    GsbmStokesConditionDirichlet(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {
    }

    /// Constructor with Id, geometry and properties.
    GsbmStokesConditionDirichlet(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {
    }

    /// Default constructor.
    GsbmStokesConditionDirichlet()
        : Condition()
    {
    }

    /// Destructor.
    ~GsbmStokesConditionDirichlet() override = default;

    /// Create with Id, pointer to geometry and pointer to properties
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GsbmStokesConditionDirichlet>(
            NewId, pGeom, pProperties
        );
    }

    /// Create with Id, nodes and pointer to properties
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GsbmStokesConditionDirichlet>(
            NewId, GetGeometry().Create(ThisNodes), pProperties
        );
    }

    ///@}
    ///@name Operations
    ///@{

    /// Clones the constitutive law of the properties
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /// Shifted Dirichlet condition of the GSBM Stokes problem on a cut facet,
    /// by Nitsche's method, in residual form: the consistency term
    /// -([[v_\alpha]], {{\sigma_h}} n-), the viscous adjoint term with
    /// \theta = ADJOINT_FACTOR (-1 skew-symmetric by default), the q term
    /// unless ADJOINT_ON_PRESSURE is false, and the penalty if PENALTY_FACTOR
    /// is positive. The stress comes from the constitutive law, the imposed
    /// velocity from the VELOCITY of the projection node, and the shift order
    /// from SHIFT_ORDER, the largest degree of the spaces by default
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

    /// Equation Ids of u_x, u_y and p of the - side, then of the + side
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /// Dofs of u_x, u_y and p of the - side, then of the + side
    void GetDofList(
        DofsVectorType& rConditionDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Check
    ///@{

    /// Checks the geometry parts, SHIFT_ORDER, the constitutive law and, with a
    /// penalty, DYNAMIC_VISCOSITY
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "GsbmStokesConditionDirichlet #" + std::to_string(Id());
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
    /// space on the - side
    void InitializeMaterial();

    /// Constitutive matrix D of the Voigt strain [xx, yy, xy], at the current
    /// strain of a facet side
    Matrix CalculateConstitutiveMatrix(
        const GeometryType& rSide,
        const ProcessInfo& rCurrentProcessInfo) const;

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
