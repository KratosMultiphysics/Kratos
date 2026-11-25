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

class KRATOS_API(IGA_APPLICATION) GsbmStokesConditionDirichletMultiplier
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GsbmStokesConditionDirichletMultiplier);

    using SizeType = std::size_t;
    using IndexType = std::size_t;

    /// Part of the geometry with the multiplier curve, after both facet sides
    static constexpr IndexType MultiplierPart = 2;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Id and geometry.
    GsbmStokesConditionDirichletMultiplier(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {
    }

    /// Constructor with Id, geometry and properties.
    GsbmStokesConditionDirichletMultiplier(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {
    }

    /// Default constructor.
    GsbmStokesConditionDirichletMultiplier()
        : Condition()
    {
    }

    /// Destructor.
    ~GsbmStokesConditionDirichletMultiplier() override = default;

    /// Create with Id, pointer to geometry and pointer to properties
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GsbmStokesConditionDirichletMultiplier>(
            NewId, pGeom, pProperties
        );
    }

    /// Create with Id, nodes and pointer to properties
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GsbmStokesConditionDirichletMultiplier>(
            NewId, GetGeometry().Create(ThisNodes), pProperties
        );
    }

    ///@}
    ///@name Operations
    ///@{

    /// Clones the constitutive law of the properties
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /// Shifted Dirichlet condition of the GSBM Stokes problem on a cut facet,
    /// in residual form: the normal velocity by Lagrange multipliers \lambda,
    /// the pressure on the skin, and the tangential one by Nitsche's method.
    /// The viscous consistency -([[v_\alpha]], {{2 \mu \epsilon(u_h)}} n-),
    /// the pressure +([[v_\alpha]] . n-, \lambda - ({{S p_h}} - {{p_h}})), the
    /// tangential adjoint term with \theta = ADJOINT_FACTOR (-1 skew-symmetric
    /// by default) and penalty if PENALTY_FACTOR is positive, along the skin's
    /// tangent, and along n- too with PENALTY_ON_NORMAL, as an augmented
    /// Lagrangian. The multiplier rows are mean-free,
    /// ([[\alpha_h]] (\mu_k - c_k), n- . ({{S(u_h)}} - g)) = 0, with
    /// c_k = m_k / \sum_j m_j and m_k the integral of \mu_k over the skin
    /// parameter, so they sum to zero and the last follows from the others.
    /// It is the pin (|[[\alpha_h]]|, \lambda - {{S p_h}}) = 0 instead, which
    /// ties \lambda to the pressure, in place of the rows' constant
    /// combination, which repeats the continuity equations of a
    /// divergence-free velocity. With PIN_MULTIPLIERS false, for a velocity
    /// that is not, the rows are
    /// ([[\alpha_h]] \mu_k, n- . ({{S(u_h)}} - g)) = 0 for every k, without
    /// the pin. \lambda is read at the closest point's normalized skin
    /// parameter, INTEGRATION_COORDINATES[0] of the projection node, which
    /// also gives VELOCITY g and LOCAL_TANGENT
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

    /// Equation Ids of u_x, u_y and p of the - side, then of the + side, then
    /// of every multiplier
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /// Dofs of u_x, u_y and p of the - side, then of the + side, then of
    /// every multiplier
    void GetDofList(
        DofsVectorType& rConditionDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Check
    ///@{

    /// Checks the geometry parts, the multiplier curve, SHIFT_ORDER, the
    /// constitutive law and, with a penalty, DYNAMIC_VISCOSITY
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "GsbmStokesConditionDirichletMultiplier #" + std::to_string(Id());
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

    /// Number of multipliers: the control points of the periodic curve
    /// without the ones that wrap around
    SizeType NumberOfMultipliers() const;

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
