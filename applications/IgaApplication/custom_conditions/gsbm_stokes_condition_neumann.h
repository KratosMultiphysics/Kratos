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

class KRATOS_API(IGA_APPLICATION) GsbmStokesConditionNeumann
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GsbmStokesConditionNeumann);

    using SizeType = std::size_t;
    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Id and geometry.
    GsbmStokesConditionNeumann(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {
    }

    /// Constructor with Id, geometry and properties.
    GsbmStokesConditionNeumann(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {
    }

    /// Default constructor.
    GsbmStokesConditionNeumann()
        : Condition()
    {
    }

    /// Destructor.
    ~GsbmStokesConditionNeumann() override = default;

    /// Create with Id, pointer to geometry and pointer to properties
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GsbmStokesConditionNeumann>(
            NewId, pGeom, pProperties
        );
    }

    /// Create with Id, nodes and pointer to properties
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GsbmStokesConditionNeumann>(
            NewId, GetGeometry().Create(ThisNodes), pProperties
        );
    }

    ///@}
    ///@name Operations
    ///@{

    /// Clones the constitutive law of the properties
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /// Shifted Neumann condition of the GSBM Stokes problem on a cut facet,
    /// in residual form: the consistency term -([[v_\alpha]], {{\sigma_h}} n-)
    /// and the shifted traction ([[v_\alpha]], (n- . n) ({{S_d^m[\sigma_h]}} n
    /// - t_N)), with the stress of the constitutive law, the prescribed
    /// traction t_N = -p_D n + NORMAL_STRESS from the PRESSURE and the
    /// optional NORMAL_STRESS of the projection node, and m = SHIFT_ORDER, the
    /// largest degree of the spaces by default
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

    /// Checks the geometry parts and the constitutive law
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "GsbmStokesConditionNeumann #" + std::to_string(Id());
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
