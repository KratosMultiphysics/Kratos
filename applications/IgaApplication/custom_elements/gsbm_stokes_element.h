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
#include "includes/element.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(IGA_APPLICATION) GsbmStokesElement
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GsbmStokesElement);

    using BaseType = Element;

    /// Geometry parts of the spaces of the fields
    static constexpr IndexType VelocityXPart = 0;
    static constexpr IndexType VelocityYPart = 1;
    static constexpr IndexType PressurePart = 2;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Id and geometry.
    GsbmStokesElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    /// Constructor with Id, geometry and properties.
    GsbmStokesElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /// Default constructor.
    GsbmStokesElement()
        : BaseType()
    {
    }

    /// Destructor.
    ~GsbmStokesElement() override = default;

    ///@}
    ///@name Operations
    ///@{

    /// Create with Id, nodes and pointer to properties
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GsbmStokesElement>(
            NewId, GetGeometry().Create(ThisNodes), pProperties
        );
    }

    /// Create with Id, pointer to geometry and pointer to properties
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GsbmStokesElement>(NewId, pGeom, pProperties);
    }

    /// Clones the constitutive law of the properties
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /// Galerkin Stokes terms at the Gauss point in residual form, with one
    /// geometry part per field (u_x, u_y and p), and PSPG if PSPG_FACTOR is
    /// given, weighted by the volume fraction; and grad-div if GRAD_DIV_FACTOR
    /// is given, not weighted. PSPG reads second derivatives, so it needs
    /// shape_function_derivatives_order = 3 in the modeler
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

    /// Equation Ids of the u_x space's control points, then u_y's, then p's
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /// Dofs of the u_x space's control points, then u_y's, then p's
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Check
    ///@{

    /// Checks the geometry parts, VOLUME_FRACTION, the constitutive law, the
    /// dofs and, with PSPG_FACTOR, DYNAMIC_VISCOSITY and KNOT_SPAN_SIZES
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "GsbmStokesElement #" + std::to_string(Id());
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

    /// Current values of the velocity dofs: u_x's control points, then u_y's
    void GetSolutionCoefficientVector(Vector& rValues) const;

    /// Voigt strain matrix [xx, yy, xy] of the velocity dofs, u_x's then u_y's
    void CalculateB(Matrix& rB) const;

    /// Stress and D of the constitutive law at the current strain B u
    void ApplyConstitutiveLaw(
        const Matrix& rB,
        ConstitutiveLaw::Parameters& rValues,
        Vector& rStressVector,
        Matrix& rD) const;

    /// Viscous term (D B u, B v) and body force (f, v), with their residual
    void AddMomentumTerms(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const array_1d<double, 3>& rBodyForce,
        const Matrix& rB,
        const Matrix& rD,
        const Vector& rStressVector,
        const double Weight) const;

    /// Pressure terms -(p, div v) and +(q, div u), with their residual
    void AddContinuityTerms(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const Matrix& rB,
        const double Weight) const;

    /// PSPG (grad q, tau (grad p - div sigma - f)), with its residual
    void AddPspgTerms(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const array_1d<double, 3>& rBodyForce,
        const Matrix& rD,
        const double Factor,
        const double Weight) const;

    /// Grad-div (beta_div div v, div u), with its residual
    void AddGradDivTerms(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const Matrix& rB,
        const double Factor,
        const double Weight) const;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
        rSerializer.save("ConstitutiveLaw", mpConstitutiveLaw);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
        rSerializer.load("ConstitutiveLaw", mpConstitutiveLaw);
    }

    ///@}
};

///@}

} // namespace Kratos
