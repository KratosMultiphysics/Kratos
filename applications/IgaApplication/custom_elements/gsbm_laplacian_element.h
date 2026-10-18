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
#include "custom_elements/laplacian_element.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(IGA_APPLICATION) GsbmLaplacianElement
    : public LaplacianElement
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GsbmLaplacianElement);

    using BaseType = LaplacianElement;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Id and geometry.
    GsbmLaplacianElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : LaplacianElement(NewId, pGeometry)
    {
    }

    /// Constructor with Id, geometry and properties.
    GsbmLaplacianElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : LaplacianElement(NewId, pGeometry, pProperties)
    {
    }

    /// Default constructor.
    GsbmLaplacianElement()
        : LaplacianElement()
    {
    }

    /// Destructor.
    ~GsbmLaplacianElement() override = default;

    ///@}
    ///@name Operations
    ///@{

    /// Create with Id, nodes and pointer to properties
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GsbmLaplacianElement>(
            NewId, GetGeometry().Create(ThisNodes), pProperties
        );
    }

    /// Create with Id, pointer to geometry and pointer to properties
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GsbmLaplacianElement>(NewId, pGeom, pProperties);
    }

    /// Left hand side of LaplacianElement, weighted by the volume fraction
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /// Right hand side of LaplacianElement, weighted by the volume fraction
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /// Mass matrix of LaplacianElement, weighted by the volume fraction
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /// Checks LaplacianElement's requirements and that VOLUME_FRACTION is set
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "GsbmLaplacianElement #" + std::to_string(Id());
    }

    ///@}

private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, LaplacianElement);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, LaplacianElement);
    }

    ///@}
};

///@}

} // namespace Kratos
