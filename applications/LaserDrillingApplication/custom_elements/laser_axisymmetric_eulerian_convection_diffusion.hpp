//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

#pragma once

// System includes


// External includes


// Project includes
#include "geometries/geometry_data.h"
#include "includes/define.h"
#include "includes/serializer.h"

// Application includes
#include "../ConvectionDiffusionApplication/custom_elements/axisymmetric_eulerian_convection_diffusion.h"

///@addtogroup LaserDrillingApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(LASER_DRILLING_APPLICATION) LaserAxisymmetricEulerianConvectionDiffusionElement
    : public AxisymmetricEulerianConvectionDiffusionElement<TDim, TNumNodes>
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of LaserAxisymmetricEulerianConvectionDiffusionElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LaserAxisymmetricEulerianConvectionDiffusionElement);

    /// Base convection-diffusion element type
    using BaseType = AxisymmetricEulerianConvectionDiffusionElement<TDim, TNumNodes>;

    /// Geometry type
    using GeometryType = typename BaseType::GeometryType;

    /// Properties type
    using PropertiesType = typename BaseType::PropertiesType;

    /// Index type
    using IndexType = typename BaseType::IndexType;

    /// Size type
    using SizeType = typename BaseType::SizeType;

    /// Vector type
    using VectorType = typename BaseType::VectorType;

    /// Matrix type
    using MatrixType = typename BaseType::MatrixType;

    /// Nodes array type
    using NodesArrayType = typename BaseType::NodesArrayType;

    /// Shape functions gradient container type
    using ShapeFunctionsGradientsType = typename GeometryType::ShapeFunctionsGradientsType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructors.
    LaserAxisymmetricEulerianConvectionDiffusionElement() : BaseType()
    {
    }

    LaserAxisymmetricEulerianConvectionDiffusionElement(
        IndexType NewId,
        typename GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {}

    LaserAxisymmetricEulerianConvectionDiffusionElement(
        IndexType NewId,
        typename GeometryType::Pointer pGeometry,
        typename PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~LaserAxisymmetricEulerianConvectionDiffusionElement(){};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<LaserAxisymmetricEulerianConvectionDiffusionElement<TDim, TNumNodes>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(
        IndexType NewId,
        typename GeometryType::Pointer pGeom,
        typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<LaserAxisymmetricEulerianConvectionDiffusionElement<TDim, TNumNodes>>(NewId, pGeom, pProperties);
    }


    ///@}
    ///@name Access
    ///@{

    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;


    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override
    {
        return "AxisymmetricConvectionDiffusion #";
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << this->Info() << this->Id();
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

    static constexpr GeometryData::IntegrationMethod mMyIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}
private:
    ///@name Static Member Variables
    ///@{

    /// Integration rule to be employed

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Friends
    ///@{


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
    }

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}
}; // Class AxisymmetricConvectionDiffusion

///@}

///@} // LaserDrillingApplication group

} // namespace Kratos.
