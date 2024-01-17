// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes


// External includes


// Project includes
#include "geometries/geometry_data.h"
#include "includes/define.h"
#include "includes/serializer.h"

// Application includes
#include "eulerian_conv_diff.h"

///@addtogroup FluidDynamicsApplication
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
class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) AxisymmetricEulerianConvectionDiffusionElement
    : public EulerianConvectionDiffusionElement<TDim, TNumNodes>
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of AxisymmetricEulerianConvectionDiffusionElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AxisymmetricEulerianConvectionDiffusionElement);

    /// Base convection-diffusion element type
    using BaseType = EulerianConvectionDiffusionElement<TDim, TNumNodes>;

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
    AxisymmetricEulerianConvectionDiffusionElement() : BaseType()
    {
    }

    AxisymmetricEulerianConvectionDiffusionElement(
        IndexType NewId,
        typename GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {}

    AxisymmetricEulerianConvectionDiffusionElement(
        IndexType NewId,
        typename GeometryType::Pointer pGeometry,
        typename PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~AxisymmetricEulerianConvectionDiffusionElement(){};

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
        return Kratos::make_intrusive<AxisymmetricEulerianConvectionDiffusionElement<TDim, TNumNodes>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(
        IndexType NewId,
        typename GeometryType::Pointer pGeom,
        typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AxisymmetricEulerianConvectionDiffusionElement<TDim, TNumNodes>>(NewId, pGeom, pProperties);
    }

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

	void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{

    int Check(const ProcessInfo &rCurrentProcessInfo) const override;

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


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    double CalculateTau(
        const typename BaseType::ElementVariables &rVariables,
        const double NormVelocity,
        const double ElementSize);

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
    static constexpr GeometryData::IntegrationMethod mIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;

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

    /**
     * @brief Calculates values at current Gauss point
     * This method calculates some values at current Gauss point
     * Note that the evaluation of these terms is already done at n+\theta
     * @param rN Current Gauss point shape functions values
     * @param rDNDX Current Gauss point shape functions gradients values
     * @param rVariables Elemental data container
     * @param rRadius Current Gauss point radius
     * @param rVelocity Current Gauss point velocity at n+\theta
     * @param rConvectionOperator Current Gauss point evaluation of the convection operator at n+\theta
     * @param rVelocityGradient Current Gauss point velocity gradient
     */
    void CalculateGaussPointData(
        const array_1d<double, TNumNodes> &rN,
        const BoundedMatrix<double, TNumNodes, TDim> &rDNDX,
        typename BaseType::ElementVariables &rVariables,
        double &rRadius,
        array_1d<double, TDim> &rVelocity,
        array_1d<double, TNumNodes> &rConvectionOperator,
        BoundedMatrix<double, TDim, TDim> &rVelocityGradient) const;

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

///@} // ConvectionDiffusionApplication group

} // namespace Kratos.
