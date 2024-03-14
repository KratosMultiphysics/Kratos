//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alessandro Franci
//                   Ruben Zorrilla
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "includes/constitutive_law.h"
#include "includes/define.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "modified_shape_functions/modified_shape_functions.h"
#include "utilities/math_utils.h"

// Application includes
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_fluid_FIC_element.h"

namespace Kratos
{

///@addtogroup PfemFluidDynamicsApplication
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

/**
 * @brief A FIC-stabilized CutPFEM element for the weakly-compressible Navier-Stokes equations
 *
 * @tparam TDim Number of dimensions
 */
template <unsigned int TDim>
class TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement : public TwoStepUpdatedLagrangianVPImplicitFluidFicElement<TDim>
{

  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement);

    ///base type:
    using BaseType = TwoStepUpdatedLagrangianVPImplicitFluidFicElement<TDim>;

    /// Node type
    using NodeType = Node;

    /// Geometry type (using with given NodeType)
    using GeometryType = Geometry<NodeType>;

    /// Definition of nodes container type, redefined from GeometryType
    using NodesArrayType = Geometry<NodeType>::PointsArrayType;

    /// Vector type for local contributions to the linear system
    using VectorType = Vector;

    /// Matrix type for local contributions to the linear system
    using MatrixType = Matrix;

    using IndexType = ::size_t;

    using SizeType = std::size_t;

    using EquationIdVectorType = std::vector<std::size_t>;

    using DofsVectorType = std::vector<Dof<double>::Pointer>;

    using DofsArrayType = PointerVectorSet<Dof<double>, IndexedObject>;

    /// Type for shape function values container
    using ShapeFunctionsType = Kratos::Vector;

    /// Type for a matrix containing the shape function gradients
    using ShapeFunctionDerivativesType = Kratos::Matrix;

    /// Type for an array of shape function gradient matrices
    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    using PropertiesType = typename BaseType::PropertiesType;

    using pPropertiesType = typename BaseType::PropertiesType::Pointer;

    using ElementalVariables = typename BaseType::ElementalVariables;

    using NodeWeakPtrVectorType = GlobalPointersVector<NodeType>;

    /// Reference type definition for constitutive laws
    using ConstitutiveLawType = ConstitutiveLaw;

    ///Pointer type for constitutive laws
    using ConstitutiveLawPointerType = ConstitutiveLawType::Pointer;

    /// Number of nodes
    static constexpr SizeType NumNodes = TDim + 1;

    /// Voigt size
    static constexpr SizeType StrainSize = TDim == 2 ? 3 : 6;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constuctor.
    TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement(IndexType NewId = 0)
      : BaseType(NewId)
    {
    }

    /// Constructor using an array of nodes.
    TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement(
        IndexType NewId,
        const NodesArrayType &ThisNodes)
        : BaseType(NewId, ThisNodes)
    {
    }

    /// Constructor using a geometry object.
    TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    /// Constuctor using geometry and properties.
    TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        pPropertiesType pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /// copy constructor
    TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement(TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement const &rOther) : BaseType(rOther)
    {
    }

    /// Destructor
    virtual ~TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement &operator=(TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement const &rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const &ThisNodes,
        pPropertiesType pProperties) const override
    {
      return Kratos::make_intrusive<TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement>(NewId, BaseType::GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Clone(
        IndexType NewId,
        NodesArrayType const &ThisNodes) const override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Elemental Data
    ///@{

    int Check(const ProcessInfo &rCurrentProcessInfo) const override;

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
      std::stringstream buffer;
      buffer << "TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement #" << BaseType::Id();
      return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
      rOStream << "TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement" << TDim << "D";
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

    void CalculateLocalMomentumEquations(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo) override;

    void CalculateGeometryData(
        ShapeFunctionDerivativesArrayType &rDN_DX,
        Matrix &rNContainer,
        Vector &rGaussWeights) override;

    void CalculateGeometryData(Vector &rGaussWeights) override;

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

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer &rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer &rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    bool IsCut() const;

    bool IsPositive() const;

    void CalculateCutGeometryData(
        ShapeFunctionDerivativesArrayType &rDNDX,
        Matrix &rN,
        Vector &rGaussWeights);

    void CalculateIntersectionGeometryData(
        ShapeFunctionDerivativesArrayType &rInterfaceDNDX,
        Matrix &rInterfaceN,
        Vector &rInterfaceGaussWeights,
        ModifiedShapeFunctions::AreaNormalsContainerType& rInterfaceUnitNormals);

    void CalculateCutGeometryData(Vector &rGaussWeights);

    void VoigtStressNormalProjection(
      const Vector& rVoigtStress,
      const array_1d<double,3>& rUnitNormal,
      array_1d<double,TDim>& rProjectedStress);

    void CalculateBMatrix(
        const Matrix& rDNDX,
        BoundedMatrix<double,StrainSize, TDim*NumNodes>& rB);

    void VoigtTransformForProduct(
        const array_1d<double,3>& rVector,
        BoundedMatrix<double, TDim, StrainSize>& rVoigtMatrix);

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
}; // Class TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement
///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TDim>
inline std::istream &operator>>(
    std::istream &rIStream,
    TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim> &rThis)
{
    return rIStream;
}

/// output stream function
template <unsigned int TDim>
inline std::ostream &operator<<(
    std::ostream &rOStream,
    const TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim> &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.
