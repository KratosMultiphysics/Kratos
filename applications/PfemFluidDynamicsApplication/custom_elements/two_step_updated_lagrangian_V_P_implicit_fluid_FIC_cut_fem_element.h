//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 0.0 $
//
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
/* #include "includes/element.h" */
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "modified_shape_functions/modified_shape_functions.h"

#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_fluid_FIC_element.h"

namespace Kratos
{

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

  /// A stabilized element for the incompressible Navier-Stokes equations.
  /**
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
    using GometryType = Geometry<NodeType>;

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

    //Constructors.

    /// Default constuctor.
    /**
       * @param NewId Index number of the new element (optional)
       */
    TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement(IndexType NewId = 0) : BaseType(NewId)
    {
    }

    /// Constructor using an array of nodes.
    /**
       * @param NewId Index of the new element
       * @param ThisNodes An array containing the nodes of the new element
       */
    TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement(IndexType NewId, const NodesArrayType &ThisNodes) : BaseType(NewId, ThisNodes)
    {
    }

    /// Constructor using a geometry object.
    /**
       * @param NewId Index of the new element
       * @param pGeometry Pointer to a geometry object
       */
    TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry)
    {
    }

    /// Constuctor using geometry and properties.
    /**
       * @param NewId Index of the new element
       * @param pGeometry Pointer to a geometry object
       * @param pProperties Pointer to the element's properties
       */
    TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement(IndexType NewId, GeometryType::Pointer pGeometry, pPropertiesType pProperties) : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /// copy constructor

    TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement(TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement const &rOther) : BaseType(rOther)
    {
    }

    /// Destructor.
    virtual ~TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Create a new element of this type
    /**
       * Returns a pointer to a new TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement element, created using given input
       * @param NewId: the ID of the new element
       * @param ThisNodes: the nodes of the new element
       * @param pProperties: the properties assigned to the new element
       * @return a Pointer to the new element
       */
    Element::Pointer Create(IndexType NewId, NodesArrayType const &ThisNodes,
                            pPropertiesType pProperties) const override
    {
      return Element::Pointer(new TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement(NewId, BaseType::GetGeometry().Create(ThisNodes), pProperties));
    }

    Element::Pointer Clone(IndexType NewId, NodesArrayType const &ThisNodes) const override;

    // void Initialize(const ProcessInfo &rCurrentProcessInfo) override;

    // /// Initializes the element and all geometric information required for the problem.
    // void InitializeSolutionStep(const ProcessInfo &rCurrentProcessInfo) override{};

    // void InitializeNonLinearIteration(const ProcessInfo &rCurrentProcessInfo) override{};

    // void CalculateLeftHandSide(MatrixType &rLeftHandSideMatrix,
    //                            const ProcessInfo &rCurrentProcessInfo) override
    // {
    //   KRATOS_TRY;
    //   KRATOS_THROW_ERROR(std::logic_error, "TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement::CalculateLeftHandSide not implemented", "");
    //   KRATOS_CATCH("");
    // }

    // void CalculateRightHandSide(VectorType &rRightHandSideVector,
    //                             const ProcessInfo &rCurrentProcessInfo) override
    // {
    //   KRATOS_TRY;
    //   KRATOS_THROW_ERROR(std::logic_error, "TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement::CalculateRightHandSide not implemented", "");
    //   KRATOS_CATCH("");
    // }

    // The following methods have different implementations depending on TDim
    /// Provides the global indices for each one of this element's local rows
    /**
       * this determines the elemental equation ID vector for all elemental
       * DOFs
       * @param rResult A vector containing the global Id of each row
       * @param rCurrentProcessInfo the current process info object (unused)
       */

    /// Returns a list of the element's Dofs
    /**
       * @param ElementalDofList the list of DOFs
       * @param rCurrentProcessInfo the current process info instance
       */

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Elemental Data
    ///@{

    /// Checks the input and that all required Kratos variables have been registered.
    /**
       * This function provides the place to perform checks on the completeness of the input.
       * It is designed to be called only once (or anyway, not often) typically at the beginning
       * of the calculations, so to verify that nothing is missing from the input
       * or that no common error is found.
       * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
       * @return 0 if no errors were found.
       */
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

    //        /// Print object's data.
    //        virtual void PrintData(std::ostream& rOStream) const;

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

    ConstitutiveLaw::Pointer mpConstitutiveLaw = nullptr;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    // /**
    //    * A constistent mass matrix is used.
    //    * @param rMassMatrix The local matrix where the result will be added.
    //    * @param rN Elemental shape functions.
    //    * @param Weight Multiplication coefficient for the matrix, typically Density times integration point weight.
    //    */

    void CalculateLocalMomentumEquations(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo) override;

    // void ComputeBoundLHSMatrix(MatrixType &BoundLHSMatrix,
    //                            const ShapeFunctionsType &rN,
    //                            const double Weight) override;

    // void ComputeBoundRHSVector(VectorType &BoundRHSVector,
    //                            const ShapeFunctionsType &rN,
    //                            const double TimeStep,
    //                            const double BoundRHSCoeffAcc,
    //                            const double BoundRHSCoeffDev) override;

    // void ComputeBoundRHSVectorComplete(VectorType &BoundRHSVector,
    //                                    const double TimeStep,
    //                                    const double BoundRHSCoeffAcc,
    //                                    const double BoundRHSCoeffDev,
    //                                    const VectorType SpatialDefRate);

    // void ComputeStabLaplacianMatrix(MatrixType &StabLaplacianMatrix,
    //                                 const ShapeFunctionDerivativesType &rShapeDeriv,
    //                                 const double Weight) override;

    // void CalculateTauFIC(double &TauOne,
    //                      double ElemSize,
    //                      const double Density,
    //                      const double Viscosity,
    //                      const ProcessInfo &rCurrentProcessInfo) override;

    // void AddStabilizationMatrixLHS(MatrixType &rLeftHandSideMatrix,
    //                                Matrix &BulkAccMatrix,
    //                                const ShapeFunctionsType &rN,
    //                                const double Weight) override;

    // void AddStabilizationNodalTermsLHS(MatrixType &rLeftHandSideMatrix,
    //                                    const double Tau,
    //                                    const double Weight,
    //                                    const ShapeFunctionDerivativesType &rDN_DX,
    //                                    const SizeType i) override;

    // void AddStabilizationNodalTermsRHS(VectorType &rRightHandSideVector,
    //                                    const double Tau,
    //                                    const double Density,
    //                                    const double Weight,
    //                                    const ShapeFunctionDerivativesType &rDN_DX,
    //                                    const SizeType i) override;

    // void CalculateLocalContinuityEqForPressure(MatrixType &rLeftHandSideMatrix,
    //                                            VectorType &rRightHandSideVector,
    //                                            const ProcessInfo &rCurrentProcessInfo) override;

    // void GetPressureAccelerationValues(Vector &rValues,
    //                                    const int Step);

    void CalculateGeometryData(ShapeFunctionDerivativesArrayType &rDN_DX,
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

    /// Assignment operator.
    TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement &operator=(TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement const &rOther);

    /* /// Copy constructor. */
    /* TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement(TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement const& rOther); */

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
  inline std::istream &operator>>(std::istream &rIStream,
                                  TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim> &rThis)
  {
    return rIStream;
  }

  /// output stream function
  template <unsigned int TDim>
  inline std::ostream &operator<<(std::ostream &rOStream,
                                  const TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim> &rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }

} // namespace Kratos.
