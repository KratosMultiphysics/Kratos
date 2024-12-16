//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:              April 2018 $
//   Revision:            $Revision:                 0.0 $
//
//

#if !defined(KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED)
#define KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

#include "pfem_fluid_dynamics_application_variables.h"

#include "custom_elements/updated_lagrangian_element.h"

#include "includes/model_part.h"

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
  class TwoStepUpdatedLagrangianElement : public UpdatedLagrangianElement<TDim>
  {

  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TwoStepUpdatedLagrangianElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TwoStepUpdatedLagrangianElement);

    typedef UpdatedLagrangianElement<TDim> BaseType;

    /// Node type (default is: Node)
    typedef Node NodeType;

    /// Geometry type (using with given NodeType)
    typedef Geometry<NodeType> GeometryType;

    /// Definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    /// Vector type for local contributions to the linear system
    typedef Vector VectorType;

    /// Matrix type for local contributions to the linear system
    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector<Dof<double>::Pointer> DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    /// Type for shape function values container
    typedef Kratos::Vector ShapeFunctionsType;

    /// Type for a matrix containing the shape function gradients
    typedef Kratos::Matrix ShapeFunctionDerivativesType;

    /// Type for an array of shape function gradient matrices
    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    /// Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;

    /// Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    typedef typename BaseType::PropertiesType PropertiesType;

    typedef typename BaseType::PropertiesType::Pointer pPropertiesType;

    typedef typename BaseType::ElementalVariables ElementalVariables;

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    TwoStepUpdatedLagrangianElement(IndexType NewId = 0) : BaseType(NewId)
    {
    }

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    TwoStepUpdatedLagrangianElement(IndexType NewId, const NodesArrayType &ThisNodes) : BaseType(NewId, ThisNodes)
    {
    }

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    TwoStepUpdatedLagrangianElement(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry)
    {
    }

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    TwoStepUpdatedLagrangianElement(IndexType NewId, GeometryType::Pointer pGeometry, pPropertiesType pProperties) : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /// copy constructor

    TwoStepUpdatedLagrangianElement(TwoStepUpdatedLagrangianElement const &rOther) : BaseType(rOther){};

    /// Destructor.
    virtual ~TwoStepUpdatedLagrangianElement()
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
     * Returns a pointer to a new TwoStepUpdatedLagrangianElement element, created using given input
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const &ThisNodes,
                            pPropertiesType pProperties) const override
    {
      return Element::Pointer(new TwoStepUpdatedLagrangianElement(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
    }

    Element::Pointer Clone(IndexType NewId, NodesArrayType const &ThisNodes) const override;

    /// Calculate the element's local contribution to the system for the current step.
    void CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                              VectorType &rRightHandSideVector,
                              const ProcessInfo &rCurrentProcessInfo) override{};

    // The following methods have different implementations depending on TDim
    /// Provides the global indices for each one of this element's local rows
    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    void EquationIdVector(EquationIdVectorType &rResult,
                          const ProcessInfo &rCurrentProcessInfo) const override;

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType &rElementalDofList,
                    const ProcessInfo &rCurrentProcessInfo) const override;

    virtual void UpdateCauchyStress(unsigned int g, const ProcessInfo &rCurrentProcessInfo) override{};

    virtual void InitializeElementalVariables(ElementalVariables &rElementalVariables) override{};

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
    int Check(const ProcessInfo &rCurrentProcessInfo) const override { return 0; };

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
      buffer << "TwoStepUpdatedLagrangianElement #" << this->Id();
      return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override {}

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

    virtual void CalculateLocalMomentumEquations(MatrixType &rLeftHandSideMatrix,
                                                 VectorType &rRightHandSideVector,
                                                 const ProcessInfo &rCurrentProcessInfo) override{};

    virtual void CalculateLocalContinuityEqForPressure(MatrixType &rLeftHandSideMatrix,
                                                       VectorType &rRightHandSideVector,
                                                       const ProcessInfo &rCurrentProcessInfo) override{};

    virtual void CalculateExplicitContinuityEquation(MatrixType &rLeftHandSideMatrix,
                                                     VectorType &rRightHandSideVector,
                                                     const ProcessInfo &rCurrentProcessInfo) override{};

    virtual double GetThetaMomentum() override
    {
      std::cout << "I SHOULD NOT ENTER HERE!" << std::endl;
      return 0.0;
    };

    virtual double GetThetaContinuity() override
    {
      std::cout << "I SHOULD NOT ENTER HERE!" << std::endl;
      return 0.0;
    };

    void CalculateMassMatrix(Matrix &rMassMatrix,
                             const ProcessInfo &rCurrentProcessInfo) override{};

    void ComputeMassMatrix(Matrix &rMassMatrix,
                           const ShapeFunctionsType &rN,
                           const double Weight,
                           double &MeanValue) override;

    void ComputeLumpedMassMatrix(Matrix &rMassMatrix,
                                 const double Weight,
                                 double &MeanValue) override;

    void AddExternalForces(Vector &rRHSVector,
                           const double Density,
                           const ShapeFunctionsType &rN,
                           const double Weight) override;

    void AddInternalForces(Vector &rRHSVector,
                           const ShapeFunctionDerivativesType &rDN_DX,
                           ElementalVariables &rElementalVariables,
                           const double Weight) override;

    void ComputeBulkMatrixLump(MatrixType &BulkMatrix,
                               const double Weight) override;

    virtual void ComputeBulkMatrixRHS(MatrixType &BulkMatrix,
                                      const double Weight) override{};

    virtual void CalcElasticPlasticCauchySplitted(
        ElementalVariables &rElementalVariables,
        const unsigned int g,
        const Vector& rN,
        const ProcessInfo &rCurrentProcessInfo,
        double &Density,
        double &DeviatoricCoeff,
        double &VolumetricCoeff) override
    {};

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
    TwoStepUpdatedLagrangianElement &operator=(TwoStepUpdatedLagrangianElement const &rOther);

    /* /// Copy constructor. */
    /* TwoStepUpdatedLagrangianElement(TwoStepUpdatedLagrangianElement const& rOther); */

    ///@}

  }; // Class TwoStepUpdatedLagrangianElement

  ///@}

  ///@name Type Definitions
  ///@{

  ///@}
  ///@name Input and output
  ///@{

  /// input stream function
  template <unsigned int TDim>
  inline std::istream &operator>>(std::istream &rIStream,
                                  TwoStepUpdatedLagrangianElement<TDim> &rThis)
  {
    return rIStream;
  }

  /// output stream function
  template <unsigned int TDim>
  inline std::ostream &operator<<(std::ostream &rOStream,
                                  const TwoStepUpdatedLagrangianElement<TDim> &rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }

} // namespace Kratos.

#endif // KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_ELEMENT  defined
