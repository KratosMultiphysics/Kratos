//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:               June 2021 $
//   Revision:            $Revision:                 0.0 $
//
//

#if !defined(KRATOS_THREE_STEP_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED)
#define KRATOS_THREE_STEP_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "includes/constitutive_law.h"

#include "custom_elements/updated_lagrangian_element.h"
#include "pfem_fluid_dynamics_application_variables.h"

#include "includes/model_part.h"
/* #include "includes/node.h" */

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
  class ThreeStepUpdatedLagrangianElement : public UpdatedLagrangianElement<TDim>
  {

  protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{
    typedef struct
    {
      unsigned int voigtsize;
      // strain state
      double DetFgrad;
      double DetFgradVel;
      double DeviatoricInvariant;
      double EquivalentStrainRate;
      double VolumetricDefRate;
      Vector SpatialDefRate;
      Vector MDGreenLagrangeMaterial;
      Matrix Fgrad;
      Matrix InvFgrad;
      Matrix FgradVel;
      Matrix InvFgradVel;
      Matrix SpatialVelocityGrad;
      Matrix ConstitutiveMatrix;
      // Stress state
      double MeanPressure;
      Vector CurrentTotalCauchyStress;
      Vector UpdatedTotalCauchyStress;
      Vector CurrentDeviatoricCauchyStress;
      Vector UpdatedDeviatoricCauchyStress;

    } ElementalVariables;

  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ThreeStepUpdatedLagrangianElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ThreeStepUpdatedLagrangianElement);

    typedef UpdatedLagrangianElement<TDim> BaseType;

    /// Node type (default is: Node<3>)
    typedef Node<3> NodeType;

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

    /* typedef Element::PropertiesType::Pointer PropertiesType::Pointer; */

    typedef Element::PropertiesType PropertiesType;

    /// Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;

    /// Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    /**
       * @param NewId Index number of the new element (optional)
       */
    ThreeStepUpdatedLagrangianElement(IndexType NewId = 0) : BaseType(NewId)
    {
    }

    /// Constructor using an array of nodes.
    /**
       * @param NewId Index of the new element
       * @param ThisNodes An array containing the nodes of the new element
       */
    ThreeStepUpdatedLagrangianElement(IndexType NewId, const NodesArrayType &ThisNodes) : BaseType(NewId, ThisNodes)
    {
    }

    /// Constructor using a geometry object.
    /**
       * @param NewId Index of the new element
       * @param pGeometry Pointer to a geometry object
       */
    ThreeStepUpdatedLagrangianElement(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry)
    {
    }

    /// Constuctor using geometry and properties.
    /**
       * @param NewId Index of the new element
       * @param pGeometry Pointer to a geometry object
       * @param pProperties Pointer to the element's properties
       */
    ThreeStepUpdatedLagrangianElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /// copy constructor

    ThreeStepUpdatedLagrangianElement(ThreeStepUpdatedLagrangianElement const &rOther) : BaseType(rOther){};

    /// Destructor.
    virtual ~ThreeStepUpdatedLagrangianElement()
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
       * Returns a pointer to a new ThreeStepUpdatedLagrangianElement element, created using given input
       * @param NewId: the ID of the new element
       * @param ThisNodes: the nodes of the new element
       * @param pProperties: the properties assigned to the new element
       * @return a Pointer to the new element
       */
    Element::Pointer Create(IndexType NewId, NodesArrayType const &ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
      return Element::Pointer(new ThreeStepUpdatedLagrangianElement(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
    }

    Element::Pointer Clone(IndexType NewId, NodesArrayType const &ThisNodes) const override;

    void Initialize(const ProcessInfo &rCurrentProcessInfo) override{};

    /// Initializes the element and all geometric information required for the problem.
    void InitializeSolutionStep(const ProcessInfo &rCurrentProcessInfo) override{};

    void InitializeNonLinearIteration(const ProcessInfo &rCurrentProcessInfo) override{};

    /// Calculate the element's local contribution to the system for the current step.
    void CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                              VectorType &rRightHandSideVector,
                              const ProcessInfo &rCurrentProcessInfo) override{};

    void CalculateLeftHandSide(MatrixType &rLeftHandSideMatrix,
                               const ProcessInfo &rCurrentProcessInfo) override
    {
      KRATOS_TRY;
      KRATOS_THROW_ERROR(std::logic_error, "ThreeStepUpdatedLagrangianElement::CalculateLeftHandSide not implemented", "");
      KRATOS_CATCH("");
    }

    void CalculateRightHandSide(VectorType &rRightHandSideVector,
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

    virtual void InitializeElementalVariables(ElementalVariables &rElementalVariables){};

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
      buffer << "ThreeStepUpdatedLagrangianElement #" << this->Id();
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
    double mMaterialDeviatoricCoefficient = 0;
    double mMaterialVolumetricCoefficient = 0;
    double mMaterialDensity = 0;

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
                           const double Weight);

    void ComputeBulkMatrixLump(MatrixType &BulkMatrix,
                               const double Weight) override;

    void ComputeBulkMatrixConsistent(MatrixType &BulkMatrix,
                                     const double Weight) override;

    void ComputeBulkMatrix(MatrixType &BulkMatrix,
                           const ShapeFunctionsType &rN,
                           const double Weight) override;

    virtual void ComputeBulkMatrixRHS(MatrixType &BulkMatrix,
                                      const double Weight) override{};

    bool CalcMechanicsUpdated(ElementalVariables &rElementalVariables,
                              const ProcessInfo &rCurrentProcessInfo,
                              const ShapeFunctionDerivativesType &rDN_DX,
                              unsigned int g);

    bool CalcStrainRate(ElementalVariables &rElementalVariables,
                        const ProcessInfo &rCurrentProcessInfo,
                        const ShapeFunctionDerivativesType &rDN_DX,
                        const double theta);

    bool CalcCompleteStrainRate(ElementalVariables &rElementalVariables,
                                const ProcessInfo &rCurrentProcessInfo,
                                const ShapeFunctionDerivativesType &rDN_DX,
                                const double theta);

    virtual void CalcElasticPlasticCauchySplitted(ElementalVariables &rElementalVariables, double TimeStep,
                                                  unsigned int g, const ProcessInfo &rCurrentProcessInfo, double &Density,
                                                  double &DeviatoricCoeff, double &VolumetricCoeff){};

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
    ThreeStepUpdatedLagrangianElement &operator=(ThreeStepUpdatedLagrangianElement const &rOther);

    /* /// Copy constructor. */
    /* ThreeStepUpdatedLagrangianElement(ThreeStepUpdatedLagrangianElement const& rOther); */

    ///@}

  }; // Class ThreeStepUpdatedLagrangianElement

  ///@}

  ///@name Type Definitions
  ///@{

  ///@}
  ///@name Input and output
  ///@{

  /// input stream function
  template <unsigned int TDim>
  inline std::istream &operator>>(std::istream &rIStream,
                                  ThreeStepUpdatedLagrangianElement<TDim> &rThis)
  {
    return rIStream;
  }

  /// output stream function
  template <unsigned int TDim>
  inline std::ostream &operator<<(std::ostream &rOStream,
                                  const ThreeStepUpdatedLagrangianElement<TDim> &rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }

} // namespace Kratos.

#endif // KRATOS_THREE_STEP_UPDATED_LAGRANGIAN_ELEMENT  defined
