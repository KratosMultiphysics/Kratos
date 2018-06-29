//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:               April 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_UPDATED_LAGRANGIAN_SEGREGATED_FLUID_ELEMENT_H_INCLUDED)
#define  KRATOS_UPDATED_LAGRANGIAN_SEGREGATED_FLUID_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/fluid_elements/fluid_element.hpp"


namespace Kratos
{
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

/// Updated Lagrangian Segregated Fluid Element for 3D and 2D geometries

/**
 * Implements a Large Displacement Lagrangian definition for fluid analysis.
 * This works for linear Triangles and Tetrahedra (base class)
 */

class KRATOS_API(PFEM_APPLICATION) UpdatedLagrangianSegregatedFluidElement : public FluidElement
{
public:

  ///@name Type Definitions
  ///@{
  ///Reference type definition for constitutive laws
  typedef ConstitutiveLaw ConstitutiveLawType;
  ///Pointer type for constitutive laws
  typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
  ///StressMeasure from constitutive laws
  typedef ConstitutiveLawType::StressMeasure StressMeasureType;
  ///Type definition for integration methods
  typedef GeometryData::IntegrationMethod IntegrationMethod;
  ///Type for size
  typedef GeometryData::SizeType SizeType;

  /// Counted pointer of UpdatedLagrangianSegregatedFluidElement
  KRATOS_CLASS_POINTER_DEFINITION( UpdatedLagrangianSegregatedFluidElement );
  ///@}

  enum StepType{VELOCITY_STEP = 0, PRESSURE_STEP = 1};

  ///@name Life Cycle
  ///@{

  /// Empty constructor needed for serialization
  UpdatedLagrangianSegregatedFluidElement();

  /// Default constructors
  UpdatedLagrangianSegregatedFluidElement(IndexType NewId, GeometryType::Pointer pGeometry);

  UpdatedLagrangianSegregatedFluidElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

  ///Copy constructor
  UpdatedLagrangianSegregatedFluidElement(UpdatedLagrangianSegregatedFluidElement const& rOther);


  /// Destructor.
  virtual ~UpdatedLagrangianSegregatedFluidElement();

  ///@}
  ///@name Operators
  ///@{

  /// Assignment operator.
  UpdatedLagrangianSegregatedFluidElement& operator=(UpdatedLagrangianSegregatedFluidElement const& rOther);

  ///@}
  ///@name Operations
  ///@{
  /**
   * Returns the currently selected integration method
   * @return current integration method selected
   */
  /**
   * creates a new total lagrangian updated element pointer
   * @param NewId: the ID of the new element
   * @param ThisNodes: the nodes of the new element
   * @param pProperties: the properties assigned to the new element
   * @return a Pointer to the new element
   */
  Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

  /**
   * clones the selected element variables, creating a new one
   * @param NewId: the ID of the new element
   * @param ThisNodes: the nodes of the new element
   * @param pProperties: the properties assigned to the new element
   * @return a Pointer to the new element
   */
  Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

  //************* STARTING - ENDING  METHODS

  /**
   * this is called for non-linear analysis at the beginning of the iteration process
   */
  void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

  /**
   * this is called for non-linear analysis at the beginning of the iteration process
   */
  void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

  //************* GETTING METHODS

  /**
   * Sets on rElementalDofList the degrees of freedom of the considered element geometry
   */
  void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override;

  /**
   * Sets on rResult the ID's of the element degrees of freedom
   */
  void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

  /**
   * Sets on rValues the nodal displacements
   */
  void GetValuesVector(Vector& rValues, int Step = 0) override;

  /**
   * Sets on rValues the nodal velocities
   */
  void GetFirstDerivativesVector(Vector& rValues, int Step = 0) override;

  /**
   * Sets on rValues the nodal accelerations
   */
  void GetSecondDerivativesVector(Vector& rValues, int Step = 0) override;

  //************************************************************************************
  //************************************************************************************
  /**
   * This function provides the place to perform checks on the completeness of the input.
   * It is designed to be called only once (or anyway, not often) typically at the beginning
   * of the calculations, so to verify that nothing is missing from the input
   * or that no common error is found.
   * @param rCurrentProcessInfo
   */
  int Check(const ProcessInfo& rCurrentProcessInfo) override;

  ///@}
  ///@name Access
  ///@{

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
    buffer << "Updated Lagrangian Fluid Element #" << Id();
    return buffer.str();
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "Updated Lagrangian Fluid Element #" << Id();
  }

  /// Print object's data.
  void PrintData(std::ostream& rOStream) const override
  {
    GetGeometry().PrintData(rOStream);
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

  StepType mStepVariable;

  ///@}
  ///@name Protected Operators
  ///@{
  ///@}
  ///@name Protected Operations
  ///@{

  /**
   * Calculate Element Kinematics
   */
  void CalculateKinematics(ElementDataType& rVariables,
                           const double& rPointNumber) override;

  /**
   * Calculate Element Jacobian
   */
  void CalculateKinetics(ElementDataType& rVariables,
                         const double& rPointNumber) override;

  /**
   * Calculation and addition of the matrices of the LHS
   */
  void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                          ElementDataType& rVariables) override;

  /**
   * Calculation and addition of the matrices of the RHS
   */
  void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
                          ElementDataType& rVariables) override;

  /**
   * Get element size from the dofs
   */
  unsigned int GetDofsSize() override;


  /**
   * Calculation of the Geometric Stiffness Matrix. Kvvg = BT * S
   */
  void CalculateAndAddKvvg(MatrixType& rLeftHandSideMatrix,
                           ElementDataType& rVariables) override;

  /**
   * Calculation of the Bulk Matrix.
   */
  void CalculateAndAddKBulk(MatrixType& rLeftHandSideMatrix,
                            ElementDataType& rVariables);

  /**
   * Calculation of the Pressure Stiffness Matrix. Kpp
   */
  void CalculateAndAddKpp(MatrixType& rLeftHandSideMatrix,
                          ElementDataType& rVariables);

  /**
   * Calculation of the Pressure Vector.
   */
  void CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
                                     ElementDataType & rVariables);

  /**
   * Set Variables of the Element to the Parameters of the Constitutive Law
   */
  void SetElementData(ElementDataType& rVariables,
                      ConstitutiveLaw::Parameters& rValues,
                      const int & rPointNumber) override;

  /**
   * Calculate stabilization factor
   */
  void CalculateStabilizationTau(ElementDataType& rVariables);

  /**
   * Get Faces with nodes in the freesurface
   */
  void GetFreeSurfaceFaces(std::vector<std::vector<SizeType> >& Faces);

  /**
   * Get Historical variables
   */
  void GetHistoricalVariables(ElementDataType& rVariables, const double& rPointNumber);


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
  ///@name Private Operators
  ///@{
  ///@}
  ///@name Private Operations
  ///@{
  ///@}
  ///@name Private  Access
  ///@{
  ///@}
  ///@}
  ///@name Serialization
  ///@{
  friend class Serializer;

  // A private default constructor necessary for serialization

  void save(Serializer& rSerializer) const override;

  void load(Serializer& rSerializer) override;


  ///@name Private Inquiry
  ///@{
  ///@}
  ///@name Un accessible methods
  ///@{
  ///@}

}; // Class UpdatedLagrangianSegregatedFluidElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_UPDATED_LAGRANGIAN_SEGREGATED_FLUID_ELEMENT_H_INCLUDED  defined
