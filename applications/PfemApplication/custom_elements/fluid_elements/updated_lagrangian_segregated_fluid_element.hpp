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
   * Called at the beginning of each solution step
   */
  void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

  /**
   * Called at the end of eahc solution step
   */
  void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

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


  /**
   * this is called during the assembling process in order
   * to calculate the elemental mass matrix
   * @param rMassMatrix: the elemental mass matrix
   * @param rCurrentProcessInfo: the current process info instance
   */
  void CalculateMassMatrix(MatrixType& rMassMatrix,
                           ProcessInfo& rCurrentProcessInfo) override;

  /**
   * this is called during the assembling process in order
   * to calculate the elemental damping matrix
   * @param rDampingMatrix: the elemental damping matrix
   * @param rCurrentProcessInfo: the current process info instance
   */
  void CalculateDampingMatrix(MatrixType& rDampingMatrix,
                              ProcessInfo& rCurrentProcessInfo) override;



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
   * Sets process information to set member variables like mStepVariable
   */
  void SetProcessInformation(const ProcessInfo& rCurrentProcessInfo) override;

  /**
   * Initialize Element General Variables
   */
  void InitializeElementData(ElementDataType & rVariables,
                             const ProcessInfo& rCurrentProcessInfo) override;

  /**
   * Calculate Integration Step Alpha
   */
  void GetStepAlpha(double& rAlpha);

  /**
   * Calculate Element Kinematics
   */
  void CalculateKinematics(ElementDataType& rVariables,
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
   * Calculation of the Geometric Stiffness Matrix. Kvvm = BT * C * B
   */
  void CalculateAndAddKvvm(MatrixType& rLeftHandSideMatrix,
                           ElementDataType& rVariables) override;

  /**
   * Calculation of the Geometric Stiffness Matrix. Kvvg = BT * S
   */
  void CalculateAndAddKvvg(MatrixType& rLeftHandSideMatrix,
                           ElementDataType& rVariables) override;


  /**
   * Calculation of the Pressure Stiffness Matrix. Kpp
   */
  void CalculateAndAddKpp(MatrixType& rLeftHandSideMatrix,
                          ElementDataType& rVariables);

  /**
   * Calculation of the Internal Forces Vector. Fi = B * sigma
   */
  void CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
                                     ElementDataType & rVariables) override;

  /**
   * Calculation of the Pressure Vector.
   */
  void CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
                                     ElementDataType & rVariables);

  /**
   * Add volumetric part to Stress Vector
   */
  void AddVolumetricPart(Vector& rStressVector, double& rMeanPressure);

  /**
   * Remove volumetric part from Stress Vector
   */
  void RemoveVolumetricPart(Vector& rStressVector, double& rMeanPressure);

  /**
   * Add volumetric part to Constitutive Matrix
   */
  void AddVolumetricPart(Matrix& rConstitutiveMatrix, double& rBulkFactor);

  /**
   * Remove volumetric part from Constitutive Matrix
   */
  void RemoveVolumetricPart(Matrix& rConstitutiveMatrix, double& rBulkFactor);

  /**
   * Calculation of the Mean value considering a Dense Matrix.
   */
  void CalculateDenseMatrixMeanValue(MatrixType& rMatrix, double& rMeanValue);

  /**
   * Calculation of the Mean value considering a Lumped Matrix.
   */
  void CalculateLumpedMatrixMeanValue(MatrixType& rMatrix, double& rMeanValue);



  /**
   * Set Variables of the Element to the Parameters of the Constitutive Law
   */
  void SetElementData(ElementDataType& rVariables,
                      ConstitutiveLaw::Parameters& rValues,
                      const int & rPointNumber) override;

  /**
   * Set Parameters for the Constitutive Law and Calculate Material Response
   */
  void CalculateMaterialResponse(ElementDataType& rVariables,
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
   * Get Face normal
   */
  void GetFaceNormal(const std::vector<SizeType>& rFace, const ElementDataType & rVariables, Vector& rNormal);

  /**
   * Get Face weight
   */
  void GetFaceWeight(const std::vector<SizeType>& rFace, const ElementDataType & rVariables, double& rWeight, double& rNormalSize);

  /**
   * Get Face normal
   */
  void GetFaceNormal(const std::vector<SizeType>& rFace, Vector& rNormal);


  /**
   * Calculation of the Volume Change of the Element
     */
  double& CalculateVolumeChange(double& rVolumeChange, ElementDataType& rVariables) override;


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
