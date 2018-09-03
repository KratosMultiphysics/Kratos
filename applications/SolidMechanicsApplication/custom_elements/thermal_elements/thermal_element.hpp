//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_THERMAL_ELEMENT_H_INCLUDED )
#define  KRATOS_THERMAL_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "includes/constitutive_law.h"


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



class KRATOS_API(SOLID_MECHANICS_APPLICATION) ThermalElement
    : public Element
{
 public:

  ///@name Type Definitions_
  ///@{
  ///Reference type definition for constitutive laws
  typedef ConstitutiveLaw ConstitutiveLawType;
  ///Pointer type for constitutive laws
  typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
  ///StressMeasure from constitutive laws
  typedef ConstitutiveLawType::StressMeasure StressMeasureType;
  ///Type definition for integration methods
  typedef GeometryData::IntegrationMethod IntegrationMethod;


  /// Counted pointer of ThermalElement
  KRATOS_CLASS_POINTER_DEFINITION(ThermalElement);


 protected:

  /**
   * Flags related to the element computation
   */

  KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR );
  KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX );

  /**
   * Parameters to be used in the Element as they are. Direct interface to Parameters Struct
   */

  struct GeneralVariables
  {

    //for axisymmetric use only
    double  CurrentRadius;
    double  ReferenceRadius;

    //general thermal variables
    double  DeltaTime;
    double  HeatCapacity;
    double  HeatConductivity;
    double  PlasticDissipation;
    double  DeltaPlasticDissipation;

    //general mechanical variables
    double  detJ;
    Vector  N;
    Matrix  DN_DX;

    //variables including all integration points
    const GeometryType::ShapeFunctionsGradientsType* pDN_De;
    GeometryType::JacobiansType J;
    GeometryType::JacobiansType j;
    const Matrix* pNcontainer;
    Matrix  DeltaPosition;

    /**
     * sets the value of a specified pointer variable
     */

    void SetShapeFunctionsGradients(const GeometryType::ShapeFunctionsGradientsType &rDN_De)
    {
      pDN_De=&rDN_De;
    };

    void SetShapeFunctions(const Matrix& rNcontainer)
    {
      pNcontainer=&rNcontainer;
    };


    /**
     * returns the value of a specified pointer variable
     */

    const GeometryType::ShapeFunctionsGradientsType& GetShapeFunctionsGradients()
    {
      return *pDN_De;
    };

    const Matrix& GetShapeFunctions()
    {
      return *pNcontainer;
    };


  };


 public:

  ///@}
  ///@name Life Cycle
  ///@{


  /// Default constructors.
  ThermalElement(IndexType NewId, GeometryType::Pointer pGeometry);
  ThermalElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

  ///Copy constructor
  ThermalElement(ThermalElement const& rOther);


  /// Destructor.
  ~ThermalElement() override;

  ///@}
  ///@name Operators
  ///@{

  /// Assignment operator.
  ThermalElement& operator=(ThermalElement const& rOther);


  ///@}
  ///@name Operations
  ///@{

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
  //Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const;


  //************* GETTING METHODS

  /**
   * Returns the currently selected integration method
   * @return current integration method selected
   */
  IntegrationMethod GetIntegrationMethod() const override;

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


  //on integration points:
  /**
   * Access for variables on Integration points.
   * This gives access to variables stored in the constitutive law on each integration point.
   * Specialisations of element.h (e.g. the TotalLagrangian) must specify the actual
   * interface to the constitutive law!
   * Note, that these functions expect a std::vector of values for the
   * specified variable type that contains a value for each integration point!
   * SetValueOnIntegrationPoints: set the values for given Variable.
   * GetValueOnIntegrationPoints: get the values for given Variable.
   */

  //SET
  /**
   * Set a double  Value on the Element Constitutive Law
   */
  void SetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

  /**
   * Set a Vector Value on the Element Constitutive Law
   */
  void SetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

  /**
   * Set a Matrix Value on the Element Constitutive Law
   */
  void SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;


  //GET:
  /**
   * Get on rVariable a double Value from the Element Constitutive Law
   */
  void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

  /**
   * Get on rVariable a Vector Value from the Element Constitutive Law
   */
  void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

  /**
   * Get on rVariable a Matrix Value from the Element Constitutive Law
   */
  void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;


  //************* STARTING - ENDING  METHODS

  /**
   * Called to initialize the element.
   * Must be called before any calculation is done
   */
  void Initialize() override;


  /**
   * Called at the beginning of each solution step
   */
  void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;


  /**
   * Called at the end of each solution step
   */
  void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;



  //************* COMPUTING  METHODS

  /**
   * this is called during the assembling process in order
   * to calculate all elemental contributions to the global system
   * matrix and the right hand side
   * @param rLeftHandSideMatrix: the elemental left hand side matrix
   * @param rRightHandSideVector: the elemental right hand side
   * @param rCurrentProcessInfo: the current process info instance
   */

  void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

  /**
   * this is called during the assembling process in order
   * to calculate the elemental right hand side vector only
   * @param rRightHandSideVector: the elemental right hand side vector
   * @param rCurrentProcessInfo: the current process info instance
   */
  void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

  /**
   * this is called during the assembling process in order
   * to calculate the elemental left hand side vector only
   * @param rLeftHandSideVector: the elemental left hand side vector
   * @param rCurrentProcessInfo: the current process info instance
   */
  void CalculateLeftHandSide (MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) override;

  /**
   * this is called during the assembling process in order
   * to calculate the elemental mass matrix
   * @param rMassMatrix: the elemental mass matrix
   * @param rCurrentProcessInfo: the current process info instance
   */
  void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override;

  /**
   * this is called during the assembling process in order
   * to calculate the elemental damping matrix
   * @param rDampingMatrix: the elemental damping matrix
   * @param rCurrentProcessInfo: the current process info instance
   */
  void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override;


  //on integration points:
  /**
   * Calculate a double Variable on the Element Constitutive Law
   */
  void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

  /**
   * Calculate a Vector Variable on the Element Constitutive Law
   */
  void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

  /**
   * Calculate a Matrix Variable on the Element Constitutive Law
   */
  void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo) override;


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

  //std::string Info() const;

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
  //      virtual String Info() const;

  /// Print information about this object.
  //      virtual void PrintInfo(std::ostream& rOStream) const;

  /// Print object's data.
  //      virtual void PrintData(std::ostream& rOStream) const;
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

  /**
   * Currently selected integration methods
   */
  IntegrationMethod mThisIntegrationMethod;

  ///@}
  ///@name Protected Operators
  ///@{
  ThermalElement() : Element()
  {
  }

  /**
   * Calculates the elemental contributions
   * \f$ K^e = w\,B^T\,D\,B \f$ and
   * \f$ r^e \f$
   */
  virtual void CalculateElementalSystem(MatrixType& rLeftHandSideMatrix,
                                        VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo,
                                        Flags & rCalculationFlags);
  ///@}
  ///@name Protected Operations
  ///@{


  /**
   * Calculation and addition of the matrices of the LHS
   */

  virtual void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                                  GeneralVariables& rVariables,
                                  double& rIntegrationWeight);

  /**
   * Calculation and addition of the vectors of the RHS
   */

  virtual void CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                  GeneralVariables& rVariables,
                                  double& rHeatSource,
                                  double& rIntegrationWeight);



  /**
   * Calculation of the Material Stiffness Matrix. Kthermal = BT * k * B
   */
  virtual void CalculateAndAddKthermal(MatrixType& rK,
                                       GeneralVariables & rVariables,
                                       double& rIntegrationWeight
                                       );

  /**
   * Calculation of the Material Stiffness Matrix. H = Ni * Nj * h
   */
  virtual void CalculateAndAddHthermal(MatrixType& rK,
                                       GeneralVariables & rVariables,
                                       double& rIntegrationWeight
                                       );


  /**
   * Calculation of the Material Stiffness Matrix. Mthermal = Ni * Nj * c * 1/Delta_t
   */
  virtual void CalculateAndAddMthermal(MatrixType& rM,
                                       GeneralVariables & rVariables,
                                       double& rIntegrationWeight
                                       );


  /**
   * Calculation of the External Forces Vector. Fe = N * t + N * b
   */
  virtual void CalculateAndAddExternalForces(GeneralVariables& rVariables,
                                             VectorType& rRightHandSideVector,
                                             double& rHeatSource,
                                             double& rIntegrationWeight
                                             );


  /**
   * Calculation of the Thermal Forces due to sigma
   */
  virtual void CalculateAndAddThermalForces(GeneralVariables & rVariables,
                                            VectorType& rRightHandSideVector,
                                            double& rIntegrationWeight
                                            );


  /**
   * Initialize System Matrices
   */
  virtual void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
                                        VectorType& rRightHandSideVector,
                                        Flags& rCalculationFlags);

  /**
   * Initialize Variables
   */
  void InitializeVariables ();

  /**
   * Calculate Element Kinematics
   */
  virtual void CalculateKinematics(GeneralVariables& rVariables,
                                   const double& rPointNumber);


  /**
   * Calculation of the Position increment
   */
  Matrix& CalculateDeltaPosition(Matrix & rDeltaPosition);


  /**
   * Calculate Variation of Thermal Properties
   */
  void CalculateThermalProperties(GeneralVariables& rVariables);



  /**
   * Calculation integration weights W
   */
  virtual double& CalculateIntegrationWeight(double & rIntegrationWeight);


  /**
   * Initialize Element General Variables
   */
  virtual void InitializeGeneralVariables(GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo);



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

  void save(Serializer& rSerializer) const override;

  void load(Serializer& rSerializer) override;

  ///@name Private Inquiry
  ///@{
  ///@}
  ///@name Un accessible methods
  ///@{
  ///@}

}; // Class ThermalElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_THERMAL_ELEMENT_H_INCLUDED  defined
