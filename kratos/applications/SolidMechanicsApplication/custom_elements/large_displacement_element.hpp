//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_LARGE_DISPLACEMENT_ELEMENT_H_INCLUDED )
#define  KRATOS_LARGE_DISPLACEMENT_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "custom_utilities/comparison_utilities.hpp"


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

/// Large Displacement Lagrangian Element for 3D and 2D geometries. (base class)

/**
 * Implements a Large Displacement Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D (base class)
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) LargeDisplacementElement
    : public Element
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

    /// Counted pointer of LargeDisplacementElement
    KRATOS_CLASS_POINTER_DEFINITION( LargeDisplacementElement );
    ///@}

protected:

    /**
     * Flags related to the element computation
     */

    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR_WITH_COMPONENTS );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX_WITH_COMPONENTS );

    /**
     * Parameters to be used in the Element as they are. Direct interface to Parameters Struct
     */

    struct GeneralVariables
    {
      private:

        //variables including all integration points
        const GeometryType::ShapeFunctionsGradientsType* pDN_De;
        const Matrix* pNcontainer;

      public:

        StressMeasureType StressMeasure;

        //for axisymmetric use only
        double  CurrentRadius;
        double  ReferenceRadius;

        //general variables for large displacement use
        double  detF;
        double  detF0;
        double  detFT;
        double  detJ;
        Vector  StrainVector;
        Vector  StressVector;
        Vector  N;
        Matrix  B;
        Matrix  F;    //F from the reference to the current configuration ( Delta F )
        Matrix  F0;   //F in the reference configuration, ( historical F )
        Matrix  FT;   //FT = F0 * F  ( total F )
        Matrix  DN_DX;
        Matrix  ConstitutiveMatrix;

        //variables including all integration points
        GeometryType::JacobiansType J;
        GeometryType::JacobiansType j;
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

        void Initialize( const unsigned int& voigt_size, 
			 const unsigned int& dimension, 
			 const unsigned int& number_of_nodes )
        {
	  StressMeasure = ConstitutiveLaw::StressMeasure_PK2;
	  //doubles
	  //radius
	  CurrentRadius = 0;
	  ReferenceRadius = 0;
	  //jacobians
	  detF  = 1;
	  detF0 = 1;
	  detFT = 1;
	  detJ  = 1;
	  //vectors
	  StrainVector.resize(voigt_size,false);
          StressVector.resize(voigt_size,false);
	  N.resize(number_of_nodes,false);	  
	  noalias(StrainVector) = ZeroVector(voigt_size);
	  noalias(StressVector) = ZeroVector(voigt_size);
	  noalias(N) = ZeroVector(number_of_nodes);
	  //matrices
	  B.resize(voigt_size, dimension*number_of_nodes,false);
	  F.resize(dimension,dimension,false);
	  F0.resize(dimension,dimension,false);
	  FT.resize(dimension,dimension,false);
	  DN_DX.resize(number_of_nodes, dimension,false);
	  ConstitutiveMatrix.resize(voigt_size, voigt_size,false);
	  DeltaPosition.resize(number_of_nodes, dimension,false);

	  noalias(B)  = ZeroMatrix(voigt_size, dimension*number_of_nodes);
	  noalias(F)  = IdentityMatrix(dimension);
	  noalias(F0) = IdentityMatrix(dimension);
	  noalias(FT) = IdentityMatrix(dimension);
	  noalias(DN_DX) = ZeroMatrix(number_of_nodes, dimension);
	  noalias(ConstitutiveMatrix) = ZeroMatrix(voigt_size, voigt_size);
	  noalias(DeltaPosition) = ZeroMatrix(number_of_nodes, dimension);
	  //others
	  J.resize(1,false);
	  j.resize(1,false);
	  noalias(J[0]) = ZeroMatrix(1,1);
	  noalias(j[0]) = ZeroMatrix(1,1);
	  //pointers
	  pDN_De = NULL;
	  pNcontainer = NULL;
	}


    };


    /**
     * This struct is used in the component wise calculation only
     * is defined here and is used to declare a member variable in the component wise elements
     * private pointers can only be accessed by means of set and get functions
     * this allows to set and not copy the local system variables
     */

    struct LocalSystemComponents
    {
    private:
      
      //for calculation local system with compacted LHS and RHS 
      MatrixType *mpLeftHandSideMatrix;
      VectorType *mpRightHandSideVector;

      //for calculation local system with LHS and RHS components 
      std::vector<MatrixType> *mpLeftHandSideMatrices;
      std::vector<VectorType> *mpRightHandSideVectors;

      //LHS variable components 
      const std::vector< Variable< MatrixType > > *mpLeftHandSideVariables;

      //RHS variable components 
      const std::vector< Variable< VectorType > > *mpRightHandSideVariables;

    
    public:

      //calculation flags
      Flags  CalculationFlags;

      /**
       * sets the value of a specified pointer variable
       */
      void SetLeftHandSideMatrix( MatrixType& rLeftHandSideMatrix ) { mpLeftHandSideMatrix = &rLeftHandSideMatrix; };
      void SetLeftHandSideMatrices( std::vector<MatrixType>& rLeftHandSideMatrices ) { mpLeftHandSideMatrices = &rLeftHandSideMatrices; };
      void SetLeftHandSideVariables(const std::vector< Variable< MatrixType > >& rLeftHandSideVariables ) { mpLeftHandSideVariables = &rLeftHandSideVariables; }; 

      void SetRightHandSideVector( VectorType& rRightHandSideVector ) { mpRightHandSideVector = &rRightHandSideVector; };
      void SetRightHandSideVectors( std::vector<VectorType>& rRightHandSideVectors ) { mpRightHandSideVectors = &rRightHandSideVectors; };
      void SetRightHandSideVariables(const std::vector< Variable< VectorType > >& rRightHandSideVariables ) { mpRightHandSideVariables = &rRightHandSideVariables; }; 

 
      /**
       * returns the value of a specified pointer variable
       */
      MatrixType& GetLeftHandSideMatrix() { return *mpLeftHandSideMatrix; };
      std::vector<MatrixType>& GetLeftHandSideMatrices() { return *mpLeftHandSideMatrices; };
      const std::vector< Variable< MatrixType > >& GetLeftHandSideVariables() { return *mpLeftHandSideVariables; }; 

      VectorType& GetRightHandSideVector() { return *mpRightHandSideVector; };
      std::vector<VectorType>& GetRightHandSideVectors() { return *mpRightHandSideVectors; };
      const std::vector< Variable< VectorType > >& GetRightHandSideVariables() { return *mpRightHandSideVariables; }; 

    };


public:


    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    LargeDisplacementElement();

    /// Default constructors
    LargeDisplacementElement(IndexType NewId, GeometryType::Pointer pGeometry);

    LargeDisplacementElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    LargeDisplacementElement(LargeDisplacementElement const& rOther);

    /// Destructor.
    virtual ~LargeDisplacementElement();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    LargeDisplacementElement& operator=(LargeDisplacementElement const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

    /**
     * clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const;


    //************* GETTING METHODS

    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
    IntegrationMethod GetIntegrationMethod() const;

    /**
     * Sets on rElementalDofList the degrees of freedom of the considered element geometry
     */
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo);

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    /**
     * Sets on rValues the nodal displacements
     */
    void GetValuesVector(Vector& rValues, int Step = 0);

    /**
     * Sets on rValues the nodal velocities
     */
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0);

    /**
     * Sets on rValues the nodal accelerations
     */
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0);



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
    virtual void SetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Set a Vector Value on the Element Constitutive Law
     */
    void SetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Set a Matrix Value on the Element Constitutive Law
     */
    void SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);

    /**
    * Set a Constitutive Law Value
    */
    void SetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
                                      std::vector<ConstitutiveLaw::Pointer>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo );


    //GET:
    /**
     * Get on rVariable a double Value from the Element Constitutive Law
     */
    virtual void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Get on rVariable a Vector Value from the Element Constitutive Law
     */
    virtual void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Get on rVariable a Matrix Value from the Element Constitutive Law
     */
    virtual void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Get a Constitutive Law Value
     */
    void GetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
                                      std::vector<ConstitutiveLaw::Pointer>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo );



    //************* STARTING - ENDING  METHODS

    /**
      * Called to initialize the element.
      * Must be called before any calculation is done
      */
    virtual void Initialize();

    /**
     * Called at the beginning of each solution step
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);

    /**
     * this is called for non-linear analysis at the beginning of the iteration process
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);

    /**
     * this is called for non-linear analysis at the beginning of the iteration process
     */
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);

    /**
     * Called at the end of eahc solution step
     */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);


    //************* COMPUTING  METHODS


    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, 
			      VectorType& rRightHandSideVector, 
			      ProcessInfo& rCurrentProcessInfo);

    /**
     * this function provides a more general interface to the element.
     * it is designed so that rLHSvariables and rRHSvariables are passed TO the element
     * thus telling what is the desired output
     * @param rLeftHandSideMatrices: container with the output left hand side matrices
     * @param rLHSVariables: paramter describing the expected LHSs
     * @param rRightHandSideVectors: container for the desired RHS output
     * @param rRHSVariables: parameter describing the expected RHSs
     */
    void CalculateLocalSystem(std::vector< MatrixType >& rLeftHandSideMatrices,
			      const std::vector< Variable< MatrixType > >& rLHSVariables,
			      std::vector< VectorType >& rRightHandSideVectors,
			      const std::vector< Variable< VectorType > >& rRHSVariables,
			      ProcessInfo& rCurrentProcessInfo);

    /**
      * this is called during the assembling process in order
      * to calculate the elemental right hand side vector only
      * @param rRightHandSideVector: the elemental right hand side vector
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateRightHandSide(VectorType& rRightHandSideVector, 
				ProcessInfo& rCurrentProcessInfo);

    /**
     * this function provides a more general interface to the element.
     * it is designed so that rRHSvariables are passed TO the element
     * thus telling what is the desired output
     * @param rRightHandSideVectors: container for the desired RHS output
     * @param rRHSVariables: parameter describing the expected RHSs
     */
    void CalculateRightHandSide(std::vector< VectorType >& rRightHandSideVectors,
				const std::vector< Variable< VectorType > >& rRHSVariables,
				ProcessInfo& rCurrentProcessInfo);

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side vector only
     * @param rLeftHandSideVector: the elemental left hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLeftHandSide (MatrixType& rLeftHandSideMatrix, 
				ProcessInfo& rCurrentProcessInfo);

    /**
     * this is called during the assembling process in order
     * to calculate the first derivatives contributions for the LHS and RHS
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateFirstDerivativesContributions(MatrixType& rLeftHandSideMatrix,
						VectorType& rRightHandSideVector,
						ProcessInfo& rCurrentProcessInfo);

    /**
     * this is called during the assembling process in order
     * to calculate the second derivatives contributions for the LHS and RHS 
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix,
						VectorType& rRightHandSideVector,
						ProcessInfo& rCurrentProcessInfo);

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix for the second derivatives constributions
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
				       ProcessInfo& rCurrentProcessInfo);


    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector for the second derivatives constributions
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
				       ProcessInfo& rCurrentProcessInfo);

    /**
      * this is called during the assembling process in order
      * to calculate the elemental mass matrix
      * @param rMassMatrix: the elemental mass matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateMassMatrix(MatrixType& rMassMatrix, 
		    ProcessInfo& rCurrentProcessInfo);

    /**
      * this is called during the assembling process in order
      * to calculate the elemental damping matrix
      * @param rDampingMatrix: the elemental damping matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateDampingMatrix(MatrixType& rDampingMatrix, 
		    ProcessInfo& rCurrentProcessInfo);


    /**
     * this function is designed to make the element to assemble an rRHS vector
     * identified by a variable rRHSVariable by assembling it to the nodes on the variable
     * rDestinationVariable.
     * @param rRHSVector: input variable containing the RHS vector to be assembled
     * @param rRHSVariable: variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable: variable in the database to which the rRHSvector will be assembled 
      * @param rCurrentProcessInfo: the current process info instance
     */      
    virtual void AddExplicitContribution(const VectorType& rRHSVector, 
					 const Variable<VectorType>& rRHSVariable, 
					 Variable<array_1d<double,3> >& rDestinationVariable, 
					 const ProcessInfo& rCurrentProcessInfo);

    //on integration points:
    /**
     * Calculate a double Variable on the Element Constitutive Law
     */
    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculate a Vector Variable on the Element Constitutive Law
     */
    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculate a Matrix Variable on the Element Constitutive Law
     */
    void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo);


    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo);


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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Large Displacement Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Large Displacement Element #" << Id();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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

    /**
     * Currently selected integration methods
     */
    IntegrationMethod mThisIntegrationMethod;

    /**
     * Container for constitutive law instances on each integration point
     */
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;


    /**
     * Finalize and Initialize label
     */
    bool mFinalizedStep;

	
    ///@}
    ///@name Protected Operators
    ///@{

 
    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Calculates the elemental contributions
     */
    virtual void CalculateElementalSystem(LocalSystemComponents& rLocalSystem,
                                          ProcessInfo& rCurrentProcessInfo);
    
    /**
     * Calculates the elemental dynamic contributions
     */
    virtual void CalculateDynamicSystem(LocalSystemComponents& rLocalSystem,
					ProcessInfo& rCurrentProcessInfo);

    /**
     * Prints element information for each gauss point (debugging purposes)
     */
    void PrintElementCalculation(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables);


    /**
     * Calculation of the tangent via perturbation of the dofs variables : testing purposes
     */
    void CalculatePerturbedLeftHandSide (MatrixType& rLeftHandSideMatrix, 
					 ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculation and addition of the matrices of the LHS
     */

    virtual void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                                    GeneralVariables& rVariables,
                                    double& rIntegrationWeight);

    /**
     * Calculation and addition of the vectors of the RHS
     */

    virtual void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
                                    GeneralVariables& rVariables,
                                    Vector& rVolumeForce,
                                    double& rIntegrationWeight);



    /**
     * Calculation and addition of the matrices of the LHS
     */

    virtual void CalculateAndAddDynamicLHS(MatrixType& rLeftHandSideMatrix, 
					   GeneralVariables& rVariables, 
					   ProcessInfo& rCurrentProcessInfo, 
					   double& rIntegrationWeight);

    /**
     * Calculation and addition of the vectors of the RHS
     */

    virtual void CalculateAndAddDynamicRHS(VectorType& rRightHandSideVector, 
					   GeneralVariables& rVariables, 
					   ProcessInfo& rCurrentProcessInfo, 
					   double& rIntegrationWeight);



    /**
     * Calculation of the Material Stiffness Matrix. Kuum = BT * C * B
     */

    virtual void CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
                                     GeneralVariables& rVariables,
                                     double& rIntegrationWeight);

    /**
     * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
     */
    virtual void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
                                     GeneralVariables& rVariables,
                                     double& rIntegrationWeight);


    /**
     * Calculation of the External Forces Vector. Fe = N * t + N * b
     */
    virtual void CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
					       GeneralVariables& rVariables,
					       Vector& rVolumeForce,
					       double& rIntegrationWeight);


    /**
      * Calculation of the Internal Forces Vector. Fi = B * sigma
      */
    virtual void CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
					       GeneralVariables & rVariables,
					       double& rIntegrationWeight);


    /**
     * Set Variables of the Element to the Parameters of the Constitutive Law
     */
    virtual void SetGeneralVariables(GeneralVariables& rVariables,
                                     ConstitutiveLaw::Parameters& rValues,
                                     const int & rPointNumber);



    /**
     * Initialize System Matrices
     */
    virtual void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          Flags& rCalculationFlags);



    /**
     * Initialize Material Properties on the Constitutive Law
     */
    void InitializeMaterial ();


    /**
     * Reset the Constitutive Law Parameters
     */
    void ResetConstitutiveLaw();


    /**
     * Clear Nodal Forces
     */
    void ClearNodalForces ();


    /**
     * Calculate Element Kinematics
     */
    virtual void CalculateKinematics(GeneralVariables& rVariables,
                                     const double& rPointNumber);


    /**
     * Calculation of the Deformation Gradient F
     */
    Matrix& CalculateDeltaPosition(Matrix & rDeltaPosition);

    /**
     * Correct Precision Errors (for rigid free movements)
     */
    void DecimalCorrection(Vector& rVector);


    /**
     * Initialize Element General Variables
     */
    virtual void InitializeGeneralVariables(GeneralVariables & rVariables, 
					    const ProcessInfo& rCurrentProcessInfo);

    /**
     * Transform Element General Variables
     */
    virtual void TransformGeneralVariables(GeneralVariables & rVariables, 
					   const double& rPointNumber);

    /**
     * Finalize Element Internal Variables
     */
    virtual void FinalizeStepVariables(GeneralVariables & rVariables, 
				       const double& rPointNumber);


    /**
     * Get the Historical Deformation Gradient to calculate after finalize the step
     */
    virtual void GetHistoricalVariables( GeneralVariables& rVariables, 
					 const double& rPointNumber );


    /**
     * Calculation of the Green Lagrange Strain Vector
     */
    virtual void CalculateGreenLagrangeStrain(const Matrix& rF,
            Vector& rStrainVector);

    /**
     * Calculation of the Almansi Strain Vector
     */
    virtual void CalculateAlmansiStrain(const Matrix& rF,
                                        Vector& rStrainVector);



    /**
     * Calculation of the Velocity Gradient
     */
    void CalculateVelocityGradient(const Matrix& rDN_DX,
                                   Matrix& rDF );


    /**
     * Calculation of the Integration Weight
     */
    virtual double& CalculateIntegrationWeight(double& rIntegrationWeight);


    /**
     * Calculation of the Total Mass of the Element
     */
    virtual double& CalculateTotalMass(double& rTotalMass, const ProcessInfo& rCurrentProcessInfo);


    /**
     * Calculation of the Volume Change of the Element
     */
    virtual double& CalculateVolumeChange(double& rVolumeChange, GeneralVariables& rVariables);

    /**
     * Calculation of the Volume Force of the Element
     */
    virtual Vector& CalculateVolumeForce(Vector& rVolumeForce, GeneralVariables& rVariables);


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

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);


    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class LargeDisplacementElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_LARGE_DISPLACEMENT_ELEMENT_H_INCLUDED  defined 
