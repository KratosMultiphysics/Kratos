//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_FORCE_LOAD_CONDITION_H_INCLUDED )
#define  KRATOS_FORCE_LOAD_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/condition.h"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"

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

/// Force Load Condition for 3D and 2D geometries. (base class)

/**
 * Implements a Force Load definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D (base class)
 */
class KRATOS_API(SOLID_MECHANICS_APPLICATION) ForceLoadCondition
    : public Condition
{
public:

    ///@name Type Definitions
    ///@{
    // Counted pointer of ForceLoadCondition
    KRATOS_CLASS_POINTER_DEFINITION( ForceLoadCondition );
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
     * Parameters to be used in the Condition as they are. 
     */

    struct GeneralVariables
    {
      private:

        //variables including all integration points
        const GeometryType::ShapeFunctionsGradientsType* pDN_De;
        const Matrix* pNcontainer;

      public:

        //for axisymmetric use only
        double  CurrentRadius;
        double  ReferenceRadius;

        //general variables
        double  DomainSize;
        double  Jacobian;
        Vector  N;
        Matrix  DN_De;
        Matrix  DeltaPosition;

        //pressure loads
        double  Pressure;
        Vector  Normal;
        Vector  Tangent1;
        Vector  Tangent2;


        //variables including all integration points
        GeometryType::JacobiansType j;
        GeometryType::JacobiansType J;

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

        void Initialize( const unsigned int& dimension,
			 const unsigned int& local_dimension,
			 const unsigned int& number_of_nodes )
        {
	  //doubles
	  //radius
	  CurrentRadius   = 0;
	  ReferenceRadius = 0;
	  //jacobians
	  DomainSize = 1;
	  Jacobian   = 1;
	  Pressure   = 0;
	  //vectors
	  N.resize(number_of_nodes,false);
 	  Normal.resize(dimension,false);
	  Tangent1.resize(dimension,false);
	  Tangent2.resize(dimension,false);
	  noalias(N) = ZeroVector(number_of_nodes);
 	  noalias(Normal) = ZeroVector(dimension);
	  noalias(Tangent1) = ZeroVector(dimension);
	  noalias(Tangent2) = ZeroVector(dimension);
	  //matrices
	  DN_De.resize(number_of_nodes, local_dimension,false);
	  noalias(DN_De) = ZeroMatrix(number_of_nodes, local_dimension);
	  DeltaPosition.resize(number_of_nodes, dimension,false);
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
     * is defined here and is used to declare a member variable in the component wise condition
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
    ForceLoadCondition();
  
    /// Default constructor.
    ForceLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    ForceLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Copy constructor
    ForceLoadCondition( ForceLoadCondition const& rOther);

    /// Destructor
    virtual ~ForceLoadCondition();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new condition pointer
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(IndexType NewId,
			      NodesArrayType const& ThisNodes,
			      PropertiesType::Pointer pProperties ) const;


    /**
     * clones the selected condition variables, creating a new one
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone(IndexType NewId, 
			     NodesArrayType const& ThisNodes) const;


    //************* STARTING - ENDING  METHODS


    /**
     * Called at the beginning of each solution step
     */
    void Initialize();

    /**
     * Called at the beginning of each solution step
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);

    /**
     * Called at the beginning of each iteration
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);


    //************* GETTING METHODS

    /**
     * Sets on rConditionDofList the degrees of freedom of the considered element geometry
     */
    void GetDofList(DofsVectorType& rConditionDofList,
		    ProcessInfo& rCurrentProcessInfo );

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     */
    void EquationIdVector(EquationIdVectorType& rResult,
			  ProcessInfo& rCurrentProcessInfo );

    /**
     * Sets on rValues the nodal displacements
     */
    void GetValuesVector(Vector& rValues,
			 int Step = 0 );

    /**
     * Sets on rValues the nodal velocities
     */
    void GetFirstDerivativesVector(Vector& rValues,
				   int Step = 0 );

    /**
     * Sets on rValues the nodal accelerations
     */
    void GetSecondDerivativesVector(Vector& rValues,
				    int Step = 0 );


    //************* COMPUTING  METHODS

    /**
     * this is called during the assembling process in order
     * to calculate all condition contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the condition left hand side matrix
     * @param rRightHandSideVector: the condition right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
			      VectorType& rRightHandSideVector,
			      ProcessInfo& rCurrentProcessInfo );


    /**
     * this function provides a more general interface to the condition.
     * it is designed so that rLHSvariables and rRHSvariables are passed TO the condition
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
      * to calculate the condition right hand side vector only
      * @param rRightHandSideVector: the condition right hand side vector
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
				ProcessInfo& rCurrentProcessInfo );


    /**
     * this function provides a more general interface to the condition.
     * it is designed so that rRHSvariables are passed TO the condition
     * thus telling what is the desired output
     * @param rRightHandSideVectors: container for the desired RHS output
     * @param rRHSVariables: parameter describing the expected RHSs
     */
    void CalculateRightHandSide(std::vector< VectorType >& rRightHandSideVectors,
				const std::vector< Variable< VectorType > >& rRHSVariables,
				ProcessInfo& rCurrentProcessInfo);

    /**
     * this is called during the assembling process in order
     * to calculate the condition left hand side matrix only
     * @param rLeftHandSideMatrix: the condition left hand side matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, 
			       ProcessInfo& rCurrentProcessInfo);

    /**
      * this is called during the assembling process in order
      * to calculate the condition mass matrix
      * @param rMassMatrix: the condition mass matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo );

    /**
      * this is called during the assembling process in order
      * to calculate the condition damping matrix
      * @param rDampingMatrix: the condition damping matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfo );


    /**
     * this function is designed to make the element to assemble an rRHS vector
     * identified by a variable rRHSVariable by assembling it to the nodes on the variable
     * rDestinationVariable.
     * @param rRHSVector: input variable containing the RHS vector to be assembled
     * @param rRHSVariable: variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable: variable in the database to which the rRHSvector will be assembled 
      * @param rCurrentProcessInfo: the current process info instance
     */      
    virtual void AddExplicitContribution(const VectorType& rRHS, 
					 const Variable<VectorType>& rRHSVariable, 
					 Variable<array_1d<double,3> >& rDestinationVariable, 
					 const ProcessInfo& rCurrentProcessInfo);


    /**
     * Get on rVariable a double Value
     */
    void GetValueOnIntegrationPoints( const Variable<double>& rVariable, 
				      std::vector<double>& rValues, 
				      const ProcessInfo& rCurrentProcessInfo );

    /**
     * Calculate a double Variable
     */
    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, 
				      std::vector<double>& rOutput, 
				      const ProcessInfo& rCurrentProcessInfo);



    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    virtual int Check( const ProcessInfo& rCurrentProcessInfo );

    ///@}
    ///@name Access
    ///@{
    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{
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
     * Energy variable for loads
     */
    double mEnergy; 

    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Clear Nodal Forces
     */
    void ClearNodalForces ();


    /**
     * Initialize System Matrices
     */

    virtual void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
					  VectorType& rRightHandSideVector,
					  Flags& rCalculationFlags);

    /**
     * Initialize General Variables
     */
    virtual void InitializeGeneralVariables(GeneralVariables& rVariables, 
					    const ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculate Condition Kinematics
     */
    virtual void CalculateKinematics(GeneralVariables& rVariables, 
				     const double& rPointNumber);

    /**
     * Calculation of the Vector Force of the Condition
     */
    virtual Vector& CalculateVectorForce(Vector& rVectorForce, GeneralVariables& rVariables);


    /**
     * Calculates the condition contributions
     */
    virtual void CalculateConditionSystem(LocalSystemComponents& rLocalSystem,
					  const ProcessInfo& rCurrentProcessInfo);


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
     * Calculation of the Load Stiffness Matrix which usually is subtracted to the global stiffness matrix
     */
    virtual void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
				     GeneralVariables& rVariables,
				     double& rIntegrationWeight);


    /**
     * Calculation of the External Forces Vector for a force or pressure vector 
     */
    virtual void CalculateAndAddExternalForces(Vector& rRightHandSideVector,
					       GeneralVariables& rVariables,
					       Vector& rVectorLoad,
					       double& rIntegrationWeight );


    /**
     * Get Node Movements for energy computation
     */
    void GetNodalDeltaMovements(Vector& rValues, const int& rNode);


    /**
     * Get Current Value, buffer 0 with FastGetSolutionStepValue
     */    
    Vector& GetCurrentValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const unsigned int& rNode);

   /**
     * Get Previous Value, buffer 1 with FastGetSolutionStepValue
     */    
    Vector& GetPreviousValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const unsigned int& rNode);


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
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);


}; // class ForceLoadCondition.

} // namespace Kratos.

#endif // KRATOS_FORCE_LOAD_CONDITION_H_INCLUDED defined 
