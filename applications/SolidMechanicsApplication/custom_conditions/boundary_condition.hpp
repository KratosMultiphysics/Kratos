//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_BOUNDARY_CONDITION_H_INCLUDED)
#define  KRATOS_BOUNDARY_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/condition.h"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"
#include "utilities/beam_math_utilities.hpp"

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

/// General Boundary Condition base type for 3D and 2D geometries.

/**
 * Implements a General definitions for a boundary neumann or mixed condition.
 * This works for arbitrary geometries in 3D and 2D (base class)
 */
class KRATOS_API(SOLID_MECHANICS_APPLICATION) BoundaryCondition
    : public Condition
{
public:

    ///@name Type Definitions
    ///@{

    typedef Variable<array_1d<double,3>>      VariableVectorType;
    typedef Variable<double>                  VariableScalarType;

    ///Type for size
    typedef GeometryData::SizeType                      SizeType;

    // Counted pointer of BoundaryCondition
    KRATOS_CLASS_POINTER_DEFINITION( BoundaryCondition );
    ///@}

protected:

    /**
     * Flags related to the element computation
     */

    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX );

    /**
     * Parameters to be used in the Condition as they are.
     */

    struct ConditionVariables
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
        double  GeometrySize;
        double  Jacobian;
        Vector  N;
        Matrix  DN_De;
        Matrix  DeltaPosition;

        //external boundary values
        double  ExternalScalarValue;
        Vector  ExternalVectorValue;

        //boundary characteristics
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
	  GeometrySize = 1;
	  Jacobian     = 1;
	  //external boundary values
	  ExternalScalarValue  = 0;
	  ExternalVectorValue.resize(dimension,false);
	  noalias(ExternalVectorValue) = ZeroVector(dimension);
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
	  J[0].resize(dimension,dimension,false);
	  j[0].resize(dimension,dimension,false);
	  noalias(J[0]) = ZeroMatrix(dimension,dimension);
	  noalias(j[0]) = ZeroMatrix(dimension,dimension);
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

    public:

      //calculation flags
      Flags  CalculationFlags;

      /**
       * sets the value of a specified pointer variable
       */
      void SetLeftHandSideMatrix( MatrixType& rLeftHandSideMatrix ) { mpLeftHandSideMatrix = &rLeftHandSideMatrix; };

      void SetRightHandSideVector( VectorType& rRightHandSideVector ) { mpRightHandSideVector = &rRightHandSideVector; };

      /**
       * returns the value of a specified pointer variable
       */
      MatrixType& GetLeftHandSideMatrix() { return *mpLeftHandSideMatrix; };

      VectorType& GetRightHandSideVector() { return *mpRightHandSideVector; };

    };


public:


    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    BoundaryCondition();

    /// Default constructor.
    BoundaryCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    BoundaryCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Copy constructor
    BoundaryCondition( BoundaryCondition const& rOther);

    /// Destructor
    ~BoundaryCondition() override;

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
			      PropertiesType::Pointer pProperties ) const override;


    /**
     * clones the selected condition variables, creating a new one
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone(IndexType NewId,
			     NodesArrayType const& ThisNodes) const override;


    //************* STARTING - ENDING  METHODS


    /**
     * Called at the beginning of each solution step
     */
    void Initialize() override;

    /**
     * Called at the beginning of each solution step
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Called at the beginning of each iteration
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;


    //************* GETTING METHODS

    /**
     * Sets on rConditionDofList the degrees of freedom of the considered element geometry
     */
    void GetDofList(DofsVectorType& rConditionDofList,
		    ProcessInfo& rCurrentProcessInfo ) override;

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     */
    void EquationIdVector(EquationIdVectorType& rResult,
			  ProcessInfo& rCurrentProcessInfo ) override;

    /**
     * Sets on rValues the nodal displacements
     */
    void GetValuesVector(Vector& rValues,
			 int Step = 0 ) override;

    /**
     * Sets on rValues the nodal velocities
     */
    void GetFirstDerivativesVector(Vector& rValues,
				   int Step = 0 ) override;

    /**
     * Sets on rValues the nodal accelerations
     */
    void GetSecondDerivativesVector(Vector& rValues,
				    int Step = 0 ) override;


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
			      ProcessInfo& rCurrentProcessInfo ) override;

    /**
      * this is called during the assembling process in order
      * to calculate the condition right hand side vector only
      * @param rRightHandSideVector: the condition right hand side vector
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
				ProcessInfo& rCurrentProcessInfo ) override;

    /**
     * this is called during the assembling process in order
     * to calculate the condition left hand side matrix only
     * @param rLeftHandSideMatrix: the condition left hand side matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
			       ProcessInfo& rCurrentProcessInfo) override;

    /**
      * this is called during the assembling process in order
      * to calculate the condition mass matrix
      * @param rMassMatrix: the condition mass matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateMassMatrix(MatrixType& rMassMatrix,
			     ProcessInfo& rCurrentProcessInfo ) override;

    /**
      * this is called during the assembling process in order
      * to calculate the condition damping matrix
      * @param rDampingMatrix: the condition damping matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateDampingMatrix(MatrixType& rDampingMatrix,
				ProcessInfo& rCurrentProcessInfo ) override;


    /**
     * this function is designed to make the element to assemble an rRHS vector
     * identified by a variable rRHSVariable by assembling it to the nodes on the variable
     * rDestinationVariable.
     * @param rRHSVector: input variable containing the RHS vector to be assembled
     * @param rRHSVariable: variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable: variable in the database to which the rRHSvector will be assembled
      * @param rCurrentProcessInfo: the current process info instance
     */
    void AddExplicitContribution(const VectorType& rRHS,
					 const Variable<VectorType>& rRHSVariable,
					 Variable<array_1d<double,3> >& rDestinationVariable,
					 const ProcessInfo& rCurrentProcessInfo) override;


    /**
     * Get on rVariable a double Value
     */
    void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
				     std::vector<double>& rValues,
				     const ProcessInfo& rCurrentProcessInfo ) override;

    /**
     * Calculate a double Variable
     */
    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
				      std::vector<double>& rOutput,
				      const ProcessInfo& rCurrentProcessInfo) override;



    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check( const ProcessInfo& rCurrentProcessInfo ) override;

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


    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Initialize Explicit Contributions
     */
    void InitializeExplicitContributions();


    /**
     * Check dof for a vector variable
     */
    virtual bool HasVariableDof(VariableVectorType& rVariable);

    /**
     * Check dof for a double variable
     */
    virtual bool HasVariableDof(VariableScalarType& rVariable);


    /**
     * Get condition size from the dofs
     */
    virtual unsigned int GetDofsSize();

    /**
     * Initialize System Matrices
     */
    virtual void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
					  VectorType& rRightHandSideVector,
					  Flags& rCalculationFlags);

    /**
     * Initialize General Variables
     */
    virtual void InitializeConditionVariables(ConditionVariables& rVariables,
					    const ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculate Condition Kinematics
     */
    virtual void CalculateKinematics(ConditionVariables& rVariables,
				     const double& rPointNumber);


    /**
     * Calculates the condition contributions
     */
    virtual void CalculateConditionSystem(LocalSystemComponents& rLocalSystem,
					  const ProcessInfo& rCurrentProcessInfo);


    /**
     * Calculation and addition of the matrices of the LHS
     */
    virtual void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                                    ConditionVariables& rVariables,
                                    double& rIntegrationWeight);

    /**
     * Calculation and addition of the vectors of the RHS
     */
    virtual void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
                                    ConditionVariables& rVariables,
                                    double& rIntegrationWeight);


    /**
     * Calculation of the Load Stiffness Matrix which usually is subtracted to the global stiffness matrix
     */
    virtual void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
				     ConditionVariables& rVariables,
				     double& rIntegrationWeight);


    /**
     * Calculation of the External Forces Vector for a force or pressure vector
     */
    virtual void CalculateAndAddExternalForces(Vector& rRightHandSideVector,
					       ConditionVariables& rVariables,
					       double& rIntegrationWeight);

    /**
     * Calculation of the External Forces Vector for a force or pressure vector
     */
    virtual double& CalculateAndAddExternalEnergy(double& rEnergy,
						  ConditionVariables& rVariables,
						  double& rIntegrationWeight,
						  const ProcessInfo& rCurrentProcessInfo);
    /**
     * Get Node Movements for energy computation
     */
    void GetNodalDeltaMovements(Vector& rValues, const int& rNode);

    /**
     * Calculation of the Position Increment
     */
    virtual Matrix& CalculateDeltaPosition(Matrix & rDeltaPosition);

    /**
     * Calculation of the Total Position Increment
     */
    virtual Matrix& CalculateTotalDeltaPosition(Matrix & rDeltaPosition);

    /**
     * Get Current Value, buffer 0 with FastGetSolutionStepValue
     */
    Vector& GetNodalCurrentValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const unsigned int& rNode);

   /**
     * Get Previous Value, buffer 1 with FastGetSolutionStepValue
     */
    Vector& GetNodalPreviousValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const unsigned int& rNode);


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

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;


}; // class BoundaryCondition.

} // namespace Kratos.

#endif // KRATOS_BOUNDARY_CONDITION_H_INCLUDED defined
