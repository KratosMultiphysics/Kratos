//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_RIGID_BODY_ELEMENT_H_INCLUDED )
#define  KRATOS_RIGID_BODY_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"
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

/// Rigid Body Element for 3D space dimension

/**
 * Nodal Variables: DISPLACEMENT, STEP_DISPLACEMENT, VELOCITY, ACCELERATION, ROTATION, STEP_ROTATION, DELTA_ROTATION, ANGULAR_VELOCITY, ANGULAR_ACCELERATION
 * Nodal Dofs: DISPLACEMENT, ROTATION
 */

class KRATOS_API(CONTACT_MECHANICS_APPLICATION) RigidBodyElement
    :public Element
{
public:

    ///@name Type Definitions
    ///@{
    ///Type definition for beam utilities
    typedef BeamMathUtils<double>                      BeamMathUtilsType;
    ///Type definition for quaternion
    typedef Quaternion<double>                            QuaternionType;
    ///Type for nodes
    typedef Node<3>                                             NodeType;
    ///Type for nodes container
    typedef PointerVectorSet<NodeType, IndexedObject> NodesContainerType;
    ///Type for size
    typedef GeometryData::SizeType                              SizeType;

    /// Counted pointer of RigidBodyElement
    KRATOS_CLASS_POINTER_DEFINITION( RigidBodyElement );

    ///@}

protected:

    /**
     * Flags related to the element computation
     */
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX );

    /**
     * Parameters to be used to store section properties
     */
    struct RigidBodyProperties
    {
      double Mass;                            // Mass of the Rigid Body
      Matrix InertiaTensor;                   // Global Inertia Tensor
    };

    /**
     * Parameters to be used in the Element as they are. Direct interface to Parameters Struct
     */
    struct ElementVariables
    {
     private:

      //variables including all integration points
      const ProcessInfo* pProcessInfo;

     public:
      
      //section properties
      RigidBodyProperties RigidBody;
      Vector VolumeForce;
      Matrix DeltaPosition;
      
      void SetProcessInfo(const ProcessInfo& rProcessInfo)
      {
        pProcessInfo=&rProcessInfo;
      }

      const ProcessInfo& GetProcessInfo()
      {
        return *pProcessInfo;
      }
      
      void Initialize(const unsigned int& dimension, const ProcessInfo& rProcessInfo)
      {
        VolumeForce.resize(dimension);
        noalias(VolumeForce) = ZeroVector(dimension);
        DeltaPosition.resize(1,dimension,false);
        noalias(DeltaPosition) = ZeroMatrix(1, dimension);
        pProcessInfo=&rProcessInfo;
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

    public:

      //calculation flags
      Flags        CalculationFlags;

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

    /// Serializer constructor
    RigidBodyElement() {};

    /// Default constructors
    RigidBodyElement(IndexType NewId, GeometryType::Pointer pGeometry);

    RigidBodyElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    RigidBodyElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties, NodesContainerType::Pointer pNodes);

    ///Copy constructor
    RigidBodyElement(RigidBodyElement const& rOther);

    /// Destructor.
    virtual ~RigidBodyElement();


    ///@}
    ///@name Operators
    ///@{


   /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;



    //************* GETTING METHODS


    /**
     * Sets on rElementalDofList the degrees of freedom of the considered element geometry
     */
    void GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo) override;

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
     * this is called for non-linear analysis at the beginning of the iteration process
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called for non-linear analysis at the beginning of the iteration process
     */
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;


    /**
     * Called at the end of eahc solution step
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
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) override;


    /**
     * this is called during the assembling process in order
     * to calculate the second derivatives contributions for the LHS and RHS
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix,
						VectorType& rRightHandSideVector,
						ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix for the second derivatives constributions
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
				       ProcessInfo& rCurrentProcessInfo) override;


    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector for the second derivatives constributions
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
				       ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental mass matrix
     * @param rMassMatrix: the elemental mass matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override;


    /**
     * this function is designed to make the element to assemble an rRHS vector
     * identified by a variable rRHSVariable by assembling it to the nodes on the variable
     * rDestinationVariable.
     * The "AddEXplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH AN ELEMENT
     * IS ALLOWED TO WRITE ON ITS NODES.
     * the caller is expected to ensure thread safety hence
     * SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
     * @param rRHSVector: input variable containing the RHS vector to be assembled
     * @param rRHSVariable: variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable: variable in the database to which the rRHSVector will be assembled
      * @param rCurrentProcessInfo: the current process info instance
     */
    void AddExplicitContribution(const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable, Variable<array_1d<double,3> >& rDestinationVariable, const ProcessInfo& rCurrentProcessInfo) override;


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
        buffer << "Rigid Body Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Rigid Body Element #" << Id();
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

    /**
     * Global to Local Quaternion for Global to Local tensor transformation SPATIAL
     */
    QuaternionType  mInitialLocalQuaternion;

    NodesContainerType::Pointer     mpNodes;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    
    /**
     * Calculates the elemental dynamic contributions
     */
    void CalculateDynamicSystem(LocalSystemComponents& rLocalSystem,
                                ProcessInfo& rCurrentProcessInfo);

    /**
     * Initialize System Matrices
     */
    virtual void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          Flags& rCalculationFlags);

    /**
     * Transform Vector Variable from Global Frame to the Spatial Local Frame
     */
    Vector& MapToInitialLocalFrame(Vector& rVariable);


    /**
     * Get Current Value, buffer 0 with FastGetSolutionStepValue
     */
    Vector& GetNodalCurrentValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const unsigned int& rNode);

    /**
     * Get Previous Value, buffer 1 with FastGetSolutionStepValue
     */
    Vector& GetNodalPreviousValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const unsigned int& rNode);


    /**
     * Calculation of the Rigid Body Properties
     */
    void CalculateRigidBodyProperties(RigidBodyProperties & rRigidBody);


    /**
      * Calculation of the Tangent Matrix
      */
    virtual void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                                    ElementVariables& rVariables);

    /**
      * Calculation of the Force Vector
      */
    virtual void CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                    ElementVariables& rVariables);
    
    /**
     * Calculation of the External Forces Vector. Fe = N * t + N * b
     */
    virtual void CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
					       ElementVariables& rVariables);

    /**
      * Calculation of the Tangent Intertia Matrix
      */
    virtual void CalculateAndAddInertiaLHS(MatrixType& rLeftHandSideMatrix,
					   ElementVariables& rVariables);

    /**
      * Calculation of the Inertial Forces Vector
      */
    virtual void CalculateAndAddInertiaRHS(VectorType& rRightHandSideVector,
					   ElementVariables& rVariables);


    /**
      * Calculation of the time integration parameters
      */
    virtual void GetTimeIntegrationParameters(double& rP0,double& rP1,double& rP2,
                                              const ProcessInfo& rCurrentProcessInfo);
    
    /**
     * Calculation of the Volume Force of the Element
     */
    virtual Vector& CalculateVolumeForce(Vector& rVolumeForce);


    /**
     * Calculation Complementary Method : Inertial Matrix Calculation Part 1
     */
    virtual void CalculateRotationLinearPartTensor(Vector& rRotationVector, Matrix& rRotationTensor);


    /**
      * Update rigid body nodes and positions
      */
    virtual void UpdateRigidBodyNodes(ProcessInfo& rCurrentProcessInfo);

    /**
     * Get element size from the dofs
     */
    virtual SizeType GetDofsSize();


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

    virtual void save(Serializer& rSerializer) const override;

    virtual void load(Serializer& rSerializer) override;

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}


}; // Class RigidBodyElement

} // namespace Kratos.
#endif //  KRATOS_RIGID_BODY_ELEMENT_H_INCLUDED defined
