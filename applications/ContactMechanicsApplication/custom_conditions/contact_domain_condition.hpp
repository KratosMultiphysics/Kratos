//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_CONTACT_DOMAIN_CONDITION_H_INCLUDED )
#define  KRATOS_CONTACT_DOMAIN_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"
#include "includes/checks.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/kratos_flags.h"

#include "custom_utilities/contact_domain_utilities.hpp"
#include "custom_utilities/mesh_data_transfer_utilities.hpp"

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


class KRATOS_API(CONTACT_MECHANICS_APPLICATION) ContactDomainCondition
    : public Condition
{
public:


    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    ///NodeType
    typedef Node < 3 > NodeType;
    ///Geometry Type
    typedef Geometry<NodeType> GeometryType;
    ///Element Type
    typedef Element::ElementType ElementType;


    ///Tensor order 1 definition
    typedef ContactDomainUtilities::PointType             PointType;
    ///SurfaceVector
    typedef ContactDomainUtilities::SurfaceVector     SurfaceVector;
    ///SurfaceScalar
    typedef ContactDomainUtilities::SurfaceScalar     SurfaceScalar;
    ///BaseLengths
    typedef ContactDomainUtilities::BaseLengths         BaseLengths;

    ///For 3D contact surfaces definition
    typedef ContactDomainUtilities::ScalarBaseType  ScalarBaseType;
    typedef ContactDomainUtilities::SurfaceBase           SurfaceBase;

protected:

    /**
     * Parameters to be used in the Condition as they are. Direct interface to Parameters Struct
     */

    typedef struct
    {
      //Geometrical surface tangent gaps:
      ScalarBaseType   CurrentGap;     //tangential gap

      //Contact constraint parameters
      double   Multiplier;            //Lagrange Multipliyer tangent
      double   Penalty;               //Penalty Parameter tangent

      //Variables of the contact domain elements
      Vector          dN_dt;      //Discrete variacion of the shape function  in the current tangent direction
      std::vector<Vector >       Tsigma;

      ScalarBaseType   CurrentTensil;
      double           GapSign;

    } ContactSurfaceParameters;



    typedef struct
    {
      ContactSurfaceParameters A;
      ContactSurfaceParameters B;

      SurfaceBase    CovariantBase;
      SurfaceBase    ContravariantBase;

      //geometrical variables
      double ReferenceArea;
      double CurrentArea;

      double FactorArea;

      double EquivalentHeigh;

      double ElementSize;

    } ContactTangentParameters;


    typedef struct
    {
        Flags           Options;               //calculation options

        //Geometrical gaps:
        SurfaceScalar   CurrentGap;             //normal and tangential gap

        //The stabilization factor or penalty factor
        SurfaceScalar   ContactFactor;

        //Friction:
        double          FrictionCoefficient;   //total friction coeffitient mu
	double          TangentialGapSign;     //sign or direction of the tangential gap

        SurfaceScalar   Multiplier;            //Lagrange Multipliyer normal and tangent
        SurfaceScalar   Penalty;               //Penalty Parameter normal and tangent


        //variables of the contact element 2D
        Vector          dN_dn;      //Discrete variacion of the shape function  in the current normal direction
        Vector          dN_dt;      //Discrete variacion of the shape function  in the current tangent direction
        Vector          dN_drn;     //Discrete variacion of the shape function  in the reference normal direction

        std::vector<Vector >       Nsigma;
        std::vector<Vector >       Tsigma;

        //Geometric variables
        SurfaceVector        CurrentSurface;

        std::vector<BaseLengths>   CurrentBase;    //Current Base Lengths variables
        std::vector<BaseLengths>   ReferenceBase;  //Reference Base Lengths variables

        //Resultant mechanical tractions
        SurfaceScalar              CurrentTensil;  //Tangential and Normal modulus of the traction vector components


        //Tangent parameters used in 3D
        ContactTangentParameters Tangent;


    } ContactParameters;



    typedef struct
    {
      ScalarBaseType   PreviousGapA;
      ScalarBaseType   PreviousGapB;

      SurfaceBase    CovariantBase;
      SurfaceBase    ContravariantBase;

    } ContactTangentVariables;


    typedef struct
    {
        double  detF;
        double  detJ;
        Vector  StrainVector;
        Vector  StressVector;
        Vector  N;
        Matrix  F;
        Matrix  DN_DX;
        Matrix  ConstitutiveMatrix;

	ContactParameters Contact;

        //Axisymmetric
        double  CurrentRadius;
        double  ReferenceRadius;

	//variables including all integration points
        const GeometryType::ShapeFunctionsGradientsType* pDN_De;
	const Matrix* pNcontainer;

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

    } ConditionVariables;


    typedef struct
    {

        //Iteration counter:
        int             IterationCounter;     //the number of the step iteration

        //The stabilization parameter and penalty parameter
        double          StabilizationFactor;
	double          PenaltyFactor;

        //Geometrical gaps:
        SurfaceScalar        PreviousGap;     //effective normal and tangential gap in previous time step configuration

	//Geometric variables
        SurfaceVector        PreStepSurface;
        SurfaceVector        ReferenceSurface;

        //Tangent variables in 3D
        ContactTangentVariables  Tangent;

        //Contact condition conectivities
        std::vector<unsigned int> nodes;
	std::vector<unsigned int> order;
        std::vector<unsigned int> slaves;

        //Resultant mechanical tractions
	PointType       TractionVector;       //Traction Vector in the reference configuration

	//Pointer Variables
        GeometryType*         mpMasterGeometry;
        ElementType*          mpMasterElement;
        Condition*            mpMasterCondition;
        NodeType*             mpMasterNode;

	/**
         * sets the value of a specified pointer variable
	 */

	void SetMasterGeometry  (GeometryType& rGeometry){ mpMasterGeometry = &rGeometry; }
        void SetMasterElement   (ElementType& rElement){ mpMasterElement = &rElement; }
        void SetMasterCondition (ConditionType& rCondition){ mpMasterCondition = &rCondition; }
        void SetMasterNode      (NodeType& rNode){ mpMasterNode = &rNode; }

	/**
         * returns the value of a specified pointer variable
         */

        GeometryType& GetMasterGeometry()   { return (*mpMasterGeometry); }
	ElementType& GetMasterElement()     { return (*mpMasterElement); }
	ConditionType& GetMasterCondition() { return (*mpMasterCondition); }
	NodeType& GetMasterNode()           { return (*mpMasterNode); }


    } ContactVariables;




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

    /// Counted pointer of ContactDomainCondition
    KRATOS_CLASS_POINTER_DEFINITION( ContactDomainCondition );


    ///@}
    ///@name Life Cycle
    ///@{


    /// Default constructors.
    ContactDomainCondition() : Condition() {};

    ContactDomainCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    ContactDomainCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    ContactDomainCondition(ContactDomainCondition const& rOther);


    /// Destructor.
    virtual ~ContactDomainCondition();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ContactDomainCondition& operator=(ContactDomainCondition const& rOther);


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
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    /**
     * clones the selected condition variables, creating a new one
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone(IndexType NewId,
			     NodesArrayType const& ThisNodes) const override;

    //************* GETTING METHODS

    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
    IntegrationMethod GetIntegrationMethod() override;

    /**
     * Sets on rConditionalDofList the degrees of freedom of the considered element geometry
     */
    void GetDofList(DofsVectorType& rConditionalDofList, ProcessInfo& rCurrentProcessInfo) override;

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
     * Set a double Value on the Condition Constitutive Law
     */
    void SetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;
    /**
     * Set a Vector Value on the Condition Constitutive Law
     */
    void SetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Set a Matrix Value on the Condition Constitutive Law
     */
    void SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    //GET:
    /**
     * Set on rVariable a double Value from the Condition Constitutive Law
     */
    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Set on rVariable a Vector Value from the Condition Constitutive Law
     */
    void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Set on rVariable a Matrix Value from the Condition Constitutive Law
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
    void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo) override;

    /**
     * this is called for non-linear analysis at the beginning of the iteration process
     */
    void InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo) override;


    /**
     * this is called for non-linear analysis at the end of the iteration process
     */
    void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo) override;


    /**
     * Called at the end of eahc solution step
     */
    void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo) override;



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
			      ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;


   /**
     * this function provides a more general interface to the condition.
     * it is designed so that rRHSvariables are passed TO the condition
     * thus telling what is the desired output
     * @param rRightHandSideVectors: container for the desired RHS output
     * @param rRHSVariables: parameter describing the expected RHSs
     */
    void CalculateRightHandSide(std::vector< VectorType >& rRightHandSideVectors,
				const std::vector< Variable< VectorType > >& rRHSVariables,
				ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side vector only
     * @param rLeftHandSideVector: the elemental left hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLeftHandSide (MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) override;


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
					 const ProcessInfo& rCurrentProcessInfo) override;


    //on integration points:
    /**
     * Calculate a double Variable on the Condition Constitutive Law
     */
    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculate a Vector Variable on the Condition Constitutive Law
     */
    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculate a Matrix Variable on the Condition Constitutive Law
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

    /**
     * Container for constitutive law instances on each integration point
     */
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    /**
     * Variables stored in the element during the computation
     */
    ContactVariables      mContactVariables;

    /**
     * Contact Domain Utilities
     */
    ContactDomainUtilities  mContactUtilities;

    ///@}
    ///@name Protected Operators
    ///@{

    /**
     * Calculation of the Contact Master Nodes and Mechanical variables
     */
    virtual void SetMasterGeometry()
	{
		KRATOS_THROW_ERROR( std::invalid_argument, "Calling base class in contact domain", "" )

	};


    /**
     * Calculate Tau stabilization or Penalty factor
     */
    virtual void CalculateContactFactor(ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_THROW_ERROR( std::invalid_argument, "Calling base class in contact domain", "" )

	};

    /**
     * Calculation of the Contact Previous Gap
     */
    virtual void CalculatePreviousGap()
	{
		KRATOS_THROW_ERROR( std::invalid_argument, "Calling base class in contact domain", "" )

	};


    /**
     * Initialize Variables
     */
    virtual void InitializeConditionVariables (ConditionVariables& rVariables,
					     const ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculates the condition contributions
     */
    virtual void CalculateConditionSystem(LocalSystemComponents& rLocalSystem,
					  ProcessInfo& rCurrentProcessInfo);
    /**
     * Clear Nodal Forces
     */
    void ClearNodalForces ();


    /**
     * Calculate Nodal Forces
     */
    void CalculateNodalForces (ProcessInfo& CurrentProcessInfo);


    /**
     * Clear Nodal Forces on Master Element
     */
    void ClearMasterElementNodalForces(ElementType& rMasterElement);


    /**
     * Set Master element information on integration points to Contact element information
     */
    void SetContactIntegrationVariable (Vector & rContactVariable,
					std::vector<Vector> & rMasterVariables,
					const unsigned int& rPointNumber);

    /**
     * Calculate Condition Kinematics
     */
    virtual void CalculateKinematics(ConditionVariables& rVariables,
				     ProcessInfo& rCurrentProcessInfo,
				     const unsigned int& rPointNumber);

    /**
     * Calculation of the Deformation Gradient F
     */
    Matrix& CalculateDeltaPosition(Matrix & rDeltaPosition);

    /**
     * Calculation of the Contact Multipliers or Penalty Factors
     */
    virtual void CalculateExplicitFactors(ConditionVariables& rVariables,
					  ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_THROW_ERROR( std::invalid_argument, "Calling base class in contact domain", "" )

	};
    /**
     * Tangent Matrix construction methods:
     */
    virtual void CalculateDomainShapeN(ConditionVariables& rVariables)
	{
		KRATOS_THROW_ERROR( std::invalid_argument, "Calling base class in contact domain", "" )

	};



    /**
     *  Parameters for friction law Relative Tangent Velocity:
     */
    virtual void CalculateRelativeVelocity(ConditionVariables& rVariables,
					   PointType & TangentVelocity,
					   ProcessInfo& rCurrentProcessInfo);

    /**
     *  Parameters for friction law Relative Tangent Displacement:
     */
    virtual void CalculateRelativeDisplacement(ConditionVariables& rVariables,
					       PointType & TangentDisplacement,
					       ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculate current tangent vector
     */
    virtual PointType & CalculateCurrentTangent(PointType &rTangent)
	{
		KRATOS_THROW_ERROR( std::invalid_argument, "Calling base class in contact domain", "" )

	};


    /**
     * Friction Parameters:
     */
    virtual void CalculateFrictionCoefficient(ConditionVariables& rVariables,
					      const PointType & TangentVelocity);


    /**
     * Calculate Integration Weight:
     */
    virtual double& CalculateIntegrationWeight(double& rIntegrationWeight);


    /**
     * Initialize System Matrices
     */
    virtual void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
					  VectorType& rRightHandSideVector,
					  Flags& rCalculationFlags);


    /**
     * Calculation of the tangent via perturbation of the dofs variables : testing purposes
     */
    void CalculatePerturbedLeftHandSide (MatrixType& rLeftHandSideMatrix,
					 ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculate LHS
     */
    virtual void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
				    ConditionVariables& rVariables,
				    double& rIntegrationWeight);

    /**
     * Calculate RHS
     */
    virtual void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
				    ConditionVariables& rVariables,
				    double& rIntegrationWeight);


    /**
     * Calculation of the Contact Stiffness Matrix
     */
    virtual void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
				     ConditionVariables& rVariables,
				     double& rIntegrationWeight);

    /**
     * Calculation of the Material Stiffness Matrix by components
     */
    virtual void CalculateContactStiffness (double &Kcont,ConditionVariables& rVariables,
					    unsigned int& ndi,unsigned int& ndj,
					    unsigned int& idir,unsigned int& jdir)
	{
		KRATOS_THROW_ERROR( std::invalid_argument, "Calling base class in contact domain", "" )

	};

    /**
     * Calculation of the Internal Forces Vector. Fi = B * sigma
     */
    virtual void CalculateAndAddContactForces(VectorType& rRightHandSideVector,
					      ConditionVariables& rVariables,
					      double& rIntegrationWeight);



    /**
     * Normal Force construction by components
     */
    virtual void CalculateNormalForce       (double &F,ConditionVariables& rVariables,
					     unsigned int& ndi,unsigned int& idir)
	{
		KRATOS_THROW_ERROR( std::invalid_argument, "Calling base class in contact domain", "" )

	};
    /**
     * Tangent Stick Force construction by components
     */
    virtual void CalculateTangentStickForce (double &F,ConditionVariables& rVariables,
					     unsigned int& ndi,unsigned int& idir)
	{
		KRATOS_THROW_ERROR( std::invalid_argument, "Calling base class in contact domain", "" )

	};
    /**
     * Tangent Slip Force construction by components
     */
    virtual void CalculateTangentSlipForce  (double &F,ConditionVariables& rVariables,
					     unsigned int& ndi,unsigned int& idir)
	{
		KRATOS_THROW_ERROR( std::invalid_argument, "Calling base class in contact domain", "" )

	};


    ///@}
    ///@name Protected Operations
    ///@{
    ///@}
    ///@name Protected  Access
    ///@{
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Protected Inquiry
    ///@{
    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:

    ///@name Private static Member Variables
    ///@{
    ///@}
    ///@name Private member Variables
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
    ///@name Serialization
    ///@{
    ///@}
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class ContactDomainCondition

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
    ContactDomainCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
    const ContactDomainCondition& rThis)
    {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
    }*/
///@}

} // namespace Kratos.
#endif // KRATOS_CONTACT_DOMAIN_CONDITION_H_INCLUDED  defined
