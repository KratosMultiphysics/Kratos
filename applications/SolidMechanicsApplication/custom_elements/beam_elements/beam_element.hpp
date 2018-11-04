//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_BEAM_ELEMENT_H_INCLUDED)
#define  KRATOS_BEAM_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/element.h"
#include "utilities/beam_math_utilities.hpp"
#include "custom_utilities/element_utilities.hpp"

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

/// Beam Element for 2D and 3D space dimensions

class KRATOS_API(SOLID_MECHANICS_APPLICATION) BeamElement
    :public Element
{
public:

    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw                         ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer     ConstitutiveLawPointerType;
    ///StressMeasure from constitutive laws
    typedef ConstitutiveLawType::StressMeasure        StressMeasureType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod           IntegrationMethod;
    ///Type definition for beam utilities
    typedef BeamMathUtils<double>                     BeamMathUtilsType;
    ///Type definition for quaternion
    typedef Quaternion<double>                           QuaternionType;
    ///Type for size
    typedef GeometryData::SizeType                             SizeType;

    /// Counted pointer of BeamElement
    KRATOS_CLASS_POINTER_DEFINITION( BeamElement );

    ///@}

protected:

    /**
     * Flags related to the element computation
     */
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX );
    KRATOS_DEFINE_LOCAL_FLAG( FINALIZED_STEP );
    /**
     * Parameters to be used to store section properties
     */

    struct SectionProperties
    {
      double Area;                            // Area or the beam section
      double Inertia_z;                       // Moment of Inertia about the local z axis, Iz local
      double Inertia_y;                       // Moment of Inertia about the local y axis, Iy local
      double Polar_Inertia;                   // Polar Moment of Inertia, measures an object's ability to resist twisting, when acted upon by differences of torque along its length.
      double Rotational_Inertia;              // Moment of Inertia about the local x axis, measures an object's resistance to changes in its rotational velocity when acted by a net resultant torque
    };

    /**
     * This is used to compute directions for this element
     */
    struct DirectorsVariables
    {
      //Current
      std::vector<Vector> Current;
      std::vector<Vector> CurrentDerivatives;

      //Initial
      std::vector<Vector> Initial;
      std::vector<Vector> InitialDerivatives;

      //Previous
      std::vector<Vector> Previous;
      std::vector<Vector> PreviousDerivatives;


      //NODAL VARIABLES
      std::vector<Matrix> InitialNode;
      std::vector<Matrix> CurrentNode;
      std::vector<Matrix> PreviousNode;

      std::vector<Matrix> CurrentNodeVelocities;
      std::vector<Matrix> PreviousNodeVelocities;
    };


    /**
     * Parameters to be used in the Element as they are. Direct interface to Parameters Struct
     */

    struct ElementData
    {
      private:

        //variables including all integration points
        const GeometryType::ShapeFunctionsGradientsType* pDN_De;
        const Matrix* pNcontainer;

        //directors variables
        DirectorsVariables* Directors;

      public:

        //integration point
        unsigned int PointNumber;

        //delta time
        double DeltaTime;

        //element length
        double  Length;
        double  detJ;

        //equilibrium point
        double Alpha;

        //general variables
        Vector  StrainVector;
        Vector  StressVector;
        Vector  N;

        Matrix  B;
        Matrix  DN_DX;
        Matrix  ConstitutiveMatrix;
        Matrix  DeltaPosition;

        //large displacement
        Vector  CurrentCurvatureVector;
        Vector  PreviousCurvatureVector;

        Vector  CurrentStepRotationVector;

        Vector  CurrentStrainResultantsVector;
        Vector  PreviousStrainResultantsVector;

        Vector  InitialAxisPositionDerivatives;
        Vector  CurrentAxisPositionDerivatives;
        Vector  PreviousAxisPositionDerivatives;


        Matrix  AlphaRotationMatrix;
        Matrix  AlphaRotationMatrixAsterisk;

        Matrix  CurrentRotationMatrix;
        Matrix  PreviousRotationMatrix;


        //variables including all integration points
        GeometryType::JacobiansType J;
        GeometryType::JacobiansType j;

        //section properties
        SectionProperties Section;

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


        void SetDirectors(DirectorsVariables& rDirectors)
        {
            Directors=&rDirectors;
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

        DirectorsVariables& GetDirectors()
        {
            return *Directors;
        };

        void Initialize( const unsigned int& voigt_size,
			 const unsigned int& dimension,
			 const unsigned int& number_of_nodes )
        {
	  //scalars
	  PointNumber = 0;
	  Length = 0;
	  detJ  = 1;
	  Alpha = 1;
          DeltaTime = 0;

	  //vectors
	  StrainVector.resize(voigt_size,false);
          StressVector.resize(voigt_size,false);
	  N.resize(number_of_nodes,false);
	  noalias(StrainVector) = ZeroVector(voigt_size);
	  noalias(StressVector) = ZeroVector(voigt_size);
	  noalias(N) = ZeroVector(number_of_nodes);
	  //matrices
	  B.resize(voigt_size, (dimension-1) * 3 * number_of_nodes,false);
	  DN_DX.resize(number_of_nodes, 1,false); //Local 1D
	  ConstitutiveMatrix.resize(voigt_size, voigt_size,false);
	  DeltaPosition.resize(number_of_nodes, dimension,false);

	  noalias(DN_DX) = ZeroMatrix(number_of_nodes, 1);
	  noalias(ConstitutiveMatrix) = ZeroMatrix(voigt_size, voigt_size);
	  noalias(DeltaPosition) = ZeroMatrix(number_of_nodes, dimension);

	  //large displacement
	  CurrentRotationMatrix.resize(dimension, dimension,false);
	  InitialAxisPositionDerivatives.resize(dimension,false);
	  CurrentAxisPositionDerivatives.resize(dimension,false);
	  PreviousAxisPositionDerivatives.resize(dimension,false);

	  noalias(CurrentRotationMatrix) = ZeroMatrix(dimension, dimension);
	  noalias(InitialAxisPositionDerivatives) = ZeroVector(dimension);
	  noalias(CurrentAxisPositionDerivatives) = ZeroVector(dimension);
	  noalias(PreviousAxisPositionDerivatives) = ZeroVector(dimension);

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
	  Directors = NULL;
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

    ///Type for element variables
    typedef ElementData                                 ElementDataType;

    ///@name Life Cycle
    ///@{

    /// Default constructors
    BeamElement(IndexType NewId, GeometryType::Pointer pGeometry);

    BeamElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    ///Copy constructor
    BeamElement(BeamElement const& rOther);

    /// Destructor.
    ~BeamElement() override;


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


    //************* GETTING METHODS

    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */

    IntegrationMethod GetIntegrationMethod() const override;

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

    //GET:
    /**
     * Get on rVariable a double Value from the Element Constitutive Law
     */
    void GetValueOnIntegrationPoints( const Variable<double>& rVariable,
				      std::vector<double>& rValues,
				      const ProcessInfo& rCurrentProcessInfo ) override;

    /**
     * Get on rVariable a Vector Value from the Element Constitutive Law
     */
    void GetValueOnIntegrationPoints( const Variable< array_1d<double, 3 > >& rVariable,
				      std::vector< array_1d<double, 3 > >& rValues,
				      const ProcessInfo& rCurrentProcessInfo ) override;


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

    //on integration points:
    /**
     * Calculate a double Variable on the Element Constitutive Law
     */
    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
				      std::vector<double>& rOutput,
				      const ProcessInfo& rCurrentProcessInfo) override;


    /**
     * Calculate a double Variable on the Element Constitutive Law
     */
    void CalculateOnIntegrationPoints( const Variable< array_1d<double, 3 > >& rVariable,
                                       std::vector< array_1d<double, 3 > >& Output,
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
        buffer << "Large Displacement Beam Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Large Displacement Beam Element #" << Id();
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
     * Currently selected integration methods
     */
    IntegrationMethod mThisIntegrationMethod;

    /**
     * Global to Local Quaternion for Global to Local tensor transformation SPATIAL
     */
    QuaternionType  mInitialLocalQuaternion;

    ///@}
    ///@name Protected Operators
    ///@{
    BeamElement() {};

    //constexpr const std::size_t& Dimension() const {return GetGeometry().WorkingSpaceDimension();}

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Increases the integration method in the "increment" order
     */
    void IncreaseIntegrationMethod(IntegrationMethod& rThisIntegrationMethod,
				   unsigned int increment) const;


    /**
     * Calculates the elemental contributions
     * \f$ K^e = w\,B^T\,D\,B \f$ and
     * \f$ r^e \f$
     */

    virtual void CalculateElementalSystem( LocalSystemComponents& rLocalSystem,
					   ProcessInfo& rCurrentProcessInfo );


    /**
     * Calculates the elemental dynamic contributions
      */
    virtual void CalculateDynamicSystem( LocalSystemComponents& rLocalSystem,
					 ProcessInfo& rCurrentProcessInfo );

    /**
     * Get element size from the dofs
     */
    virtual unsigned int GetDofsSize();

    /**
     * Initialize System Matrices
     */
    void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
				  VectorType& rRightHandSideVector,
				  Flags& rCalculationFlags);


    /**
     * Transform Vector Variable from Global Frame to the Spatial Local Frame
     */
    Vector& MapToInitialLocalFrame(Vector& rVariable, unsigned int PointNumber = 0);

    /**
     * Transform Vector Variable form Spatial Frame to Global Frame
     */
    virtual void MapLocalToGlobal(ElementDataType& rVariables, Matrix& rVariable);

    /**
     * Transform Vector Variable form Spatial Frame to Global Frame
     */
    virtual void MapLocalToGlobal(ElementDataType& rVariables, VectorType& rVector);


    /**
     * Get Current Value in the integration point in the Local Reference configuration, buffer 0 with FastGetSolutionStepValue
     */
    Vector& GetLocalCurrentValue(const Variable<array_1d<double,3> >&rVariable,
				 Vector& rValue,
				 const Vector& rN);


    /**
     * Get Reference Value in the integration point in the Local Reference configuration, buffer 1 with FastGetSolutionStepValue
     */
    Vector& GetLocalPreviousValue(const Variable<array_1d<double,3> >&rVariable,
				  Vector& rValue,
				  const Vector& rN);


    /**
     * Get Current Value, buffer 0 with FastGetSolutionStepValue
     */
    Vector& GetNodalCurrentValue(const Variable<array_1d<double,3> >&rVariable,
				 Vector& rValue,
				 const unsigned int& rNode);

    /**
     * Get Reference Value, buffer 1 with FastGetSolutionStepValue
     */
    Vector& GetNodalPreviousValue(const Variable<array_1d<double,3> >&rVariable,
				  Vector& rValue,
				  const unsigned int& rNode);

    /**
     * Initialize Element General Variables
     */
    virtual void InitializeElementData(ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo);


    /**
     * Calculation of the Section Properties
     */
    void CalculateSectionProperties(SectionProperties & rSection);


    /**
     * Calculate Element Kinematics
     */
    virtual void CalculateKinematics(ElementDataType& rVariables,
                                     const unsigned int& rPointNumber);

    /**
     * Calculate Element Constitutive Matrix
     */
    virtual void CalculateConstitutiveMatrix(ElementDataType& rVariables);


    /**
     * Calculate Element Material Constitutive Matrix
     */
    void CalculateMaterialConstitutiveMatrix(Matrix& rConstitutiveMatrix, ElementDataType& rVariables);


    /**
     * Calculate Element Stress Resultants and Couples
     */
    virtual void CalculateStressResultants(ElementDataType& rVariables, const unsigned int& rPointNumber);

    /**
     * Calculation and addition of the matrices of the LHS
     */
    virtual void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, double& rIntegrationWeight);

    /**
     * Calculation and addition of the vectors of the RHS
     */
    virtual void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, Vector& VolumeForce, double& rIntegrationWeight);

   /**
     * Calculation of the Integration Weight
     */
    virtual double& CalculateIntegrationWeight(double& rIntegrationWeight);

    /**
     * Calculation of the Material Stiffness Matrix. Kuum = BT * C * B
     */

    virtual void CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
                                     ElementDataType& rVariables,
                                     double& rIntegrationWeight);



    /**
     * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
     */
    virtual void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
                                     ElementDataType& rVariables,
                                     double& rIntegrationWeight);

    /**
     * Calculation of the External Forces Vector. Fe = N * t + N * b
     */
    virtual void CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
					       ElementDataType& rVariables,
					       Vector& rVolumeForce,
					       double& rIntegrationWeight);


    /**
      * Calculation of the Tangent Intertia Matrix
      */
    virtual void CalculateAndAddInertiaLHS(MatrixType& rLeftHandSideMatrix,
					   ElementDataType& rVariables,
					   ProcessInfo& rCurrentProcessInfo,
					   double& rIntegrationWeight);

    /**
      * Calculation of the Inertial Forces Vector
      */
    virtual void CalculateAndAddInertiaRHS(VectorType& rRightHandSideVector,
					   ElementDataType& rVariables,
					   ProcessInfo& rCurrentProcessInfo,
					   double& rIntegrationWeight);

    /**
      * Calculation of the Internal Forces Vector. Fi = B * sigma
      */
    virtual void CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
					       ElementDataType & rVariables,
					       double& rIntegrationWeight);


    /**
     * Calculation of the Rotation tensor
     */
    void CalculateLocalAxesMatrix(Matrix& rRotationMatrix);


    /**
     * Calculation of the Volume Force of the Element
     */
    virtual Vector& CalculateVolumeForce(Vector& rVolumeForce, const Vector& rN);

    /**
     * Calculation of Element Mass
     */
    double& CalculateTotalMass( SectionProperties& Section, double& rTotalMass );

    /**
     * Calculation Complementary Method : Inertial Matrix Calculation Part 2
     */
    void CalculateInertiaDyadic(SectionProperties& rSection, Matrix& rInertiaDyadic);


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


}; // Class BeamElement

} // namespace Kratos.
#endif //  KRATOS_BEAM_ELEMENT_H_INCLUDED defined
