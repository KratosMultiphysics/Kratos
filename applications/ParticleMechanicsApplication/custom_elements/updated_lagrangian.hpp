//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//


#if !defined(KRATOS_UPDATED_LAGRANGIAN_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
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

/// Large Displacement Lagrangian Element for 3D and 2D geometries. (base class)

/**
 * Implements a Large Displacement Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D (base class)
 */

class UpdatedLagrangian
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
    KRATOS_CLASS_POINTER_DEFINITION( UpdatedLagrangian );
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
        //const GeometryType::ShapeFunctionsGradientsType* pDN_De;
        const Matrix* pDN_De;
        const Vector* pNcontainer;

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
        Matrix  F;
        Matrix  FT;
        Matrix  F0;
        Matrix  DN_DX;
        Matrix  DN_De;
        Matrix  ConstitutiveMatrix;

        //variables including all integration points
        //GeometryType::JacobiansType J;
        //GeometryType::JacobiansType j;
        Matrix J;
        Matrix j;
        Matrix  DeltaPosition;
        Matrix CurrentDisp;
        Matrix PreviousDisp;

        /**
         * sets the value of a specified pointer variable
         */
        void SetShapeFunctionsGradients(const Matrix &rDN_De)
        {
            pDN_De=&rDN_De;
        };

        void SetShapeFunctions(const Vector& rNcontainer)
        {
            pNcontainer=&rNcontainer;
        };


        /**
         * returns the value of a specified pointer variable
         */
        const Matrix& GetShapeFunctionsGradients()
        {
            return *pDN_De;
        };

        const Vector& GetShapeFunctions()
        {
            return *pNcontainer;
        };


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
        void SetLeftHandSideMatrix( MatrixType& rLeftHandSideMatrix )
        {
            mpLeftHandSideMatrix = &rLeftHandSideMatrix;
        };
        void SetLeftHandSideMatrices( std::vector<MatrixType>& rLeftHandSideMatrices )
        {
            mpLeftHandSideMatrices = &rLeftHandSideMatrices;
        };
        void SetLeftHandSideVariables(const std::vector< Variable< MatrixType > >& rLeftHandSideVariables )
        {
            mpLeftHandSideVariables = &rLeftHandSideVariables;
        };

        void SetRightHandSideVector( VectorType& rRightHandSideVector )
        {
            mpRightHandSideVector = &rRightHandSideVector;
        };
        void SetRightHandSideVectors( std::vector<VectorType>& rRightHandSideVectors )
        {
            mpRightHandSideVectors = &rRightHandSideVectors;
        };
        void SetRightHandSideVariables(const std::vector< Variable< VectorType > >& rRightHandSideVariables )
        {
            mpRightHandSideVariables = &rRightHandSideVariables;
        };


        /**
         * returns the value of a specified pointer variable
         */
        MatrixType& GetLeftHandSideMatrix()
        {
            return *mpLeftHandSideMatrix;
        };
        std::vector<MatrixType>& GetLeftHandSideMatrices()
        {
            return *mpLeftHandSideMatrices;
        };
        const std::vector< Variable< MatrixType > >& GetLeftHandSideVariables()
        {
            return *mpLeftHandSideVariables;
        };

        VectorType& GetRightHandSideVector()
        {
            return *mpRightHandSideVector;
        };
        std::vector<VectorType>& GetRightHandSideVectors()
        {
            return *mpRightHandSideVectors;
        };
        const std::vector< Variable< VectorType > >& GetRightHandSideVariables()
        {
            return *mpRightHandSideVariables;
        };

    };


public:


    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    UpdatedLagrangian();


    /// Default constructors
    UpdatedLagrangian(IndexType NewId, GeometryType::Pointer pGeometry);

    UpdatedLagrangian(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    UpdatedLagrangian(UpdatedLagrangian const& rOther);

    /// Destructor.
    ~UpdatedLagrangian() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    UpdatedLagrangian& operator=(UpdatedLagrangian const& rOther);

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
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    /**
     * clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;


    //************* GETTING METHODS

    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
    //IntegrationMethod GetIntegrationMethod() const;

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
    //virtual void SetValueOnIntegrationPoints(const Variable<double>& rVariable, double& rValues, ProcessInfo& rCurrentProcessInfo);

    /**
     * Set a Vector Value on the Element Constitutive Law
     */
    //void SetValueOnIntegrationPoints(const Variable<Vector>& rVariable, Vector& rValues, ProcessInfo& rCurrentProcessInfo);

    /**
     * Set a Matrix Value on the Element Constitutive Law
     */
    //void SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, Matrix& rValues, ProcessInfo& rCurrentProcessInfo);

    /**
    * Set a Constitutive Law Value
    */
    //void SetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
    //ConstitutiveLaw::Pointer& rValues,
    //ProcessInfo& rCurrentProcessInfo );


    //GET:
    /**
     * Get on rVariable a double Value from the Element Constitutive Law
     */
    //virtual void GetValueOnIntegrationPoints(const Variable<double>& rVariable, double& rValues, ProcessInfo& rCurrentProcessInfo);

    /**
     * Get on rVariable a Vector Value from the Element Constitutive Law
     */
    //virtual void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, Vector& rValues, ProcessInfo& rCurrentProcessInfo);

    /**
     * Get on rVariable a Matrix Value from the Element Constitutive Law
     */
    //virtual void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, Matrix& rValues, ProcessInfo& rCurrentProcessInfo);

    /**
     * Get a Constitutive Law Value
     */
    //void GetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
    //ConstitutiveLaw::Pointer& rValues,
    //ProcessInfo& rCurrentProcessInfo );



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

    void Calculate(const Variable<double>& rVariable,
                   double& Output,
                   const ProcessInfo& rCurrentProcessInfo) override;

    void Calculate(const Variable<array_1d<double, 3 > >& rVariable,
                   array_1d<double, 3 > & Output,
                   const ProcessInfo& rCurrentProcessInfo) override;

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

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override;

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
                              ProcessInfo& rCurrentProcessInfo) override;

    /**
      * this is called during the assembling process in order
      * to calculate the elemental right hand side vector only
      * @param rRightHandSideVector: the elemental right hand side vector
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this function provides a more general interface to the element.
     * it is designed so that rRHSvariables are passed TO the element
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
    void CalculateLeftHandSide (MatrixType& rLeftHandSideMatrix,
                                ProcessInfo& rCurrentProcessInfo) override;

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


    /**
     * this function is designed to make the element to assemble an rRHS vector
     * identified by a variable rRHSVariable by assembling it to the nodes on the variable
     * rDestinationVariable.
     * @param rRHSVector: input variable containing the RHS vector to be assembled
     * @param rRHSVariable: variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable: variable in the database to which the rRHSvector will be assembled
      * @param rCurrentProcessInfo: the current process info instance
     */
    //virtual void AddExplicitContribution(const VectorType& rRHSVector,
    //const Variable<VectorType>& rRHSVariable,
    //Variable<array_1d<double,3> >& rDestinationVariable,
    //const ProcessInfo& rCurrentProcessInfo);

    //on integration points:
    /**
     * Calculate a double Variable on the Element Constitutive Law
     */
    //void CalculateOnIntegrationPoints(const Variable<double>& rVariable, double& rOutput, ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculate a Vector Variable on the Element Constitutive Law
     */
    //void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, Vector& rOutput, ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculate a Matrix Variable on the Element Constitutive Law
     */
    //void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, Matrix& rOutput, ProcessInfo& rCurrentProcessInfo);




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
        buffer << "MPM Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MPM Element #" << Id();
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
     * Container for historical total elastic deformation measure F0 = dx/dX
     */
    Matrix mDeformationGradientF0;
    /**
     * Container for the total deformation gradient determinants
     */
    double mDeterminantF0;
    /**
     * Container for historical inverse of Jacobian at reference configuration invJ0
     */
    Matrix mInverseJ0;
    Matrix mInverseJ;
    /**
     * Container for the total Jacobian determinants
     */
    double mDeterminantJ0;

    /**
     * Currently selected integration methods
     */
    //IntegrationMethod mThisIntegrationMethod;

    /**
     * Container for constitutive law instances on each integration point
     */
    ConstitutiveLaw::Pointer mConstitutiveLawVector;


    /**
     * Finalize and Initialize label
     */
    bool mFinalizedStep;


    ///@}
    ///@name Protected Operators
    ///@{

    /**
     * Calculates the elemental contributions
     * \f$ K^e = w\,B^T\,D\,B \f$ and
     * \f$ r^e \f$
     */
    virtual void CalculateElementalSystem(LocalSystemComponents& rLocalSystem,
                                          ProcessInfo& rCurrentProcessInfo);
    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Calculation and addition of the matrices of the LHS
     */

    virtual void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                                    GeneralVariables& rVariables,
                                    const double& rIntegrationWeight);

    /**
     * Calculation and addition of the vectors of the RHS
     */

    virtual void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
                                    GeneralVariables& rVariables,
                                    Vector& rVolumeForce,
                                    const double& rIntegrationWeight);


    /**
     * Calculation of the Material Stiffness Matrix. Kuum = BT * C * B
     */

    virtual void CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
                                     GeneralVariables& rVariables,
                                     const double& rIntegrationWeight);

    /**
     * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
     */
    virtual void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
                                     GeneralVariables& rVariables,
                                     const double& rIntegrationWeight);


    /**
     * Calculation of the External Forces Vector. Fe = N * t + N * b
     */
    virtual void CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
            GeneralVariables& rVariables,
            Vector& rVolumeForce,
            const double& rIntegrationWeight);


    /**
      * Calculation of the Internal Forces Vector. Fi = B * sigma
      */
    virtual void CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
            GeneralVariables & rVariables,
            const double& rIntegrationWeight);


    /**
     * Set Variables of the Element to the Parameters of the Constitutive Law
     */
    virtual void SetGeneralVariables(GeneralVariables& rVariables,
                                     ConstitutiveLaw::Parameters& rValues);


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
    void ResetConstitutiveLaw() override;


    /**
     * Clear Nodal Forces
     */
    void ClearNodalForces ();
    /**
     * Clear Nodal Displacement Velocity and Acceleration
     */

    /**
     * Calculate Element Kinematics
     */
    virtual void CalculateKinematics(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo);



    /**
     * Calculation of the Current Displacement
     */
    Matrix& CalculateCurrentDisp(Matrix & rCurrentDisp, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Correct Precision Errors (for rigid free movements)
     */
    void DecimalCorrection(Vector& rVector);


    /**
     * Initialize Element General Variables
     */
    virtual void InitializeGeneralVariables(GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo);


    /**
      * Finalize Element Internal Variables
      */
    virtual void FinalizeStepVariables(GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo);

    /**
      * Update the position of the MP or Gauss point when Finalize Element Internal Variables is called
      */

    virtual void UpdateGaussPoint(GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Get the Historical Deformation Gradient to calculate after finalize the step
     */
    virtual void GetHistoricalVariables( GeneralVariables& rVariables);


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



    // /**
    //* Calculation of the Velocity Gradient
    //*/
    //void CalculateVelocityGradient(const Matrix& rDN_DX,
    //Matrix& rDF );

    /**
     * Calculation of the Deformation Matrix  BL
     */
    virtual void CalculateDeformationMatrix(Matrix& rB,
                                            Matrix& rF,
                                            Matrix& rDN_DX);

    /**
     * Calculation of the Integration Weight
     */
    virtual double& CalculateIntegrationWeight(double& rIntegrationWeight);

    /**
     * Calculate Jacobian in a given point
     */
    virtual Matrix& MPMJacobian(Matrix& rResult, const array_1d<double,3>& rPoint);

    /**
     * Calculate Jacobian in a given point and given a delta position
     */
    virtual Matrix& MPMJacobianDelta(Matrix& rResult, array_1d<double,3>& rPoint, Matrix& rDeltaPosition);

    /**
     * Calculate Shape Function Values in a given point
     */

    virtual Vector& MPMShapeFunctionPointValues(Vector& rResult, array_1d<double,3>& rPoint);

    /**
     * Calculate Shape Function grandient local Values in a given point in 3 dimension
     */
    virtual Matrix& MPMShapeFunctionsLocalGradients(Matrix& rResult);
    /**
     * Calculate local coordinated of a given point in 3 dimension
     */
    virtual Vector& MPMLocalCoordinates(Vector& rResult, array_1d<double,3>& rPoint);

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

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;


    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class UpdatedLagrangian

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_UPDATED_LAGRANGIAN_H_INCLUDED  defined
