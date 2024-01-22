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

class MPMUpdatedLagrangian
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

    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    /// Counted pointer of LargeDisplacementElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( MPMUpdatedLagrangian );
    ///@}

protected:

    /**
     * Flags related to the element computation
     */

    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR_WITH_COMPONENTS );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX_WITH_COMPONENTS );

    struct MaterialPointVariables
    {
    public:
        // Material Point Position
        CoordinatesArrayType xg;
        // MP_MASS
        double mass;
        // MP_DENSITY
        double density;
        // MP_VOLUME
        double volume;

        // MP_DISPLACEMENT
        array_1d<double, 3> displacement;
        // MP_VELOCITY
        array_1d<double, 3> velocity;
        // MP_ACCELERATION
        array_1d<double, 3> acceleration;

        // MP_VOLUME_ACCELERATION
        array_1d<double, 3> volume_acceleration;

        // MP_CAUCHY_STRESS_VECTOR
        Vector cauchy_stress_vector;
        // MP_ALMANSI_STRAIN_VECTOR
        Vector almansi_strain_vector;

        // MP_DELTA_PLASTIC_STRAIN
        double delta_plastic_strain;
        // MP_DELTA_PLASTIC_VOLUMETRIC_STRAIN
        double delta_plastic_volumetric_strain;
        // MP_DELTA_PLASTIC_DEVIATORIC_STRAIN
        double delta_plastic_deviatoric_strain;
        // MP_EQUIVALENT_PLASTIC_STRAIN
        double equivalent_plastic_strain;
        // MP_ACCUMULATED_PLASTIC_VOLUMETRIC_STRAIN
        double accumulated_plastic_volumetric_strain;
        // MP_ACCUMULATED_PLASTIC_DEVIATORIC_STRAIN
        double accumulated_plastic_deviatoric_strain;

        explicit MaterialPointVariables()
        {
            // MP_MASS
            mass = 1.0;
            // MP_DENSITY
            density = 1.0;
            // MP_VOLUME
            volume = 1.0;

            // MP_DELTA_PLASTIC_STRAIN
            delta_plastic_strain = 1.0;
            // MP_DELTA_PLASTIC_VOLUMETRIC_STRAIN
            delta_plastic_volumetric_strain = 1.0;
            // MP_DELTA_PLASTIC_DEVIATORIC_STRAIN
            delta_plastic_deviatoric_strain = 1.0;
            // MP_EQUIVALENT_PLASTIC_STRAIN
            equivalent_plastic_strain = 1.0;
            // MP_ACCUMULATED_PLASTIC_VOLUMETRIC_STRAIN
            accumulated_plastic_volumetric_strain = 1.0;
            // MP_ACCUMULATED_PLASTIC_DEVIATORIC_STRAIN
            accumulated_plastic_deviatoric_strain = 1.0;
        }

    private:

        ///@}
        ///@name Serialization
        ///@{
        friend class Serializer;

        void save( Serializer& rSerializer ) const
        {
            rSerializer.save("xg",xg);
            rSerializer.save("mass",mass);
            rSerializer.save("density",density);
            rSerializer.save("volume",volume);
            rSerializer.save("displacement",displacement);
            rSerializer.save("velocity",velocity);
            rSerializer.save("acceleration",acceleration);
            rSerializer.save("volume_acceleration",volume_acceleration);
            rSerializer.save("cauchy_stress_vector",cauchy_stress_vector);
            rSerializer.save("almansi_strain_vector",almansi_strain_vector);
            rSerializer.save("delta_plastic_strain",delta_plastic_strain);
            rSerializer.save("delta_plastic_volumetric_strain",delta_plastic_volumetric_strain);
            rSerializer.save("delta_plastic_deviatoric_strain",delta_plastic_deviatoric_strain);
            rSerializer.save("equivalent_plastic_strain",equivalent_plastic_strain);
            rSerializer.save("accumulated_plastic_volumetric_strain",accumulated_plastic_volumetric_strain);
            rSerializer.save("accumulated_plastic_deviatoric_strain",accumulated_plastic_deviatoric_strain);
        }

        void load( Serializer& rSerializer )
        {
            rSerializer.load("xg",xg);
            rSerializer.load("mass",mass);
            rSerializer.load("density",density);
            rSerializer.load("volume",volume);
            rSerializer.load("displacement",displacement);
            rSerializer.load("velocity",velocity);
            rSerializer.load("acceleration",acceleration);
            rSerializer.load("volume_acceleration",volume_acceleration);
            rSerializer.load("cauchy_stress_vector",cauchy_stress_vector);
            rSerializer.load("almansi_strain_vector",almansi_strain_vector);
            rSerializer.load("delta_plastic_strain",delta_plastic_strain);
            rSerializer.load("delta_plastic_volumetric_strain",delta_plastic_volumetric_strain);
            rSerializer.load("delta_plastic_deviatoric_strain",delta_plastic_deviatoric_strain);
            rSerializer.load("equivalent_plastic_strain",equivalent_plastic_strain);
            rSerializer.load("accumulated_plastic_volumetric_strain",accumulated_plastic_volumetric_strain);
            rSerializer.load("accumulated_plastic_deviatoric_strain",accumulated_plastic_deviatoric_strain);
        }
        ///@}
    };

    /**
     * Parameters to be used in the Element as they are. Direct interface to Parameters Struct
     */

    struct GeneralVariables
    {
    public:

        StressMeasureType StressMeasure;

        // For axisymmetric use only
        double  CurrentRadius;
        double  ReferenceRadius;

        // General variables for large displacement use
        double  detF;
        double  detF0;
        double  detFT;
        Vector  StrainVector;
        Vector  StressVector;
        Matrix  B;
        Matrix  F;
        Matrix  FT;
        Matrix  F0;
        Matrix  DN_DX;
        Matrix  ConstitutiveMatrix;

        // Variables including all integration points
        Matrix CurrentDisp;
    };

public:


    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    MPMUpdatedLagrangian();


    /// Default constructors
    MPMUpdatedLagrangian(IndexType NewId, GeometryType::Pointer pGeometry);

    MPMUpdatedLagrangian(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    MPMUpdatedLagrangian(MPMUpdatedLagrangian const& rOther);

    /// Destructor.
    ~MPMUpdatedLagrangian() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    MPMUpdatedLagrangian& operator=(MPMUpdatedLagrangian const& rOther);

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

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

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
     * Sets on rElementalDofList the degrees of freedom of the considered element geometry
     */
    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     */
    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * Sets on rValues the nodal displacements
     */
    void GetValuesVector(Vector& rValues, int Step = 0) const override;

    /**
     * Sets on rValues the nodal velocities
     */
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override;

    /**
     * Sets on rValues the nodal accelerations
     */
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const  override;

    //************* STARTING - ENDING  METHODS

    /**
      * Called to initialize the element.
      * Must be called before any calculation is done
      */
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Called at the beginning of each solution step
     */
    virtual void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Called at the end of eahc solution step
     */
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;


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
                              const ProcessInfo& rCurrentProcessInfo) override;

    /**
      * this is called during the assembling process in order
      * to calculate the elemental right hand side vector only
      * @param rRightHandSideVector: the elemental right hand side vector
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side vector only
     * @param rLeftHandSideVector: the elemental left hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLeftHandSide (MatrixType& rLeftHandSideMatrix,
                                const ProcessInfo& rCurrentProcessInfo) override;

    /**
      * this is called during the assembling process in order
      * to calculate the elemental mass matrix
      * @param rMassMatrix: the elemental mass matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateMassMatrix(MatrixType& rMassMatrix,
                             const ProcessInfo& rCurrentProcessInfo) override;

    /**
      * this is called during the assembling process in order
      * to calculate the elemental damping matrix
      * @param rDampingMatrix: the elemental damping matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateDampingMatrix(MatrixType& rDampingMatrix,
                                const ProcessInfo& rCurrentProcessInfo) override;


    void AddExplicitContribution(const VectorType& rRHSVector,
                                 const Variable<VectorType>& rRHSVariable,
                                 const Variable<array_1d<double, 3> >& rDestinationVariable,
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
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;


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
    ///@name Access Get Values
    ///@{

    void CalculateOnIntegrationPoints(const Variable<bool>& rVariable,
        std::vector<bool>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<int>& rVariable,
        std::vector<int>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
        std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access Set Values
    ///@{

    void SetValuesOnIntegrationPoints(const Variable<int>& rVariable,
        const std::vector<int>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void SetValuesOnIntegrationPoints(const Variable<double>& rVariable,
        const std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void SetValuesOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
        const std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable,
        const std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{

    MaterialPointVariables mMP;

    /**
     * Container for historical total elastic deformation measure F0 = dx/dX
     */
    Matrix mDeformationGradientF0;
    /**
     * Container for the total deformation gradient determinants
     */
    double mDeterminantF0;

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

    virtual SizeType GetNumberOfDofs() {
        return GetGeometry().WorkingSpaceDimension();
    }

    /**
     * Calculates the elemental contributions
     * \f$ K^e = w\,B^T\,D\,B \f$ and
     * \f$ r^e \f$
     */
    virtual void CalculateElementalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag);

    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Calculation and addition of the matrices of the LHS
     */

    virtual void CalculateAndAddLHS(
        MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight,
        const ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculation and addition of the vectors of the RHS
     */

    virtual void CalculateAndAddRHS(
        VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        Vector& rVolumeForce,
        const double& rIntegrationWeight,
        const ProcessInfo& rCurrentProcessInfo);


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
                                     const double& rIntegrationWeight,
                                     const bool IsAxisymmetric = false);


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

    /// Calculation of the Explicit Stresses from velocity gradient.
    virtual void CalculateExplicitStresses(const ProcessInfo& rCurrentProcessInfo,
        GeneralVariables& rVariables);


    /**
     * Set Variables of the Element to the Parameters of the Constitutive Law
     */
    virtual void SetGeneralVariables(GeneralVariables& rVariables,
                                     ConstitutiveLaw::Parameters& rValues, const Vector& rN);


    /**
     * Initialize Material Properties on the Constitutive Law
     */
    virtual void InitializeMaterial (const ProcessInfo& rCurrentProcessInfo);


    /**
     * Reset the Constitutive Law Parameters
     */
    void ResetConstitutiveLaw() override;


    /**
     * Clear Nodal Forces
     */
    void ClearNodalForces ();

    /**
     * Calculate Element Kinematics
     */
    virtual void CalculateKinematics(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);


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


    /**
     * Calculation of the Deformation Matrix  BL
     */
    virtual void CalculateDeformationMatrix(Matrix& rB,
                                            const Matrix& rDN_DX,
                                            const Matrix& rN,
                                            const bool IsAxisymmetric = false);

    /**
     * Calculation of the Integration Weight
     */
    virtual double& CalculateIntegrationWeight(double& rIntegrationWeight);


    /**
     * Calculation of the Volume Change of the Element
     */
    virtual double& CalculateVolumeChange(double& rVolumeChange, GeneralVariables& rVariables);


    /// Calculation of the Deformation Gradient F
    void CalculateDeformationGradient(const Matrix& rDN_DX, Matrix& rF, Matrix& rDeltaPosition,
        const bool IsAxisymmetric = false);

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

}; // Class MPMUpdatedLagrangian

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_UPDATED_LAGRANGIAN_H_INCLUDED  defined
