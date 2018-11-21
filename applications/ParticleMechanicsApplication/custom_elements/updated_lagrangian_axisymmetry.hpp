//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					    Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//


#if !defined(KRATOS_UPDATED_LAGRANGIAN_AXISYMMETRY_H_INCLUDED )
#define      KRATOS_UPDATED_LAGRANGIAN_AXISYMMETRY_H_INCLUDED

// System includes

// External includes
#include "custom_elements/updated_lagrangian.hpp"


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

class UpdatedLagrangianAxisymmetry
    : public UpdatedLagrangian
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
    KRATOS_CLASS_POINTER_DEFINITION( UpdatedLagrangianAxisymmetry );
    ///@}

    /**
     * Parameters to be used in the Element as they are. Direct interface to Parameters Struct
     */



public:


    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    UpdatedLagrangianAxisymmetry();


    /// Default constructors
    UpdatedLagrangianAxisymmetry(IndexType NewId, GeometryType::Pointer pGeometry);

    UpdatedLagrangianAxisymmetry(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    UpdatedLagrangianAxisymmetry(UpdatedLagrangianAxisymmetry const& rOther);

    /// Destructor.
    virtual ~UpdatedLagrangianAxisymmetry();

    ///@}
    ///@name Operators
    ///@{

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
    //virtual void Initialize();

    /**
     * Called at the beginning of each solution step
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    //************* COMPUTING  METHODS

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
     * Container for historical total elastic deformation measure F0 = dx/dX
     */
    //Matrix mDeformationGradientF0;
    /**
     * Container for the total deformation gradient determinants
     */
    //double mDeterminantF0;
    /**
     * Container for historical inverse of Jacobian at reference configuration invJ0
     */
    //Matrix mInverseJ0;
    //Matrix mInverseJ;
    /**
     * Container for the total Jacobian determinants
     */
    //double mDeterminantJ0;

    /**
     * Currently selected integration methods
     */
    //IntegrationMethod mThisIntegrationMethod;

    /**
     * Container for constitutive law instances on each integration point
     */
    //ConstitutiveLaw::Pointer mConstitutiveLawVector;


    /**
     * Finalize and Initialize label
     */
    //bool mFinalizedStep;


    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Calculation and addition of the matrices of the LHS
     */

    //virtual void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                                    //GeneralVariables& rVariables,
                                    //double& rIntegrationWeight);

    /**
     * Calculation and addition of the vectors of the RHS
     */

    //virtual void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
                                    //GeneralVariables& rVariables,
                                    //Vector& rVolumeForce,
                                    //double& rIntegrationWeight);


    /**
     * Calculation of the Material Stiffness Matrix. Kuum = BT * C * B
     */

    //virtual void CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
                                     //GeneralVariables& rVariables,
                                     //double& rIntegrationWeight);

    /**
     * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
     */
    void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
                                     GeneralVariables& rVariables,
                                     double& rIntegrationWeight) override;

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
    void CalculateKinematics(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo) override;


    virtual double CalculateRadius(Vector& N, GeometryType& Geom, std::string ThisConfiguration);

     /**
     * Initialize Element General Variables
     */
    void Initialize() override;

    /**
     * Initialize Element General Variables
     */
    void InitializeGeneralVariables(GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculation of the Green Lagrange Strain Vector
     */
    void CalculateGreenLagrangeStrain(const Matrix& rF,
            Vector& rStrainVector) override;

    /**
     * Calculation of the Almansi Strain Vector
     */
    void CalculateAlmansiStrain(const Matrix& rF, Vector& rStrainVector) override;



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
                                            Matrix& rDN_DX,
                                            Vector& rN);

    /**
     * Calculate Jacobian in a given point
     */
    Matrix& MPMJacobian(Matrix& rResult, array_1d<double,3>& rPoint) override;

    /**
     * Calculate Jacobian in a given point and given a delta position
     */
    Matrix& MPMJacobianDelta(Matrix& rResult, array_1d<double,3>& rPoint, Matrix& rDeltaPosition) override;

    /**
     * Calculate Shape Function Values in a given point
     */

    Vector& MPMShapeFunctionPointValues(Vector& rResult, array_1d<double,3>& rPoint) override;

    /**
     * Calculate Shape Function grandient local Values in a given point in 3 dimension
     */
    Matrix& MPMShapeFunctionsLocalGradients(Matrix& rResult) override;

    /**
     * Calculation of the Deformation Gradient F
     */
    virtual void CalculateDeformationGradient(const Matrix& rDN_DX,
                                      Matrix& rF,
                                      Matrix& rDeltaPosition,
                                      const double & rCurrentRadius,
                                      const double & rReferenceRadius);

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
