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

    /**
     * Called at the beginning of each solution step
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

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

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
     */
    void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
                                     GeneralVariables& rVariables,
                                     const double& rIntegrationWeight) override;

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
    Matrix& MPMJacobian(Matrix& rResult, const array_1d<double,3>& rPoint) override;

    /**
     * Calculate Jacobian in a given point and given a delta position
     */
    Matrix& MPMJacobianDelta(Matrix& rResult, const array_1d<double,3>& rPoint, const Matrix& rDeltaPosition) override;

    /**
     * Calculate Shape Function Values in a given point
     */

    Vector& MPMShapeFunctionPointValues(Vector& rResult, const array_1d<double,3>& rPoint) override;

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
