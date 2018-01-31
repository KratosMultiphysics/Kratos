// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_BASE_SOLID_ELEMENT_H_INCLUDED )
#define  KRATOS_BASE_SOLID_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/integration_utilities.h"
#include "structural_mechanics_application_variables.h"

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
    
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION)  BaseSolidElement
    : public Element
{
protected:
    /**
     * Internal variables used in the kinematic calculations
     */
    struct KinematicVariables
    {
        Vector  N;
        Matrix  B;
        double  detF;
        Matrix  F;
        double  detJ0;
        Matrix  J0;
        Matrix  InvJ0;
        Matrix  DN_DX;
        
        /**
         * The default constructor
         * @param StrainSize: The size of the strain vector in Voigt notation
         * @param Dimension: The size of the strain vector in Voigt notation
         * @param NumberOfNodes: The size of the strain vector in Voigt notation
         */
        KinematicVariables( 
            const unsigned int& StrainSize, 
            const unsigned int& Dimension, 
            const unsigned int& NumberOfNodes 
            )
        {
            detF = 1.0;
            detJ0 = 1.0;
            N = ZeroVector(NumberOfNodes);
            B = ZeroMatrix(StrainSize, Dimension * NumberOfNodes);
            F = IdentityMatrix(Dimension);
            DN_DX = ZeroMatrix(NumberOfNodes, Dimension);
            J0 = ZeroMatrix(Dimension, Dimension);
            InvJ0 = ZeroMatrix(Dimension, Dimension);
        }
    };
    
    /**
     * Internal variables used in the kinematic calculations
     */
    struct ConstitutiveVariables
    {
        Vector StrainVector;
        Vector StressVector;
        Matrix D;
        
        /**
         * The default constructor
         * @param StrainSize: The size of the strain vector in Voigt notation
         */
        ConstitutiveVariables(const unsigned int& StrainSize)
        {
            StrainVector = ZeroVector(StrainSize);
            StressVector = ZeroVector(StrainSize);
            D = ZeroMatrix(StrainSize, StrainSize);
        }
    };
public:

    ///@name Type Definitions
    ///@{
    
    // Counted pointer of BaseSolidElement
    KRATOS_CLASS_POINTER_DEFINITION( BaseSolidElement );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    BaseSolidElement()
    {};

    // Constructor using an array of nodes
    BaseSolidElement( IndexType NewId, GeometryType::Pointer pGeometry ):Element(NewId,pGeometry)
    {};

    // Constructor using an array of nodes with properties
    BaseSolidElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ):Element(NewId,pGeometry,pProperties)
    {};

    // Destructor
    ~BaseSolidElement() override
    {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Called to initialize the element.
     * Must be called before any calculation is done
     */
    void Initialize() override;

    /**
      * This resets the constitutive law
      */
    void ResetConstitutiveLaw() override;

    /**
     * Called at the beginning of each solution step
     * @param rCurrentProcessInfo: the current process info instance
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * This is called for non-linear analysis at the beginning of the iteration process
     * @param rCurrentProcessInfo: the current process info instance
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * This is called for non-linear analysis at the beginning of the iteration process
     * @param rCurrentProcessInfo: the current process info instance
     */
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;
    
    /**
     * Called at the end of eahc solution step
     * @param rCurrentProcessInfo: the current process info instance
     */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;
    
    /**
     * Sets on rResult the ID's of the element degrees of freedom
     * @param rResult: The vector containing the equation id
     * @param rCurrentProcessInfo: The current process info instance
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo 
        ) override;

    /**
     * Sets on rElementalDofList the degrees of freedom of the considered element geometry
     * @param rElementalDofList: The vector containing the dof of the element
     * @param rCurrentProcessInfo: The current process info instance
     */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo 
        ) override;

    /**
     * Sets on rValues the nodal displacements
     * @param rValues: The values of displacements
     * @param Step: The step to be computed
     */
    void GetValuesVector(
        Vector& rValues,
        int Step = 0 
        ) override;
    
    /**
     * Sets on rValues the nodal velocities
     * @param rValues: The values of velocities
     * @param Step: The step to be computed
     */
    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step = 0 
        ) override;

    /**
     * Sets on rValues the nodal accelerations
     * @param rValues: The values of accelerations
     * @param Step: The step to be computed
     */
    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step = 0 
        ) override;

    /**
     * This function provides a more general interface to the element. 
     * It is designed so that rLHSvariables and rRHSvariables are passed to the element thus telling what is the desired output
     * @param rLeftHandSideMatrices: container with the output left hand side matrices
     * @param rLHSVariables: paramter describing the expected LHSs
     * @param rRightHandSideVectors: container for the desired RHS output
     * @param rRHSVariables: parameter describing the expected RHSs
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix, 
        VectorType& rRightHandSideVector, 
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * This is called during the assembling process in order to calculate the elemental right hand side vector only
      * @param rRightHandSideVector: the elemental right hand side vector
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector, 
        ProcessInfo& rCurrentProcessInfo
        ) override;
        
    /**
      * This is called during the assembling process in order to calculate the elemental mass matrix
      * @param rMassMatrix: the elemental mass matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo 
        ) override;

    /**
      * This is called during the assembling process in order
      * to calculate the elemental damping matrix
      * @param rDampingMatrix: the elemental damping matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfo 
        ) override;

    /**
     * Calculate a double Variable on the Element Constitutive Law
     * @param rVariable: The variable we want to get
     * @param rOutput: The values obtained int the integration points
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable, 
        std::vector<double>& rOutput, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;
        
    /**
     * Calculate a double array_1d on the Element Constitutive Law
     * @param rVariable: The variable we want to get
     * @param rOutput: The values obtained int the integration points
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3>>& rVariable, 
        std::vector<array_1d<double, 3>>& rOutput, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * Calculate a Vector Variable on the Element Constitutive Law
     * @param rVariable: The variable we want to get
     * @param rOutput: The values obtained int the integration points
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable, 
        std::vector<Vector>& rOutput, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * Calculate a Matrix Variable on the Element Constitutive Law
     * @param rVariable: The variable we want to get
     * @param rOutput: The values obtained int the integration points
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<Matrix >& rVariable, 
        std::vector< Matrix >& rOutput, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;

     /**
      * Set a double Value on the Element Constitutive Law
      * @param rVariable: The variable we want to set
      * @param rValues: The values to set in the integration points
      * @param rCurrentProcessInfo: the current process info instance
      */
    void SetValueOnIntegrationPoints(
        const Variable<double>& rVariable, 
        std::vector<double>& rValues, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;
    
     /**
      * Set a Vector Value on the Element Constitutive Law
      * @param rVariable: The variable we want to set
      * @param rValues: The values to set in the integration points
      * @param rCurrentProcessInfo: the current process info instance
      */
    void SetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable, 
        std::vector<Vector>& rValues, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;

     /**
      * Set a Matrix Value on the Element Constitutive Law
      * @param rVariable: The variable we want to set
      * @param rValues: The values to set in the integration points
      * @param rCurrentProcessInfo: the current process info instance
      */
    void SetValueOnIntegrationPoints(
        const Variable<Matrix>& rVariable, 
        std::vector<Matrix>& rValues, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * Get on rVariable a double Value from the Element Constitutive Law
     * @param rVariable: The variable we want to get
     * @param rValues: The results in the integration points
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<double>& rVariable, 
        std::vector<double>& rValues, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * Get on rVariable a array_1d Value from the Element Constitutive Law
     * @param rVariable: The variable we want to get
     * @param rValues: The results in the integration points
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<array_1d<double, 3>>& rVariable, 
        std::vector<array_1d<double, 3>>& rValues, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;
        
    /**
     * Get on rVariable a Vector Value from the Element Constitutive Law
     * @param rVariable: The variable we want to get
     * @param rValues: The results in the integration points
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable, 
        std::vector<Vector>& rValues, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * Get on rVariable a Matrix Value from the Element Constitutive Law
     * @param rVariable: The variable we want to get
     * @param rValues: The results in the integration points
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<Matrix>& rVariable, 
        std::vector<Matrix>& rValues, 
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * This method is not defined yet...
     * @param rCurrentProcessInfo: the current process info instance
     */
    void Calculate(
        const Variable<double>& rVariable, 
        double& Output,
        const ProcessInfo& rCurrentProcessInfo
        ) override;
    
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
    
protected:
    
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{
    
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector; // The vector containing the constitutive laws

    ///@}
    ///@name Protected Operators
    ///@{
    
    ///@}
    ///@name Protected Operations
    ///@{
    
    /**
     * It initializes the material
     */
    virtual void InitializeMaterial();
    
    /**
     * Gives the StressMeasure used
     */
    virtual ConstitutiveLaw::StressMeasure GetStressMeasure() const;
    
    /**
     * This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix: The LHS
     * @param rRightHandSideVector: The RHS
     * @param rCurrentProcessInfo: The current process info instance
     * @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag: The flag to set if compute the RHS
     */
    virtual void CalculateAll(
        MatrixType& rLeftHandSideMatrix, 
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
        );
    
    /**
     * This functions updates the kinematics variables
     * @param rThisKinematicVariables: The kinematic variables to be calculated 
     * @param PointNumber: The integration point considered
     * @param IntegrationPoints: The list of integration points
     */ 
    virtual void CalculateKinematicVariables(
        KinematicVariables& rThisKinematicVariables, 
        const unsigned int PointNumber,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints
        );
        
    /**
     * This functions updates the constitutive variables
     * @param rThisKinematicVariables: The kinematic variables to be calculated 
     * @param rThisConstitutiveVariables: The constitutive variables
     * @param rValues: The CL parameters
     * @param PointNumber: The integration point considered
     * @param IntegrationPoints: The list of integration points
     * @param ThisStressMeasure: The stress measure considered
     * @param Displacements: The displacements vector
     */ 
    virtual void CalculateConstitutiveVariables(
        KinematicVariables& rThisKinematicVariables, 
        ConstitutiveVariables& rThisConstitutiveVariables, 
        ConstitutiveLaw::Parameters& rValues,
        const unsigned int PointNumber,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure,
        const Vector Displacements
        );
    
    /**
     * This methods gives us a matrix with the increment of displacement
     * @return DeltaDisplacement: The matrix containing the increment of displacements
     */
    Matrix CalculateDeltaDisplacement(Matrix& DeltaDisplacement);
    
    /**
     * This functions calculate the derivatives in the reference frame
     */ 
    virtual double CalculateDerivativesOnReferenceConfiguration(
        Matrix& J0, 
        Matrix& InvJ0, 
        Matrix& DN_DX, 
        const unsigned int PointNumber,
        IntegrationMethod ThisIntegrationMethod
        );
    
    /**
     * This functions calculate the derivatives in the current frame
     */ 
    double CalculateDerivativesOnCurrentConfiguration(
        Matrix& J, 
        Matrix& InvJ, 
        Matrix& DN_DX, 
        const unsigned int PointNumber,
        IntegrationMethod ThisIntegrationMethod
        );
    
    /**
     * This function computes the body force
     * @param rGaussPoint: The gauss point considered
     */ 
    Vector GetBodyForce(
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const unsigned int PointNumber
        ) const;
    
    /**
     * Calculation of the Material Stiffness Matrix. Km = B^T * D *B
     */
    void CalculateAndAddKm(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& B,
        const Matrix& D,
        const double IntegrationWeight
        );

    /**
     * Calculation of the Geometric Stiffness Matrix. Kg = dB * S
     */
    void CalculateAndAddKg(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& DN_DX,
        const Vector& StressVector,
        const double IntegrationWeight
        );
    
    /**
     * Calculation of the RHS
     */
    void CalculateAndAddResidualVector(
        VectorType& rRightHandSideVector,
        const KinematicVariables& rThisKinematicVariables,
        const ProcessInfo& CurrentProcessInfo,
        const Vector& BodyForce,
        const Vector& StressVector,
        const double IntegrationWeight
        );
    
    /**
     * This function add the external force contribution
     */ 
    void CalculateAndAddExtForceContribution(
        const Vector& N,
        const ProcessInfo& CurrentProcessInfo,
        const Vector& BodyForce,
        VectorType& mResidualVector,
        const double IntegrationWeight
        );

    /**
     * This functions computes the integration weight to consider
     * @param ThisIntegrationMethod: The array containing the integration points
     * @param PointNumber: The id of the integration point considered
     * @param detJ: The determinant of the jacobian of the element
     */
    virtual double GetIntegrationWeight(
        const GeometryType::IntegrationPointsArrayType& ThisIntegrationMethod,
        const unsigned int PointNumber,
        const double detJ
        );
    
    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

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

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    }

}; // class BaseSolidElement.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_BASE_SOLID_ELEMENT_H_INCLUDED  defined 
