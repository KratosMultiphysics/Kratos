// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_BASE_SOLID_ELEMENT_H_INCLUDED )
#define  KRATOS_BASE_SOLID_ELEMENT_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
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
        Vector  StrainVector;
        Vector  StressVector;
        Vector  N;
        Matrix  B;
        double  detF;
        Matrix  F;
        double  detJ;
        double  detJ0;
        Matrix  J0;
        Matrix  InvJ0;
        Matrix  DN_DX;
        Matrix  D;
        
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
            detJ0 = 1.0;
            StrainVector = ZeroVector(StrainSize);
            StressVector = ZeroVector(StrainSize);
            N = ZeroVector(NumberOfNodes);
            B = ZeroMatrix(StrainSize, Dimension * NumberOfNodes);
            F = IdentityMatrix(Dimension);
            DN_DX = ZeroMatrix(NumberOfNodes, Dimension);
            D = ZeroMatrix(StrainSize, StrainSize);
            J0 = ZeroMatrix(Dimension, Dimension);
            InvJ0 = ZeroMatrix(Dimension, Dimension);
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
    void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo) override;

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
    void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo) override;
    
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
     * @param rValues: The CL parameters
     * @param PointNumber: The integration point considered
     * @param Displacements: The displacements vector
     */ 
    virtual void UpdateKinematics(
        KinematicVariables& rThisKinematicVariables, 
        ConstitutiveLaw::Parameters& rValues,
        const unsigned int PointNumber,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const Vector Displacements
        );
    
    /**
     * This functions calculate the derivatives in the reference frame
     */ 
    double CalculateDerivativesOnReference(
        Matrix& J0, 
        Matrix& InvJ0, 
        Matrix& DN_DX, 
        const unsigned int PointNumber,
        IntegrationMethod ThisIntegrationMethod
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
