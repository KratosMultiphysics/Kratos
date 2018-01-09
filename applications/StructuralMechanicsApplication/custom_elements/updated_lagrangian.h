// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//


#if !defined(KRATOS_UPDATED_LAGRANGIAN_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_H_INCLUDED


// System includes


// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/base_solid_element.h"
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

/// Updated Lagrangian element for 2D and 3D geometries.

/**
 * Implements a total Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 2D and 3D
 */

class UpdatedLagrangian
    : public BaseSolidElement
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

    /// Counted pointer of UpdatedLagrangian
    KRATOS_CLASS_POINTER_DEFINITION(UpdatedLagrangian);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    UpdatedLagrangian(IndexType NewId, GeometryType::Pointer pGeometry);
    UpdatedLagrangian(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    ~UpdatedLagrangian() override;

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
     * Called at the beginning of each solution step
     * @param rCurrentProcessInfo: the current process info instance
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;
    
    /**
     * Called at the end of eahc solution step
     * @param rCurrentProcessInfo: the current process info instance
     */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;
    
    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
    //TODO: ADD THE OTHER CREATE FUNCTION
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

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

    /* Historical total elastic deformation measure */
    bool mF0Computed;           // To avoid computing more than once the historical total elastic deformation measure 
    std::vector<double> mDetF0; // The historical total elastic deformation measure determinant
    std::vector<Matrix> mF0;    // The historical total elastic deformation measure
    
    ///@}
    ///@name Protected Operators
    ///@{
    
    UpdatedLagrangian() : BaseSolidElement()
    {
    }

    /**
     * Gives the StressMeasure used
     */
    ConstitutiveLaw::StressMeasure GetStressMeasure() const override;
    
    /**
     * It updates the historical database
     * @param rThisKinematicVariables: The kinematic variables to be calculated 
     * @param PointNumber: The integration point considered
     */ 
    void UpdateHistoricalDatabase(
        KinematicVariables& rThisKinematicVariables,
        const unsigned int PointNumber
        );
        
    /**
     * This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix: The LHS
     * @param rRightHandSideVector: The RHS
     * @param rCurrentProcessInfo: The current process info instance
     * @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag: The flag to set if compute the RHS
     */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix, 
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
        ) override;
        
    /**
     * This functions updates the kinematics variables
     * @param rThisKinematicVariables: The kinematic variables to be calculated 
     * @param PointNumber: The integration point considered
     */ 
    void CalculateKinematicVariables(
        KinematicVariables& rThisKinematicVariables,
        const unsigned int PointNumber,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints
        ) override;
        
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
    void CalculateConstitutiveVariables(
        KinematicVariables& rThisKinematicVariables, 
        ConstitutiveVariables& rThisConstitutiveVariables, 
        ConstitutiveLaw::Parameters& rValues,
        const unsigned int PointNumber,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure,
        const Vector Displacements = ZeroVector(1)
        ) override;
    
    /**
     * This functions calculate the derivatives in the reference frame
     */ 
    double CalculateDerivativesOnReferenceConfiguration(
        Matrix& J0, 
        Matrix& InvJ0, 
        Matrix& DN_DX, 
        const unsigned int PointNumber,
        IntegrationMethod ThisIntegrationMethod
        ) override;
    
    ///@}
    ///@name Protected Operations
    ///@{
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

    void CalculateBodyForces(
        Vector& BodyForce,
        const ProcessInfo& CurrentProcessInfo
        );

    void InitializeVariables();
    
    void CalculateB(
        Matrix& B,
        const Matrix& DN_DX,
        const unsigned int StrainSize,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const unsigned int PointNumber
        );
    
    /**
     * It returns the reference configuration deformation gradient determinant
     * @param PointNumber: The integration point considered
     * @return The reference configuration deformation gradient determinant
     */
    double ReferenceConfigurationDeformationGradientDeterminant(const unsigned PointNumber) const;
    
    /**
     * It returns the reference configuration deformation gradient
     * @param PointNumber: The integration point considered
     * @return The reference configuration deformation gradient
     */
    Matrix ReferenceConfigurationDeformationGradient(const unsigned PointNumber) const;
    
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
    /// Assignment operator.
    //UpdatedLagrangian& operator=(const UpdatedLagrangian& rOther);
    /// Copy constructor.
    //UpdatedLagrangian(const UpdatedLagrangian& rOther);
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
