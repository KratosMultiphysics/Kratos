// KRATOS
//
//  License:		 BSD License
//					 license: ../../license.txt
//
//  Main authors:    Ramsubramanian Pazhanisamy
//                   Ricky Aristio
//


#if !defined(KRATOS_ACOUSTIC_ELEMENT_H_INCLUDED )
#define  KRATOS_ACOUSTIC_ELEMENT_H_INCLUDED

// System includes


// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "mor_application_variables.h"
#include "utilities/integration_utilities.h"
#include "utilities/geometrical_sensitivity_utility.h"
#include "includes/variables.h"


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

/**
 * @class AcousticElement
 * @ingroup MORApplication
 * @brief Acoustic element
 * @details
 * @author
 * @author
 */

class AcousticElement
    : public Element
{
public:
    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    //typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    //typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// The base element type
    typedef Element BaseType;

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    /// Counted pointer of AcousticElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AcousticElement);

    ///@}
    ///@name Life Cycle
    ///@{
    
    /// Default constructor.
    AcousticElement(IndexType NewId, GeometryType::Pointer pGeometry);
    AcousticElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    // Copy constructor
    AcousticElement(AcousticElement const& rOther)
        :BaseType(rOther)
    {};

    /// Destructor.
    ~AcousticElement() override;

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param pGeom The pointer to the geometry of the element
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param ThisNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief It creates a new element pointer and clones the previous element data
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone (
        IndexType NewId,
        NodesArrayType const& rThisNodes
        ) const override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The current process info instance
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
        buffer << "Acoustic Element #" << Id() << "\n";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Acoustic Element #" << Id() << "\n";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }


    /**
     * @brief Sets on rResult the ID's of the element degrees of freedom
     * @param rResult The vector containing the equation id
     * @param rCurrentProcessInfo The current process info instance
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
     * @param rElementalDofList The vector containing the dof of the element
     * @param rCurrentProcessInfo The current process info instance
     */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
        ) override;


    ///@}
    ///@name Friends
    ///@{
    ///@}

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
        Vector  Acoustic_pressure;
        Vector  PressureGradient;
    

        /**
         * The default constructor
         * @param StrainSize The size of the strain vector in Voigt notation
         * @param Dimension The problem dimension: 2D or 3D
         * @param NumberOfNodes The number of nodes in the element
         */
        KinematicVariables(
            const SizeType PressureGradientSize,
            const SizeType Dimension,
            const SizeType NumberOfNodes
            )
        {
            detF = 1.0;
            detJ0 = 1.0;
            N = ZeroVector(NumberOfNodes);
            B = ZeroMatrix(Dimension, NumberOfNodes);
            F = IdentityMatrix(Dimension);
            DN_DX = ZeroMatrix(NumberOfNodes, Dimension); 
            J0 = ZeroMatrix(Dimension, Dimension);
            InvJ0 = ZeroMatrix(Dimension, Dimension);
            Acoustic_pressure = ZeroVector(NumberOfNodes);
            PressureGradient = ZeroVector(Dimension);
        }
    };

    // /**
    //  * Internal variables used in the kinematic calculations
    //  */
      struct ConstitutiveVariables
    {

        Matrix D;

        ConstitutiveVariables()
        {
              D = ZeroMatrix(1, 1);
        }
    };


    /**
     * @brief This functions calculate the derivatives in the reference frame
     * @param rJ0 The jacobian in the reference configuration
     * @param rInvJ0 The inverse of the jacobian in the reference configuration
     * @param rDN_DX The gradient derivative of the shape function
     * @param PointNumber The id of the integration point considered
     * @param ThisIntegrationMethod The integration method considered
     * @return The determinant of the jacobian in the reference configuration
     */
    virtual double CalculateDerivativesOnReferenceConfiguration(
        Matrix& rJ0,
        Matrix& rInvJ0,
        Matrix& rDN_DX,
        const IndexType PointNumber,
        IntegrationMethod ThisIntegrationMethod
        ) const;
        
    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    AcousticElement() : Element()
    {
    }


    /**
     * @brief This method returns if the element provides the strain
     */
    bool UseElementProvidedStrain() const;

    /**
     * @brief This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix The LHS
     * @param rRightHandSideVector The RHS
     * @param rCurrentProcessInfo The current process info instance
     * @param CalculateStiffnessMatrixFlag The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag The flag to set if compute the RHS
     */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
        ) ;


    /**
     * @brief This function provides a more general interface to the element.
     * @details It is designed so that rLHSvariables and rRHSvariables are passed to the element thus telling what is the desired output
     * @param rLeftHandSideMatrix container with the output Left Hand Side matrix
     * @param rRightHandSideVector container for the desired RHS output
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This is called during the assembling process in order to calculate the elemental left hand side matrix only
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief This is called during the assembling process in order to calculate the elemental right hand side vector only
      * @param rRightHandSideVector the elemental right hand side vector
      * @param rCurrentProcessInfo the current process info instance
      */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This functions updates the kinematics variables
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param PointNumber The integration point considered
     */
    void CalculateKinematicVariables(
        KinematicVariables& rThisKinematicVariables,
        const IndexType PointNumber,
        const GeometryType::IntegrationMethod& rIntegrationMethod
        ) ;

    /**
     * Calculation of the Deformation Matrix B
     * @param rB The deformation matrix
     * @param rDN_DX The derivatives of the shape functions
     * @param IntegrationPoints The array containing the integration points
     * @param PointNumber The integration point considered
     */
    virtual void CalculateB(
        Matrix& rB,
        const Matrix& rDN_DX,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const IndexType PointNumber
        ) const;

    double GetBodyForce( 
        const GeometryType::IntegrationPointsArrayType& rIntegrationPoints, 
        const IndexType PointNumber);

    /**
     * @brief Calculation of the Material Stiffness Matrix. Km = B^T * D *B
     * @param rLeftHandSideMatrix The local LHS of the element
     * @param B The deformationmmatrix
     * @param D The constitutive matrix
     * @param IntegrationWeight The integration weight of the corresponding Gauss point
     */
    void CalculateAndAddKm(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& B,
    const double IntegrationWeight
    ) const;


    /**
      * @brief This is called during the assembling process in order to calculate the elemental mass matrix
      * @param rMassMatrix The elemental mass matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This method computes directly the lumped mass vector
     * @param rMassMatrix The lumped mass vector
     */
    void CalculateLumpedMassVector(VectorType& rMassVector) const;

    // void CalculateAndAddMm(
    // MatrixType& rLeftHandSideMatrix,
    // const Vector& N,
    // const double IntegrationWeight,
    // const Matrix& D
    // ) const;

    void CalculateAndAddResidualVector(
    VectorType& rRightHandSideVector,
    const KinematicVariables& rThisKinematicVariables,
    const ProcessInfo& rCurrentProcessInfo,
    const double rBodyForce,
    const Vector& rPressureGradientVector,
    const double IntegrationWeight
    ) const;

    double GetIntegrationWeight(
    const GeometryType::IntegrationPointsArrayType& rThisIntegrationPoints,
    const IndexType PointNumber,
    const double detJ
    ) const;

    void CalculateAndAddExtForceContribution(
    const Vector& rN,
    const ProcessInfo& rCurrentProcessInfo,
    const double rBodyForce,
    VectorType& rRightHandSideVector,
    const double Weight
    ) const;


    /**
     * @brief This functions updates the data structure passed to the CL
     * @param rThisConstitutiveVariables The constitutive variables
     * @param rValues The CL parameters
     */
    void SetConstitutiveVariables(
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues
        ) ;

    /**
     * @brief This functions updates the constitutive variables
     * @param rThisConstitutiveVariables The constitutive variables
     * @param rValues The CL parameters
      */
    void CalculateConstitutiveVariables(
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues
        ) ;
        

    /**
     * @brief Calculation of the equivalent deformation gradient
     * @param rF The deformation gradient F
     * @param StrainVector The strain tensor (Voigt notation)
     
    */

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

    // void save(Serializer& rSerializer) const override;

    // void load(Serializer& rSerializer) override;

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    /// Assignment operator.
    //AcousticElement& operator=(const AcousticElement& rOther);
    /// Copy constructor.
    //AcousticElement(const AcousticElement& rOther);
    ///@}

}; // Class AcousticElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_ACOUSTIC_ELEMENT_H_INCLUDED  defined
