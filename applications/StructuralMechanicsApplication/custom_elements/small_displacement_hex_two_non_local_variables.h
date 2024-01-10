// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Athira Vadakkekkara
//


#if !defined(KRATOS_SMALL_DISPLACEMENT_HEX_TWO_NON_LOCAL_VARIABLES_H_INCLUDED )
#define  KRATOS_SMALL_DISPLACEMENT_HEX_TWO_NON_LOCAL_VARIABLES_H_INCLUDED

// System includes


// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/base_solid_element.h"
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
 * @class SmallDisplacementHexTwoNonLocalVariables
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details Equivalent strains in tension and compression Et and Ec are new DOFS and are defined as new kratos variables. In the other nonlocal  models with Eq strain or damage as the nonlocal variable,
 * additional kratos variable(string) -NONLOCAL_VARIABLE_NAME has to be passed in material properties- to set either eq strain or damage as the
 * nonlocal variable. To avoid conflicts, this is removed in this class and instead two separate kratos variables- equivalent strain in tension and compression are defined
 * and set as the DOFS. So the material law should also take this in consideration.
 * @author Athira Vadakkekkara
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SmallDisplacementHexTwoNonLocalVariables
    : public BaseSolidElement
{
public:

    /**
     * Internal variables (material tangents) used in the kinematic calculations of non local formulation
     * DuNL
     * DNLu
     * DNLNL
     */
    struct NonLocalConstitutiveVariables
    {
        Vector DuNL1;
        Vector DuNL2;
        Vector DNL1u;
        Vector DNL2u;
        Vector Local_Variables_GP;
        Vector NonLocal_Variables_GP;
        double DNL1NL1;
        double DNL2NL2;

        /**
         * The default constructor
         * @param StrainSize The size of the strain vector in Voigt notation
         */
        NonLocalConstitutiveVariables(const SizeType StrainSize)
        {
            if (DuNL1.size() != StrainSize)
                DuNL1.resize(StrainSize);
            if (DuNL2.size() != StrainSize)
                DuNL2.resize(StrainSize);
            if (DNL1u.size() != StrainSize)
                DNL1u.resize(StrainSize);
            if (DNL2u.size() != StrainSize)
                DNL2u.resize(StrainSize);

            noalias(DuNL1)        = ZeroVector(StrainSize);
            noalias(DuNL2)        = ZeroVector(StrainSize);
            noalias(DNL1u)        = ZeroVector(StrainSize);
            noalias(DNL2u)        = ZeroVector(StrainSize);
            Local_Variables_GP    = ZeroVector(2);
            NonLocal_Variables_GP = ZeroVector(2);
            DNL1NL1               = 0.0;
            DNL2NL2               = 0.0;
        }
    };
    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// The base element type
    typedef BaseSolidElement BaseType;

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    /// Counted pointer of SmallDisplacementHexTwoNonLocalVariables
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SmallDisplacementHexTwoNonLocalVariables);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SmallDisplacementHexTwoNonLocalVariables(IndexType NewId, GeometryType::Pointer pGeometry);
    SmallDisplacementHexTwoNonLocalVariables(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    // Copy constructor
    SmallDisplacementHexTwoNonLocalVariables(SmallDisplacementHexTwoNonLocalVariables const& rOther)
        :BaseType(rOther)
    {};

    /// Destructor.
    ~SmallDisplacementHexTwoNonLocalVariables() override;

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

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
     * @brief Sets on rResult the ID's of the element degrees of freedom
     * @param rResult The vector containing the equation id
     * @param rCurrentProcessInfo The current process info instance
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /**
     * @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
     * @param rElementalDofList The vector containing the dof of the element
     * @param rCurrentProcessInfo The current process info instance
     */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /**
     * @brief Sets on rValues the nodal displacements
     * @param rValues The values of displacements
     * @param Step The step to be computed
     */
    void GetValuesVector(
        Vector& rValues,
        int Step = 0
        ) const override;

    /**
     * @brief Sets on rValues the nodal velocities
     * @param rValues The values of velocities
     * @param Step The step to be computed
     */
    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step = 0
        ) const override;

    /**
     * @brief Sets on rValues the nodal accelerations
     * @param rValues The values of accelerations
     * @param Step The step to be computed
     */
    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step = 0
        ) const override;

    /**
    * Called at the end of eahc solution step
    * @param rCurrentProcessInfo the current process info instance
    */
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;
    /**
     * @brief Calculate a Vector Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{

    const Parameters GetSpecifications() const override;

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Small Displacement Solid Element #" << Id() << "\nConstitutive law: " << BaseType::mConstitutiveLawVector[0]->Info();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Small Displacement Solid Element #" << Id() << "\nConstitutive law: " << BaseType::mConstitutiveLawVector[0]->Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
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

    ///@}
    ///@name Protected Operators
    ///@{

    SmallDisplacementHexTwoNonLocalVariables() : BaseSolidElement()
    {
    }

     /**
     * @brief This method returns if the element provides the strain
     */
    bool UseElementProvidedStrain() const override;

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
        ) override;

    /**
     * @brief This functions calculates non local constitutive variables
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param PointNumber The integration point considered
     */
    void CalculateNonLocalConstitutiveVariables(
        KinematicVariables& rThisKinematicVariables,
        NonLocalConstitutiveVariables& rThisNonLocalConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const Variable<Vector>& rThisVariable,
        const IndexType PointNumber
        );

    /**
     * This functions updates the constitutive variables
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param rThisConstitutiveVariables The constitutive variables
     * @param rThisNonLocalConstitutiveVariables the non local constitutive variables
     * @param rValues The CL parameters
     * @param PointNumber The integration point considered
     * @param IntegrationPoints The list of integration points
     * @param ThisStressMeasure The stress measure considered
     */
    void CalculateAllConstitutiveVariables(
        KinematicVariables& rThisKinematicVariables,
        ConstitutiveVariables& rThisConstitutiveVariables,
        NonLocalConstitutiveVariables& rThisNonLocalConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const IndexType PointNumber,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
        );


    /**
     * @brief This functions updates the data structure passed to the CL
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param rThisConstitutiveVariables The constitutive variables
     * @param rValues The CL parameters
     * @param PointNumber The integration point considered
     * @param IntegrationPoints The list of integration points
     */
    void SetConstitutiveVariables(
        KinematicVariables& rThisKinematicVariables,
        ConstitutiveVariables& rThisConstitutiveVariables,
        NonLocalConstitutiveVariables& rThisNonLocalConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const IndexType PointNumber,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints
        ) ;

    /**
     * @brief This functions assembles D(Duu) from rThisConstitutiveVariables and DuNL, DNLu, DNLNL from rThisNonLocalConstitutiveVariables to get D
     * @param rThisConstitutiveVariables The constitutive variables
     * @param rThisNonLocalConstitutiveVariables The non local constitutive variables
     */
    void SetNonLocalConstitutiveVariables(
        const Matrix& rConstitutiveMatrix,
        ConstitutiveVariables& rThisConstitutiveVariables,
        NonLocalConstitutiveVariables& rThisNonLocalConstitutiveVariables
        ) const;

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

    /**
     * @brief Calculation of the equivalent deformation gradient
     * @param rF The deformation gradient F
     * @param StrainVector The strain tensor (Voigt notation)
     */
    virtual void ComputeEquivalentF(
        Matrix& rF,
        const Vector& StrainVector
        ) const;

    /**
     * @brief Calculation of element stiffness KuNL
     * @param B The deformation matrix
     * @param DuNL material tangent
     * @param N shape function values
     */
    void CalculateAndAddKuNL(
        Matrix& rStiffnessMatrixKuNL1,
        Matrix& rStiffnessMatrixKuNL2,
        NonLocalConstitutiveVariables& rThisNonLocalConstitutiveVariables,
        const Matrix& B,
        const Vector& N,
        const double IntegrationWeight
        ) const;

    /**
     * @brief Calculation of element stiffness KNLu
     * @param B The deformation matrix
     * @param DNLu material tangent
     * @param N shape function values
     */
    void CalculateAndAddKNLu(
        Matrix& rStiffnessMatrixKNL1u,
        Matrix& rStiffnessMatrixKNL2u,
        NonLocalConstitutiveVariables& rThisNonLocalConstitutiveVariables,
        const Matrix& B,
        const Vector& N,
        const double IntegrationWeight
        ) const;

    /**
     * @brief Calculation of element stiffness KNLNL
     * @param B The deformation matrix
     * @param DNLNL material tangent
     * @param N shape function values
     */
    void CalculateAndAddKNLNL(
        Matrix& rStiffnessMatrixKNL1NL1,
        Matrix& rStiffnessMatrixKNL2NL2,
        NonLocalConstitutiveVariables& rThisNonLocalConstitutiveVariables,
        const Matrix& DN_DX,
        const Vector& N,
        const double IntegrationWeight
        ) const;

    /**
     * @brief Assembling stiffness matrix and residual vector
     * @param LeftHandSideMatrix assembled element stiffness matrix
     * @param RightHandSideVector assembled residual vector
     */
    void AssembleRHSAndLHS(
        MatrixType& LeftHandSideMatrix,
        VectorType& RightHandSideVector,
        const bool StiffnessMatrixFlag,
        const bool ResidualVectorFlag,
        const Matrix& Kuu,
        const Matrix& KuNL1,
        const Matrix& KuNL2,
        const Matrix& KNL1u,
        const Matrix& KNL2u,
        const Matrix& KNL1NL1,
        const Matrix& KNL2NL2,
        const Vector& rFu,
        const Vector& rFNL1,
        const Vector& rFNL2
        );

    /**
     * @brief Assembling stiffness matrix and residual vector
     * @param LeftHandSideMatrix assembled element stiffness matrix
     * @param RightHandSideVector assembled residual vector
     */
    void CalculateAndAddResidualForceVector(
        VectorType& rRightHandSideVector,
        const KinematicVariables& rThisKinematicVariables,
        const ProcessInfo& rCurrentProcessInfo,
        const array_1d<double, 3>& rBodyForce,
        const Vector& rStressVector,
        const double IntegrationWeight) const;

    /**
     * @brief adding non local components of right hand side vector
     * @param rResidualNonLocalVector
     */
    void CalculateAndAddResidualNonLocalVector(
        Vector& rResidualNonLocalVector1,
        Vector& rResidualNonLocalVector2,
        const KinematicVariables& rThisKinematicVariables,
        const NonLocalConstitutiveVariables& rThisNonLocalConstitutiveVariables,
        const ProcessInfo& rCurrentProcessInfo,
        const double IntegrationWeight);
    /**
     * @brief
     * @param rKuu
     */
    void CalculateAndAddKuu(
        MatrixType& rKuu,
        const Matrix& B,
        const Matrix& D,
        const double IntegrationWeight) const;

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

    Matrix mConstitutiveMatrix;
    const Variable<Vector>*mpEquivalentStrainsTC;


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
    /// Assignment operator.
    //SmallDisplacementHexTwoNonLocalVariables& operator=(const SmallDisplacementHexTwoNonLocalVariables& rOther);
    /// Copy constructor.
    //SmallDisplacementHexTwoNonLocalVariables(const SmallDisplacementHexTwoNonLocalVariables& rOther);
    ///@}

}; // Class SmallDisplacementHexTwoNonLocalVariables

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_SMALL_DISPLACEMENT_H_INCLUDED  defined
