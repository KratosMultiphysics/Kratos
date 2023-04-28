// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//


#if !defined(KRATOS_SMALL_DISPLACEMENT_NON_LOCAL_HEX_H_INCLUDED )
#define  KRATOS_SMALL_DISPLACEMENT_NON_LOCAL_HEX_H_INCLUDED

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
 * @class SmallDisplacementNonLocalHex
 * @ingroup StructuralMechanicsApplication
 * @brief Small displacement element for 2D and 3D geometries.
 * @details Implements a small displacement definition for structural analysis. This works for arbitrary geometries in 2D and 3D
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SmallDisplacementNonLocalHex
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
        Vector DuNL;
        Vector DNLu;
        double DNLNL;
        double Local_Variable_GP;
        double NonLocal_Variable_GP;
        
        /**
         * The default constructor
         * @param StrainSize The size of the strain vector in Voigt notation
         */
        NonLocalConstitutiveVariables(const SizeType StrainSize)
        {
            if (DuNL.size() != StrainSize)
                DuNL.resize(StrainSize);
            if (DNLu.size() != StrainSize)
                DNLu.resize(StrainSize);

            noalias(DuNL)        = ZeroVector(StrainSize);
            noalias(DNLu)        = ZeroVector(StrainSize);
            DNLNL                = 0.0;
            Local_Variable_GP    = 0.0;
            NonLocal_Variable_GP = 0.0;
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

    /// Counted pointer of SmallDisplacementNonLocalHex
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SmallDisplacementNonLocalHex);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SmallDisplacementNonLocalHex(IndexType NewId, GeometryType::Pointer pGeometry);
    SmallDisplacementNonLocalHex(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    // Copy constructor
    SmallDisplacementNonLocalHex(SmallDisplacementNonLocalHex const& rOther)
        :BaseType(rOther)
    {};

    /// Destructor.
    ~SmallDisplacementNonLocalHex() override;

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

    SmallDisplacementNonLocalHex() : BaseSolidElement()
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
        const Variable<double>& rThisVariable,
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
    void CalculateAndAddKmud(
        Matrix& rStiffnessMatrixKud,
        const Matrix& B,
        const Vector& DuNL,
        const Vector& N,
        const double IntegrationWeight
        ) const;

    /**
     * @brief Calculation of element stiffness KNLu
     * @param B The deformation matrix
     * @param DNLu material tangent
     * @param N shape function values
     */
    void CalculateAndAddKmdu(
        Matrix& rStiffnessMatrixKdu,
        const Matrix& B,
        const Vector& DNLu,
        const Vector& N,
        const double IntegrationWeight
        ) const;

    /**
     * @brief Calculation of element stiffness KNLNL
     * @param B The deformation matrix
     * @param DNLNL material tangent
     * @param N shape function values
     */
    void CalculateAndAddKmdd(
        Matrix& rStiffnessMatrixKdd,
        const Matrix& DN_DX,
        const double& DNLNL,
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
        const Matrix& Kud,
        const Matrix& Kdu,
        const Matrix& Kdd,
        const Vector& rFu,
        const Vector& rFNL 
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
        Vector& rResidualNonLocalVector, 
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
    //SmallDisplacementNonLocalHex& operator=(const SmallDisplacementNonLocalHex& rOther);
    /// Copy constructor.
    //SmallDisplacementNonLocalHex(const SmallDisplacementNonLocalHex& rOther);
    ///@}

}; // Class SmallDisplacementNonLocalHex

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_SMALL_DISPLACEMENT_H_INCLUDED  defined
