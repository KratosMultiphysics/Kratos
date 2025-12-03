//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andrea Gorgi
//

#pragma once


// System includes
#include "includes/define.h"
#include "includes/condition.h"

// External includes

// Project includes
#include "iga_application_variables.h"
#include "includes/constitutive_law.h"

namespace Kratos
{
/// Condition for penalty support condition
class KRATOS_API(IGA_APPLICATION) SupportSolidCondition
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer definition of SupportSolidCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SupportSolidCondition);

    /// Size types
    using SizeType = std::size_t;
    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{
    
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /// Constructor with Id and geometry
    SupportSolidCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {};

    /// Constructor with Id, geometry and property
    SupportSolidCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {};

    /// Default constructor
    SupportSolidCondition() : Condition()
    {};

    /// Destructor
    virtual ~SupportSolidCondition() override
    {};

    ///@}
    ///@name Life Cycle
    ///@{

    /// Create with Id, pointer to geometry and pointer to property
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive<SupportSolidCondition>(
            NewId, pGeom, pProperties);
    };

    /// Create with Id, pointer to geometry and pointer to property
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive<SupportSolidCondition>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    };

    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief This is called during the assembling process in order
    *        to calculate the condition right hand side matrix
    * @param rLeftHandSideMatrix the condition right hand side matrix
    * @param rCurrentProcessInfo the current process info
    */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
    * @brief This is called during the assembling process in order
    *        to calculate the condition left hand side matrix
    * @param rLeftHandSideMatrix the condition left hand side matrix
    * @param rCurrentProcessInfo the current process info
    */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This function provides a more general interface to the element.
     * @details It is designed so that rLHSvariables and rRHSvariables are
     *          passed to the element thus telling what is the desired output
     * @param rLeftHandSideMatrix container with the output Left Hand Side matrix
     * @param rRightHandSideVector container for the desired RHS output
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

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
    * @brief Sets on rConditionDofList the degrees of freedom of the considered element geometry
    * @param rElementalDofList The vector containing the dof of the element
    * @param rCurrentProcessInfo The current process info instance
    */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const override;

    /**
     * @brief Get the solution coefficient at the previous time step in the two-dimensional case.
     * 
     * @param rValues solution coefficients at the previous time step
     */
    void GetSolutionCoefficientVector(
        Vector& rValues) const;


    /**
     * @brief Calculate the B matrix for the element in the two-dimensional case.
     * 
     * @param rB B matrix to be calculated
     * @param r_DN_DX The shape function derivatives in the global coordinate system
     */
    void CalculateB(
        Matrix& rB,
        Matrix& r_DN_DX) const;

    /**
     * @brief Fills the delta position matrix as CurrentPosition - InitialPosition.
     * @param rGeometry Geometry providing nodal coordinates.
     * @param rDeltaPosition Output matrix sized [num_nodes x 3].
     */
    void CalculateDeltaPositionMatrix(
        const GeometryType& rGeometry,
        Matrix& rDeltaPosition);
    ///@}
    ///@name Check
    ///@{

    /// Performs check if Penalty factor is provided.
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "\"SupportSolidCondition\" #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "\"SupportSolidCondition\" #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }


    void AnalyticalConstitutiveMatrix(Matrix& rD, Vector& rStrain)
    {
        double nu = this->GetProperties().GetValue(POISSON_RATIO);
        double E = this->GetProperties().GetValue(YOUNG_MODULUS);
        double sig_y = this->GetProperties().GetValue(YIELD_STRESS);

        double eps_xx =rStrain[0];
        double gamma_xy = rStrain[2];
        double eps_yy = rStrain[1];

        rD(0,0) = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            E*(nu - 1)/((nu + 1)*(2*nu - 1))
         )
         : (
            (-E*nu*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2) + (1.0/3.0)*E*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2) - 4*sig_y*std::pow(2*eps_xx - eps_yy, 2)*(nu + 1)*(2*nu - 1)*std::sqrt(27*std::pow(gamma_xy, 2) + 6*std::pow(eps_xx - 2*eps_yy, 2) + 6*std::pow(eps_xx + eps_yy, 2) + 6*std::pow(2*eps_xx - eps_yy, 2)) + (4.0/3.0)*std::sqrt(3)*sig_y*(nu + 1)*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 3.0/2.0))/((nu + 1)*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2))
         ));
         rD(0,1) = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            -E*nu/((nu + 1)*(2*nu - 1))
         )
         : (
            (-E*nu*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2) + (1.0/3.0)*E*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2) + 4*sig_y*(eps_xx - 2*eps_yy)*(2*eps_xx - eps_yy)*(nu + 1)*(2*nu - 1)*std::sqrt(27*std::pow(gamma_xy, 2) + 6*std::pow(eps_xx - 2*eps_yy, 2) + 6*std::pow(eps_xx + eps_yy, 2) + 6*std::pow(2*eps_xx - eps_yy, 2)) - 2.0/3.0*std::sqrt(3)*sig_y*(nu + 1)*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 3.0/2.0))/((nu + 1)*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2))
         ));
         rD(0,2) = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            0
         )
         : (
            2*gamma_xy*sig_y*(-2*eps_xx + eps_yy)/std::pow(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2), 3.0/2.0)
         ));
         rD(1,0) = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            -E*nu/((nu + 1)*(2*nu - 1))
         )
         : (
            (-E*nu*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2) + (1.0/3.0)*E*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2) + 4*sig_y*(eps_xx - 2*eps_yy)*(2*eps_xx - eps_yy)*(nu + 1)*(2*nu - 1)*std::sqrt(27*std::pow(gamma_xy, 2) + 6*std::pow(eps_xx - 2*eps_yy, 2) + 6*std::pow(eps_xx + eps_yy, 2) + 6*std::pow(2*eps_xx - eps_yy, 2)) - 2.0/3.0*std::sqrt(3)*sig_y*(nu + 1)*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 3.0/2.0))/((nu + 1)*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2))
         ));
         rD(1,1) = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            E*(nu - 1)/((nu + 1)*(2*nu - 1))
         )
         : (
            (-E*nu*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2) + (1.0/3.0)*E*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2) - 4*sig_y*std::pow(eps_xx - 2*eps_yy, 2)*(nu + 1)*(2*nu - 1)*std::sqrt(27*std::pow(gamma_xy, 2) + 6*std::pow(eps_xx - 2*eps_yy, 2) + 6*std::pow(eps_xx + eps_yy, 2) + 6*std::pow(2*eps_xx - eps_yy, 2)) + (4.0/3.0)*std::sqrt(3)*sig_y*(nu + 1)*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 3.0/2.0))/((nu + 1)*(2*nu - 1)*std::pow(9*std::pow(gamma_xy, 2) + 2*std::pow(eps_xx - 2*eps_yy, 2) + 2*std::pow(eps_xx + eps_yy, 2) + 2*std::pow(2*eps_xx - eps_yy, 2), 2))
         ));
         rD(1,2) = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            0
         )
         : (
            2*gamma_xy*sig_y*(eps_xx - 2*eps_yy)/std::pow(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2), 3.0/2.0)
         ));
         rD(2,0) = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            0
         )
         : (
            2*gamma_xy*sig_y*(-2*eps_xx + eps_yy)/std::pow(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2), 3.0/2.0)
         ));
         rD(2,1) = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            0
         )
         : (
            2*gamma_xy*sig_y*(eps_xx - 2*eps_yy)/std::pow(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2), 3.0/2.0)
         ));
         rD(2,2) = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            (1.0/2.0)*E/(nu + 1)
         )
         : (
            4*sig_y*(std::pow(eps_xx, 2) - eps_xx*eps_yy + std::pow(eps_yy, 2))/std::pow(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2), 3.0/2.0)
         ));
    }

    void AnalyticalStress(Vector& rStressVector, Vector& rStrain)
    {
        double nu = this->GetProperties().GetValue(POISSON_RATIO);
        double E = this->GetProperties().GetValue(YOUNG_MODULUS);
        double sig_y = this->GetProperties().GetValue(YIELD_STRESS);

        double eps_xx =rStrain[0];
        double gamma_xy = rStrain[2];
        double eps_yy = rStrain[1];
        
        rStressVector[0] = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            E*(eps_xx*nu - eps_xx - eps_yy*nu)/(2*std::pow(nu, 2) + nu - 1)
         )
         : (
            (1.0/3.0)*(-E*eps_xx*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2)) - E*eps_yy*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2)) + 8*eps_xx*nu*sig_y - 4*eps_xx*sig_y - 4*eps_yy*nu*sig_y + 2*eps_yy*sig_y)/((2*nu - 1)*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2)))
         ));
         rStressVector[1] = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            E*(-eps_xx*nu + eps_yy*nu - eps_yy)/(2*std::pow(nu, 2) + nu - 1)
         )
         : (
            (1.0/3.0)*(-E*eps_xx*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2)) - E*eps_yy*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2)) - 4*eps_xx*nu*sig_y + 2*eps_xx*sig_y + 8*eps_yy*nu*sig_y - 4*eps_yy*sig_y)/((2*nu - 1)*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2)))
         ));
         rStressVector[2] = ((sig_y >= (1.0/2.0)*E*std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))/(nu + 1)) ? (
            (1.0/2.0)*E*gamma_xy/(nu + 1)
         )
         : (
            gamma_xy*sig_y/std::sqrt(4*std::pow(eps_xx, 2) - 4*eps_xx*eps_yy + 4*std::pow(eps_yy, 2) + 3*std::pow(gamma_xy, 2))
         ));
    }

    ///@}

protected:


/**
 * Internal variables used in the constitutive calculations
 */
struct ConstitutiveVariables
{
    ConstitutiveLaw::StrainVectorType StrainVector;
    ConstitutiveLaw::StressVectorType StressVector;
    ConstitutiveLaw::VoigtSizeMatrixType D;

    /**
     * The default constructor
     * @param StrainSize The size of the strain vector in Voigt notation
     */
    ConstitutiveVariables(const SizeType StrainSize)
    {
        if (StrainVector.size() != StrainSize)
            StrainVector.resize(StrainSize);

        if (StressVector.size() != StrainSize)
            StressVector.resize(StrainSize);

        if (D.size1() != StrainSize || D.size2() != StrainSize)
            D.resize(StrainSize, StrainSize);

        noalias(StrainVector) = ZeroVector(StrainSize);
        noalias(StressVector) = ZeroVector(StrainSize);
        noalias(D)            = ZeroMatrix(StrainSize, StrainSize);
    }
};

///@name Protected static Member Variables
///@{
void InitializeMaterial();

void InitializeMemberVariables();

void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

/**
 * @brief Compute the constitutive law response for the given strain vector.
 * 
 * @param matSize 
 * @param rStrain 
 * @param rValues 
 * @param rConstitutiVariables 
 */
void ApplyConstitutiveLaw(
        SizeType matSize, 
        Vector& rStrain, 
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveVariables& rConstitutiVariables);

//@}
///@name Protected member Variables
///@{
ConstitutiveLaw::Pointer mpConstitutiveLaw; /// The pointer containing the constitutive laws
unsigned int mDim; /// The dimension of the condition 
double mCharacteristicGeometryLength = 0.0; /// Characteristic length used by the constitutive law
Matrix mOldConstitutiveMatrix; /// Constitutive matrix cached at InitializeSolutionStep
Matrix mDBOld; /// Stored product D*B computed at InitializeSolutionStep

///@}

private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
    }

    ///@}

}; // Class SupportPenaltyLaplacianCondition

}  // namespace Kratos.
