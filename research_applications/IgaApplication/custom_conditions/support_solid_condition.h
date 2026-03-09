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
