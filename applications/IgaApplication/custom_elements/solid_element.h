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

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/math_utils.h"
#include "includes/variables.h"

// Application includes
#include "iga_application_variables.h"

namespace Kratos
{
////@name Kratos Globals
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
class KRATOS_API(IGA_APPLICATION) SolidElement : public Element
{
protected:
/**
    /**
     * Internal variables used in the kinematic calculations
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

public:

    ///@name Type Definitions
    ///@{

    /// Counted pointer of SolidElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SolidElement);

    // static constexpr std::size_t NumNodes = TDim + 1;

    ///@}
    ///@name Life Cycle
    ///@{

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /// Default constructor.
    SolidElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry);

    SolidElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~SolidElement();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override;


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
    * @brief This is called during the assembling process in order
    *        to calculate the condition left hand side matrix
    * @param rLeftHandSideMatrix the condition left hand side matrix
    * @param rCurrentProcessInfo the current process info
    */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix, 
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
    * @brief This is called during the assembling process in order
    *        to calculate the condition right hand side matrix
    * @param rLeftHandSideMatrix the condition right hand side matrix
    * @param rCurrentProcessInfo the current process info
    */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector, 
        const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(
        EquationIdVectorType& rResult, 
        const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(
        DofsVectorType& ElementalDofList, 
        const ProcessInfo& CurrentProcessInfo) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    SolidElement() : Element()
    {
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{
    
    IntegrationMethod GetIntegrationMethod() const override;

    ///@}
    ///@name Input and output
    ///@{


    ///@}
    ///@name Friends
    ///@{
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "SolidElement #" << Id();
        return buffer.str();
    }

    ///@}
/**
    * @brief Calculate a double Variable on the Element Constitutive Law
    * @param rVariable The variable we want to get
    * @param rValues The values obtained int the integration points
    * @param rCurrentProcessInfo the current process info instance
    */
    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
    * @brief Calculate a Vector Variable on the Element Constitutive Law
    * @param rVariable The variable we want to get
    * @param rValues The values obtained int the integration points
    * @param rCurrentProcessInfo the current process info instance
    */
    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3>>& rVariable,
        std::vector<array_1d<double, 3>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    void CalculateB(
        Matrix& rB,
        Matrix& r_DN_DX) const;

    void GetValuesVector(
        Vector& rValues) const;

protected:
    ///@name Protected static Member Variables
    ///@{

    void InitializeMaterial();


    ///@}
    ///@name Protected member Variables
    ///@{

    ConstitutiveLaw::Pointer mpConstitutiveLaw; /// The pointer containing the constitutive laws


    ///@}
    ///@name Protected Operators
    ///@{


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

    // Protected default constructor necessary for serialization

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    // std::vector<std::size_t> GetSurrogateFacesIds();

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class SolidElement

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.
