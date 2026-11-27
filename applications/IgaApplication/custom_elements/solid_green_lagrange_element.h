//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Thomas Lillieskold, Andrea Gorgi 
//

#pragma once


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/math_utils.h"
#include "includes/variables.h"
#include "custom_elements/solid_element.h"

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
class KRATOS_API(IGA_APPLICATION) SolidGreenLagrangeElement : public SolidElement
{
protected:
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

    /// Counted pointer of SolidGreenLagrangeElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SolidGreenLagrangeElement);

    // static constexpr std::size_t NumNodes = TDim + 1;

    ///@}
    ///@name Life Cycle
    ///@{

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /// Default constructor.
    SolidGreenLagrangeElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry);

    SolidGreenLagrangeElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties);
    
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

    /// Destructor.
    virtual ~SolidGreenLagrangeElement();

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

    SolidGreenLagrangeElement() : SolidElement()
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
        buffer << "SolidGreenLagrangeElement #" << Id();
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
     * @brief Get the solution coefficient at the previous time step in the two-dimensional case.
     * 
     * @param rValues solution coefficients at the previous time step
     */
    void GetSolutionCoefficientVector(
        Vector& rValues) const;

    
protected:
    /** @brief Calculate the deformation gradient for the element.
     * 
     * @param rDN_DX The shape function derivatives in the global coordinate system
     * @return The deformation gradient
     */
    Matrix CalculateDeformationGradient(const Matrix& rDN_DX) const;

    /** @brief Calculate the B matrix for the element.
     * 
     * @param rDeformationGradient The deformation gradient
     * @param rDN_DX The shape function derivatives in the global coordinate system
     * @param rB The B matrix to be calculated
     */
    void CalculateBMatrix(
        const Matrix& rDeformationGradient,
        const Matrix& rDN_DX,
        Matrix& rB) const;

    /** @brief Add the material contribution to the tangent matrix. */
    void CalculateAndAddKm(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& rB,
        const Matrix& rConstitutiveMatrix,
        const double IntegrationWeight) const;

    /** @brief Add the geometric contribution to the tangent matrix. */
    void CalculateAndAddKg(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& rDN_DX,
        const Vector& rStressVector,
        const double IntegrationWeight) const;

    ///@name Protected static Member Variables
    ///@{

    void InitializeMaterial();

    /** @brief Compute the constitutive response from the deformation gradient. */
    void ApplyConstitutiveLaw(
        const Matrix& rDeformationGradient,
        const double DeterminantF,
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveVariables& rConstitutiveVariables,
        const bool ComputeConstitutiveTensor,
        const ConstitutiveLaw::StressMeasure StressMeasure);


    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{
    
    /**
     * @brief Calculate the initial Jacobian matrix for the element.
     * 
     * @param rGeometry 
     * @param rJacobian 
     */
    void CalculateInitialJacobian(
        const GeometryType& rGeometry, Matrix& rJacobian) const
    {
        GeometryType::JacobiansType J0;
        rGeometry.Jacobian(J0,this->GetIntegrationMethod());

        switch (mDim) {
            case 2:
            {
                rJacobian.resize(2,2);
                rJacobian(0,0) = J0[0](0,0);
                rJacobian(0,1) = J0[0](0,1);
                rJacobian(1,0) = J0[0](1,0);
                rJacobian(1,1) = J0[0](1,1);
                return;
            }
            case 3:
            {
                rJacobian.resize(3,3);
                rJacobian(0,0) = J0[0](0,0);
                rJacobian(0,1) = J0[0](0,1);
                rJacobian(0,2) = J0[0](0,2);
                rJacobian(1,0) = J0[0](1,0);
                rJacobian(1,1) = J0[0](1,1);
                rJacobian(1,2) = J0[0](1,2);
                rJacobian(2,0) = J0[0](2,0);
                rJacobian(2,1) = J0[0](2,1);
                rJacobian(2,2) = J0[0](2,2);
                return;
            }
            default:
                KRATOS_ERROR << "Dimension not supported: " << mDim << std::endl;
        }
    }


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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SolidElement);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SolidElement);
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

}; // Class SolidGreenLagrangeElement

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.
