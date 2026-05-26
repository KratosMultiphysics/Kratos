//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Ricky Aristio
//                   Maram Alkhlaifat
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/math_utils.h"

// Application includes
#include "iga_application_variables.h"

namespace Kratos
{
///@}
///@name Kratos Classes
///@{
/**
 * @class Shell6pElement
 * @ingroup IgaApplication
 * @brief This class defines a Reissner–Mindlin shell element for isogeometric analysis based on degenerated solid theory.
 * Theory behind this implementation can be found in Benson et al. "Isogeometric shell analysis: The Reissner–Mindlin shell" (2010),
 * DOI:10.1016/j.cma.2009.05.011 
 * 
 * @details The element uses a curvilinear geometry representation from the IGA framework and employing six kinematic parameters per control point (three translational and three rotational degrees of freedom). 
 * The present implementation is modified by the formulations in Benson et al, in terms of using displacements as the primal variable, instead of velocities.
 * Moreover, it is enhanced for geometrically nonlinear analysis by applying co-rotational method.
 */

class Shell6pElement
    : public Element
{
protected:

    /// @brief Internal structs
    struct KinematicVariables
    {
        array_1d<double, 3> MetricCovariant;
        array_1d<double, 3> CurvatureCovariant;
        array_1d<double, 3> BaseVector1;
        array_1d<double, 3> BaseVector2;
        array_1d<double, 3> NormalVector;
        array_1d<double, 3> NormalVectorTilde;
        double DifferentialArea;
        KinematicVariables(std::size_t Dimension);
    };

    struct ConstitutiveVariables
    {
        Vector StrainVector;
        Vector StressVector;
        Matrix ConstitutiveMatrix;
        ConstitutiveVariables(std::size_t StrainSize);
    };

public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of Shell6pElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(Shell6pElement);

    using IndexType = std::size_t;
    using GeometryType = Geometry<Node>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor using an array of nodes
    Shell6pElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry);

    /// Constructor using an array of nodes with properties
    Shell6pElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties);

    /// Default constructor necessary for serialization
    Shell6pElement() = default;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Create with Id, pointer to geometry and pointer to property
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override;

    /// Create with Id, pointer to geometry and pointer to property
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override;

    ///@}
    ///@name Operations
    ///@{

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

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

    /**
    * @brief This is called during the assembling process in order to calculate the elemental mass matrix
    * @param rMassMatrix The elemental mass matrix
    * @param rCurrentProcessInfo The current process info instance
    */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
    * @brief This is called during the assembling process in order to calculate the elemental damping matrix
    * @param rDampingMatrix The elemental damping matrix
    * @param rCurrentProcessInfo The current process info instance
    */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

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

    void GetValuesVector(
        Vector& rValues,
        int Step) const override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Check
    ///@{

    /**
    * This function provides the place to perform checks on the completeness of the input.
    * It is designed to be called only once (or anyway, not often) typically at the beginning
    * of the calculations, so to verify that nothing is missing from the input
    * or that no common error is found.
    * @param rCurrentProcessInfo
    */
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Shell6pElement #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Shell6pElement #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    // Components of the metric coefficient tensor on the contravariant basis
    std::vector<array_1d<double, 3>> m_A_ab_covariant_vector;
    // Components of the curvature coefficient tensor on the contravariant basis
    std::vector<array_1d<double, 3>> m_B_ab_covariant_vector;

    // Determinant of the geometrical Jacobian.
    Vector m_dA_vector;

    /* Transformation the strain tensor from the curvilinear system
    *  to the local cartesian in voigt notation including a 2 in the
    *  shear part. */
    std::vector<Matrix> m_T_vector;

    /// The vector containing the constitutive laws for all integration points.
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    // curvilinear coordinate zeta (theta3)
    double mZeta;

    /// @brief Informations regarding the Gauss-quadrature in thickness direction
    struct GaussQuadratureThickness
    {
        unsigned int num_GP_thickness;
        Vector integration_weight_thickness;
        Vector zeta;

        // The default constructor
        GaussQuadratureThickness(){}
        // constructor
        GaussQuadratureThickness(const unsigned int& rNumGPThickness)
        {
            num_GP_thickness = rNumGPThickness;
            integration_weight_thickness = ZeroVector(rNumGPThickness);
            zeta = ZeroVector(rNumGPThickness);

            if (rNumGPThickness == 2)
            {
                integration_weight_thickness(0) = 1.0;
                zeta(0) = -sqrt(1.0 / 3.0);
                integration_weight_thickness(1) = 1.0;
                zeta(1) = sqrt(1.0 / 3.0);
            }
            else if (rNumGPThickness == 3)
            {
                integration_weight_thickness(0) = 5.0 / 9.0;
                zeta(0) = -sqrt(3.0 / 5.0);
                integration_weight_thickness(1) = 8.0/9.0;
                zeta(1) = 0.0;
                integration_weight_thickness(2) = 5.0 / 9.0;
                zeta(2) = sqrt(3.0 / 5.0);
            }
            else
            {
                KRATOS_ERROR << "Desired number of Gauss-Points unlogical or not implemented. You can choose 3 Gauss-Points." << std::endl;
            }   
        }
    };

    // Specified the number of Gauss-Points over the thickness 
    GaussQuadratureThickness mGaussIntegrationThickness = GaussQuadratureThickness(2);

    ///@}
    ///@name Operations
    ///@{

    /// Calculates LHS and RHS dependent on flags
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    );

    /// Initialize Operations
    void InitializeMaterial();

    void CalculateKinematics(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables) const;

    // Compute the transformation matrix from local to global cartesian coordinates
    void CalculateTransformationFromLocalToGlobalCartesian(
        const KinematicVariables& rKinematicVariables,
        Matrix& rTransformationMatrix) const; 
    
    void CalculateBOperator(
        const IndexType IntegrationPointIndex,
        Matrix& rB,
        double zeta,
        Matrix& J_inv,
        Matrix& dn,
        const KinematicVariables& rActualKinematic) const;

    void CalculateBGeometric(
        const IndexType IntegrationPointIndex,
        Matrix& rB,
        double zeta,
        Matrix& J_inv,
        Matrix& dn,
        const KinematicVariables& rActualKinematic) const;

    void CalculateNormalVectorDerivatives(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables,
        Matrix& DerivativeNormalMatrix) const;

    void CalculateBDrilling(
        const IndexType IntegrationPointIndex,
        Matrix& rBd,
        Matrix& J_inv,
        const KinematicVariables& rActualKinematic) const;

    void CalculateStressMatrix(
        array_1d<double, 6> stress_vector,
        Matrix& stress_matrix
    ) const;

    /**
    * This functions updates the constitutive variables
    * @param rActualMetric: The actual metric
    * @param rThisConstitutiveVariables: The constitutive variables to be calculated
    * @param rValues: The CL parameters
    * @param ThisStressMeasure: The stress measure considered
    */
    void CalculateConstitutiveVariables(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rActualMetric,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    ) const;

    /**
     * @brief This method gets a value directly from the CL
     * @details Avoids code repetition
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained at the integration points
     * @tparam TType The variable type
     */
    template<class TType>
    void GetValueOnConstitutiveLaw(
        const Variable<TType>& rVariable,
        std::vector<TType>& rOutput
        )
    {
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

        for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {
            mConstitutiveLawVector[point_number]->GetValue(rVariable, rOutput[point_number]);
        }
    }

    ///@}
    ///@name Geometrical Functions
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
        rSerializer.save("A_ab_covariant_vector", m_A_ab_covariant_vector);
        rSerializer.save("B_ab_covariant_vector", m_B_ab_covariant_vector);
        rSerializer.save("dA_vector", m_dA_vector);
        rSerializer.save("T_vector", m_T_vector);
        rSerializer.save("constitutive_law_vector", mConstitutiveLawVector);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
        rSerializer.load("A_ab_covariant_vector", m_A_ab_covariant_vector);
        rSerializer.load("B_ab_covariant_vector", m_B_ab_covariant_vector);
        rSerializer.load("dA_vector", m_dA_vector);
        rSerializer.load("T_vector", m_T_vector);
        rSerializer.load("constitutive_law_vector", mConstitutiveLawVector);
    }

    ///@}

};     // Class Shell6pElement
///@}

}  // namespace Kratos.