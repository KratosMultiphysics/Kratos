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
///@name Kratos Classes
///@{
/// Short class definition.
/** Reissner–Mindlin shell element for isogeometric B-Rep analysis based on degenerated solid theory.
    The element uses a curvilinear geometry representation from the IGA framework.
    Employing six kinematic parameters per control point (three translational and three rotational degrees of freedom). 

    The present implementation is inspired by the formulations in Benson et al. and Du et al. , but does not follow any 
    single reference verbatim; the kinematic assumptions, strain measures and numerical integration are adapted to the 
    Kratos IGA framework and to geometrically nonlinear analysis of thin to moderately thick CAD-based shell structures.  

    For further details, see:
     [1] Benson, D.J., Bazilevs, Y., Hsu, M.C., & Hughes, T.J.R. (2010).
     "Isogeometric shell analysis: The Reissner–Mindlin shell."
     Computer Methods in Applied Mechanics and Engineering, 199, 276-289

     [2] Du, X., Li, J., Wang, W., Zhao, G., Liu, Y., & Zhang, P. (2024).
     "Isogeometric Shape Optimization of Reissner–Mindlin Shell with Analytical Sensitivity 
      and Application to Cellular Sandwich Structures."
     Composite Structures, Elsevier.
*/

class Shell6pBbarElement
    : public Element
{
protected:

    /// @brief Internal structs
    struct KinematicVariables
    {
        array_1d<double, 3> a_ab_covariant;
        array_1d<double, 3> b_ab_covariant;
        array_1d<double, 3> a1;
        array_1d<double, 3> a2;
        array_1d<double, 3> a3;
        array_1d<double, 3> a3_tilde;
        double dA;
        KinematicVariables(std::size_t Dimension);
    };


    struct ConstitutiveVariables
    {
        Vector StrainVector;
        Vector StressVector;
        Matrix ConstitutiveMatrix;
        ConstitutiveVariables(std::size_t StrainSize);
    };


    struct SecondVariations
    {
        Matrix B11;
        Matrix B22;
        Matrix B12;
        SecondVariations(const int& mat_size);
    }; 

public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of Shell6pBbarElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(Shell6pBbarElement);

    using IndexType = std::size_t;
    using GeometryType = Geometry<Node>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor using an array of nodes
    Shell6pBbarElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {};

    /// Constructor using an array of nodes with properties
    Shell6pBbarElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {};

    /// Default constructor necessary for serialization
    Shell6pBbarElement() = default;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Create with Id, pointer to geometry and pointer to property
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive<Shell6pBbarElement>(
            NewId, pGeom, pProperties);
    };

    /// Create with Id, pointer to geometry and pointer to property
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive< Shell6pBbarElement >(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    };

    ///@}
    ///@name Obligatory KRATOS Operations
    ///@{

    /**
    * @brief This is called during the assembling process in order
    *        to calculate the condition right hand side matrix
    * @param rLeftHandSideMatrix the condition right hand side matrix
    * @param rCurrentProcessInfo the current process info
    */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const std::size_t number_of_nodes = GetGeometry().size();
        const std::size_t mat_size = number_of_nodes * 6;

        if (rRightHandSideVector.size() != mat_size)
            rRightHandSideVector.resize(mat_size);
        noalias(rRightHandSideVector) = ZeroVector(mat_size);

        MatrixType left_hand_side_matrix;

        CalculateAll(left_hand_side_matrix, rRightHandSideVector,
            rCurrentProcessInfo, false, true);
    }

    /**
    * @brief This is called during the assembling process in order
    *        to calculate the condition left hand side matrix
    * @param rLeftHandSideMatrix the condition left hand side matrix
    * @param rCurrentProcessInfo the current process info
    */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const std::size_t mat_size = GetGeometry().size() * 6;
        
        if (rLeftHandSideMatrix.size1() != mat_size)
            rLeftHandSideMatrix.resize(mat_size, mat_size);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

        VectorType right_hand_side_vector;
        CalculateAll(rLeftHandSideMatrix, right_hand_side_vector,
            rCurrentProcessInfo, true, false);
    }

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
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const std::size_t mat_size = GetGeometry().size() * 6;

        if (rRightHandSideVector.size() != mat_size)
            rRightHandSideVector.resize(mat_size);
        noalias(rRightHandSideVector) = ZeroVector(mat_size);

        if (rLeftHandSideMatrix.size1() != mat_size)
            rLeftHandSideMatrix.resize(mat_size, mat_size);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
            rCurrentProcessInfo, true, true);
    }

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

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

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

    ///@}
    ///@name Base Class Operations
    ///@{

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void GetValuesVector(
        Vector& rValues,
        int Step) const override;

    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step) const override;

    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step) const override;

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
        buffer << "Shell6pBbarElement #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Shell6pBbarElement #" << Id();
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

    // Computes transformation
    void CalculateTransformation(
        const KinematicVariables& rKinematicVariables,
        Matrix& rT) const;
    
    // Computes transformation for the stress tensor 
    void CalculateTransformationFromCovariantToCartesian(
        const KinematicVariables& rKinematicVariables,
        Matrix& rTCovToCar) const;

    void CalculateB(
        const IndexType IntegrationPointIndex,
        Matrix& rB,
        double zeta,
        Matrix& DN_De_Jn,
        Matrix& J_inv,
        Matrix& dn,
        const KinematicVariables& rActualKinematic) const;
     // Average strain through the thickness-B
    void CalculateMID(
        const IndexType IntegrationPointIndex,
        Matrix& rMID,
        Matrix& DN_De_Jn,
        Matrix& J_inv,
        Matrix& dn,
        const KinematicVariables& rActualKinematic) const;
    void CalculateN_A(
        const IndexType IntegrationPointIndex,
        Matrix& rN_sigma_A,
        const KinematicVariables& rActualKinematic) const;
    void CalculateN_B(
        const IndexType IntegrationPointIndex,
        Matrix& rN_sigma_B,
        const KinematicVariables& rActualKinematic) const;
    void CalculateM(
        const double IntegrationWeight,
        const double IntegrationWeight_zeta,
        Matrix& rN_sigma_A,
        Matrix& rN_sigma_B,
        Matrix& rM) const;
    void CalculateMID_bar(
        const double IntegrationWeight,
        const double IntegrationWeight_zeta,
        Matrix& rMID_bar,
        Matrix& rMID,
        Matrix& rN_sigma_A,
        Matrix& rN_sigma_B,
        Matrix& rM) const;

    void CalculateBGeometric(
        const IndexType IntegrationPointIndex,
        Matrix& rB,
        double zeta,
        Matrix& DN_De_Jn,
        Matrix& J_inv,
        Matrix& dn,
        const KinematicVariables& rActualKinematic) const;

    void CalculateJn(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables,
        double zeta,
        Matrix& DN_De_Jn,
        Matrix& J_inv,
        Matrix& dn,
        double& area) const;

    void CalculateBDrill(
        const IndexType IntegrationPointIndex,
        Matrix& rBd,
        Matrix& DN_De_Jn,
        const KinematicVariables& rActualKinematic) const;

    void CalculateStressMatrix(
        array_1d<double, 6> stress_vector,
        Matrix& stress_matrix
    ) const;

    void CalculateSecondVariationStrainCurvature(
        const IndexType IntegrationPointIndex,
        SecondVariations& rSecondVariationsStrain,
        SecondVariations& rSecondVariationsCurvature,
        const KinematicVariables& rActualKinematic) const;

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

     void CalculateAndAddK(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& rKm,
        const Matrix& rKd,                                                                                                               
        const double IntegrationWeight,
        const double IntegrationWeight_zeta) const;

     void CalculateAndAddKm(
        MatrixType& rKm,
        const Matrix& rB,
        const Matrix& rD,
        const Matrix& rMID
        const Matrix& rMID_bar  ) const;
    
     void CalculateAndAddKmBd(                                             
        MatrixType& rKd,
        const MatrixType& rBd) const;

     void CalculateAndAddNonlinearKm(
        Matrix& rLeftHandSideMatrix,
        const Matrix& rB,
        const Matrix& rD,
        const double IntegrationWeight,
        const double IntegrationWeight_zeta) const;

    // Calculation of the PK2 stress
    void CalculatePK2Stress(
        const IndexType IntegrationPointIndex,
        array_1d<double, 3>& rPK2MembraneStressCartesian,
        array_1d<double, 3>& rPK2BendingStressCartesian,
        const ProcessInfo& rCurrentProcessInfo) const;

    // Calculation of the Cauchy stress by transforming the PK2 stress
    void CalculateCauchyStress(
        const IndexType IntegrationPointIndex,
        array_1d<double, 3>& rCauchyMembraneStressesCartesian, 
        array_1d<double, 3>& rCauchyBendingStressesCartesian, 
        const ProcessInfo& rCurrentProcessInfo) const;

    // Calculation of the shear force, shear force = derivative of moment
    void CalculateShearForce(
        const IndexType IntegrationPointIndex,
        array_1d<double, 2>& rq, 
        const ProcessInfo& rCurrentProcessInfo) const;

    void CalculateDerivativeOfCurvatureInitial(
        const IndexType IntegrationPointIndex,
        array_1d<double, 3>& rDCurvature_D1,
        array_1d<double, 3>& rDCurvature_D2,
        const Matrix& rHessian) const;

    void CalculateDerivativeOfCurvatureActual(
        const IndexType IntegrationPointIndex,
        array_1d<double, 3>& rDCurvature_D1,
        array_1d<double, 3>& rDCurvature_D2,
        const Matrix& rHessian,
        const KinematicVariables& rKinematicVariables) const;

    void CalculateDerivativeTransformationMatrices(
        const IndexType IntegrationPointIndex,
        std::vector<Matrix>& rDQ_Dalpha_init,
        std::vector<Matrix>& rDTransCartToCov_Dalpha_init,
        const Matrix& rHessian) const;

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

    // void CalculateHessian(
    //     Matrix& Hessian,
    //     const Matrix& rDDN_DDe) const;

    // void CalculateSecondDerivativesOfBaseVectors(
    //     const Matrix& rDDDN_DDDe,
    //     array_1d<double, 3>& rDDa1_DD11,
    //     array_1d<double, 3>& rDDa1_DD12,
    //     array_1d<double, 3>& rDDa2_DD21,
    //     array_1d<double, 3>& rDDa2_DD22) const;

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

};     // Class Shell6pBbarElement
///@}

}  // namespace Kratos.