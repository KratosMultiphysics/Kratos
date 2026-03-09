#if !defined(KRATOS_SHELL_5p_HIERARCHIC_ELEMENT_H_INCLUDED)
#define  KRATOS_SHELL_5p_HIERARCHIC_ELEMENT_H_INCLUDED


// System includes
#include "includes/variables.h"
#include "includes/element.h"

// External includes

// Project includes

// Application includes
#include "iga_application_variables.h"

namespace Kratos
{
///@name Kratos Classes
///@{
/// 5 parameter shell element.
/* Isogeometric hierarchic Reissner-Mindlin shell element parameterized by 5 parameters (5p). */
class Shell5pHierarchicElement
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of Shell5pHierarchicElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(Shell5pHierarchicElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor using an array of nodes
    Shell5pHierarchicElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {};

    /// Constructor using an array of nodes with properties
    Shell5pHierarchicElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {};

    /// default constructor
    Shell5pHierarchicElement() : Element() {};

    /// Destructor.
    virtual ~Shell5pHierarchicElement() override
    {};

    ///@}
    ///@name Create
    ///@{

    /// Create with GeometryType::Pointer
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive<Shell5pHierarchicElement>(
            NewId, pGeom, pProperties);
    };

    /// Create with NodesArrayType
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        KRATOS_ERROR << "Can not construct Shell5pHierarchicElement with NodesArrayType, because "
            << "it is based on high order background with QuadraturePointGeometries. "
            << "QuadraturePointGeometries would loose their evaluated shape function within this constructor."
            << std::endl;
    };

    ///@}
    ///@name Operations
    ///@{

    /* Called to initialize the element.
     * Must be called before any calculation is done
     */
    void Initialize(
        const ProcessInfo& rCurrentProcessInfo) override;

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
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    );

    /// @brief Sets on rResult the ID's of the element degrees of freedom
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const override;

    /// returns dof list
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const override;

    /// Calls CalculateAll without computing StiffnessMatrix
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        // Calculation flags
        const bool CalculateStiffnessMatrixFlag = false;
        const bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();

        SizeType number_of_control_points = GetGeometry().size();
        SizeType mat_size = number_of_control_points * 5;

        //resizing as needed the RHS
        if (rRightHandSideVector.size() != mat_size)
            rRightHandSideVector.resize(mat_size);
        rRightHandSideVector = ZeroVector(mat_size);

        CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }

    /// Calls CalculateAll without computing ResidualVector
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        // Calculation flags
        const bool CalculateStiffnessMatrixFlag = true;
        const bool CalculateResidualVectorFlag = false;
        VectorType temp = Vector();

        SizeType number_of_control_points = GetGeometry().size();
        SizeType mat_size = number_of_control_points * 5;

        //resizing as needed the LHS
        if (rLeftHandSideMatrix.size1() != mat_size || rLeftHandSideMatrix.size2() != mat_size)
            rLeftHandSideMatrix.resize(mat_size, mat_size);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

        CalculateAll(rLeftHandSideMatrix, temp, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }

    /// Calls CalculateAll
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        //calculation flags
        const bool CalculateStiffnessMatrixFlag = true;
        const bool CalculateResidualVectorFlag = true;

        SizeType number_of_control_points = GetGeometry().size();
        SizeType mat_size = number_of_control_points * 5;

        //resizing as needed the LHS
        if (rLeftHandSideMatrix.size1() != mat_size || rLeftHandSideMatrix.size2() != mat_size)
            rLeftHandSideMatrix.resize(mat_size, mat_size);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

        //resizing as needed the RHS
        if (rRightHandSideVector.size() != mat_size)
            rRightHandSideVector.resize(mat_size);
        rRightHandSideVector = ZeroVector(mat_size);

        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }

    /// @brief Stress recovery
    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(
        const ProcessInfo& rCurrentProcessInfo) const override;

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Hierarchic 5p Shell #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Hierarchic 5p Shell #" << Id();

    }

    ///@}

private:
    ///@name Private Member Variables
    ///@{

    /// The vector containing the constitutive laws
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    // curvilinear coordinate zeta (theta3)
    double mZeta;

    // transformation matrix from contravariant curvilinear base to local Cartesian base both in the initial configuration
    Matrix mInitialTransConToCar;

    /// Internal variables used for metric transformation
    struct MetricVariables
    {
        Vector a_ab; // covariant metric
        Vector a_ab_con; // contravariant metric
        Vector curvature; //
        Matrix J; //Jacobian
        Vector a1; //base vector 1
        Vector a2; //base vector 2
        Vector a3_kirchhoff_love; //base vector 3
        Vector a3_kirchhoff_love_tilde; // unnormalized base vector 3, in Kiendl (2011) a_3_tilde
        double dA; //differential area
        Vector a1_con;  // contravariant base vector 1
        Vector a2_con;  // contravariant base vector 2
        Vector Da1_D1;  // derivative of base vector 1 w.r.t. theta1
        Vector Da1_D2;  // derivative of base vector 1 w.r.t. theta2
        Vector Da2_D2;  // derivative of base vector 2 w.r.t. theta2
        Matrix H; //Hessian (second derivative of cartesian coordinates w.r.t. curvilinear coordinates)

        /* The default constructor
         * @param rWorkingSpaceDimension: The size of working space dimension
         * @param rStrainSize: The size of the StrainVector
         */
        MetricVariables(const unsigned int& rWorkingSpaceDimension = 3, const unsigned int& rStrainSize = 5)
        {
            a_ab = ZeroVector(rWorkingSpaceDimension);
            a_ab_con = ZeroVector(rWorkingSpaceDimension);

            curvature = ZeroVector(rWorkingSpaceDimension);

            J = ZeroMatrix(rWorkingSpaceDimension, 2);

            a1 = ZeroVector(rWorkingSpaceDimension);
            a2 = ZeroVector(rWorkingSpaceDimension);
            a3_kirchhoff_love = ZeroVector(rWorkingSpaceDimension);
            a3_kirchhoff_love_tilde = ZeroVector(rWorkingSpaceDimension);

            dA = 1.0;

            a1_con = ZeroVector(rWorkingSpaceDimension);
            a2_con = ZeroVector(rWorkingSpaceDimension);

            Da1_D1 = ZeroVector(rWorkingSpaceDimension);
            Da1_D2 = ZeroVector(rWorkingSpaceDimension);
            Da2_D2 = ZeroVector(rWorkingSpaceDimension);

            H = ZeroMatrix(rWorkingSpaceDimension, rWorkingSpaceDimension);
        }
    };

    /// Internal variables used in the constitutive equations
    struct ConstitutiveVariables
    {
        Vector E; //strain
        Vector S; //stress
        Matrix D; //constitutive matrix

        /* The default constructor
         * @param StrainSize: The size of the strain vector in Voigt notation
         */
        ConstitutiveVariables(const unsigned int& rStrainSize)
        {
            E = ZeroVector(rStrainSize);
            S = ZeroVector(rStrainSize);
            D = ZeroMatrix(rStrainSize, rStrainSize);
        }
    };

    /// Internal variables used in the constitutive equations
    struct SecondVariations
    {
        Matrix B11;
        Matrix B22;
        Matrix B12;
        Matrix B23;
        Matrix B13;

        /// default constructor
        SecondVariations(const unsigned int& mat_size)
        {
            B11 = ZeroMatrix(mat_size, mat_size);
            B22 = ZeroMatrix(mat_size, mat_size);
            B12 = ZeroMatrix(mat_size, mat_size);
            B23 = ZeroMatrix(mat_size, mat_size);
            B13 = ZeroMatrix(mat_size, mat_size);
        }

        /// operator for addition (+)
        SecondVariations operator+ (const SecondVariations& rSecondVariations)
        {
            KRATOS_TRY

            KRATOS_ERROR_IF(B11.size1() != rSecondVariations.B11.size1()) << "Addition of SecondVariations of different size." << std::endl;
            
            unsigned int mat_size = B11.size1();
            SecondVariations second_variations(mat_size);
            second_variations.B11 = B11 + rSecondVariations.B11;
            second_variations.B22 = B22 + rSecondVariations.B22;
            second_variations.B12 = B12 + rSecondVariations.B12;
            second_variations.B23 = B23 + rSecondVariations.B23;
            second_variations.B13 = B13 + rSecondVariations.B13;

            return second_variations;

            KRATOS_CATCH("")
        }
    };

    MetricVariables mInitialMetric = MetricVariables(3, 5);

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

            if (rNumGPThickness == 3)
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

    // here the number of Gauss-Points over the thickness is specified
    GaussQuadratureThickness mGaussIntegrationThickness = GaussQuadratureThickness(3);

    ///@}
    ///@name Operations
    ///@{

    /* @brief Calculation of the Stiffness Matrix
     * @detail Km = integration_weight * B^T * D *B
     * @param B = B matrix
     * @param D = material stiffness matrix
     */
    void CalculateAndAddKm(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& B,
        const Matrix& D,
        const double& rIntegrationWeight );

    /* @brief The method calculates and adds the non-linear part of the stiffness matrix
     * @param SD = stress
     */
    void CalculateAndAddNonlinearKm(
        Matrix& rLeftHandSideMatrix,
        const SecondVariations& SecondVariationsStrain,
        const Vector& SD,
        const double& rIntegrationWeight);

    /// @brief The method calculates the metric
    void CalculateMetric( MetricVariables& rMetric);
    
    /* @brief The method calculates the transformation matrix from the contravariant to the local Cartesian base both in the initial configuration
     * @param rG1_con = first initial contravariant base vector
     * @param rG2_con = second initial contravariant base vector 
     */
    void CalculateTransformationMatrixInitialTransConToCar(
        const array_1d<double, 3>& rG1_con, 
        const array_1d<double, 3>& rG2_con);

    /* @brief The method calculates the transformation matrix from the local Cartesian to the covariant base both in the initial configuration
     * @param rInitialTransCarToCov = transformation matrix from the local Cartesian to the covariant base both in the initial configuration
     */
    void CalculateTransformationMatrixInitialTransCarToCov(Matrix& rInitialTransCarToCov);

    /* @brief The method calculates the transformation matrix from the contravariant to the local Cartesian base both in the actual configuration
     * @param rg1 = first actual covariant base vector
     * @param rg2 = second actual covariant base vector
     * @param rg3 = third actual covariant base vector
     * @param ra2_con = second actual covariant base vector of the mid-surface
     * @param ra3_kirchhoff_love = third actual covariant Kirchhoff-Love base vector of the mid-surface
     */
    void CalculateTransformationMatrixActualTransCovToCar(
        Matrix& rActualTransCovToCar,
        const Vector& rg1,
        const Vector& rg2,
        const Vector& rg3,
        const Vector& ra2_con,
        const Vector& ra3_kirchhoff_love);

    /* @brief Function determines the values of the shear dofs w_1 and w_2 and calculates the shear difference vector as well as the derivatives of itself and its components  
     * @detail Reissner-Mindlin shell with hierarchic rotations (Oesterle 2018)
     * @param rw = shear difference vector
     * @param rDw_D1 = derivative of the shear difference vector w.r.t. theta1
     * @param rDw_D2 = derivative of the shear difference vector w.r.t. theta2
     * @param rw_alpha = the components w_1 and w_2 of the shear difference vector (shear dofs)
     * @param rDw_alpha_Dbeta = derivative of the components w_1 resp. w_2 w.r.t theta1 and theta2
     */
    void CalculateShearDifferenceVector(
        array_1d<double, 3>& rw,
        array_1d<double, 3>& rDw_D1,
        array_1d<double, 3>& rDw_D2,
        array_1d<double, 2>& rw_alpha,
        Matrix& rDw_alpha_Dbeta,
        const MetricVariables& rActualMetric,
        IndexType IntegrationPointIndex = 0);
    
    /* @brief Calculation of the base vectors of the shell body (in contrast to the mid-surface) for the initial configuration
     * @detail A linearized metric (g_alpha = a_alpha + zeta * Da3_Dalpha) is assumed. Base vector a refer to the mid-span and base vector g refer to the shell body.
     * @param rG1 = first covariant base vector of the shell body in the initial configuration
     * @param rG2 = second covariant base vector of the shell body in the initial configuration
     * @param rG1_con = first contravariant base vector of the shell body in the initial configuration
     * @param rG2_con = second contravariant base vector of the shell body in the initial configuration
     */
    void CalculateInitialBaseVectorsLinearised(
        array_1d<double, 3>& rG1,
        array_1d<double, 3>& rG2,
        array_1d<double, 3>& rG1_con,
        array_1d<double, 3>& rG2_con);

    /* @brief Calculation of the base vectors of the shell body (in contrast to the mid-surface) for the actual configuration
     * @detail A linearized metric (g_alpha = a_alpha + zeta * Da3_Dalpha) is assumed.
     * @param rw = shear difference vector
     * @param rDw_D1 = derivative of the shear difference vector w.r.t. theta1
     * @param rDw_D2 = derivative of the shear difference vector w.r.t. theta2
     * @param rg1 = first covariant base vector of the shell body in the actual configuration
     * @param rg2 = second covariant base vector of the shell body in the actual configuration
     * @param rg3 = third covariant base vector of the shell body in the actual configuration
     */
    void CalculateActualBaseVectorsLinearised(
        const MetricVariables& rActualMetric,
        const Vector& rw,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        array_1d<double, 3>& rg1,
        array_1d<double, 3>& rg2,
        array_1d<double, 3>& rg3);

    /* @brief Calculates deformation gradient F for a Gauss point
     * @param rG1, rG2 = base vectors of the shell body of the reference configuration (G3=A3)
     * @param rg1, rg2, rg3 = base vectors of the shell body of the actual configuration
     * @param rF = deformation gradient
     * @param rdetF = determinant of deformation gradient
     */
    void CalculateDeformationGradient(
        const array_1d<double, 3> rG1,
        const array_1d<double, 3> rG2,
        const array_1d<double, 3> rg1,
        const array_1d<double, 3> rg2,
        const array_1d<double, 3> rg3,
        Matrix& rF,
        double& rdetF);

    /* @brief This functions updates the constitutive variables.
     * @param rActualMetric = actual metric
     * @param rw = shear difference vector
     * @param rDw_D1 = derivative of shear difference vector w.r.t. theta1
     * @param rDw_D2 = derivative of shear difference vector w.r.t. theta2
     * @param rThisConstitutiveVariables = The constitutive variables to be calculated
     * @param rValues = The constitutive law parameters
     * @param ThisStressMeasure = The stress measure considered
     */
    void CalculateConstitutiveVariables(
        const MetricVariables& rActualMetric,
        const Vector& rw,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure);

    /* @brief This method calculates the strain according to the Kirchhoff-Love theory
     * @param ra_ab = covariant metric
     */
    void CalculateStrain(
        array_1d<double, 5>& rStrainVector,
        const Vector& ra_ab,
        const Vector& rCurvature);

    /* @brief This method calculates the additional strain components according to the Reissner-Mindlin theory
     * @param rw = shear difference vector
     * @param rDw_D1, rDw_D2 = derivative of the shear difference vector w.r.t. theta1 respecitvely theta2
     * @param ra1, ra2 = first repecitvely second covariant base vector of the actual configuration
     */
    void CalculateStrainReissnerMindlin(
        array_1d<double, 5>& rStrainVectorReissnerMindlin,
        const Vector& rw,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        const Vector& ra1,
        const Vector& ra2);

    /// @brief This method performs a transformation of a curvilinear strain vector with Voigt size 5 to a Cartesian strain vector with Voigt size 6
    void TransformationCurvilinearStrainSize5ToCartesianStrainSize6(
        const Vector& rCurvilinearStrain,
        Vector& rCartesianStrain);

    /* @brief The method calculates the B matrix of the Kirchhoff-Love contributions
     * @param rB = B matrix
     */
    void CalculateB(
        Matrix& rB,
        const MetricVariables& rActualMetric,
        IndexType IntegrationPointIndex = 0);

    /// @brief The method calculates the second variations according to Kirchhoff-Love
    void CalculateSecondVariations(
        SecondVariations& rSecondVariations,
        const MetricVariables& rActualMetric,
        IndexType IntegrationPointIndex = 0);

    /* @brief The method calculates the additional terms of the B matrix and the second variations according to Reissner-Mindlin
     * @param rw = hierarchic shear difference vector
     * @param rDw_D1 = derivative of the shear difference vector w.r.t. theta1
     * @param rDw_D2 = derivative of the shear difference vector w.r.t. theta2
     * @param rw_alpha = components of the shear difference vector
     * @param rDw_alpha_Dbeta = derivative of the component alpha of the shear difference vector w.r.t. theta(beta) which means theta1 or theta2
     */
    void CalculateVariationsReissnerMindlin(
        Matrix& rB,
        SecondVariations& rSecondVariations,
        const Vector& rw,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        const Vector& rw_alpha,
        const Matrix& rDw_alpha_Dbeta,
        const MetricVariables& rActualMetric,
        const bool& rCalculateStiffnessMatrixFlag,
        IndexType IntegrationPointIndex = 0);

    ///@}
    ///@name Protected Operations
    ///@{

    /// @brief It initializes the material
    virtual void InitializeMaterial();

    /// Computes the Hessian of this problem
    void CalculateHessian(
        Matrix& Hessian,
        const Matrix& rDDN_DDe);

    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
    }

    ///@}

}; // Class Shell5pHierarchicElement

///@}

}  // namespace Kratos.

#endif // KRATOS_SHELL_5p_HIERARCHIC_ELEMENT_H_INCLUDED  defined
