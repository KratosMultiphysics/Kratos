#if !defined(KRATOS_SHELL_3P_ELEMENT_H_INCLUDED )
#define  KRATOS_SHELL_3P_ELEMENT_H_INCLUDED


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
/** Kirchhoff-Love Shell. Optimized for Isogeometric Analysis by Kiendl et al. .
*/
class Shell3pElement
    : public Element
{
protected:

    /// Internal variables used for metric transformation
    struct KinematicVariables
    {
        // covariant metric
        array_1d<double, 3> a_ab_covariant;
        array_1d<double, 3> b_ab_covariant;

        //base vector 1
        array_1d<double, 3> a1;
        //base vector 2
        array_1d<double, 3> a2;
        //base vector 3 normalized
        array_1d<double, 3> a3;
        //not-normalized base vector 3
        array_1d<double, 3> a3_tilde;

        //differential area
        double dA;

        /**
        * The default constructor
        * @param Dimension: The size of working space dimension
        */
        KinematicVariables(SizeType Dimension)
        {
            noalias(a_ab_covariant) = ZeroVector(Dimension);
            noalias(b_ab_covariant) = ZeroVector(Dimension);

            noalias(a1) = ZeroVector(Dimension);
            noalias(a2) = ZeroVector(Dimension);
            noalias(a3) = ZeroVector(Dimension);

            noalias(a3_tilde) = ZeroVector(Dimension);

            dA = 1.0;
        }
    };

    /**
    * Internal variables used in the constitutive equations
    */
    struct ConstitutiveVariables
    {
        Vector StrainVector;
        Vector StressVector;
        Matrix ConstitutiveMatrix;

        /**
        * @param StrainSize: The size of the strain vector in Voigt notation
        */
        ConstitutiveVariables(SizeType StrainSize)
        {
            StrainVector = ZeroVector(StrainSize);
            StressVector = ZeroVector(StrainSize);
            ConstitutiveMatrix = ZeroMatrix(StrainSize, StrainSize);
        }
    };

    /**
    * Internal variables used in the constitutive equations
    */
    struct SecondVariations
    {
        Matrix B11;
        Matrix B22;
        Matrix B12;

        /**
        * The default constructor
        * @param StrainSize: The size of the strain vector in Voigt notation
        */
        SecondVariations(const int& mat_size)
        {
            B11 = ZeroMatrix(mat_size, mat_size);
            B22 = ZeroMatrix(mat_size, mat_size);
            B12 = ZeroMatrix(mat_size, mat_size);
        }
    };

public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of Shell3pElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(Shell3pElement);

    /// Size types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    // GometryType
    typedef Geometry<Node<3>> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor using an array of nodes
    Shell3pElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {};

    /// Constructor using an array of nodes with properties
    Shell3pElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {};

    /// Default constructor necessary for serialization
    Shell3pElement()
        : Element()
    {};

    /// Destructor.
    virtual ~Shell3pElement() override
    {};

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
        return Kratos::make_intrusive<Shell3pElement>(
            NewId, pGeom, pProperties);
    };

    /// Create with Id, pointer to geometry and pointer to property
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive< Shell3pElement >(
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
        ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType mat_size = number_of_nodes * 3;

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
        ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType mat_size = number_of_nodes * 3;

        VectorType right_hand_side_vector;

        if (rLeftHandSideMatrix.size1() != mat_size)
            rLeftHandSideMatrix.resize(mat_size, mat_size);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

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
        ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType mat_size = number_of_nodes * 3;

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
    * @brief Sets on rResult the ID's of the element degrees of freedom
    * @param rResult The vector containing the equation id
    * @param rCurrentProcessInfo The current process info instance
    */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
    * @brief Sets on rConditionDofList the degrees of freedom of the considered element geometry
    * @param rElementalDofList The vector containing the dof of the element
    * @param rCurrentProcessInfo The current process info instance
    */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
    ) override;

    ///@}
    ///@name Base Class Operations
    ///@{

    void Initialize() override;

    void GetValuesVector(
        Vector& rValues,
        int Step) override;

    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step) override;

    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step) override;

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
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "KLElement #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "KLElement #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {
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

    ///@}
    ///@name Operations
    ///@{

    /// Calculates LHS and RHS dependent on flags
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    );

    /// Initialize Operations
    void InitializeMaterial();

    void CalculateKinematics(
        IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables);

    // Computes transformation
    void CalculateTransformation(
        const KinematicVariables& rKinematicVariables,
        Matrix& rT);

    void CalculateBMembrane(
        IndexType IntegrationPointIndex,
        Matrix& rB,
        const KinematicVariables& rActualKinematic);

    void CalculateBCurvature(
        IndexType IntegrationPointIndex,
        Matrix& rB,
        const KinematicVariables& rActualKinematic);

    void CalculateSecondVariationStrainCurvature(
        IndexType IntegrationPointIndex,
        SecondVariations& rSecondVariationsStrain,
        SecondVariations& rSecondVariationsCurvature,
        const KinematicVariables& rActualKinematic);

    /**
    * This functions updates the constitutive variables
    * @param rActualMetric: The actual metric
    * @param rThisConstitutiveVariables: The constitutive variables to be calculated
    * @param rValues: The CL parameters
    * @param ThisStressMeasure: The stress measure considered
    */
    void CalculateConstitutiveVariables(
        IndexType IntegrationPointIndex,
        KinematicVariables& rActualMetric,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    );

    inline void CalculateAndAddKm(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& B,
        const Matrix& D,
        const double IntegrationWeight);

    inline void CalculateAndAddNonlinearKm(
        Matrix& rLeftHandSideMatrix,
        const SecondVariations& rSecondVariationsStrain,
        const Vector& rSD,
        const double IntegrationWeight);

    ///@}
    ///@name Geometrical Functions
    ///@{

    void CalculateHessian(
        Matrix& Hessian,
        const Matrix& rDDN_DDe);

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void load(Serializer& rSerializer) const override
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

};     // Class Shell3pElement
///@}

}  // namespace Kratos.

#endif // KRATOS_MESHLESS_SHELL_3P_ELEMENT_H_INCLUDED  defined