#if !defined(KRATOS_BEAM_THICK_ELEMENT_2D_H_INCLUDED )
#define  KRATOS_BEAM_THICK_ELEMENT_2D_H_INCLUDED


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
/** Non-linear 2D Timoshenko Beam. Optimized for Isogeometric Analysis.
*/
class KRATOS_API(IGA_APPLICATION) BeamThickElement2D
    : public Element
{
protected:

    /// Internal variables used for metric transformation
    struct KinematicVariables
    {
        // covariant metric
        double a_11_covariant;
        double b_11_covariant;

        //base vector 1
        array_1d<double, 3> a1;
        //base vector 2 normalized
        array_1d<double, 3> a2;
        //not-normalized base vector 2
        array_1d<double, 3> a2_tilde;

        //beta
        double beta;
        //first derivative of beta
        double beta_deriv;


        //differential length
        double dL;

        /**
        * The default constructor
        * @param Dimension: The size of working space dimension
        */
        KinematicVariables(SizeType Dimension)
        {
            a_11_covariant = 0.0;
            b_11_covariant = 0.0;

            beta = 0.0;
            beta_deriv = 0.0;

            noalias(a1) = ZeroVector(Dimension);
            noalias(a2) = ZeroVector(Dimension);
            noalias(a2_tilde) = ZeroVector(Dimension);

            dL = 1.0;
        }
    };

    /**
    * Internal variables used in the constitutive equations
    */
    struct ConstitutiveVariables
    {
        double StrainValue;
        double StressValue;
        double ConstitutiveValue;

        ConstitutiveVariables()
        {
            StrainValue = 0.0;
            StressValue = 0.0;
            ConstitutiveValue = 0.0;
        }
    };

public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of BeamThickElement2D
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(BeamThickElement2D);

    /// Size types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    // GometryType
    typedef Geometry<Node> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor using an array of nodes
    BeamThickElement2D(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {};

    /// Constructor using an array of nodes with properties
    BeamThickElement2D(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {};

    /// Default constructor necessary for serialization
    BeamThickElement2D()
        : Element()
    {};

    /// Destructor.
    virtual ~BeamThickElement2D() = default;

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
        return Kratos::make_intrusive<BeamThickElement2D>(
            NewId, pGeom, pProperties);
    };

    /// Create with Id, pointer to geometry and pointer to property
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive< BeamThickElement2D >(
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
        const ProcessInfo& rCurrentProcessInfo) override
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
        const ProcessInfo& rCurrentProcessInfo) override
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
        buffer << "Kirchhoff-Love BeamThickElement2D #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Kirchhoff-Love BeamThickElement2D #" << Id();
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
    std::vector<double> m_A_11_covariant_vector;
    // Components of the curvature coefficient tensor on the contravariant basis
    std::vector<double> m_B_11_covariant_vector;
    // Determinant of the geometrical Jacobian.
    Vector m_dL_vector;

    /// The vector containing the constitutive laws for all integration points.
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

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
    ) const;

    /// Initialize Operations
    void InitializeMaterial();

    void CalculateKinematics(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables) const;

    void CalculateBMembrane(
        const IndexType IntegrationPointIndex,
        Vector& rB,
        const KinematicVariables& rActualKinematic) const;

    void CalculateBCurvature(
        const IndexType IntegrationPointIndex,
        Vector& rB,
        const KinematicVariables& rActualKinematic) const;
    
    void CalculateBShear(
        const IndexType IntegrationPointIndex,
        Vector& rB,
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
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature,
        ConstitutiveVariables& rThisConstitutiveVariablesShear,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    ) const;

    inline void CalculateAndAddKm(
        MatrixType& rLeftHandSideMatrix,
        const Vector& rB,
        const double& rD,
        const double IntegrationWeight) const;

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

    void CalculateHessian(
        array_1d<double, 3>& Hessian,
        const Matrix& rDDN_DDe) const;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
        rSerializer.save("A_11_covariant", m_A_11_covariant_vector);
        rSerializer.save("B_11_covariant", m_B_11_covariant_vector);
        rSerializer.save("dA_vector", m_dL_vector);
        rSerializer.save("constitutive_law_vector", mConstitutiveLawVector);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
        rSerializer.load("A_11_covariant", m_A_11_covariant_vector);
        rSerializer.load("B_11_covariant", m_B_11_covariant_vector);
        rSerializer.load("dA_vector", m_dL_vector);
        rSerializer.load("constitutive_law_vector", mConstitutiveLawVector);
    }

    ///@}

};     // Class BeamThickElement2D
///@}

}  // namespace Kratos.

#endif // KRATOS__BEAM_THICK_ELEMENT_2D_H_INCLUDED  defined
