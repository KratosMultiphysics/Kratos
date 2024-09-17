#if !defined(KRATOS_SHELL_5P_ELEMENT_H_INCLUDED )
#define  KRATOS_SHELL_5P_ELEMENT_H_INCLUDED


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
/** Reissner-Mindlin shell using Green-Lagrangian strains and formulated using stress resultants
*/
class KRATOS_API(IGA_APPLICATION) Shell5pElement final
    : public Element
{
protected:

    /// Internal variables used for metric transformation
    struct KinematicVariables
    {
        // covariant metric
        array_1d<double, 3> metricChange;
        array_1d<double, 3> curvature;
        array_1d<double, 2> transShear; //shear

        //base vectors current and reference
        BoundedVector<double, 3> a1;
        BoundedVector<double, 3> a2;
        BoundedVector<double, 3> A1;
        BoundedVector<double, 3> A2;

        BoundedVector<double, 3> dud1;
        BoundedVector<double, 3> dud2;
        //director
        BoundedVector<double, 3> t;
        //array_1d<double, 3> t;

        //director derivatives wrt physical space
        array_1d<double, 3> dtd1;
        array_1d<double, 3> dtd2;

        //differential area
        double dA;

        void setZero()
        {
            metricChange = ZeroVector(3);
            curvature = ZeroVector(3);
            transShear = ZeroVector(2);

            a1 = ZeroVector(3);
            a2 = ZeroVector(3);
            A1 = ZeroVector(3);
            A2 = ZeroVector(3);
            t = ZeroVector(3);
            dtd1 = ZeroVector(3);
            dtd2 = ZeroVector(3);
            dud1 = ZeroVector(3);
            dud2 = ZeroVector(3);

            dA = 0;
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

    using Matrix3d = BoundedMatrix<double, 3, 3>;
    using Matrix2d = BoundedMatrix<double, 2, 2>;
    using Matrix32d = BoundedMatrix<double, 3, 2>;
    using Matrix23d = BoundedMatrix<double, 2, 3>;

    /**
    * Internal variables used in the material and geometric stiffness
    */
    struct VariationVariables
    {
        Matrix3d P;
        Matrix3d Q1;
        Matrix3d Q2;
        Matrix3d S1;
        Matrix3d S2;
        Matrix3d Chi11;
        Matrix3d Chi12Chi21;
        Matrix3d Chi22;
    };

public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of Shell5pElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(Shell5pElement);

    /// Size types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    // GometryType
    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor using an array of nodes
    Shell5pElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {};

    /// Constructor using an array of nodes with properties
    Shell5pElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {};

    /// Default constructor necessary for serialization
    Shell5pElement()
        : Element()
    {};

    /// Destructor.
    virtual ~Shell5pElement() final
    {};

    ///@}
    ///@name Life Cycle
    ///@{

    /// Create with Id, pointer to geometry and pointer to property
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const final
    {
        return Kratos::make_intrusive<Shell5pElement>(
            NewId, pGeom, pProperties);
    };

    /// Create with Id, pointer to geometry and pointer to property
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const final
    {
        return Kratos::make_intrusive< Shell5pElement >(
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
        const ProcessInfo& rCurrentProcessInfo) final
    {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType mat_size = number_of_nodes * 5;

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
        const ProcessInfo& rCurrentProcessInfo) final
    {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType mat_size = number_of_nodes * 5;

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
        const ProcessInfo& rCurrentProcessInfo) final
    {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType mat_size = number_of_nodes * 5;

        if (rRightHandSideVector.size() != mat_size)
            rRightHandSideVector.resize(mat_size);
        noalias(rRightHandSideVector) = ZeroVector(mat_size);

        if (rLeftHandSideMatrix.size1() != mat_size || rLeftHandSideMatrix.size2() != mat_size)
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
        const ProcessInfo& rCurrentProcessInfo
    ) const final;

    /**
    * @brief Sets on rConditionDofList the degrees of freedom of the considered element geometry
    * @param rElementalDofList The vector containing the dof of the element
    * @param rCurrentProcessInfo The current process info instance
    */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const final;

    ///@}
    ///@name Base Class Operations
    ///@{

    void Initialize(const ProcessInfo& rCurrentProcessInfo) final;

    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) final;

    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) final;

    void GetValuesVector(
        Vector& rValues,
        int Step) const final;

    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step) const final;

    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step) const final;

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
    int Check(const ProcessInfo& rCurrentProcessInfo) const final {
        auto& r_geometry = GetGeometry();
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            KRATOS_ERROR_IF_NOT(r_geometry[i].Has(DIRECTOR))
                << "DIRECTOR not provided at node #" << r_geometry[i].Id() << std::endl;
        }

        return 0;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const final
    {
        std::stringstream buffer;
        buffer << "RMElement #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const final
    {
        rOStream << "RMElement #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    // Components of the metric coefficient tensor on the contravariant basis
    //std::vector<array_1d<double, 3>> m_A_ab_covariant_vector;
    // Components of the curvature coefficient tensor on the contravariant basis
    std::vector<array_1d<double, 3>> reference_Curvature;
    // Components of the shear coefficient tensor on the contravariant basis
    std::vector<array_1d<double, 2>> reference_TransShear;

    // Determinant of the geometrical Jacobian.
    Vector m_dA_vector;

    Kratos::Variable<Vector>::Type const& (NodeType::*m_GetValueFunctor)(const Kratos::Variable<Vector >&) const = &NodeType::GetValue;
    Kratos::Variable<array_1d<double, 3>>::Type const& (NodeType::*m_FastGetSolutionStepValueFunctor)(const Kratos::Variable<array_1d<double,3> >&) const = &NodeType::FastGetSolutionStepValue;
    NodeType::PointType::CoordinatesArrayType const& (NodeType::* m_GetCoordinatesFunctor)() const = &NodeType::Coordinates;
    NodeType::PointType const& (NodeType::* m_GetInitialPositionFunctor)() const = &NodeType::GetInitialPosition;

    // Transformed curvilinear derivatives into cartesian derivatives/
    std::vector<Matrix> m_cart_deriv;

    //The St. Venant Kirchhoff 5P Material Tangent
    BoundedMatrix<double, 8, 8> mC;

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

    std::pair< Shell5pElement::KinematicVariables, Shell5pElement::VariationVariables>
        CalculateKinematics(const IndexType IntegrationPointIndex) const;

    // Computes the cartesian derivatives from curvilinear ones
    Matrix CalculateCartesianDerivatives(
        const IndexType IntegrationPointIndex);

    Matrix CalculateStrainDisplacementOperator(
        const IndexType IntegrationPointIndex,
        const KinematicVariables& rActualKinematic,
        const VariationVariables& rVariations) const ;

    Matrix CalculateGeometricStiffness(
        const IndexType IntegrationPointIndex,
        const KinematicVariables& rActKin,
        const VariationVariables& ractVar,
        const ConstitutiveVariables& rThisConstitutiveVariables) const;

    /**
    * This functions updates the constitutive variables
    * @param rActualMetric: The actual metric
    * @param rThisConstitutiveVariables: The constitutive variables to be calculated
    * @param rValues: The CL parameters
    * @param ThisStressMeasure: The stress measure considered
    */
    void CalculateConstitutiveVariables(
        const IndexType IntegrationPointIndex,
        const KinematicVariables& rActualMetric,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    );

    /// Helper
    void CalculateSVKMaterialTangent();

    template< typename ContainerType, typename NodeFunctor, typename ...Args>
    BoundedVector<double, 3> InterpolateNodalVariable(const ContainerType& vec, const NodeFunctor& funct, const Args&... args) const;

    public:
    static BoundedMatrix<double, 3, 2> TangentSpaceFromStereographicProjection(const array_1d<double, 3 >& director);
    private:
    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const final
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
        rSerializer.save("reference_Curvature", reference_Curvature);
        rSerializer.save("reference_TransShear", reference_TransShear);
        rSerializer.save("dA_vector", m_dA_vector);
        rSerializer.save("cart_deriv", m_cart_deriv);
    }

    void load(Serializer& rSerializer) final
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
        rSerializer.load("curvature", reference_Curvature);
        rSerializer.load("reference_TransShear", reference_TransShear);
        rSerializer.load("dA_vector", m_dA_vector);
        rSerializer.save("cart_deriv", m_cart_deriv);
    }

    ///@}

};     // Class Shell5pElement
///@}

}  // namespace Kratos.

#endif // KRATOS_MESHLESS_SHELL_5P_ELEMENT_H_INCLUDED  defined