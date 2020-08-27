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
class Shell5pElement
    : public Element
{
protected:

    /// Internal variables used for metric transformation
    struct KinematicVariables
    {
        // covariant metric
        array_1d<double, 3> metric;
        array_1d<double, 3> curvature;
        array_1d<double, 2> transShear; //shear 

        //base vector 1
        BoundedVector<double, 3> a1;
        BoundedVector<double, 3> a2;

       //array_1d<double, 3> a1;
      // array_1d<double, 3> a2;
        //director
        BoundedVector<double, 3> t;
        //array_1d<double, 3> t;

        //director derivatives wrt physical space
        array_1d<double, 3> dtd1;
        array_1d<double, 3> dtd2;

        //differential area
        double dA;

        /**
        * The default constructor
        * @param Dimension: The size of working space dimension
        */
        KinematicVariables()
        {
            noalias(metric) = ZeroVector(3);
            noalias(curvature) = ZeroVector(3);
            noalias(transShear) = ZeroVector(2);

            noalias(a1) = ZeroVector(3);
            noalias(a2) = ZeroVector(3);
            noalias(t) = ZeroVector(3);
            noalias(dtd1) = ZeroVector(3);
            noalias(dtd2) = ZeroVector(3);

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
    * Internal variables used in the material and geometric stiffness
    */
    struct VariationVariables
    {
        Matrix P;
        Matrix Q1;
        Matrix Q2;
        Matrix S1;
        Matrix S2;
        Matrix Chi11;
        Matrix Chi12;
        Matrix Chi21;
        Matrix Chi22;
        Matrix WI;
        Matrix WJ;

        /**
        * The default constructor
        */
        VariationVariables()
        {
            P = ZeroMatrix(3, 3);
            Q1 = ZeroMatrix(3, 3);
            Q2 = ZeroMatrix(3, 3);
            S1 = ZeroMatrix(3, 3);
            S2 = ZeroMatrix(3, 3);
            Chi11 = ZeroMatrix(3, 3);
            Chi12 = ZeroMatrix(3, 3);
            Chi21 = ZeroMatrix(3, 3);
            Chi22 = ZeroMatrix(3, 3);
            WI = ZeroMatrix(3, 3);
            WJ = ZeroMatrix(3, 3);
        }
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
    typedef Geometry<Node<3>> GeometryType;

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
    virtual ~Shell5pElement() override
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
        return Kratos::make_intrusive<Shell5pElement>(
            NewId, pGeom, pProperties);
    };

    /// Create with Id, pointer to geometry and pointer to property
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
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
        ProcessInfo& rCurrentProcessInfo) override
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
        ProcessInfo& rCurrentProcessInfo) override
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
        ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType mat_size = number_of_nodes * 5;

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
        buffer << "RMElement #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "RMElement #" << Id();
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
    //std::vector<array_1d<double, 3>> m_A_ab_covariant_vector;
    // Components of the curvature coefficient tensor on the contravariant basis
    std::vector<array_1d<double, 3>> reference_Curvature;
    // Components of the shear coefficient tensor on the contravariant basis
    std::vector<array_1d<double, 2>> reference_TransShear;

    // Shape functions at all integration points
    Matrix m_N;
    // Determinant of the geometrical Jacobian.
    Vector m_dA_vector;

    // Transformed curvilinear derivatives into cartesian derivatives/
    std::vector<Matrix> m_cart_deriv;

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
        KinematicVariables& rKinematicVariables, VariationVariables& rVariationVariables);

    // Computes the cartesian derivatives from curvilinear ones
    Matrix CalculateCartesianDerivatives(
        IndexType IntegrationPointIndex,
        const KinematicVariables& rKinematicVariables);

    Matrix CalculateStrainDisplacementOperator(
        IndexType IntegrationPointIndex,
        const KinematicVariables& rActualKinematic,
        const VariationVariables& rVariations);

    Matrix CalculateGeometricStiffness(
        IndexType IntegrationPointIndex,
        const KinematicVariables& rActKin,
        const VariationVariables& ractVar,
        const ConstitutiveVariables& rThisConstitutiveVariables);

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
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    );


    //void Shell5pElement::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo);
    //void Shell5pElement::constructReferenceDirectorL2FitSystem(MatrixType& rLeftHandSideMatrix, MatrixType& rRightHandSideMatrix);
  //  inline void CalculateAndAddKm(
 //       MatrixType& rLeftHandSideMatrix,
  //      const Matrix& B,
  //      const Matrix& D,
   //     const double IntegrationWeight);

  //  inline void CalculateAndAddNonlinearKm(
   //     Matrix& rLeftHandSideMatrix,
   //     const SecondVariations& rSecondVariationsStrain,
  //      const Vector& rSD,
   //     const double IntegrationWeight);


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
        rSerializer.save("reference_Curvature", reference_Curvature);
        rSerializer.save("reference_TransShear", reference_TransShear);
        rSerializer.save("dA_vector", m_dA_vector);
        rSerializer.save("cart_deriv", m_cart_deriv);
        rSerializer.save("constitutive_law_vector", mConstitutiveLawVector);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
        rSerializer.load("curvature", reference_Curvature);
        rSerializer.load("reference_TransShear", reference_TransShear);
        rSerializer.load("dA_vector", m_dA_vector);
        rSerializer.save("cart_deriv", m_cart_deriv);
        rSerializer.load("constitutive_law_vector", mConstitutiveLawVector);
    }

    ///@}

};     // Class Shell5pElement
///@}

}  // namespace Kratos.

#endif // KRATOS_MESHLESS_SHELL_5P_ELEMENT_H_INCLUDED  defined