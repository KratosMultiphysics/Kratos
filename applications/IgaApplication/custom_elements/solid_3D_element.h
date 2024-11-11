#if !defined(KRATOS_SOLID_3D_ELEMENT_H_INCLUDED )
#define  KRATOS_SOLID_3D_ELEMENT_H_INCLUDED


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
///@name Kratos Classes
///@{
/// Short class definition.
/** Kirchhoff-Love Shell. Optimized for Isogeometric Analysis by Kiendl et al. .
*/
class Solid3DElement
    : public Element
{
protected:
/**
     * Internal variables used in the kinematic calculations
     */
    struct KinematicVariables
    {
        Vector  N;
        Matrix  B;
        double  detF;
        Matrix  F;
        double  detJ0;
        Matrix  J0;
        Matrix  InvJ0;
        Matrix  DN_DX;
        Vector Displacements;

        /**
         * The default constructor
         * @param StrainSize The size of the strain vector in Voigt notation
         * @param Dimension The problem dimension: 2D or 3D
         * @param NumberOfNodes The number of nodes in the element
         */
        KinematicVariables(
            const SizeType StrainSize,
            const SizeType Dimension,
            const SizeType NumberOfNodes
            )
        {
            detF = 1.0;
            detJ0 = 1.0;
            N = ZeroVector(NumberOfNodes);
            B = ZeroMatrix(StrainSize, Dimension * NumberOfNodes);
            F = IdentityMatrix(Dimension);
            DN_DX = ZeroMatrix(NumberOfNodes, Dimension);
            J0 = ZeroMatrix(Dimension, Dimension);
            InvJ0 = ZeroMatrix(Dimension, Dimension);
            Displacements = ZeroVector(Dimension * NumberOfNodes);
        }
    };

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

    /**
     * @brief This method returns if the element provides the strain
     */
    bool UseElementProvidedStrain() const;
public:

    ///@name Type Definitions
    ///@{

    /// Counted pointer of Solid3DElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(Solid3DElement);

    // static constexpr std::size_t NumNodes = TDim + 1;

    ///@}
    ///@name Life Cycle
    ///@{

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /// Default constructor.
    Solid3DElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry);

    Solid3DElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~Solid3DElement();

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


    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    Solid3DElement() : Element()
    {
    }

    // New functions NICO
    const virtual GeometryType::IntegrationPointsArrayType  IntegrationPoints() const 
    {
        return GetGeometry().IntegrationPoints();
    }

    const virtual GeometryType::IntegrationPointsArrayType  IntegrationPoints(IntegrationMethod ThisMethod) const
    {
        return GetGeometry().IntegrationPoints(ThisMethod);
    }

    const Matrix& ShapeFunctionsValues(IntegrationMethod ThisMethod) const
    {
        return GetGeometry().ShapeFunctionsValues(ThisMethod);
    }

    // void CalculateRightHandSide(
    //     VectorType& rRightHandSideVector,
    //     const ProcessInfo& rCurrentProcessInfo) override
    // {
    //     const SizeType number_of_nodes = GetGeometry().size();
    //     const SizeType mat_size = number_of_nodes * 3;

    //     if (rRightHandSideVector.size() != mat_size)
    //         rRightHandSideVector.resize(mat_size);
    //     noalias(rRightHandSideVector) = ZeroVector(mat_size);

    //     MatrixType left_hand_side_matrix;

    //     KRATOS_WATCH('RHS1')
        
    //     // CalculateAll(left_hand_side_matrix, rRightHandSideVector,
    //     //     rCurrentProcessInfo, false, true);
    //     MatrixType temp(0,0);
    //     CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
    // }


    // void CalculateLeftHandSide(
    //     MatrixType& rLeftHandSideMatrix,
    //     const ProcessInfo& rCurrentProcessInfo) override
    // {
    //     KRATOS_WATCH('LHS0')
    //     const SizeType number_of_nodes = GetGeometry().size();
    //     const SizeType mat_size = number_of_nodes * 3;

    //     VectorType right_hand_side_vector;

    //     if (rLeftHandSideMatrix.size1() != mat_size)
    //         rLeftHandSideMatrix.resize(mat_size, mat_size);
    //     noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    //     // CalculateAll(rLeftHandSideMatrix, right_hand_side_vector,
    //     //     rCurrentProcessInfo, true, false);
        
    //     KRATOS_WATCH('LHS1')

    //     VectorType temp(0);
    //     CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
    // }


    // void CalculateLocalSystem(
    //     MatrixType& rLeftHandSideMatrix,
    //     VectorType& rRightHandSideVector,
    //     const ProcessInfo& rCurrentProcessInfo) override
    // {
    //     const SizeType number_of_nodes = GetGeometry().size();
    //     const SizeType mat_size = number_of_nodes * 3;

    //     if (rRightHandSideVector.size() != mat_size)
    //         rRightHandSideVector.resize(mat_size);
    //     noalias(rRightHandSideVector) = ZeroVector(mat_size);

    //     if (rLeftHandSideMatrix.size1() != mat_size)
    //         rLeftHandSideMatrix.resize(mat_size, mat_size);
    //     noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    //     CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
    //         rCurrentProcessInfo, true, true);
    // }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{
    
    IntegrationMethod GetIntegrationMethod() const override;

    /**
     * @brief Returns the used integration method
     * @return default integration method of the used Geometry
     */
    // IntegrationMethod GetIntegrationMethod() const override
    // {
    //     return mThisIntegrationMethod;
    // }

    // std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector; /// The vector containing the constitutive laws


    // const GeometryType::IntegrationPointsArrayType  IntegrationPoints() const 
    // {
    //     return GetGeometry().IntegrationPoints();
    // }
    ///@}
    ///@name Input and output
    ///@{


    ///@}
    ///@name Friends
    ///@{


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
    
    // void Solid3DElement::CalculateOnIntegrationPoints(
    //     const Variable<Vector>& rVariable,
    //     std::vector<Vector>& rOutput,
    //     const ProcessInfo& rCurrentProcessInfo
    // ) override;
    /**
     * @brief This function computes the body force
     * @param IntegrationPoints The array containing the integration points
     * @param PointNumber The id of the integration point considered
     * @return The vector of body forces
     */
    // virtual array_1d<double, 3> GetBodyForce(
    //     const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    //     const IndexType PointNumber
    //     ) const;

    /**
     * @brief This functions updates the kinematics variables
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param PointNumber The integration point considered
     */
    void CalculateKinematicVariables(
        KinematicVariables& rThisKinematicVariables,
        const IndexType PointNumber,
        const GeometryType::IntegrationMethod& rIntegrationMethod
        ) ;


    void CalculateB(
        Matrix& rB,
        Matrix& r_DN_DX) const;

    void GetValuesVector(
        Vector& rValues) const;        

protected:
    ///@name Protected static Member Variables
    ///@{

    void InitializeMaterial();

    /**
     * @brief Gives the StressMeasure used
     */
    ConstitutiveLaw::StressMeasure GetStressMeasure() const;

    ///@}
    ///@name Protected member Variables
    ///@{

    ConstitutiveLaw::Pointer mpConstitutiveLaw; /// The pointer containing the constitutive laws

    IntegrationMethod mThisIntegrationMethod; /// Currently selected integration methods

    /**
     * @brief This functions updates the constitutive variables
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param rThisConstitutiveVariables The constitutive variables
     * @param rValues The CL parameters
     * @param PointNumber The integration point considered
     * @param IntegrationPoints The list of integration points
     * @param ThisStressMeasure The stress measure considered
     */
    void CalculateConstitutiveVariables(
        KinematicVariables& rThisKinematicVariables,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure = ConstitutiveLaw::StressMeasure_PK2
        );


        /**
     * @brief This functions updates the data structure passed to the CL
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param rThisConstitutiveVariables The constitutive variables
     * @param rValues The CL parameters
     * @param PointNumber The integration point considered
     * @param IntegrationPoints The list of integration points
     */
    void SetConstitutiveVariables(
        KinematicVariables& rThisKinematicVariables,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints
        );


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

    
    // /// Calculates LHS and RHS dependent on flags
    // void CalculateAll(
    //     MatrixType& rLeftHandSideMatrix,
    //     VectorType& rRightHandSideVector,
    //     const ProcessInfo& rCurrentProcessInfo,
    //     const bool CalculateStiffnessMatrixFlag,
    //     const bool CalculateResidualVectorFlag
    // ) const;


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

    /// Assignment operator.
    //Solid3DElement& operator=(const Solid3DElement& rOther);

    /// Copy constructor.
    //Solid3DElement(const Solid3DElement& rOther);

    ///@}

}; // Class Solid3DElement

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    Solid3DElement& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const Solid3DElement& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);
      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_SOLID_3D_ELEMENT_H_INCLUDED  defined