
#if !defined(KRATOS_IGA_SHELL_3P_ELEMENT_H_INCLUDED)
#define  KRATOS_IGA_SHELL_3P_ELEMENT_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_elements/base_discrete_element.h"

// Application includes
#include "iga_application_variables.h"

namespace Kratos
{
///@name Kratos Classes
///@{
/// Short class definition.
/** 
 * @class IgaShell3pElement
 * @ingroup IGAApplication
 * @brief Reissner-Mindlin 5-Parameter hierarchic shell element by Echter et al. (2013)
    */
class IgaShell3pElement
    : public BaseDiscreteElement
{
public:
    ///@name Type Definitions
    ///@{
    /// Counted pointer of IgaShell3pElement
    KRATOS_CLASS_POINTER_DEFINITION(IgaShell3pElement);
    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
    // Constructor using an array of nodes
    IgaShell3pElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseDiscreteElement(NewId, pGeometry)
    {};
    // Constructor using an array of nodes with properties
    IgaShell3pElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseDiscreteElement(NewId, pGeometry, pProperties)
    {};

    // default constructor necessary for serialization
    IgaShell3pElement() : BaseDiscreteElement() 
    {};

    /// Destructor.
    virtual ~IgaShell3pElement() override
    {};

    /**
    * @brief Creates a new element
    * @param NewId The Id of the new created element
    * @param pGeom The pointer to the geometry of the element
    * @param pProperties The pointer to property
    * @return The pointer to the created element
    */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override
    {
        
        KRATOS_WATCH("CreateStart");
        
        return Kratos::make_shared<IgaShell3pElement>(
            NewId, pGeom, pProperties);
    };

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {

        KRATOS_WATCH("CreateStart");
        
        return Kratos::make_shared<IgaShell3pElement>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    };

    ///@}
    ///@name Operations
    ///@{

    /**
    * Called to initialize the element.
    * Must be called before any calculation is done
    */
    void Initialize() override;

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
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    ) override;
    
    // function not checked yet, just copied from surface_base, needed? (ML)
    /**
    * Calculation of the Material Stiffness Matrix. Km = B^T * D *B
    */
    void CalculateAndAddKm(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& B,
        const Matrix& D,
        const double IntegrationWeight) {KRATOS_WATCH("CalculateAndAddKm");};

    // function not checked yet, just copied from surface_base, needed? (ML)
    void CalculateAndAddNonlinearKm(
        Matrix& rLeftHandSideMatrix,
        // const SecondVariations& SecondVariationsStrain,
        const Vector& SD,
        const double& rIntegrationWeight) {KRATOS_WATCH("CalculateAndAddNonlinearKm");};
   
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
    * @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
    * @param rElementalDofList The vector containing the dof of the element
    * @param rCurrentProcessInfo The current process info instance
    */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
    * This function provides the place to perform checks on the completeness of the input.
    * It is designed to be called only once (or anyway, not often) typically at the beginning
    * of the calculations, so to verify that nothing is missing from the input
    * or that no common error is found.
    * @param rCurrentProcessInfo
    */
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    std::string Info() const override
    {
        
        KRATOS_WATCH("Check");

        std::stringstream buffer;
        buffer << "IgaShell3pElement #" << Id();
        KRATOS_WATCH(GetValue(SHAPE_FUNCTION_VALUES));
        KRATOS_WATCH(GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES));
        KRATOS_WATCH(GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES));
        KRATOS_WATCH(GetValue(INTEGRATION_WEIGHT));
        return buffer.str();
    };

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        
        KRATOS_WATCH("PrintInfo");
        
        rOStream << "IgaShell3pElement #" << Id();
    };

    // function not checked yet, just copied from surface_base, needed? (ML)
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo) override {KRATOS_WATCH("CalculateMassMatrix");};

    ///@}

protected:

private:
    ///@name Static Member Variables
    ///@{

    ///@}
	///@name Member Variables
	///@{
    /**
        * Internal variables used for metric transformation
        */
    struct MetricVariables
    {
        Vector gab; // covariant metric (only in-plane, g_11, g_22, g_12)       // not used for sure (ML)
        Vector gab_con; // contravariant metric (only in-plane)     // not used for sure (ML)
        Vector curvature; //
        Matrix J; //Jacobian
        double  detJ;     // not used for sure (ML)
        Vector g1, g2, g3; //base vector 1, 2, 3
        // derivatives of base vector w.r.t. to convective coordinates theta1 and theta2
        Vector Dg1_D1, Dg1_D2, Dg2_D1, Dg2_D2, Dg3_D1, Dg3_D2;
        double dA; //differential area
        Matrix H; //Hessian
        Matrix Q; //Transformation matrix Q from contravariant to cartesian basis
        Matrix T; //Transformation matrix T from contravariant to local cartesian basis
        
        /**
            * The default constructor
            * @param Dimension: The size of working space dimension
            */
        MetricVariables(const unsigned int& Dimension)
        {
            gab = ZeroVector(Dimension);
            gab_con = ZeroVector(Dimension);

            curvature = ZeroVector(Dimension);

            J = ZeroMatrix(Dimension, Dimension);
            detJ = 1.0;

            g1 = ZeroVector(Dimension);
            g2 = ZeroVector(Dimension);
            g3 = ZeroVector(Dimension);

            Dg1_D1 = ZeroVector(Dimension);
            Dg1_D2 = ZeroVector(Dimension);
            Dg2_D1 = ZeroVector(Dimension);
            Dg2_D2 = ZeroVector(Dimension);
            Dg3_D1 = ZeroVector(Dimension);
            Dg3_D2 = ZeroVector(Dimension);

            dA = 1.0;

            Matrix H = ZeroMatrix(Dimension, Dimension);
            Matrix Q = ZeroMatrix(Dimension, Dimension);
            Matrix T = ZeroMatrix(Dimension, Dimension);
        }
    };

    MetricVariables m_initial_metric = MetricVariables(3);
    
    // rotation vector and its derivatives
    Vector m_Phi, m_DPhi_D1, m_DPhi_D2;
    // rotation angle
    double m_phi1, m_phi2;

    // first derivative of the mid-surface displacement v w.r.t. theta1 (convective coordinate) 
    Vector m_Dv_D1;
    // first derivative of the mid-surface displacement v w.r.t. theta2 (convective coordinate) 
    Vector m_Dv_D2;        

    /**
    * Internal variables used in the constitutive equations
    */
    struct ConstitutiveVariables        // should be checked before usage (ML)
    {
        Vector E; //strain
        Vector S; //stress
        Matrix D; //constitutive matrix

        /**
        * @brief The default constructor
        * @param StrainSize: The size of the strain vector in Voigt notation
        */
        ConstitutiveVariables(const unsigned int& StrainSize)
        {
            E = ZeroVector(StrainSize);
            S = ZeroVector(StrainSize);
            D = ZeroMatrix(StrainSize, StrainSize);
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


    ///@{        
    ///@name Operations
    ///@{
   
    /**
        * @brief This function calculates all metric variables
        */
    void CalculateMetric(MetricVariables& rMetric);

    /**
     * @brief This function calculates the rotation vector and its derivatives
     * @details The rotation vector is needed to compute the difference
     * between the director in the undeformed and deformed configuration
     * @param rActualMetric: The actual metric
     */
    void CalculateRotationVector(
        MetricVariables& rAcutalMetric);

    /**
        * @brief This functions updates the constitutive variables
        * @param rActualMetric: The actual metric
        * @param rThisConstitutiveVariables: The constitutive variables to be calculated
        * @param rValues: The CL parameters
        * @param ThisStressMeasure: The stress measure considered
        */
    void CalculateConstitutiveVariables(
        MetricVariables& rActualMetric,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure);
    
    /** 
     * @brief This function computes the membrane strains (strains at mid-surface resp. constant part of strains)
     * @param rStrainVector: container to save the calculated membrane strain
     */
    void CalculateStrain(
        Vector& rStrainVector);
    
    /** 
     * @brief This function computes the strains due to bending (strains depending on theta3 (thickness direction))
     * @param rCurvatureVector: container to save the calculated strain due to bending
     */
    void CalculateCurvature(
        Vector& rCurvatureVector,
        const MetricVariables& rActualMetric);
    
    // function not checked yet, just copied from surface_base, needed? (ML)
    void CalculateBMembrane(        
        Matrix& rB,
        const MetricVariables& rMetric,
        const Matrix& rDN_De,
        const unsigned int& rNumberOfNodes,
        const unsigned int& rMatSize);

    // function not checked yet, just copied from surface_base, needed? (ML)
    void CalculateBCurvature(
        Matrix& rB,
        const MetricVariables& metric) {KRATOS_WATCH("CalculateBCurvature");};

    ///@}

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
    };

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
    };

    ///@}

};     // Class IgaShell3pElement
///@}

}  // namespace Kratos.

#endif // KRATOS_IGA_SHELL_5P_ELEMENT_H_INCLUDED  defined