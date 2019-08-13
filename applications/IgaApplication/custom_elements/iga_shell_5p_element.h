#if !defined(KRATOS_IGA_SHELL_5p_ELEMENT_H_INCLUDED)
#define  KRATOS_IGA_SHELL_5p_ELEMENT_H_INCLUDED


// System includes
#include "includes/variables.h"

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
/** 3D Shell with hierarchical shear vector (7p). Optimized for Isogeometric Analysis by Echter et al..
*/
class IgaShell5pElement
    : public BaseDiscreteElement
{
public:
    ///@name Type Definitions
    ///@{
    /// Counted pointer of IgaShell5pElement
    KRATOS_CLASS_POINTER_DEFINITION(IgaShell5pElement);
    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
    // Constructor using an array of nodes
    IgaShell5pElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseDiscreteElement(NewId, pGeometry)
    {};
    // Constructor using an array of nodes with properties
    IgaShell5pElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseDiscreteElement(NewId, pGeometry, pProperties)
    {};

    // default constructor necessary for serialization
    IgaShell5pElement() : BaseDiscreteElement() {};

    /// Destructor.
    virtual ~IgaShell5pElement() override
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
        return Kratos::make_shared<IgaShell5pElement>(
            NewId, pGeom, pProperties);
    };

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_shared< IgaShell5pElement >(
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
        std::stringstream buffer;
        buffer << "KLElement #" << Id();
        KRATOS_WATCH(GetValue(SHAPE_FUNCTION_VALUES));
        KRATOS_WATCH(GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES));
        KRATOS_WATCH(GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES));
        KRATOS_WATCH(GetValue(INTEGRATION_WEIGHT));
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "KLElement #" << Id();

    }

    ///@}

protected:

private:
    ///@name Static Member Variables
    ///@{
        
    ///@}
	///@name Member Variables
	///@{
    // curvilinear coordinate zeta (theta3)
    double mZeta;

    /**
        * Internal variables used for metric transformation
        */
    struct MetricVariables
    {
        Vector a_ab; // covariant metric
        Vector a_ab_con; // contravariant metric
        Vector curvature; //
        Matrix J; //Jacobian
        double detJ;       // not used (ML)
        Vector a1; //base vector 1
        Vector a2; //base vector 2
        Vector a3_KL; //base vector 3
        Vector a3_KL_tilde; // unnormalized base vector 3, in Kiendl (2011) a_3_tilde
        double dA; //differential area
        Vector a1_con;  // contravariant base vector 1
        Vector a2_con;  // contravariant base vector 2
        Vector Da1_D1;  // derivative of base vector 1 w.r.t. theta1
        Vector Da1_D2;  // derivative of base vector 1 w.r.t. theta2
        Vector Da2_D2;  // derivative of base vector 2 w.r.t. theta2
        Matrix H; //Hessian (second derivative of cartesian coordinates w.r.t. curvilinear coordinates)
        Matrix Q; //Transformation matrix Q from contravariant to local Cartesian basis (only for strains!!!)
        Matrix TransCartToCov; // Transformation matrix from local Cartesian to covariant basis
        Matrix TransCovToCart; // Transformation matrix from covariant to local Cartesian basis

        /**
        * The default constructor
        * @param rWorkingSpaceDimension: The size of working space dimension
        * @param rStrainSize: The size of the StrainVector
        */
        MetricVariables(const unsigned int& rWorkingSpaceDimension = 3, const unsigned int& rStrainSize = 5)
        {
            a_ab = ZeroVector(rWorkingSpaceDimension);
            a_ab_con = ZeroVector(rWorkingSpaceDimension);

            curvature = ZeroVector(rWorkingSpaceDimension);

            J = ZeroMatrix(rWorkingSpaceDimension, 2);
            detJ = 1.0;

            a1 = ZeroVector(rWorkingSpaceDimension);
            a2 = ZeroVector(rWorkingSpaceDimension);
            a3_KL = ZeroVector(rWorkingSpaceDimension);
            a3_KL_tilde = ZeroVector(rWorkingSpaceDimension);

            dA = 1.0;

            a1_con = ZeroVector(rWorkingSpaceDimension);
            a2_con = ZeroVector(rWorkingSpaceDimension);

            Da1_D1 = ZeroVector(rWorkingSpaceDimension);
            Da1_D2 = ZeroVector(rWorkingSpaceDimension);
            Da2_D2 = ZeroVector(rWorkingSpaceDimension);

            H = ZeroMatrix(rWorkingSpaceDimension, rWorkingSpaceDimension);
            Q = ZeroMatrix(rStrainSize, rStrainSize);
            TransCartToCov = ZeroMatrix(rStrainSize, rStrainSize);
            TransCovToCart = ZeroMatrix(rStrainSize, rStrainSize);
        }
    };

    /**
    * Internal variables used in the constitutive equations
    */
    struct ConstitutiveVariables
    {
        Vector E; //strain
        Vector S; //stress
        Matrix D; //constitutive matrix

        /**
        * The default constructor
        * @param StrainSize: The size of the strain vector in Voigt notation
        */
        ConstitutiveVariables(const unsigned int& rStrainSize)
        {
            E = ZeroVector(rStrainSize);
            S = ZeroVector(rStrainSize);
            D = ZeroMatrix(rStrainSize, rStrainSize);
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
        Matrix B23;
        Matrix B13;

        /**
        * The default constructor
        */
        SecondVariations(const unsigned int& mat_size)
        {
            B11 = ZeroMatrix(mat_size, mat_size);
            B22 = ZeroMatrix(mat_size, mat_size);
            B12 = ZeroMatrix(mat_size, mat_size);
            B23 = ZeroMatrix(mat_size, mat_size);
            B13 = ZeroMatrix(mat_size, mat_size);
        }

        /**
         * operator for addition (+)
         */
        SecondVariations operator+ (const SecondVariations& rSecondVariations)
        {
            KRATOS_TRY

            if (B11.size1() != rSecondVariations.B11.size1()){
                KRATOS_WATCH("Addition of SecondVariations of different size.")     // ML
                KRATOS_ERROR << "Addition of SecondVariations of different size." << std::endl;
            }
            
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

    /**
     * @brief Informations regarding the Gauss-quadrature in thickness direction
     */
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
                KRATOS_WATCH("Desired number of Gauss-Points unlogical or not implemented. You can choose 3 Gauss-Points.")     // ML
                KRATOS_ERROR << "Desired number of Gauss-Points unlogical or not implemented. You can choose 3 Gauss-Points." << std::endl;
            }
            
        }
    };

    // here the number of Gauss-Points over the thickness can be determined
    GaussQuadratureThickness mGaussQuadratureThickness = GaussQuadratureThickness(3);
    ///@}
    ///@name Operations
    ///@{
    /**
        * Calculation of the Material Stiffness Matrix. Km = B^T * D *B
        */
    void CalculateAndAddKm(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& B,
        const Matrix& D,
        const double& rIntegrationWeight );

    /**
     * @brief The method calculates and adds the non-linear part of the stiffness matrix
     * @param SecondVariationsStrain are seperated into a membrane and a curvature part
     * @param SD = stress (bending or membrane stress respectively)
     */
    void CalculateAndAddNonlinearKm(
        Matrix& rLeftHandSideMatrix,
        const SecondVariations& SecondVariationsStrain,
        const Vector& SD,
        const double& rIntegrationWeight);

    void CalculateMetric( MetricVariables& rMetric);
    
    /**
     * @brief Function determines the values of the shear dofs w_1 and w_2 and calculates the shear difference vector
     * @detail Reissner-Mindlin shell with hierarchic rotations (Oesterle 2018)
     * @param rw = shear difference vector
     */
    void CalculateShearDifferenceVector(
        array_1d<double, 3>& rw,
        array_1d<double, 3>& rDw_D1,
        array_1d<double, 3>& rDw_D2,
        array_1d<double, 2>& rw_alpha,
        Matrix& rDw_alpha_Dbeta,
        const MetricVariables& rActualMetric);
    
    /**
     * @brief Calculation of the base vectors of the shell body (in contrast to the mid-surface) for the initial configuration
     * @detail A linearized metric (g_alpha = a_alpha + zeta * Da3_Dalpha) is assumed
     */
    void CalculateInitialBaseVectorsGLinearized(
        array_1d<double, 3>& rG1,
        array_1d<double, 3>& rG2);

    /**
     * @brief Calculation of the base vectors of the shell body (in contrast to the mid-surface) for the actual configuration
     * @detail A linearized metric (g_alpha = a_alpha + zeta * Da3_Dalpha) is assumed
     * @param rw = shear difference vector
     */
    void CalculateActualBaseVectorsgLinearized(
        const MetricVariables& rActualMetric,
        const Vector& rw,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        array_1d<double, 3>& rg1,
        array_1d<double, 3>& rg2,
        array_1d<double, 3>& rg3);

    /**
     * @brief Calculates deformation gradient F for a Gauss point
     * @param rG1, rG2 = base vectors of the shell body of the reference configuration (G3=A3)
     * @param rg1, rg2, rg3 = base vectors of the shell body of the actual configuration
     */
    void CalculateDeformationGradient(
        const array_1d<double, 3> rG1,
        const array_1d<double, 3> rG2,
        const array_1d<double, 3> rg1,
        const array_1d<double, 3> rg2,
        const array_1d<double, 3> rg3,
        Matrix& rF,
        double& rdetF);

    /**
    * This functions updates the constitutive variables
    * @param rActualMetric: The actual metric
    * @param rThisConstitutiveVariables: The constitutive variables to be calculated
    * @param rValues: The CL parameters
    * @param ThisStressMeasure: The stress measure considered
    */
    void CalculateConstitutiveVariables(
        const MetricVariables& rActualMetric,
        const Vector& rw,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure);

    void CalculateStrain(
        array_1d<double, 5>& rStrainVector,
        const Vector& rgab,
        const Vector& rCurvature);

    void CalculateStrainRM(
        array_1d<double, 5>& rStrainVectorRM,
        const Vector& rw,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        const Vector& rg1,
        const Vector& rg2);

    void TransformationCurvilinearStrainSize5ToCartesianStrainSize6(
        const Vector& rCurvilinearStrain,
        Vector& rCartesianStrain);

	void CalculateB(
		Matrix& rB,
		const MetricVariables& rMetric);
    
    void CalculateSecondVariations(
        SecondVariations& rSecondVariations,
        const MetricVariables& rMetric);

    void CalculateVariationsRM(        
        Matrix& rB,
        SecondVariations& rSecondVariations,
        const Vector& rw,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        const Vector& rw_alpha,
        const Matrix& rDw_alpha_Dbeta,
        const MetricVariables& rActualMetric,
        const bool& rCalculateStiffnessMatrixFlag);
 
    /**
     * @brief Stress recovery
     */
    void Calculate(
        const Variable<double>& rVariable,
        double& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;
    
    unsigned int mcount = 0.0;
    ///@}

    ///@}
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

};     // Class IgaShell5pElement
///@}

}  // namespace Kratos.

#endif // KRATOS_IGA_SHELL_5p_ELEMENT_H_INCLUDED  defined