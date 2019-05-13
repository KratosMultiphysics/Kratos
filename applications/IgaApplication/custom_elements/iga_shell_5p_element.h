#if !defined(KRATOS_IGA_SHELL_5p_ELEMENT_H_INCLUDED)
#define  KRATOS_IGA_SHELL_5p_ELEMENT_H_INCLUDED


// System includes
#include "includes/variables.h"

// External includes

// Project includes
#include "custom_elements/surface_base_discrete_element.h"

// Application includes
#include "iga_application_variables.h"

namespace Kratos
{
///@name Kratos Classes
///@{
/// Short class definition.
/** Kirchhoff-Love Shell. Optimized for Isogeometric Analysis by Kiendl et al. .
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
    /**
        * Internal variables used for metric transformation
        */
    struct MetricVariables
    {
        Vector gab; // covariant metric
        Vector gab_con; // contravariant metric
        Vector curvature; //
        Matrix J; //Jacobian
        double  detJ;
        Vector g1; //base vector 1
        Vector g2; //base vector 2
        Vector g3; //base vector 3
        Vector g3_notnorm; // unnormalized base vector 3, in Kiendl (2011) a_3_tilde
        double dA; //differential area
        Matrix H; //Hessian (second derivative of cartesian coordinates w.r.t. curvilinear coordinates)
        Matrix Q; //Transformation matrix Q from contravariant to cartesian basis

        /**
        * The default constructor
        * @param rWorkingSpaceDimension: The size of working space dimension
        * @param rStrainSize: The size of the StrainVector
        */
        MetricVariables(const unsigned int& rWorkingSpaceDimension = 3, const unsigned int& rStrainSize = 5)
        {
            gab = ZeroVector(rWorkingSpaceDimension);
            gab_con = ZeroVector(rWorkingSpaceDimension);

            curvature = ZeroVector(rWorkingSpaceDimension);

            J = ZeroMatrix(rWorkingSpaceDimension, rWorkingSpaceDimension);
            detJ = 1.0;

            g1 = ZeroVector(rWorkingSpaceDimension);
            g2 = ZeroVector(rWorkingSpaceDimension);
            g3 = ZeroVector(rWorkingSpaceDimension);
            g3_notnorm = ZeroVector(rWorkingSpaceDimension);

            dA = 1.0;

            H = ZeroMatrix(rWorkingSpaceDimension, rWorkingSpaceDimension);
            Q = ZeroMatrix(rStrainSize, rStrainSize);
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
        * @param StrainSize: The size of the strain vector in Voigt notation
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
        const double IntegrationWeight );

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

    void CalculateMetric( MetricVariables& rMetric );
    
    /**
     * @brief Function determines the values of the shear dofs w_1 and w_2 and calculates the shear difference vector
     * @detail Reissner-Mindlin shell with hierarchic rotations (Oesterle 2018)
     */
    void CalculateShearDifferenceVector(
        Vector& rShearDifferenceVector,
        Vector& rDw_D1,
        Vector& rDw_D2,
        Vector& rw_alpha,
        Matrix& rDw_alpha_Dbeta,
        const MetricVariables& rActualMetric);

    /**
    * This functions updates the constitutive variables
    * @param rActualMetric: The actual metric
    * @param rThisConstitutiveVariables: The constitutive variables to be calculated
    * @param rValues: The CL parameters
    * @param ThisStressMeasure: The stress measure considered
    */
    void CalculateConstitutiveVariables(
        const MetricVariables& rActualMetric,
        const Vector& rShearDifferenceVector,
        const Vector& rw_alpha,
        const Matrix& rDw_alpha_Dbeta,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure);

    void CalculateStrain(
        Vector& rStrainVector,
        const Vector& rgab);

    void CalculateCurvature(
        Vector& rCurvatureVector,
        const Vector& rCurvature);

    void CalculateStrainRM(
        Vector& rStrainVectorRM,
        const Vector& rShearDifferenceVector,
        const Vector& rg1,
        const Vector& rg2);

    void CalculateCurvatureRM(
        Vector& rCurvatureVectorRM,
        const Vector& rDw_D1,
        const Vector& rDw_D2,
        const Vector& rg1,
        const Vector& rg2);

    void TransformationCurvilinearStrainSize5ToCartesianStrainSize6(
        const Vector rCurvilinearStrain,
        Vector rCartesianStrain);

	void CalculateBMembrane(
		Matrix& rB,
		const MetricVariables& rMetric);
    
	void CalculateBCurvature(
		Matrix& rB,
		const MetricVariables& rMetric);

    void CalculateSecondVariationStrainCurvature(
        SecondVariations& rSecondVariationsMembrane,
        SecondVariations& rSecondVariationsCurvature,
        const MetricVariables& rMetric);

    void CalculateVariationsRM(        
        Matrix& rBMembraneRM,
        Matrix& rBCurvatureRM,
        SecondVariations& rSecondVariationsMembraneRM,
        SecondVariations& rSecondVariationsCurvatureRM,
        const Vector& rShearDifferenceVector,
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
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo);
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