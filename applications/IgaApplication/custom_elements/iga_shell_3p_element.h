#if !defined(KRATOS_IGA_SHELL_3p_ELEMENT_H_INCLUDED)
#define  KRATOS_IGA_SHELL_3p_ELEMENT_H_INCLUDED


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
/** Kirchhoff-Love Shell. Optimized for Isogeometric Analysis by Kiendl et al. .
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
    IgaShell3pElement() : BaseDiscreteElement() {};

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
        return Kratos::make_shared<IgaShell3pElement>(
            NewId, pGeom, pProperties);
    };

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_shared< IgaShell3pElement >(
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
        Vector g3_tilde; // unnormalized base vector 3, in Kiendl (2011) a_3_tilde
        double dA; //differential area
        Vector DDg1_DD11;  // second derivative of base vector 1 w.r.t. theta1
        Vector DDg1_DD12;  // second derivative of base vector 1 w.r.t. theta1 and theta2
        Vector DDg2_DD21;  // second derivative of base vector 2 w.r.t. theta1 and theta2
        Vector DDg2_DD22;  // second derivative of base vector 2 w.r.t. theta2
        Vector Dcurvature_D1;   // derivative of curvature w.r.t. theta1
        Vector Dcurvature_D2;   // derivative of curvature w.r.t. theta2
        Matrix H; //Hessian
        Matrix Q; //Transformation matrix Q from contravariant to local Cartesian basis (only for strains!!!)
        Matrix TransCartToCov; // Transformation matrix from local Cartesian to covariant basis
        Matrix TransCovToCart; // Transformation matrix from covariant to local Cartesian basis

        /**
        * The default constructor
        * @param Dimension: The size of working space dimension
        */
        MetricVariables(const unsigned int& Dimension)
        {
            gab = ZeroVector(Dimension);
            gab_con = ZeroVector(Dimension);

            curvature = ZeroVector(Dimension);

            J = ZeroMatrix(Dimension, 2);
            detJ = 1.0;

            g1 = ZeroVector(Dimension);
            g2 = ZeroVector(Dimension);
            g3 = ZeroVector(Dimension);

            dA = 1.0;

            DDg1_DD11 = ZeroVector(Dimension);
            DDg1_DD12 = ZeroVector(Dimension);
            DDg2_DD21 = ZeroVector(Dimension);
            DDg2_DD22 = ZeroVector(Dimension);

            Dcurvature_D1 = ZeroVector(Dimension);
            Dcurvature_D2 = ZeroVector(Dimension);

            H = ZeroMatrix(3, 3);
            Q = ZeroMatrix(3, 3);
            TransCartToCov = ZeroMatrix(3, 3);
            TransCovToCart = ZeroMatrix(3, 3);
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

    MetricVariables mInitialMetric = MetricVariables(3);
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

    void CalculateAndAddNonlinearKm(
        Matrix& rLeftHandSideMatrix,
        const SecondVariations& SecondVariationsStrain,
        const Vector& SD,
        const double& rIntegrationWeight);

    void CalculateMetric( MetricVariables& rMetric );

    /**
    * This functions updates the constitutive variables
    * @param rActualMetric: The actual metric
    * @param rThisConstitutiveVariables: The constitutive variables to be calculated
    * @param rValues: The CL parameters
    * @param ThisStressMeasure: The stress measure considered
    */
    void CalculateConstitutiveVariables(
        const MetricVariables& rActualMetric,
        ConstitutiveVariables& rThisConstitutiveVariablesMembrane,
        ConstitutiveVariables& rThisConstitutiveVariablesCurvature,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    );

    void CalculateStrain(
        Vector& rStrainVector,
        const Vector& rgab);

    void CalculateCurvature(
        Vector& rCurvatureVector,
        const Vector& rCurvature);

	void CalculateBMembrane(
		Matrix& rB,
		const MetricVariables& rMetric);
    
	void CalculateBCurvature(
		Matrix& rB,
		const MetricVariables& rMetric);

    void CalculateSecondVariationStrainCurvature(
        SecondVariations& rSecondVariationsStrain,
        SecondVariations& rSecondVariationsCurvature,
        const MetricVariables& rMetric);
    
    /**
     * @brief Stress recovery
     */
    void Calculate(
        const Variable<double>& rVariable,
        double& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Calculation of the shear force, shear force = derivative of moment
     * @detail based on Carat ElementShell_NURBS_KL.cpp
     */
    void CalculateShearForce(
        array_1d<double, 2>& rq, 
        const MetricVariables& rActualMetric,
        const ConstitutiveVariables& rConstitutiveVariablesCurvature);

    void CalculateDerivativeTransformationMatrices(
        std::vector<Matrix>& rDQ_Dalpha_init,
        std::vector<Matrix>& rDTransCartToCov_Dalpha_init);
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

};     // Class IgaShell3pElement
///@}

}  // namespace Kratos.

#endif // KRATOS_IGA_SHELL_3p_ELEMENT_H_INCLUDED  defined