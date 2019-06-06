#if !defined(KRATOS_IGA_SHELL_5P_ELEMENT_H_INCLUDED )
#define  KRATOS_IGA_SHELL_5P_ELEMENT_H_INCLUDED


// System includes
#include "utilities/math_utils.h"

// External includes

// Project includes

// Application includes
#include "iga_application_variables.h"

#include "custom_elements/iga_base_element.h"
#include "custom_utilities/geometry_utilities/iga_geometry_utilities.h"



namespace Kratos
{
///@name Kratos Classes
///@{
/// Short class definition.
/** Kirchhoff-Love Shell. Developed for Isogeometric Analysis by Kiendl et al.
*   (https://doi.org/10.1016/j.cma.2009.08.013)
*   The formulation requires a continuity of C >= 2.
*/
class IgaShell3pElement
    : public IgaBaseElement
{
private:
    /**
    * Internal variables used for metric transformation
    */
    struct MetricVariables
    {
        array_1d<double, 3> a_ab_covariant; // covariant metric
        array_1d<double, 3> b_ab_covariant; //
        Matrix J; //Jacobian
        array_1d<double, 3> a1; //base vector 1
        array_1d<double, 3> a2; //base vector 2
        array_1d<double, 3> a3; //base vector 3 normalized
        array_1d<double, 3> a3_tilde; //not-normalized base vector 3
        double dA; //differential area
        Matrix H; //Hessian

        /**
        * The default constructor
        * @param Dimension: The size of working space dimension
        */
        MetricVariables(const unsigned int& Dimension)
        {
            a_ab_covariant = ZeroVector(Dimension);
            b_ab_covariant = ZeroVector(Dimension);

            J = ZeroMatrix(Dimension, Dimension);

            a1 = ZeroVector(Dimension);
            a2 = ZeroVector(Dimension);
            a3 = ZeroVector(Dimension);

            a3_tilde = ZeroVector(Dimension);

            dA = 1.0;

            H = ZeroMatrix(3, 3);
        }
    };

    /**
    * Internal variables used in the constitutive equations.
    * Contains of strain vector E, stress vector S and constitutive Matrix D.
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
        : IgaBaseElement(NewId, pGeometry)
    {};

    // Constructor using an array of nodes with properties
    IgaShell3pElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : IgaBaseElement(NewId, pGeometry, pProperties)
    {};

    // default constructor necessary for serialization
    IgaShell3pElement() : IgaBaseElement() {};

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
        return Kratos::make_intrusive<IgaShell3pElement>(
            NewId, pGeom, pProperties);
    };

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive< IgaShell3pElement >(
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

    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo
    );

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
        buffer << "IgaShell3pElement #" << Id();
        KRATOS_WATCH(GetValue(SHAPE_FUNCTION_VALUES));
        KRATOS_WATCH(GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES));
        KRATOS_WATCH(GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES));
        KRATOS_WATCH(GetValue(INTEGRATION_WEIGHT));
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IgaShell3pElement #" << Id();
    }

    ///@}

protected:

private:
    // Components of the metric coefficient tensor on the contravariant basis
    array_1d<double, 3> mA_ab_covariant;
    // Components of the curvature coefficient tensor on the contravariant basis
    array_1d<double, 3> mB_ab_covariant;
    // Determinant of the geometrical Jacobian.
    double mdetJ;
    // Transformation the strain tensor from the curvilinear system
    // to the local cartesian in voigt notation including a 2 in the 
    // shear part.
    Matrix mT;

    ConstitutiveLaw::Pointer mConstitutiveLaw;

    ///@name Static Member Variables
    ///@{
    ///@name Operations
    ///@{

    /* Initilaizes the constitutive law and sets the member 
    *  variable mConstitutiveLaw
    */
    void InitializeMaterial();

    void CalculateBMembrane(
        Matrix& rB,
        const MetricVariables& rMetric);

    void CalculateSecondVariationStrainMembrane(
        SecondVariations& rSecondVariationsStrain,
        const MetricVariables& rMetric);

    void CalculateBBending(
        Matrix& rB,
        const MetricVariables& rMetric);

    void CalculateSecondVariationStrainAndBending(
        SecondVariations& rSecondVariationsStrain,
        SecondVariations& rSecondVariationsCurvature,
        const MetricVariables& rMetric);

    void CalculateMetrics(
        MetricVariables& metric );

    void CalculateTransformation(
        const MetricVariables& metric,
        Matrix& T);

    void CalculateAndAddKm(
        Matrix& rLeftHandSideMatrix,
        const Matrix& B,
        const Matrix& D,
        const double IntegrationWeight);

    void CalculateAndAddNonlinearKm(
        Matrix& rLeftHandSideMatrix,
        const SecondVariations& SecondVariationsStrain,
        const Vector& S,
        const double& rIntegrationWeight);

    /**
    * This functions updates the constitutive variables
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
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    );
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

#endif // KRATOS_MESHLESS_IGA_SHELL_5P_ELEMENT_H_INCLUDED  defined