#if !defined(KRATOS_SURFACE_BASE_DISCRETE_ELEMENT_H_INCLUDED )
#define  KRATOS_SURFACE_BASE_DISCRETE_ELEMENT_H_INCLUDED


// System includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"

// External includes

// Project includes
#include "custom_elements/base_discrete_element.h"

namespace Kratos
{
///@name Kratos Classes
///@{
/// Short class definition.
/** Discrete surface element deals as base class for thin walled structures.
KRATOS_API(IGA_STRUCTURAL_MECHANICS_APPLICATION)
*/
class  SurfaceBaseDiscreteElement
    : public BaseDiscreteElement
{
protected:
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

            dA = 1.0;

            Matrix H = ZeroMatrix(3, 3);
            Matrix Q = ZeroMatrix(3, 3);
            Matrix T = ZeroMatrix(3, 3);
        }
    };

    /**
    * Internal variables used in the constitutive equations
    */
    struct ConstitutiveVariables
    {
        Vector StrainVector;
        Vector StressVector;
        Matrix DMembrane;

        /**
        * The default constructor
        * @param StrainSize: The size of the strain vector in Voigt notation
        */
        ConstitutiveVariables(const unsigned int& StrainSize)
        {
            StrainVector = ZeroVector(StrainSize);
            StressVector = ZeroVector(StrainSize);
            DMembrane = ZeroMatrix(StrainSize, StrainSize);
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


    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "SurfaceBaseDiscreteElement #" << Id();
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SurfaceBaseDiscreteElement #" << Id();
    }

public:
    ///@name Type Definitions
    ///@{
    /// Counted pointer of SurfaceBaseDiscreteElement
    KRATOS_CLASS_POINTER_DEFINITION(SurfaceBaseDiscreteElement);
    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
     // Constructor using an array of nodes
    SurfaceBaseDiscreteElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseDiscreteElement(NewId, pGeometry)
    {};
     // Constructor using an array of nodes with properties
    SurfaceBaseDiscreteElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
        : BaseDiscreteElement(NewId, pGeometry, pProperties)
    {};

    SurfaceBaseDiscreteElement() : BaseDiscreteElement()
    {};

    /// Destructor.
    virtual ~SurfaceBaseDiscreteElement() override
    {};

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_ERROR << "Trying to create a \"BaseDiscreteElement\"" << std::endl;
    };

    ///@}
    ///@name Operations
    ///@{

    /**
    * Called to initialize the element.
    * Must be called before any calculation is done
    */
    void Initialize() override;

    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo
    ) override;

    ///@}
protected:

    ConstitutiveLaw::Pointer mConstitutiveLaw;
    MetricVariables mInitialMetric = MetricVariables(3);

    ///@name Operations
    ///@{

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
    );

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

    /**
    * Hessian calculates the Hessian for the given system with the
    shape function derivatives DN_De.
    *
    * @param[in] DN_De derivatives of shape functions.
    * @param[out] Hessian calculated Hessian. Is always of size 3x3.
    */
    void Hessian(
        Matrix& Hessian,
        const Matrix& DDN_DDe,
        const int rDimension = 3);


    /**
    * Is called to compute the respective metric depending on the deformation.
    * @param metric: the current metric
    */
    virtual void CalculateMetric(
        MetricVariables& metric);

    /**
    * This functions updates the constitutive variables
    * @param rActual
    Metric: The actual metric
    * @param rThisConstitutiveVariables: The constitutive variables to be calculated
    * @param rValues: The CL parameters
    * @param ThisStressMeasure: The stress measure considered
    */
    virtual void CalculateConstitutiveVariables(
        MetricVariables& rActualMetric,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    );

    void CalculateStrain(
        Vector& StrainVector,
        Vector& gab,
        Vector& gab0);

    void CalculateCurvature(
        Vector& CurvatureVector,
        Vector& bv,
        Vector& bv_ref);

    void CalculateBMembrane(
        Matrix& rB,
        const MetricVariables& metric);

    void CalculateBCurvature(
        Matrix& rB,
        const MetricVariables& metric);

    void CalculateSecondVariationStrainMembrane(
        SecondVariations& rSecondVariationsStrain,
        const MetricVariables& rMetric);

    void CalculateSecondVariationStrainCurvature(
        Matrix& Strain_in_Q_coordinates11,
        Matrix& Strain_in_Q_coordinates22,
        Matrix& Strain_in_Q_coordinates12,
        Matrix& Curvature_in_Q_coordinates11,
        Matrix& Curvature_in_Q_coordinates22,
        Matrix& Curvature_in_Q_coordinates12,
        const MetricVariables& rMetric);
    ///@}

private:
    ///@name Operations
    ///@{

    ///@}

    ///@name Static Member Variables

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
        rSerializer.save("ConstitutiveLaw", mConstitutiveLaw);
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
        rSerializer.load("ConstitutiveLaw", mConstitutiveLaw);
    }
    ///@}
};     // Class SurfaceBaseDiscreteElement
///@}
}  // namespace Kratos.

#endif // KRATOS_SURFACE_BASE_DISCRETE_ELEMENT_H_INCLUDED  defined