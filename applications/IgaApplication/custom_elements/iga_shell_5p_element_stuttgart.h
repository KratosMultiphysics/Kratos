#if !defined(KRATOS_IGA_SHELL_5p_ELEMENT_STUTTGART_H_INCLUDED)
#define  KRATOS_IGA_SHELL_5p_ELEMENT_STUTTGART_H_INCLUDED


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
/** Reissner-Mindlin Shell with hierarchical shear vector. Optimized for Isogeometric Analysis by Oesterle et al..
*/
class IgaShell5pElementStuttgart
    : public BaseDiscreteElement
{
public:
    ///@name Type Definitions
    ///@{
    /// Counted pointer of IgaShell5pElementStuttgart
    KRATOS_CLASS_POINTER_DEFINITION(IgaShell5pElementStuttgart);
    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
    // Constructor using an array of nodes
    IgaShell5pElementStuttgart(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseDiscreteElement(NewId, pGeometry)
    {};
    // Constructor using an array of nodes with properties
    IgaShell5pElementStuttgart(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseDiscreteElement(NewId, pGeometry, pProperties)
    {};

    // default constructor necessary for serialization
    IgaShell5pElementStuttgart() : BaseDiscreteElement() {};

    /// Destructor.
    virtual ~IgaShell5pElementStuttgart() override
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
        return Kratos::make_shared<IgaShell5pElementStuttgart>(
            NewId, pGeom, pProperties);
    };

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_shared< IgaShell5pElementStuttgart >(
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
        Vector gab; // covariant metric
        Vector gab_con; // contravariant metric
        Vector curvature; //
        Matrix J; //Jacobian
        double  detJ;       // not used (ML)
        Vector a1; //base vector 1 of the mid-surface
        Vector a2; //base vector 2 of the mid-surface
        Vector a3_KL; //base vector 3 of the mid-surface
        Vector a3_KL_tilde; // unnormalized base vector 3 of the mid-surface, in Kiendl (2011) a_3_tilde
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

            a1 = ZeroVector(rWorkingSpaceDimension);
            a2 = ZeroVector(rWorkingSpaceDimension);
            a3_KL = ZeroVector(rWorkingSpaceDimension);
            a3_KL_tilde = ZeroVector(rWorkingSpaceDimension);

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

    MetricVariables mInitialMetric = MetricVariables(3, 5);
    ///@}
    ///@name Operations
    ///@{
    void CalculateMetric(MetricVariables& rMetric);
    
    void calc_G_ref(
        array_1d<double, 3>&      G1,            ///< erster Basisvektor (o)
        array_1d<double, 3>&      G2,            ///< zweiter Basisvektor (o)
        array_1d<double, 3>&      G3             ///< dritter Basisvektor (o)
        );

    void calc_g_act_linearisiert(
        const MetricVariables& rActualMetric,
        array_1d<double, 3>&      g1,            ///< erster Basisvektor (o)
        array_1d<double, 3>&      g2,            ///< zweiter Basisvektor (o)
        array_1d<double, 3>&      g3             ///< dritter Basisvektor (o)
        );

    void boperator_nln_linearisiert(
        Matrix&              bop,                     ///< B-Operator (o)
        array_1d<double, 5>&              Egl,                     ///< Green-Lagrange Verzerrungen (o)
        const Vector&               funct,                   ///< Ansatzfunktionen ausgewertet an xi, eta (i)
        const Matrix&               deriv,                   ///< erste Ableitungen der Ansatzfunktionen (i)
        const Matrix&               s_deriv,                 ///< zweite Ableitungen der Ansatzfunktionen (i)
        const MetricVariables& rActualMetric
        );
    
    void kgeom_linearisiert(
        Matrix&              IKg,                     ///< Integrand des geometrischen Steifigkeitsmatrix (o)
        const Vector&               S,                       ///< Zweite Piola-Kirchhoff-Spannungen (i)
        const Vector&               funct,                   ///< Ansatzfunktionen ausgewertet an xi, eta (i)
        const Matrix&               deriv,                   ///< erste Ableitungen der Ansatzfunktionen (i)
        const Matrix&               s_deriv,                 ///< zweite Ableitungen der Ansatzfunktionen (i)
        const MetricVariables& rActualMetric
        );

    void TransformationCurvilinearStrainSize5ToCartesianStrainSize6(
        const Vector rCurvilinearStrain,
        Vector& rCartesianStrain);
        
    /**
     * @brief Stress recovery
     */
    void Calculate(
        const Variable<double>& rVariable,
        double& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;
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

};     // Class IgaShell5pElementStuttgart
///@}

}  // namespace Kratos.

#endif // KRATOS_IGA_SHELL_5p_ELEMENT_STUTTGART_H_INCLUDED  defined