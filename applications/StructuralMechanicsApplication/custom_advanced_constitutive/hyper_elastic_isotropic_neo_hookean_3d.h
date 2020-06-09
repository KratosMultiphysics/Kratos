// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined (KRATOS_HYPER_ELASTIC_ISOTROPIC_NEO_HOOKEAN_3D_LAW_H_INCLUDED)
#define  KRATOS_HYPER_ELASTIC_ISOTROPIC_NEO_HOOKEAN_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
/**
 * @class HyperElasticIsotropicNeoHookean3D
 * @ingroup StructuralMechanicsApplication
 * @brief This law defines an hyperelastic material according to the NeoHookean formulation for 3D  cases
 * @details A neo-Hookean solid is a hyperelastic material model, similar to Hooke's law, that can be used for predicting the nonlinear stress-strain behavior of materials undergoing large deformations. The model was proposed by Ronald Rivlin in 1948. In contrast to linear elastic materials, the stress-strain curve of a neo-Hookean material is not linear. Instead, the relationship between applied stress and strain is initially linear, but at a certain point the stress-strain curve will plateau. The neo-Hookean model does not account for the dissipative release of energy as heat while straining the material and perfect elasticity is assumed at all stages of deformation. he neo-Hookean model is based on the statistical thermodynamics of cross-linked polymer chains and is usable for plastics and rubber-like substances.
 * More info https://en.wikipedia.org/wiki/Neo-Hookean_solid
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) HyperElasticIsotropicNeoHookean3D
    : public ConstitutiveLaw
{
public:

    ///@name Type Definitions
    ///@{

    /// The definition of the process info
    typedef ProcessInfo ProcessInfoType;

    /// The definition of the CL base  class
    typedef ConstitutiveLaw    BaseType;

    /// The definition of the size type
    typedef std::size_t        SizeType;

    /// The definition of the index type
    typedef std::size_t       IndexType;

    /// Static definition of the dimension
    static constexpr SizeType Dimension = 3;

    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = 6;

    /// Pointer definition of HyperElasticIsotropicNeoHookean3D
    KRATOS_CLASS_POINTER_DEFINITION( HyperElasticIsotropicNeoHookean3D );

    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    HyperElasticIsotropicNeoHookean3D();

    /**
     * @brief Copy constructor.
     */
    HyperElasticIsotropicNeoHookean3D (const HyperElasticIsotropicNeoHookean3D& rOther);

    /**
     * @brief Destructor.
     */
    ~HyperElasticIsotropicNeoHookean3D() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Clone operation
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * @brief This function is designed to be called once to check compatibility with element
     * @param rFeatures The Features of the law
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * @brief Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
        return Dimension;
    };

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return VoigtSize;
    };

    /**
     * @brief Returns the expected strain measure of this constitutive law (by default Green-Lagrange)
     * @return the expected strain measure
     */
    StrainMeasure GetStrainMeasure() override
    {
        return StrainMeasure_GreenLagrange;
    }

    /**
     * @brief Returns the stress measure of this constitutive law (by default 2st Piola-Kirchhoff stress in voigt notation)
     * @return the expected stress measure
     */
    StressMeasure GetStressMeasure() override
    {
        return StressMeasure_PK2;
    }

    /**
     * @brief Computes the material response: PK1 stresses and algorithmic ConstitutiveMatrix
     * @param rValues The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters & rValues) override;

    /**
     * @brief Computes the material response: PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK2 (ConstitutiveLaw::Parameters & rValues) override;

    /**
     * @brief Computes the material response: Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters & rValues) override;

    /**
     * @brief Computes the material response: Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters & rValues) override;

    /**
     * @brief Initialize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK2 (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void InitializeMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues) override;

    /**
      * @brief Updates the material response: Cauchy stresses and Internal Variables
      * @param rValues The Internalvalues of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponsePK1 (ConstitutiveLaw::Parameters & rValues) override;

    /**
      * @brief Updates the material response: Cauchy stresses and Internal Variables
      * @param rValues The Internalvalues of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponsePK2 (ConstitutiveLaw::Parameters & rValues) override;

    /**
      * @brief Updates the material response: Cauchy stresses and Internal Variables
      * @param rValues The Internalvalues of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponseKirchhoff (ConstitutiveLaw::Parameters & rValues)  override;

    /**
      * @brief Updates the material response: Cauchy stresses and Internal Variables
      * @param rValues The Internalvalues of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponseCauchy (ConstitutiveLaw::Parameters & rValues) override;

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresInitializeMaterialResponse() override
    {
        return false;
    }

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresFinalizeMaterialResponse() override
    {
        return false;
    }

    /**
     * @brief It calculates the value of a specified variable (double)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    double& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<double>& rThisVariable,
        double& rValue
        ) override;

    /**
     * @brief It calculates the value of a specified variable (Vector)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Vector& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Vector>& rThisVariable,
        Vector& rValue
        ) override;

    /**
     * @brief It calculates the value of a specified variable (Matrix)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Matrix& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue
        ) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rMaterialProperties The properties of the material
     * @param rElementGeometry The geometry of the element
     * @param rCurrentProcessInfo The current process info instance
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief It calculates the constitutive matrix C (PK2)
     * @param rConstitutiveMatrix The constitutive matrix
     * @param rInverseCTensor The inverse right Cauchy-Green tensor
     * @param DeterminantF The determinant of the deformation gradient
     * @param LameLambda First Lame parameter
     * @param LameMu Second Lame parameter
     */
    virtual void CalculateConstitutiveMatrixPK2(
        Matrix& rConstitutiveMatrix,
        const Matrix& rInverseCTensor,
        const double DeterminantF,
        const double LameLambda,
        const double LameMu
        );

    /**
     * @brief It calculates the constitutive matrix C (Kirchoff)
     * @param rConstitutiveMatrix The constitutive matrix
     * @param DeterminantF The determinant of the deformation gradient
     * @param LameLambda First Lame parameter
     * @param LameMu Second Lame parameter
     */
    virtual void CalculateConstitutiveMatrixKirchhoff(
        Matrix& rConstitutiveMatrix,
        const double DeterminantF,
        const double LameLambda,
        const double LameMu
        );

    /**
     * @brief It calculates the PK2 stress vector
     * @param rInvCTensor The inverse of the right Cauchy-Green tensor
     * @param rStressVector The stress vector in Voigt notation
     * @param DeterminantF The determinant of the deformation gradient
     * @param LameLambda First Lame parameter
     * @param LameMu Second Lame parameter
     */
    virtual void CalculatePK2Stress(
        const Matrix& rInvCTensor,
        Vector& rStressVector,
        const double DeterminantF,
        const double LameLambda,
        const double LameMu
        );

    /**
     * @brief It calculates the Kirchoff stress vector
     * @param rBTensor The left Cauchy-Green tensor
     * @param rStressVector The stress vector in Voigt notation
     * @param DeterminantF The determinant of the deformation gradient
     * @param LameLambda First Lame parameter
     * @param LameMu Second Lame parameter
     */
    virtual void CalculateKirchhoffStress(
        const Matrix& rBTensor,
        Vector& rStressVector,
        const double DeterminantF,
        const double LameLambda,
        const double LameMu
        );

    /**
     * @brief It calculates the strain vector
     * @param rValues The Internalvalues of the law
     * @param rStrainVector The strain vector in Voigt notation
     */
    virtual void CalculateGreenLagrangianStrain(
        ConstitutiveLaw::Parameters& rValues,
        Vector& rStrainVector
        );

    /**
     * @brief Calculates the Almansi strains
     * @param rValues The Internalvalues of the law
     * @param rStrainVector The strain vector in Voigt notation
     */
    virtual void CalculateAlmansiStrain(
        ConstitutiveLaw::Parameters& rValues,
        Vector& rStrainVector
        );

    ///@}
    ///@name Private Operations
    ///@{
    ///@}

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
    }


}; // Class HyperElasticIsotropicNeoHookean3D
}  // namespace Kratos.
#endif // KRATOS_HYPER_ELASTIC_ISOTROPIC_NEO_HOOKEAN_3D_LAW_H_INCLUDED  defined
