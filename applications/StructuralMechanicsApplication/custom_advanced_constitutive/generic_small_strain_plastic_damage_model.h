// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Sergio Jimenez
//
//

#if !defined(KRATOS_GENERIC_SMALL_STRAIN_PLASTIC_DAMAGE_MODEL_H_INCLUDED)
#define KRATOS_GENERIC_SMALL_STRAIN_PLASTIC_DAMAGE_MODEL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/elastic_isotropic_3d.h"
#include "custom_constitutive/linear_plane_strain.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

// The size type definition
typedef std::size_t SizeType;

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
 * @class GenericSmallStrainPlasticDamageModel
 * @ingroup StructuralMechanicsApplication
 * @brief This class is the base class which define the Plastic Damage model developed by Luccioni B. and Oller S.
 * @details This class considers a constitutive law integrator for the damage and another one for the plasticity process
 * @tparam TPlasticityIntegratorType and TDamageIntegratorType The constitutive law integrators considered
 * @author Alejandro Cornejo & Sergio Jimenez
 */
template <class TPlasticityIntegratorType, class TDamageIntegratorType>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) GenericSmallStrainPlasticDamageModel
    : public std::conditional<TPlasticityIntegratorType::VoigtSize == 6, ElasticIsotropic3D, LinearPlaneStrain >::type
{
public:
    ///@name Type Definitions
    ///@{

    // The index type definition
    typedef std::size_t IndexType;

    /// The define the working dimension size, already defined in the integrator
    static constexpr SizeType Dimension = TPlasticityIntegratorType::Dimension;

    /// The define the Voigt size, already defined in the  integrator
    static constexpr SizeType VoigtSize = TPlasticityIntegratorType::VoigtSize;

    /// Definition of the base class
    typedef typename std::conditional<VoigtSize == 6, ElasticIsotropic3D, LinearPlaneStrain >::type BaseType;

    /// The definition of the Voigt array type
    typedef array_1d<double, VoigtSize> BoundedArrayType;

    /// The definition of the bounded matrix type
    typedef BoundedMatrix<double, Dimension, Dimension> BoundedMatrixType;

    /// Counted pointer of GenericSmallStrainPlasticDamageModel
    KRATOS_CLASS_POINTER_DEFINITION(GenericSmallStrainPlasticDamageModel);

    /// The node definition
    typedef Node<3> NodeType;

    /// The geometry definition
    typedef Geometry<NodeType> GeometryType;

    /// The machine precision tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    struct PlasticDamageParameters {
        array_1d<double, VoigtSize> PlasticityFFLux;
        array_1d<double, VoigtSize> PlasticityGFLux;
        array_1d<double, VoigtSize> DamageYieldFLux;
        double DamageIndicator = 0.0;
        double PlasticityIndicator = 0.0;
        array_1d<double, VoigtSize> PlasticStrain;
        array_1d<double, VoigtSize> StrainVector;
        array_1d<double, VoigtSize> StressVector;
        double DamageIncrement = 0.0;
        double PlasticConsistencyIncrement = 0.0;
        double UniaxialStressPlasticity = 0.0;
        double UniaxialStressDamage = 0.0;
        double HardeningParameterDamage = 0.0;
        double DamageDissipationIncrement = 0.0;
        array_1d<double, VoigtSize> PlasticStrainIncrement;
        double CharacteristicLength = 0.0;
        double Damage = 0.0;
        double PlasticDissipation = 0.0;
        double DamageDissipation = 0.0;
        double DamageThreshold = 0.0;
        double PlasticityThreshold = 0.0;
        double PlasticDenominator = 0.0;
        double UndamagedFreeEnergy = 0.0;
    };
    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Default constructor.
    */
    GenericSmallStrainPlasticDamageModel()
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<GenericSmallStrainPlasticDamageModel<TPlasticityIntegratorType, TDamageIntegratorType>>(*this);
    }

    /**
    * Copy constructor.
    */
    GenericSmallStrainPlasticDamageModel(const GenericSmallStrainPlasticDamageModel &rOther)
        : BaseType(rOther),
          mPlasticDissipation(rOther.mPlasticDissipation),
          mThresholdPlasticity(rOther.mThresholdPlasticity),
          mPlasticStrain(rOther.mPlasticStrain),
		  mThresholdDamage(rOther.mThresholdDamage),
          mDamage(rOther.mDamage)
    {
    }

    /**
    * Destructor.
    */
    ~GenericSmallStrainPlasticDamageModel() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Computes the material response in terms of Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief This is to be called at the very beginning of the calculation (e.g. from InitializeElement) in order to initialize all relevant attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void InitializeMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues
        ) override;

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
     * @brief Finalize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Finalize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues) override;
    /**
     * Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<double> &rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Vector)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Vector> &rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Matrix)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Matrix> &rThisVariable) override;

    /**
     * @brief Sets the value of a specified variable (double)
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<double> &rThisVariable,
        const double& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<Vector> &rThisVariable,
        const Vector& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Returns the value of a specified variable (double)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& GetValue(
        const Variable<double> &rThisVariable,
        double& rValue
        ) override;

    /**
     * @brief Returns the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector& GetValue(
        const Variable<Vector> &rThisVariable,
        Vector& rValue
        ) override;

    /**
     * @brief Returns the value of a specified variable (matrix)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Matrix& GetValue(
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue
        ) override;

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresInitializeMaterialResponse() override
    {
        return true;
    }

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresFinalizeMaterialResponse() override
    {
        return true;
    }

    /**
     * @brief Returns the value of a specified variable (double)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    double& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<double>& rThisVariable,
        double& rValue) override;

    /**
     * @brief Returns the value of a specified variable (vector)
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
     * @brief Returns the value of a specified variable (matrix)
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
     * of the calculations, so to verify that nothing is missing from the input or that no common error is found.
     * @param rMaterialProperties The properties of the material
     * @param rElementGeometry The geometry of the element
     * @param rCurrentProcessInfo The current process info instance
     * @return 0 if OK, 1 otherwise
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This method works as the damage integrator im the isotropic damage CL and computes the associated parameters.
     * @details It is used to call the damage integrator and apply the particularities of this CL.
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rStrainVector The strain vector
     * @param rUniaxialStress The uniaxial stress of the damage model
     * @param rDamageThreshold The maximum uniaxial stress achieved previously
     * @param rDamageDissipation The internal variable of energy dissipation due to damage
     * @param rConstitutiveMatrix The elastic constitutive matrix
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     * @param rDamageFlux The derivative of the yield surface used for damage
     * @param rPlasticStrain The plastic component of the strain
     * @param Damage The internal variable of the damage model
     * @param DamageIncrement The increment of the internal variable of the damage model at this time step
     * @param UndamagedFreeEnergy The undamaged free energy
     * @param rHardd Hardening parameter for the damage model
     * @param rDamageDissipationIncrement Increment of the internal variable of energy dissipation due to damage
     */
    double CalculateDamageParameters(
        PlasticDamageParameters& rParameters,
        const Matrix& rElasticMatrix,
        ConstitutiveLaw::Parameters& rValues);

    /**
     * @brief This method computes the tensile/compressive indicators
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rTensileIndicatorFactor The tensile indicator
     * @param rCompressionIndicatorFactor The compressive indicator
     * @param rSumPrincipalStresses The sum of the principal stresses
     * @param rPrincipalStresses The principal stresses
     */
    void CalculateIndicatorsFactors(
        const array_1d<double, 6>& rPredictiveStressVector,
        double& rTensileIndicatorFactor,
        double& rCompressionIndicatorFactor,
        double& rSumPrincipalStresses,
        array_1d<double, 3>& rPrincipalStresses);

    /**
     * @brief This method checks the value of some of the variables used in the CL and guarantees that their value is between 0 and 1.
     * @param rInternalVariable Internal variable of the CL (damage, rDamageDissipationIncrement and rDamageDissipation)
     */
    void CheckInternalVariable(
        double& rInternalVariable);

    /**
     * @brief This method computes the increments for the damage internal variable and for the plastic consistency parameter.
     * @param rFluxDamageYield The derivative of the yield surface used for damage
     * @param rStrainVector The strain vector
     * @param Damage The internal variable of the damage model
     * @param rPlasticityFlux The derivative of the yield surface used for plasticity
     * @param rPlasticityGFlux The derivative of the potential used for plasticity
     * @param rElasticMatrix The elastic constitutive matrix
     * @param DamageIndicator The difference between the uniaxial stress and the damage threshold
     * @param PlasticityIndicator The difference between the uniaxial stress and the plasticity threshold
     * @param rPlasticStrain The plastic component of the strain
     * @param rDamageIncrement The increment of the internal variable of the damage model at this time step
     * @param rPlasticConsistencyIncrement The increment of the internal plastic consistency variable of the plasticity model at this time step
     * @param PlasticDenominator The plasticity numerical value to obtain the pastic consistency factor
     * @param UniaxialStressPlast The uniaxial stress of the plasticity model
     * @param Hardd Hardening parameter for the damage model
     * @param rDamageDissipationIncrement Increment of the internal variable of energy dissipation due to damage
     */
    void CalculateIncrementsPlasticDamageCase(
        PlasticDamageParameters& rParameters,
        const Matrix& rElasticMatrix);

    /**
     * @brief This method works as the damage integrator im the isotropic damage CL and computes the associated parameters.
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rStrainVector The strain vector
     * @param rUniaxialStress The uniaxial stress of the plasticity model
     * @param rThreshold The maximum uniaxial stress achieved previously
     * @param rPlasticDenominator The plasticity numerical value to obtain the plastic consistency factor
     * @param rFflux The derivative of the yield surface used for plasticity
     * @param rGflux The derivative of the potential used for plasticity
     * @param rPlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param rPlasticStrainIncrement The increment of the plastic strain
     * @param rConstitutiveMatrix The elastic constitutive matrix
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     * @param rPlasticStrain The plastic component of the strain
     * @param Damage The internal variable of the damage model
     */
    double CalculatePlasticParameters(
        PlasticDamageParameters& rParameters,
        const Matrix& rConstitutiveMatrix,
        ConstitutiveLaw::Parameters& rValues);

    /**
     * @brief This method computes the plastic denominator needed to obtain the plastic consistency factor
     * to compute the plastic consistency factor
     * @param rFflux The derivative of the yield surface
     * @param rGflux The derivative of the plastic potential
     * @param rConstitutiveMatrix The elastic constitutive matrix
     * @param rHardeningParameter The hardening parameter needed for the plasticity algorithm
     * @param rPlasticDenominator The plasticity numerical value to obtain the plastic consistency factor
     */
    void CalculatePlasticDenominator(
        const array_1d<double, VoigtSize>& rFFlux,
        const array_1d<double, VoigtSize>& rGFlux,
        const Matrix& rConstitutiveMatrix,
        double& rHardeningParameter,
        const double Damage,
        double& rPlasticDenominator);

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    double& GetThresholdPlasticity() { return mThresholdPlasticity; }
    double& GetPlasticDissipation() { return mPlasticDissipation; }
    Vector& GetPlasticStrain() { return mPlasticStrain; }

    void SetThresholdPlasticity(const double ThresholdPlasticity) { mThresholdPlasticity = ThresholdPlasticity; }
    void SetPlasticDissipation(const double PlasticDissipation) { mPlasticDissipation = PlasticDissipation; }
    void SetPlasticStrain(const array_1d<double, VoigtSize>& rPlasticStrain) { mPlasticStrain = rPlasticStrain; }

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

    ///@}
  private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    // Converged values
    double mPlasticDissipation = 0.0;
    double mThresholdPlasticity = 0.0;
    Vector mPlasticStrain = ZeroVector(VoigtSize);
	double mThresholdDamage = 0.0;
	double mDamage = 0.0;
	double mDamageDissipation = 0.0;

    // only for printing
    double mUniaxialStress = 0.0;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method computes the tangent tensor
     * @param rValues The constitutive law parameters and flags
     */
    void CalculateTangentTensor(ConstitutiveLaw::Parameters &rValues);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    // Serialization

    friend class Serializer;

    void save(Serializer &rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.save("PlasticDissipation", mPlasticDissipation);
        rSerializer.save("ThresholdPlasticity", mThresholdPlasticity);
        rSerializer.save("PlasticStrain", mPlasticStrain);
        rSerializer.save("ThresholdDamage", mThresholdDamage);
        rSerializer.save("Damage", mDamage);
        rSerializer.save("DamageDissipation", mDamageDissipation);
    }

    void load(Serializer &rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.load("PlasticDissipation", mPlasticDissipation);
        rSerializer.load("ThresholdPlasticity", mThresholdPlasticity);
        rSerializer.load("PlasticStrain", mPlasticStrain);
        rSerializer.load("ThresholdDamage", mThresholdDamage);
        rSerializer.load("Damage", mDamage);
        rSerializer.load("DamageDissipation", mDamageDissipation);
    }

    ///@}

}; // Class GenericYieldSurface

} // namespace Kratos
#endif
