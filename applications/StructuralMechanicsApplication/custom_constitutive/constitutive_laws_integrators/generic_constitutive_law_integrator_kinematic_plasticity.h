// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#if !defined(KRATOS_GENERIC_CONSTITUTIVE_LAW_INTEGRATOR_KINEMATIC_PLASTICITY_H_INCLUDED)
#define KRATOS_GENERIC_CONSTITUTIVE_LAW_INTEGRATOR_KINEMATIC_PLASTICITY_H_INCLUDED

// System includes

// Project includes
#include "includes/define.h"
#include "includes/checks.h"
#include "includes/properties.h"
#include "utilities/math_utils.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

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
 * @class GenericConstitutiveLawIntegratorKinematicPlasticity
 * @ingroup StructuralMechanicsApplication
 * @brief This object integrates the predictive stress using the plasticity theory by means of
 * linear/exponential softening or hardening + softening evolution laws
 * @details The definitions of these classes is completely static, the derivation is done in a static way
 * @tparam TYieldSurfaceType The yield surface considered
 * The plasticity integrator requires the definition of the following properties:
 * - SOFTENING_TYPE: The softening behaviour considered (linear, exponential,etc...)
 * - HARDENING_CURVE: The type of considered hardening curve (linear, exponential, pure plastic, etc...)
 * - MAXIMUM_STRESS: The maximum stress that defines the exponential hardening
 * - MAXIMUM_STRESS_POSITION: The maximum stress position that defines the exponential hardening
 * - FRACTURE_ENERGY: A fracture energy-based function is used to describe strength degradation in post-peak regime
 * - YOUNG_MODULUS: It defines the relationship between stress (force per unit area) and strain (proportional deformation) in a material in the linear elasticity regime of a uniaxial deformation.
 * - YIELD_STRESS: Yield stress is the amount of stress that an object needs to experience for it to be permanently deformed. Does not require to be defined simmetrically, one YIELD_STRESS_COMPRESSION and other YIELD_STRESS_TENSION can be defined for not symmetric cases
 * @author Alejandro Cornejo & Lucia Barbu
 */
template<class TYieldSurfaceType>
class GenericConstitutiveLawIntegratorKinematicPlasticity
{
  public:
    ///@name Type Definitions
    ///@{
      
    /// The machine precision tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    /// Definition of index
    typedef std::size_t IndexType;

    /// The type of yield surface
    typedef TYieldSurfaceType YieldSurfaceType;

    /// The define the working dimension size, already defined in the yield surface
    static constexpr SizeType Dimension = YieldSurfaceType::Dimension;

    /// The define the Voigt size, already defined in the yield surface
    static constexpr SizeType VoigtSize = YieldSurfaceType::VoigtSize;
    
    /// The type of plastic potential
    typedef typename YieldSurfaceType::PlasticPotentialType PlasticPotentialType;

    /// Counted pointer of GenericConstitutiveLawIntegratorKinematicPlasticity
    KRATOS_CLASS_POINTER_DEFINITION(GenericConstitutiveLawIntegratorKinematicPlasticity);

    ///@}
    ///@name  Enum's
    ///@{

    enum class HardeningCurveType
    {
        LinearSoftening = 0,
        ExponentialSoftening = 1,
        InitialHardeningExponentialSoftening = 2,
        PerfectPlasticity = 3,
        CurveFittingHardening = 4
    };

    enum class KinematicHardeningType
    {
        LinearKinematicHardening = 0,
        AmstrongFrederickKinematicHardening = 1,
        AraujoVoyiadjisKinematicHardening = 2
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor
    GenericConstitutiveLawIntegratorKinematicPlasticity()
    {
    }

    /// Copy constructor
    GenericConstitutiveLawIntegratorKinematicPlasticity(GenericConstitutiveLawIntegratorKinematicPlasticity const &rOther)
    {
    }

    /// Assignment operator
    GenericConstitutiveLawIntegratorKinematicPlasticity &operator=(GenericConstitutiveLawIntegratorKinematicPlasticity const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~GenericConstitutiveLawIntegratorKinematicPlasticity()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method integrates the predictive stress vector with the CL using differents evolution laws using the backward euler scheme
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rStrainVector The equivalent strain vector of that integration point
     * @param rUniaxialStress The equivalent uniaxial stress
     * @param rThreshold The maximum uniaxial stress of the linear behaviour
     * @param rPlasticDenominator The plasticity numerical value to obtain the pastic consistency factor
     * @param rFflux The derivative of the yield surface
     * @param rGflux The derivative of the plastic potential
     * @param rPlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param rPlasticStrainIncrement The increment of plastic strain of this time step
     * @param rConstitutiveMatrix The elastic constitutive matrix
     * @param rPlasticStrain The elastic constitutive matrix
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void IntegrateStressVector(
        array_1d<double, VoigtSize>& rPredictiveStressVector,
        Vector& rStrainVector,
        double& rUniaxialStress,
        double& rThreshold,
        double& rPlasticDenominator,
        array_1d<double, VoigtSize>& rFflux,
        array_1d<double, VoigtSize>& rGflux,
        double& rPlasticDissipation,
        array_1d<double, VoigtSize>& rPlasticStrainIncrement,
        const Matrix& rConstitutiveMatrix,
        Vector& rPlasticStrain,
        ConstitutiveLaw::Parameters& rValues,
        const double CharacteristicLength,
        Vector& rBackStressVector,
        const Vector& rPreviousStressVector
        )
    {
        bool is_converged = false;
        IndexType iteration = 0, max_iter = 100;
        array_1d<double, VoigtSize> delta_sigma;
        double plastic_consistency_factor_increment, F;

        // Backward Euler
        while (is_converged == false && iteration <= max_iter) {
            F = rUniaxialStress - rThreshold;
            plastic_consistency_factor_increment = F * rPlasticDenominator;
            //if (plastic_consistency_factor_increment < 0.0) plastic_consistency_factor_increment = 0.0; // NOTE: It could be useful, maybe
            noalias(rPlasticStrainIncrement) = plastic_consistency_factor_increment * rGflux;
            noalias(rPlasticStrain) += rPlasticStrainIncrement;
            noalias(delta_sigma) = prod(rConstitutiveMatrix, rPlasticStrainIncrement);

            noalias(rPredictiveStressVector) -= delta_sigma;
            CalculateAndSubstractBackStress(rPredictiveStressVector, rValues, rPreviousStressVector,
                                            rPlasticStrainIncrement, rBackStressVector);
            
            CalculatePlasticParameters(rPredictiveStressVector, rStrainVector, rUniaxialStress, rThreshold,
                                       rPlasticDenominator, rFflux, rGflux, rPlasticDissipation, rPlasticStrainIncrement,
                                       rConstitutiveMatrix, rValues, CharacteristicLength, rPlasticStrain, rBackStressVector);

            F = rUniaxialStress - rThreshold;

            if (std::abs(F) <= std::abs(1.0e-4 * rThreshold)) { // Has converged
                is_converged = true;
            } else {
                iteration++;
            }
        }
        KRATOS_WARNING_IF("Backward Euler Plasticity", iteration > max_iter) << "Maximum number of iterations in plasticity loop reached..." << std::endl;
    }

    /**
     * @brief This method calculates all the plastic parameters required for the integration of the PredictiveStressVector
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rStrainVector The equivalent strain vector of that integration point
     * @param rUniaxialStress The equivalent uniaxial stress
     * @param rThreshold The maximum uniaxial stress of the linear behaviour
     * @param rPlasticDenominator The plasticity numerical value to obtain the pastic consistency factor
     * @param rFflux The derivative of the yield surface
     * @param rGflux The derivative of the plastic potential
     * @param rPlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param rPlasticStrainIncrement The increment of plastic strain of this time step
     * @param rConstitutiveMatrix The elastic constitutive matrix
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculatePlasticParameters(
        array_1d<double, VoigtSize>& rPredictiveStressVector,
        Vector& rStrainVector,
        double& rUniaxialStress,
        double& rThreshold,
        double& rPlasticDenominator,
        array_1d<double, VoigtSize>& rFflux,
        array_1d<double, VoigtSize>& rGflux,
        double& rPlasticDissipation,
        array_1d<double, VoigtSize>& rPlasticStrainIncrement,
        const Matrix& rConstitutiveMatrix,
        ConstitutiveLaw::Parameters& rValues,
        const double CharacteristicLength,
        const Vector& rPlasticStrain,
        const Vector& rBackStressVector
        )
    {
        array_1d<double, VoigtSize> deviator = ZeroVector(6);
        array_1d<double, VoigtSize> h_capa = ZeroVector(6);
        double J2, tensile_indicator_factor, compression_indicator_factor, slope, hardening_parameter, equivalent_plastic_strain;

        YieldSurfaceType::CalculateEquivalentStress( rPredictiveStressVector, rStrainVector, rUniaxialStress, rValues);
        const double I1 = rPredictiveStressVector[0] + rPredictiveStressVector[1] + rPredictiveStressVector[2];
        ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(rPredictiveStressVector, I1, deviator, J2);
        CalculateFFluxVector(rPredictiveStressVector, deviator, J2, rFflux, rValues);
        CalculateGFluxVector(rPredictiveStressVector, deviator, J2, rGflux, rValues);
        CalculateIndicatorsFactors(rPredictiveStressVector, tensile_indicator_factor,compression_indicator_factor);
        CalculatePlasticDissipation(rPredictiveStressVector, tensile_indicator_factor,compression_indicator_factor, rPlasticStrainIncrement,rPlasticDissipation, h_capa, rValues, CharacteristicLength);
        CalculateEquivalentPlasticStrain(rPredictiveStressVector, rUniaxialStress, rPlasticStrain, tensile_indicator_factor, rValues, equivalent_plastic_strain);
        CalculateEquivalentStressThreshold( rPlasticDissipation, tensile_indicator_factor,compression_indicator_factor, rThreshold, slope, rValues, equivalent_plastic_strain);
        CalculateHardeningParameter(rFflux, slope, h_capa, hardening_parameter);
        CalculatePlasticDenominator(rFflux, rGflux, rConstitutiveMatrix, hardening_parameter, rPlasticDenominator, rBackStressVector, rValues);
    }

    /**
     * @brief This method calculates the derivative of the yield surface
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rDeviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the deviatoric part of the stress vector
     * @param rFFluxVector The derivative of the yield surface
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateFFluxVector(
        const array_1d<double, VoigtSize>& rPredictiveStressVector,
        const array_1d<double, VoigtSize>& rDeviator,
        const double J2,
        array_1d<double, VoigtSize>& rFFluxVector,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        YieldSurfaceType::CalculateYieldSurfaceDerivative(rPredictiveStressVector, rDeviator, J2, rFFluxVector, rValues);
    }

    /**
     * @brief This method computes the back stress for the kinematic plasticity
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rKinematicParameters The oarameters required for the models
     * @param rPreviousStressVector The previous stress vector
     * @param rPlasticStrainIncrement The plastric strain increment of this step
     * @param rBackStressVector The back stress for the kinematic plasticity
     */
    static void CalculateAndSubstractBackStress(
        array_1d<double, VoigtSize>& rPredictiveStressVector,
        ConstitutiveLaw::Parameters& rValues,
        const Vector& rPreviousStressVector,
        const Vector& rPlasticStrainIncrement,
        Vector& rBackStressVector
        )
    {
        const Vector& kinematic_parameters = rValues.GetMaterialProperties()[KINEMATIC_PLASTICITY_PARAMETERS];
        const int kinematic_hardening_type = rValues.GetMaterialProperties()[KINEMATIC_HARDENING_TYPE];

        switch (static_cast<KinematicHardeningType>(kinematic_hardening_type))
        {
            double pDot, denominator, dot_product_dp;
            case KinematicHardeningType::LinearKinematicHardening:
                KRATOS_ERROR_IF(kinematic_parameters.size() == 0) << "Kinematic Parameters not defined..." << std::endl;
                rBackStressVector += 2.0 / 3.0 * kinematic_parameters[0] * rPlasticStrainIncrement;
                break;

            case KinematicHardeningType::AmstrongFrederickKinematicHardening:
                KRATOS_ERROR_IF(kinematic_parameters.size() < 2) << "Kinematic Parameters not defined..." << std::endl;
                dot_product_dp = 0.0;
                for (int i = 0; i < rPlasticStrainIncrement.size(); ++i) {
                    dot_product_dp += rPlasticStrainIncrement[i]*rPlasticStrainIncrement[i];
                }
                pDot = std::sqrt(2.0 / 3.0 * dot_product_dp);
                denominator = 1.0 + (kinematic_parameters[1] * pDot);
                rBackStressVector += (2.0 / 3.0 * kinematic_parameters[0] * rPlasticStrainIncrement) / denominator;
                //KRATOS_WATCH(rBackStressVector)
				break;

            case KinematicHardeningType::AraujoVoyiadjisKinematicHardening:
                KRATOS_ERROR_IF(kinematic_parameters.size() != 3) << "Kinematic Parameters not defined..." << std::endl;
                dot_product_dp = 0.0;
                for (int i = 0; i < rPlasticStrainIncrement.size(); ++i) {
                    dot_product_dp += rPlasticStrainIncrement[i]*rPlasticStrainIncrement[i];
                }
                pDot = std::sqrt(2.0 / 3.0 * dot_product_dp);
                denominator = 1.0 + (kinematic_parameters[1] * pDot);
                if (pDot != 0.0) {
                    rBackStressVector += (2.0 / 3.0 * kinematic_parameters[0] * rPlasticStrainIncrement) / denominator;
                } else {
                    const Vector& delta_stress = rPredictiveStressVector - rPreviousStressVector;
                    rBackStressVector += ((2.0 / 3.0 * kinematic_parameters[0] * rPlasticStrainIncrement) +
                                         kinematic_parameters[2] * delta_stress) / denominator;
                }
                break;

            default:
                KRATOS_ERROR << " The Kinematic hardening type of plasticity is not set or wrong..." << kinematic_hardening_type << std::endl;
                break;
        }
        rPredictiveStressVector -= rBackStressVector;
    }

    /**
     * @brief This method calculates the derivative of the plastic potential
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rDeviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the deviatoric part of the stress vector
     * @param rGFluxVector The derivative of the yield surface
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateGFluxVector(
        const array_1d<double, VoigtSize>& rPredictiveStressVector,
        const array_1d<double, VoigtSize>& rDeviator,
        const double J2,
        array_1d<double, VoigtSize>& rGFluxVector,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        YieldSurfaceType::CalculatePlasticPotentialDerivative(rPredictiveStressVector, rDeviator, J2, rGFluxVector, rValues);
    }

    /**
     * @brief This method computes the tensile/compressive indicators
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rTensileIndicatorFactor The tensile indicator
     * @param rCompressionIndicatorFactor The compressive indicator
     */
    static void CalculateIndicatorsFactors(
        const array_1d<double, VoigtSize>& rPredictiveStressVector,
        double& rTensileIndicatorFactor,
        double& rCompressionIndicatorFactor
        )
    {
        // We do an initial check
        if (norm_2(rPredictiveStressVector) < 1.0e-8) {
            rTensileIndicatorFactor = 0.5;
            rCompressionIndicatorFactor = 0.5;
            return;
        }

        // We proceed as usual
        array_1d<double, Dimension> principal_stresses = ZeroVector(Dimension);
        ConstitutiveLawUtilities<VoigtSize>::CalculatePrincipalStresses(principal_stresses, rPredictiveStressVector);

        double suma = 0.0, sumb = 0.0, sumc = 0.0;
        double aux_sa;

        for (IndexType i = 0; i < Dimension; ++i) {
            aux_sa = std::abs(principal_stresses[i]);
            suma += aux_sa;
            sumb += 0.5 * (principal_stresses[i] + aux_sa);
            sumc += 0.5 * (-principal_stresses[i] + aux_sa);
        }

        if (std::abs(suma) > tolerance) {
            rTensileIndicatorFactor = sumb / suma;
            rCompressionIndicatorFactor = sumc / suma;
        } else {
            rTensileIndicatorFactor = sumb;
            rCompressionIndicatorFactor = sumc;
        }

        // Final check
        if ((std::abs(rTensileIndicatorFactor) + std::abs(rCompressionIndicatorFactor)) < tolerance) {
            rTensileIndicatorFactor = 0.5;
            rCompressionIndicatorFactor = 0.5;
            return;
        }
    }

    /**
     * @brief This method computes the plastic dissipation of the plasticity model
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param TensileIndicatorFactor The tensile indicator
     * @param CompressionIndicatorFactor The compressive indicator
     * @param PlasticStrainIncrement The increment of plastic strain of this time step
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param rHCapa The slope of the PlasticDiss-Threshold curve
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculatePlasticDissipation(
        const array_1d<double, VoigtSize>& rPredictiveStressVector,
        const double TensileIndicatorFactor,
        const double CompressionIndicatorFactor,
        const Vector& PlasticStrainInc,
        double& rPlasticDissipation,
        array_1d<double, VoigtSize>& rHCapa,
        ConstitutiveLaw::Parameters& rValues,
        const double CharacteristicLength
        )
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();

        const double young_modulus = r_material_properties[YOUNG_MODULUS];
        const bool has_symmetric_yield_stress = r_material_properties.Has(YIELD_STRESS);
        const double yield_compression = has_symmetric_yield_stress ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_COMPRESSION];
        const double yield_tension = has_symmetric_yield_stress ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_TENSION];
        const double n = yield_compression / yield_tension;
        const double fracture_energy_tension = r_material_properties[FRACTURE_ENERGY]; // Frac energy in tension
        const double fracture_energy_comprression = r_material_properties[FRACTURE_ENERGY] * std::pow(n, 2); // Frac energy in compression

        const double characteristic_fracture_energy_tension = fracture_energy_tension / CharacteristicLength;
        const double characteristic_fracture_energy_compression = fracture_energy_comprression / CharacteristicLength;

        const double hlim = 2.0 * young_modulus * characteristic_fracture_energy_compression / (std::pow(yield_compression, 2));
        KRATOS_ERROR_IF(CharacteristicLength > hlim) << "The Fracture Energy is to low: " << characteristic_fracture_energy_compression << std::endl;

        double constant0 = 0.0, constant1 = 0.0, dplastic_dissipation = 0.0;
        if (characteristic_fracture_energy_tension > 0.000001) {
            constant0 = TensileIndicatorFactor / characteristic_fracture_energy_tension;
            constant1 = CompressionIndicatorFactor / characteristic_fracture_energy_compression;
        }
        const double constant = constant0 + constant1;

        for (IndexType i = 0; i < VoigtSize; ++i) {
            rHCapa[i] = constant * rPredictiveStressVector[i];
            dplastic_dissipation += rHCapa[i] * PlasticStrainInc[i];
        }

        if (dplastic_dissipation < 0.0 || dplastic_dissipation > 1.0)
            dplastic_dissipation = 0.0;

        rPlasticDissipation += dplastic_dissipation;
        if (rPlasticDissipation >= 1.0)
            rPlasticDissipation = 0.9999;
        else if (rPlasticDissipation < 0.0)
            rPlasticDissipation = 0.0;

        // We add a check
        KRATOS_DEBUG_ERROR_IF(std::isnan(rPlasticDissipation)) << "rPlasticDissipation is nan" << std::endl;
    }

    /**
     * @brief This method computes the uniaxial threshold that differentiates the elastic-plastic behaviour
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param TensileIndicatorFactor The tensile indicator
     * @param CompressionIndicatorFactor The compressive indicator
     * @param rEquivalentStressThreshold The maximum uniaxial stress of the linear behaviour
     * @param rSlope The slope of the PlasticDiss-Threshold curve
     * @param rValues Parameters of the constitutive law
     * @param EquivalentPlasticStrain The equivalent plastic strain
     */
    static void CalculateEquivalentStressThreshold(
        const double PlasticDissipation,
        const double TensileIndicatorFactor,
        const double CompressionIndicatorFactor,
        double& rEquivalentStressThreshold,
        double& rSlope,
        ConstitutiveLaw::Parameters& rValues,
        const double EquivalentPlasticStrain
        )
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();
        const int curve_type = r_material_properties[HARDENING_CURVE];
        BoundedVector<double, 2> slopes, eq_thresholds;

        for (IndexType i = 0; i < 2; ++i) { // i:0 Tension ; i:1 compression
            switch (static_cast<HardeningCurveType>(curve_type))
            {
              case HardeningCurveType::LinearSoftening:
                CalculateEquivalentStressThresholdHardeningCurveLinearSoftening(
                    PlasticDissipation, TensileIndicatorFactor,
                    CompressionIndicatorFactor, eq_thresholds[i], slopes[i],
                    rValues);
                break;

              case HardeningCurveType::ExponentialSoftening:
                CalculateEquivalentStressThresholdHardeningCurveExponentialSoftening(
                    PlasticDissipation, TensileIndicatorFactor,
                    CompressionIndicatorFactor, eq_thresholds[i], slopes[i],
                    rValues);
                break;

              case HardeningCurveType::InitialHardeningExponentialSoftening:
                CalculateEquivalentStressThresholdHardeningCurveInitialHardeningExponentialSoftening(
                    PlasticDissipation, TensileIndicatorFactor,
                    CompressionIndicatorFactor, eq_thresholds[i], slopes[i],
                    rValues);
                break;

              case HardeningCurveType::PerfectPlasticity:
                CalculateEquivalentStressThresholdHardeningCurvePerfectPlasticity(
                    PlasticDissipation, TensileIndicatorFactor,
                    CompressionIndicatorFactor, eq_thresholds[i], slopes[i],
                    rValues);
                break;

              case HardeningCurveType::CurveFittingHardening: // Only for VonMises + Tresca!
                CalculateEquivalentStressThresholdCurveFittingHardening(
                    PlasticDissipation, TensileIndicatorFactor,
                    CompressionIndicatorFactor, eq_thresholds[i], slopes[i],
                    rValues, EquivalentPlasticStrain);
                break;

                // Add more cases...

              default:
                KRATOS_ERROR << " The HARDENING_CURVE of plasticity is not set or wrong..." << curve_type << std::endl;
                break;
            }
        }

        rEquivalentStressThreshold = TensileIndicatorFactor * eq_thresholds[0] + CompressionIndicatorFactor * eq_thresholds[1];
        rSlope = rEquivalentStressThreshold * ((TensileIndicatorFactor * slopes[0] / eq_thresholds[0]) + (CompressionIndicatorFactor * slopes[1] / eq_thresholds[1]));
        KRATOS_DEBUG_ERROR_IF(rEquivalentStressThreshold < tolerance) << "Threshold set to zero. TensileIndicatorFactor: " << TensileIndicatorFactor << " eq_thresholds[0]: " << eq_thresholds[0] << " CompressionIndicatorFactor: " << CompressionIndicatorFactor << " eq_thresholds[1]:" << eq_thresholds[1] << std::endl;
    }

    /**
     * @brief This method computes the uniaxial threshold using a linear softening
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param TensileIndicatorFactor The tensile indicator
     * @param CompressionIndicatorFactor The compressive indicator
     * param rEquivalentStressThreshold The maximum uniaxial stress of the linear behaviour
     * @param rSlope The slope of the PlasticDiss-Threshold curve
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateEquivalentStressThresholdHardeningCurveLinearSoftening(
        const double PlasticDissipation,
        const double TensileIndicatorFactor,
        const double CompressionIndicatorFactor,
        double& rEquivalentStressThreshold,
        double& rSlope,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        double initial_threshold;
        GetInitialUniaxialThreshold(rValues, initial_threshold);

        rEquivalentStressThreshold = initial_threshold * std::sqrt(1.0 - PlasticDissipation);
        rSlope = -0.5 * (std::pow(initial_threshold, 2.0) / (rEquivalentStressThreshold));
    }

    /**
     * @brief This method computes the uniaxial threshold using a exponential softening
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param TensileIndicatorFactor The tensile indicator
     * @param CompressionIndicatorFactor The compressive indicator
     * param rEquivalentStressThreshold The maximum uniaxial stress of the linear behaviour
     * @param rSlope The slope of the PlasticDiss-Threshold curve
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateEquivalentStressThresholdHardeningCurveExponentialSoftening(
        const double PlasticDissipation,
        const double TensileIndicatorFactor,
        const double CompressionIndicatorFactor,
        double& rEquivalentStressThreshold,
        double& rSlope,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        double initial_threshold;
        GetInitialUniaxialThreshold(rValues, initial_threshold);

        rEquivalentStressThreshold = initial_threshold * (1.0 - PlasticDissipation);
        rSlope = -0.5 * initial_threshold;
    }

    /**
     * @brief This method computes the uniaxial threshold using a hardening-softening law
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param TensileIndicatorFactor The tensile indicator
     * @param CompressionIndicatorFactor The compressive indicator
     * param rEquivalentStressThreshold The maximum uniaxial stress of the linear behaviour
     * @param rSlope The slope of the PlasticDiss-Threshold curve
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateEquivalentStressThresholdHardeningCurveInitialHardeningExponentialSoftening(
        const double PlasticDissipation,
        const double TensileIndicatorFactor,
        const double CompressionIndicatorFactor,
        double& rEquivalentStressThreshold,
        double& rSlope,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();

        double initial_threshold;
        GetInitialUniaxialThreshold(rValues, initial_threshold);
        const double ultimate_stress = r_material_properties[MAXIMUM_STRESS];              // sikpi
        const double max_stress_position = r_material_properties[MAXIMUM_STRESS_POSITION]; // cappi

        if (PlasticDissipation < 1.0) {
            const double ro = std::sqrt(1.0 - initial_threshold / ultimate_stress);
            double alpha = std::log((1.0 - (1.0 - ro) * (1.0 - ro)) / ((3.0 - ro) * (1.0 + ro) * max_stress_position));
            alpha = std::exp(alpha / (1.0 - max_stress_position));
            const double phi = std::pow((1.0 - ro), 2.0) + ((3.0 - ro) * (1.0 + ro) * PlasticDissipation * (std::pow(alpha, (1.0 - PlasticDissipation))));

            rEquivalentStressThreshold = ultimate_stress * (2.0 * std::sqrt(phi) - phi);
            rSlope = ultimate_stress * ((1.0 / std::sqrt(phi)) - 1.0) * (3.0 - ro) * (1.0 + ro) * (std::pow(alpha, (1.0 - PlasticDissipation))) *
                     (1.0 - std::log(alpha) * PlasticDissipation);
        } else {
            KRATOS_ERROR << "PlasticDissipation > 1.0 " << PlasticDissipation << std::endl;
        }
    }

    /**
     * @brief This method computes the uniaxial threshold using a perfect plasticity law
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param TensileIndicatorFactor The tensile indicator
     * @param CompressionIndicatorFactor The compressive indicator
     * param rEquivalentStressThreshold The maximum uniaxial stress of the linear behaviour
     * @param rSlope The slope of the PlasticDiss-Threshold curve
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateEquivalentStressThresholdHardeningCurvePerfectPlasticity(
        const double PlasticDissipation,
        const double TensileIndicatorFactor,
        const double CompressionIndicatorFactor,
        double& rEquivalentStressThreshold,
        double& rSlope,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        double initial_threshold;
        GetInitialUniaxialThreshold(rValues, initial_threshold);

        rEquivalentStressThreshold = initial_threshold;
        rSlope = -0.5 * initial_threshold;
    }

    /**
     * @brief This method computes the uniaxial threshold using a perfect plasticity law
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param TensileIndicatorFactor The tensile indicator
     * @param CompressionIndicatorFactor The compressive indicator
     * param rEquivalentStressThreshold The maximum uniaxial stress of the linear behaviour
     * @param rSlope The slope of the PlasticDiss-Threshold curve
     * @param rValues Parameters of the constitutive law
     * @param rPlasticStrain The Plastic Strain internal variable
     */
	static void CalculateEquivalentStressThresholdCurveFittingHardening(
        const double PlasticDissipation,
        const double TensileIndicatorFactor,
        const double CompressionIndicatorFactor,
        double& rEquivalentStressThreshold,
        double& rSlope,
        ConstitutiveLaw::Parameters& rValues,
        const double EquivalentPlasticStrain      
        )
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();
        const Vector& curve_fitting_parameters = r_material_properties[CURVE_FITTING_PARAMETERS];
        const Vector& plastic_strain_indicators = r_material_properties[PLASTIC_STRAIN_INDICATORS];
        const double fracture_energy = r_material_properties[FRACTURE_ENERGY];

        const SizeType order_polinomial = curve_fitting_parameters.size();
        const double plastic_strain_indicator_1 = plastic_strain_indicators[0];
        const double plastic_strain_indicator_2 = plastic_strain_indicators[1];

        // Compute the initial and the final stresses
        double stress_indicator_1 = curve_fitting_parameters[0];
        double dS_dEp_slope = 0.0;

        for (IndexType i = 1; i < order_polinomial; ++i) {
            stress_indicator_1 += curve_fitting_parameters[i] * (std::pow(plastic_strain_indicator_1, i));
            dS_dEp_slope += i * curve_fitting_parameters[i] * std::pow(plastic_strain_indicator_1, i - 1);
        }

        double dKp_dEp = stress_indicator_1 / fracture_energy;
        double dS_dEp = 0.0; // initial slope is zero

        const double stress_indicator_2 =  stress_indicator_1 +  dS_dEp * (plastic_strain_indicator_2 - plastic_strain_indicator_1);

        // Compute volumetric fracture energies of each region
        double Gt1 = 0.0;
        for (IndexType i = 0; i < order_polinomial; ++i) {
            Gt1 += curve_fitting_parameters[i] * (std::pow(plastic_strain_indicator_1, i + 1)) / (i + 1);
        }
        const double Gt2 = (stress_indicator_1 * stress_indicator_2) * (plastic_strain_indicator_2 - plastic_strain_indicator_1) * 0.5;
        const double Gt3 = fracture_energy - Gt2 - Gt1;

        KRATOS_ERROR_IF(Gt3 < 0.0) << "Fracture energy too low in CurveFittingHardening of plasticity..."  << std::endl;

        // Compute segment threshold
        const double segment_threshold = (Gt2 + Gt1) / fracture_energy;

        if (PlasticDissipation <= segment_threshold) {
            const double Eps = EquivalentPlasticStrain;

            if (EquivalentPlasticStrain < plastic_strain_indicator_1) { // Polinomial region
                double S_Ep = curve_fitting_parameters[0];
                double dS_dEp = 0.0;
                for (IndexType i = 1; i < order_polinomial; ++i) {
                    S_Ep += curve_fitting_parameters[i] * std::pow(Eps, i);
                    dS_dEp += i *  curve_fitting_parameters[i] * std::pow(Eps, i - 1);
                }
                dKp_dEp = S_Ep / fracture_energy;

                rEquivalentStressThreshold = S_Ep;
                rSlope = dS_dEp / dKp_dEp;
            } else { // Linear region
                const double alpha =  std::pow(stress_indicator_1, 2);
                const double beta = (std::pow(stress_indicator_2, 2) - alpha) / (plastic_strain_indicator_2 - plastic_strain_indicator_1);

                const double S_Ep = std::sqrt(alpha + beta * (Eps - plastic_strain_indicator_1));

                dS_dEp = 0.5 * (beta) / S_Ep;
                dKp_dEp = S_Ep / fracture_energy;
                
                rEquivalentStressThreshold = S_Ep;
                rSlope = dS_dEp / dKp_dEp;
            }
        } else { // Exponential softening
            const double Eps = EquivalentPlasticStrain;
            const double alpha = std::pow(stress_indicator_1, 2);
            const double beta = (std::pow(stress_indicator_2, 2) - alpha) / (plastic_strain_indicator_2 - plastic_strain_indicator_1);

            const double S_Ep = std::sqrt(alpha + beta * (Eps - plastic_strain_indicator_1));
            const double plastic_dissipation_region_3 = PlasticDissipation - segment_threshold;

            const double beta2 = 1.5 * S_Ep / Gt3;
            const double alpha2 = std::sqrt((plastic_dissipation_region_3 * 2.0 * beta2 * fracture_energy / S_Ep) + 1.0);
            rEquivalentStressThreshold = S_Ep * alpha2 * (2.0 - alpha2);
            rSlope = 2.0 * beta2 * fracture_energy * (1.0 / alpha2 - 1.0);
        }
    }

    /**
     * @brief This method returns the equivalent plastic strain
     * @param rThreshold The uniaxial stress threshold
     * @param rValues Parameters of the constitutive law
     * @param rStressVector The stress vector
     * @param r0 The tensile indicator
     * @param rEquivalentPlasticStrain The equivalent plastic strain
     * @param rPlasticStrain The plastic strain vector
     */
    static void CalculateEquivalentPlasticStrain(
        const Vector& rStressVector, 
        const double UniaxialStress,
        const Vector& rPlasticStrain,
        const double r0,
        ConstitutiveLaw::Parameters& rValues,
        double& rEquivalentPlasticStrain
        ) 
    {
        double scalar_product = 0.0;
        for (IndexType i = 0; i < rPlasticStrain.size(); ++i) {
            scalar_product += rStressVector[i] * rPlasticStrain[i];
        }

        /*Since this law is used for Von Mises and Tresca, no
        scaling is necessary, even though the needed params are available*/
        rEquivalentPlasticStrain = scalar_product / UniaxialStress;
    }

    /**
     * @brief This method returns the initial uniaxial stress threshold
     * @param rThreshold The uniaxial stress threshold
     * @param rValues Parameters of the constitutive law
     */
    static void GetInitialUniaxialThreshold(ConstitutiveLaw::Parameters& rValues, double& rThreshold)
    {
        TYieldSurfaceType::GetInitialUniaxialThreshold(rValues, rThreshold);
    }

    /**
     * @brief This method computes hardening parameter needed for the algorithm
     * @param rGFlux The derivative of the plastic potential
     * @param SlopeThreshold The slope of the PlasticDiss-Threshold curve
     * @param rHardeningParameter The hardening parameter needed for the algorithm
     * @param rSlope The slope of the PlasticDiss-Threshold curve
     */
    static void CalculateHardeningParameter(
        const array_1d<double, VoigtSize>& rGFlux,
        const double SlopeThreshold,
        const array_1d<double, VoigtSize>& rHCapa,
        double& rHardeningParameter
        )
    {
        rHardeningParameter = -SlopeThreshold;
        double aux = 0.0;

        for (IndexType i = 0; i < VoigtSize; ++i) {
            aux += rHCapa[i] * rGFlux[i];
        }
        if (aux != 0.0)
            rHardeningParameter *= aux;
    }

    /**
     * @brief This method computes the plastic denominator needed
     * to compute the plastic consistency factor
     * @param rFflux The derivative of the yield surface
     * @param rGflux The derivative of the plastic potential
     * @param rConstitutiveMatrix The elastic constitutive matrix
     * @param rHardeningParameter The hardening parameter needed for the algorithm
     * @param rPlasticDenominator The plasticity numerical value to obtain the pastic consistency factor
     */
    static void CalculatePlasticDenominator(
        const array_1d<double, VoigtSize>& rFFlux,
        const array_1d<double, VoigtSize>& rGFlux,
        const Matrix& rConstitutiveMatrix,
        double& rHardeningParameter,
        double& rPlasticDenominator,
        const Vector& rBackStressVector,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        const Vector& kinematic_parameters = rValues.GetMaterialProperties()[KINEMATIC_PLASTICITY_PARAMETERS];
        const int kinematic_hardening_type = rValues.GetMaterialProperties()[KINEMATIC_HARDENING_TYPE];

        const array_1d<double, VoigtSize> delta_vector = prod(rGFlux, rConstitutiveMatrix);
        double A1 = 0.0;
        for (IndexType i = 0; i < VoigtSize; ++i) {
            A1 += rFFlux[i] * delta_vector[i];
        }
        A1 *= (1.0 - kinematic_parameters[2]);

        double dot_fflux_gflux = 0.0, A2;
        for (IndexType i = 0; i < VoigtSize; ++i) {
            dot_fflux_gflux += rFFlux[i] * rGFlux[i];
        }
        const double two_thirds = 2.0 / 3.0;
		double dot_fflux_backstress = 0.0;
        switch (static_cast<KinematicHardeningType>(kinematic_hardening_type))
        {
            case KinematicHardeningType::LinearKinematicHardening:
                A2 = two_thirds * kinematic_parameters[0] * dot_fflux_gflux;
                break;

            case KinematicHardeningType::AmstrongFrederickKinematicHardening:
                A2 = two_thirds * kinematic_parameters[0] * dot_fflux_gflux;

                for (IndexType i = 0; i < VoigtSize; ++i) {
                    dot_fflux_backstress += rFFlux[i] * rBackStressVector[i];
                }
                A2 -= kinematic_parameters[1] * dot_fflux_backstress * std::sqrt(two_thirds * std::abs(dot_fflux_backstress));
                break;

            default:
                KRATOS_ERROR << " The Kinematic hardening type of plasticity is not set or wrong..." << kinematic_hardening_type << std::endl;
                break;
        }

        const double A3 = rHardeningParameter;
        rPlasticDenominator = 1.0 / (A1 + A2 + A3);
		rPlasticDenominator *= (1.0 - kinematic_parameters[2]);
    }

    /**
     * @brief This method defines in the CL integrator
     * @return 0 if OK, 1 otherwise
     */
    static int Check(const Properties& rMaterialProperties)
    {
        // Checking is registered
        KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS);
        KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS);
        KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS_TENSION);
        KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS_COMPRESSION);
        KRATOS_CHECK_VARIABLE_KEY(HARDENING_CURVE);
        KRATOS_CHECK_VARIABLE_KEY(MAXIMUM_STRESS);
        KRATOS_CHECK_VARIABLE_KEY(MAXIMUM_STRESS_POSITION);
        KRATOS_CHECK_VARIABLE_KEY(FRACTURE_ENERGY);

        // Checking is defined
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS)) << "HARDENING_CURVE is not a defined value" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(HARDENING_CURVE)) << "HARDENING_CURVE is not a defined value" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(FRACTURE_ENERGY)) << "FRACTURE_ENERGY is not a defined value" << std::endl;

        // Checking curves
        const int curve_type = rMaterialProperties[HARDENING_CURVE];
        if (static_cast<HardeningCurveType>(curve_type) == HardeningCurveType::InitialHardeningExponentialSoftening) {
            KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(MAXIMUM_STRESS)) << "MAXIMUM_STRESS is not a defined value" << std::endl;
            KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(MAXIMUM_STRESS_POSITION)) << "MAXIMUM_STRESS_POSITION is not a defined value" << std::endl;
        } else if (static_cast<HardeningCurveType>(curve_type) == HardeningCurveType::CurveFittingHardening) {
            KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(CURVE_FITTING_PARAMETERS)) << "CURVE_FITTING_PARAMETERS is not a defined value" << std::endl;
            KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(PLASTIC_STRAIN_INDICATORS)) << "PLASTIC_STRAIN_INDICATORS is not a defined value" << std::endl;
        }
        
        if (!rMaterialProperties.Has(YIELD_STRESS)) {
            KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YIELD_STRESS_TENSION)) << "YIELD_STRESS_TENSION is not a defined value" << std::endl;
            KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YIELD_STRESS_COMPRESSION)) << "YIELD_STRESS_COMPRESSION is not a defined value" << std::endl;
            
            const double yield_compression = rMaterialProperties[YIELD_STRESS_COMPRESSION];
            const double yield_tension = rMaterialProperties[YIELD_STRESS_TENSION];

            KRATOS_ERROR_IF(yield_compression < tolerance) << "Yield stress in compression almost zero or negative, include YIELD_STRESS_COMPRESSION in definition";
            KRATOS_ERROR_IF(yield_tension < tolerance) << "Yield stress in tension almost zero or negative, include YIELD_STRESS_TENSION in definition";
        } else {
            const double yield_stress = rMaterialProperties[YIELD_STRESS];

            KRATOS_ERROR_IF(yield_stress < tolerance) << "Yield stress almost zero or negative, include YIELD_STRESS in definition";
        }

        return TYieldSurfaceType::Check(rMaterialProperties);
    }

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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class GenericYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
#endif
