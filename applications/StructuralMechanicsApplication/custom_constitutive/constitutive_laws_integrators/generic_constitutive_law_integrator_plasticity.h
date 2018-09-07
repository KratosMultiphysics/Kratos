// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Lucia Barbu
//

#if !defined(KRATOS_GENERIC_CONSTITUTIVE_LAW_INTEGRATOR_PLASTICITY_H_INCLUDED)
#define KRATOS_GENERIC_CONSTITUTIVE_LAW_INTEGRATOR_PLASTICITY_H_INCLUDED

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
 * @class GenericConstitutiveLawIntegratorPlasticity
 * @ingroup StructuralMechanicsApplication
 * @brief: This object integrates the predictive stress using the plasticity theory by means of
 * linear/exponential softening or hardening+softening evolution laws
 * @details
 * @tparam TYieldSurfaceType
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <class TYieldSurfaceType>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) GenericConstitutiveLawIntegratorPlasticity
{
  public:
    ///@name Type Definitions
    ///@{
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    /// Definition of index
    typedef std::size_t IndexType;

    /// Definition of size type
    typedef std::size_t SizeType;

    /// The type of yield surface
    typedef TYieldSurfaceType YieldSurfaceType;

    /// The type of plastic potential
    typedef typename YieldSurfaceType::PlasticPotentialType PlasticPotentialType;

    /// Counted pointer of GenericConstitutiveLawIntegratorPlasticity
    KRATOS_CLASS_POINTER_DEFINITION(GenericConstitutiveLawIntegratorPlasticity);

    ///@}
    ///@name  Enum's
    ///@{

    enum class HardeningCurveType
    {
        LinearSoftening = 0,
        ExponentialSoftening = 1,
        InitialHardeningExponentialSoftening = 2,
        PerfectPlasticity = 3
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor
    GenericConstitutiveLawIntegratorPlasticity()
    {
    }

    /// Copy constructor
    GenericConstitutiveLawIntegratorPlasticity(GenericConstitutiveLawIntegratorPlasticity const &rOther)
    {
    }

    /// Assignment operator
    GenericConstitutiveLawIntegratorPlasticity &operator=(GenericConstitutiveLawIntegratorPlasticity const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~GenericConstitutiveLawIntegratorPlasticity()
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
     * @param rC The elastic constitutive matrix
     * @param rPlasticStrain The elastic constitutive matrix
     * @param rMaterialProperties The material properties
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void IntegrateStressVector(
        Vector& rPredictiveStressVector,
        Vector& rStrainVector,
        double& rUniaxialStress,
        double& rThreshold,
        double& rPlasticDenominator,
        Vector& rFflux,
        Vector& rGflux,
        double& rPlasticDissipation,
        Vector& rPlasticStrainIncrement,
        const Matrix& rC,
        Vector& rPlasticStrain,
        const Properties& rMaterialProperties,
        const double CharacteristicLength
        )
    {
        bool is_converged = false;
        IndexType iteration = 0, max_iter = 9000;
        BoundedVector<double, 6> delta_sigma;
        double plastic_consistency_factor_increment, F;

        // Backward Euler
        while (is_converged == false && iteration <= max_iter) {
            F = rUniaxialStress - rThreshold;
            plastic_consistency_factor_increment = F * rPlasticDenominator;
            //if (plastic_consistency_factor_increment < 0.0) plastic_consistency_factor_increment = 0.0;
            noalias(rPlasticStrainIncrement) = plastic_consistency_factor_increment * rGflux;
            noalias(rPlasticStrain) += rPlasticStrainIncrement;
            noalias(delta_sigma) = prod(rC, rPlasticStrainIncrement);
            noalias(rPredictiveStressVector) -= delta_sigma;

            CalculatePlasticParameters(rPredictiveStressVector, rStrainVector, rUniaxialStress, rThreshold,
                                       rPlasticDenominator, rFflux, rGflux, rPlasticDissipation, rPlasticStrainIncrement,
                                       rC, rMaterialProperties, CharacteristicLength);

            F = rUniaxialStress - rThreshold;

            if (std::abs(F) <= std::abs(1.0e-4 * rThreshold)) { // Has converged
                is_converged = true;
            } else {
                iteration++;
            }
        }
        KRATOS_WARNING_IF("Backward Euler Plasticity", iteration == max_iter) << "Maximum number of iterations in plasticity loop reached..." << std::endl;
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
     * @param rC The elastic constitutive matrix
     * @param rPlasticStrain The elastic constitutive matrix
     * @param rMaterialProperties The material properties
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculatePlasticParameters(
        Vector& rPredictiveStressVector,
        Vector& rStrainVector,
        double& rUniaxialStress,
        double& rThreshold,
        double& rPlasticDenominator,
        Vector& rFflux,
        Vector& rGflux,
        double& rPlasticDissipation,
        Vector& rPlasticStrainIncrement,
        const Matrix& rC,
        const Properties& rMaterialProperties,
        const double CharacteristicLength
        )
    {
        Vector deviator = ZeroVector(6);
        Vector h_capa = ZeroVector(6);
        double J2, r0, r1, slope, hardening_parameter;

        YieldSurfaceType::CalculateEquivalentStress(rPredictiveStressVector, rStrainVector, rUniaxialStress, rMaterialProperties);
        const double I1 = rPredictiveStressVector[0] + rPredictiveStressVector[1] + rPredictiveStressVector[2];
        ConstitutiveLawUtilities::CalculateJ2Invariant(rPredictiveStressVector, I1, deviator, J2);
        CalculateFFluxVector(rPredictiveStressVector, deviator, J2, rFflux, rMaterialProperties);
        CalculateGFluxVector(rPredictiveStressVector, deviator, J2, rGflux, rMaterialProperties);
        CalculateRFactors(rPredictiveStressVector, r0, r1);
        CalculatePlasticDissipation(rPredictiveStressVector, r0, r1, rPlasticStrainIncrement, rPlasticDissipation, h_capa, rMaterialProperties, CharacteristicLength);
        CalculateEquivalentStressThreshold(rPlasticDissipation, r0, r1, rThreshold, slope, rMaterialProperties);
        CalculateHardeningParameter(rFflux, slope, h_capa, hardening_parameter);
        CalculatePlasticDenominator(rFflux, rGflux, rC, hardening_parameter, rPlasticDenominator);
    }

    /**
     * @brief This method calculates the derivative of the yield surface
     * @param StressVector The stress vector
     * @param Deviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the deviatoric part of the stress vector
     * @param FFluxVector The derivative of the yield surface
     * @param rMaterialProperties The material properties
     */
    static void CalculateFFluxVector(
        const Vector& StressVector,
        const Vector& Deviator,
        const double J2,
        Vector& FFluxVector,
        const Properties& rMaterialProperties
        )
    {
        YieldSurfaceType::CalculateYieldSurfaceDerivative(StressVector, Deviator, J2,
                                                          FFluxVector, rMaterialProperties);
    }

    /**
     * @brief This method calculates the derivative of the plastic potential
     * @param StressVector The stress vector
     * @param Deviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the deviatoric part of the stress vector
     * @param GFluxVector The derivative of the yield surface
     * @param rMaterialProperties The material properties
     */
    static void CalculateGFluxVector(
        const Vector& StressVector,
        const Vector& Deviator,
        const double J2,
        Vector& GFluxVector,
        const Properties& rMaterialProperties
        )
    {
        YieldSurfaceType::CalculatePlasticPotentialDerivative(StressVector, Deviator, J2,
                                                              GFluxVector, rMaterialProperties);
    }

    /**
     * @brief This method computes the tensile/compressive indicators
     * @param rStressVector The stress vector
     * @param r0 The tensile indicator
     * @param r1 The compressive indicator
     */
    static void CalculateRFactors(
        const Vector& rStressVector,
        double& r0,
        double& r1
        )
    {
        // We do an initial check
        if (norm_2(rStressVector) < 1.0e-8) {
            r0 = 0.5;
            r1 = 0.5;
            return;
        }

        // We proceed as usual
        Vector principal_stresses = ZeroVector(3);
        ConstitutiveLawUtilities::CalculatePrincipalStresses(principal_stresses, rStressVector);

        double suma = 0.0, sumb = 0.0, sumc = 0.0;
        double aux_sa;

        for (IndexType i = 0; i < 3; ++i) {
            aux_sa = std::abs(principal_stresses[i]);
            suma += aux_sa;
            sumb += 0.5 * (principal_stresses[i] + aux_sa);
            sumc += 0.5 * (-principal_stresses[i] + aux_sa);
        }

        if (std::abs(suma) > tolerance) {
            r0 = sumb / suma;
            r1 = sumc / suma;
        } else {
            r0 = sumb;
            r1 = sumc;
        }

        // Final check
        if ((std::abs(r0) + std::abs(r1)) < tolerance) {
            r0 = 0.5;
            r1 = 0.5;
            return;
        }
    }

    /**
     * @brief This method computes the plastic dissipation of the plasticity model
     * @param StressVector The stress vector
     * @param r0 The tensile indicator
     * @param r1 The compressive indicator
     * @param PlasticStrainIncrement The increment of plastic strain of this time step
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param rHCapa The slope of the PlasticDiss-Threshold curve
     * @param rMaterialProperties The material properties
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculatePlasticDissipation(
        const Vector& StressVector,
        const double r0,
        const double r1,
        const Vector& PlasticStrainInc,
        double& rPlasticDissipation,
        Vector& rHCapa,
        const Properties& rMaterialProperties,
        const double CharacteristicLength
        )
    {
        const double Young = rMaterialProperties[YOUNG_MODULUS];
        const double yield_compression = rMaterialProperties[YIELD_STRESS_COMPRESSION];
        const double yield_tension = rMaterialProperties[YIELD_STRESS_TENSION];
        const double n = yield_compression / yield_tension;
        const double Gf = rMaterialProperties[FRACTURE_ENERGY];                   // Frac energy in tension
        const double Gfc = rMaterialProperties[FRACTURE_ENERGY] * std::pow(n, 2); // Frac energy in compression

        const double gf = Gf / CharacteristicLength;
        const double gfc = Gfc / CharacteristicLength;

        const double hlim = 2.0 * Young * gfc / (std::pow(yield_compression, 2));
        KRATOS_ERROR_IF(CharacteristicLength > hlim) << "The Fracture Energy is to low: " << gfc << std::endl;

        double constant0 = 0.0, constant1 = 0.0, dplastic_dissipation = 0.0;
        if (gf > 0.000001) {
            constant0 = r0 / gf;
            constant1 = r1 / gfc;
        }
        const double constant = constant0 + constant1;

        for (IndexType i = 0; i < 6; ++i) {
            rHCapa[i] = constant * StressVector[i];
            dplastic_dissipation += rHCapa[i] * PlasticStrainInc[i];
        }

        if (dplastic_dissipation<0.0 || dplastic_dissipation> 1.0)
            dplastic_dissipation = 0.0;

        rPlasticDissipation += dplastic_dissipation;
        if (rPlasticDissipation >= 1.0)
            rPlasticDissipation = 0.9999;
    }

    /**
     * @brief This method computes the uniaxial threshold that differentiates the elastic-plastic behaviour
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param r0 The tensile indicator
     * @param r1 The compressive indicator
     * param rEquivalentStressThreshold The maximum uniaxial stress of the linear behaviour
     * @param rSlope The slope of the PlasticDiss-Threshold curve
     * @param rMaterialProperties The material properties
     */
    static void CalculateEquivalentStressThreshold(
        const double PlasticDissipation,
        const double r0,
        const double r1,
        double& rEquivalentStressThreshold,
        double& rSlope,
        const Properties& rMaterialProperties
        )
    {
        const int curve_type = rMaterialProperties[HARDENING_CURVE];
        const double yield_comp = rMaterialProperties[YIELD_STRESS_COMPRESSION];
        const double yield_tension = rMaterialProperties[YIELD_STRESS_TENSION];
        const double n = yield_comp / yield_tension;

        BoundedVector<double, 2> Gf, slopes, eq_thresholds;

        Gf[0] = rMaterialProperties[FRACTURE_ENERGY];
        Gf[1] = std::pow(n, 2) * Gf[0];

        for (IndexType i = 0; i < 2; ++i) { // i:0 Tension ; i:1 compression
            switch (static_cast<HardeningCurveType>(curve_type))
            {
            case HardeningCurveType::LinearSoftening:
                CalculateEquivalentStressThresholdHardeningCurveLinearSoftening(PlasticDissipation, r0, r1, eq_thresholds[i], slopes[i], rMaterialProperties);
                break;

            case HardeningCurveType::ExponentialSoftening:
                CalculateEquivalentStressThresholdHardeningCurveExponentialSoftening(PlasticDissipation, r0, r1, eq_thresholds[i], slopes[i], rMaterialProperties);
                break;

            case HardeningCurveType::InitialHardeningExponentialSoftening:
                CalculateEquivalentStressThresholdHardeningCurveInitialHardeningExponentialSoftening(PlasticDissipation, r0, r1, eq_thresholds[i], slopes[i], rMaterialProperties);
                break;

            case HardeningCurveType::PerfectPlasticity:
                CalculateEquivalentStressThresholdHardeningCurvePerfectPlasticity(PlasticDissipation, r0, r1, eq_thresholds[i], slopes[i], rMaterialProperties);
                break;

                // Add more cases...

            default:
                KRATOS_ERROR << " The HARDENING_CURVE of plasticity is not set or wrong..." << curve_type << std::endl;
                break;
            }
        }
        rEquivalentStressThreshold = r0 * eq_thresholds[0] + r1 * eq_thresholds[1];
        rSlope = rEquivalentStressThreshold * ((r0 * slopes[0] / eq_thresholds[0]) + (r1 * slopes[1] / eq_thresholds[1]));
        KRATOS_DEBUG_ERROR_IF(rEquivalentStressThreshold < tolerance) << "Threshold set to zero. r0: " << r0 << " eq_thresholds[0]: " << eq_thresholds[0] << " r1: " << r1 << " eq_thresholds[1]:" << eq_thresholds[1] << std::endl;
    }

    /**
     * @brief This method computes the uniaxial threshold using a linear softening
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param r0 The tensile indicator
     * @param r1 The compressive indicator
     * param rEquivalentStressThreshold The maximum uniaxial stress of the linear behaviour
     * @param rSlope The slope of the PlasticDiss-Threshold curve
     * @param rMaterialProperties The material properties
     */
    static void CalculateEquivalentStressThresholdHardeningCurveLinearSoftening(
        const double PlasticDissipation,
        const double r0,
        const double r1,
        double& rEquivalentStressThreshold,
        double& rSlope,
        const Properties& rMaterialProperties
        )
    {
        double initial_threshold;
        GetInitialUniaxialThreshold(rMaterialProperties, initial_threshold);

        rEquivalentStressThreshold = initial_threshold * std::sqrt(1.0 - PlasticDissipation);
        rSlope = -0.5 * (std::pow(initial_threshold, 2.0) / (rEquivalentStressThreshold));
    }

    /**
     * @brief This method computes the uniaxial threshold using a exponential softening
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param r0 The tensile indicator
     * @param r1 The compressive indicator
     * param rEquivalentStressThreshold The maximum uniaxial stress of the linear behaviour
     * @param rSlope The slope of the PlasticDiss-Threshold curve
     * @param rMaterialProperties The material properties
     */
    static void CalculateEquivalentStressThresholdHardeningCurveExponentialSoftening(
        const double PlasticDissipation,
        const double r0,
        const double r1,
        double& rEquivalentStressThreshold,
        double& rSlope,
        const Properties& rMaterialProperties)
    {
        double initial_threshold;
        GetInitialUniaxialThreshold(rMaterialProperties, initial_threshold);

        rEquivalentStressThreshold = initial_threshold * (1.0 - PlasticDissipation);
        rSlope = -0.5 * initial_threshold;
    }

    /**
     * @brief This method computes the uniaxial threshold using a hardening-softening law
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param r0 The tensile indicator
     * @param r1 The compressive indicator
     * param rEquivalentStressThreshold The maximum uniaxial stress of the linear behaviour
     * @param rSlope The slope of the PlasticDiss-Threshold curve
     * @param rMaterialProperties The material properties
     */
    static void CalculateEquivalentStressThresholdHardeningCurveInitialHardeningExponentialSoftening(
        const double PlasticDissipation,
        const double r0,
        const double r1,
        double& rEquivalentStressThreshold,
        double& rSlope,
        const Properties& rMaterialProperties
        )
    {
        double initial_threshold;
        GetInitialUniaxialThreshold(rMaterialProperties, initial_threshold);
        const double ultimate_stress = rMaterialProperties[MAXIMUM_STRESS];              // sikpi
        const double max_stress_position = rMaterialProperties[MAXIMUM_STRESS_POSITION]; // cappi

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
     * @param r0 The tensile indicator
     * @param r1 The compressive indicator
     * param rEquivalentStressThreshold The maximum uniaxial stress of the linear behaviour
     * @param rSlope The slope of the PlasticDiss-Threshold curve
     * @param rMaterialProperties The material properties
     */
    static void CalculateEquivalentStressThresholdHardeningCurvePerfectPlasticity(
        const double PlasticDissipation,
        const double r0,
        const double r1,
        double& rEquivalentStressThreshold,
        double& rSlope,
        const Properties& rMaterialProperties
        )
    {
        double initial_threshold;
        GetInitialUniaxialThreshold(rMaterialProperties, initial_threshold);

        rEquivalentStressThreshold = initial_threshold;
        rSlope = -0.5 * initial_threshold;
    }

    /**
     * @brief This method returns the initial uniaxial stress threshold
     * @param rThreshold The uniaxial stress threshold
     * @param rMaterialProperties The material properties
     */
    static void GetInitialUniaxialThreshold(const Properties& rMaterialProperties, double& rThreshold)
    {
        TYieldSurfaceType::GetInitialUniaxialThreshold(rMaterialProperties, rThreshold);
    }

    /**
     * @brief This method computes hardening parameter needed for the algorithm
     * @param Gflux The derivative of the plastic potential
     * @param SlopeThreshold The slope of the PlasticDiss-Threshold curve
     * @param rHardeningParameter The hardening parameter needed for the algorithm
     * @param rSlope The slope of the PlasticDiss-Threshold curve
     */
    static void CalculateHardeningParameter(
        const Vector& GFlux,
        const double SlopeThreshold,
        const Vector& HCapa,
        double& rHardeningParameter
        )
    {
        rHardeningParameter = -SlopeThreshold;
        double aux = 0.0;

        for (IndexType i = 0; i < 6; ++i) {
            aux += HCapa[i] * GFlux[i];
        }
        if (aux != 0.0)
            rHardeningParameter *= aux;
    }

    /**
     * @brief This method computes the plastic denominator needed
     * to compute the plastic consistency factor
     * @param rFflux The derivative of the yield surface
     * @param rGflux The derivative of the plastic potential
     * @param rC The elastic constitutive matrix
     * @param rHardeningParameter The hardening parameter needed for the algorithm
     * @param rPlasticDenominator The plasticity numerical value to obtain the pastic consistency factor
     */
    static void CalculatePlasticDenominator(
        const Vector& rFFlux,
        const Vector& rGFlux,
        const Matrix& rC,
        double& rHardeningParameter,
        double& rPlasticDenominator
        )
    {
        //const Vector delta_vector = prod(C, GFlux);
        const Vector delta_vector = prod(rGFlux, rC);
        double A1 = 0.0;

        for (IndexType i = 0; i < 6; ++i) {
            A1 += rFFlux[i] * delta_vector[i];
        }

        const double A2 = 0.0; // Only for isotropic hard
        const double A3 = rHardeningParameter;
        rPlasticDenominator = 1.0 / (A1 + A2 + A3);
    }

    /**
     * @brief This method defines in the CL integrator
     * @return 0 if OK, 1 otherwise
     */
    static int Check(const Properties& rMaterialProperties)
    {
        KRATOS_CHECK_VARIABLE_KEY(HARDENING_CURVE);
        KRATOS_CHECK_VARIABLE_KEY(MAXIMUM_STRESS);
        KRATOS_CHECK_VARIABLE_KEY(MAXIMUM_STRESS_POSITION);
        KRATOS_CHECK_VARIABLE_KEY(FRACTURE_ENERGY);

        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(HARDENING_CURVE)) << "HARDENING_CURVE is not a defined value" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(MAXIMUM_STRESS)) << "MAXIMUM_STRESS is not a defined value" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(MAXIMUM_STRESS_POSITION)) << "MAXIMUM_STRESS_POSITION is not a defined value" << std::endl;
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(FRACTURE_ENERGY)) << "FRACTURE_ENERGY is not a defined value" << std::endl;

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

    // Serialization

    friend class Serializer;

    void save(Serializer &rSerializer) const
    {
    }

    void load(Serializer &rSerializer)
    {
    }

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
