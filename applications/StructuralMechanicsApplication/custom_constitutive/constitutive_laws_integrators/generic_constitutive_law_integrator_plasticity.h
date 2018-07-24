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
     * @param PredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param StrainVector The equivalent strain vector of that integration point
     * @param UniaxialStress The equivalent uniaxial stress 
     * @param Threshold The maximum uniaxial stress of the linear behaviour
     * @param PlasticDenominator The plasticity numerical value to obtain the pastic consistency factor
     * @param Fflux The derivative of the yield surface
     * @param Gflux The derivative of the plastic potential
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param PlasticStrainIncrement The increment of plastic strain of this time step
     * @param C The elastic constitutive matrix
     * @param PlasticStrain The elastic constitutive matrix
     * @param rMaterialProperties The material properties
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void IntegrateStressVector(
        Vector& rPredictiveStressVector,
        Vector& StrainVector,
        double& UniaxialStress,
        double& Threshold,
        double& PlasticDenominator,
        Vector& Fflux,
        Vector& Gflux,
        double& PlasticDissipation,
        Vector& PlasticStrainIncrement,
        const Matrix& C,
        Vector& PlasticStrain,
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
            F = UniaxialStress - Threshold;
            plastic_consistency_factor_increment = F * PlasticDenominator;
            //if (plastic_consistency_factor_increment < 0.0) plastic_consistency_factor_increment = 0.0;
            noalias(PlasticStrainIncrement) = plastic_consistency_factor_increment * Gflux;
            noalias(PlasticStrain) += PlasticStrainIncrement;
            noalias(delta_sigma) = prod(C, PlasticStrainIncrement);
            noalias(rPredictiveStressVector) -= delta_sigma;

            CalculatePlasticParameters(rPredictiveStressVector, StrainVector, UniaxialStress, Threshold,
                                       PlasticDenominator, Fflux, Gflux, PlasticDissipation, PlasticStrainIncrement,
                                       C, rMaterialProperties, CharacteristicLength);

            F = UniaxialStress - Threshold;

            if (std::abs(F) <= std::abs(1.0e-4 * Threshold)) { // Has converged
                is_converged = true;
            } else {
                iteration++;
            }
        }
        KRATOS_WARNING_IF("Backward Euler Plasticity", iteration == max_iter) << "Maximum number of iterations in plasticity loop reached..." << std::endl;
    }

    /**
     * @brief This method calculates all the plastic parameters required for the integration of the PredictiveStressVector
     * @param PredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param StrainVector The equivalent strain vector of that integration point
     * @param UniaxialStress The equivalent uniaxial stress 
     * @param Threshold The maximum uniaxial stress of the linear behaviour
     * @param PlasticDenominator The plasticity numerical value to obtain the pastic consistency factor
     * @param Fflux The derivative of the yield surface
     * @param Gflux The derivative of the plastic potential
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param PlasticStrainIncrement The increment of plastic strain of this time step
     * @param C The elastic constitutive matrix
     * @param PlasticStrain The elastic constitutive matrix
     * @param rMaterialProperties The material properties
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculatePlasticParameters(
        Vector& PredictiveStressVector,
        Vector& StrainVector,
        double& UniaxialStress,
        double& Threshold,
        double& PlasticDenominator,
        Vector& Fflux,
        Vector& Gflux,
        double& PlasticDissipation,
        Vector& PlasticStrainIncrement,
        const Matrix& C,
        const Properties& rMaterialProperties,
        const double CharacteristicLength
        )
    {
        Vector deviator = ZeroVector(6);
        Vector h_capa = ZeroVector(6);
        double J2, r0, r1, Slope, HardParam;

        YieldSurfaceType::CalculateEquivalentStress(PredictiveStressVector, StrainVector,
                                                    UniaxialStress, rMaterialProperties);
        const double I1 = PredictiveStressVector[0] + PredictiveStressVector[1] + PredictiveStressVector[2];
        ConstitutiveLawUtilities::CalculateJ2Invariant(PredictiveStressVector, I1, deviator, J2);
        CalculateFFluxVector(PredictiveStressVector, deviator, J2, Fflux, rMaterialProperties);
        CalculateGFluxVector(PredictiveStressVector, deviator, J2, Gflux, rMaterialProperties);
        CalculateRFactors(PredictiveStressVector, r0, r1);
        CalculatePlasticDissipation(PredictiveStressVector, r0, r1, PlasticStrainIncrement,
                                    PlasticDissipation, h_capa, rMaterialProperties, CharacteristicLength);
        CalculateEquivalentStressThreshold(PlasticDissipation, r0,
                                           r1, Threshold, Slope, rMaterialProperties);
        CalculateHardeningParameter(Fflux, Slope, h_capa, HardParam);
        CalculatePlasticDenominator(Fflux, Gflux, C, HardParam, PlasticDenominator);
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
     * @param StressVector The stress vector 
     * @param r0 The tensile indicator
     * @param r1 The compressive indicator
     */
    static void CalculateRFactors(
        const Vector& StressVector,
        double& r0,
        double& r1
        )
    {
        Vector principal_stresses = ZeroVector(3);
        ConstitutiveLawUtilities::CalculatePrincipalStresses(principal_stresses, StressVector);

        double suma = 0.0, sumb = 0.0, sumc = 0.0;
        Vector SA = ZeroVector(3), SB = ZeroVector(3), SC = ZeroVector(3);

        for (IndexType i = 0; i < 3; i++) {
            SA[i] = std::abs(principal_stresses[i]);
            SB[i] = 0.5 * (principal_stresses[i] + SA[i]);
            SC[i] = 0.5 * (-principal_stresses[i] + SA[i]);

            suma += SA[i];
            sumb += SB[i];
            sumc += SC[i];
        }

        if (std::abs(suma) > tolerance) {
            r0 = sumb / suma;
            r1 = sumc / suma;
        } else {
            r0 = sumb;
            r1 = sumc;
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

        for (IndexType i = 0; i < 6; i++) {
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

        for (IndexType i = 0; i < 2; i++) { // i:0 Tension ; i:1 compression
            switch (static_cast<HardeningCurveType>(curve_type))
            {
            case HardeningCurveType::LinearSoftening:
                CalculateEqStressThresholdHardCurve1(PlasticDissipation, r0, r1,
                                                     eq_thresholds[i], slopes[i], rMaterialProperties);
                break;

            case HardeningCurveType::ExponentialSoftening:
                CalculateEqStressThresholdHardCurve2(PlasticDissipation, r0, r1,
                                                     eq_thresholds[i], slopes[i], rMaterialProperties);
                break;

            case HardeningCurveType::InitialHardeningExponentialSoftening:
                CalculateEqStressThresholdHardCurve3(PlasticDissipation, r0, r1,
                                                     eq_thresholds[i], slopes[i], rMaterialProperties);
                break;

            case HardeningCurveType::PerfectPlasticity:
                CalculateEqStressThresholdHardCurve4(PlasticDissipation, r0, r1,
                                                     eq_thresholds[i], slopes[i], rMaterialProperties);
                break;

                // Add more cases...

            default:
                KRATOS_ERROR << " The HARDENING_CURVE of plasticity is not set or wrong..." << curve_type << std::endl;
                break;
            }
        }
        rEquivalentStressThreshold = r0 * eq_thresholds[0] + r1 * eq_thresholds[1];
        rSlope = rEquivalentStressThreshold * ((r0 * slopes[0] / eq_thresholds[0]) + (r1 * slopes[1] / eq_thresholds[1]));
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
    static void CalculateEqStressThresholdHardCurve1(
        const double PlasticDissipation,
        const double r0,
        const double r1,
        double& rEquivalentStressThreshold,
        double& rSlope,
        const Properties& rMaterialProperties
        )
    {
        //const double initial_threshold = rMaterialProperties[YIELD_STRESS_COMPRESSION];
        double initial_threshold;
        TYieldSurfaceType::GetInitialUniaxialThreshold(rMaterialProperties, initial_threshold);

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
    static void CalculateEqStressThresholdHardCurve2(
        const double PlasticDissipation,
        const double r0,
        const double r1,
        double& rEquivalentStressThreshold,
        double& rSlope,
        const Properties& rMaterialProperties)
    {
//         const double initial_threshold = rMaterialProperties[YIELD_STRESS_COMPRESSION];
        double initial_threshold;
        TYieldSurfaceType::GetInitialUniaxialThreshold(rMaterialProperties, initial_threshold);

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
    static void CalculateEqStressThresholdHardCurve3(
        const double PlasticDissipation,
        const double r0,
        const double r1,
        double& rEquivalentStressThreshold,
        double& rSlope,
        const Properties& rMaterialProperties
        )
    {
//         const double initial_threshold = rMaterialProperties[YIELD_STRESS_COMPRESSION];  // sikma
        double initial_threshold;
        TYieldSurfaceType::GetInitialUniaxialThreshold(rMaterialProperties, initial_threshold);
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
    static void CalculateEqStressThresholdHardCurve4(
        const double PlasticDissipation,
        const double r0,
        const double r1,
        double& rEquivalentStressThreshold,
        double& rSlope,
        const Properties& rMaterialProperties
        )
    {
        double initial_threshold;
        TYieldSurfaceType::GetInitialUniaxialThreshold(rMaterialProperties, initial_threshold);

        rEquivalentStressThreshold = initial_threshold;
        rSlope = -0.5 * initial_threshold;
    }
    /**
     * @brief This method computes hardening parameter needed for the algorithm
     * @param Gflux The derivative of the plastic potential
     * @param SlopeThreshold The slope of the PlasticDiss-Threshold curve
     * @param rHardParameter The hardening parameter needed for the algorithm
     * @param rSlope The slope of the PlasticDiss-Threshold curve
     */
    static void CalculateHardeningParameter(
        const Vector& GFlux,
        const double SlopeThreshold,
        const Vector& HCapa,
        double& rHardParameter
        )
    {
        rHardParameter = -SlopeThreshold;
        double aux = 0.0;

        for (IndexType i = 0; i < 6; i++) {
            aux += HCapa[i] * GFlux[i];
        }
        if (aux != 0.0)
            rHardParameter *= aux;
    }

    /**
     * @brief This method computes the plastic denominator needed 
     * to compute the plastic consistency factor
     * @param Fflux The derivative of the yield surface
     * @param Gflux The derivative of the plastic potential
     * @param C The elastic constitutive matrix
     * @param rHardParameter The hardening parameter needed for the algorithm
     * @param PlasticDenominator The plasticity numerical value to obtain the pastic consistency factor
     */
    static void CalculatePlasticDenominator(
        const Vector& FFlux,
        const Vector& GFlux,
        const Matrix& C,
        double& rHardParameter,
        double& PlasticDenominator
        )
    {
        //const Vector delta_vector = prod(C, GFlux);
        const Vector delta_vector = prod(GFlux, C);
        double A1 = 0.0;

        for (IndexType i = 0; i < 6; i++) {
            A1 += FFlux[i] * delta_vector[i];
        }

        const double A2 = 0.0; // Only for isotropic hard
        const double A3 = rHardParameter;
        PlasticDenominator = 1.0 / (A1 + A2 + A3);
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
