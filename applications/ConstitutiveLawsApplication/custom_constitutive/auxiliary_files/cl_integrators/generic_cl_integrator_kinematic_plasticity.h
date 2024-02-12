// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#pragma once

// System includes

// Project includes
#include "includes/define.h"
#include "includes/checks.h"
#include "includes/properties.h"
#include "utilities/math_utils.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "constitutive_laws_application_variables.h"
#include "generic_cl_integrator_plasticity.h"

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

    /// The definition of the Voigt array type
    typedef array_1d<double, VoigtSize> BoundedArrayType;

    /// The definition of the bounded matrix type
    typedef BoundedMatrix<double, Dimension, Dimension> BoundedMatrixType;

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
        ArmstrongFrederickKinematicHardening = 1,
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
     * @param rYieldSurfaceDerivative The derivative of the yield surface
     * @param rDerivativePlasticPotential The derivative of the plastic potential
     * @param rPlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param rPlasticStrainIncrement The increment of plastic strain of this time step
     * @param rConstitutiveMatrix The elastic constitutive matrix
     * @param rPlasticStrain The elastic constitutive matrix
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     * @param rBackStressVector The so-called BackStressVector used for the kinematic hardening
     * @param rPreviousStressVector The previous converged stress vector
     */
    static void IntegrateStressVector(
        BoundedArrayType& rPredictiveStressVector,
        Vector& rStrainVector,
        double& rUniaxialStress,
        double& rThreshold,
        double& rPlasticDenominator,
        BoundedArrayType& rYieldSurfaceDerivative,
        BoundedArrayType& rDerivativePlasticPotential,
        double& rPlasticDissipation,
        BoundedArrayType& rPlasticStrainIncrement,
        Matrix& rConstitutiveMatrix,
        Vector& rPlasticStrain,
        ConstitutiveLaw::Parameters& rValues,
        const double CharacteristicLength,
        Vector& rBackStressVector,
        const Vector& rPreviousStressVector
        )
    {
        // Material properties
        const Properties& r_material_properties = rValues.GetMaterialProperties();

        // Defining some variables
        bool is_converged = false;
        IndexType iteration = 0, max_iter = r_material_properties.Has(MAX_NUMBER_NL_CL_ITERATIONS) ? r_material_properties.GetValue(MAX_NUMBER_NL_CL_ITERATIONS) : 100;
        BoundedArrayType delta_sigma;
        double plastic_consistency_factor_increment, threshold_indicator;
        BoundedArrayType kin_hard_stress_vector;
        const bool analytic_tangent = (r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) && r_material_properties[TANGENT_OPERATOR_ESTIMATION] == 0) ? true : false;

        // Backward Euler
        while (!is_converged && iteration <= max_iter) {
            threshold_indicator = rUniaxialStress - rThreshold;
            plastic_consistency_factor_increment = threshold_indicator * rPlasticDenominator;
            plastic_consistency_factor_increment = (plastic_consistency_factor_increment < 0.0) ? 0.0 : plastic_consistency_factor_increment;
            noalias(rPlasticStrainIncrement) = plastic_consistency_factor_increment * rDerivativePlasticPotential;
            noalias(rPlasticStrain) += rPlasticStrainIncrement;
            noalias(delta_sigma) = prod(rConstitutiveMatrix, rPlasticStrainIncrement);
            noalias(rPredictiveStressVector) -= delta_sigma;
            CalculateBackStress(rPredictiveStressVector, rValues, rPreviousStressVector, rPlasticStrainIncrement, rBackStressVector);
            noalias(kin_hard_stress_vector) = rPredictiveStressVector - rBackStressVector;
            threshold_indicator = CalculatePlasticParameters(kin_hard_stress_vector, rStrainVector, rUniaxialStress, rThreshold, rPlasticDenominator, rYieldSurfaceDerivative, rDerivativePlasticPotential, rPlasticDissipation, rPlasticStrainIncrement,rConstitutiveMatrix, rValues, CharacteristicLength, rPlasticStrain, rBackStressVector);


            if (std::abs(threshold_indicator) <= std::abs(1.0e-4 * rThreshold)) { // Has converged
                is_converged = true;
            } else {
                iteration++;
            }
        }
        if (analytic_tangent) {
            Matrix tangent_tensor;
            tangent_tensor.resize(VoigtSize, VoigtSize, false);
            CalculateTangentMatrix(tangent_tensor, rConstitutiveMatrix, rYieldSurfaceDerivative, rDerivativePlasticPotential, rPlasticDenominator);
            noalias(rConstitutiveMatrix) = tangent_tensor;
        }
        KRATOS_WARNING_IF("GenericConstitutiveLawIntegratorKinematicPlasticity", iteration > max_iter) << "Maximum number of iterations in plasticity loop reached..." << std::endl;
    }

    /**
     * @brief This method calculates the analytical tangent tensor
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rDeviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the deviatoric part of the stress vector
     * @param rFFluxVector The derivative of the yield surface
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateTangentMatrix(
        Matrix& rTangent,
        const Matrix& rElasticMatrix,
        const array_1d<double, VoigtSize>& rFFluxVector,
        const array_1d<double, VoigtSize>& rGFluxVector,
        const double Denominator
        )
    {
        rTangent = rElasticMatrix - outer_prod(Vector(prod(rElasticMatrix, rGFluxVector)), Vector(prod(rElasticMatrix, rFFluxVector))) * Denominator;
    }


    /**
     * @brief This method calculates all the plastic parameters required for the integration of the PredictiveStressVector
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rStrainVector The equivalent strain vector of that integration point
     * @param rUniaxialStress The equivalent uniaxial stress
     * @param rThreshold The maximum uniaxial stress of the linear behaviour
     * @param rPlasticDenominator The plasticity numerical value to obtain the pastic consistency factor
     * @param rYieldSurfaceDerivative The derivative of the yield surface
     * @param rDerivativePlasticPotential The derivative of the plastic potential
     * @param rPlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param rPlasticStrainIncrement The increment of plastic strain of this time step
     * @param rConstitutiveMatrix The elastic constitutive matrix
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static double CalculatePlasticParameters(
        BoundedArrayType& rPredictiveStressVector,
        Vector& rStrainVector,
        double& rUniaxialStress,
        double& rThreshold,
        double& rPlasticDenominator,
        BoundedArrayType& rYieldSurfaceDerivative,
        BoundedArrayType& rDerivativePlasticPotential,
        double& rPlasticDissipation,
        BoundedArrayType& rPlasticStrainIncrement,
        const Matrix& rConstitutiveMatrix,
        ConstitutiveLaw::Parameters& rValues,
        const double CharacteristicLength,
        const Vector& rPlasticStrain,
        const Vector& rBackStressVector
        )
    {
        BoundedArrayType deviator = ZeroVector(VoigtSize);
        BoundedArrayType h_capa = ZeroVector(VoigtSize);
        double J2, tensile_indicator_factor, compression_indicator_factor, slope, hardening_parameter, equivalent_plastic_strain;

        YieldSurfaceType::CalculateEquivalentStress( rPredictiveStressVector, rStrainVector, rUniaxialStress, rValues);
        const double I1 = rPredictiveStressVector[0] + rPredictiveStressVector[1] + rPredictiveStressVector[2];
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(rPredictiveStressVector, I1, deviator, J2);
        CalculateDerivativeYieldSurface(rPredictiveStressVector, deviator, J2, rYieldSurfaceDerivative, rValues);
        CalculateDerivativePlasticPotential(rPredictiveStressVector, deviator, J2, rDerivativePlasticPotential, rValues);
        CalculateIndicatorsFactors(rPredictiveStressVector, tensile_indicator_factor,compression_indicator_factor);
        CalculatePlasticDissipation(rPredictiveStressVector, tensile_indicator_factor,compression_indicator_factor, rPlasticStrainIncrement,rPlasticDissipation, h_capa, rValues, CharacteristicLength);
        CalculateEquivalentPlasticStrain(rPredictiveStressVector, rUniaxialStress, rPlasticStrain, tensile_indicator_factor, rValues, equivalent_plastic_strain);
        CalculateEquivalentStressThreshold(rPlasticDissipation, tensile_indicator_factor,compression_indicator_factor, rThreshold, slope, rValues, equivalent_plastic_strain, CharacteristicLength);
        CalculateHardeningParameter(rDerivativePlasticPotential, slope, h_capa, hardening_parameter);
        CalculatePlasticDenominator(rYieldSurfaceDerivative, rDerivativePlasticPotential, rConstitutiveMatrix, hardening_parameter, rPlasticDenominator, rBackStressVector, rValues);

        return rUniaxialStress - rThreshold;
    }

    /**
     * @brief This method calculates the derivative of the yield surface
     * @param rPredictiveStressVector The predictive stress vector S = C : (E-Ep)
     * @param rDeviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the deviatoric part of the stress vector
     * @param rDerivativeYieldSurface The derivative of the yield surface
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateDerivativeYieldSurface(
        const BoundedArrayType& rPredictiveStressVector,
        const BoundedArrayType& rDeviator,
        const double J2,
        BoundedArrayType& rDerivativeYieldSurface,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        YieldSurfaceType::CalculateYieldSurfaceDerivative(rPredictiveStressVector, rDeviator, J2, rDerivativeYieldSurface, rValues);
    }

    /**
     * @brief This method computes the back stress for the kinematic plasticity
     * This method has 3 different ways of computing this back-stress:
     * Linear hardening, Armstrong-Frederick  and Araujo-Voyiadjis.
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rKinematicParameters The parameters required for the models
     * @param rPreviousStressVector The previous stress vector
     * @param rPlasticStrainIncrement The plastic strain increment of this iteration
     * @param rBackStressVector The back-stress vector for the kinematic plasticity
     */
    static void CalculateBackStress(
        BoundedArrayType& rPredictiveStressVector,
        ConstitutiveLaw::Parameters& rValues,
        const Vector& rPreviousStressVector,
        const Vector& rPlasticStrainIncrement,
        Vector& rBackStressVector
        )
    {
        const Vector& r_kinematic_parameters = rValues.GetMaterialProperties()[KINEMATIC_PLASTICITY_PARAMETERS];
        const unsigned int kinematic_hardening_type = rValues.GetMaterialProperties()[KINEMATIC_HARDENING_TYPE];

        switch (static_cast<KinematicHardeningType>(kinematic_hardening_type))
        {
            double pDot, denominator, dot_product_dp;
            case KinematicHardeningType::LinearKinematicHardening:
                KRATOS_ERROR_IF(r_kinematic_parameters.size() == 0) << "Kinematic Parameters not defined..." << std::endl;
                rBackStressVector += 2.0 / 3.0 * r_kinematic_parameters[0] * rPlasticStrainIncrement;
                break;

            case KinematicHardeningType::ArmstrongFrederickKinematicHardening:
                KRATOS_ERROR_IF(r_kinematic_parameters.size() < 2) << "Kinematic Parameters not defined..." << std::endl;
                dot_product_dp = 0.0;
                for (IndexType i = 0; i < rPlasticStrainIncrement.size(); ++i) {
                    dot_product_dp += rPlasticStrainIncrement[i] * rPlasticStrainIncrement[i];
                }
                pDot = std::sqrt(2.0 / 3.0 * dot_product_dp);
                denominator = 1.0 + (r_kinematic_parameters[1] * pDot);
                rBackStressVector = (rBackStressVector + ((2.0 / 3.0 * r_kinematic_parameters[0]) * rPlasticStrainIncrement)) / denominator;
                break;

            case KinematicHardeningType::AraujoVoyiadjisKinematicHardening:
                KRATOS_ERROR_IF(r_kinematic_parameters.size() != 3) << "Kinematic Parameters not defined..." << std::endl;
                dot_product_dp = 0.0;
                for (IndexType i = 0; i < rPlasticStrainIncrement.size(); ++i) {
                    dot_product_dp += rPlasticStrainIncrement[i] * rPlasticStrainIncrement[i];
                }
                pDot = std::sqrt(2.0 / 3.0 * dot_product_dp);
                denominator = 1.0 + (r_kinematic_parameters[1] * pDot);
                if (pDot > tolerance) {
                    rBackStressVector = (rBackStressVector + ((2.0 / 3.0 * r_kinematic_parameters[0]) * rPlasticStrainIncrement)) / denominator;
                } else {
                    const Vector& r_delta_stress = rPredictiveStressVector - rPreviousStressVector;
                    rBackStressVector = (rBackStressVector + ((2.0 / 3.0 * r_kinematic_parameters[0]) * rPlasticStrainIncrement) +
                                         r_kinematic_parameters[2] * r_delta_stress) / denominator;
                }
                break;
            default:
                KRATOS_ERROR << " The Kinematic hardening type of plasticity is not set or wrong..." << kinematic_hardening_type << std::endl;
                break;
        }
    }

    /**
     * @brief This method calculates the derivative of the plastic potential
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rDeviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the deviatoric part of the stress vector
     * @param rDerivativePlasticPotential The derivative of the yield surface
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateDerivativePlasticPotential(
        const BoundedArrayType& rPredictiveStressVector,
        const BoundedArrayType& rDeviator,
        const double J2,
        BoundedArrayType& rDerivativePlasticPotential,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        YieldSurfaceType::CalculatePlasticPotentialDerivative(rPredictiveStressVector, rDeviator, J2, rDerivativePlasticPotential, rValues);
    }

    /**
     * @brief This method computes the tensile/compressive indicators
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rTensileIndicatorFactor The tensile indicator
     * @param rCompressionIndicatorFactor The compressive indicator
     */
    static void CalculateIndicatorsFactors(
        const BoundedArrayType& rPredictiveStressVector,
        double& rTensileIndicatorFactor,
        double& rCompressionIndicatorFactor
        )
    {
        GenericConstitutiveLawIntegratorPlasticity<YieldSurfaceType>::CalculateIndicatorsFactors(
            rPredictiveStressVector,rTensileIndicatorFactor,rCompressionIndicatorFactor);
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
        const BoundedArrayType& rPredictiveStressVector,
        const double TensileIndicatorFactor,
        const double CompressionIndicatorFactor,
        const Vector& PlasticStrainInc,
        double& rPlasticDissipation,
        BoundedArrayType& rHCapa,
        ConstitutiveLaw::Parameters& rValues,
        const double CharacteristicLength
        )
    {
        GenericConstitutiveLawIntegratorPlasticity<YieldSurfaceType>::CalculatePlasticDissipation(
            rPredictiveStressVector,TensileIndicatorFactor,CompressionIndicatorFactor,
            PlasticStrainInc,rPlasticDissipation,rHCapa,rValues,CharacteristicLength);
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
     * @param CharacteristicLength Characteristic length of the finite element
     */
    static void CalculateEquivalentStressThreshold(
        const double PlasticDissipation,
        const double TensileIndicatorFactor,
        const double CompressionIndicatorFactor,
        double& rEquivalentStressThreshold,
        double& rSlope,
        ConstitutiveLaw::Parameters& rValues,
        const double EquivalentPlasticStrain,
        const double CharacteristicLength
        )
    {
        GenericConstitutiveLawIntegratorPlasticity<YieldSurfaceType>::CalculateEquivalentStressThreshold(
            PlasticDissipation,TensileIndicatorFactor,CompressionIndicatorFactor,rEquivalentStressThreshold,
            rSlope,rValues,EquivalentPlasticStrain,CharacteristicLength);
    }

    /**
     * @brief This method computes the uniaxial threshold using a linear softening
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param TensileIndicatorFactor The tensile indicator
     * @param CompressionIndicatorFactor The compressive indicator
     * @param rEquivalentStressThreshold The maximum uniaxial stress of the linear behaviour
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
        GenericConstitutiveLawIntegratorPlasticity<YieldSurfaceType>::CalculateEquivalentStressThresholdHardeningCurveLinearSoftening(
            PlasticDissipation,TensileIndicatorFactor,CompressionIndicatorFactor,rEquivalentStressThreshold,rSlope,rValues);
    }

    /**
     * @brief This method computes the uniaxial threshold using a exponential softening
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param TensileIndicatorFactor The tensile indicator
     * @param CompressionIndicatorFactor The compressive indicator
     * @param rEquivalentStressThreshold The maximum uniaxial stress of the linear behaviour
     * @param rSlope The slope of the PlasticDiss-Threshold curve
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateEquivalentStressThresholdHardeningCurveExponentialSoftening(
        const double PlasticDissipation,
        const double TensileIndicatorFactor,
        const double CompressionIndicatorFactor,
        double& rEquivalentStressThreshold,
        double& rSlope,
        ConstitutiveLaw::Parameters& rValues,
        const double CharacteristicLength
        )
    {
        GenericConstitutiveLawIntegratorPlasticity<YieldSurfaceType>::CalculateEquivalentStressThresholdHardeningCurveExponentialSoftening(
            PlasticDissipation,TensileIndicatorFactor,CompressionIndicatorFactor,rEquivalentStressThreshold,rSlope,rValues,CharacteristicLength);
    }

    /**
     * @brief This method computes the uniaxial threshold using a hardening-softening law
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param TensileIndicatorFactor The tensile indicator
     * @param CompressionIndicatorFactor The compressive indicator
     * @param rEquivalentStressThreshold The maximum uniaxial stress of the linear behaviour
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
        GenericConstitutiveLawIntegratorPlasticity<YieldSurfaceType>::CalculateEquivalentStressThresholdHardeningCurveInitialHardeningExponentialSoftening(
            PlasticDissipation,TensileIndicatorFactor,CompressionIndicatorFactor,rEquivalentStressThreshold,rSlope,rValues);
    }

    /**
     * @brief This method computes the uniaxial threshold using a perfect plasticity law
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param TensileIndicatorFactor The tensile indicator
     * @param CompressionIndicatorFactor The compressive indicator
     * @param rEquivalentStressThreshold The maximum uniaxial stress of the linear behaviour
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
        GenericConstitutiveLawIntegratorPlasticity<YieldSurfaceType>::CalculateEquivalentStressThresholdHardeningCurvePerfectPlasticity(
            PlasticDissipation,TensileIndicatorFactor,CompressionIndicatorFactor,rEquivalentStressThreshold,rSlope,rValues);
    }

    /**
     * @brief This method computes the uniaxial threshold using a perfect plasticity law
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param TensileIndicatorFactor The tensile indicator
     * @param CompressionIndicatorFactor The compressive indicator
     * @param rEquivalentStressThreshold The maximum uniaxial stress of the linear behaviour
     * @param rSlope The slope of the PlasticDiss-Threshold curve
     * @param rValues Parameters of the constitutive law
     * @param EquivalentPlasticStrain The Plastic Strain internal variable
     * @param CharacteristicLength Characteristic length of the finite element
     */
    static void CalculateEquivalentStressThresholdCurveFittingHardening(
        const double PlasticDissipation,
        const double TensileIndicatorFactor,
        const double CompressionIndicatorFactor,
        double& rEquivalentStressThreshold,
        double& rSlope,
        ConstitutiveLaw::Parameters& rValues,
        const double EquivalentPlasticStrain,
        const double CharacteristicLength
        )
    {
        GenericConstitutiveLawIntegratorPlasticity<YieldSurfaceType>::CalculateEquivalentStressThresholdCurveFittingHardening(
            PlasticDissipation,TensileIndicatorFactor,CompressionIndicatorFactor,rEquivalentStressThreshold,rSlope,rValues,
            EquivalentPlasticStrain,CharacteristicLength);
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
        GenericConstitutiveLawIntegratorPlasticity<YieldSurfaceType>::CalculateEquivalentPlasticStrain(
            rStressVector,UniaxialStress,rPlasticStrain,r0,rValues,rEquivalentPlasticStrain);
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
        const BoundedArrayType& rGFlux,
        const double SlopeThreshold,
        const BoundedArrayType& rHCapa,
        double& rHardeningParameter
        )
    {
        GenericConstitutiveLawIntegratorPlasticity<YieldSurfaceType>::CalculateHardeningParameter(
            rGFlux,SlopeThreshold,rHCapa,rHardeningParameter);
    }

    /**
     * @brief This method computes the plastic denominator needed
     * to compute the plastic consistency factor
     * @param rYieldSurfaceDerivative The derivative of the yield surface
     * @param rDerivativePlasticPotential The derivative of the plastic potential
     * @param rConstitutiveMatrix The elastic constitutive matrix
     * @param rHardeningParameter The hardening parameter needed for the algorithm
     * @param rPlasticDenominator The plasticity numerical value to obtain the pastic consistency factor
     */
    static void CalculatePlasticDenominator(
        const BoundedArrayType& rFFlux,
        const BoundedArrayType& rGFlux,
        const Matrix& rConstitutiveMatrix,
        double& rHardeningParameter,
        double& rPlasticDenominator,
        const Vector& rBackStressVector,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        const Vector& r_kinematic_parameters = rValues.GetMaterialProperties()[KINEMATIC_PLASTICITY_PARAMETERS];
        const int kinematic_hardening_type = rValues.GetMaterialProperties()[KINEMATIC_HARDENING_TYPE];

        const BoundedArrayType delta_vector = prod(rGFlux, rConstitutiveMatrix);
        double A1 = 0.0;
        for (IndexType i = 0; i < VoigtSize; ++i) {
            A1 += rFFlux[i] * delta_vector[i];
        }
        if (r_kinematic_parameters.size() == 3) {
            A1 *= (1.0 - r_kinematic_parameters[2]);
        } // Araujo case with 3 params

        double dot_fflux_gflux = 0.0, A2;
        for (IndexType i = 0; i < VoigtSize; ++i) {
            dot_fflux_gflux += rFFlux[i] * rGFlux[i];
        }
        const double two_thirds = 2.0 / 3.0;
        double dot_fflux_backstress = 0.0, dot_gflux_gflux = 0.0;
        switch (static_cast<KinematicHardeningType>(kinematic_hardening_type))
        {
            case KinematicHardeningType::LinearKinematicHardening:
                A2 = two_thirds * r_kinematic_parameters[0] * dot_fflux_gflux;
                break;

            case KinematicHardeningType::ArmstrongFrederickKinematicHardening:
                A2 = two_thirds * r_kinematic_parameters[0] * dot_fflux_gflux;
                for (IndexType i = 0; i < VoigtSize; ++i) {
                    dot_fflux_backstress += rFFlux[i] * rBackStressVector[i];
                }
                for (IndexType i = 0; i < VoigtSize; ++i) {
                    dot_gflux_gflux += rGFlux[i] * rGFlux[i];
                }
                A2 -= r_kinematic_parameters[1] * dot_fflux_backstress * std::sqrt(two_thirds * dot_gflux_gflux);
                break;

            case KinematicHardeningType::AraujoVoyiadjisKinematicHardening:
                A2 = two_thirds * r_kinematic_parameters[0] * dot_fflux_gflux;
                for (IndexType i = 0; i < VoigtSize; ++i) {
                    dot_fflux_backstress += rFFlux[i] * rBackStressVector[i];
                }
                for (IndexType i = 0; i < VoigtSize; ++i) {
                    dot_gflux_gflux += rGFlux[i] * rGFlux[i];
                }
                A2 -= r_kinematic_parameters[1] * dot_fflux_backstress * std::sqrt(two_thirds * dot_gflux_gflux);
                break;

            default:
                KRATOS_ERROR << " The Kinematic hardening type of plasticity is not set or wrong..." << kinematic_hardening_type << std::endl;
                break;
        }

        const double A3 = rHardeningParameter;
        rPlasticDenominator = 1.0 / (A1 + A2 + A3);

        if (r_kinematic_parameters.size() == 3) {
           rPlasticDenominator *= (1.0 - r_kinematic_parameters[2]);
        } // Araujo case with 3 params
    }

    /**
     * @brief This method defines in the CL integrator
     * @return 0 if OK, 1 otherwise
     */
    static int Check(const Properties& rMaterialProperties)
    {
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
