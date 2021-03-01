// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Alejandro Cornejo
//                   Lucia Barbu
//

#if !defined(KRATOS_GENERIC_FINITE_STRAIN_CONSTITUTIVE_LAW_INTEGRATOR_KINEMATIC_PLASTICITY_H_INCLUDED)
#define KRATOS_GENERIC_FINITE_STRAIN_CONSTITUTIVE_LAW_INTEGRATOR_KINEMATIC_PLASTICITY_H_INCLUDED

// System includes

// Project includes
#include "includes/define.h"
#include "includes/checks.h"
#include "includes/properties.h"
#include "utilities/math_utils.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_kinematic_plasticity.h"

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
 * @class GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity
 * @ingroup StructuralMechanicsApplication
 * @brief This object integrates the predictive stress using the plasticity theory by means of
 * linear/exponential softening or hardening + softening evolution laws
 * @details The crucial difference between the discretisation of the large strain problem and the infinitesimal one lies in the numerical approximation of the plastic flow equation. The structure of the plastic flow equation makes algorithms based on exponential map integrators ideal for numerical approximation. (COMPUTATIONAL METHODS FOR PLASTICITY THEORY AND APPLICATIONS. EA de Souza Neto,D Perić, DRJ Owen pag. 616).
 * The definitions of these classes is completely static, the derivation is done in a static way
 * @tparam TYieldSurfaceType The yield surface considered
 * The plasticity integrator requires the definition of the following properties:
 * - SOFTENING_TYPE: The softening behaviour considered (linear, exponential,etc...)
 * - HARDENING_CURVE: The type of considered hardening curve (linear, exponential, pure plastic, etc...)
 * - MAXIMUM_STRESS: The maximum stress that defines the exponential hardening
 * - MAXIMUM_STRESS_POSITION: The maximum stress position that defines the exponential hardening
 * - FRACTURE_ENERGY: A fracture energy-based function is used to describe strength degradation in post-peak regime
 * - YOUNG_MODULUS: It defines the relationship between stress (force per unit area) and strain (proportional deformation) in a material in the linear elasticity regime of a uniaxial deformation.
 * - YIELD_STRESS: Yield stress is the amount of stress that an object needs to experience for it to be permanently deformed. Does not require to be defined simmetrically, one YIELD_STRESS_COMPRESSION and other YIELD_STRESS_TENSION can be defined for not symmetric cases
 * @author Vicente Mataix Ferrandiz
 * @author Alejandro Cornejo
 * @author Lucia Barbu
 */
template<class TYieldSurfaceType>
class GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity
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

    /// Small strain plasticity law integrator
    typedef GenericConstitutiveLawIntegratorKinematicPlasticity<YieldSurfaceType> SmallStrainIntegratorType;

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

    /// Definition of the tolerance for convergence
    static constexpr double ConvergenceTolerance = 1.0e-4;

    /// Counted pointer of GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity
    KRATOS_CLASS_POINTER_DEFINITION(GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity);

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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor
    GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity()
    {
    }

    /// Copy constructor
    GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity(GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity const &rOther)
    {
    }

    /// Assignment operator
    GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity &operator=(GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~GenericFiniteStrainConstitutiveLawIntegratorKinematicPlasticity()
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
     * @param rConstitutiveLaw The constitutive law necessary to compute the sigma
     * @param rPredictiveStressVector The predictive stress vector
     * @param rUniaxialStress The equivalent uniaxial stress
     * @param rThreshold The maximum uniaxial stress of the linear behaviour
     * @param rPlasticDenominator The plasticity numerical value to obtain the pastic consistency factor
     * @param rYieldSurfaceDerivative The derivative of the yield surface
     * @param rPlasticPotentialDerivative The derivative of the plastic potential
     * @param rPlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param rPlasticDeformationGradient The plastic deformation gradient
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     * @param rBackStressVector The so-called BackStressVector used for the kinematic hardening
     */
    static void IntegrateStressVector(
        ConstitutiveLaw& rConstitutiveLaw,
        const Variable<Vector>& rStrainVariable,
        const Variable<Vector>& rStressVariable,
        BoundedArrayType& rPredictiveStressVector,
        Vector& rStrainVector,
        double& rUniaxialStress,
        double& rThreshold,
        double& rPlasticDenominator,
        BoundedArrayType& rYieldSurfaceDerivative,
        BoundedArrayType& rPlasticPotentialDerivative,
        double& rPlasticDissipation,
        Matrix& rPlasticDeformationGradient,
        const Matrix& rPreviousDeformationGradient,
        const Matrix& rConstitutiveMatrix,
        ConstitutiveLaw::Parameters& rValues,
        const double CharacteristicLength,
        Vector& rBackStressVector
        )
    {
        // We do a backup of the parameters
        const ConstitutiveLaw::Parameters values_backup = rValues;

        // Material properties
        const Properties& r_material_properties = rValues.GetMaterialProperties();

        // Some values initialization
        IndexType iteration = 0, max_iter = r_material_properties.Has(MAX_NUMBER_NL_CL_ITERATIONS) ? r_material_properties.GetValue(MAX_NUMBER_NL_CL_ITERATIONS) : 100;
        double plastic_consistency_factor_increment = 0.0;
        double threshold_indicator = rUniaxialStress - rThreshold;
        Vector delta_sigma = ZeroVector(VoigtSize);
        Vector plastic_strain = ZeroVector(VoigtSize);
        Vector delta_plastic_strain = ZeroVector(VoigtSize);

        // Initialize the pastic deformation gradient increment
        double aux_det;
        Vector aux_vector(VoigtSize);
        Matrix Re(Dimension, Dimension), Ue(Dimension, Dimension);
        Matrix inverse_previous_deformation_gradient_increment(Dimension, Dimension);
        const Matrix previous_plastic_deformation_gradient = rPlasticDeformationGradient;
        Matrix plastic_deformation_gradient_increment(Dimension, Dimension);
        Matrix inverse_plastic_deformation_gradient(Dimension, Dimension);

        // Predictive deformation gradient
        Matrix predictive_deformation_gradient = rValues.GetDeformationGradientF();

        // Calculate the incremental deformation gradient
        MathUtils<double>::InvertMatrix(rPreviousDeformationGradient, inverse_previous_deformation_gradient_increment, aux_det);
        const Matrix incremental_deformation_gradient = prod(predictive_deformation_gradient, inverse_previous_deformation_gradient_increment);

        // The inverse of the plastic deformation gradient
        MathUtils<double>::InvertMatrix(previous_plastic_deformation_gradient, inverse_plastic_deformation_gradient, aux_det);

        // We compute Fen auxiliarly
        const Matrix previous_elastic_deformation_matrix = prod(rPreviousDeformationGradient, inverse_plastic_deformation_gradient);

        // We compute Fen+1 trial auxiliarly
        const Matrix trial_elastic_deformation_matrix = prod(incremental_deformation_gradient, previous_elastic_deformation_matrix);

        // We compute Fe auxiliarly
        Matrix elastic_deformation_matrix = previous_elastic_deformation_matrix;

        // We compute the previous stress vector
        rValues.SetDeterminantF(MathUtils<double>::Det(rPreviousDeformationGradient));
        rValues.SetDeformationGradientF(rPreviousDeformationGradient);
        Vector previous_stress_vector(VoigtSize);
        rConstitutiveLaw.CalculateValue(rValues, rStressVariable, previous_stress_vector);

        // Definition of the kinematic hardening stress vector
        BoundedArrayType kin_hard_stress_vector;

        // Backward Euler
        while (iteration <= max_iter) {
            // With this we can compute the polar decomposition in order to compute the Ren+1
            ConstitutiveLawUtilities<VoigtSize>::PolarDecomposition(elastic_deformation_matrix, Re, Ue);

            noalias(elastic_deformation_matrix) = ConstitutiveLawUtilities<VoigtSize>::CalculateDirectElasticDeformationGradient(trial_elastic_deformation_matrix, rPlasticPotentialDerivative, plastic_consistency_factor_increment, Re);

            // We calculate the plastic consistency factor
            plastic_consistency_factor_increment = threshold_indicator * rPlasticDenominator;
            if (plastic_consistency_factor_increment < 0.0) plastic_consistency_factor_increment = 0.0; // NOTE: It could be useful, maybe

//             noalias(plastic_deformation_gradient_increment) = ConstitutiveLawUtilities<VoigtSize>::CalculateLinearPlasticDeformationGradientIncrement(rPlasticPotentialDerivative, plastic_consistency_factor_increment);
            noalias(plastic_deformation_gradient_increment) = ConstitutiveLawUtilities<VoigtSize>::CalculateDirectPlasticDeformationGradientIncrement(rPlasticPotentialDerivative, plastic_consistency_factor_increment, Re);
//             noalias(plastic_deformation_gradient_increment) = ConstitutiveLawUtilities<VoigtSize>::CalculateExponentialPlasticDeformationGradientIncrement(rPlasticPotentialDerivative, plastic_consistency_factor_increment, Re);

            // We check that the increment is not a zero matrix or to big
            const double norm_increment = norm_frobenius(plastic_deformation_gradient_increment);
            if (norm_increment < 1.0e-8
                || norm_increment > 1.0e2
                ) {
                Vector& r_strain_vector = rValues.GetStrainVector();
                rConstitutiveLaw.CalculateValue(rValues, rStrainVariable, r_strain_vector);
                rConstitutiveLaw.CalculateValue(rValues, rStressVariable, aux_vector);
                noalias(rPredictiveStressVector) = aux_vector;
                break;
            }

            // The increment of the deformation is not added but multiplied in finite strain
            noalias(rPlasticDeformationGradient) = prod(plastic_deformation_gradient_increment, previous_plastic_deformation_gradient);
            rValues.SetDeterminantF(MathUtils<double>::Det(plastic_deformation_gradient_increment));
            rValues.SetDeformationGradientF(plastic_deformation_gradient_increment);
            rConstitutiveLaw.CalculateValue(rValues, rStrainVariable, delta_plastic_strain);

            // We compute the plastic strain
            rValues.SetDeterminantF(MathUtils<double>::Det(rPlasticDeformationGradient));
            rValues.SetDeformationGradientF(rPlasticDeformationGradient);
            rConstitutiveLaw.CalculateValue(rValues, rStrainVariable, plastic_strain);

            // We compute the new predictive stress vector
            MathUtils<double>::InvertMatrix(plastic_deformation_gradient_increment, inverse_plastic_deformation_gradient, aux_det);
            predictive_deformation_gradient = prod(inverse_plastic_deformation_gradient, predictive_deformation_gradient);
            rValues.SetDeterminantF(MathUtils<double>::Det(predictive_deformation_gradient));
            rValues.SetDeformationGradientF(predictive_deformation_gradient);
            rConstitutiveLaw.CalculateValue(rValues, rStressVariable, aux_vector);
            noalias(rPredictiveStressVector) = aux_vector;

            // Calculate back-stress
            SmallStrainIntegratorType::CalculateBackStress(rPredictiveStressVector, rValues, previous_stress_vector,
                                                           delta_plastic_strain, rBackStressVector);

            // Calculate plastic parameters
            noalias(kin_hard_stress_vector) = rPredictiveStressVector - rBackStressVector;
            threshold_indicator = CalculatePlasticParameters(kin_hard_stress_vector, rStrainVector, rUniaxialStress, rThreshold,
                                                             rPlasticDenominator, rYieldSurfaceDerivative, rPlasticPotentialDerivative,
                                                             rPlasticDissipation, delta_plastic_strain, rConstitutiveMatrix, rValues, CharacteristicLength, plastic_strain, rBackStressVector);

            if (threshold_indicator <= std::abs(1.0e-4 * rThreshold)) { // Has converged
                break;
            } else {
                ++iteration;
            }
        }

        rValues = values_backup;

        if (iteration > max_iter) {
            KRATOS_WARNING_FIRST_N("Backward Euler Plasticity", 20) << "Maximum number of iterations in plasticity loop reached..." << std::endl;
        }
    }

    /**
     * @brief This method calculates all the plastic parameters required for the integration of the PredictiveStressVector
     * @param rPredictiveStressVector The predictive stress vector
     * @param rUniaxialStress The equivalent uniaxial stress
     * @param rThreshold The maximum uniaxial stress of the linear behaviour
     * @param rPlasticDenominator The plasticity numerical value to obtain the pastic consistency factor
     * @param rYieldSurfaceDerivative The derivative of the yield surface
     * @param rDerivativePlasticPotential The derivative of the plastic potential
     * @param rPlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param rPlasticStrainIncrement The increment of plastic strain of this time step
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     * @return It returns the threshold indicator
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
        const Vector& rPlasticStrainIncrement,
        const Matrix& rConstitutiveMatrix,
        ConstitutiveLaw::Parameters& rValues,
        const double CharacteristicLength,
        const Vector& rPlasticStrain,
        const Vector& rBackStressVector
        )
    {
        BoundedArrayType deviator = ZeroVector(6);
        BoundedArrayType h_capa = ZeroVector(6);
        double J2, tensile_indicator_factor, compression_indicator_factor, slope, hardening_parameter, equivalent_plastic_strain;

        YieldSurfaceType::CalculateEquivalentStress(rPredictiveStressVector, rStrainVector, rUniaxialStress, rValues);
        const double I1 = rPredictiveStressVector[0] + rPredictiveStressVector[1] + rPredictiveStressVector[2];
        ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(rPredictiveStressVector, I1, deviator, J2);
        CalculateDerivativeYieldSurface(rPredictiveStressVector, deviator, J2, rYieldSurfaceDerivative, rValues);
        CalculateDerivativePlasticPotential(rPredictiveStressVector, deviator, J2, rDerivativePlasticPotential, rValues);
        CalculateIndicatorsFactors(rPredictiveStressVector, tensile_indicator_factor, compression_indicator_factor);
        CalculatePlasticDissipation(rPredictiveStressVector, tensile_indicator_factor, compression_indicator_factor, rPlasticStrainIncrement, rPlasticDissipation, h_capa, rValues, CharacteristicLength);
        CalculateEquivalentPlasticStrain(rPredictiveStressVector, rUniaxialStress, rPlasticStrain, tensile_indicator_factor, rValues, equivalent_plastic_strain);
        CalculateEquivalentStressThreshold(rPlasticDissipation, tensile_indicator_factor, compression_indicator_factor, rThreshold, slope, equivalent_plastic_strain, rValues, CharacteristicLength);
        CalculateHardeningParameter(rDerivativePlasticPotential, slope, h_capa, hardening_parameter);
        CalculatePlasticDenominator(rYieldSurfaceDerivative, rDerivativePlasticPotential, rConstitutiveMatrix, hardening_parameter, rPlasticDenominator, rBackStressVector, rValues);

        return rUniaxialStress - rThreshold;
    }

    /**
     * @brief This method calculates the derivative of the yield surface
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
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
     * @brief This method calculates the derivative of the plastic potential
     * @param rDeviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the deviatoric part of the stress vector
     * @param rDerivativePlasticPotential The derivative of the plastic potential
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateDerivativePlasticPotential(
        BoundedArrayType& rPredictiveStressVector,
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
        SmallStrainIntegratorType::CalculateIndicatorsFactors(rPredictiveStressVector, rTensileIndicatorFactor, rCompressionIndicatorFactor);
    }

    /**
     * @brief This method computes the plastic dissipation of the plasticity model
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param TensileIndicatorFactor The tensile indicator
     * @param CompressionIndicatorFactor The compressive indicator
     * @param rPlasticStrainIncrement The increment of plastic strain of this time step
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param rHCapa The slope of the PlasticDiss-Threshold curve
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculatePlasticDissipation(
        const BoundedArrayType& rPredictiveStressVector,
        const double TensileIndicatorFactor,
        const double CompressionIndicatorFactor,
        const Vector& rPlasticStrainIncrement,
        double& rPlasticDissipation,
        BoundedArrayType& rHCapa,
        ConstitutiveLaw::Parameters& rValues,
        const double CharacteristicLength
        )
    {
        SmallStrainIntegratorType::CalculatePlasticDissipation(rPredictiveStressVector, TensileIndicatorFactor, CompressionIndicatorFactor, rPlasticStrainIncrement, rPlasticDissipation, rHCapa, rValues, CharacteristicLength);
    }

    /**
     * @brief This method computes the uniaxial threshold that differentiates the elastic-plastic behaviour
     * @param PlasticDissipation The internal variable of energy dissipation due to plasticity
     * @param TensileIndicatorFactor The tensile indicator
     * @param CompressionIndicatorFactor The compressive indicator
     * @param rEquivalentStressThreshold The maximum uniaxial stress of the linear behaviour
     * @param rSlope The slope of the PlasticDiss-Threshold curve
     * @param EquivalentPlasticStrain The equivalent plastic strain
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateEquivalentStressThreshold(
        const double PlasticDissipation,
        const double TensileIndicatorFactor,
        const double CompressionIndicatorFactor,
        double& rEquivalentStressThreshold,
        double& rSlope,
        const double EquivalentPlasticStrain,
        ConstitutiveLaw::Parameters& rValues,
        const double CharacteristicLength
        )
    {
        SmallStrainIntegratorType::CalculateEquivalentStressThreshold(PlasticDissipation, TensileIndicatorFactor, CompressionIndicatorFactor, rEquivalentStressThreshold, rSlope, rValues, EquivalentPlasticStrain, CharacteristicLength);
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
        SmallStrainIntegratorType::CalculateEquivalentPlasticStrain(rStressVector, UniaxialStress, rPlasticStrain, r0, rValues, rEquivalentPlasticStrain);
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
     * @param rDerivativePlasticPotential The derivative of the plastic potential
     * @param SlopeThreshold The slope of the PlasticDiss-Threshold curve
     * @param rHardeningParameter The hardening parameter needed for the algorithm
     */
    static void CalculateHardeningParameter(
        const BoundedArrayType& rDerivativePlasticPotential,
        const double SlopeThreshold,
        const BoundedArrayType& rHCapa,
        double& rHardeningParameter
        )
    {
        SmallStrainIntegratorType::CalculateHardeningParameter(rDerivativePlasticPotential, SlopeThreshold, rHCapa, rHardeningParameter);
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
        const array_1d<double, VoigtSize>& rYieldSurfaceDerivative,
        const array_1d<double, VoigtSize>& rDerivativePlasticPotential,
        const Matrix& rConstitutiveMatrix,
        double& rHardeningParameter,
        double& rPlasticDenominator,
        const Vector& rBackStressVector,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        SmallStrainIntegratorType::CalculatePlasticDenominator(rYieldSurfaceDerivative, rDerivativePlasticPotential, rConstitutiveMatrix, rHardeningParameter, rPlasticDenominator, rBackStressVector, rValues);
    }

    /**
     * @brief This method defines in the CL integrator
     * @return 0 if OK, 1 otherwise
     */
    static int Check(const Properties& rMaterialProperties)
    {
        return SmallStrainIntegratorType::Check(rMaterialProperties);
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
