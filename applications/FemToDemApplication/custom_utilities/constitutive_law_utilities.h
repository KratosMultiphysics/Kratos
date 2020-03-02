//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//                   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_CONSTITUTIVE_LAW_UTILITIES)
#define KRATOS_CONSTITUTIVE_LAW_UTILITIES

// System includes

// External includes

// Project includes
#include "includes/ublas_interface.h"
#include "includes/node.h"
#include "geometries/geometry.h"
#include "includes/constitutive_law.h"
#include "fem_to_dem_application_variables.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

/// The size type definition
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
 * @class ConstitutiveLawUtilities
 * @ingroup StructuralMechanicsApplication / FemToDemApplication
 * @brief This class includes several utilities necessaries for the computation of the constitutive law
 * @details The methods are static, so it can be called without constructing the class
 * @tparam TVoigtSize The number of components on the Voigt notation
 * @author Alejandro Cornejo
 * @author Vicente Mataix Ferrandiz
 */
template <SizeType TVoigtSize = 6>
class ConstitutiveLawUtilities
{
  public:
    ///@name Type definitions
    ///@{

    /// The index type definition
    typedef std::size_t IndexType;

    /// We define the dimension
    static constexpr SizeType Dimension = TVoigtSize == 6 ? 3 : 2;

    /// We define the Voigt size
    static constexpr SizeType VoigtSize = TVoigtSize;

    /// The matrix type definition
    typedef Matrix MatrixType;

    /// the vector type definition
    typedef Vector VectorType;

    /// The definition of the bounded vector type
    typedef array_1d<double, VoigtSize> BoundedVectorType;

    /// The definition of the bounded matrix type
    typedef BoundedMatrix<double, Dimension, Dimension> BoundedMatrixType;

    /// Node type definition
    typedef Node<3> NodeType;

    /// Geometry definitions
    typedef Geometry<NodeType> GeometryType;

    /// The zero tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method computes the first invariant from a given stress vector
     * @param rStressVector The stress vector on Voigt notation
     * @param rI1 The first invariant
     */
    static void CalculateI1Invariant(
        const BoundedVectorType& rStressVector,
        double& rI1
        );

    /**
     * @brief This method computes the second invariant from a given stress vector
     * @param rStressVector The stress vector on Voigt notation
     * @param rI2 The second invariant
     * @todo Adapt for 2D dimension
     */
    static void CalculateI2Invariant(
        const BoundedVectorType& rStressVector,
        double& rI2
        );

    /**
     * @brief This method computes the third invariant from a given stress vector
     * @param rStressVector The stress vector on Voigt notation
     * @param rI3 The third invariant
     * @todo Adapt for 2D dimension
     */
    static void CalculateI3Invariant(
        const BoundedVectorType& rStressVector,
        double& rI3
        );

    /**
     * @brief This method computes the second invariant of J
     * @param rStressVector The stress vector on Voigt notation
     * @param I1 The first invariant
     * @param rDeviator The deviator of the stress
     * @param rJ2 The second invariant of J
     */
    static void CalculateJ2Invariant(
        const BoundedVectorType& rStressVector,
        const double I1,
        BoundedVectorType& rDeviator,
        double& rJ2
        );

    /**
     * @brief This method computes the third invariant of J
     * @param rDeviator The deviator of the stress
     * @param rJ3 The third invariant of J
     */
    static void CalculateJ3Invariant(
        const BoundedVectorType& rDeviator,
        double& rJ3
        );

    /**
     * @brief This method computes the first vector
     * @param rFirstVector The first vector
     */
    static void CalculateFirstVector(BoundedVectorType& rFirstVector);

    /**
     * @brief This method computes the second vector
     * @param rDeviator The deviator of the stress
     * @param J2 The resultant J2 stress
     * @param rSecondVector The second vector
     */
    static void CalculateSecondVector(
        const BoundedVectorType& rDeviator,
        const double J2,
        BoundedVectorType& rSecondVector
        );

    /**
     * @brief This method computes the third vector
     * @param rDeviator The deviator of the stress
     * @param J2 The resultant J2 stress
     * @param rThirdVector The third vector
     * @todo Adapt for 2D dimension
     */
    static void CalculateThirdVector(
        const BoundedVectorType& rDeviator,
        const double J2,
        BoundedVectorType& rThirdVector
        );

    /**
     * @brief This method computes the lode angle
     * @param J2 The resultant J2 stress
     * @param J3 The resultant J3 stress
     * @param rLodeAngle The lode angle
     */
    static void CalculateLodeAngle(
        const double J2,
        const double J3,
        double& rLodeAngle
        );

    /**
     * @brief This method computes the principal stresses vector
     * @details http://www.continuummechanics.org/principalstress.html
     * @param rPrincipalStressVector The vector of principal stresses
     * @param rStressVector The vector of stresses
     * @todo Adapt for 2D dimension
     */
    static void CalculatePrincipalStresses(
        array_1d<double, Dimension>& rPrincipalStressVector,
        const BoundedVectorType& rStressVector
        );

    /**
     * @brief This method computes the principal stresses vector
     * @details Using Cardano formula and renormalizing (TODO)
     * @param rPrincipalStressVector The vector of principal stresses
     * @param rStressVector The vector of stresses
     * @todo Adapt for 2D dimension
     */
    static void CalculatePrincipalStressesWithCardano(
        array_1d<double, Dimension>& rPrincipalStressVector,
        const BoundedVectorType& rStressVector
        );

    /**
     * @brief This method the uniaxial equivalent stress of Huber-VonMises
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rStrainVector The StrainVector vector
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateEquivalentStressHuberVonMises(
        const BoundedVectorType& rPredictiveStressVector,
        const Vector& rStrainVector,
        double& rEquivalentStress,
        ConstitutiveLaw::Parameters& rValues)
    {
        double I1, J2;
        array_1d<double, VoigtSize> deviator = ZeroVector(VoigtSize);
        ConstitutiveLawUtilities<VoigtSize>::CalculateI1Invariant(rPredictiveStressVector, I1);
        ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(rPredictiveStressVector, I1, deviator, J2);
        rEquivalentStress = std::sqrt(3.0 * J2);
    }

    /**
     * @brief This method the uniaxial equivalent stress of ModifiedMohrCoulomb
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rStrainVector The StrainVector vector
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateEquivalentStressModifiedMohrCoulomb(
        const BoundedVectorType& rPredictiveStressVector,
        const Vector& rStrainVector,
        double& rEquivalentStress,
        ConstitutiveLaw::Parameters& rValues)
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();
        const double yield_compression = r_material_properties[YIELD_STRESS_C];
        const double yield_tension = r_material_properties[YIELD_STRESS_T];
        double friction_angle = r_material_properties[INTERNAL_FRICTION_ANGLE] * Globals::Pi / 180.0; // In radians!

        // Check input variables
        if (friction_angle < tolerance) {
            friction_angle = 32.0 * Globals::Pi / 180.0;
            KRATOS_WARNING("ModifiedMohrCoulombYieldSurface") << "Friction Angle not defined, assumed equal to 32 deg " << std::endl;
        }

        double theta;
        const double R = std::abs(yield_compression / yield_tension);
        const double Rmorh = std::pow(std::tan((Globals::Pi / 4.0) + friction_angle / 2.0), 2);
        const double alpha_r = R / Rmorh;
        const double sin_phi = std::sin(friction_angle);

        double I1, J2, J3;
        array_1d<double, TVoigtSize> deviator = ZeroVector(TVoigtSize);
        ConstitutiveLawUtilities<TVoigtSize>::CalculateI1Invariant(rPredictiveStressVector, I1);
        ConstitutiveLawUtilities<TVoigtSize>::CalculateJ2Invariant(rPredictiveStressVector, I1, deviator, J2);
        ConstitutiveLawUtilities<TVoigtSize>::CalculateJ3Invariant(deviator, J3);

        const double K1 = 0.5 * (1.0 + alpha_r) - 0.5 * (1.0 - alpha_r) * sin_phi;
        const double K2 = 0.5 * (1.0 + alpha_r) - 0.5 * (1.0 - alpha_r) / sin_phi;
        const double K3 = 0.5 * (1.0 + alpha_r) * sin_phi - 0.5 * (1.0 - alpha_r);

        // Check Modified Mohr-Coulomb criterion
        if (std::abs(I1) < tolerance) {
            rEquivalentStress = 0.0;
        } else {
            ConstitutiveLawUtilities<TVoigtSize>::CalculateLodeAngle(J2, J3, theta);
            rEquivalentStress = (2.0 * std::tan(Globals::Pi * 0.25 + friction_angle * 0.5) / std::cos(friction_angle)) * ((I1 * K3 / 3.0) +
                        std::sqrt(J2) * (K1 * std::cos(theta) - K2 * std::sin(theta) * sin_phi / std::sqrt(3.0)));
        }
    }

    /**
     * @brief This method the uniaxial equivalent stress of MohrCoulomb
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rStrainVector The StrainVector vector
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateEquivalentStressMohrCoulomb(
        const BoundedVectorType& rPredictiveStressVector,
        const Vector& rStrainVector,
        double& rEquivalentStress,
        ConstitutiveLaw::Parameters& rValues)
    {
        double I1, J2, J3, lode_angle;
        array_1d<double, TVoigtSize> deviator = ZeroVector(TVoigtSize);

        ConstitutiveLawUtilities<TVoigtSize>::CalculateI1Invariant(rPredictiveStressVector, I1);
        ConstitutiveLawUtilities<TVoigtSize>::CalculateJ2Invariant(rPredictiveStressVector, I1, deviator, J2);
        ConstitutiveLawUtilities<TVoigtSize>::CalculateJ3Invariant(deviator, J3);
        ConstitutiveLawUtilities<TVoigtSize>::CalculateLodeAngle(J2, J3, lode_angle);

        const Properties& r_material_properties = rValues.GetMaterialProperties();
        const double friction_angle = r_material_properties[INTERNAL_FRICTION_ANGLE] * Globals::Pi / 180.0;

        rEquivalentStress = (std::cos(lode_angle) - std::sin(lode_angle) * std::sin(friction_angle) / std::sqrt(3.0)) * std::sqrt(J2) +
            I1 * std::sin(friction_angle) / 3.0;
    }

    /**
     * @brief This method the uniaxial equivalent stress of Rankine
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rStrainVector The StrainVector vector
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateEquivalentStressRankine(
        const BoundedVectorType& rPredictiveStressVector,
        const Vector& rStrainVector,
        double& rEquivalentStress,
        ConstitutiveLaw::Parameters& rValues)
    {
        array_1d<double, Dimension> principal_stress_vector = ZeroVector(Dimension);
        ConstitutiveLawUtilities<TVoigtSize>::CalculatePrincipalStresses(principal_stress_vector, rPredictiveStressVector);
        // The rEquivalentStress is the maximum principal stress
        if (Dimension == 3) 
            rEquivalentStress = std::max(std::max(principal_stress_vector[0], principal_stress_vector[1]), principal_stress_vector[2]);
        else // 2D
            rEquivalentStress = std::max(principal_stress_vector[0], principal_stress_vector[1]);
    }

    /**
     * @brief This method the uniaxial equivalent stress of Tresca
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rStrainVector The StrainVector vector
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateEquivalentStressTresca(
        const BoundedVectorType& rPredictiveStressVector,
        const Vector& rStrainVector,
        double& rEquivalentStress,
        ConstitutiveLaw::Parameters& rValues)
    {
        double I1, J2, J3, lode_angle;
        array_1d<double, TVoigtSize> deviator = ZeroVector(TVoigtSize);
        ConstitutiveLawUtilities<TVoigtSize>::CalculateI1Invariant(rPredictiveStressVector, I1);
        ConstitutiveLawUtilities<TVoigtSize>::CalculateJ2Invariant(rPredictiveStressVector, I1, deviator, J2);
        ConstitutiveLawUtilities<TVoigtSize>::CalculateJ3Invariant(deviator, J3);
        ConstitutiveLawUtilities<TVoigtSize>::CalculateLodeAngle(J2, J3, lode_angle);
        rEquivalentStress = 2.0 * std::cos(lode_angle) * std::sqrt(J2);
    }

    /**
     * @brief This method the uniaxial equivalent stress of SimoJu
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rStrainVector The StrainVector vector
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateEquivalentStressSimoJu(
        const BoundedVectorType& rPredictiveStressVector,
        const Vector& rStrainVector,
        double& rEquivalentStress,
        ConstitutiveLaw::Parameters& rValues)
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();
        // It compares with fc / sqrt(E)
        array_1d<double, Dimension> principal_stress_vector;
        ConstitutiveLawUtilities<VoigtSize>::CalculatePrincipalStresses(principal_stress_vector, rPredictiveStressVector);
        const double yield_compression = r_material_properties[YIELD_STRESS_C];
        const double yield_tension = r_material_properties[YIELD_STRESS_T];
        const double n = std::abs(yield_compression / yield_tension);

        double sum_a = 0.0, sum_b = 0.0, sum_c = 0.0, ere0, ere1;
        for (std::size_t cont = 0; cont < 2; ++cont) {
            sum_a += std::abs(principal_stress_vector[cont]);
            sum_b += 0.5 * (principal_stress_vector[cont] + std::abs(principal_stress_vector[cont]));
            sum_c += 0.5 * (-principal_stress_vector[cont] + std::abs(principal_stress_vector[cont]));
        }
        ere0 = sum_b / sum_a;
        ere1 = sum_c / sum_a;

        double auxf = 0.0;
        for (std::size_t cont = 0; cont < VoigtSize; ++cont) {
            auxf += rStrainVector[cont] * rPredictiveStressVector[cont]; // E:S
        }
        rEquivalentStress = std::sqrt(auxf);
        rEquivalentStress *= (ere0 * n + ere1);
    }

    /**
     * @brief This method the uniaxial equivalent stress of Drucker-Prager
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rStrainVector The StrainVector vector
     * @param rValues Parameters of the constitutive law
     */
    static void CalculateEquivalentStressDruckerPrager(
        const BoundedVectorType& rPredictiveStressVector,
        const Vector& rStrainVector,
        double& rEquivalentStress,
        ConstitutiveLaw::Parameters& rValues)
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();

        double friction_angle = r_material_properties[INTERNAL_FRICTION_ANGLE] * Globals::Pi / 180.0; // In radians!
        const double sin_phi = std::sin(friction_angle);
        const double root_3 = std::sqrt(3.0);

        // Check input variables
        if (friction_angle < tolerance) {
            friction_angle = 32.0 * Globals::Pi / 180.0;
            KRATOS_WARNING("DruckerPragerYieldSurface") << "Friction Angle not defined, assumed equal to 32 " << std::endl;
        }

        double I1, J2;
        ConstitutiveLawUtilities<TVoigtSize>::CalculateI1Invariant(rPredictiveStressVector, I1);
        array_1d<double, TVoigtSize> deviator = ZeroVector(TVoigtSize);
        ConstitutiveLawUtilities<TVoigtSize>::CalculateJ2Invariant(rPredictiveStressVector, I1, deviator, J2);

        if (std::abs(I1) < tolerance) {
            rEquivalentStress = 0.0;
        } else {
            const double CFL = -root_3 * (3.0 - sin_phi) / (3.0 * sin_phi - 3.0);
            const double TEN0 = 2.0 * I1 * sin_phi / (root_3 * (3.0 - sin_phi)) + std::sqrt(J2);
            rEquivalentStress = std::abs(CFL * TEN0);
        } 
    }

    /**
     * @brief This method returns the initial uniaxial stress threshold
     * for VonMises
     * @param rThreshold The uniaxial stress threshold
     * @param rValues Parameters of the constitutive law
     */
    static void GetInitialUniaxialThresholdHuberVonMises(
        ConstitutiveLaw::Parameters& rValues,
        double& rThreshold)
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();
        const double yield_tension = r_material_properties[YIELD_STRESS_T];
        rThreshold = std::abs(yield_tension);
    }

    /**
     * @brief This method returns the initial uniaxial stress threshold
     * for ModifiedMohrCoulomb
     * @param rThreshold The uniaxial stress threshold
     * @param rValues Parameters of the constitutive law
     */
    static void GetInitialUniaxialThresholdModifiedMohrCoulomb(
        ConstitutiveLaw::Parameters& rValues,
        double& rThreshold)
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();
        const double yield_compression = r_material_properties[YIELD_STRESS_C];
        rThreshold = std::abs(yield_compression);
    }

    /**
     * @brief This method returns the initial uniaxial stress threshold
     * for MohrCoulomb
     * @param rThreshold The uniaxial stress threshold
     * @param rValues Parameters of the constitutive law
     */
    static void GetInitialUniaxialThresholdMohrCoulomb(
        ConstitutiveLaw::Parameters& rValues,
        double& rThreshold)
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();
        const double cohesion = r_material_properties[COHESION];
        const double friction_angle = r_material_properties[INTERNAL_FRICTION_ANGLE] * Globals::Pi / 180.0;
        rThreshold = cohesion * std::cos(friction_angle);
    }

    /**
     * @brief This method returns the initial uniaxial stress threshold
     * for Rankine
     * @param rThreshold The uniaxial stress threshold
     * @param rValues Parameters of the constitutive law
     */
    static void GetInitialUniaxialThresholdRankine(
        ConstitutiveLaw::Parameters& rValues,
        double& rThreshold)
    {
        GetInitialUniaxialThresholdHuberVonMises(rValues, rThreshold);
    }

    /**
     * @brief This method returns the initial uniaxial stress threshold
     * for Tresca
     * @param rThreshold The uniaxial stress threshold
     * @param rValues Parameters of the constitutive law
     */
    static void GetInitialUniaxialThresholdTresca(
        ConstitutiveLaw::Parameters& rValues,
        double& rThreshold)
    {
        GetInitialUniaxialThresholdHuberVonMises(rValues, rThreshold);
    }

    /**
     * @brief This method returns the initial uniaxial stress threshold
     * for SimoJu
     * @param rThreshold The uniaxial stress threshold
     * @param rValues Parameters of the constitutive law
     */
    static void GetInitialUniaxialThresholdSimoJu(
        ConstitutiveLaw::Parameters& rValues,
        double& rThreshold)
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();
        const double yield_compression = r_material_properties[YIELD_STRESS_C];
        rThreshold = std::abs(yield_compression / std::sqrt(r_material_properties[YOUNG_MODULUS]));
    }

    /**
     * @brief This method returns the initial uniaxial stress threshold
     * for SimoJu
     * @param rThreshold The uniaxial stress threshold
     * @param rValues Parameters of the constitutive law
     */
    static void GetInitialUniaxialThresholdDruckerPrager(
        ConstitutiveLaw::Parameters& rValues,
        double& rThreshold)
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();
        const double yield_tension = r_material_properties[YIELD_STRESS_T];
        const double friction_angle = r_material_properties[INTERNAL_FRICTION_ANGLE] * Globals::Pi / 180.0; // In radians!
        const double sin_phi = std::sin(friction_angle);
        rThreshold = std::abs(yield_tension * (3.0 + sin_phi) / (3.0 * sin_phi - 3.0));
    }

    /**
     * @brief This method returns the damage parameter needed in the exp/linear expressions of damage
     * @param rAParameter The damage parameter
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateDamageParameterHuberVonMises(
        ConstitutiveLaw::Parameters& rValues,
        double& rAParameter,
        const double CharacteristicLength)
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();
        const double fracture_energy = r_material_properties[FRAC_ENERGY_T];
        const double young_modulus = r_material_properties[YOUNG_MODULUS];
        const double yield_tension = r_material_properties[YIELD_STRESS_T];
        rAParameter = 1.00 / (fracture_energy * young_modulus / (CharacteristicLength * std::pow(yield_tension, 2)) - 0.5);
        KRATOS_ERROR_IF(rAParameter < 0.0) << "Fracture energy is too low, increase FRACTURE_ENERGY..." << std::endl;
    }

    /**
     * @brief This method returns the damage parameter needed in the exp/linear expressions of damage
     * @param rAParameter The damage parameter
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateDamageParameterModifiedMohrCoulomb(
        ConstitutiveLaw::Parameters& rValues,
        double& rAParameter,
        const double CharacteristicLength)
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();
        const double fracture_energy = r_material_properties[FRAC_ENERGY_T];
        const double young_modulus = r_material_properties[YOUNG_MODULUS];
        const double yield_compression = r_material_properties[YIELD_STRESS_C];
        const double yield_tension = r_material_properties[YIELD_STRESS_T];
        const double n = yield_compression / yield_tension;
        rAParameter = 1.00 / (fracture_energy * n * n * young_modulus / (CharacteristicLength * std::pow(yield_compression, 2)) - 0.5);
        KRATOS_ERROR_IF(rAParameter < 0.0) << "Fracture energy is too low, increase FRACTURE_ENERGY..." << std::endl;
    }

    /**
     * @brief This method returns the damage parameter needed in the exp/linear expressions of damage
     * @param rAParameter The damage parameter
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateDamageParameterMohrCoulomb(
        ConstitutiveLaw::Parameters& rValues,
        double& rAParameter,
        const double CharacteristicLength)
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();
        const double fracture_energy = r_material_properties[FRAC_ENERGY_T];
        const double young_modulus = r_material_properties[YOUNG_MODULUS];
        const double cohesion = r_material_properties[COHESION];
        rAParameter = 1.00 / (fracture_energy * young_modulus / (CharacteristicLength * std::pow(cohesion, 2)) - 0.5);
        KRATOS_ERROR_IF(rAParameter < 0.0) << "Fracture energy is too low, increase FRACTURE_ENERGY..." << std::endl;
    }

    /**
     * @brief This method returns the damage parameter needed in the exp/linear expressions of damage
     * @param rAParameter The damage parameter
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateDamageParameterRankine(
        ConstitutiveLaw::Parameters& rValues,
        double& rAParameter,
        const double CharacteristicLength)
    {
        CalculateDamageParameterHuberVonMises(rValues, rAParameter, CharacteristicLength);
    }

    /**
     * @brief This method returns the damage parameter needed in the exp/linear expressions of damage
     * @param rAParameter The damage parameter
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateDamageParameterTresca(
        ConstitutiveLaw::Parameters& rValues,
        double& rAParameter,
        const double CharacteristicLength)
    {
        CalculateDamageParameterHuberVonMises(rValues, rAParameter, CharacteristicLength);
    }

    /**
     * @brief This method returns the damage parameter needed in the exp/linear expressions of damage
     * @param rAParameter The damage parameter
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateDamageParameterSimoJu(
        ConstitutiveLaw::Parameters& rValues,
        double& rAParameter,
        const double CharacteristicLength)
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();
        const double fracture_energy = r_material_properties[FRAC_ENERGY_T];
        const double yield_compression = r_material_properties[YIELD_STRESS_C];
        const double yield_tension = r_material_properties[YIELD_STRESS_T];
        const double n = yield_compression / yield_tension;
        rAParameter = 1.0 / (fracture_energy * n * n / (CharacteristicLength * std::pow(yield_compression, 2)) - 0.5);
        KRATOS_ERROR_IF(rAParameter < 0.0) << "Fracture energy is too low, increase FRACTURE_ENERGY..." << std::endl;
    }

    /**
     * @brief This method returns the damage parameter needed in the exp/linear expressions of damage
     * @param rAParameter The damage parameter
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateDamageParameterDruckerPrager(
        ConstitutiveLaw::Parameters& rValues,
        double& rAParameter,
        const double CharacteristicLength)
    {
        CalculateDamageParameterModifiedMohrCoulomb(rValues, rAParameter, CharacteristicLength); 
    }

    /**
     * @brief This method returns max value over a vector
     * @param rValues The values
     */
    static double GetMaxValue(const Vector& rValues)
    {
        double aux = 0.0;
        for (IndexType i = 0; i < rValues.size(); ++i) {
            if (aux < rValues[i]) aux = rValues[i];
        }
        return aux;
    }

    /**
     * @brief This method returns max abs value over a vector
     * @param rValues The values
     */
    static double GetMaxAbsValue(const Vector& rArrayValues)
    {
        const SizeType dimension = rArrayValues.size();

        IndexType counter = 0;
        double aux = 0.0;
        for (IndexType i = 0; i < dimension; ++i) {
            if (std::abs(rArrayValues[i]) > aux) {
                aux = std::abs(rArrayValues[i]);
                ++counter;
            }
        }
        return aux;
    }

    /**
     * @brief This method returns min abs value over a vector
     * @param rValues The values
     */
    static double GetMinAbsValue(const Vector& rArrayValues)
    {
        const SizeType dimension = rArrayValues.size();
        IndexType counter = 0;
        double aux = std::numeric_limits<double>::max();
        for (IndexType i = 0; i < dimension; ++i) {
            if (std::abs(rArrayValues[i]) < aux) {
                aux = std::abs(rArrayValues[i]);
                ++counter;
            }
        }
        return aux;
    }

    /**
     * @brief This method returns the 2 max values of a vector
     * @param rValues The values
     */
    static void Get2MaxValues(
        Vector& rMaxValues, 
        const double a, 
        const double b, 
        const double c
        )
    {
        rMaxValues.resize(2);
        Vector V;
        V.resize(3);
        V[0] = a;
        V[1] = b;
        V[2] = c;
        const int n = 3;

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n - 1; j++) {
                if (V[j] > V[j + 1]) {
                    double aux = V[j];
                    V[j] = V[j + 1];
                    V[j + 1] = aux;
                }
            }
        }
        rMaxValues[0] = V[2];
        rMaxValues[1] = V[1];
    }

    /**
     * @brief Calculation of the Hencky strain vector (true strain, natural strain, logarithmic strain)
     * @details See https://en.wikipedia.org/wiki/Finite_strain_theory#Seth%E2%80%93Hill_family_of_generalized_strain_tensors
     * @param rCauchyTensor The right Cauchy tensor
     * @param rStrainVector The Hencky strain vector
     */
    static void CalculateHenckyStrain(
        const MatrixType& rCauchyTensor,
        VectorType& rStrainVector);

    /**
     * @brief Calculation of the Biot strain vector
     * @details See https://en.wikipedia.org/wiki/Finite_strain_theory#Seth%E2%80%93Hill_family_of_generalized_strain_tensors
     * @param rCauchyTensor The right Cauchy tensor
     * @param rStrainVector The Biot strain vector
     */
    static void CalculateBiotStrain(
        const MatrixType& rCauchyTensor,
        VectorType& rStrainVector);

    /**
     * @brief Calculation of the Almansi strain vector
     * @details See https://en.wikipedia.org/wiki/Finite_strain_theory#Seth%E2%80%93Hill_family_of_generalized_strain_tensors
     * @param rLeftCauchyTensor The left Cauchy tensor
     * @param rStrainVector The Almansi strain vector
     */
    static void CalculateAlmansiStrain(
        const MatrixType& rLeftCauchyTensor,
        VectorType& rStrainVector);

    /**
     * @brief Calculation of the Green-Lagrange strain vector
     * @details See https://en.wikipedia.org/wiki/Finite_strain_theory#Seth%E2%80%93Hill_family_of_generalized_strain_tensors
     * @param rCauchyTensor The right Cauchy tensor
     * @param rStrainVector The Green-Lagrange strain vector
     */
    static void CalculateGreenLagrangianStrain(
        const MatrixType& rCauchyTensor,
        VectorType& rStrainVector);
  private:

}; // class ConstitutiveLawUtilities
} // namespace Kratos
#endif /* KRATOS_CONSTITUTIVE_LAW_UTILITIES defined */
