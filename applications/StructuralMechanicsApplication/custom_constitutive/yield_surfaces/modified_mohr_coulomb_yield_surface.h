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

#if !defined(KRATOS_MODIFIED_MOHR_COULOMB_YIELD_SURFACE_H_INCLUDED)
#define  KRATOS_MODIFIED_MOHR_COULOMB_YIELD_SURFACE_H_INCLUDED

// System includes

// Project includes
#include "custom_constitutive/yield_surfaces/generic_yield_surface.h"

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
 * @class ModifiedMohrCoulombYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @tparam TPlasticPotentialType The plastic potential considered
 * @tparam TVoigtSize The number of components on the Voigt notation
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <class TPlasticPotentialType>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ModifiedMohrCoulombYieldSurface
{
public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    typedef TPlasticPotentialType PlasticPotentialType;

    /// Counted pointer of ModifiedMohrCoulombYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION( ModifiedMohrCoulombYieldSurface);

    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    ModifiedMohrCoulombYieldSurface()
    {
    }

    /// Copy constructor
    ModifiedMohrCoulombYieldSurface(ModifiedMohrCoulombYieldSurface const& rOther)
    {
    }

    /// Assignment operator
    ModifiedMohrCoulombYieldSurface& operator=(ModifiedMohrCoulombYieldSurface const& rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~ModifiedMohrCoulombYieldSurface() {};

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method the uniaxial equivalent stress
     * @param StressVector The stress vector 
     * @param StrainVector The StrainVector vector
     * @param rMaterialProperties The material properties
     */
    static void CalculateEquivalentStress(  
        const Vector& StressVector,
        const Vector& StrainVector, 
        double& rEqStress, 
        const Properties& rMaterialProperties
    )
    {
        const double YieldCompression = rMaterialProperties[YIELD_STRESS_COMPRESSION];
        const double YieldTension = rMaterialProperties[YIELD_STRESS_TENSION];
        double FrictionAngle = rMaterialProperties[INTERNAL_FRICTION_ANGLE] * Globals::Pi / 180.0; // In radians!
		
        // Check input variables
        if (FrictionAngle < tolerance) {
            FrictionAngle = 32.0 * Globals::Pi / 180.0;
            KRATOS_WARNING("ModifiedMohrCoulombYieldSurface") << "Friction Angle not defined, assumed equal to 32 deg " << std::endl;
        }
        KRATOS_ERROR_IF(YieldCompression < tolerance) << " ERROR: Yield stress in compression not defined, include YIELD_STRESS_COMPRESSION in .mdpa ";
        KRATOS_ERROR_IF(YieldTension < tolerance) << " ERROR: Yield stress in tension not defined, include YIELD_STRESS_TENSION in .mdpa ";

        double K1, K2, K3, Rmorh, R, alpha_r, theta;
        R = std::abs(YieldCompression / YieldTension);
        Rmorh = std::pow(std::tan((Globals::Pi / 4.0) + FrictionAngle / 2.0), 2);
        alpha_r = R / Rmorh;
        double sinphi = std::sin(FrictionAngle);
		
        double I1, J2, J3; 
        ConstitutiveLawUtilities::CalculateI1Invariant(StressVector, I1);
        Vector Deviator = ZeroVector(6);
        ConstitutiveLawUtilities::CalculateJ2Invariant(StressVector, I1, Deviator, J2);
        ConstitutiveLawUtilities::CalculateJ3Invariant(Deviator, J3);

        K1 = 0.5*(1.0 + alpha_r) - 0.5*(1.0 - alpha_r)*sinphi;
        K2 = 0.5*(1.0 + alpha_r) - 0.5*(1.0 - alpha_r) / sinphi;
        K3 = 0.5*(1.0 + alpha_r)*sinphi - 0.5*(1.0 - alpha_r);

        // Check Modified Mohr-Coulomb criterion
        if (I1 == 0.0) {
            rEqStress = 0.0;
        } else {
            ConstitutiveLawUtilities::CalculateLodeAngle(J2, J3, theta);
            rEqStress = (2.0*std::tan(Globals::Pi*0.25 + FrictionAngle*0.5) / std::cos(FrictionAngle))*((I1*K3 / 3.0) +
                        std::sqrt(J2)*(K1*std::cos(theta) - K2*std::sin(theta)*sinphi / std::sqrt(3.0)));
        }
    }

    /**
     * @brief This method returns the initial uniaxial stress threshold
     * @param rThreshold The uniaxial stress threshold
     * @param rMaterialProperties The material properties
     */
    static void GetInitialUniaxialThreshold(const Properties& rMaterialProperties, double& rThreshold)
    {
        rThreshold = std::abs(rMaterialProperties[YIELD_STRESS_COMPRESSION]);
    }

    /**
     * @brief This method calculates the derivative of the plastic potential DG/DS
     * @param StressVector The stress vector 
     * @param Deviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the Deviator 
     * @param rg The derivative of the plastic potential
     * @param rMaterialProperties The material properties
     */
    static void CalculatePlasticPotentialDerivative(
        const Vector& StressVector,
        const Vector& Deviator,
        const double J2, 
        Vector& GFlux,
        const Properties& rMaterialProperties
    )
    {
        TPlasticPotentialType::CalculatePlasticPotentialDerivative(StressVector, Deviator, J2, GFlux, rMaterialProperties);
    }

    /**
     * @brief This method returns the damage parameter needed in the exp/linear expressions of damage
     * @param AParameter The damage parameter
     * @param rMaterialProperties The material properties
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateDamageParameter(
        const Properties& rMaterialProperties, 
        double& AParameter, 
        const double CharacteristicLength
        )
    {
        const double Gf = rMaterialProperties[FRACTURE_ENERGY];
        const double E  = rMaterialProperties[YOUNG_MODULUS];
        const double sigma_c = rMaterialProperties[YIELD_STRESS_COMPRESSION];
        const double sigma_t = rMaterialProperties[YIELD_STRESS_TENSION];
        const double n = sigma_c / sigma_t;
		
        if (rMaterialProperties[SOFTENING_TYPE] == static_cast<int>(SofteningType::Exponential)) {
            AParameter = 1.00 / (Gf*n*n*E / (CharacteristicLength * std::pow(sigma_c, 2)) - 0.5);
        } else { // linear
            AParameter = - std::pow(sigma_c, 2) / (2.0*E*Gf*n*n / CharacteristicLength);
        }
    }

    /**
     * @brief This  script  calculates  the derivatives  of the Yield Surf
    according   to   NAYAK-ZIENKIEWICZ   paper International
    journal for numerical methods in engineering vol 113-135 1972.
     As:            DF/DS = c1*V1 + c2*V2 + c3*V3
     * @param StressVector The stress vector 
     * @param Deviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the Deviator 
     * @param rFFlux The derivative of the yield surface
     * @param rMaterialProperties The material properties
     */
    static void CalculateYieldSurfaceDerivative(
        const Vector& StressVector, 
        const Vector& Deviator,
        const double J2, 
        Vector& rFFlux,
        const Properties& rMaterialProperties
    )
    {
        Vector FirstVector, SecondVector, ThirdVector;

        ConstitutiveLawUtilities::CalculateFirstVector(FirstVector);
        ConstitutiveLawUtilities::CalculateSecondVector(Deviator, J2, SecondVector);
        ConstitutiveLawUtilities::CalculateThirdVector(Deviator, J2, ThirdVector);

        double J3, LodeAngle;
        ConstitutiveLawUtilities::CalculateJ3Invariant(Deviator, J3);
        ConstitutiveLawUtilities::CalculateLodeAngle(J2, J3, LodeAngle);

        const double Checker = std::abs(LodeAngle*57.29577951308);

        double c1, c2, c3;
        const double FrictionAngle = rMaterialProperties[INTERNAL_FRICTION_ANGLE];
        const double SinPhi    = std::sin(FrictionAngle);
        const double CosPhi    = std::cos(FrictionAngle);
        const double SinTheta  = std::sin(LodeAngle);
        const double CosTheta  = std::cos(LodeAngle);
        const double Cos3Theta = std::cos(3.0*LodeAngle);
        const double TanTheta  = std::tan(LodeAngle);
        const double Tan3Theta = std::tan(3.0*LodeAngle);
        const double Root3     = std::sqrt(3.0);

        const double ComprYield = rMaterialProperties[YIELD_STRESS_COMPRESSION];
        const double TensiYield = rMaterialProperties[YIELD_STRESS_TENSION];
        const double n = ComprYield / TensiYield;

        const double Dilatancy = rMaterialProperties[DILATANCY_ANGLE];
        const double AnglePhi = (Globals::Pi * 0.25) + Dilatancy * 0.5;
        const double alpha = n / (std::tan(AnglePhi) * std::tan(AnglePhi));

        const double CFL = 2.0 * std::tan(AnglePhi) / CosPhi;

        const double K1 = 0.5*(1.0 + alpha) - 0.5*(1.0 - alpha)*SinPhi;
        const double K2 = 0.5*(1.0 + alpha) - 0.5*(1.0 - alpha) / SinPhi;
        const double K3 = 0.5*(1.0 + alpha)*SinPhi - 0.5*(1.0 - alpha);

        if (SinPhi != 0.0) c1 = CFL * K3 / 3.0;
        else c1 = 0.0; // check

        if (Checker < 29.0) {
            c2 = CosTheta * CFL * (K1*(1+TanTheta*Tan3Theta) + K2*SinPhi*(Tan3Theta-TanTheta) / Root3);
            c3 = CFL*(K1*Root3*SinTheta + K2*SinPhi*CosTheta) / (2.0*J2*Cos3Theta);
        } else {
            c3 = 0.0;
            double Aux = 1.0;
            if (LodeAngle > 0.0) Aux = -1.0;
            c2 = 0.5*CFL*(K1*Root3 + Aux*K2*SinPhi/Root3);
        }

        noalias(rFFlux) = c1*FirstVector + c2*SecondVector + c3*ThirdVector;
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

    void save(Serializer& rSerializer) const
    {
    }

    void load(Serializer& rSerializer)
    {
    }

    ///@}

}; // Class ModifiedMohrCoulombYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.
#endif
