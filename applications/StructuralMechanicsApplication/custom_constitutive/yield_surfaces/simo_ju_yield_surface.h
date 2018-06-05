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

#if !defined(KRATOS_SIMO_JU_YIELD_SURFACE_H_INCLUDED)
#define  KRATOS_SIMO_JU_YIELD_SURFACE_H_INCLUDED

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
 * @class SimoJuYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @tparam TPlasticPotentialType The plastic potential considered
 * @tparam TVoigtSize The number of components on the Voigt notation
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <class TPlasticPotentialType , std::size_t TVoigtSize>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SimoJuYieldSurface
{
public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    typedef TPlasticPotentialType PlasticPotentialType;

    /// Counted pointer of SimoJuYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(SimoJuYieldSurface);

    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    SimoJuYieldSurface()
    {
    }

    /// Copy constructor
    SimoJuYieldSurface(SimoJuYieldSurface const& rOther)
    {
    }

    /// Assignment operator
    SimoJuYieldSurface& operator=(SimoJuYieldSurface const& rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~SimoJuYieldSurface() {};

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{

    static void CalculateEquivalentStress(  
        const Vector& StressVector,
        const Vector& StrainVector, 
        double& rEqStress, 
        const Properties& rMaterialProperties
    )
    {   // It compares with fc / sqrt(E)
        Vector PrincipalStressVector;
        ConstitutiveLawUtilities::CalculatePrincipalStresses(PrincipalStressVector, StressVector);

        double sigma_t, sigma_c, n;
        sigma_t = rMaterialProperties[YIELD_STRESS_TENSION];
        sigma_c = rMaterialProperties[YIELD_STRESS_COMPRESSION];
        n = std::abs(sigma_c / sigma_t);

        double SumA, SumB, SumC, ere0, ere1;
        for (std::size_t cont = 0;cont < 2;cont++) {
            SumA += std::abs(PrincipalStressVector[cont]);
            SumB += 0.5*(PrincipalStressVector[cont]  + std::abs(PrincipalStressVector[cont]));
            SumC += 0.5*(-PrincipalStressVector[cont] + std::abs(PrincipalStressVector[cont]));
        }
        ere0 = SumB / SumA;
        ere1 = SumC / SumA;

        // Check SimoJu criterion
        if (std::abs(StrainVector[0]) > tolerance) {
            rEqStress = 0;
        } else {
            double auxf = 0.0;
            for (std::size_t cont = 0; cont < 6; cont++) {
                auxf += StrainVector[cont] * StressVector[cont];  // E*S
            }
            rEqStress = std::sqrt(auxf);
            rEqStress *= (ere0*n + ere1);
        }
    }

    static void GetInitialUniaxialThreshold(const Properties& rMaterialProperties, double& rThreshold)
    {
        rThreshold = std::abs(rMaterialProperties[YIELD_STRESS_COMPRESSION] / std::sqrt(rMaterialProperties[YOUNG_MODULUS]));
    }

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

         if (rMaterialProperties[SOFTENING_TYPE] == static_cast<std::size_t>(SofteningType::Exponential)) {
            AParameter = 1.00 / (n*n*Gf*E / (CharacteristicLength * std::pow(sigma_c, 2)) - 0.5);
        } else {
            
        }
    }

    // Computes dG/dS
    static void CalculatePlasticPotentialDerivative(
        const Vector& StressVector,
        const Vector& Deviator,
        const double J2, 
        Vector& rg,
        const Properties& rMaterialProperties
    )
    {
        TPlasticPotentialType::CalculatePlasticPotentialDerivative(StressVector, Deviator, J2, rg, rMaterialProperties);
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

}; // Class SimoJuYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.
#endif
