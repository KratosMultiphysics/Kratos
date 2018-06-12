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


#if !defined(KRATOS_DRUCKER_PRAGER_PLASTIC_POTENTIAL_H_INCLUDED)
#define  KRATOS_DRUCKER_PRAGER_PLASTIC_POTENTIAL_H_INCLUDED

// System includes

// Project includes
#include "custom_constitutive/plastic_potentials/generic_plastic_potential.h"

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
 * @class DruckerPragerPlasticPotential
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @author Alejandro Cornejo & Lucia Barbu
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) DruckerPragerPlasticPotential
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of DruckerPragerPlasticPotential
    KRATOS_CLASS_POINTER_DEFINITION(DruckerPragerPlasticPotential);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    DruckerPragerPlasticPotential()
    {
    }

    /// Copy constructor
    DruckerPragerPlasticPotential(DruckerPragerPlasticPotential const& rOther)
    {
    }

    /// Assignment operator
    DruckerPragerPlasticPotential& operator=(DruckerPragerPlasticPotential const& rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~DruckerPragerPlasticPotential() {};


    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /*
    This  script  calculates  the derivatives  of the potential
    functions according to NAYAK-ZIENKIEWICZ paper International
    journal for numerical methods in engineering vol 113-135 1972.
    As:            DG/DS = c1*V1 + c2*V2 + c3*V3
    */

    static void CalculatePlasticPotentialDerivative(
        const Vector& StressVector, 
        const Vector& Deviator,
        const double J2, 
        Vector& rGFlux,
        const Properties& rMaterialProperties
    )
    {
        Vector FirstVector, SecondVector, ThirdVector;

        ConstitutiveLawUtilities::CalculateFirstVector(FirstVector);
        ConstitutiveLawUtilities::CalculateSecondVector(Deviator, J2, SecondVector);
        ConstitutiveLawUtilities::CalculateThirdVector(Deviator, J2, ThirdVector);

        double c1, c2, c3;
        c3 = 0.0;

        const double Dilatancy = rMaterialProperties[DILATANCY_ANGLE] * Globals::Pi / 180.0;
        const double SinDil    = std::sin(Dilatancy);
        const double Root3     = std::sqrt(3.0);

        const double CFL = -Root3*(3.0-SinDil) / (3.0*SinDil-3.0);
        c1 = CFL*2.0*SinDil / (Root3*(3.0-SinDil));
        c2 = CFL;

        noalias(rGFlux) = c1*FirstVector + c2*SecondVector + c3*ThirdVector;
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

}; // Class GenericYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.
#endif
