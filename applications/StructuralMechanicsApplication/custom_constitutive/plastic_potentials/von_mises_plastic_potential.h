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


#if !defined(KRATOS_VON_MISES_PLASTIC_POTENTIAL_H_INCLUDED)
#define  KRATOS_VON_MISES_PLASTIC_POTENTIAL_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "utilities/math_utils.h"

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
 * @class VonMisesPlasticPotential
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @author Alejandro Cornejo & Lucia Barbu
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) VonMisesPlasticPotential
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of VonMisesPlasticPotential
    KRATOS_CLASS_POINTER_DEFINITION(VonMisesPlasticPotential);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    VonMisesPlasticPotential()
    {
    }

    /// Copy constructor
    VonMisesPlasticPotential(VonMisesPlasticPotential const& rOther)
    {
    }

    /// Assignment operator
    VonMisesPlasticPotential& operator=(VonMisesPlasticPotential const& rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~VonMisesPlasticPotential() {};


    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /*
    This script calculates the derivatives of the yield
    functions according to NAYAK-ZIENKIEWICZ paper. International
    journal for numerical methods in engineering vol 113-135 1972.
    As:            DF/DS = c1*V1 + c2*V2 + c3*V3
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

        CalculateFirstVector(FirstVector);
        CalculateSecondVector(Deviator, J2, SecondVector);
        CalculateThirdVector(Deviator, J2, ThirdVector);

        double c1, c2, c3;
        c1 = 0.0;
        c2 = std::sqrt(3.0);
        c3 = 0.0;

        noalias(rGFlux) = c1*FirstVector + c2*SecondVector + c3*ThirdVector;
    }

    static void CalculateFirstVector(Vector& FirstVector)
    {
        FirstVector = ZeroVector(6);
        FirstVector[0] = 1.0;
        FirstVector[1] = 1.0;
        FirstVector[2] = 1.0;

    }

    static void CalculateSecondVector(
        const Vector Deviator, 
        const double J2, 
        Vector& SecondVector
    )
    {
        const double twosqrtJ2 = 2.0*std::sqrt(J2);
        for (int i = 0; i < 6; i++)
        {
            SecondVector[i] = Deviator[i] / (twosqrtJ2);
        }

        SecondVector[3] *= 2.0;
        SecondVector[4] *= 2.0;
        SecondVector[5] *= 2.0;
    }

    static void CalculateThirdVector(
        const Vector Deviator, 
        const double J2, 
        Vector& ThirdVector
    )
    {
        ThirdVector.resize(6);
        const double J2thirds = J2 / 3.0;

        ThirdVector[0] = Deviator[1]*Deviator[2] - Deviator[4]*Deviator[4] + J2thirds;
        ThirdVector[1] = Deviator[0]*Deviator[2] - Deviator[5]*Deviator[5] + J2thirds;
        ThirdVector[2] = Deviator[0]*Deviator[1] - Deviator[3]*Deviator[3] + J2thirds;
        ThirdVector[3] = 2.0*(Deviator[4]*Deviator[5] - Deviator[3]*Deviator[2]);
        ThirdVector[4] = 2.0*(Deviator[3]*Deviator[4] - Deviator[1]*Deviator[5]);
        ThirdVector[5] = 2.0*(Deviator[5]*Deviator[3] - Deviator[0]*Deviator[4]);
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
