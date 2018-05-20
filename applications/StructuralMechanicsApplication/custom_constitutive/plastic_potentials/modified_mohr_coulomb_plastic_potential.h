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


#if !defined(KRATOS_MODIFIED_MOHR_COULOMB_PLASTIC_POTENTIAL_H_INCLUDED)
#define  KRATOS_MODIFIED_MOHR_COULOMB_PLASTIC_POTENTIAL_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "utilities/math_utils.h"
#include "includes/global_variables.h"


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
 * @class ModifiedMohrCoulombPlasticPotential
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @author Alejandro Cornejo
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ModifiedMohrCoulombPlasticPotential
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of ModifiedMohrCoulombPlasticPotential
    KRATOS_CLASS_POINTER_DEFINITION(ModifiedMohrCoulombPlasticPotential);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Initialization constructor.
    ModifiedMohrCoulombPlasticPotential()
    {
    }

    /// Copy constructor
    ModifiedMohrCoulombPlasticPotential(ModifiedMohrCoulombPlasticPotential const& rOther)
    {
    }

    /// Assignment operator
    ModifiedMohrCoulombPlasticPotential& operator=(ModifiedMohrCoulombPlasticPotential const& rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~ModifiedMohrCoulombPlasticPotential() {};


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

        CalculateFirstVector(FirstVector);
        CalculateSecondVector(Deviator, J2, SecondVector);
        CalculateThirdVector(Deviator, J2, ThirdVector);

        double J3, LodeAngle;
        CalculateJ3Invariant(Deviator, J3);
        CalculateLodeAngle(J2, J3, LodeAngle);

        const double Checker = std::abs(LodeAngle*57.29577951308);

        double c1, c2, c3;
        const double Dilatancy = rMaterialProperties[DILATANCY_ANGLE];
        const double SinDil    = std::sin(Dilatancy);
        const double CosDil    = std::cos(Dilatancy);
        const double SinTheta  = std::sin(LodeAngle);
        const double CosTheta  = std::cos(LodeAngle);
        const double Cos3Theta = std::cos(3.0*LodeAngle);
        const double TanTheta  = std::tan(LodeAngle);
        const double Tan3Theta = std::tan(3.0*LodeAngle);
        const double Root3     = std::sqrt(3.0);

        const double ComprYield = rMaterialProperties[YIELD_STRESS_C];
		const double TensiYield = rMaterialProperties[YIELD_STRESS_T];
        const double n = ComprYield / TensiYield;

        const double AnglePhi = (Globals::Pi * 0.25) + Dilatancy * 0.5;
        const double alpha = n / (std::tan(AnglePhi) * std::tan(AnglePhi));

        const double CFL = 2.0 * std::tan(AnglePhi) / CosDil;

		const double K1 = 0.5*(1 + alpha) - 0.5*(1 - alpha)*SinDil;
		const double K2 = 0.5*(1 + alpha) - 0.5*(1 - alpha) / SinDil;
		const double K3 = 0.5*(1 + alpha)*SinDil - 0.5*(1 - alpha);

        if (SinDil != 0.0) c1 = CFL * K3 / 3.0;
        else c1 = 0.0; // check

        if (Checker < 29.0)
        {
            c2 = CosTheta * CFL * (K1*(1+TanTheta*Tan3Theta) + K2*SinDil*(Tan3Theta-TanTheta) / Root3);
            c3 = CFL*(K1*Root3*SinTheta + K2*SinDil*CosTheta) / (2.0*J2*Cos3Theta);
        }
        else
        {
            c3 = 0.0;
            double Aux = 1.0;
            if (LodeAngle > 0.0) Aux = -1.0;
            c2 = 0.5*CFL*(K1*Root3 + Aux*K2*SinDil/Root3);
        }

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

    static void CalculateLodeAngle(const double J2, const double J3, double& LodeAngle)
    {
		double sint3 = (-3.0*std::sqrt(3.0)*J3) / (2.0*J2*std::sqrt(J2));
		if (sint3 < -0.95) sint3 = -1;
		if (sint3 > 0.95)  sint3 = 1; 
		LodeAngle = std::asin(sint3) / 3.0;
    }

    static void CalculateJ3Invariant(const Vector& Deviator, double& rJ3)
    {
        rJ3 = Deviator[0]*(Deviator[1]*Deviator[2] - Deviator[4]*Deviator[4])  +
			Deviator[3]*(-Deviator[3]*Deviator[2]  + Deviator[5]*Deviator[4])  +
			Deviator[5]*(Deviator[3]*Deviator[4] - Deviator[5]*Deviator[1]);
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
