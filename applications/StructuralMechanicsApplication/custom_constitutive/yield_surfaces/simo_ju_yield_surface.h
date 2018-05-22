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

#if !defined(KRATOS_SIMO_JU_YIELD_SURFACE_H_INCLUDED)
#define  KRATOS_SIMO_JU_YIELD_SURFACE_H_INCLUDED

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
 * @class SimoJuYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @tparam TPlasticPotentialType 
 * @author Alejandro Cornejo
 */
template <class TPlasticPotentialType , class TVoigtSize>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SimoJuYieldSurface
{
public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    typedef typename TPlasticPotentialType PlasticPotentialType;

    /// Counted pointer of SimoJuYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(SimoJuYieldSurface);

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
		CalculatePrincipalStresses(PrincipalStressVector, StressVector);

		double sigma_t, sigma_c, n;
		sigma_t = rMaterialProperties[YIELD_STRESS_T];
		sigma_c = rMaterialProperties[YIELD_STRESS_C];
		n = abs(sigma_c / sigma_t);

		double SumA, SumB, SumC, ere0, ere1;
		for (int cont = 0;cont < 2;cont++)
		{
			SumA += std::abs(PrincipalStressVector[cont]);
			SumB += 0.5*(PrincipalStressVector[cont]  + std::abs(PrincipalStressVector[cont]));
			SumC += 0.5*(-PrincipalStressVector[cont] + std::abs(PrincipalStressVector[cont]));
		}
		ere0 = SumB / SumA;
		ere1 = SumC / SumA;

		// Check SimoJu criterion
		if (StrainVector[0] == 0) rEqStress = 0; 
		else
		{
			double auxf = 0.0;
			for (int cont = 0; cont < 6; cont++)
			{
				auxf += StrainVector[cont] * StressVector[cont];  // E*S
			}
			rEqStress = std::sqrt(auxf);
			rEqStress *= (ere0*n + ere1);
		}
    }

    static void GetInitialUniaxialThreshold(const Properties& rMaterialProperties, double& rThreshold)
    {
        rThreshold = std::abs(rMaterialProperties[YIELD_STRESS_C] / std::sqrt(rMaterialProperties[YOUNG_MODULUS]));
    }

    static void CalculateDamageParameter(
        const Properties& rMaterialProperties, 
        double& AParameter, 
        const double CharacteristicLength
    )
    {
        const double Gf = rMaterialProperties[FRACTURE_ENERGY];
        const double E  = rMaterialProperties[YOUNG_MODULUS];
        const double sigma_c = rMaterialProperties[YIELD_STRESS_C];
        const double sigma_t = rMaterialProperties[YIELD_STRESS_T];
        const double n = sigma_c / sigma_t;

        if (rMaterialProperties[SOFTENING_TYPE] == Exponential)
        {
            AParameter = 1.00 / (n*n*Gt*E / (CharacteristicLength * std::pow(sigma_c, 2)) - 0.5);
        }
        else
        {
            
        }
    }

    static void CalculateI1Invariant(const Vector& StressVector, double& rI1)
    {
        rI1 = StressVector[0] + StressVector[1] + StressVector[2];
    }

    static void CalculateI2Invariant(const Vector& StressVector, double& rI2)
    {
        rI2 = (StressVector[0] + StressVector[2])*StressVector[1] + StressVector[0]*StressVector[2] +
            - StressVector[3]*StressVector[3] - StressVector[4]*StressVector[4] - StressVector[5]*StressVector[5];
    }

    static void CalculateI3Invariant(const Vector& StressVector, double& rI3)
    {
        rI3 = (StressVector[1]*StressVector[2] - StressVector[4]*StressVector[4])*StressVector[0] -
            StressVector[1]*StressVector[5]*StressVector[5] - StressVector[2]*StressVector[3]*StressVector[3] +
            2.0*StressVector[3]*StressVector[4]*StressVector[5];
    }

    static void CalculateJ2Invariant(const Vector& StressVector, const double& I1, Vector& rDeviator, double& rJ2)
    {
        if (TVoigtSize == 6)
        {
            rDeviator = StressVector;
            double Pmean = I1 / 3.0;

            rDeviator[0] -= Pmean;
            rDeviator[1] -= Pmean;
            rDeviator[2] -= Pmean;

            rJ2 = 0.5*(rDeviator[0]*rDeviator[0] + rDeviator[1]*rDeviator[1] + rDeviator[2]*rDeviator[2]) +
                (rDeviator[3]*rDeviator[3] + rDeviator[4]*rDeviator[4] + rDeviator[5]*rDeviator[5]);
        }
        else
        {
            // 2d
        }

    }

    static void CalculateJ3Invariant(const Vector& Deviator, double& rJ3)
    {
        rJ3 = Deviator[0]*(Deviator[1]*Deviator[2] - Deviator[4]*Deviator[4])  +
			Deviator[3]*(-Deviator[3]*Deviator[2]  + Deviator[5]*Deviator[4])  +
			Deviator[5]*(Deviator[3]*Deviator[4] - Deviator[5]*Deviator[1]);
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

    static void CalculateLodeAngle(const double J2, const double J3, double& LodeAngle)
    {
		const double sint3 = (-3.0*std::sqrt(3.0)*J3) / (2.0*J2*std::sqrt(J2));
		if (sint3 < -0.95) sint3 = -1;
		if (sint3 > 0.95)  sint3 = 1; 
		LodeAngle = asin(sint3) / 3.0;
    }

	static void CalculatePrincipalStresses(Vector& rPrincipalStressVector, const Vector StressVector)
	{
		rPrincipalStressVector.resize(3);
		double I1, I2, I3, phi, Num, Denom, II1;
		CalculateI1Invariant(StressVector, I1);
		CalculateI2Invariant(StressVector, I2);
		CalculateI3Invariant(StressVector, I3);
		II1 = I1*I1;

		Num = (2.0*II1 - 9.0*I2)*I1 + 27.0*I3;
		Denom = (II1 - 3.0*I2);

		if (Denom != 0.0)
		{
			phi = Num / (2.0*Denom*std::sqrt(Denom));

			if (std::abs(phi) > 1.0)
			{
				if (phi > 0.0) phi = 1.0;
				else phi = -1.0;
			}

			double acosphi = std::acos(phi);
			phi = acosphi / 3.0;

			double aux1 = 0.666666666666667*sqrt(II1 - 3.0*I2);
			double aux2 = I1 / 3.0;

			rPrincipalStressVector[0] = aux2 + aux1*cos(phi);
			rPrincipalStressVector[1] = aux2 + aux1*cos(phi - 2.09439510239);
			rPrincipalStressVector[2] = aux2 + aux1*cos(phi - 4.18879020478);
		}
		else 
		{
			rPrincipalStressVector = ZeroVector(3);
		}
		
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
