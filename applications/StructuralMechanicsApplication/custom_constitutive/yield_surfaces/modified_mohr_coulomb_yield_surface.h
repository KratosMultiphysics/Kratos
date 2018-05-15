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

#if !defined(KRATOS_MODIFIED_MOHR_COULOMB_YIELD_SURFACE_H_INCLUDED)
#define  KRATOS_MODIFIED_MOHR_COULOMB_YIELD_SURFACE_H_INCLUDED

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
 * @class ModifiedMohrCoulombYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @tparam TPlasticPotentialType 
 * @author Alejandro Cornejo
 */
template <class TPlasticPotentialType , class TVoigtSize>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ModifiedMohrCoulombYieldSurface
{
public:
    ///@name Type Definitions
    ///@{

    /// The type of potential plasticity
    typedef typename TPlasticPotentialType PlasticPotentialType;

    /// Counted pointer of ModifiedMohrCoulombYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION( ModifiedMohrCoulombYieldSurface );

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

    static void CalculateEquivalentStress(  
        const Vector& StressVector, 
        double& rEqStress, 
        const Properties& rMaterialProperties
    )
    {       
        double sigma_c = 0.0, sigma_t = 0.0, friction_angle = 0.0, E = 0.0, Gt = 0.0;
		sigma_c = rMaterialProperties[YIELD_STRESS_C];
		sigma_t = rMaterialProperties[YIELD_STRESS_T];
		friction_angle = rMaterialProperties[INTERNAL_FRICTION_ANGLE] * 3.14159265359 / 180.0; // In radians!
		E  = rMaterialProperties[YOUNG_MODULUS];
		Gt = rMaterialProperties[FRAC_ENERGY_T];

		// Check input variables 
		if (friction_angle < 1e-24) { friction_angle = 32 * 3.14159 / 180; std::cout << "Friction Angle not defined, assumed equal to 32 deg " << std::endl; }
		if (sigma_c < 1e-24) { KRATOS_ERROR << " ERROR: Yield stress in compression not defined, include YIELD_STRESS_C in .mdpa "; }
		if (sigma_t < 1e-24) { KRATOS_ERROR << " ERROR: Yield stress in tension not defined, include YIELD_STRESS_T in .mdpa "; }
		if (Gt < 1e-24) { KRATOS_ERROR << " ERROR: Fracture Energy not defined in the model part, include FRAC_ENERGY_T in .mdpa "; }

		double K1, K2, K3, Rmorh, R, alpha_r, c_max, theta, c_threshold;
		R = std::abs(sigma_c / sigma_t);
		Rmorh = pow(tan((3.14159265359 / 4) + friction_angle / 2), 2);
		alpha_r = R / Rmorh;
		c_max = std::abs(sigma_c);
		double sinphi = sin(friction_angle);

		double I1 = 0.0, J2 = 0.0, J3 = 0.0;
        CalculateI1Invariant(StressVector, I1);
		Vector Deviator = ZeroVector(6);
        CalculateJ2Invariant(StressVector, I1, Deviator, J2);
		CalculateJ3Invariant(Deviator, J3);
		K1 = 0.5*(1 + alpha_r) - 0.5*(1 - alpha_r)*sinphi;
		K2 = 0.5*(1 + alpha_r) - 0.5*(1 - alpha_r) / sinphi;
		K3 = 0.5*(1 + alpha_r)*sinphi - 0.5*(1 - alpha_r);

		double rEqStress = 0.0; 
		// Check Modified Mohr-Coulomb criterion
		if (I1 == 0.0)  rEqStress = 0.0; 
		else
		{
			theta = Calculate_Theta_Angle(J2, J3);
			rEqStress = (2.0*tan(3.14159265359*0.25 + friction_angle*0.5) / cos(friction_angle))*((I1*K3 / 3.0) + 
                sqrt(J2)*(K1*cos(theta) - K2*sin(theta)*sinphi / sqrt(3.0)));
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
        Vector& rg)
    {
        TPlasticPotentialType::CalculatePlasticPotentialDerivative(StressVector, Deviator, J2, rg);
    }

    static void CalculateLodeAngle(const double J2, const double J3, double& LodeAngle)
    {
        double sint3;
		sint3 = (-3.0*sqrt(3)*J3) / (2.0 * J2*sqrt(J2));
		if (sint3 < -0.95) { sint3 = -1;}
		if (sint3 > 0.95)  { sint3 = 1; }
		LodeAngle = asin(sint3) / 3.0;
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
