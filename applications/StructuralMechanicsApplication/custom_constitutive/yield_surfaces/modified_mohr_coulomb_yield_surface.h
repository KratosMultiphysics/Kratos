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
 * @class ModifiedMohrCoulombYieldSurface
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @tparam TPlasticPotentialType 
 * @author Alejandro Cornejo & Lucia Barbu
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
    template <class TDim>
    static void CalculateEquivalentStress(  
        const Vector& StressVector,
        const Vector& StrainVector, 
        double& rEqStress, 
        const Properties& rMaterialProperties
    )
    {
		const double YieldCompression = rMaterialProperties[YIELD_STRESS_COMPRESSION];
		const double YieldTension = rMaterialProperties[YIELD_STRESS_TENSION];
		const double FrictionAngle = rMaterialProperties[FRICTION_ANGLE] * Globals::Pi / 180.0; // In radians!

		// Check input variables 
        double tol = std::numeric_limits<double>::epsilon();
		if (FrictionAngle < tol) { FrictionAngle = 32 * Globals::Pi / 180; std::cout << "Friction Angle not defined, assumed equal to 32 deg " << std::endl; }
		if (YieldCompression < tol) { KRATOS_ERROR << " ERROR: Yield stress in compression not defined, include YIELD_STRESS_COMPRESSION in .mdpa "; }
		if (YieldTension < tol) { KRATOS_ERROR << " ERROR: Yield stress in tension not defined, include YIELD_STRESS_TENSION in .mdpa "; }

		double K1, K2, K3, Rmorh, R, alpha_r, theta;
		R = std::abs(YieldCompression / YieldTension);
		Rmorh = std::pow(tan((Globals::Pi / 4.0) + FrictionAngle / 2.0), 2);
		alpha_r = R / Rmorh;
		double sinphi = std::sin(FrictionAngle);

		double I1, J2, J3;
        CalculateI1Invariant(StressVector, I1);
		Vector Deviator = ZeroVector(6);
        CalculateJ2Invariant(StressVector, I1, Deviator, J2);
		CalculateJ3Invariant(Deviator, J3);

		K1 = 0.5*(1 + alpha_r) - 0.5*(1 - alpha_r)*sinphi;
		K2 = 0.5*(1 + alpha_r) - 0.5*(1 - alpha_r) / sinphi;
		K3 = 0.5*(1 + alpha_r)*sinphi - 0.5*(1 - alpha_r);

		double rEqStress; 
		// Check Modified Mohr-Coulomb criterion
		if (I1 == 0.0)  rEqStress = 0.0; 
		else
		{
			CalculateLodeAngle(J2, J3, theta);
			rEqStress = (2.0*std::tan(Globals::Pi*0.25 + FrictionAngle*0.5) / std::cos(FrictionAngle))*((I1*K3 / 3.0) + 
                std::sqrt(J2)*(K1*std::cos(theta) - K2*std::sin(theta)*sinphi / std::sqrt(3.0)));
		}
    }

    static void GetInitialUniaxialThreshold(const Properties& rMaterialProperties, double& rThreshold)
    {
        rThreshold = std::abs(rMaterialProperties[YIELD_STRESS_COMPRESSION]);
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
        rDeviator = StressVector;
        const double Pmean = I1 / 3.0;

        rDeviator[0] -= Pmean;
        rDeviator[1] -= Pmean;
        rDeviator[2] -= Pmean;

        rJ2 = 0.5*(rDeviator[0]*rDeviator[0] + rDeviator[1]*rDeviator[1] + rDeviator[2]*rDeviator[2]) +
            (rDeviator[3]*rDeviator[3] + rDeviator[4]*rDeviator[4] + rDeviator[5]*rDeviator[5]);
        
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
        Vector& GFlux,
        const Properties& rMaterialProperties
    )
    {
        TPlasticPotentialType::CalculatePlasticPotentialDerivative(StressVector, Deviator, J2, GFlux, rMaterialProperties);
    }

    static void CalculateLodeAngle(const double J2, const double J3, double& LodeAngle)
    {
		const double sint3 = (-3.0*std::sqrt(3.0)*J3) / (2.0*J2*std::sqrt(J2));
		if (sint3 < -0.95) sint3 = -1;
		if (sint3 > 0.95)  sint3 = 1; 
		LodeAngle = std::asin(sint3) / 3.0;
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

        if (rMaterialProperties[SOFTENING_TYPE] == "Exponential")
        {
            AParameter = 1.00 / (n*n*Gt*E / (CharacteristicLength * std::pow(sigma_c, 2)) - 0.5);
        }
        else
        {
            
        }
    }

    /*
    This  script  calculates  the derivatives  of the Yield Surf
    according   to   NAYAK-ZIENKIEWICZ   paper International
    journal for numerical methods in engineering vol 113-135 1972.
    As:            DF/DS = c1*V1 + c2*V2 + c3*V3
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

        CalculateFirstVector(FirstVector);
        CalculateSecondVector(Deviator, J2, SecondVector);
        CalculateThirdVector(Deviator, J2, ThirdVector);

        double J3, LodeAngle;
        CalculateJ3Invariant(Deviator, J3);
        CalculateLodeAngle(J2, J3, LodeAngle);

        const double Checker = std::abs(LodeAngle*57.29577951308);

        double c1, c2, c3;
        const double FrictionAngle = rMaterialProperties[FRICTION_ANGLE];
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

        const double AnglePhi = (Globals::Pi * 0.25) + Dilatancy * 0.5;
        const double alpha = n / (std::tan(AnglePhi) * std::tan(AnglePhi));

        const double CFL = 2.0 * std::tan(AnglePhi) / CosPhi;

		const double K1 = 0.5*(1 + alpha) - 0.5*(1 - alpha)*SinPhi;
		const double K2 = 0.5*(1 + alpha) - 0.5*(1 - alpha) / SinPhi;
		const double K3 = 0.5*(1 + alpha)*SinPhi - 0.5*(1 - alpha);

        if (SinPhi != 0.0) c1 = CFL * K3 / 3.0;
        else c1 = 0.0; // check

        if (Checker < 29.0)
        {
            c2 = CosTheta * CFL * (K1*(1+TanTheta*Tan3Theta) + K2*SinPhi*(Tan3Theta-TanTheta) / Root3);
            c3 = CFL*(K1*Root3*SinTheta + K2*SinPhi*CosTheta) / (2.0*J2*Cos3Theta);
        }
        else
        {
            c3 = 0.0;
            double Aux = 1.0;
            if (LodeAngle > 0.0) Aux = -1.0;
            c2 = 0.5*CFL*(K1*Root3 + Aux*K2*SinPhi/Root3);
        }

        noalias(rFFlux) = c1*FirstVector + c2*SecondVector + c3*ThirdVector;
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
    //template class CalculateI1Invariant<3>;
    //template class CalculateI2Invariant<3>;
    //template class CalculateI3Invariant<3>;
    //template class CalculateJ2Invariant<3>;
    //template class CalculateJ3Invariant<3>;

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
