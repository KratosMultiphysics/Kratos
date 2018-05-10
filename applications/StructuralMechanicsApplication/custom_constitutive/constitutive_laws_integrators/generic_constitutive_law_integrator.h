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

#if !defined(KRATOS_GENERIC_CONSTITUTIVE_LAW_INTEGRATOR_H_INCLUDED)
#define  KRATOS_GENERIC_CONSTITUTIVE_LAW_INTEGRATOR_H_INCLUDED

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
 * @class GenericConstitutiveLawIntegrator
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @tparam TYieldSurfaceType
 * @author Alejandro Cornejo
 */
template <class TYieldSurfaceType, class TVoigtSize>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) GenericConstitutiveLawIntegrator
{
public:
    ///@name Type Definitions
    ///@{

    struct PlasticParameters
    {
        Vector PredictiveStressVector; 
        double UniaxialStress;
        double PlasticDenominator; 
        Vector Fflux;
        Vector Gflux;
        double PlasticDissipation;
        Vector PlasticStrainIncrement;
        Vector PlasticStrain;

        void Initialize()
        {

        }
	};

    /// The type of yield surface
    typedef typename TYieldSurfaceType YieldSurfaceType;

    /// The type of plastic potential 
    typedef typename YieldSurfaceType::TPlasticPotentialType PlasticPotentialType;

    /// Counted pointer of GenericConstitutiveLawIntegrator
    KRATOS_CLASS_POINTER_DEFINITION(GenericConstitutiveLawIntegrator);

    /// Initialization constructor
    GenericConstitutiveLawIntegrator()
    {
    }

    /// Copy constructor
    GenericConstitutiveLawIntegrator(GenericConstitutiveLawIntegrator const& rOther)
    {
    }

    /// Assignment operator
    GenericConstitutiveLawIntegrator& operator=(GenericConstitutiveLawIntegrator const& rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~GenericConstitutiveLawIntegrator()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    static void IntegrateStressVector(
        Vector& PredictiveStressVector, 
        double& UniaxialStress, 
        double& Threshold, 
        double& PlasticDenominator, 
        Vector& Fflux, Vector& Gflux, 
        double& PlasticDissipation, 
        Vector& PlasticStrainIncrement,  
        const Matrix& C, 
        Vector& PlasticStrain, 
        const Properties& rMaterialProperties
    )  
    { 
        bool is_converged = false; 
        int iteration = 0, max_iter = 9000; 
        BoundedVector<double, TVoigtSize> DSigma, DS; 
        double PlasticConsistencyFactorIncrement = 0.0;     // Lambda 

        // Backward Euler  
        while (is_converged == false && iteration <= max_iter) 
        { 
            PlasticConsistencyFactorIncrement = UniaxialStress * PlasticDenominator; 
            if (PlasticConsistencyFactorIncrement < 0.0) PlasticConsistencyFactorIncrement = 0.0; 

            noalias(PlasticStrainIncrement) = PlasticConsistencyFactorIncrement * Gflux; // check 
            noalias(PlasticStrain) += PlasticStrainIncrement; 
            noalias(DS) = prod(C, PlasticStrainIncrement); 
            noalias(DSigma) -= DS; 
            noalias(PredictiveStressVector) -= DSigma; 

            CalculatePlasticParameters(PredictiveStressVector, UniaxialStress, Threshold, 
                PlasticDenominator, Fflux, Gflux, PlasticDissipation, PlasticStrainIncrement, 
                C, rMaterialProperties); 

            double F = UniaxialStress - Threshold; 

            if (F < std::abs(1.0e-8 * Threshold)) // Has converged 
            { 
                is_converged = true; 
            }
            else iteration++ 
        } 
        // if (iteration == max_iter) KRATOS_ERROR << "Reached max iterations inside the Plasticity loop" << std::endl; 
    }

    static void CalculatePlasticParameters(
        Vector& PredictiveStressVector, 
        double& UniaxialStress, 
        double& Threshold,
        double& PlasticDenominator, 
        Vector& Fflux, 
        Vector& Gflux, 
        double& PlasticDissipation, 
        Vector& PlasticStrainIncrement,
        const Matrix& C, 
        const Properties& rMaterialProperties
    )
    {
        BoundedVector<double, TVoigtSize> Deviator = ZeroVector(TVoigtSize); // TODO -> poner 2d o 3d?
        Vector HCapa = ZeroVector(6);
        double J2 = 0.0, r0 = 0.0, r1 = 0.0, Slope = 0.0, HardParam = 0.0;

        YieldSurfaceType::CalculateEquivalentStress(PredictiveStressVector, UniaxialStress, rMaterialProperties);
        CalculateDeviatorVector(PredictiveStressVector, Deviator, J2);
        CalculateFFluxVector(PredictiveStressVector, Deviator, J2, Fflux);
        CalculateGFluxVector(PredictiveStressVector, Deviator, J2, Gflux);
        CalculateRFactors(PredictiveStressVector, r0, r1);
        CalculatePlasticDissipation(PredictiveStressVector, r0,
            r1, PlasticStrainIncrement, PlasticDissipation, HCapa);
        CalculateEquivalentStressThreshold(PlasticDissipation, r0,
            r1, Threshold, Slope, rMaterialProperties);
        CalculateHardeningParameter(Fflux, Slope, HCapa, HardParam); // FFlux or GFlux????
        CalculatePlasticDenominator(Fflux, C, HardParam, PlasticDenominator)

    }

    // DF/DS
    static void CalculateFFluxVector(
        const Vector& StressVector, 
        const Vector& Deviator,
        const double J2, 
        Vector& FFluxVector
    )
    {

    }

    // DG/DS
    static void CalculateGFluxVector(
        const Vector& StressVector, 
        const Vector& Deviator,
        const double J2, 
        Vector& GFluxVector
    )
    {

    }
    
    // Calculates the McAully factors 
    static void CalculateRFactors(
        const Vector& StressVector, 
        double& r0, 
        double& r1
    )
    {
        Vector PrincipalStresses = ZeroVector(3);
		CalculatePrincipalStresses(PrincipalStresses, StressVector);

		double suma = 0.0, sumb = 0.0, sumc = 0.0;
		Vector SA = ZeroVector(3) , SB = ZeroVector(3), SC = ZeroVector(3);

		for (int i = 0; i < 3; i++)
		{
			SA[i] = abs(PrincipalStresses[i]);
			SB[i] = 0.5*(PrincipalStresses[i]  + SA[i]);
			SC[i] = 0.5*(-PrincipalStresses[i] + SA[i]);

			suma += SA[i];
			sumb += SB[i];
			sumc += SC[i];
		}
		if (suma != 0.0)
		{
			r0 = sumb/suma;
			r1 = sumc/suma;
		}
		else
		{
			r0 = sumb;
			r1 = sumc;
        }
    }

    // Calculates Plastic Dissipation
    static void CalculatePlasticDissipation(
        const Vector& StressVector, 
        const double r0,
        const double r1, 
        const Vector& PlasticStrainInc, 
        double& rCapap, 
        Vector& HCapa
    )
    {

    }

    // Calculates the stress threshold 
    static void CalculateEquivalentStressThreshold(
        const double Capap, 
        const double r0,
        const double r1, 
        double& rEquivalentStressThreshold, 
        double& rSlope,
        const Properties& rMaterialProperties
    )
    {

    }

    static void CalculateHardeningParameter(
        const Vector& FluxVector, 
        const double SlopeThreshold,
        const Vector& HCapa, 
        double& rHardParameter
    ) // todo which Flux=??????
    {

    }

    static void CalculatePlasticDenominator(
        const Vector& FluxVector, 
        const Matrix& C,
        const double HardParam, 
        double& PlasticDenominator
    )
    {

    }

    static void CalculatePrincipalStresses(
        Vector& rPrincipalStressVector, 
        const Vector StressVector
    )
	{
		rPrincipalStressVector.resize(3);
		double I1, I2, I3, phi, Num, Denom, II1;
		I1 = YieldSurfaceType::CalculateI1Invariant(StressVector, I1);
		I2 = YieldSurfaceType::CalculateI2Invariant(StressVector, I2);
		I3 = YieldSurfaceType::CalculateI3Invariant(StressVector, I3);
		II1 = I1*I1;

		Num = (2.0*II1 - 9.0*I2)*I1 + 27.0*I3;
		Denom = (II1 - 3.0*I2);

		if (Denom != 0.0)
		{
			phi = Num / (2.0*Denom*sqrt(Denom));

			if (std::abs(phi) > 1.0)
			{
				if (phi > 0.0) phi = 1.0;
				else phi = -1.0;
			}

			double acosphi = acos(phi);
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
