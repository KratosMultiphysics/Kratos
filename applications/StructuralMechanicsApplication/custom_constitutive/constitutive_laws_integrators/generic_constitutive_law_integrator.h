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
template <class TYieldSurfaceType>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) GenericConstitutiveLawIntegrator
{
public:
    ///@name Type Definitions
    ///@{

    struct PlasticParameters
    {
        Vector PredictiveStressVector; /// This bka
        double UniaxialStress;
        double Kp;
        double PlasticDenominator; 
        Vector Fflux;
        Vector Gflux;
        double PlasticDissipation;
        Vector PlasticStrainIncrement;
        Vector PlasticStrain;

        void Initialize()
        {
            
        }
    }

    /// The type of yield surface
    typedef TYieldSurfaceType YieldSurfaceType;

    /// The type of potential plasticity
    typedef YieldSurfaceType::TPlasticPotentialType PlasticPotentialType;

    /// Counted pointer of GenericConstitutiveLawIntegrator
    KRATOS_CLASS_POINTER_DEFINITION(GenericConstitutiveLawIntegrator);

    /// Initialization constructor.
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
        double& Kp,
        double& PlasticDenominator, 
        Vector& Fflux, 
        Vector& Gflux, 
        double& Capap, 
        Vector& PlasticStrainIncrement,
        const Matrix& C,
        Vector& PlasticStrain
        )
    {

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
        Vector Deviator = ZeroVector(6); // TODO -> poner 2d o 3d?
        Vector HCapa = ZeroVector(6);
        double J2 = 0.0, r0 = 0.0, r1 = 0.0, Slope = 0.0, HardParam = 0.0;

        YieldSurfaceType::CalculateEquivalentStress(PredictiveStressVector, UniaxialStress, rMaterialProperties);
        this->CalculateDeviatorVector(PredictiveStressVector, Deviator, J2);
        this->CalculateFFluxVector(PredictiveStressVector, Deviator, J2, Fflux);
        this->CalculateGFluxVector(PredictiveStressVector, Deviator, J2, Gflux);
        this->CalculateRFactors(PredictiveStressVector, r0, r1);
        this->CalculatePlasticDissipation(PredictiveStressVector, r0,
            r1, PlasticStrainIncrement, PlasticDissipation, HCapa);
        this->CalculateEquivalentStressThreshold(PlasticDissipation, r0,
            r1, Threshold, Slope, rMaterialProperties);
        this->CalculateHardeningParameter(Fflux, Slope, HCapa, HardParam); // FFlux or GFlux????
        this->CalculatePlasticDenominator(Fflux, C, HardParam, PlasticDenominator)

    }

    // DF/DS
    static void CalculateFFluxVector(const Vector& StressVector, const Vector& Deviator,
        const double J2, Vector& FFluxVector)
    {

    }

    // DG/DS
    static void CalculateGFluxVector(const Vector& StressVector, const Vector& Deviator,
        const double J2, Vector& GFluxVector)
    {

    }

    static void CalculateRFactors(const Vector& StressVector, double& r0, double& r1)
    {

    }

    // Calculates Capap
    static void CalculatePlasticDissipation(const Vector& StressVector, const double r0,
        const double r1, const Vector& PlasticStrainInc, double& rCapap, Vector& HCapa)
    {

    }

    // Calculates Kp
    static void CalculateEquivalentStressThreshold(const double Capap, const double r0,
        const double r1, double& rEquivalentStressThreshold, double& rSlope)
    {

    }

    static void CalculateHardeningParameter(const Vector& FluxVector, const double SlopeThreshold,
        const Vector& HCapa, double& rHardParameter) // todo which Flux=??????
    {

    }

    static void CalculatePlasticDenominator(const Vector& FluxVector, const Matrix& C,
        const double HardParam, double& PlasticDenominator)
    {

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
