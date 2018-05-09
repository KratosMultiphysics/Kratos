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
 * @author Alejandro Cornejo
 */
template <class  YieldSurfaceType, class PlasticPotentialType>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) GenericConstitutiveLawIntegrator
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of GenericConstitutiveLawIntegrator
    KRATOS_CLASS_POINTER_DEFINITION(GenericConstitutiveLawIntegrator);

    /// Initialization constructor.
    GenericConstitutiveLawIntegrator()
    {
        //mpYieldSurface = YieldSurfaceType().Clone();
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

//     /// Clone
//     GenericConstitutiveLawIntegrator::Pointer Clone() const
//     {
//         GenericConstitutiveLawIntegrator<class YieldSurfaceType>::Pointer p_clone(new GenericConstitutiveLawIntegrator<class YieldSurfaceType>(*this));
//         return p_clone;
//     }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    static void IntegrateStressVector(Vector& PredictiveStressVector, double& UniaxialStress, double& Kp,
        double& PlasticDenominator, Vector& Fflux, Vector& Gflux, double& Capap, Vector& PlasticStrainIncrement,
        const Matrix& C, Vector& PlasticStrain)
    {

    }

    static void CalculatePlasticParameters(Vector& PredictiveStressVector, double& UniaxialStress, double& Kp,
        double& PlasticDenominator, Vector& Fflux, Vector& Gflux, double& Capap, Vector& PlasticStrainIncrement,
        const Matrix& C)
    {

    }

    // DF/DS
    static void CalculateFFluxVector(const Vector& StressVector, const Vector& Deviator,
        const double& J2, Vector& FFluxVector)
    {

    }

    // DG/DS
    static void CalculateGFluxVector(const Vector& StressVector, const Vector& Deviator,
        const double& J2, Vector& GFluxVector)
    {

    }

    static void CalculateRFactors(const Vector& StressVector, double& r0, double& r1)
    {

    }

    // Calculates Capap
    static void CalculatePlasticDissipation(const Vector& StressVector, const double& r0,
        const double& r1, const Vector& PlasticStrainInc, double& rCapap, Vector& HCapa)
    {

    }

    // Calculates Kp
    static void CalculateEquivalentStressThreshold(const double& Capap, const double& r0,
        const double& r1, double& rEquivalentStressThreshold, double& rSlope)
    {

    }

    static void CalculateHardeningParameter(const Vector& FluxVector, const double& SlopeThreshold,
        const Vector& HCapa, double& rHardParameter) // todo which Flux=??????
    {

    }

    static void CalculatePlasticDenominator(const Vector& FluxVector, const Matrix& C,
        const double& HardParam, double& PlasticDenominator)
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

//     typename YieldSurfaceType::Pointer mpYieldSurface;

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
//         rSerializer.save("YieldSurface", mpYieldSurface);
    }

    void load(Serializer& rSerializer)
    {
//         rSerializer.load("YieldSurface", mpYieldSurface);
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
