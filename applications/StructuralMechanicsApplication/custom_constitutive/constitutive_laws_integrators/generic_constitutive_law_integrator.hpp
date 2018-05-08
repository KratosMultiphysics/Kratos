//
//   Project Name:        KratosStructuralMechanicsApplication $
//   Created by:          $Author:            A.Cornejo        $
//   Last modified by:    $Co-Author:                          $
//   Date:                $Date:                April 2018     $
//   Revision:            $Revision:                  0.0      $
//
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
    template <class  YieldSurfaceType, class PlasticPotentialType>
    class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) GenericConstitutiveLawIntegrator
    {

        public:
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

            /// Clone
            // GenericConstitutiveLawIntegrator::Pointer Clone() const 
            // {
            //     GenericConstitutiveLawIntegrator<class YieldSurfaceType>::Pointer p_clone(new GenericConstitutiveLawIntegrator<class YieldSurfaceType>(*this));
            //     return p_clone;
            // }

            // ***************************************************************************
            // ***************************************************************************

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

        protected:

            //typename YieldSurfaceType::Pointer mpYieldSurface;


    };
}
#endif