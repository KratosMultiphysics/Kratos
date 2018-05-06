//
//   Project Name:        KratosStructuralMechanicsApplication $
//   Created by:          $Author:            A.Cornejo        $
//   Last modified by:    $Co-Author:                          $
//   Date:                $Date:                April 2018     $
//   Revision:            $Revision:                  0.0      $
//
//

#if !defined(KRATOS_GENERIC_YIELD_SURFACE_H_INCLUDED)
#define  KRATOS_GENERIC_YIELD_SURFACE_H_INCLUDED

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

    template <class PlasticPotentialType>
    class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) GenericYieldSurface
    {

        public:

            KRATOS_CLASS_POINTER_DEFINITION( GenericYieldSurface );

            /// Initialization constructor.
            GenericYieldSurface()
            {
                mpPlasticPotential = PlasticPotentialType().Clone();
            }

            /// Copy constructor
            GenericYieldSurface(GenericYieldSurface const& rOther)
            {
                mPlasticPotential = rOther.mPlasticPotential;
            }

            /// Assignment operator
            GenericYieldSurface& operator=(GenericYieldSurface const& rOther)
            {
                mPlasticPotential = rOther.mPlasticPotential;
                return *this;
            }

            /// Destructor
            virtual ~GenericYieldSurface() {};

            /// Clone
            // GenericYieldSurface::Pointer Clone() const override
            // {
            //     GenericYieldSurface<class PlasticPotentialType>::Pointer p_clone(new GenericYieldSurface<class PlasticPotentialType>(*this));
            //     return p_clone;
            // }


            // ***************************************************************************
            // ***************************************************************************

            static void CalculateEquivalentStress(const Vector& StressVector, double& rEqStress)
            {
                // Implement for each yield surf
            }

            void CalculateI1Invariant(const Vector& StressVector, double& rI1)
            {
                rI1 = StressVector[0] + StressVector[1] + StressVector[2];
            }

            void CalculateI2Invariant(const Vector& StressVector, double& rI2)
            {
                rI2 = (StressVector[0] + StressVector[2])*StressVector[1] + StressVector[0]*StressVector[2] +
                    - StressVector[3]*StressVector[3] - StressVector[4]*StressVector[4] - StressVector[5]*StressVector[5];
            }

            void CalculateI3Invariant(const Vector& StressVector, double& rI3)
            {
                rI3 = (StressVector[1]*StressVector[2] - StressVector[4]*StressVector[4])*StressVector[0] -
                    StressVector[1]*StressVector[5]*StressVector[5] - StressVector[2]*StressVector[3]*StressVector[3] +
                    2.0*StressVector[3]*StressVector[4]*StressVector[5];
            }

            void CalculateJ2Invariant(const Vector& StressVector, const double& I1, Vector& rDeviator, double& rJ2)
            {
                rDeviator = StressVector;
                double Pmean = I1 / 3.0;

                rDeviator[0] -= Pmean;
                rDeviator[1] -= Pmean;
                rDeviator[2] -= Pmean;

                rJ2 = 0.5*(rDeviator[0]*rDeviator[0] + rDeviator[1]*rDeviator[1] + rDeviator[2]*rDeviator[2]) +
                    (rDeviator[3]*rDeviator[3] + rDeviator[4]*rDeviator[4] + rDeviator[5]*rDeviator[5]);
            }

            // Computes dG/dS
            void CalculatePlasticPotentialDerivative(const Vector& StressVector,const Vector& Deviator,const double& J2, Vector& rg)
            {
                //mPlasticPotential->CalculatePlasticPotentialDerivative(StressVector,Deviator,J2,rg);
                PlasticPotentialType::CalculatePlasticPotentialDerivative(StressVector,Deviator,J2,rg);
            }

		protected:

			//typename PlasticPotentialType::Pointer mpPlasticPotential;

    };

}

#endif