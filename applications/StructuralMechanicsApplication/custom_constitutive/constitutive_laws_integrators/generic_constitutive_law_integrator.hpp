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
    template <typename YieldSurfaceType>
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

            static void IntegrateStressVector() // While Loop plasticity
            {

            }


        protected:

            //typename YieldSurfaceType::Pointer mpYieldSurface;


    };
}
#endif