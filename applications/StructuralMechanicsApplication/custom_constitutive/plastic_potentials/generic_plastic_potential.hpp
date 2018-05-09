// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo

#if !defined(KRATOS_GENERIC_PLASTIC_POTENTIAL_H_INCLUDED)
#define  KRATOS_GENERIC_PLASTIC_POTENTIAL_H_INCLUDED

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
    class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) GenericPlasticPotential
    {

        public:

            KRATOS_CLASS_POINTER_DEFINITION(GenericPlasticPotential);

            /// Initialization constructor.
            GenericPlasticPotential()
            {
            }

            /// Copy constructor
            GenericPlasticPotential(GenericPlasticPotential const& rOther)
            {
            }

            /// Assignment operator
            GenericPlasticPotential& operator=(GenericPlasticPotential const& rOther)
            {
                return *this;
            }

            /// Destructor
            virtual ~GenericPlasticPotential() {};

            /// Clone
            // GenericPlasticPotential::Pointer Clone() const 
            // {
            //     GenericPlasticPotential::Pointer p_clone(new GenericPlasticPotential(*this));
            //     return p_clone;
            // }
            
            // ***************************************************************************
            // ***************************************************************************

            static void CalculatePlasticPotentialDerivative(const Vector& StressVector, const Vector& Deviator,const double& J2, Vector& rg)
            {

            }


    };
}
#endif