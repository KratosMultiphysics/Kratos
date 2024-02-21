//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//                   Danilo Cavalcanti
//

#pragma once

// System includes
#include <iostream>

// Project includes
#include "includes/define.h"
#include "custom_saturation/saturation_law.hpp"
#include "custom_saturation/brooksandcorey_law.hpp"
#include "custom_saturation/vangenuchten_law.hpp"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos
{

class SaturationLawWrapper
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(SaturationLawWrapper);

    static Kratos::unique_ptr<SaturationLaw> Clone(const Properties& rMaterialProperties)
     {
         if (rMaterialProperties.Has(SATURATION_LAW_NAME))
         {
            const std::string &SaturationLawName = rMaterialProperties[SATURATION_LAW_NAME];
            
            if (SaturationLawName == "BrooksAndCoreyLaw")
                return Kratos::make_unique<BrooksAndCoreyLaw>();

            if (SaturationLawName == "VanGenuchtenLaw")
                return Kratos::make_unique<VanGenuchtenLaw>();

            KRATOS_ERROR << "Undefined SATURATION_LAW_NAME name" << SaturationLawName << std::endl;

            return nullptr;
         }

         // The default is Brooks and Corey Law
         return Kratos::make_unique<BrooksAndCoreyLaw>();

     }

}; // Class SaturationLawWrapper

}  // namespace Kratos.