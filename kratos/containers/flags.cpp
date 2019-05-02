//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "containers/flags.h"
#include "includes/serializer.h"

namespace Kratos
{

    void Flags::save(Serializer& rSerializer) const
    {
        rSerializer.save("IsDefined",  mIsDefined);
        rSerializer.save("Flags",  mFlags);
    }

    void Flags::load(Serializer& rSerializer)
    {
        rSerializer.load("IsDefined",  mIsDefined);
        rSerializer.load("Flags",  mFlags);
    }

  
  
}  // namespace Kratos.


