//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Angel Celigueta
//

// System includes

// External includes

// Project includes
#include "includes/gid_io.h"

namespace Kratos
{
    int GidIOBase::msLiveInstances = 0;
    template class GidIO<GidGaussPointsContainer, GidMeshContainer>;
}  // namespace Kratos.


