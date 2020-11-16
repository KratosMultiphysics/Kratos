//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//


// System includes


// External includes


// Project includes
#include "mapping_matrix_builder.h"


namespace Kratos
{




///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class MappingMatrixBuilder< false >; // non-distributed version
template class MappingMatrixBuilder< true >; // distributed version

}  // namespace Kratos.


