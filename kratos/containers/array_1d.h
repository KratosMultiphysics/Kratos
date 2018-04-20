//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//                    
//


#if	!defined(KRATOS_ARRAY_1D_H_INCLUDED	)
#define	 KRATOS_ARRAY_1D_H_INCLUDED



// System includes

// External	includes

// Project includes
#include "includes/amatrix_interface.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type	Definitions
///@{
    
template<typename TDataType, std::size_t TSize1> 
using array_1d = AMatrix::Matrix<TDataType, TSize1, 1>;

///@}


}  // namespace	Kratos.

#endif // KRATOS_ARRAY_1D_H_INCLUDED  defined 
