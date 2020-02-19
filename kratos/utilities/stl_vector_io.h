//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//
//

#if !defined(KRATOS_STL_VECTOR_IO_H_INCLUDED )
#define  KRATOS_STL_VECTOR_IO_H_INCLUDED

// System includes
#include <iostream>
#include <vector>

namespace Kratos {

template<class T>
std::ostream& operator<<(std::ostream& rOStream, const std::vector<T>& rVec) {

    std::size_t vector_size = rVec.size();

    rOStream << "[";
    if(vector_size>0) rOStream << rVec[0];
    if(vector_size>1) {
        for(std::size_t i = 1; i < vector_size; i++)
            rOStream<<", "<<rVec[i];
    }
    rOStream << "]";

    return rOStream;
}

} //namespace Kratos


#endif // KRATOS_STL_VECTOR_IO_H_INCLUDED  defined
