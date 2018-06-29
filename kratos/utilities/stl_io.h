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

#if !defined(KRATOS_KRATOS_STL_IO_H_INCLUDED )
#define  KRATOS_KRATOS_STL_IO_H_INCLUDED

// System includes
#include <algorithm>
#include <iterator>
#include <iostream>
#include <vector>

// Project includes
#include "includes/define.h"

// Std::vecotr << operator

namespace Kratos {

// Std::vector << operator
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> & data) {

    std::size_t vector_size = data.size();

    os << "[";
    if(vector_size>0) os << data[0];
    if(vector_size>1) {
        for(std::size_t i = 1; i < vector_size; i++) {
            os<<", "<<data[i];
        }
    }
    os << "]";

    return os;
}

} //namespace Kratos


#endif // KRATOS_KRATOS_STL_IO_H_INCLUDED  defined
