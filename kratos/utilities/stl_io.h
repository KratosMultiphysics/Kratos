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
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> & data) {

    std::cout << "[";
    std::copy(data.begin(), data.end(), std::ostream_iterator<T>(std::cout, ", "));
    std::cout << "[";
    
    return os;
}

#endif // KRATOS_KRATOS_STL_IO_H_INCLUDED  defined
