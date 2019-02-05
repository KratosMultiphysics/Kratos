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
#include <iostream>
#include <vector>
#include <set>
#include <memory>

// Project includes
#include "includes/define.h"

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

template<class T>
std::ostream& operator <<(std::ostream& rOStream, const std::set<T>& rSet) {

  const std::size_t set_size = rSet.size();

  rOStream << "[";
  if(set_size>0) rOStream << *(rSet.begin());
  if(set_size>1) {
    for(auto it(std::next(rSet.begin(),1)); it!=rSet.end(); ++it)
      rOStream<<", "<<*it;
  }
  rOStream << "]";

  return rOStream;
}

template<class T>
std::ostream& operator <<(std::ostream& rOStream, const std::weak_ptr<T>& rData) {

  if(!rData.expired())
    rOStream << *rData.lock().get();
  else
    rOStream <<" expired weak_ptr ";

  return rOStream;
}


} //namespace Kratos


#endif // KRATOS_KRATOS_STL_IO_H_INCLUDED  defined
