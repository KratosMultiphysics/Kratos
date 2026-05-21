//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Pooyan Dadvand
//

#pragma once

// System includes
#include <iostream>
#include <utility>
#include <memory>

// External includes
#include "intrusive_ptr/intrusive_ptr.hpp"

namespace Kratos {

template<class T>
using shared_ptr = std::shared_ptr<T>;

template<class T>
using weak_ptr = std::weak_ptr<T>;

template<class T>
using unique_ptr = std::unique_ptr<T>;

template<typename C, typename...Args>
intrusive_ptr<C> make_intrusive(Args &&...args) {
    return intrusive_ptr<C>(new C(std::forward<Args>(args)...));
}
template<typename C, typename...Args>
shared_ptr<C> make_shared(Args &&...args) {
    return std::make_shared<C>(std::forward<Args>(args)...);
}

template<typename C, typename...Args>
unique_ptr<C> make_unique(Args &&...args) {
    return std::make_unique<C>(std::forward<Args>(args)...);
}

template<class T>
std::ostream& operator <<(std::ostream& rOStream, const Kratos::weak_ptr<T>& rData) {

  if(!rData.expired())
    rOStream << *rData.lock().get();
  else
    rOStream <<" expired weak_ptr ";

  return rOStream;
}

template<class T>
std::ostream& operator <<(std::ostream& rOStream, const Kratos::intrusive_ptr<T>& rData) {
  rOStream << *rData.get();
  return rOStream;
}

} // namespace Kratos


#define KRATOS_CLASS_POINTER_DEFINITION(a) typedef Kratos::shared_ptr<a > Pointer; \
typedef Kratos::shared_ptr<a > SharedPointer; \
typedef Kratos::weak_ptr<a > WeakPointer; \
typedef Kratos::unique_ptr<a > UniquePointer

namespace Kratos {
template< class T > class GlobalPointer;
}

#define KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(a) typedef typename Kratos::intrusive_ptr<a > Pointer; \
typedef Kratos::GlobalPointer<a > WeakPointer; \
typedef Kratos::unique_ptr<a > UniquePointer; \
typename a::Pointer shared_from_this(){ return a::Pointer(this); }
