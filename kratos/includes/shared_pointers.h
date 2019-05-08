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
//                   Pooyan Dadvand
//

#if !defined(KRATOS_MEMORY_H_INCLUDED )
#define  KRATOS_MEMORY_H_INCLUDED

/* System includes */
#include <utility>

/* External includes */
#include <memory>
#include "intrusive_ptr/intrusive_ptr.hpp"

namespace Kratos {

template<class T>
using intrusive_ptr = std::intrusive_ptr<T>; 

template<class T>
using intrusive_weak_ptr = std::intrusive_weak_ptr<T>; 

template<class T>
using shared_ptr = std::shared_ptr<T>;

template<class T>
using weak_ptr = std::weak_ptr<T>;

template<class T>
using unique_ptr = std::unique_ptr<T>;

template<typename C, typename...Args>
intrusive_ptr<C> make_intrusive(Args &&...args) {
    return std::make_intrusive<C>(std::forward<Args>(args)...);

}
template<typename C, typename...Args>
shared_ptr<C> make_shared(Args &&...args) {
    return std::make_shared<C>(std::forward<Args>(args)...);

}

template<typename C, typename...Args>
unique_ptr<C> make_unique(Args &&...args) {
    // Note: std::make_unique is C++14, this can be updated once we upgrade from C++11
    return unique_ptr<C>(new C(std::forward<Args>(args)...));
}

template<typename C, typename...Args>
shared_ptr<C> static_pointer_cast(Args &&...args) {
    return std::static_pointer_cast<C>(std::forward<Args>(args)...);
}

template<typename C, typename...Args>
shared_ptr<C> dynamic_pointer_cast(Args &&...args) {
    return std::dynamic_pointer_cast<C>(std::forward<Args>(args)...);
}

template<typename C, typename...Args>
shared_ptr<C> const_pointer_cast(Args &&...args) {
    return std::const_pointer_cast<C>(std::forward<Args>(args)...);
}

// template<typename C, typename...Args>
// shared_ptr<C> reinterpret_pointer_cast(Args &&...args) {
//     return std::reinterpret_pointer_cast<C>(std::forward<Args>(args)...);
// }
} // namespace Kratos


#define KRATOS_CLASS_POINTER_DEFINITION(a) typedef Kratos::shared_ptr<a > Pointer; \
typedef Kratos::shared_ptr<a > SharedPointer; \
typedef Kratos::weak_ptr<a > WeakPointer; \
typedef Kratos::unique_ptr<a > UniquePointer

#endif /* KRATOS_MEMORY_H_INCLUDED  defined */