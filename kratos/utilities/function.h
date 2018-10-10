//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    
//
//

#if !defined(KRATOS_FUNCTION_H_INCLUDED)
#define  KRATOS_FUNCTION_H_INCLUDED

// System includes
#include <functional>

// Project includes
#include "utilities/function_fwd.h"

namespace Kratos
{

///@name Kratos Classes
///@{

template <class T>
class Function : public std::function<T>
{
    public:
    Function() : std::function<T>() {}

    template <class F>
    Function(F f) : std::function<T>(f) {}

    private:
    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
    }

    void load(Serializer& rSerializer)
    {
    }
};

template <class T>
std::ostream& operator<<(std::ostream& os, const Function<T>& f)
{
    return os;
}

}  // namespace Kratos.

#endif // KRATOS_FUNCTION_H_INCLUDED  defined 
