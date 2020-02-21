//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_REGISTERED_COMPONENT_LOOKUP_H_INCLUDED)
#define KRATOS_REGISTERED_COMPONENT_LOOKUP_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_components.h"

// Application includes

namespace Kratos
{

/// A class for calling a functor on a variable specified by it's registered name.
/**
 * Example:
 *
 * template <typename TVariable>
 * class CopyVariableFunctor
 * {
 *    public:
 *    void operator()(TVariable const& rVariable,
 *                    DataValueContainer const& rSrc)
 *                    DataValueContainer& rDest)
 *  {
 *      rDest[rVariable] = rSrc[rVariable];
 *  }
 * };
 *
 * RegisteredComponentLookup<Variable<double>>("PRESSURE").Execute<CopyVariableFunctor>(rSrc, rDest);
 */
template <typename ...TVariables> 
class RegisteredComponentLookup
{
    template <template<typename T> class TFunctor, typename ...Targs>
    struct FunctorWrapper
    {
        template <typename TVariable>
        void Execute(std::string const& rName, bool& found, Targs&... args)
        {
            if (!found && KratosComponents<TVariable>::Has(rName))
            {
                const TVariable& rVariable = KratosComponents<TVariable>::Get(rName);
                TFunctor<TVariable>()(rVariable, args...);
                found = true;
            }
        }
    };

public:
    explicit RegisteredComponentLookup(std::string const& rName) : mName(rName)
    {}

  template <template<typename T> class TFunctor, typename ...Targs>
  void Execute(Targs&... args)
  {
      bool found = false;
      int dummy[sizeof...(TVariables)] = {(FunctorWrapper<TFunctor, Targs...>().template Execute<TVariables>(mName, found, args...), 0)...};
      ignore_unused_variable_warning(dummy);
      KRATOS_ERROR_IF(!found) << "Variable \"" << mName << "\" was not found.\n";
  }

  private:
  std::string mName;

  void ignore_unused_variable_warning(int[]) {}
};

} // namespace Kratos.

#endif // KRATOS_REGISTERED_COMPONENT_LOOKUP_H_INCLUDED defined
