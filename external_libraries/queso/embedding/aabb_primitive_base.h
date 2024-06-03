// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef AABB_PRIMITIVE_BASE_INCLUDE_H
#define AABB_PRIMITIVE_BASE_INCLUDE_H

//// STL includes
#include <cstddef>
#include <array>
//// Project includes
#include "queso/includes/define.hpp"

namespace queso {

///@name QuESo Classes
///@{

/// Forward declaration
class AABB_primitive;

/**
 * @class  AABB_primitive
 * @author Manuel Messmer
 * @brief  Base class for aabb primitives. Derived classes must override intersect().
*/
class AABB_primitive_base
{
public:
    ///@}
    ///@name Life cycle
    ///@{

    // Destructor
    virtual ~AABB_primitive_base() = default;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns true, if AABB intersect with this object. Interface for AABB tree.
    /// @param aabb
    /// @return bool
    virtual bool intersect(const AABB_primitive &aabb) const = 0;

    ///@}
}; // End AABB_primitive_base class
///@} // End QuESo classes

} // End namespace queso

#endif // AABB_PRIMITIVE_BASE_INCLUDE_H
