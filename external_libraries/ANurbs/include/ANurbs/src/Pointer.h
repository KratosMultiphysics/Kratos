#pragma once

#include <memory>

namespace ANurbs {

template <typename T>
using Pointer = std::shared_ptr<T>;

template <typename T>
using Unique = std::unique_ptr<T>;

template <typename T, typename... TArgs>
Unique<T>
New(
    TArgs&&... args)
{
    return Unique<T>(new T(std::forward<TArgs>(args)...));
}

} // namespace ANurbs
