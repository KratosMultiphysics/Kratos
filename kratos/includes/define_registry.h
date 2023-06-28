//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos Roig
//                   Ruben Zorrilla
//

#pragma once

// System includes
#include <stdexcept>
#include <sstream>


// External includes

// Project includes
#include "includes/registry.h"
#include "includes/registry_item.h"

#define KRATOS_REGISTRY_NAME_(A,B) A##B
#define KRATOS_REGISTRY_NAME(A,B) KRATOS_REGISTRY_NAME_(A,B)

/**
 * @brief Macro to register class prototypes
 * This macro creates a static bool in the class (which value is not relevant and is supposed to be never accessed)
 * This macro must be placed within the private section of the class (e.g. Process, Controller, DerivedProcess, etc)
 * The macro forces a function to be evaluated in the declaration of the class
 * First it checks if the provided key has been already registered (note that it checks the key but not the object)
 * to then create a dispacher function that will initialize the registry object the first time it is accessed
 * NOTE 1: The return type is intentionally the HasItem of the added item (and not a simple bool) in order to avoid the 
 * compiler to statically analize the function.
 * NOTE 2: Also note that it is explicitly needed to define the dispatcher as TFunctionType or, in other words, putting 
 * the dispatcher function directly inside the std::move would break current code as the compiler will deduce the type
 * of the lambda function to something that is not an std::function.
 */
#define KRATOS_REGISTRY_ADD_PROTOTYPE(NAME, X)                                            \
    static inline bool KRATOS_REGISTRY_NAME(_is_registered_, __LINE__) = []() -> bool {   \
        using TFunctionType = std::function<std::shared_ptr<X>()>;                        \
        std::string key_name = NAME + std::string(".") + std::string(#X);                 \
        if (!Registry::HasItem(key_name))                                                 \
        {                                                                                 \
            auto &r_item = Registry::AddItem<RegistryItem>(key_name);                     \
            TFunctionType dispatcher = [](){return std::make_shared<X>();};               \
            r_item.AddItem<TFunctionType>("Prototype", std::move(dispatcher));            \
        }                                                                                 \
        return Registry::HasItem(key_name);                                               \
    }();

#define KRATOS_REGISTRY_ADD_TEMPLATE_PROTOTYPE(NAME, X, ...)                                            \
    static inline bool KRATOS_REGISTRY_NAME(_is_registered_, __LINE__) = []() -> bool {                 \
        using TFunctionType = std::function<std::shared_ptr<X<__VA_ARGS__>>()>;                         \
        std::string key_name = NAME + std::string(".") + std::string(#X) + "<" + #__VA_ARGS__ + ">";    \
        if (!Registry::HasItem(key_name))                                                               \
        {                                                                                               \
            auto &r_item = Registry::AddItem<RegistryItem>(key_name);                                   \
            TFunctionType dispatcher = [](){return std::make_shared<X<__VA_ARGS__>>();};                \
            r_item.AddItem<TFunctionType>("Prototype", std::move(dispatcher));                          \
        } else {                                                                                        \
        }                                                                                               \
        return Registry::HasItem(key_name);                                                             \
    }();
