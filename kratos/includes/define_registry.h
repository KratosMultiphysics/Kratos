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

// Based on https://stackoverflow.com/questions/8709340
#define KRATOS_REGISTRY_ADD_BASE() \
    bool __kratos_registry_base = true;

namespace Kratos::_Registry 
{
    template<class T, class R>
    T base_of(R T::*);

    template<class T>
    class has_registry_base 
    {
        private:
            typedef char ok_type;
            typedef long ko_type;
            template<typename U> static ok_type test(decltype(&U::__kratos_registry_base));
            template<typename U> static ko_type test(...);
        public:
            static constexpr bool value = sizeof(test<T>(0)) == sizeof(ok_type);
    };

    template<class T, class TRitem>
    void _registry_add(TRitem & r_item, std::string key, std::string name, std::string value)  
    {
        using TFunctionType = std::function<std::shared_ptr<T>()>;

        TFunctionType dispatcher = [](){return std::make_shared<T>();};
        r_item.template AddItem<TFunctionType>(value, std::move(dispatcher));
    };

    template<class T, class TRitem>
    typename std::enable_if<has_registry_base<T>::value == true, void>::type _registry_add_base(TRitem & r_item, std::string key, std::string name, std::string value) 
    {
        using TBase = decltype(base_of(&T::__kratos_registry_base));

        _registry_add<TBase>(r_item, key, name, value);
    };

    template<class T, class TRitem>
    typename std::enable_if<has_registry_base<T>::value == false, void>::type _registry_add_base(TRitem & r_item, std::string key, std::string name, std::string value) {};
} // namespace Kratos::_Registry

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
// #define KRATOS_REGISTRY_ADD_PROTOTYPE(NAME, X)                                            \
//     static inline bool KRATOS_REGISTRY_NAME(_is_registered_, __LINE__) = []() -> bool {   \
//         using TType = X;                                                                  \
//         using TBase = decltype(base_of(&X::__kratos_registry_base));                      \
//         using TFunctionType = std::function<std::shared_ptr<TType>()>;                    \
//         using TFunctionBase = std::function<std::shared_ptr<TBase>()>;                    \
//         std::string key_name = NAME + std::string(".") + std::string(#X);                 \
//         if (!Registry::HasItem(key_name))                                                 \
//         {                                                                                 \
//             auto &r_item = Registry::AddItem<RegistryItem>(key_name);                     \
//             TFunctionType pdispatcher = [](){return std::make_shared<TType>();};          \
//             TFunctionBase bdispatcher = [](){return std::make_shared<TBase>();};          \
//             r_item.AddItem<TFunctionType>("Prototype", std::move(pdispatcher));           \
//             r_item.AddItem<TFunctionType>("Base",      std::move(bdispatcher));           \
//         }                                                                                 \
//         return Registry::HasItem(key_name);                                               \
//     }();

#define KRATOS_REGISTRY_ADD_PROTOTYPE(NAME, X)                                              \
    static inline bool KRATOS_REGISTRY_NAME(_is_registered_, __LINE__) = []() -> bool {     \
        std::string key_name = NAME + std::string(".") + std::string(#X);                   \
        if (!Registry::HasItem(key_name))                                                   \
        {                                                                                   \
            auto &r_item = Registry::AddItem<RegistryItem>(key_name);                       \
            Kratos::_Registry::_registry_add<X>(r_item, NAME, #X, "Prototype");             \
            Kratos::_Registry::_registry_add_base<X>(r_item, NAME, #X, "Base");             \
        }                                                                                   \
        return Registry::HasItem(key_name);                                                 \
    }();

#define KRATOS_REGISTRY_ADD_TEMPLATE_PROTOTYPE(NAME, X, ...)                                            \
    static inline bool KRATOS_REGISTRY_NAME(_is_registered_, __LINE__) = []() -> bool {                 \
        using TFunctionType = std::function<std::shared_ptr<X<__VA_ARGS__>>()>;                         \
        std::string key_name = NAME + std::string(".") + std::string(#X) + "<" + Kratos::Registry::RegistryTemplateToString(__VA_ARGS__) + ">";    \
        if (!Registry::HasItem(key_name))                                                               \
        {                                                                                               \
            auto &r_item = Registry::AddItem<RegistryItem>(key_name);                                   \
            TFunctionType dispatcher = [](){return std::make_shared<X<__VA_ARGS__>>();};                \
            r_item.AddItem<TFunctionType>("Prototype", std::move(dispatcher));                          \
        } else {                                                                                        \
        }                                                                                               \
        return Registry::HasItem(key_name);                                                             \
    }();
