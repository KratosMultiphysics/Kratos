//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Carlos A. Roig
//

// System includes

// External includes

// Project includes

/* This is only valid for C++20. Uploaded here as a prof of concept */

namespace Kratos {

    template<typename TDataType, typename TInterfaceType, typename TInterfaceFunctor>
    void AddToInterfaceFold(TInterfaceType& rInterface, TInterfaceFunctor & rFunctor) {
        rFunctor.template operator()<TDataType>(rInterface);
    }

    template<typename... args, typename TInterfaceType, typename TInterfaceFunctor>
    static void AddToInterface(TInterfaceType& rInterface, TInterfaceFunctor && rFunctor) {
        (AddToInterfaceFold<args, TInterfaceType, TInterfaceFunctor>(rInterface, rFunctor), ...);
    }

} // namespace Kratos

