//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Máté Kelemen
//

// --- HDF5 Includes ---
#include "add_custom_predicates_to_python.h"
#include "custom_utilities/piped_model_predicate.h"


namespace Kratos::Python
{


void AddCustomPredicatesToPython(pybind11::module& rModule)
{
    #define KRATOS_DEFINE_PIPED_PREDICATE_BINDINGS(NAME)    \
        pybind11::class_<HDF5::NAME,                        \
                         HDF5::NAME::Pointer,               \
                         ModelPredicate>(rModule, #NAME)    \
            .def(pybind11::init<>())                        \
            .def(pybind11::init<const Parameters&>())       \
            .def("__call__", &HDF5::NAME::operator())

    KRATOS_DEFINE_PIPED_PREDICATE_BINDINGS(TimeIntervalPredicate);

    KRATOS_DEFINE_PIPED_PREDICATE_BINDINGS(StepIntervalPredicate);

    KRATOS_DEFINE_PIPED_PREDICATE_BINDINGS(PeriodicTimeIntervalPredicate);

    KRATOS_DEFINE_PIPED_PREDICATE_BINDINGS(PeriodicStepIntervalPredicate);

    #undef KRATOS_DEFINE_PIPED_PREDICATE_BINDINGS
}


} // namespace Kratos::Python
