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

// --- Core Includes ---
#include "add_predicates_to_python.h"
#include "utilities/model_predicate.h"
#include "includes/define.h"


namespace Kratos::Python
{


void AddPredicatesToPython(pybind11::module& rModule)
{
    pybind11::class_<ModelPredicate, ModelPredicate::Pointer>(rModule, "ModelPredicate")
        .def("__call__", &ModelPredicate::operator())
        ;
}


} // namespace Kratos::Python
