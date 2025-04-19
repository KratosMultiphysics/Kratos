// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Mate Kelemen
//

// Project includes
#include "custom_python/add_custom_constraints_to_python.hpp" // AddCustomConstraintsToPython
#include "custom_constraints/link_constraint.hpp" // LinkConstraint


namespace Kratos::Python {


void AddCustomConstraintsToPython(pybind11::module& rModule)
{
    pybind11::class_<LinkConstraint, LinkConstraint::Pointer, MasterSlaveConstraint>(rModule, "LinkConstraint")
        .def(pybind11::init<const LinkConstraint::IndexType,Node&,Node&,const std::size_t,bool>())
        ;
}


} // namespace Kratos::Python
