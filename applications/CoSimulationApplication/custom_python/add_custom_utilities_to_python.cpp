// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//					 license: CoSimulationApplication/license.txt
//
//  Main authors:
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/feti_dynamic_coupling_utilities.h"

namespace Kratos{

namespace Python{

    void  AddCustomUtilitiesToPython(pybind11::module& m)
    {
        pybind11::class_< FetiDynamicCouplingUtilities>(m, "FetiDynamicCouplingUtilities")
            .def(pybind11::init<ModelPart&, ModelPart&, Parameters>())
            .def("SetOriginAndDestinationDomainsWithInterfaceModelParts",
                &FetiDynamicCouplingUtilities::SetOriginAndDestinationDomainsWithInterfaceModelParts)
            .def("SetEffectiveStiffnessMatrix",
                pybind11::overload_cast<CompressedMatrix&, IndexType>(&FetiDynamicCouplingUtilities::SetEffectiveStiffnessMatrix))
            .def("SetEffectiveStiffnessMatrix",
                pybind11::overload_cast<IndexType>(&FetiDynamicCouplingUtilities::SetEffectiveStiffnessMatrix))
            .def("EquilibrateDomains",
                &FetiDynamicCouplingUtilities::EquilibrateDomains)
            .def("SetOriginInitialKinematics",
                &FetiDynamicCouplingUtilities::SetOriginInitialKinematics)
            .def("SetMappingMatrix",
                &FetiDynamicCouplingUtilities::SetMappingMatrix)
            .def("SetLinearSolver",
                &FetiDynamicCouplingUtilities::SetLinearSolver)
            ;
    }

}  // namespace Python.
} // Namespace Kratos

