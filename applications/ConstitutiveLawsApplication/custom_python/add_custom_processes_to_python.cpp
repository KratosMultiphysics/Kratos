// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Sergio Jimenez/Alejandro Cornejo
//

// System includes

// External includes


// Project includes
#include "custom_python/add_custom_processes_to_python.h"
#include "constitutive_laws_application_variables.h"
#include "includes/model_part.h"

//Processes
#include "custom_processes/advance_in_time_high_cycle_fatigue_process.h"

namespace Kratos {
    namespace Python {
        void  AddCustomProcessesToPython(pybind11::module& m)
        {
            namespace py = pybind11;

            /// Processes
            py::class_<AdvanceInTimeHighCycleFatigueProcess, AdvanceInTimeHighCycleFatigueProcess::Pointer, Process>(m, "AdvanceInTimeHighCycleFatigueProcess")
                .def(py::init< ModelPart&, Parameters >())
                ;
        }
    } // namespace Python.
} // Namespace Kratos

