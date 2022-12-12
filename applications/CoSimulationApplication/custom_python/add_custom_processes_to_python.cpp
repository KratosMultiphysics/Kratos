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
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_processes/data_transfer_3D_1D_process.h"

namespace Kratos::Python{

    void  AddCustomProcessesToPython(pybind11::module& m)
    {
        pybind11::class_< DataTransfer3D1DProcess>(m, "DataTransfer3D1DProcess")
            .def_static("From3Dto1DDataTransfer", &DataTransfer3D1DProcess::From3Dto1DDataTransfer)
            .def_static("From1Dto3DDataTransfer", &DataTransfer3D1DProcess::From1Dto3DDataTransfer)
            ;
    }

}  // namespace Kratos::Python.
