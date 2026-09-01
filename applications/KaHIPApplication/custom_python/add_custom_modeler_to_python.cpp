//     __ __      __  __________  ___                ___            __  _
//    / //_/___ _/ / / /  _/ __ \/   |  ____  ____  / (_)________ _/ /_(_)___  ____
//   / ,< / __ `/ /_/ // // /_/ / /| | / __ \/ __ \/ / / ___/ __ `/ __/ / __ \/ __
//  / /| / /_/ / __  // // ____/ ___ |/ /_/ / /_/ / / / /__/ /_/ / /_/ / / /_/ / / /
// /_/ |_\__,_/_/ /_/___/_/   /_/  |_/ .___/ .___/_/_/\___/\__,_/\__/_/\____/_/ /_/
//                                  /_/   /_/
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// Project includes
#include "custom_python/add_custom_modeler_to_python.h"
#include "custom_modeler/kahip_partitioning_modeler.h"

namespace Kratos::Python {

void AddCustomModelerToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // ── KaHIPPartitioningModeler ─────────────────────────────────────────────
    // Full-workflow modeler: reads input mesh, partitions with KaHIP, and loads
    // the local partition into a ModelPart in a single SetupModelPart() call.
    py::class_<
        KaHIPPartitioningModeler,
        KaHIPPartitioningModeler::Pointer,
        Modeler>(m, "KaHIPPartitioningModeler")
        .def(py::init<Model&, Parameters>(),
             py::arg("model"), py::arg("parameters"))
        .def("SetupModelPart", &KaHIPPartitioningModeler::SetupModelPart)
        .def("GetDefaultParameters", &KaHIPPartitioningModeler::GetDefaultParameters)
        .def("Info", &KaHIPPartitioningModeler::Info)
        .def("__str__", [](const KaHIPPartitioningModeler& self) {
            return self.Info();
        })
        ;
}

} // namespace Kratos::Python
