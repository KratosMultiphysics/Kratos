// KRATOS 
// _____   __               __  __      _ _   _             
//|  __ \ / _|             |  \/  |    | | | (_)            
//| |__) | |_ ___ _ __ ___ | \  / | ___| | |_ _ _ __   __ _ 
//|  ___/|  _/ _ \ '_ ` _ \| |\/| |/ _ \ | __| | '_ \ / _` |
//| |    | ||  __/ | | | | | |  | |  __/ | |_| | | | | (_| |
//|_|    |_| \___|_| |_| |_|_|  |_|\___|_|\__|_|_| |_|\__, |
//                                                     __/ |
//                                                    |___/ APPLICATION
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Julio Marti
//



// System includes

// External includes
#include <pybind11/pybind11.h>


// Project includes
//#include "includes/define.h"
#include "includes/define_python.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "custom_utilities/face_heat_distribution.h"
#include "custom_utilities/streamline.h"
#include "custom_utilities/heat_source.h"
namespace Kratos
{

namespace Python
{


void  AddCustomUtilitiesToPython(pybind11::module& m)
{
namespace py = pybind11;



 py::class_<Streamline < 3 > >(m,"Streamline").def(py::init<>())
   .def("SubSteppingElementbasedSI", &Streamline < 3 > ::SubSteppingElementbasedSI)
   .def("RungeKutta4ElementbasedSI", &Streamline < 3 > ::RungeKutta4ElementbasedSI)
   .def("RungeKutta4KernelbasedSI", &Streamline < 3 > ::RungeKutta4KernelbasedSI)
   .def("CheckInvertElement", &Streamline < 3 > ::CheckInvertElement)
   ;


 py::class_<FaceHeatFlux < 3 > >(m,"FaceHeatFlux").def(py::init<>())
   .def("FaceHeatFluxDistribution", &FaceHeatFlux < 3 > ::FaceHeatFluxDistribution)
   .def("FlameDistribution", &FaceHeatFlux < 3 > ::FlameDistribution)
   ;

 py::class_<HeatSource < 3 > >(m,"HeatSource").def(py::init<>())
   .def("Heat_Source", &HeatSource < 3 > ::Heat_Source)
   ;


}

}  // namespace Python.

} // Namespace Kratos
