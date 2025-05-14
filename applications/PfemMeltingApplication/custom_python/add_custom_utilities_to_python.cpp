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
#include "custom_utilities/pfem2_utilities.h"
#include "custom_utilities/pfemmelting_apply_bc_process.h"
#include "custom_utilities/trial.h"

namespace Kratos
{

namespace Python
{


void  AddCustomUtilitiesToPython(pybind11::module& m)
{
namespace py = pybind11;



 py::class_<Streamline < 3 > >(m,"Streamline").def(py::init<>())
   .def("MovingParticlesN", &Streamline < 3 > ::MovingParticlesN)
   .def("MovingParticles", &Streamline < 3 > ::MovingParticles)
   .def("RungeKutta4ElementbasedSI", &Streamline < 3 > ::RungeKutta4ElementbasedSI)
   ;


/* py::class_<Streamline < 3 > >(m,"Streamline").def(py::init<>())
   .def("SubSteppingElementbasedSI", &Streamline < 3 > ::SubSteppingElementbasedSI)
   .def("RungeKutta4ElementbasedSI", &Streamline < 3 > ::RungeKutta4ElementbasedSI)
   .def("RungeKutta4KernelbasedSI", &Streamline < 3 > ::RungeKutta4KernelbasedSI)
   .def("CheckInvertElement", &Streamline < 3 > ::CheckInvertElement)
   .def("CalculateVolume", &Streamline < 3 > ::CalculateVolume)
   .def("MovingParticles", &Streamline < 3 > ::MovingParticles)
   .def("MovingParticlesN", &Streamline < 3 > ::MovingParticles)
   .def("Aux", &Streamline < 3 > ::MovingParticles)
   ;
*/

 py::class_<FaceHeatFlux < 3 > >(m,"FaceHeatFlux").def(py::init<>())
   .def("FaceHeatFluxDistribution", &FaceHeatFlux < 3 > ::FaceHeatFluxDistribution)
   .def("FlameDistribution", &FaceHeatFlux < 3 > ::FlameDistribution)
   ;

 py::class_<HeatSource < 3 > >(m,"HeatSource").def(py::init<>())
   .def("Heat_Source", &HeatSource < 3 > ::Heat_Source)
   ;

 py::class_<PfemMeltingApplyBCProcess, PfemMeltingApplyBCProcess::Pointer, Process>(m, "PfemMeltingApplyBCProcess").def(py::init<ModelPart &>());
 
 
 py::class_<Pfem2Utils < 2 > >(m,"Pfem2Utils").def(py::init<>())
   .def("MarkExcessivelyCloseNodes", &Pfem2Utils < 2 > ::MarkExcessivelyCloseNodes)
   .def("MarkNodesCloseToWall", &Pfem2Utils < 2 > ::MarkNodesCloseToWall)
   .def("MoveLonelyNodes", &Pfem2Utils < 2 > ::MoveLonelyNodes)
   
   ;
   
   
   

}

}  // namespace Python.

} // Namespace Kratos
