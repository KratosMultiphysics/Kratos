/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
mossaiby@yahoo.com
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"

#include "custom_python/add_utilities_to_python.h"

// #if defined(__APPLE__) || defined(__MACOSX)
// 	#include <OpenCL/cl.h>
// #else
// 	#include <CL/cl.h>
// #endif

#include "includes/model_part.h"

#include "custom_utilities/opencl_interface.h"
#include "custom_utilities/opencl_edge_data.h"
#include "custom_utilities/opencl_pure_convection_edgebased.h"
#include "custom_utilities/opencl_fluid_solver.h"
#include "custom_utilities/renumber_by_neighbour_count_util.h"


namespace Kratos
{

namespace Python
{
    void  AddUtilitiesToPython()
    {

        using namespace boost::python;

          class_< OpenCL::DeviceGroup, boost::noncopyable > ("OpenCLDeviceGroup", init<cl_device_type, bool> ())
						  .def("AddCLSearchPath", &OpenCL::DeviceGroup::AddCLSearchPath)
						;

          class_< OpenCLMatrixContainer,  boost::noncopyable >   ("OpenCLMatrixContainer3D", init< OpenCL::DeviceGroup& >() )
                          .def("ConstructCSRVector",&OpenCLMatrixContainer::ConstructCSRVector)
                          .def("BuildCSRData",&OpenCLMatrixContainer::BuildCSRData)
                          .def("Clear",&OpenCLMatrixContainer::Clear)
                        ;

          class_< OpenCLPureConvectionEdgeBased3D,  boost::noncopyable >   ("OpenCLPureConvectionEdgeBased3D", init< OpenCLMatrixContainer&, ModelPart& >() )
                          .def("Initialize", &OpenCLPureConvectionEdgeBased3D::Initialize)
                          .def("Solve", &OpenCLPureConvectionEdgeBased3D::Solve)
						  .def("ComputeTimeStep", &OpenCLPureConvectionEdgeBased3D::ComputeTimeStep)
                        ;

          class_< OpenCLFluidSolver3D,  boost::noncopyable >   ("OpenCLFluidSolver3D", init< OpenCLMatrixContainer&, ModelPart&, const double, const double, const Vector, bool, double, double, double, double, bool >() )
						  .def("Initialize", &OpenCLFluidSolver3D::Initialize)
						  .def("ComputeTimeStep", &OpenCLFluidSolver3D::ComputeTimeStep)
						  .def("SolveStep1", &OpenCLFluidSolver3D::SolveStep1)
						  .def("SolveStep2", &OpenCLFluidSolver3D::SolveStep2)
						  .def("SolveStep3", &OpenCLFluidSolver3D::SolveStep3)
						  .def("ComputeTimeStep", &OpenCLFluidSolver3D::ComputeTimeStep)
						  .def("CalculateNormals", &OpenCLFluidSolver3D::CalculateNormals)
						  .def("UpdateFixedVelocityValues", &OpenCLFluidSolver3D::UpdateFixedVelocityValues)
                                                  .def("ActivateWallResistance", &OpenCLFluidSolver3D::ActivateWallResistance)
                        ;

          enum_<cl_device_type>("cl_device_type")
                          .value("CL_DEVICE_TYPE_CPU", CL_DEVICE_TYPE_CPU)
                          .value("CL_DEVICE_TYPE_GPU", CL_DEVICE_TYPE_GPU)
                          .value("CL_DEVICE_TYPE_ALL", CL_DEVICE_TYPE_ALL)
                        ;

          class_< RenumberByNeighbourCountUtil,  boost::noncopyable >   ("RenumberByNeighbourCountUtil", init< >() )
                          .def("Renumber",&RenumberByNeighbourCountUtil::Renumber)
                        ;


  }

}  // namespace Python.

} // Namespace Kratos

