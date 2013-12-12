/*
==============================================================================
KratosTestApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
 
//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: G.Casas$
//   Date:                $Date: 2011-6-27 16:07:33 $
//   Revision:            $Revision: 1.2 $
//
//

// System includes 

// External includes 
#include <boost/python.hpp>

// Project includes

//#include "includes/model_part.h" -S
//#include "custom_python/add_custom_utilities_to_python.h" -S
#include "add_custom_utilities_to_python.h" //+S
//#include "custom_python/add_custom_utilities_to_python.h" -S
//#include "custom_utilities/create_and_destroy.h"

#include "custom_utilities/custom_functions.h"
#include "custom_utilities/binbased_DEM_fluid_coupled_mapping.h" //S
#include "custom_utilities/volume_averaging_tool.h"


namespace Kratos{

namespace Python{
    
typedef ModelPart::NodesContainerType::iterator PointIterator;
typedef std::vector<array_1d<double, 3 > > ComponentVectorType;
typedef std::vector<array_1d<double, 3 > >::iterator ComponentIteratorType;

void  AddCustomUtilitiesToPython(){
using namespace boost::python;
    
    class_<CustomFunctionsCalculator, boost::noncopyable >
        ("CustomFunctionsCalculator", init<>())
        .def("CalculatePressureGradient", &CustomFunctionsCalculator::CalculatePressureGradient)
        .def("AssessStationarity", &CustomFunctionsCalculator::AssessStationarity)
        ; 
    
    class_<BinBasedDEMFluidCoupledMapping < 2 > >("BinBasedDEMFluidCoupledMapping2D", init<double, int, int>())
        .def("InterpolateFromFluidMesh", &BinBasedDEMFluidCoupledMapping < 2 > ::InterpolateFromFluidMesh)
        .def("InterpolateFromDEMMesh", &BinBasedDEMFluidCoupledMapping < 2 > ::InterpolateFromDEMMesh)
        .def("ComputePostProcessResults", &BinBasedDEMFluidCoupledMapping < 2 > ::ComputePostProcessResults)
        ;
 
    class_<BinBasedDEMFluidCoupledMapping < 3 > >("BinBasedDEMFluidCoupledMapping3D", init<double, int>())
        .def("InterpolateFromFluidMesh", &BinBasedDEMFluidCoupledMapping < 3 > ::InterpolateFromFluidMesh)
        .def("InterpolateFromNewestFluidMesh", &BinBasedDEMFluidCoupledMapping < 3 > ::InterpolateFromNewestFluidMesh)
        .def("InterpolateFromDEMMesh", &BinBasedDEMFluidCoupledMapping < 3 > ::InterpolateFromDEMMesh)
        .def("ComputePostProcessResults", &BinBasedDEMFluidCoupledMapping < 3 > ::ComputePostProcessResults)
        ;
    }

}  // namespace Python.

} // Namespace Kratos
