/*
==============================================================================
KratosR1IncompressibleFluidApplication 
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
//   Last modified by:    $Author: jmarti $
//   Date:                $Date: 2009-01-23 14:34:00 $
//   Revision:            $Revision: 1.1 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"



#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/utilities.h" 
#include "custom_utilities/q_utilities.h"

namespace Kratos
{
  
  namespace Python
  {
    
    
    
    void GenerateModelPart(Utils& Utils,ModelPart& origin_model_part,ModelPart& destination_model_part,unsigned int domain_size )
    {
      if(domain_size == 2)
	{
	  Utils.GenerateModelPart(origin_model_part, destination_model_part,
				  KratosComponents<Element>::Get("ConvDiff2D"),//
				  KratosComponents<Condition>::Get("ThermalFace2D")	); //("ConvDiffQ2D" ConvDiff2D
	}
      else if(domain_size == 3)
	{
			Utils.GenerateModelPart(origin_model_part, destination_model_part,
						KratosComponents<Element>::Get("ConvDiff3D"),
						KratosComponents<Condition>::Get("ThermalFace3D")	); 
	}
    }
    
    
    
    
    void  AddCustomUtilitiesToPython()
    {
      using namespace boost::python;
      class_<Utils>("Utils", init<>())
	.def("GenerateModelPart",GenerateModelPart)
	.def("ApplyInitialTemperature",&Utils::ApplyInitialTemperature)
	.def("FindFluidLevel",&Utils::FindFluidLevel)
	;
      
      
      class_<qUtils>("qUtils", init<>())
	.def("EstimateDeltaTime",&qUtils::EstimateDeltaTime)
	.def("IdentifyFluidNodes",&qUtils::IdentifyFluidNodes)
	.def("QuasiLagrangianMove",&qUtils::QuasiLagrangianMove)
	.def("MarkExcessivelyCloseNodes",&qUtils::MarkExcessivelyCloseNodes)
	.def("Predict",&qUtils::Predict)
	.def("EstimateDeltaTime", &qUtils::EstimateDeltaTime)
        .def("CalculateNodalMass",&qUtils::CalculateNodalMass)
	.def("MoveLonelyNodes",&qUtils::MoveLonelyNodes)
        .def("CalculateVolume",&qUtils::CalculateVolume)
	.def("Return",&qUtils::Return)
	.def("ReduceTimeStep",&qUtils::ReduceTimeStep)
	;
    }
    
  }  // namespace Python.
  
} // Namespace Kratos

