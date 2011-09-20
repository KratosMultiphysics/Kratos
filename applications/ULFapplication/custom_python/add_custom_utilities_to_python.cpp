/*
==============================================================================
KratosULFApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).
 
Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Pawel Ryzhakov
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
//   Last modified by:    $Author: anonymous $
//   Date:                $Date: 2009-01-15 14:50:24 $
//   Revision:            $Revision: 1.6 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/ulf_utilities.h"
#include "custom_utilities/nist_utilities.h"

namespace Kratos
{
	
namespace Python
{
	void GenerateModelPart(NistUtils& NistUtils,ModelPart& origin_model_part,ModelPart& destination_model_part,unsigned int domain_size )
	{
		if(domain_size == 2)
		{
			NistUtils.GenerateModelPart(origin_model_part, destination_model_part,
				KratosComponents<Element>::Get("ConvDiff2D"),
				KratosComponents<Condition>::Get("ThermalFace2D")	); 
		}
		else if(domain_size == 3)
		{
			NistUtils.GenerateModelPart(origin_model_part, destination_model_part,
				KratosComponents<Element>::Get("ConvDiff3D"),
				KratosComponents<Condition>::Get("ThermalFace3D")	); 
		}
	}
	
  void  AddCustomUtilitiesToPython()
  {
	using namespace boost::python;

	  class_<UlfUtils>("UlfUtils", init<>())
		.def("ApplyBoundaryConditions",&UlfUtils::ApplyBoundaryConditions)
		.def("ApplyMinimalPressureConditions",&UlfUtils::ApplyMinimalPressureConditions)
		.def("EstimateDeltaTime",&UlfUtils::EstimateDeltaTime)
		.def("MarkOuterNodes",&UlfUtils::MarkOuterNodes)
		.def("Predict",&UlfUtils::Predict)
		.def("MoveLonelyNodes",&UlfUtils::MoveLonelyNodes)
		.def("CalculateVolume",&UlfUtils::CalculateVolume)
		.def("CalculateNodalArea",&UlfUtils::CalculateNodalArea)
		.def("IdentifyFluidNodes",&UlfUtils::IdentifyFluidNodes)
		.def("ReduceTimeStep",&UlfUtils::ReduceTimeStep)
		.def("MarkExcessivelyCloseNodes",&UlfUtils::MarkExcessivelyCloseNodes)
		.def("NodalIncrementalPressureCalculation",&UlfUtils::NodalIncrementalPressureCalculation)
		.def("SaveNodalArea",&UlfUtils::SaveNodalArea)
		.def("SaveLagrangianInlet", &UlfUtils::SaveLagrangianInlet)
		.def("InjectNodesAtInlet", &UlfUtils::InjectNodesAtInlet)
		.def("MoveInletNodes", &UlfUtils::MoveInletNodes)
		.def("MarkNodesCloseToWall", &UlfUtils::MarkNodesCloseToWall)
		.def("MarkNodesTouchingWall", &UlfUtils::MarkNodesTouchingWall)
		.def("MarkNodesCloseToWallForBladder", &UlfUtils::MarkNodesCloseToWallForBladder)
		.def("MarkNodesCloseToFS", &UlfUtils::MarkNodesCloseToFS)
		.def ("MarkLonelyNodesForErasing", &UlfUtils::MarkLonelyNodesForErasing)
		.def ("RestoreLagInletNodalH", &UlfUtils::RestoreNodalHAtLagInlet)
		.def ("SetLagInletNodalH", &UlfUtils::SetNodalHAtLagInlet)
		.def ("DeleteFreeSurfaceNodesBladder", &UlfUtils::DeleteFreeSurfaceNodesBladder)		
		 ;
	  
	  class_<NistUtils>("NistUtils", init<>())
		.def("GenerateModelPart",GenerateModelPart)
		.def("ApplyInitialTemperature",&NistUtils::ApplyInitialTemperature)
		.def("FindFluidLevel",&NistUtils::FindFluidLevel)
		;	
	
  }
	
}  // namespace Python.

} // Namespace Kratos

