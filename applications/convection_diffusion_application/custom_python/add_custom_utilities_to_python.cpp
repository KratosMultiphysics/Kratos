/*
==============================================================================
KratosConvectionDiffusionApplication 
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
//   Last modified by:    $Author: antonia $
//   Date:                $Date: 2008-03-11 14:06:15 $
//   Revision:            $Revision: 1.4 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/face_heat_utilities.h"
#include "custom_utilities/pure_convection_tools.h"
#include "custom_utilities/pure_convection_CrankN_tools.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
//#include "custom_utilities/convection_diffusion_settings.h"

namespace Kratos
{
	
namespace Python
{

void GenerateModelPart(FaceHeatUtilities& FaceHeatUtilities,ModelPart& origin_model_part,ModelPart& destination_model_part,unsigned int domain_size )
	{
		if(domain_size == 2)
		{
			FaceHeatUtilities.GenerateModelPart(origin_model_part, destination_model_part, KratosComponents<Element>::Get("ConvDiff2D"),KratosComponents<Condition>::Get("ThermalFace2D")	); 
		}
		else if(domain_size == 3)
		{
			FaceHeatUtilities.GenerateModelPart(origin_model_part, destination_model_part,KratosComponents<Element>::Get("ConvDiff3D"),KratosComponents<Condition>::Get("ThermalFace3D")	); 
		}
	}
	
  void  AddCustomUtilitiesToPython()
  {
	using namespace boost::python;


	  class_<FaceHeatUtilities>("FaceHeatUtilities", init<>())
		.def("ApplyFaceHeat",&FaceHeatUtilities::ApplyFaceHeat)
		.def("ConditionModelPart",&FaceHeatUtilities::ConditionModelPart)
		.def("GenerateModelPart",GenerateModelPart)
		;
	  
	  typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
	  typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
	  typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

	  class_< PureConvectionUtilities< 2, SparseSpaceType, LinearSolverType >,  boost::noncopyable >	("PureConvectionUtilities2D", init<	>() )
			  .def("ConstructSystem",&PureConvectionUtilities< 2, SparseSpaceType, LinearSolverType >::ConstructSystem)
			  .def("CalculateProjection",&PureConvectionUtilities< 2, SparseSpaceType, LinearSolverType >::CalculateProjection)
			  .def("ConvectScalarVar",&PureConvectionUtilities< 2, SparseSpaceType, LinearSolverType >::ConvectScalarVar)
			  .def("ClearSystem",&PureConvectionUtilities< 2, SparseSpaceType, LinearSolverType >::ClearSystem)
			  ;	
	    
	  class_< PureConvectionUtilities< 3, SparseSpaceType, LinearSolverType >,  boost::noncopyable >	("PureConvectionUtilities3D", init<	>() )
			  .def("ConstructSystem",&PureConvectionUtilities< 3, SparseSpaceType, LinearSolverType >::ConstructSystem)
			  .def("CalculateProjection",&PureConvectionUtilities< 3, SparseSpaceType, LinearSolverType >::CalculateProjection)
			  .def("ConvectScalarVar",&PureConvectionUtilities< 3, SparseSpaceType, LinearSolverType >::ConvectScalarVar)
			  .def("ClearSystem",&PureConvectionUtilities< 3, SparseSpaceType, LinearSolverType >::ClearSystem)
			  ;	  

	class_< PureConvectionCrankNUtilities< 2, SparseSpaceType, LinearSolverType >,  boost::noncopyable >	("PureConvectionCrankNUtilities2D", init<	>() )
			  .def("ConstructSystem",&PureConvectionCrankNUtilities< 2, SparseSpaceType, LinearSolverType >::ConstructSystem)
			  .def("CalculateProjection",&PureConvectionCrankNUtilities< 2, SparseSpaceType, LinearSolverType >::CalculateProjection)
			  .def("ConvectScalarVar",&PureConvectionCrankNUtilities< 2, SparseSpaceType, LinearSolverType >::ConvectScalarVar)
			  .def("ClearSystem",&PureConvectionCrankNUtilities< 2, SparseSpaceType, LinearSolverType >::ClearSystem)
			  ;	
	    
	  class_< PureConvectionCrankNUtilities< 3, SparseSpaceType, LinearSolverType >,  boost::noncopyable >	("PureConvectionCrankNUtilities3D", init<	>() )
			  .def("ConstructSystem",&PureConvectionCrankNUtilities< 3, SparseSpaceType, LinearSolverType >::ConstructSystem)
			  .def("CalculateProjection",&PureConvectionCrankNUtilities< 3, SparseSpaceType, LinearSolverType >::CalculateProjection)
			  .def("ConvectScalarVar",&PureConvectionCrankNUtilities< 3, SparseSpaceType, LinearSolverType >::ConvectScalarVar)
			  .def("ClearSystem",&PureConvectionCrankNUtilities< 3, SparseSpaceType, LinearSolverType >::ClearSystem)
			  ;
			  

            

  }
	
}  // namespace Python.

} // Namespace Kratos

