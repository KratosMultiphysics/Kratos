/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
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
 
//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2007-10-24 08:20:31 $
//   Revision:            $Revision: 1.5 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "python/add_utilities_to_python.h"
#include "utilities/variable_utils.h" 
#include "utilities/normal_calculation_utils.h"
#include "utilities/body_normal_calculation_utils.h"
#include "utilities/body_distance_calculation_utils.h"
#include "utilities/signed_distance_calculation_utils.h"
#include "utilities/parallel_levelset_distance_calculator.h"
#include "utilities/openmp_utils.h"

// #include "utilities/signed_distance_calculator_bin_based.h"
#include "utilities/divide_elem_utils.h"
#include "utilities/timer.h"

//#include "spatial_containers/bounding_box.h"
#include "utilities/bounding_box_utilities.h"


namespace Kratos
{
	
namespace Python
{
 
  void  AddUtilitiesToPython()
  {
	using namespace boost::python;

	  class_<VariableUtils>("VariableUtils", init<>())
		.def("SaveVectorVar",&VariableUtils::SaveVectorVar)
		.def("SaveScalarVar",&VariableUtils::SaveScalarVar)
		.def("SelectNodeList",&VariableUtils::SelectNodeList)
		.def("CopyVectorVar",&VariableUtils::CopyVectorVar)
		.def("CopyScalarVar",&VariableUtils::CopyScalarVar)
		.def("SetToZero_VectorVar",&VariableUtils::SetToZero_VectorVar)
		.def("SetToZero_ScalarVar",&VariableUtils::SetToZero_ScalarVar)
		.def("SetToZero_ScalarVar",&VariableUtils::SetToZero_ScalarVar)
		.def("SetToZero_VelocityVectorVar",&VariableUtils::SetToZero_VelocityVectorVar)
		.def("CheckVariableExists",&VariableUtils::SetToZero_VelocityVectorVar)
		;

	  class_<NormalCalculationUtils>("NormalCalculationUtils", init<>())
		.def("CalculateOnSimplex",&NormalCalculationUtils::CalculateOnSimplex)
		;

	  class_<BodyNormalCalculationUtils>("BodyNormalCalculationUtils", init<>())
		.def("CalculateBodyNormals",&BodyNormalCalculationUtils::CalculateBodyNormals)
		;

	  class_<BodyDistanceCalculationUtils>("BodyDistanceCalculationUtils", init<>())
			  .def("CalculateDistances2D",&BodyDistanceCalculationUtils::CalculateDistances<2>)
			  .def("CalculateDistances3D",&BodyDistanceCalculationUtils::CalculateDistances<3>)
			  ;

	  class_<SignedDistanceCalculationUtils<2> >("SignedDistanceCalculationUtils2D", init<>())
			  .def("CalculateDistances",&SignedDistanceCalculationUtils<2>::CalculateDistances )
                          .def("FindMaximumEdgeSize",&SignedDistanceCalculationUtils<2>::FindMaximumEdgeSize )
			  ;

	  class_<SignedDistanceCalculationUtils<3> >("SignedDistanceCalculationUtils3D", init<>())
			  .def("CalculateDistances",&SignedDistanceCalculationUtils<3>::CalculateDistances )
                          .def("FindMaximumEdgeSize",&SignedDistanceCalculationUtils<3>::FindMaximumEdgeSize )
			  ;

	  class_<ParallelDistanceCalculator<2>, boost::noncopyable >("ParallelDistanceCalculator2D", init<>())
			  .def("CalculateDistances",&ParallelDistanceCalculator<2>::CalculateDistances )
                          .def("FindMaximumEdgeSize",&ParallelDistanceCalculator<2>::FindMaximumEdgeSize )
			  ;

	  class_<ParallelDistanceCalculator<3>, boost::noncopyable >("ParallelDistanceCalculator3D", init<>())
			  .def("CalculateDistances",&ParallelDistanceCalculator<3>::CalculateDistances )
                          .def("FindMaximumEdgeSize",&ParallelDistanceCalculator<3>::FindMaximumEdgeSize )
 			  ;
			  
			  
// 	  class_<SignedDistanceCalculationBinBased<2> >("SignedDistanceCalculationBinBased2D", init<>())
// 			  .def("CalculateDistances",&SignedDistanceCalculationBinBased<2>::CalculateDistances )
//                           .def("FindMaximumEdgeSize",&SignedDistanceCalculationBinBased<2>::FindMaximumEdgeSize )
// 			  ;
// 
// 	  class_<SignedDistanceCalculationBinBased<3> >("SignedDistanceCalculationBinBased3D", init<>())
// 			  .def("CalculateDistances",&SignedDistanceCalculationBinBased<3>::CalculateDistances )
//                           .def("FindMaximumEdgeSize",&SignedDistanceCalculationBinBased<3>::FindMaximumEdgeSize )
// 			  ;

	  class_<DivideElemUtils>("DivideElemUtils", init<>())
		.def("DivideElement_2D",&DivideElemUtils::DivideElement_2D)
		;

	  class_<Timer>("Timer", init<>())
		.add_property("PrintOnScreen", &Timer::GetPrintOnScreen,&Timer::SetPrintOnScreen)
		.def("Start", &Timer::Start)
		.def("Stop", &Timer::Stop)
		.staticmethod("Start")
		.staticmethod("Stop")
// 	    .def("PrintTimingInformation",Timer::PrintTimingInformation)
	    .def(self_ns::str(self))
		;
 



          class_<BoundingBoxUtilities>("BoundingBoxUtilities", init<ModelPart&, const unsigned int& >() )
                    .def("Test", &BoundingBoxUtilities::Test)
                    ;

         
//           class_<SplitElements, boost::noncopyable >
//                     ("SplitElements", init<ModelPart&, int >() )
//                     .def("Split", &SplitElements::Split)
//                     ;


// 	  def("PrintTimingInformation",Timer::PrintTimingInformation);

          class_<OpenMPUtils>("OpenMPUtils", init<>() )
                  .def("SetNumThreads", &OpenMPUtils::SetNumThreads)
                  .staticmethod("SetNumThreads")
                  ;

  }
	
}  // namespace Python.

} // Namespace Kratos

