//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// External includes 

// Project includes
#include "includes/node.h"
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"

//Application includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_strategies/custom_builders_and_solvers/residual_based_builder_and_solver.hpp"

//Utilities
#include "custom_utilities/boundary_normals_calculation_utilities.hpp"
#include "custom_utilities/modeler_utilities.hpp"
#include "custom_utilities/contact_domain_utilities.hpp"

namespace Kratos
{
	
  namespace Python
  {
    
    void  AddCustomUtilitiesToPython()
    {

      using namespace boost::python;

      
      typedef UblasSpace<double, CompressedMatrix, Vector>    SparseSpaceType;
      typedef UblasSpace<double, Matrix, Vector>               LocalSpaceType;
      typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

      typedef Scheme< SparseSpaceType, LocalSpaceType >            SchemeType;
      typedef SchemeType::Pointer                           SchemePointerType;

      typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType >                 BuilderAndSolverType;
      typedef BuilderAndSolverType::Pointer        BuilderAndSolverPointerType;


      //***************DOMAIN SET**************//
      class_< ModelerUtilities, boost::noncopyable > ("ModelerUtilities", init<>())
	.def("SetDomainLabels",&ModelerUtilities::SetDomainLabels)
	;


      //***************NORMALS**************//

      // This is required to recognize the different overloads 
      typedef  void (BoundaryNormalsCalculationUtilities::*CalcBoundShrinkage)(ModelPart&,int);
      typedef  void (BoundaryNormalsCalculationUtilities::*CalcBoundClassical)(ModelPart&); 

      CalcBoundShrinkage      CalcBoundNormals_Shrinkage     = &BoundaryNormalsCalculationUtilities::CalculateBoundaryNormals;
      CalcBoundClassical      CalcBoundNormals_Classical     = &BoundaryNormalsCalculationUtilities::CalculateBoundaryNormals;
      
      class_<BoundaryNormalsCalculationUtilities > ("BoundaryNormalsCalculation", init<>())
	.def("CalculateBoundaryNormals", CalcBoundNormals_Shrinkage)
	.def("CalculateBoundaryUnitNormals", CalcBoundNormals_Classical)
	;



    }

  }  // namespace Python.

} // Namespace Kratos

