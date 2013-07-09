//   
//   Project Name:        KratosSolidMechanicsApplication $      
//   Last modified by:    $Author:            JMCarbonell $ 
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// External includes 

// Project includes
#include "includes/node.h"
#include "includes/define.h"
#include "processes/process.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

//Application includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_strategies/custom_builders_and_solvers/residual_based_builder_and_solver.hpp"


namespace Kratos
{
	
  namespace Python
  {

    // //Boundary Skin Generator
    // void GenerateSkin(ModelPart& model_part,char* ConditionName,unsigned int dimension,unsigned int preserve )
    // {
    //   GenerateBoundarySkin(model_part,KratosComponents<Condition>::Get(ConditionName),dimension,preserve);
    // }
    

  	
    void  AddCustomUtilitiesToPython()
    {

      using namespace boost::python;

      
      typedef UblasSpace<double, CompressedMatrix, Vector>    SparseSpaceType;
      typedef UblasSpace<double, Matrix, Vector>               LocalSpaceType;
      typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

      typedef Scheme< SparseSpaceType, LocalSpaceType >            SchemeType;
      typedef typename SchemeType::Pointer                  SchemePointerType;

      typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
      typedef typename BuilderAndSolverType::Pointer                   BuilderAndSolverPointerType;
 
      typedef Process                                         ProcessBaseType;




 

    }

  }  // namespace Python.

} // Namespace Kratos

