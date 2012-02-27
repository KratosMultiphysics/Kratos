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
//   Last modified by:    $Author: janosch $
//   Date:                $Date: 2008-03-05 14:35:19 $
//   Revision:            $Revision: 1.8 $
//
//
// #define KRATOS_CG_SOLVER_H_EXCLUDED

// System includes 


// External includes 
#include <boost/python.hpp>


// Project includes 
#include "includes/define.h"
#include "add_vector_to_python.h" 
//#include "add_mapped_vector_to_python.h"
#include "add_matrix_to_python.h" 
#include "add_geometries_to_python.h" 
//#include "add_quadratures_to_python.h" 
#include "add_containers_to_python.h" 
#include "add_processes_to_python.h" 
#include "add_model_part_to_python.h" 
#include "add_io_to_python.h" 
#include "add_mesh_to_python.h" 
#include "add_kernel_to_python.h" 
#include "add_kratos_application_to_python.h" 
//#include "add_equation_systems_to_python.h"
#include "add_linear_solvers_to_python.h" 
#include "add_process_info_to_python.h"
#include "add_constitutive_law_to_python.h"
#include "add_serializer_to_python.h"
#include "add_table_to_python.h"

//#include "add_sparse_vector_to_python.h"
#include "pointer_vector_set_python_interface.h"
#include "solution_step_variable_indexing_python.h"

#include "add_linear_solvers_to_python.h" 
#include "add_strategies_to_python.h" 
#include "add_utilities_to_python.h" 

#include "add_parallel_strategies_to_python.h"
#include "add_parallel_linear_solvers_to_python.h"

#include "add_matrix_market_interface_to_python.h"

namespace Kratos
{ 

namespace Python 
{ 

  char const* greet()
  {
    return "Hello, I am Kratos new release :-)";
  }
  
  using namespace boost::python;
  
  BOOST_PYTHON_MODULE(Kratos)
  {
    AddVectorToPython();     
//    AddSparseVectorToPython();
    AddMatrixToPython();
    AddBandedMatrixToPython();
    AddTriangularMatrixToPython();
    AddSymmetricMatrixToPython();
#if defined KRATOS_ADD_HERMITIAN_MATRIX_INTERFACE
    AddHermitianMatrixToPython();
#endif
    AddIdentityMatrixToPython();
    AddZeroMatrixToPython();
    AddScalarMatrixToPython();
    AddSparseMatrixToPython();
    AddCompressedMatrixToPython();
#if defined KRATOS_ADD_COORDINATE_MATRIX_INTERFACE
    AddCoordinateMatrixToPython();
#endif
    AddPointsToPython();
  //  AddQuadraturesToPython();
  //  AddIntegrationPointsToPython();
    AddContainersToPython();
	AddProcessesToPython();
	AddIOToPython();
	AddModelPartToPython();
	AddNodeToPython();
	AddPropertiesToPython();
	AddMeshToPython();
	AddKernelToPython();
	AddKratosApplicationToPython();
//	AddEquationSystemsToPython();
	AddLinearSolversToPython();
	AddStrategiesToPython();
	AddUtilitiesToPython();
	AddProcessInfoToPython();
	AddConstitutiveLawToPython();
        AddSerializerToPython();
        AddTableToPython();

	AddParallelStrategiesToPython();
	AddParallelLinearSolversToPython();
	AddMatrixMarketInterfaceToPython();

    def("Hello", greet);
  }
  
  
}  // namespace Python.
  
}  // namespace Kratos.


