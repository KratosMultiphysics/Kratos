//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

// #define KRATOS_CG_SOLVER_H_EXCLUDED

// System includes


// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/kratos_version.h"
#include "add_vector_to_python.h"
//#include "add_mapped_vector_to_python.h"
#include "add_matrix_to_python.h"
#include "add_points_to_python.h"
#include "add_geometries_to_python.h"
//#include "add_quadratures_to_python.h"
#include "add_containers_to_python.h"
#include "add_processes_to_python.h"
#include "add_model_part_to_python.h"
#include "add_model_to_python.h"
#include "add_io_to_python.h"
#include "add_mesh_to_python.h"
#include "add_modeler_to_python.h"
#include "add_kernel_to_python.h"
#include "add_kratos_application_to_python.h"

#include "add_linear_solvers_to_python.h"
#include "add_amgcl_solver_to_python.h"
#include "add_process_info_to_python.h"
#include "add_constitutive_law_to_python.h"
#include "add_serializer_to_python.h"
#include "add_table_to_python.h"
//#include "add_sparse_vector_to_python.h"
#include "pointer_vector_set_python_interface.h"
#include "solution_step_variable_indexing_python.h"

#include "add_strategies_to_python.h"
#include "add_utilities_to_python.h"

// #include "add_parallel_strategies_to_python.h"
// #include "add_parallel_linear_solvers_to_python.h"

#include "add_matrix_market_interface_to_python.h"

#include "add_search_strategies_to_python.h"
#include "add_kratos_parameters_to_python.h"
#include "add_testing_to_python.h"
#include "add_logger_to_python.h"
#include "add_memory_info_to_python.h"

namespace Kratos
{

namespace Python
{

char const* greet()
{
	std::stringstream header;
	header << "Hello, I am Kratos Multi-Physics " << KRATOS_VERSION <<" ;-)";
    return header.str().c_str();
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
    AddModelToPython();
    AddNodeToPython();
    AddPropertiesToPython();
    AddMeshToPython();
    AddModelerToPython();
    AddKernelToPython();
    AddKratosApplicationToPython();
    AddLinearSolversToPython();
    AddAMGCLSolverToPython();
    AddStrategiesToPython();
    AddUtilitiesToPython();
    AddProcessInfoToPython();
    AddConstitutiveLawToPython();
    AddSerializerToPython();
    AddTableToPython();
    AddGeometriesToPython();
//     AddParallelStrategiesToPython(); 
//     AddParallelLinearSolversToPython();
    AddMatrixMarketInterfaceToPython();
    AddKratosParametersToPython();

    AddSearchStrategiesToPython();
    AddTestingToPython();
    AddLoggerToPython();
    AddMemoryInfoToPython();

    def("Hello", greet);
}


}  // namespace Python.

}  // namespace Kratos.
