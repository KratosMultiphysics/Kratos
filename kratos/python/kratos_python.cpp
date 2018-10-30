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
#include <pybind11/pybind11.h>



// Project includes
#include "includes/define_python.h"
#include "includes/kratos_version.h"
#include "add_vector_to_python.h"
#include "add_model_part_to_python.h"
#include "add_kernel_to_python.h"
#include "add_deprecated_variables_to_python.h"
#include "add_c2c_variables_to_python.h"
#include "add_cfd_variables_to_python.h"
#include "add_mesh_moving_variables_to_python.h"
#include "add_mapping_variables_to_python.h"
#include "add_dem_variables_to_python.h"
#include "add_fsi_variables_to_python.h"
#include "add_mat_variables_to_python.h"
#include "add_legacy_structural_app_vars_to_python.h"
#include "add_containers_to_python.h"
#include "add_matrix_to_python.h"
#include "add_points_to_python.h"
#include "add_geometries_to_python.h"
#include "add_containers_to_python.h"
#include "add_processes_to_python.h"
#include "add_model_to_python.h"
#include "add_io_to_python.h"
#include "add_mesh_to_python.h"
#include "add_modeler_to_python.h"
#include "add_kratos_application_to_python.h"
#include "add_linear_solvers_to_python.h"
#include "add_factories_to_python.h"
#include "add_amgcl_solver_to_python.h"
#include "add_process_info_to_python.h"
#include "add_constitutive_law_to_python.h"
#include "add_serializer_to_python.h"
#include "add_table_to_python.h"
#include "add_strategies_to_python.h"
#include "add_utilities_to_python.h"
#include "add_matrix_market_interface_to_python.h"
#include "add_search_strategies_to_python.h"
#include "add_kratos_parameters_to_python.h"
#include "add_testing_to_python.h"
#include "add_logger_to_python.h"
#include "add_memory_info_to_python.h"
#include "add_constraint_to_python.h"
#include "add_response_functions_to_python.h"
#include "add_data_communicator_to_python.h"

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

PYBIND11_MODULE(Kratos, m)
{
    namespace py = pybind11;

    AddVectorToPython(m);
    AddMatrixToPython(m);
    AddPointsToPython(m);
    AddKernelToPython(m);
    AddContainersToPython(m);
    AddModelPartToPython(m);
    AddNodeToPython(m);
    AddPropertiesToPython(m);
    AddMeshToPython(m);

    AddCFDVariablesToPython(m);
    AddDEMVariablesToPython(m);
    AddALEVariablesToPython(m);
    AddC2CVariablesToPython(m);
    AddFSIVariablesToPython(m);
    AddLegacyStructuralAppVarsToPython(m);
    AddMappingVariablesToPython(m);
    AddMATVariablesToPython(m);
    AddDeprecatedVariablesToPython(m);



//     //  AddQuadraturesToPython();
//     //  AddIntegrationPointsToPython();

     AddProcessesToPython(m);
     AddIOToPython(m);
     AddModelToPython(m);

     AddModelerToPython(m);

     AddKratosApplicationToPython(m);
     AddLinearSolversToPython(m);
     AddFactoriesToPython(m);
     AddAMGCLSolverToPython(m);
     AddStrategiesToPython(m);
     AddUtilitiesToPython(m);
     AddProcessInfoToPython(m);
     AddConstitutiveLawToPython(m);
     AddSerializerToPython(m);
     AddTableToPython(m);
     AddGeometriesToPython(m);

    AddMatrixMarketInterfaceToPython(m);
    AddKratosParametersToPython(m);
//
    AddSearchStrategiesToPython(m);
     AddTestingToPython(m);
     AddLoggerToPython(m); //TO BE SPOKEN WITH POOYAN
     AddConstraintToPython(m);
     AddResponseFunctionsToPython(m);
     AddDataCommunicatorToPython(m);

     //     AddBandedMatrixToPython();
     //     AddTriangularMatrixToPython();
     //     AddSymmetricMatrixToPython();
     //     AddIdentityMatrixToPython();
     //     AddZeroMatrixToPython();
     //     AddScalarMatrixToPython();
     //     AddSparseMatrixToPython();
     //     AddCompressedMatrixToPython();
     // #if defined KRATOS_ADD_COORDINATE_MATRIX_INTERFACE
     //     AddCoordinateMatrixToPython();
     // #endif

     m.def("Hello", greet);
}


}  // namespace Python.

}  // namespace Kratos.
