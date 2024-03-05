//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

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
#include "add_cfd_variables_to_python.h"
#include "add_mesh_moving_variables_to_python.h"
#include "add_mapping_variables_to_python.h"
#include "add_dem_variables_to_python.h"
#include "add_fsi_variables_to_python.h"
#include "add_mat_variables_to_python.h"
#include "add_containers_to_python.h"
#include "add_matrix_to_python.h"
#include "add_quaternion_to_python.h"
#include "add_points_to_python.h"
#include "add_predicates_to_python.h"
#include "add_geometries_to_python.h"
#include "add_bounding_box_to_python.h"
#include "add_containers_to_python.h"
#include "add_controllers_to_python.h"
#include "add_operations_to_python.h"
#include "add_processes_to_python.h"
#include "add_properties_to_python.h"
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
#include "add_convergence_accelerators_to_python.h"
#include "add_geometrical_utilities_to_python.h"
#include "add_other_utilities_to_python.h"
#include "add_variable_utils_to_python.h"
#include "add_matrix_market_interface_to_python.h"
#include "add_search_strategies_to_python.h"
#include "add_kratos_parameters_to_python.h"
#include "add_testing_to_python.h"
#include "add_logger_to_python.h"
#include "add_memory_info_to_python.h"
#include "add_constraint_to_python.h"
#include "add_response_functions_to_python.h"
#include "add_communicator_to_python.h"
#include "add_data_communicator_to_python.h"
#include "add_parallel_environment_to_python.h"
#include "add_global_pointers_to_python.h"
#include "add_dofs_to_python.h"
#include "add_mapper_to_python.h"
#include "add_sparse_matrices_to_python.h"
#include "add_registry_to_python.h"
#include "add_container_expression_to_python.h"
#include "add_accessors_to_python.h"
#include "add_globals_to_python.h"
#include "add_geometry_data_to_python.h"

namespace Kratos::Python
{

std::string Hello()
{
    std::stringstream header;
    header << "Hello, I am Kratos Multi-Physics " << GetVersionString() << "for" << GetOSName() << " ;-)\n";
    return header.str();
}

PYBIND11_MODULE(Kratos, m)
{
    namespace py = pybind11;

    AddVectorToPython(m);
    AddMatrixToPython(m);
    AddPointsToPython(m);
    AddPredicatesToPython(m);
    AddKernelToPython(m);
    AddContainersToPython(m);
    AddModelPartToPython(m);
    AddDofsToPython(m);
    AddNodeToPython(m);
    AddPropertiesToPython(m);
    AddAccessorsToPython(m);
    AddMeshToPython(m);
    AddQuaternionToPython(m);

    AddCFDVariablesToPython(m);
    AddDEMVariablesToPython(m);
    AddALEVariablesToPython(m);
    AddFSIVariablesToPython(m);
    AddMappingVariablesToPython(m);
    AddMATVariablesToPython(m);
    AddDeprecatedVariablesToPython(m);
    AddGlobalPointersToPython(m);

    AddOperationsToPython(m);
    AddControllersToPython(m);
    AddProcessesToPython(m);
    AddIOToPython(m);
    AddModelToPython(m);

    AddModelerToPython(m);

    AddKratosApplicationToPython(m);
    AddLinearSolversToPython(m);
    AddFactoriesToPython(m);
    AddAMGCLSolverToPython(m);
    AddStrategiesToPython(m);
    AddConvergenceAcceleratorsToPython(m);
    AddGeometricalUtilitiesToPython(m);
    AddOtherUtilitiesToPython(m);
    AddVariableUtilsToPython(m);
    AddProcessInfoToPython(m);
    AddConstitutiveLawToPython(m);
    AddSerializerToPython(m);
    AddTableToPython(m);
    AddGeometryDataToPython(m);
    AddGeometriesToPython(m);
    AddBoundingBoxToPython(m);

    AddMatrixMarketInterfaceToPython(m);
    AddKratosParametersToPython(m);

    AddSearchStrategiesToPython(m);
    AddTestingToPython(m);
    AddLoggerToPython(m);
    AddMemoryInfoToPython(m);
    AddConstraintToPython(m);
    AddResponseFunctionsToPython(m);
    AddCommunicatorToPython(m);
    AddDataCommunicatorToPython(m);
    AddParallelEnvironmentToPython(m);
    AddMapperToPython(m);
    AddSparseMatricesToPython(m);
    AddRegistryToPython(m);
    AddContainerExpressionToPython(m);
    AddGlobalsToPython(m);

    m.def("Hello", Hello);
}

}  // namespace Kratos::Python.
