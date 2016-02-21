// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
#include "add_io_to_python.h"
#include "add_mesh_to_python.h"
#include "add_modeler_to_python.h"
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

#include "add_search_strategies_to_python.h"
#include "add_kratos_parameters_to_python.h"

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
    AddNodeToPython();
    AddPropertiesToPython();
    AddMeshToPython();
    AddModelerToPython();
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
    AddGeometriesToPython();
    AddParallelStrategiesToPython();
    AddParallelLinearSolversToPython();
    AddMatrixMarketInterfaceToPython();
    AddKratosParametersToPython();
    
    AddSearchStrategiesToPython();

    def("Hello", greet);
}


}  // namespace Python.

}  // namespace Kratos.


