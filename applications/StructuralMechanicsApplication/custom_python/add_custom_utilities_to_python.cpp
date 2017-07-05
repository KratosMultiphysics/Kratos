// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

//Utilities
#include "custom_utilities/sprism_neighbours.hpp"
#include "custom_utilities/eigenvector_to_solution_step_variable_transfer_utility.hpp"

//Processes
#include "custom_processes/apply_multi_point_constraints_process.h"


namespace Kratos
{
namespace Python
{

inline
void TransferEigenvector1(
        EigenvectorToSolutionStepVariableTransferUtility& rThisUtil,
        ModelPart& rModelPart,
        int iEigenMode)
{
    rThisUtil.Transfer(rModelPart,iEigenMode);
}

inline
void TransferEigenvector2(
        EigenvectorToSolutionStepVariableTransferUtility& rThisUtil,
        ModelPart& rModelPart,
        int iEigenMode,
        int step)
{
    rThisUtil.Transfer(rModelPart,iEigenMode,step);
}

void  AddCustomUtilitiesToPython()
{
    using namespace boost::python;

//     typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
//     typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
//     typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    
    class_<SprismNeighbours>("SprismNeighbours", init<ModelPart&>())
    .def("Execute",&SprismNeighbours::Execute)
    .def("ClearNeighbours",&SprismNeighbours::ClearNeighbours)
    ;

    class_<EigenvectorToSolutionStepVariableTransferUtility>(
                "EigenvectorToSolutionStepVariableTransferUtility")
    .def("Transfer",TransferEigenvector1)
    .def("Transfer",TransferEigenvector2)
    ;

    /// Processes
    class_<ApplyMultipointConstraintsProcess, boost::noncopyable, bases<Process>>("ApplyMultipointConstraintsProcess", init<ModelPart&>())
    .def(init< ModelPart&, Parameters& >())
	.def("AddMasterSlaveRelation", &ApplyMultipointConstraintsProcess::AddMasterSlaveRelationWithNodesAndVariableComponents)
    .def("AddMasterSlaveRelation", &ApplyMultipointConstraintsProcess::AddMasterSlaveRelationWithNodeIdsAndVariableComponents)
	.def("AddMasterSlaveRelation", &ApplyMultipointConstraintsProcess::AddMasterSlaveRelationWithNodesAndVariable)
    .def("AddMasterSlaveRelation", &ApplyMultipointConstraintsProcess::AddMasterSlaveRelationWithNodeIdsAndVariable)
    .def("SetActive", &ApplyMultipointConstraintsProcess::SetActive)      
    .def("PrintData", &ApplyMultipointConstraintsProcess::PrintData);

}

}  // namespace Python.  

} // Namespace Kratos

