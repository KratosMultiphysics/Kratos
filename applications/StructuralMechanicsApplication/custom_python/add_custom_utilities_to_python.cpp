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
#include "custom_utilities/formulate_constraint_utility.hpp"

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

    class_<MpcData, boost::noncopyable>( "MpcData", init<>() )
	.def("AddSlave",&MpcData::AddSlave)
	.def("AddMaster",&MpcData::AddMaster)
    .def(self_ns::str(self))
    ;
    class_<Variable<MpcData> , bases<Flags>, boost::noncopyable >("MpcData", no_init)
    ;

    class_<AddMultiPointConstraint>("AddMultiPointConstraint")
	.def("ApplyConstraint", &AddMultiPointConstraint::ApplyConstraint);

}

}  // namespace Python.

} // Namespace Kratos

