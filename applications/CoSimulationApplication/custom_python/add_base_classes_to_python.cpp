//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//

// System includes

// External includes
#include <boost/python.hpp>
#include <string>

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "custom_python/add_base_classes_to_python.h"
#include "base_classes/base_co_simulation_solver_io.h"
#include "base_classes/base_co_simulation_solver.h"
#include "base_classes/base_co_simulation_coupling_strategy.h"
/*#include "custom_base_classes/base_co_simulation_convergence_acceleration_scheme.h"
#include "custom_base_classes/base_co_simulation_convergence_criterion.h"*/

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

namespace Kratos
{

namespace Python
{

void AddCustomBaseClassesToPython()
{
    using namespace boost::python;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

    typedef CoSimulationBaseCouplingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> CoSimulationBaseCouplingStrategyType;
    typedef typename CoSimulationBaseSolver::Pointer CoSimBaseClassPointerType;

    //********************************************************************
    //********************CoSimulationIo**********************************
    //********************************************************************
    class_<CoSimulationBaseIo,
           boost::noncopyable>("CoSimulationBaseIo", init<>())
           .def("ExportData", &CoSimulationBaseIo::ExportData)
           .def("ImportData", &CoSimulationBaseIo::ImportData)
           .def("ExportMesh", &CoSimulationBaseIo::ExportMesh)
           .def("ImportMesh", &CoSimulationBaseIo::ImportMesh) 
           .def("MakeDataAvailable", &CoSimulationBaseIo::MakeDataAvailable)
           .def("MakeMeshAvailable", &CoSimulationBaseIo::MakeDataAvailable);


    //********************************************************************
    //********************CoSimulationSolver*************************
    //********************************************************************
    class_<CoSimulationBaseSolver,boost::noncopyable>("CoSimulationBaseSolver", init<std::string>())
           .def("Name", &CoSimulationBaseSolver::Name)
           .def("SetIo", &CoSimulationBaseSolver::SetIo)
           .def("ImportData", &CoSimulationBaseSolver::ImportData)
           .def("ExportData", &CoSimulationBaseSolver::ExportData)
           .def("MakeDataAvailable", &CoSimulationBaseSolver::MakeDataAvailable)
           .def("ExportMesh", &CoSimulationBaseSolver::ExportMesh)
           .def("ImportMesh", &CoSimulationBaseSolver::ImportMesh) 
           .def("MakeMeshAvailable", &CoSimulationBaseSolver::MakeDataAvailable)
           .def("Initialize", &CoSimulationBaseSolver::ImportMesh) 
           .def("InitializeTimeStep", &CoSimulationBaseSolver::ImportMesh) 
           .def("SolveTimeStep", &CoSimulationBaseSolver::ImportMesh) 
           .def("FinalizeTimeStep", &CoSimulationBaseSolver::ImportMesh)            
           .def("Finalize", &CoSimulationBaseSolver::ImportMesh);

   //********************************************************************
    //********************CoSimulationCouplingStrategy********************
    //******************************************************************** 
    class_<CoSimulationBaseCouplingStrategyType, bases< CoSimulationBaseSolver >, 
           boost::noncopyable>("CoSimulationBaseCouplingStrategy", init<std::string, CoSimBaseClassPointerType , CoSimBaseClassPointerType>())
           .def(init<std::string>());

    //********************************************************************
    //********************CoSimulationRelaxationSchemes*******************
    //********************************************************************
/*     class_<CoSimulationBaseConvergenceAccelerationScheme,
           boost::noncopyable>("CoSimulationBaseConvergenceAccelerationScheme", init<>());

 */
    //********************************************************************
    //********************CoSimulationConvergenceCriterion****************
    //********************************************************************
/*     class_<CoSimulationBaseConvergenceCriterion,
           boost::noncopyable>("CoSimulationBaseConvergenceCriterion", init<double, double>())
           .def("IsConverged",&CoSimulationBaseConvergenceCriterion::IsConverged); */
}

} // namespace Python.

} // Namespace Kratos
