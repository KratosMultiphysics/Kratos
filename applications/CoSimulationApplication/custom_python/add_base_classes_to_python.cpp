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
#include "base_classes/base_co_simulation_application_io.h"
#include "base_classes/base_co_simulation_application.h"
/*#include "custom_base_classes/base_co_simulation_coupling_strategy.h"
#include "custom_base_classes/base_co_simulation_convergence_acceleration_scheme.h"
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

    /*typedef CoSimulationBaseCouplingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> CoSimulationBaseCouplingStrategyType;
    typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> SolvingStrategyType;*/

    //********************************************************************
    //********************CoSimulationIo**********************************
    //********************************************************************
    class_<CoSimulationBaseIo,
           boost::noncopyable>("CoSimulationBaseIo", init<>())
           .def("ExportData", &CoSimulationBaseIo::ExportData)
           .def("ImportData", &CoSimulationBaseIo::ImportData)
           .def("ExportMesh", &CoSimulationBaseIo::ExportMesh)
           .def("ImportMesh", &CoSimulationBaseIo::ImportMesh) 
           .def("MakeDataAvailable", &CoSimulationBaseIo::MakeDataAvailable);


    //********************************************************************
    //********************CoSimulationApplication*************************
    //********************************************************************
    class_<CoSimulationBaseApplication,boost::noncopyable>("CoSimulationBaseApplication", init<std::string>())
           .def("Name", &CoSimulationBaseApplication::Name)
           .def("SetIo", &CoSimulationBaseApplication::SetIo)
           .def("ImportData", &CoSimulationBaseApplication::ImportData)
           .def("ExportData", &CoSimulationBaseApplication::ExportData)
           .def("ExportMesh", &CoSimulationBaseApplication::ExportMesh)
           .def("ImportMesh", &CoSimulationBaseApplication::ImportMesh) 
           .def("MakeDataAvailable", &CoSimulationBaseApplication::MakeDataAvailable);
           
    ;

   /* //********************************************************************
    //********************CoSimulationCouplingStrategy********************
    //********************************************************************
    class_<CoSimulationBaseCouplingStrategyType,
           bases<CoSimulationBaseClassType>,
           boost::noncopyable>("CoSimulationBaseCouplingStrategy", init<CoSimulationBaseApplicationType::Pointer , CoSimulationBaseApplicationType::Pointer, CoSimulationBaseConvergenceCriterion::Pointer, CoSimulationBaseConvergenceAccelerationScheme::Pointer>())
           ;

    //********************************************************************
    //********************CoSimulationRelaxationSchemes*******************
    //********************************************************************
    class_<CoSimulationBaseConvergenceAccelerationScheme,
           boost::noncopyable>("CoSimulationBaseConvergenceAccelerationScheme", init<>());


    //********************************************************************
    //********************CoSimulationConvergenceCriterion****************
    //********************************************************************
    class_<CoSimulationBaseConvergenceCriterion,
           boost::noncopyable>("CoSimulationBaseConvergenceCriterion", init<double, double>())
           .def("IsConverged",&CoSimulationBaseConvergenceCriterion::IsConverged);   */        

}

} // namespace Python.

} // Namespace Kratos
