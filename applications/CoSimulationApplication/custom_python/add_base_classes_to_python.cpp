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

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "custom_python/add_base_classes_to_python.h"
#include "custom_base_classes/base_co_simulation_class.h"
#include "custom_base_classes/base_co_simulation_application_io.h"
#include "custom_base_classes/base_co_simulation_application.h"
#include "custom_base_classes/base_co_simulation_coupling_strategy.h"
#include "custom_base_classes/base_co_simulation_convergence_acceleration_scheme.h"

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
    typedef CoSimulationBaseClass<SparseSpaceType, LocalSpaceType, LinearSolverType> CoSimulationBaseClassType;
    typedef CoSimulationBaseApplication<SparseSpaceType, LocalSpaceType, LinearSolverType> CoSimulationBaseApplicationType;
    typedef CoSimulationBaseCouplingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> CoSimulationBaseCouplingStrategyType;
    typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> SolvingStrategyType;

    //********************************************************************
    //********************CoSimulationBaseClass***************************
    //********************************************************************
    class_<CoSimulationBaseClassType,
           bases<SolvingStrategyType>,
           boost::noncopyable>("CoSimulationBaseClass", init<ModelPart &>());    

    //********************************************************************
    //********************CoSimulationIo**********************************
    //********************************************************************
    class_<CoSimulationBaseIo,
           boost::noncopyable>("CoSimulationBaseIo", init<Parameters>());

    //********************************************************************
    //********************CoSimulationApplication*************************
    //********************************************************************
    class_<CoSimulationBaseApplicationType,
           bases<CoSimulationBaseClassType>,
           boost::noncopyable>("CoSimulationBaseApplication", init<CoSimulationBaseIo &, Parameters>());

    //********************************************************************
    //********************CoSimulationCouplingStrategy********************
    //********************************************************************
    class_<CoSimulationBaseCouplingStrategyType,
           bases<CoSimulationBaseClassType>,
           boost::noncopyable>("CoSimulationBaseCouplingStrategy", init<CoSimulationBaseApplicationType::Pointer , CoSimulationBaseApplicationType::Pointer, CoSimulationBaseConvergenceAccelerationScheme& >())
           ;

    //********************************************************************
    //********************CoSimulationRelaxationSchemes*******************
    //********************************************************************
    class_<CoSimulationBaseConvergenceAccelerationScheme,
           boost::noncopyable>("CoSimulationBaseConvergenceAccelerationScheme", init<>());
}

} // namespace Python.

} // Namespace Kratos
