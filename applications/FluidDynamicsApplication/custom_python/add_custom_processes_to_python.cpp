//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"

#include "custom_processes/spalart_allmaras_turbulence_model.h"
#include "custom_processes/Boundary_Windkessel_model.h"
#include "custom_processes/stokes_initialization_process.h"
#include "custom_processes/distance_modification_process.h"
#include "custom_processes/boussinesq_force_process.h"
#include "custom_processes/embedded_nodes_initialization_process.h"
#include "custom_processes/embedded_postprocess_process.h"
#include "custom_processes/move_rotor_process.h"
#include "spaces/ublas_space.h"

#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "linear_solvers/linear_solver.h"

namespace Kratos
{

namespace Python
{

void AddCustomProcessesToPython()
{
    using namespace boost::python;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    class_< SpalartAllmarasTurbulenceModel< SparseSpaceType, LocalSpaceType, LinearSolverType >, bases<Process>, boost::noncopyable >
    ("SpalartAllmarasTurbulenceModel", init < ModelPart&, LinearSolverType::Pointer, unsigned int, double, unsigned int, bool, unsigned int>())
    .def("ActivateDES", &SpalartAllmarasTurbulenceModel< SparseSpaceType, LocalSpaceType, LinearSolverType >::ActivateDES)
    .def("AdaptForFractionalStep", &SpalartAllmarasTurbulenceModel< SparseSpaceType, LocalSpaceType, LinearSolverType >::AdaptForFractionalStep)
    ;

    class_< StokesInitializationProcess< SparseSpaceType, LocalSpaceType, LinearSolverType >, bases<Process>, boost::noncopyable >
    ("StokesInitializationProcess",init<ModelPart::Pointer, LinearSolverType::Pointer, unsigned int, const Kratos::Variable<int>& >())
    .def("SetConditions",&StokesInitializationProcess<SparseSpaceType, LocalSpaceType, LinearSolverType>::SetConditions)
    ;

    class_< BoussinesqForceProcess, bases<Process>, boost::noncopyable >
    ("BoussinesqForceProcess",init<ModelPart::Pointer, Parameters& >())
    ;

    class_< WindkesselModel, bases<Process>, boost::noncopyable >
    ("WindkesselModel", init < ModelPart&>())
    ;

    class_< DistanceModificationProcess, bases<Process>, boost::noncopyable >
    ("DistanceModificationProcess",init < ModelPart&, const bool, const bool, const bool >())
    .def(init< ModelPart&, Parameters& >())
    ;

    class_< EmbeddedNodesInitializationProcess, bases<Process>, boost::noncopyable >
    ("EmbeddedNodesInitializationProcess",init < ModelPart&, unsigned int >())
    ;


    class_< EmbeddedPostprocessProcess, bases<Process>, boost::noncopyable >
    ("EmbeddedPostprocessProcess",init < ModelPart& >())
    ;

    class_< MoveRotorProcess, bases<Process>, boost::noncopyable >
    ("MoveRotorProcess",init < ModelPart&, const double, const double, const double, const double, const double, const double >())
    .def(init< ModelPart&, Parameters& >())
    ;

}

} // namespace Python.

} // Namespace Kratos
