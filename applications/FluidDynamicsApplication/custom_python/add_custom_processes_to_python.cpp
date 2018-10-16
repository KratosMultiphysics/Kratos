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
#ifdef KRATOS_USE_AMATRIX
#include "boost/numeric/ublas/matrix.hpp" // for the sparse space dense vector
#endif // KRATOS_USE_AMATRIX

// Project includes
#include "containers/model.h"
#include "includes/define_python.h"
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
#include "custom_processes/embedded_skin_visualization_process.h"
#include "custom_processes/integration_point_statistics_process.h"
#include "custom_processes/move_rotor_process.h"
#include "spaces/ublas_space.h"

#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "linear_solvers/linear_solver.h"

namespace Kratos
{

namespace Python
{

void AddCustomProcessesToPython(pybind11::module& m)
{
    using namespace pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double> > SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    class_<SpalartAllmarasTurbulenceModel< SparseSpaceType, LocalSpaceType, LinearSolverType >, SpalartAllmarasTurbulenceModel< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer, Process>
    (m,"SpalartAllmarasTurbulenceModel")
    .def(init < ModelPart&, LinearSolverType::Pointer, unsigned int, double, unsigned int, bool, unsigned int>())
    .def("ActivateDES", &SpalartAllmarasTurbulenceModel< SparseSpaceType, LocalSpaceType, LinearSolverType >::ActivateDES)
    .def("AdaptForFractionalStep", &SpalartAllmarasTurbulenceModel< SparseSpaceType, LocalSpaceType, LinearSolverType >::AdaptForFractionalStep)
    ;

    class_<StokesInitializationProcess< SparseSpaceType, LocalSpaceType, LinearSolverType >, StokesInitializationProcess< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer, Process>
    (m,"StokesInitializationProcess")
    .def(init<ModelPart&, LinearSolverType::Pointer, unsigned int, const Kratos::Variable<int>& >())
    .def("SetConditions",&StokesInitializationProcess<SparseSpaceType, LocalSpaceType, LinearSolverType>::SetConditions)
    ;

    class_<BoussinesqForceProcess, BoussinesqForceProcess::Pointer, Process>
    (m,"BoussinesqForceProcess")
    .def(init<ModelPart&, Parameters& >())
    ;

    class_<WindkesselModel, WindkesselModel::Pointer, Process>
    (m,"WindkesselModel")
    .def(init < ModelPart&>())
    ;

    class_<DistanceModificationProcess, DistanceModificationProcess::Pointer, Process>
    (m,"DistanceModificationProcess")
    .def(init < ModelPart&, const double, const double, const bool, const bool, const bool >())
    .def(init< ModelPart&, Parameters& >())
    ;

    class_<EmbeddedNodesInitializationProcess, EmbeddedNodesInitializationProcess::Pointer, Process>
    (m,"EmbeddedNodesInitializationProcess")
    .def(init < ModelPart&, unsigned int >())
    ;

    class_<EmbeddedPostprocessProcess, EmbeddedPostprocessProcess::Pointer, Process>
    (m,"EmbeddedPostprocessProcess")
    .def(init < ModelPart& >())
    ;

    class_<EmbeddedSkinVisualizationProcess, EmbeddedSkinVisualizationProcess::Pointer, Process>
    (m,"EmbeddedSkinVisualizationProcess")
    .def(init <
        ModelPart&,
        ModelPart&,
        const std::vector<Variable <double> >,
        const std::vector<Variable< array_1d<double, 3> > >,
        const std::vector<VariableComponent<VectorComponentAdaptor< array_1d< double, 3> > > >,
        std::string,
        const bool >())
    .def(init< ModelPart&, ModelPart&, Parameters& >())
    ;

    class_<IntegrationPointStatisticsProcess, IntegrationPointStatisticsProcess::Pointer, Process>
    (m, "IntegrationPointStatisticsProcess")
    .def(init<Model&, Parameters::Pointer>())
    ;

    class_<MoveRotorProcess, MoveRotorProcess::Pointer, Process>
    (m,"MoveRotorProcess")
    .def(init < ModelPart&, const double, const double, const double, const double, const double, const unsigned int >())
    .def(init< ModelPart&, Parameters& >())
    ;

}

} // namespace Python.

} // Namespace Kratos
