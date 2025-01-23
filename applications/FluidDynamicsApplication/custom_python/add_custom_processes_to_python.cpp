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

// Project includes
#include "containers/model.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "includes/define_python.h"
#include "includes/model_part.h"
#include "processes/process.h"

#include "custom_processes/apply_compressible_navier_stokes_boundary_conditions_process.h"
#include "custom_processes/Boundary_Windkessel_model.h"
#include "custom_processes/boussinesq_force_process.h"
#include "custom_processes/calulate_levelset_consistent_nodal_gradient_process.h"
#include "custom_processes/compute_pressure_coefficient_process.h"
#include "custom_processes/distance_modification_process.h"
#include "custom_processes/distance_smoothing_process.h"
#include "custom_processes/embedded_nodes_initialization_process.h"
#include "custom_processes/embedded_postprocess_process.h"
#include "custom_processes/embedded_skin_visualization_process.h"
#include "custom_processes/integration_point_statistics_process.h"
#include "custom_processes/mass_conservation_check_process.h"
#include "custom_processes/two_fluids_inlet_process.h"
#include "custom_processes/shock_capturing_entropy_viscosity_process.h"
#include "custom_processes/shock_capturing_physics_based_process.h"
#include "custom_processes/spalart_allmaras_turbulence_model.h"
#include "custom_processes/stokes_initialization_process.h"
#include "custom_processes/compute_y_plus_process.h"
#include "custom_processes/navier_stokes_vectorial_fractional_convection_process.h"

#include "spaces/ublas_space.h"

#include "linear_solvers/linear_solver.h"
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"


namespace Kratos
{

namespace Python
{

void AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double> > SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    py::class_<SpalartAllmarasTurbulenceModel< SparseSpaceType, LocalSpaceType, LinearSolverType >, SpalartAllmarasTurbulenceModel< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer, Process>
    (m,"SpalartAllmarasTurbulenceModel")
    .def(py::init < ModelPart&, LinearSolverType::Pointer, unsigned int, double, unsigned int, bool, unsigned int>())
    .def("ActivateDES", &SpalartAllmarasTurbulenceModel< SparseSpaceType, LocalSpaceType, LinearSolverType >::ActivateDES)
    .def("AdaptForFractionalStep", &SpalartAllmarasTurbulenceModel< SparseSpaceType, LocalSpaceType, LinearSolverType >::AdaptForFractionalStep)
    ;

    py::class_<StokesInitializationProcess< SparseSpaceType, LocalSpaceType, LinearSolverType >, StokesInitializationProcess< SparseSpaceType, LocalSpaceType, LinearSolverType >::Pointer, Process>
    (m,"StokesInitializationProcess")
    .def(py::init<ModelPart&, LinearSolverType::Pointer, unsigned int, const Kratos::Variable<int>& >())
    .def("SetConditions",&StokesInitializationProcess<SparseSpaceType, LocalSpaceType, LinearSolverType>::SetConditions)
    ;

    py::class_<BoussinesqForceProcess, BoussinesqForceProcess::Pointer, Process>
    (m,"BoussinesqForceProcess")
    .def(py::init<ModelPart&, Parameters& >())
    ;

    py::class_<WindkesselModel, WindkesselModel::Pointer, Process>
    (m,"WindkesselModel")
    .def(py::init < ModelPart&>())
    ;

    py::class_<DistanceModificationProcess, DistanceModificationProcess::Pointer, Process>
    (m,"DistanceModificationProcess")
    .def(py::init < ModelPart&, const double, const double, const bool, const bool, const bool >())
    .def(py::init< ModelPart&, Parameters& >())
    .def(py::init< Model&, Parameters& >())
    ;

    py::class_<EmbeddedNodesInitializationProcess, EmbeddedNodesInitializationProcess::Pointer, Process>
    (m,"EmbeddedNodesInitializationProcess")
    .def(py::init < ModelPart&, unsigned int >())
    ;

    py::class_<EmbeddedPostprocessProcess, EmbeddedPostprocessProcess::Pointer, Process>
    (m,"EmbeddedPostprocessProcess")
    .def(py::init < ModelPart& >())
    ;

    py::class_<EmbeddedSkinVisualizationProcess, EmbeddedSkinVisualizationProcess::Pointer, Process>
    (m,"EmbeddedSkinVisualizationProcess")
    .def(py::init <
        ModelPart&,
        ModelPart&,
        const std::vector<const Variable <double>* >,
        const std::vector<const Variable< array_1d<double, 3> >* >,
        const std::vector<const Variable <double>* >,
        const std::vector<const Variable< array_1d<double, 3> >* >,
        const EmbeddedSkinVisualizationProcess::LevelSetType&,
        const EmbeddedSkinVisualizationProcess::ShapeFunctionsType&,
        const bool >())
    .def(py::init< Model&, Parameters >())
    .def(py::init< ModelPart&, ModelPart&, Parameters >())
    ;

    py::class_<IntegrationPointStatisticsProcess, IntegrationPointStatisticsProcess::Pointer, Process>
    (m, "IntegrationPointStatisticsProcess")
    .def(py::init<Model&, Parameters::Pointer>())
    ;

    py::class_<MassConservationCheckProcess, MassConservationCheckProcess::Pointer, Process>
    (m,"MassConservationCheckProcess")
    .def(py::init < ModelPart&, const bool, const int, const bool, const std::string >())
    .def(py::init< ModelPart&, Parameters& >())
    .def("Initialize", &MassConservationCheckProcess::Initialize)
    .def("ExecuteInTimeStep", &MassConservationCheckProcess::ExecuteInTimeStep)
    .def("ComputePositiveVolume", &MassConservationCheckProcess::ComputePositiveVolume)
    .def("ComputeNegativeVolume", &MassConservationCheckProcess::ComputeNegativeVolume)
    .def("ComputeInterfaceArea", &MassConservationCheckProcess::ComputeInterfaceArea)
    .def("ComputeFlowOverBoundary", &MassConservationCheckProcess::ComputeFlowOverBoundary)
    ;

    py::class_<ShockCapturingPhysicsBasedProcess, ShockCapturingPhysicsBasedProcess::Pointer, Process>
    (m, "ShockCapturingPhysicsBasedProcess")
    .def(py::init < Model&, Parameters >())
    .def(py::init < ModelPart&, Parameters >())
    ;

    py::class_<ShockCapturingEntropyViscosityProcess, ShockCapturingEntropyViscosityProcess::Pointer, Process>
    (m, "ShockCapturingEntropyViscosityProcess")
    .def(py::init < Model&, Parameters >())
    .def(py::init < ModelPart&, Parameters >())
    ;

    py::class_<TwoFluidsInletProcess, TwoFluidsInletProcess::Pointer, Process>
    (m,"TwoFluidsInletProcess")
    .def(py::init< ModelPart&, Parameters&, Process::Pointer >())
    .def("SmoothDistanceField", &TwoFluidsInletProcess::SmoothDistanceField)
    ;

    py::class_<DistanceSmoothingProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType>, DistanceSmoothingProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType>::Pointer, Process>(m,"DistanceSmoothingProcess2D")
    .def(py::init< ModelPart&, LinearSolverType::Pointer >())
    .def(py::init< ModelPart&, Parameters >())
    .def(py::init< Model&, Parameters >())
    ;

    py::class_<DistanceSmoothingProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType>, DistanceSmoothingProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType>::Pointer, Process>(m,"DistanceSmoothingProcess3D")
    .def(py::init< ModelPart&, LinearSolverType::Pointer >())
    .def(py::init< ModelPart&, Parameters >())
    .def(py::init< Model&, Parameters >())
    ;

    py::class_<CalulateLevelsetConsistentNodalGradientProcess, CalulateLevelsetConsistentNodalGradientProcess::Pointer, Process>(m,"CalulateLevelsetConsistentNodalGradientProcess")
    .def(py::init< ModelPart& >())
    .def(py::init< ModelPart&, Parameters >())
    .def(py::init< Model&, Parameters >())
    ;

    py::class_<ApplyCompressibleNavierStokesBoundaryConditionsProcess, ApplyCompressibleNavierStokesBoundaryConditionsProcess::Pointer, Process>(m, "ApplyCompressibleNavierStokesBoundaryConditionsProcess")
    .def(py::init<Model&, Parameters>())
    ;

    py::class_<ComputePressureCoefficientProcess, ComputePressureCoefficientProcess::Pointer, Process>(m, "ComputePressureCoefficientProcess")
    .def(py::init<Model&, Parameters>())
    ;

    py::class_<ComputeYPlusProcess, ComputeYPlusProcess::Pointer, Process>(m, "ComputeYPlusProcess")
    .def(py::init<Model&, Parameters>())
    ;
    py::class_<NavierStokesVectorialFractionalConvectionProcess<2, SparseSpaceType, LocalSpaceType, LinearSolverType>, NavierStokesVectorialFractionalConvectionProcess<2, SparseSpaceType, LocalSpaceType, LinearSolverType>::Pointer, Process>(m, "NavierStokesVectorialFractionalConvectionProcess2D")
        .def(py::init<Model &, LinearSolverType::Pointer, Parameters>())
        .def(py::init<ModelPart &, LinearSolverType::Pointer, Parameters>());

    py::class_<NavierStokesVectorialFractionalConvectionProcess<3, SparseSpaceType, LocalSpaceType, LinearSolverType>, NavierStokesVectorialFractionalConvectionProcess<3, SparseSpaceType, LocalSpaceType, LinearSolverType>::Pointer, Process>(m, "NavierStokesVectorialFractionalConvectionProcess3D")
        .def(py::init<Model &, LinearSolverType::Pointer, Parameters>())
        .def(py::init<ModelPart &, LinearSolverType::Pointer, Parameters>());
}

} // namespace Python.

} // Namespace Kratos
