// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "spaces/ublas_space.h"

//Utilities
#include "custom_utilities/rayleigh_damping_coefficients_utilities.h"
#include "custom_utilities/explicit_integration_utilities.h"
#include "custom_utilities/project_vector_on_surface_utility.h"
#include "custom_utilities/perturb_geometry_sparse_utility.h"
#include "custom_utilities/perturb_geometry_subgrid_utility.h"
#include "custom_utilities/section_properties_utility.h"


namespace Kratos::Python {

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    // Base types
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverTypeSparse;
    typedef LinearSolverTypeSparse::Pointer LinearSolverPointerTypeSparse;

    typedef LinearSolver<LocalSpaceType, LocalSpaceType > LinearSolverTypeDense;
    typedef LinearSolverTypeDense::Pointer LinearSolverPointerTypeDense;

    // RayleighDampingCoefficientsUtilities
    m.def("ComputeDampingCoefficients",&RayleighDampingCoefficientsUtilities::ComputeDampingCoefficients);

    // ExplicitIntegrationUtilities
    m.def("CalculateDeltaTime",&ExplicitIntegrationUtilities::CalculateDeltaTime);

    py::class_<ProjectVectorOnSurfaceUtility>(m,"ProjectVectorOnSurfaceUtility")
        .def_static("Execute",&ProjectVectorOnSurfaceUtility::Execute);

    py::class_<PerturbGeometrySparseUtility, PerturbGeometrySparseUtility::Pointer>(m,"PerturbGeometrySparseUtility")
        .def(py::init<ModelPart&,LinearSolverPointerTypeSparse, Parameters>())
        .def("CreateRandomFieldVectors", &PerturbGeometrySparseUtility::CreateRandomFieldVectors)
        .def("ApplyRandomFieldVectorsToGeometry", &PerturbGeometrySparseUtility::ApplyRandomFieldVectorsToGeometry)
        ;

    py::class_<PerturbGeometrySubgridUtility, PerturbGeometrySubgridUtility::Pointer>(m,"PerturbGeometrySubgridUtility")
        .def(py::init<ModelPart&, LinearSolverPointerTypeDense, Parameters>())
        .def("CreateRandomFieldVectors", &PerturbGeometrySubgridUtility::CreateRandomFieldVectors)
        .def("ApplyRandomFieldVectorsToGeometry", &PerturbGeometrySubgridUtility::ApplyRandomFieldVectorsToGeometry)
        ;

    py::class_<SectionPropertiesUtility>(m, "SectionPropertiesUtility")
    .def_static("InterpretSections",
        &SectionPropertiesUtility::InterpretSections);
}

}  // namespace Kratos::Python.

