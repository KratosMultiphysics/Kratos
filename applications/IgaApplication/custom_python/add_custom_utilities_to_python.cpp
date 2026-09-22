/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
*/

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "iga_application_variables.h"

#include "spaces/ublas_space.h"
#include "custom_utilities/director_utilities.h"
#include "custom_utilities/iga_flags.h"
#include "custom_utilities/compute_interface_traction_shell_3p.h"
#include "custom_utilities/iga_sbm_mls_extension_operator_utility.h"
#include "custom_utilities/iga_sbm_taylor_extension_operator_utility.h"
#include "custom_utilities/iga_sbm_domain_classification_utility.h"


namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(
    pybind11::module& m)
{
    pybind11::class_< DirectorUtilities >(m, "DirectorUtilities")
        .def(pybind11::init<ModelPart&, Parameters>())
        .def("ComputeDirectors",
            &DirectorUtilities::ComputeDirectors)
        ;

    pybind11::class_< IgaFlags > iga_flags = pybind11::class_< IgaFlags >(m, "IgaFlags")
        .def(pybind11::init<>())
        ;

    iga_flags.attr("FIX_DISPLACEMENT_X") = IgaFlags::FIX_DISPLACEMENT_X;
    iga_flags.attr("FIX_DISPLACEMENT_Y") = IgaFlags::FIX_DISPLACEMENT_Y;
    iga_flags.attr("FIX_DISPLACEMENT_Z") = IgaFlags::FIX_DISPLACEMENT_Z;
    iga_flags.attr("FIX_ROTATION_X") = IgaFlags::FIX_ROTATION_X;
    iga_flags.attr("FIX_ROTATION_Y") = IgaFlags::FIX_ROTATION_Y;
    iga_flags.attr("FIX_ROTATION_Z") = IgaFlags::FIX_ROTATION_Z;

    pybind11::class_< ComputeInterfaceTractionShell3pUtility >(m, "ComputeInterfaceTractionShell3pUtility")
        .def_static(
            "ComputeAndSetInterfaceTraction",
            &ComputeInterfaceTractionShell3pUtility::ComputeAndSetInterfaceTraction)
        ;

    pybind11::class_< IgaSbmMlsExtensionOperatorUtility::ExtensionOperatorResult >(m, "IgaSbmMlsExtensionOperatorResult")
        .def(pybind11::init<>())
        .def_readonly("NodeIds", &IgaSbmMlsExtensionOperatorUtility::ExtensionOperatorResult::NodeIds)
        .def_readonly("Weights", &IgaSbmMlsExtensionOperatorUtility::ExtensionOperatorResult::Weights)
        .def_readonly("GradientWeights", &IgaSbmMlsExtensionOperatorUtility::ExtensionOperatorResult::GradientWeights)
        ;

    pybind11::class_< IgaSbmMlsExtensionOperatorUtility >(m, "IgaSbmMlsExtensionOperatorUtility")
        .def_static(
            "FindNClosestPoints",
            &IgaSbmMlsExtensionOperatorUtility::FindNClosestPoints)
        .def_static(
            "ComputeExtensionOperator",
            &IgaSbmMlsExtensionOperatorUtility::ComputeExtensionOperator)
        .def_static(
            "ClassifyActiveElements",
            &IgaSbmMlsExtensionOperatorUtility::ClassifyActiveElements)
        .def_static(
            "ComputeAssembledExtensionOperator",
            &IgaSbmMlsExtensionOperatorUtility::ComputeAssembledExtensionOperator)
        .def_static(
            "PrecomputeAndStoreCouplingExtensionOperators",
            &IgaSbmMlsExtensionOperatorUtility::PrecomputeAndStoreCouplingExtensionOperators)
        .def_static(
            "PrecomputeAndStoreCouplingShiftSources",
            &IgaSbmMlsExtensionOperatorUtility::PrecomputeAndStoreCouplingShiftSources)
        .def_static(
            "GetStoredDofNodeIds",
            &IgaSbmMlsExtensionOperatorUtility::GetStoredDofNodeIds)
        .def_static(
            "GetDifferentialArea",
            &IgaSbmMlsExtensionOperatorUtility::GetDifferentialArea)
        ;

    pybind11::class_< IgaSbmTaylorExtensionOperatorUtility >(m, "IgaSbmTaylorExtensionOperatorUtility")
        .def_static(
            "PrecomputeAndStoreDualCouplingTaylorData",
            &IgaSbmTaylorExtensionOperatorUtility::PrecomputeAndStoreDualCouplingTaylorData,
            pybind11::arg("rInterfaceConditions"),
            pybind11::arg("rMasterSurrogateConditions"),
            pybind11::arg("rSlaveSurrogateConditions"))
        .def_static(
            "ReplaceWithInterfaceConditions",
            &IgaSbmTaylorExtensionOperatorUtility::ReplaceWithInterfaceConditions,
            pybind11::arg("rCouplingModelPart"),
            pybind11::arg("StartId"))
        ;

    pybind11::class_< IgaSbmDomainClassificationUtility::ClassificationResult >(m, "IgaSbmDomainClassificationResult")
        .def(pybind11::init<>())
        .def_readonly("SpansU", &IgaSbmDomainClassificationUtility::ClassificationResult::SpansU)
        .def_readonly("SpansV", &IgaSbmDomainClassificationUtility::ClassificationResult::SpansV)
        .def_readonly("Classification", &IgaSbmDomainClassificationUtility::ClassificationResult::Classification)
        ;

    pybind11::class_< IgaSbmDomainClassificationUtility::SurrogateBoundarySegment >(m, "IgaSbmSurrogateBoundarySegment")
        .def(pybind11::init<>())
        .def_readonly("Start", &IgaSbmDomainClassificationUtility::SurrogateBoundarySegment::Start)
        .def_readonly("End", &IgaSbmDomainClassificationUtility::SurrogateBoundarySegment::End)
        .def_readonly("ActiveSpanIndexU", &IgaSbmDomainClassificationUtility::SurrogateBoundarySegment::ActiveSpanIndexU)
        .def_readonly("ActiveSpanIndexV", &IgaSbmDomainClassificationUtility::SurrogateBoundarySegment::ActiveSpanIndexV)
        ;

    pybind11::class_< IgaSbmDomainClassificationUtility::SurrogateBoundaryQuadraturePointInfo >(m, "IgaSbmSurrogateBoundaryQuadraturePointInfo")
        .def(pybind11::init<>())
        .def_readonly("ParametricPosition", &IgaSbmDomainClassificationUtility::SurrogateBoundaryQuadraturePointInfo::ParametricPosition)
        .def_readonly("PhysicalPosition", &IgaSbmDomainClassificationUtility::SurrogateBoundaryQuadraturePointInfo::PhysicalPosition)
        .def_readonly("NodeIds", &IgaSbmDomainClassificationUtility::SurrogateBoundaryQuadraturePointInfo::NodeIds)
        .def_readonly("ShapeFunctionValues", &IgaSbmDomainClassificationUtility::SurrogateBoundaryQuadraturePointInfo::ShapeFunctionValues)
        ;

    pybind11::class_< IgaSbmDomainClassificationUtility >(m, "IgaSbmDomainClassificationUtility")
        .def_static(
            "ClassifyKnotSpansFromGeometry",
            &IgaSbmDomainClassificationUtility::ClassifyKnotSpansFromGeometry)
        .def_static(
            "ComputeSurrogateBoundaryFromClassification",
            &IgaSbmDomainClassificationUtility::ComputeSurrogateBoundaryFromClassification)
        .def_static(
            "ValidateNonEmptyActiveDomainFromClassification",
            &IgaSbmDomainClassificationUtility::ValidateNonEmptyActiveDomainFromClassification)
        .def_static(
            "GetElementParametricPosition",
            &IgaSbmDomainClassificationUtility::GetElementParametricPosition)
        .def_static(
            "CreateSurrogateBoundaryQuadraturePoints",
            &IgaSbmDomainClassificationUtility::CreateSurrogateBoundaryQuadraturePoints)
        .def_static(
            "CreateSurrogateBoundaryConditions",
            &IgaSbmDomainClassificationUtility::CreateSurrogateBoundaryConditions)
        ;
}

} // namespace Python
} // Namespace Kratos
