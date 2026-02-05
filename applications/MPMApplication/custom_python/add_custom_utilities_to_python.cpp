//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//
//


// System includes

// External includes

// Project includes
#include "custom_utilities/mpm_energy_calculation_utility.h"
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/material_point_search_utility.h"
#include "custom_utilities/material_point_generator_utility.cpp"
#include "custom_utilities/mpm_energy_calculation_utility.cpp"


namespace Kratos::Python{

    void SearchElementAccordingToDimension(
        ModelPart& rBackgroundGridModelPart,
        ModelPart& rMPMModelPart,
        const std::size_t MaxNumberOfResults,
        const double Tolerance)
    {
        const auto dimension = rBackgroundGridModelPart.GetProcessInfo()[DOMAIN_SIZE];
        if (dimension == 2) MPMSearchElementUtility::SearchElement<2>(rBackgroundGridModelPart, rMPMModelPart, MaxNumberOfResults, Tolerance);
        else if (dimension == 3) MPMSearchElementUtility::SearchElement<3>(rBackgroundGridModelPart, rMPMModelPart, MaxNumberOfResults, Tolerance);
    }

    void GenerateMaterialPointElementAccordingToDimension(
        ModelPart& rBackgroundGridModelPart,
        ModelPart& rInitialModelPart,
        ModelPart& rMPMModelPart,
        bool IsMixedFormulation)
    {
        const auto dimension = rBackgroundGridModelPart.GetProcessInfo()[DOMAIN_SIZE];
        if (dimension == 2) MaterialPointGeneratorUtility::GenerateMaterialPointElement<2>(
            rBackgroundGridModelPart, rInitialModelPart, rMPMModelPart, IsMixedFormulation);
        else if (dimension == 3) MaterialPointGeneratorUtility::GenerateMaterialPointElement<3>(
            rBackgroundGridModelPart, rInitialModelPart, rMPMModelPart, IsMixedFormulation);
    }

    void GenerateMaterialPointConditionAccordingToDimension(
        ModelPart& rBackgroundGridModelPart,
        ModelPart& rInitialModelPart,
        ModelPart& rMPMModelPart)
    {
        const auto dimension = rBackgroundGridModelPart.GetProcessInfo()[DOMAIN_SIZE];
        if (dimension == 2) MaterialPointGeneratorUtility::GenerateMaterialPointCondition<2>(
            rBackgroundGridModelPart, rInitialModelPart, rMPMModelPart);
        else if (dimension == 3) MaterialPointGeneratorUtility::GenerateMaterialPointCondition<3>(
            rBackgroundGridModelPart, rInitialModelPart, rMPMModelPart);
    }

    void  AddCustomUtilitiesToPython(pybind11::module& m)
    {

        namespace py = pybind11;

        m.def("SearchElement", SearchElementAccordingToDimension);
        m.def("GenerateMaterialPointElement", GenerateMaterialPointElementAccordingToDimension);
        m.def("GenerateMaterialPointCondition", GenerateMaterialPointConditionAccordingToDimension);
        m.def("GenerateLagrangeNodes", MaterialPointGeneratorUtility::GenerateLagrangeNodes);

        // Calculate energy utility
        py::class_< MPMEnergyCalculationUtility> (m,"EnergyCalculationUtility")
            .def(py::init<>())
            .def_static("CalculatePotentialEnergy", py::overload_cast<Element&>(&MPMEnergyCalculationUtility::CalculatePotentialEnergy), py::arg("element"))
            .def_static("CalculatePotentialEnergy", py::overload_cast<ModelPart&>(&MPMEnergyCalculationUtility::CalculatePotentialEnergy), py::arg("model_part"))
            .def_static("CalculateStrainEnergy", py::overload_cast<Element&>(&MPMEnergyCalculationUtility::CalculateStrainEnergy), py::arg("element"))
            .def_static("CalculateStrainEnergy", py::overload_cast<ModelPart&>(&MPMEnergyCalculationUtility::CalculateStrainEnergy), py::arg("model_part"))
            .def_static("CalculateKineticEnergy", py::overload_cast<Element&>(&MPMEnergyCalculationUtility::CalculateKineticEnergy), py::arg("element"))
            .def_static("CalculateKineticEnergy", py::overload_cast<ModelPart&>(&MPMEnergyCalculationUtility::CalculateKineticEnergy), py::arg("model_part"))
            .def_static("CalculateTotalEnergy", py::overload_cast<Element&>(&MPMEnergyCalculationUtility::CalculateTotalEnergy), py::arg("element"))
            .def_static("CalculateTotalEnergy", py::overload_cast<ModelPart&>(&MPMEnergyCalculationUtility::CalculateTotalEnergy), py::arg("model_part"))
            .def_static("CalculateAllEnergies", [](ModelPart& model_part) -> py::tuple {
                    double potential, kinetic, strain, total;
                    MPMEnergyCalculationUtility().CalculateAllEnergies(model_part, potential, kinetic, strain, total);
                    return py::make_tuple(potential, kinetic, strain, total);
                }, py::arg("model_part"))
            ;
    }

}  // namespace Kratos::Python.
