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
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/material_point_search_utility.h"
#include "custom_utilities/material_point_generator_utility.cpp"
#include "custom_utilities/mpm_energy_calculation_utility.h"
#include "custom_utilities/mpm_volume_sum_utility.h"


namespace Kratos{
namespace Python{

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
        m.def("SearchElement", SearchElementAccordingToDimension, pybind11::arg("BackgroundGridModelPart"), pybind11::arg("MPMModelPart"), pybind11::arg("MaxNumberOfResults"), pybind11::arg("Tolerance"));
        m.def("GenerateMaterialPointElement", GenerateMaterialPointElementAccordingToDimension, pybind11::arg("BackgroundGridModelPart"),  pybind11::arg("InitialModelPart"), pybind11::arg("MPMModelPart"), pybind11::arg("IsMixedFormulation"));
        m.def("GenerateMaterialPointCondition", GenerateMaterialPointConditionAccordingToDimension, pybind11::arg("BackgroundGridModelPart"),  pybind11::arg("InitialModelPart"), pybind11::arg("MPMModelPart"));
        m.def("CalculateKineticEnergy", pybind11::overload_cast<ModelPart&>(&MPMEnergyCalculationUtility::CalculateKineticEnergy), pybind11::arg("MPMModelPart"));
        m.def("CalculateStrainEnergy", pybind11::overload_cast<ModelPart&>(&MPMEnergyCalculationUtility::CalculateStrainEnergy), pybind11::arg("MPMModelPart"));
        m.def("CalculatePotentialEnergy", pybind11::overload_cast<ModelPart&>(&MPMEnergyCalculationUtility::CalculatePotentialEnergy), pybind11::arg("MPMModelPart"));
        m.def("CalculateTotalEnergy", pybind11::overload_cast<ModelPart&>(&MPMEnergyCalculationUtility::CalculateTotalEnergy), pybind11::arg("MPMModelPart"));
        m.def("CalculateTotalEnergy", pybind11::overload_cast<ModelPart&>(&MPMEnergyCalculationUtility::CalculateTotalEnergy), pybind11::arg("MPMModelPart"));
        m.def("CalculateTotalMPVolume", pybind11::overload_cast<const ModelPart&>(&MPMVolumeSumUtility::AddModelPartMPMVolumeIntoGrid), pybind11::arg("MPMModelPart"));
    }

}  // namespace Python.
} // Namespace Kratos

