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
#include "custom_utilities/brute_force_material_point_locator.h"
#include "custom_utilities/mpm_energy_calculation_utility.h"
#include "custom_utilities/mapping_utilities/mpm_base_particle_mapping_utility.hpp"
#include "custom_utilities/mapping_utilities/mpm_flip_particle_mapping_utility.hpp"
#include "custom_utilities/mapping_utilities/mpm_pic_particle_mapping_utility.hpp"
#include "custom_utilities/mapping_utilities/mpm_tpic_particle_mapping_utility.hpp"


namespace Kratos::Python{

    class MPMBaseParticleMappingUtilityTrampoline : public MPMBaseParticleMappingUtility
    {
    public:
        //Inherit the constructors

        void P2GMomentum(Element& rElement, Node& rNode, const double& rN_i) override
        {
            using ReturnType = void;
            using BaseType = MPMBaseParticleMappingUtility;
            PYBIND11_OVERRIDE_PURE(
                ReturnType,
                BaseType,
                P2GMomentum,
                rElement,
                rNode,
                rN_i);
        }
        void P2GInertia(Element& rElement, Node& rNode, const double& rN_i) override
        {
            using ReturnType = void;
            using BaseType = MPMBaseParticleMappingUtility;
            PYBIND11_OVERRIDE_PURE(
                ReturnType,
                BaseType,
                P2GInertia,
                rElement,
                rNode,
                rN_i);
        }
        void G2PVelocity(Element& rElement, const array_1d<double, 3>& rNewMPAcceleration) override
        {
            using ReturnType = void;
            using BaseType = MPMBaseParticleMappingUtility;
            PYBIND11_OVERRIDE_PURE(
                ReturnType,
                BaseType,
                G2PVelocity,
                rElement,
                rNewMPAcceleration);
        }
    }; // class ControllerTrampoline

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


        // // MPM Residual Based Newton Raphson Strategy Type
        // pybind11::class_< MPMFlipParticleMappingUtility,typename MPMResidualBasedNewtonRaphsonStrategyType::Pointer, BaseSolvingStrategyType >(m,"MPMResidualBasedNewtonRaphsonStrategy")
        //     .def(pybind11::init< ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, int, bool, bool, bool >() )
        //     .def(pybind11::init< ModelPart&, BaseSchemeType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, int, bool, bool, bool >() )
        //     ;
        pybind11::class_<MPMBaseParticleMappingUtility, MPMBaseParticleMappingUtility::Pointer, MPMBaseParticleMappingUtilityTrampoline>(m, "MPMBaseParticleMappingUtility")
            .def("ResetBackgroundGrid", &MPMBaseParticleMappingUtility::ResetBackgroundGrid)
            .def("RunP2GMapping", &MPMBaseParticleMappingUtility::RunP2GMapping)
            .def("RunG2PMapping", &MPMBaseParticleMappingUtility::RunG2PMapping)
            ;

        pybind11::class_<MPMFlipParticleMappingUtility, MPMFlipParticleMappingUtility::Pointer, MPMBaseParticleMappingUtility>(m, "MPMFlipParticleMappingUtility")
            .def(pybind11::init<ModelPart&, ModelPart&, const unsigned int>(), pybind11::arg("material_point_model_part"), pybind11::arg("grid_model_part"), pybind11::arg("echo_level"))
            .def("Initialize", &MPMFlipParticleMappingUtility::Initialize)
            .def("RunP2GMapping", &MPMFlipParticleMappingUtility::RunP2GMapping)
            .def("RunG2PMapping", &MPMFlipParticleMappingUtility::RunG2PMapping)
            ;

        pybind11::class_<MPMPicParticleMappingUtility, MPMPicParticleMappingUtility::Pointer, MPMBaseParticleMappingUtility>(m, "MPMPicParticleMappingUtility")
            .def(pybind11::init<ModelPart&, ModelPart&, const unsigned int>(), pybind11::arg("material_point_model_part"), pybind11::arg("grid_model_part"), pybind11::arg("echo_level"))
            .def("Initialize", &MPMPicParticleMappingUtility::Initialize)
            .def("RunP2GMapping", &MPMPicParticleMappingUtility::RunP2GMapping)
            .def("RunG2PMapping", &MPMPicParticleMappingUtility::RunG2PMapping)
            ;

        pybind11::class_<MPMTpicParticleMappingUtility, MPMTpicParticleMappingUtility::Pointer, MPMBaseParticleMappingUtility>(m, "MPMTpicParticleMappingUtility")
            .def(pybind11::init<ModelPart&, ModelPart&, const unsigned int>(), pybind11::arg("material_point_model_part"), pybind11::arg("grid_model_part"), pybind11::arg("echo_level"))
            .def("Initialize", &MPMTpicParticleMappingUtility::Initialize)
            .def("RunP2GMapping", &MPMTpicParticleMappingUtility::RunP2GMapping)
            .def("RunG2PMapping", &MPMTpicParticleMappingUtility::RunG2PMapping)
            ;

        // Brute force material point (element/condition) locator
        py::class_<BruteForceMaterialPointLocator> (m, "BruteForceMaterialPointLocator")
            .def(py::init<ModelPart& >())
            .def("FindElement", &BruteForceMaterialPointLocator::FindElement, py::arg("point"), py::arg("abs_tolerance"))
            .def("FindCondition", &BruteForceMaterialPointLocator::FindCondition, py::arg("point"), py::arg("tolerance"))
            ;

        // Calculate energy utility
        py::class_< MPMEnergyCalculationUtility> (m,"EnergyCalculationUtility")
            .def(py::init<>())
            .def_static("CalculatePotentialEnergy", py::overload_cast<Element&>(&MPMEnergyCalculationUtility::CalculatePotentialEnergy), py::arg("element"))
            .def_static("CalculatePotentialEnergy", py::overload_cast<ModelPart&>(&MPMEnergyCalculationUtility::CalculatePotentialEnergy), py::arg("mpm_model_part"))
            .def_static("CalculateStrainEnergy", py::overload_cast<Element&>(&MPMEnergyCalculationUtility::CalculateStrainEnergy), py::arg("element"))
            .def_static("CalculateStrainEnergy", py::overload_cast<ModelPart&>(&MPMEnergyCalculationUtility::CalculateStrainEnergy), py::arg("mpm_model_part"))
            .def_static("CalculateKineticEnergy", py::overload_cast<Element&>(&MPMEnergyCalculationUtility::CalculateKineticEnergy), py::arg("element"))
            .def_static("CalculateKineticEnergy", py::overload_cast<ModelPart&>(&MPMEnergyCalculationUtility::CalculateKineticEnergy), py::arg("mpm_model_part"))
            .def_static("CalculateAllEnergies", py::overload_cast<Element&>(&MPMEnergyCalculationUtility::CalculateAllEnergies), py::arg("element"))
            .def_static("CalculateAllEnergies", py::overload_cast<ModelPart&>(&MPMEnergyCalculationUtility::CalculateAllEnergies), py::arg("mpm_model_part"))
            ;
    }

}  // namespace Kratos::Python.
