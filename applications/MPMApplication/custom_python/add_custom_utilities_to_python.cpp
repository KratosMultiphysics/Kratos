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
#include "custom_utilities/mapping_utilities/mpm_base_particle_mapping_utility.hpp"
#include "custom_utilities/mapping_utilities/mpm_flip_particle_mapping_utility.hpp"


namespace Kratos{
namespace Python{

    class MPMBaseParticleMappingUtilityTrampoline : public MPMBaseParticleMappingUtility
    {
    public:
        //Inherit the constructors

        void P2GMomentum(Element& rElement, Node& rNode, const double& rN_i,  const ProcessInfo& rCurrentProcessInfo) override
        {
            using ReturnType = void;
            using BaseType = MPMBaseParticleMappingUtility;
            PYBIND11_OVERRIDE_PURE(
                ReturnType,
                BaseType,
                P2GMomentum,
                rElement,
                rNode,
                rN_i,
                rCurrentProcessInfo);
        }
        void P2GInertia(Element& rElement, Node& rNode, const double& rN_i,  const ProcessInfo& rCurrentProcessInfo) override
        {
            using ReturnType = void;
            using BaseType = MPMBaseParticleMappingUtility;
            PYBIND11_OVERRIDE_PURE(
                ReturnType,
                BaseType,
                P2GInertia,
                rElement,
                rNode,
                rN_i,
                rCurrentProcessInfo);
        }
        void G2PVelocity(Element& rElement, const array_1d<double, 3>& rNewMPAcceleration, const ProcessInfo& rCurrentProcessInfo) override
        {
            using ReturnType = void;
            using BaseType = MPMBaseParticleMappingUtility;
            PYBIND11_OVERRIDE_PURE(
                ReturnType,
                BaseType,
                G2PVelocity,
                rElement,
                rNewMPAcceleration,
                rCurrentProcessInfo);
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

    }

}  // namespace Python.
} // Namespace Kratos

