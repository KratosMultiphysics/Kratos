//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

// External includes
//NOTE: those two includes should go first in the include list of this file
#include "includes/define_python.h"
#include "includes/define.h"

#include <pybind11/pybind11.h>

// Project includes

#include "includes/model_part.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/analytic_tools/analytic_model_part_filler.h"
#include "custom_utilities/analytic_tools/analytic_particle_watcher.h"
#include "custom_utilities/analytic_tools/analytic_face_watcher.h"
#include "custom_utilities/create_and_destroy.h"
#include "custom_utilities/calculate_global_physical_properties.h"
#include "custom_utilities/pre_utilities.h"
#include "custom_utilities/post_utilities.h"
#include "custom_utilities/search_utilities.h"
#include "custom_utilities/omp_dem_search.h"
#include "custom_utilities/dem_fem_search.h"
#include "custom_utilities/dem_fem_utilities.h"
#include "custom_utilities/inlet.h"
#include "custom_utilities/fast_filling_creator.h"
#include "custom_utilities/force_based_inlet.h"
#include "custom_utilities/reorder_consecutive_from_given_ids_model_part_io.h"
#include "custom_utilities/AuxiliaryUtilities.h"
#include "custom_utilities/excavator_utility.h"
#include "custom_utilities/analytic_tools/particles_history_watcher.h"
#include "custom_utilities/move_mesh_utility.h"
#include "custom_utilities/stationarity_checker.h"
#include "custom_utilities/multiaxial_control_module_generalized_2d_utilities.hpp"
#include "custom_utilities/random_variable.h"
#include "custom_utilities/piecewise_linear_random_variable.h"
#include "custom_utilities/discrete_random_variable.h"
#include "custom_utilities/parallel_bond_utilities.h"


namespace Kratos {

namespace Python {

typedef ModelPart::NodesContainerType::iterator      PointIterator;
typedef std::vector<array_1d<double, 3 > >           ComponentVectorType;
typedef std::vector<array_1d<double, 3 > >::iterator ComponentIteratorType;
typedef SpatialSearch::NodesContainerType            NodesArrayType;

pybind11::list Aux_MeasureTopHeight(PreUtilities& ThisPreUtils, ModelPart& rModelPart)
{
    double subtotal = 0.0;
    double weight = 0.0;
    ThisPreUtils.MeasureTopHeight(rModelPart,subtotal,weight);

    // Copy output to a Python list
    pybind11::list Out;

    Out.append( subtotal );
    Out.append( weight );
    return Out;
}

pybind11::list Aux_MeasureBotHeight(PreUtilities& ThisPreUtils, ModelPart& rModelPart)
{
    double subtotal = 0.0;
    double weight = 0.0;
    ThisPreUtils.MeasureBotHeight(rModelPart,subtotal,weight);

    // Copy output to a Python list
    pybind11::list Out;

    Out.append( subtotal );
    Out.append( weight );
    return Out;
}

Element::Pointer CreateSphericParticle1(ParticleCreatorDestructor& r_creator_destructor,
                                                ModelPart& r_modelpart,
                                                int r_Elem_Id,
                                                const array_1d<double, 3 >& coordinates,
                                                Properties::Pointer r_params,
                                                const double radius,
                                                const Element& r_reference_element) {

    return r_creator_destructor.CreateSphericParticle(r_modelpart, r_Elem_Id, coordinates, r_params, radius, r_reference_element);
}

Element::Pointer CreateSphericParticle2(ParticleCreatorDestructor& r_creator_destructor,
                                              ModelPart& r_modelpart,
                                              int r_Elem_Id,
                                              Node ::Pointer reference_node,
                                              Properties::Pointer r_params,
                                              const double radius,
                                              const Element& r_reference_element) {

    return r_creator_destructor.CreateSphericParticle(r_modelpart, r_Elem_Id, reference_node, r_params, radius, r_reference_element);
}

Element::Pointer CreateSphericParticle3(ParticleCreatorDestructor& r_creator_destructor,
                                              ModelPart& r_modelpart,
                                              int r_Elem_Id,
                                              Node ::Pointer reference_node,
                                              Properties::Pointer r_params,
                                              const double radius,
                                              const std::string& element_name) {

    return r_creator_destructor.CreateSphericParticle(r_modelpart, r_Elem_Id, reference_node, r_params, radius, element_name);
}

Element::Pointer CreateSphericParticle4(ParticleCreatorDestructor& r_creator_destructor,
                                              ModelPart& r_modelpart,
                                              Node ::Pointer reference_node,
                                              Properties::Pointer r_params,
                                              const double radius,
                                              const std::string& element_name) {

    return r_creator_destructor.CreateSphericParticle(r_modelpart, reference_node, r_params, radius, element_name);
}

Element::Pointer CreateSphericParticle5(ParticleCreatorDestructor& r_creator_destructor,
                                              ModelPart& r_modelpart,
                                              int r_Elem_Id,
                                              const array_1d<double, 3 >& coordinates,
                                              Properties::Pointer r_params,
                                              const double radius,
                                              const std::string& element_name) {

    return r_creator_destructor.CreateSphericParticle(r_modelpart, r_Elem_Id, coordinates, r_params, radius, element_name);
}

Element::Pointer CreateSphericParticle6(ParticleCreatorDestructor& r_creator_destructor,
                                              ModelPart& r_modelpart,
                                              const array_1d<double, 3 >& coordinates,
                                              Properties::Pointer r_params,
                                              const double radius,
                                              const std::string& element_name) {

    return r_creator_destructor.CreateSphericParticle(r_modelpart, coordinates, r_params, radius, element_name);
}

void CreatePropertiesProxies1(PropertiesProxiesManager& r_properties_proxy_manager, ModelPart& r_modelpart) {
    r_properties_proxy_manager.CreatePropertiesProxies(r_modelpart);
}

void CreatePropertiesProxies2(PropertiesProxiesManager& r_properties_proxy_manager, ModelPart& r_modelpart, ModelPart& r_inlet_modelpart, ModelPart& r_clusters_modelpart) {
    r_properties_proxy_manager.CreatePropertiesProxies(r_modelpart, r_inlet_modelpart, r_clusters_modelpart);
}

namespace py = pybind11;

void AddCustomUtilitiesToPython(pybind11::module& m) {

    py::class_<ParticleCreatorDestructor, ParticleCreatorDestructor::Pointer>(m, "ParticleCreatorDestructor")
        .def(py::init<>())
        .def(py::init<Parameters>())
        .def(py::init<AnalyticWatcher::Pointer>())
        .def(py::init<AnalyticWatcher::Pointer, const Parameters&>())
        .def("CalculateSurroundingBoundingBox", &ParticleCreatorDestructor::CalculateSurroundingBoundingBox)
        .def("UpdateSurroundingBoundingBox", &ParticleCreatorDestructor::UpdateSurroundingBoundingBox)
        .def("MarkParticlesForErasingGivenBoundingBox", &ParticleCreatorDestructor::MarkParticlesForErasingGivenBoundingBox<SphericParticle>)
        .def("MarkParticlesForErasingGivenBoundingBox", &ParticleCreatorDestructor::MarkParticlesForErasingGivenBoundingBox<Cluster3D>)
        .def("MarkParticlesForErasingGivenScalarVariableValue", &ParticleCreatorDestructor::MarkParticlesForErasingGivenScalarVariableValue)
        .def("MarkParticlesForErasingGivenVectorVariableModulus", &ParticleCreatorDestructor::MarkParticlesForErasingGivenVectorVariableModulus)
        .def("MarkParticlesForErasingGivenCylinder", &ParticleCreatorDestructor::MarkParticlesForErasingGivenCylinder)
        .def("GetHighNode", &ParticleCreatorDestructor::GetHighNode)
        .def("GetLowNode", &ParticleCreatorDestructor::GetLowNode)
        .def("GetDiameter", &ParticleCreatorDestructor::GetDiameter)
        .def("SetHighNode", &ParticleCreatorDestructor::SetHighNode)
        .def("SetLowNode", &ParticleCreatorDestructor::SetLowNode)
        .def("SetMaxNodeId", &ParticleCreatorDestructor::SetMaxNodeId)
        .def("DestroyParticlesOutsideBoundingBox", &ParticleCreatorDestructor::DestroyParticlesOutsideBoundingBox<SphericParticle>)
        .def("DestroyParticlesOutsideBoundingBox", &ParticleCreatorDestructor::DestroyParticlesOutsideBoundingBox<Cluster3D>)
        .def("DestroyContactElementsOutsideBoundingBox", &ParticleCreatorDestructor::DestroyContactElementsOutsideBoundingBox)
        .def("FindMaxNodeIdInModelPart", &ParticleCreatorDestructor::FindMaxNodeIdInModelPart)
        .def("FindMaxElementIdInModelPart", &ParticleCreatorDestructor::FindMaxElementIdInModelPart)
        .def("FindMaxConditionIdInModelPart", &ParticleCreatorDestructor::FindMaxConditionIdInModelPart)
        .def("RenumberElementIdsFromGivenValue", &ParticleCreatorDestructor::RenumberElementIdsFromGivenValue)
        .def("CreateSphericParticle", CreateSphericParticle1)
        .def("CreateSphericParticle", CreateSphericParticle2)
        .def("CreateSphericParticle", CreateSphericParticle3)
        .def("CreateSphericParticle", CreateSphericParticle4)
        .def("CreateSphericParticle", CreateSphericParticle5)
        .def("CreateSphericParticle", CreateSphericParticle6)
        .def("DestroyMarkedParticles", &ParticleCreatorDestructor::DestroyMarkedParticles)
        .def("MarkContactElementsForErasing", &ParticleCreatorDestructor::MarkContactElementsForErasing)
        .def("MarkContactElementsForErasingContinuum", &ParticleCreatorDestructor::MarkContactElementsForErasingContinuum)
        .def("DestroyContactElements", &ParticleCreatorDestructor::DestroyContactElements)
        .def("MarkIsolatedParticlesForErasing", &ParticleCreatorDestructor::MarkIsolatedParticlesForErasing)
        ;

    py::class_<DEM_Inlet, DEM_Inlet::Pointer>(m, "DEM_Inlet")
        .def(py::init<ModelPart&>())
        .def(py::init<ModelPart&, const int>())
        .def(py::init<ModelPart&, Parameters&, const int>())
        .def("CreateElementsFromInletMesh", &DEM_Inlet::CreateElementsFromInletMesh)
        .def("InitializeDEM_Inlet", &DEM_Inlet::InitializeDEM_Inlet
            ,py::arg("model_part")
            ,py::arg("creator_destructor")
            ,py::arg("using_strategy_for_continuum") = false
            )
        .def("GetTotalNumberOfParticlesInjectedSoFar", &DEM_Inlet::GetTotalNumberOfParticlesInjectedSoFar)
        .def("GetTotalMassInjectedSoFar", &DEM_Inlet::GetTotalMassInjectedSoFar)
        .def("GetMaxRadius", &DEM_Inlet::GetMaxRadius)
        ;

    py::class_<Fast_Filling_Creator, Fast_Filling_Creator::Pointer>(m, "Fast_Filling_Creator")
        .def(py::init<const int>())
        .def(py::init<Parameters&, const int>())
        .def("GetRandomParticleRadius", &Fast_Filling_Creator::GetRandomParticleRadius, py::arg("creator_destructor"))
        .def("CheckHasIndentationOrNot", &Fast_Filling_Creator::CheckHasIndentationOrNot)
        ;

    py::class_<DEM_Force_Based_Inlet, DEM_Force_Based_Inlet::Pointer, DEM_Inlet>(m, "DEM_Force_Based_Inlet")
        .def(py::init<ModelPart&, array_1d<double, 3>, const int>())
        .def(py::init<ModelPart&, array_1d<double, 3>>())
        ;

    py::class_<SphericElementGlobalPhysicsCalculator, SphericElementGlobalPhysicsCalculator::Pointer >(m, "SphericElementGlobalPhysicsCalculator")
        .def(py::init<ModelPart&>())
        .def("CalculateTotalVolume", &SphericElementGlobalPhysicsCalculator::CalculateTotalVolume)
        .def("CalculatePorosityWithinSphere", &SphericElementGlobalPhysicsCalculator::CalculatePorosityWithinSphere)
        .def("CalculateTotalMass", &SphericElementGlobalPhysicsCalculator::CalculateTotalMass)
        .def("CalculateMaxNodalVariable", &SphericElementGlobalPhysicsCalculator::CalculateMaxNodalVariable)
        .def("CalculateMinNodalVariable", &SphericElementGlobalPhysicsCalculator::CalculateMinNodalVariable)
        .def("CalculateD50", &SphericElementGlobalPhysicsCalculator::CalculateD50)
        .def("CalculateCenterOfMass", &SphericElementGlobalPhysicsCalculator::CalculateCenterOfMass)
        .def("GetInitialCenterOfMass", &SphericElementGlobalPhysicsCalculator::GetInitialCenterOfMass)
        .def("CalculateTranslationalKinematicEnergy", &SphericElementGlobalPhysicsCalculator::CalculateTranslationalKinematicEnergy)
        .def("CalculateRotationalKinematicEnergy", &SphericElementGlobalPhysicsCalculator::CalculateRotationalKinematicEnergy)
        .def("CalculateElasticEnergy", &SphericElementGlobalPhysicsCalculator::CalculateElasticEnergy)
        .def("CalculateInelasticFrictionalEnergy", &SphericElementGlobalPhysicsCalculator::CalculateInelasticFrictionalEnergy)
        .def("CalculateInelasticViscodampingEnergy", &SphericElementGlobalPhysicsCalculator::CalculateInelasticViscodampingEnergy)
        .def("CalculateInelasticRollingResistanceEnergy", &SphericElementGlobalPhysicsCalculator::CalculateInelasticRollingResistanceEnergy)
        .def("CalculateGravitationalPotentialEnergy", &SphericElementGlobalPhysicsCalculator::CalculateGravitationalPotentialEnergy)
        .def("CalculateTotalMomentum", &SphericElementGlobalPhysicsCalculator::CalculateTotalMomentum)
        .def("CalulateTotalAngularMomentum", &SphericElementGlobalPhysicsCalculator::CalulateTotalAngularMomentum)
        .def("CalculateSumOfInternalForces", &SphericElementGlobalPhysicsCalculator::CalculateSumOfInternalForces)
        .def("CalculateParticleNumberTimesMaxNormalBallToBallForceTimesRadius", &SphericElementGlobalPhysicsCalculator::CalculateParticleNumberTimesMaxNormalBallToBallForceTimesRadius)
        ;
    
    py::class_<ContactElementGlobalPhysicsCalculator, ContactElementGlobalPhysicsCalculator::Pointer >(m, "ContactElementGlobalPhysicsCalculator")
        .def(py::init<>())
        .def("CalculateTotalStressTensor", &ContactElementGlobalPhysicsCalculator::CalculateTotalStressTensor)
        .def("CalculateTotalStressTensorWithinSphere", &ContactElementGlobalPhysicsCalculator::CalculateTotalStressTensorWithinSphere)
        .def("CalculateFabricTensorWithinSphere", &ContactElementGlobalPhysicsCalculator::CalculateFabricTensorWithinSphere)
        .def("CalculateAveragedCoordinationNumberWithinSphere", &ContactElementGlobalPhysicsCalculator::CalculateAveragedCoordinationNumberWithinSphere)
        .def("CalculateUnbalancedForceWithinSphere", &ContactElementGlobalPhysicsCalculator::CalculateUnbalancedForceWithinSphere)
        ;

    void (DemSearchUtilities::*SearchNodeNeigboursDistancesMM)(ModelPart&,ModelPart&,const double&,const Variable<double>&) = &DemSearchUtilities::SearchNodeNeigboursDistances<Variable<double> >;
    void (DemSearchUtilities::*SearchNodeNeigboursDistancesML)(NodesArrayType&,ModelPart&,const double&,const Variable<double>&) = &DemSearchUtilities::SearchNodeNeigboursDistances<Variable<double> >;
    void (DemSearchUtilities::*SearchNodeNeigboursDistancesLM)(ModelPart&,NodesArrayType&,const double&,const Variable<double>&) = &DemSearchUtilities::SearchNodeNeigboursDistances<Variable<double> >;
    void (DemSearchUtilities::*SearchNodeNeigboursDistancesLL)(NodesArrayType&,NodesArrayType&,const double&,const Variable<double>&) = &DemSearchUtilities::SearchNodeNeigboursDistances<Variable<double> >;

    py::class_<DemSearchUtilities, DemSearchUtilities::Pointer>(m, "DemSearchUtilities")
        .def(py::init<SpatialSearch::Pointer>())
        .def("SearchNodeNeighboursDistances", SearchNodeNeigboursDistancesMM)
        .def("SearchNodeNeighboursDistances", SearchNodeNeigboursDistancesML)
        .def("SearchNodeNeighboursDistances", SearchNodeNeigboursDistancesLM)
        .def("SearchNodeNeighboursDistances", SearchNodeNeigboursDistancesLL)
        ;

    py::class_<AnalyticModelPartFiller, AnalyticModelPartFiller::Pointer>(m, "AnalyticModelPartFiller")
        .def(py::init<>())
        .def("FillAnalyticModelPartGivenFractionOfParticlesToTransform", &AnalyticModelPartFiller::FillAnalyticModelPartGivenFractionOfParticlesToTransform
            ,py::arg("fraction_of_particles_to_convert")
            ,py::arg("spheres_model_part")
            ,py::arg("particle_creator_destructor")
            ,py::arg("analytic_sub_model_part_name") = ""
            )
        ;

    py::class_<AnalyticParticleWatcher, AnalyticParticleWatcher::Pointer>(m, "AnalyticParticleWatcher")
        .def(py::init<>())
        .def("MakeMeasurements", &AnalyticParticleWatcher::MakeMeasurements)
        //.def("GetTimeStepsData", &AnalyticParticleWatcher::GetTimeStepsData)
        //.def("GetParticleData", &AnalyticParticleWatcher::GetParticleData)
        //.def("GetAllParticlesData", &AnalyticParticleWatcher::GetAllParticlesData)
        .def("SetNodalMaxImpactVelocities", &AnalyticParticleWatcher::SetNodalMaxImpactVelocities)
        .def("SetNodalMaxFaceImpactVelocities", &AnalyticParticleWatcher::SetNodalMaxFaceImpactVelocities)
        ;

    py::class_<std::list<int>>(m, "IntList")
        .def(py::init<>())
        //.def("clear", &std::list<int>::clear)
        //.def("pop_back", &std::list<int>::pop_back)
        //.def("__len__", [](const std::list<int> &v) { return v.size(); })
        //.def("__iter__", [](std::list<int> &v) {
        //return make_iterator(v.begin(), v.end());
        //}
        ;

    py::class_<std::list<double>>(m, "DoubleList")
        .def(py::init<>())
        //.def("clear", &std::list<int>::clear)
        //.def("pop_back", &std::list<int>::pop_back)
        //.def("__len__", [](const std::list<int> &v) { return v.size(); })
        //.def("__iter__", [](std::list<int> &v) {
        //return make_iterator(v.begin(), v.end());
        //}
        ;

    py::class_<AnalyticFaceWatcher, AnalyticFaceWatcher::Pointer>(m, "AnalyticFaceWatcher")
        .def(py::init<ModelPart& >())
        .def("ClearData", &AnalyticFaceWatcher::ClearData)
        .def("MakeMeasurements", &AnalyticFaceWatcher::MakeMeasurements)
        //.def("GetTimeStepsData", &AnalyticFaceWatcher::GetTimeStepsData)
        //.def("GetFaceData", &AnalyticFaceWatcher::GetFaceData)
        //.def("GetAllFacesData", &AnalyticFaceWatcher::GetAllFacesData)
        .def("GetTotalFlux", &AnalyticFaceWatcher::GetTotalFlux)
        ;

    py::class_<DEM_FEM_Search, DEM_FEM_Search::Pointer>(m, "DEM_FEM_Search")
        .def(py::init<>())
        .def("GetBBHighPoint", &DEM_FEM_Search::GetBBHighPoint)
        .def("GetBBLowPoint", &DEM_FEM_Search::GetBBLowPoint)
        ;

    py::class_<PreUtilities, PreUtilities::Pointer >(m, "PreUtilities")
        .def(py::init<>())
        .def(py::init<ModelPart&>())
        .def("MeasureTopHeigh", Aux_MeasureTopHeight)
        .def("MeasureBotHeigh", Aux_MeasureBotHeight)
        .def("SetClusterInformationInProperties", &PreUtilities::SetClusterInformationInProperties)
        .def("CreateCartesianSpecimenMdpa", &PreUtilities::CreateCartesianSpecimenMdpa)
        .def("BreakBondUtility", &PreUtilities::BreakBondUtility)
        .def("FillAnalyticSubModelPartUtility", &PreUtilities::FillAnalyticSubModelPartUtility)
        .def("MarkToEraseParticlesOutsideRadius", &PreUtilities::MarkToEraseParticlesOutsideRadius)
        .def("MarkToEraseParticlesOutsideBoundary", &PreUtilities::MarkToEraseParticlesOutsideBoundary)
        .def("MarkToEraseParticlesOutsideRadiusForGettingCylinder", &PreUtilities::MarkToEraseParticlesOutsideRadiusForGettingCylinder)
        .def("ApplyConcentricForceOnParticles", &PreUtilities::ApplyConcentricForceOnParticles)
        .def("ResetSkinParticles", &PreUtilities::ResetSkinParticles)
        .def("SetSkinParticlesInnerCircularBoundary", &PreUtilities::SetSkinParticlesInnerCircularBoundary)
        .def("SetSkinParticlesOuterCircularBoundary", &PreUtilities::SetSkinParticlesOuterCircularBoundary)
        .def("SetSkinParticlesOuterSquaredBoundary", &PreUtilities::SetSkinParticlesOuterSquaredBoundary)
        .def("PrintNumberOfNeighboursHistogram", &PreUtilities::PrintNumberOfNeighboursHistogram)
        ;

    py::class_<PostUtilities, PostUtilities::Pointer>(m, "PostUtilities")
        .def(py::init<>())
        .def("VelocityTrap", &PostUtilities::VelocityTrap)
        .def("AddModelPartToModelPart", &PostUtilities::AddModelPartToModelPart)
        .def("AddSpheresNotBelongingToClustersToMixModelPart", &PostUtilities::AddSpheresNotBelongingToClustersToMixModelPart)
        .def("QuasiStaticAdimensionalNumber", &PostUtilities::QuasiStaticAdimensionalNumber)
        .def("IntegrationOfForces", &PostUtilities::IntegrationOfForces)
        .def("IntegrationOfElasticForces", &PostUtilities::IntegrationOfElasticForces)
        .def("ComputePoisson", &PostUtilities::ComputePoisson)
        .def("ComputePoisson2D", &PostUtilities::ComputePoisson2D)
        .def("ComputeEulerAngles", &PostUtilities::ComputeEulerAngles)
        ;

    py::class_<ParallelBondUtilities, ParallelBondUtilities::Pointer>(m, "ParallelBondUtilities")
        .def(py::init<>())
        .def("SetCurrentIndentationAsAReferenceInParallelBonds", &ParallelBondUtilities::SetCurrentIndentationAsAReferenceInParallelBonds)
        .def("SetCurrentIndentationAsAReferenceInParallelBondsForPBM", &ParallelBondUtilities::SetCurrentIndentationAsAReferenceInParallelBondsForPBM)
        ;

    py::class_<DEMFEMUtilities, DEMFEMUtilities::Pointer>(m, "DEMFEMUtilities")
        .def(py::init<>())
        .def("MoveAllMeshes", &DEMFEMUtilities::MoveAllMeshes)
        .def("CreateRigidFacesFromAllElements", &DEMFEMUtilities::CreateRigidFacesFromAllElements)
        ;

    py::class_<ReorderConsecutiveFromGivenIdsModelPartIO, ReorderConsecutiveFromGivenIdsModelPartIO::Pointer, ReorderConsecutiveModelPartIO>(m, "ReorderConsecutiveFromGivenIdsModelPartIO")
        .def(py::init<std::string const& >())
        .def(py::init<std::string const&, const int, const int, const int>())
        .def(py::init<std::string const&, const int, const int, const int, const Flags>())
        ;

    py::class_<AuxiliaryUtilities, AuxiliaryUtilities::Pointer>(m, "AuxiliaryUtilities")
        .def(py::init<>())
        .def("ComputeAverageZStressFor2D", &AuxiliaryUtilities::ComputeAverageZStressFor2D)
        .def("UpdateTimeInOneModelPart", &AuxiliaryUtilities::UpdateTimeInOneModelPart)
        ;

    py::class_<PropertiesProxiesManager, PropertiesProxiesManager::Pointer>(m, "PropertiesProxiesManager")
        .def(py::init<>())
        .def("CreatePropertiesProxies", CreatePropertiesProxies1)
        .def("CreatePropertiesProxies", CreatePropertiesProxies2)
        ;

    py::class_<ExcavatorUtility, ExcavatorUtility::Pointer >(m, "ExcavatorUtility")
        .def(py::init<ModelPart&, const double, const double, const double, const double, const double, const double, const double, const double, const double, const double, const double, const double, const double>())
        .def("ExecuteInitializeSolutionStep", &ExcavatorUtility::ExecuteInitializeSolutionStep)
        ;

    py::class_<AnalyticWatcher, AnalyticWatcher::Pointer>(m, "AnalyticWatcher")
        .def(py::init<>())
        ;

    py::class_<ParticlesHistoryWatcher, ParticlesHistoryWatcher::Pointer, AnalyticWatcher>(m, "ParticlesHistoryWatcher")
        .def(py::init<>())
        .def("GetNewParticlesData", &ParticlesHistoryWatcher::GetNewParticlesData)
        ;

    py::class_<MoveMeshUtility, MoveMeshUtility::Pointer>(m, "MoveMeshUtility")
        .def(py::init<>())
        .def("MoveDemMesh", &MoveMeshUtility::MoveDemMesh)
        ;

    py::class_<StationarityChecker, StationarityChecker::Pointer>(m, "StationarityChecker")
        .def(py::init<>())
        .def("CheckIfItsTimeToChangeGravity", &StationarityChecker::CheckIfItsTimeToChangeGravity)
        .def("CheckIfVariableIsNullInModelPart", &StationarityChecker::CheckIfVariableIsNullInModelPart)
        ;

    py::class_<MultiaxialControlModuleGeneralized2DUtilities, MultiaxialControlModuleGeneralized2DUtilities::Pointer>(m, "MultiaxialControlModuleGeneralized2DUtilities")
        .def(py::init<ModelPart&,ModelPart&,Parameters&>())
        .def("ExecuteInitialize", &MultiaxialControlModuleGeneralized2DUtilities::ExecuteInitialize)
        .def("ExecuteInitializeSolutionStep", &MultiaxialControlModuleGeneralized2DUtilities::ExecuteInitializeSolutionStep)
        .def("ExecuteFinalizeSolutionStep", &MultiaxialControlModuleGeneralized2DUtilities::ExecuteFinalizeSolutionStep)
        ;

    py::class_<RandomVariable, RandomVariable::Pointer>(m, "RandomVariable")
        .def(py::init<const Parameters>())
        .def("GetSupport", &RandomVariable::GetSupport)
        ;

    py::class_<PiecewiseLinearRandomVariable, PiecewiseLinearRandomVariable::Pointer, RandomVariable>(m, "PiecewiseLinearRandomVariable")
        .def(py::init<const Parameters>())
        .def(py::init<const Parameters, const int>())
        .def("Sample", &PiecewiseLinearRandomVariable::Sample)
        .def("ProbabilityDensity", &PiecewiseLinearRandomVariable::ProbabilityDensity)
        .def("GetMean", &PiecewiseLinearRandomVariable::GetMean)
        ;

    py::class_<DiscreteRandomVariable, DiscreteRandomVariable::Pointer, RandomVariable>(m, "DiscreteRandomVariable")
        .def(py::init<const Parameters>())
        .def(py::init<const Parameters, const int>())
        .def("Sample", &DiscreteRandomVariable::Sample)
        .def("ProbabilityDensity", &DiscreteRandomVariable::ProbabilityDensity)
        .def("GetMean", &DiscreteRandomVariable::GetMean)
        ;
    }

/*ModelPart::NodesContainerType::Pointer ModelPartGetNodes1(ModelPart& rModelPart)
{
	return rModelPart.pNodes();
}*/

}  // namespace Python

} // Namespace Kratos
