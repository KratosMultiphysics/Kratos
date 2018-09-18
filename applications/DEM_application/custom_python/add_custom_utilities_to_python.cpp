//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

// External includes
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
#include "custom_utilities/benchmark_utilities.h"
#include "custom_utilities/inlet.h"
#include "custom_utilities/force_based_inlet.h"
#include "custom_utilities/reorder_consecutive_from_given_ids_model_part_io.h"
#include "custom_utilities/AuxiliaryUtilities.h"
#include "custom_utilities/excavator_utility.h"
#include "custom_utilities/analytic_tools/particles_history_watcher.h"

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
                                              Node < 3 > ::Pointer reference_node,
                                              Properties::Pointer r_params,
                                              const double radius,
                                              const Element& r_reference_element) {

    return r_creator_destructor.CreateSphericParticle(r_modelpart, r_Elem_Id, reference_node, r_params, radius, r_reference_element);
}

Element::Pointer CreateSphericParticle3(ParticleCreatorDestructor& r_creator_destructor,
                                              ModelPart& r_modelpart,
                                              int r_Elem_Id,
                                              Node < 3 > ::Pointer reference_node,
                                              Properties::Pointer r_params,
                                              const double radius,
                                              const std::string& element_name) {

    return r_creator_destructor.CreateSphericParticle(r_modelpart, r_Elem_Id, reference_node, r_params, radius, element_name);
}

Element::Pointer CreateSphericParticle4(ParticleCreatorDestructor& r_creator_destructor,
                                              ModelPart& r_modelpart,
                                              Node < 3 > ::Pointer reference_node,
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

using namespace pybind11;

void AddCustomUtilitiesToPython(pybind11::module& m) {

    class_<ParticleCreatorDestructor, ParticleCreatorDestructor::Pointer>(m, "ParticleCreatorDestructor")
        .def(init<>())
        .def(init<AnalyticWatcher::Pointer>())
        .def("CalculateSurroundingBoundingBox", &ParticleCreatorDestructor::CalculateSurroundingBoundingBox)
        .def("MarkParticlesForErasingGivenBoundingBox", &ParticleCreatorDestructor::MarkParticlesForErasingGivenBoundingBox)
        .def("MarkParticlesForErasingGivenScalarVariableValue", &ParticleCreatorDestructor::MarkParticlesForErasingGivenScalarVariableValue)
        .def("MarkParticlesForErasingGivenVectorVariableModulus", &ParticleCreatorDestructor::MarkParticlesForErasingGivenVectorVariableModulus)
        .def("GetHighNode", &ParticleCreatorDestructor::GetHighNode)
        .def("GetLowNode", &ParticleCreatorDestructor::GetLowNode)
        .def("GetDiameter", &ParticleCreatorDestructor::GetDiameter)
        .def("SetHighNode", &ParticleCreatorDestructor::SetHighNode)
        .def("SetLowNode", &ParticleCreatorDestructor::SetLowNode)
        .def("SetMaxNodeId", &ParticleCreatorDestructor::SetMaxNodeId)
        .def("DestroyParticlesOutsideBoundingBox", &ParticleCreatorDestructor::DestroyParticlesOutsideBoundingBox)
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
        ;

    class_<DEM_Inlet, DEM_Inlet::Pointer>(m, "DEM_Inlet")
        .def(init<ModelPart&>())
        .def("CreateElementsFromInletMesh", &DEM_Inlet::CreateElementsFromInletMesh)
        .def("InitializeDEM_Inlet", &DEM_Inlet::InitializeDEM_Inlet
            ,arg("model_part")
            ,arg("creator_destructor")
            ,arg("using_strategy_for_continuum") = false
            )
        .def("GetNumberOfParticlesInjectedSoFar", &DEM_Inlet::CreateElementsFromInletMesh)
        ;

    class_<DEM_Force_Based_Inlet, DEM_Force_Based_Inlet::Pointer, DEM_Inlet>(m, "DEM_Force_Based_Inlet")
        .def(init<ModelPart&, array_1d<double, 3>>())
        ;

    class_<SphericElementGlobalPhysicsCalculator, SphericElementGlobalPhysicsCalculator::Pointer >(m, "SphericElementGlobalPhysicsCalculator")
        .def(init<ModelPart&>())
        .def("CalculateTotalVolume", &SphericElementGlobalPhysicsCalculator::CalculateTotalVolume)
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
        .def("CalculateGravitationalPotentialEnergy", &SphericElementGlobalPhysicsCalculator::CalculateGravitationalPotentialEnergy)
        .def("CalculateTotalMomentum", &SphericElementGlobalPhysicsCalculator::CalculateTotalMomentum)
        .def("CalulateTotalAngularMomentum", &SphericElementGlobalPhysicsCalculator::CalulateTotalAngularMomentum)
        .def("CalculateSumOfInternalForces", &SphericElementGlobalPhysicsCalculator::CalculateSumOfInternalForces)
        .def("SetGravity", &SphericElementGlobalPhysicsCalculator::SetGravity)
        ;

    void (DemSearchUtilities::*SearchNodeNeigboursDistancesMM)(ModelPart&,ModelPart&,const double&,const Variable<double>&) = &DemSearchUtilities::SearchNodeNeigboursDistances<Variable<double> >;
    void (DemSearchUtilities::*SearchNodeNeigboursDistancesML)(NodesArrayType&,ModelPart&,const double&,const Variable<double>&) = &DemSearchUtilities::SearchNodeNeigboursDistances<Variable<double> >;
    void (DemSearchUtilities::*SearchNodeNeigboursDistancesLM)(ModelPart&,NodesArrayType&,const double&,const Variable<double>&) = &DemSearchUtilities::SearchNodeNeigboursDistances<Variable<double> >;
    void (DemSearchUtilities::*SearchNodeNeigboursDistancesLL)(NodesArrayType&,NodesArrayType&,const double&,const Variable<double>&) = &DemSearchUtilities::SearchNodeNeigboursDistances<Variable<double> >;

    class_<DemSearchUtilities, DemSearchUtilities::Pointer>(m, "DemSearchUtilities")
        .def(init<SpatialSearch::Pointer>())
        .def("SearchNodeNeighboursDistances", SearchNodeNeigboursDistancesMM)
        .def("SearchNodeNeighboursDistances", SearchNodeNeigboursDistancesML)
        .def("SearchNodeNeighboursDistances", SearchNodeNeigboursDistancesLM)
        .def("SearchNodeNeighboursDistances", SearchNodeNeigboursDistancesLL)
        ;

    class_<AnalyticModelPartFiller, AnalyticModelPartFiller::Pointer>(m, "AnalyticModelPartFiller")
        .def(init<>())
        .def("FillAnalyticModelPartGivenFractionOfParticlesToTransform", &AnalyticModelPartFiller::FillAnalyticModelPartGivenFractionOfParticlesToTransform
            ,arg("fraction_of_particles_to_convert")
            ,arg("spheres_model_part")
            ,arg("particle_creator_destructor")
            ,arg("analytic_sub_model_part_name") = ""
            )
        ;

    class_<AnalyticParticleWatcher, AnalyticParticleWatcher::Pointer>(m, "AnalyticParticleWatcher")
        .def(init<>())
        .def("MakeMeasurements", &AnalyticParticleWatcher::MakeMeasurements)
        //.def("GetTimeStepsData", &AnalyticParticleWatcher::GetTimeStepsData)
        //.def("GetParticleData", &AnalyticParticleWatcher::GetParticleData)
        //.def("GetAllParticlesData", &AnalyticParticleWatcher::GetAllParticlesData)
        .def("SetNodalMaxImpactVelocities", &AnalyticParticleWatcher::SetNodalMaxImpactVelocities)
        .def("SetNodalMaxFaceImpactVelocities", &AnalyticParticleWatcher::SetNodalMaxFaceImpactVelocities)
        ;


    class_<std::list<int>>(m, "IntList")
        .def(init<>())
        //.def("clear", &std::list<int>::clear)
        //.def("pop_back", &std::list<int>::pop_back)
        //.def("__len__", [](const std::list<int> &v) { return v.size(); })
        //.def("__iter__", [](std::list<int> &v) {
        //return make_iterator(v.begin(), v.end());
        //}
        ;

    class_<std::list<double>>(m, "DoubleList")
        .def(init<>())
        //.def("clear", &std::list<int>::clear)
        //.def("pop_back", &std::list<int>::pop_back)
        //.def("__len__", [](const std::list<int> &v) { return v.size(); })
        //.def("__iter__", [](std::list<int> &v) {
        //return make_iterator(v.begin(), v.end());
        //}
        ;


    class_<AnalyticFaceWatcher, AnalyticFaceWatcher::Pointer>(m, "AnalyticFaceWatcher")
        .def(init < ModelPart& >())
        .def("ClearData", &AnalyticFaceWatcher::ClearData)
        .def("MakeMeasurements", &AnalyticFaceWatcher::MakeMeasurements)
        //.def("GetTimeStepsData", &AnalyticFaceWatcher::GetTimeStepsData)
        //.def("GetFaceData", &AnalyticFaceWatcher::GetFaceData)
        //.def("GetAllFacesData", &AnalyticFaceWatcher::GetAllFacesData)
        .def("GetTotalFlux", &AnalyticFaceWatcher::GetTotalFlux)
        ;

    class_<DEM_FEM_Search, DEM_FEM_Search::Pointer>(m, "DEM_FEM_Search")
        .def(init<>())
        .def("GetBBHighPoint", &DEM_FEM_Search::GetBBHighPoint)
        .def("GetBBLowPoint", &DEM_FEM_Search::GetBBLowPoint)
        ;

    class_<PreUtilities, PreUtilities::Pointer >(m, "PreUtilities")
        .def(init<>())
        .def(init<ModelPart&>())
        .def("MeasureTopHeigh", Aux_MeasureTopHeight)
        .def("MeasureBotHeigh", Aux_MeasureBotHeight)
        .def("SetClusterInformationInProperties", &PreUtilities::SetClusterInformationInProperties)
        .def("CreateCartesianSpecimenMdpa", &PreUtilities::CreateCartesianSpecimenMdpa)
        .def("BreakBondUtility", &PreUtilities::BreakBondUtility)
        .def("FillAnalyticSubModelPartUtility", &PreUtilities::FillAnalyticSubModelPartUtility)
        ;

    class_<PostUtilities, PostUtilities::Pointer>(m, "PostUtilities")
        .def(init<>())
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

    class_<DEMFEMUtilities, DEMFEMUtilities::Pointer>(m, "DEMFEMUtilities")
        .def(init<>())
        .def("MoveAllMeshes", &DEMFEMUtilities::MoveAllMeshes)
        .def("MoveAllMeshesUsingATable", &DEMFEMUtilities::MoveAllMeshesUsingATable)
        .def("CreateRigidFacesFromAllElements", &DEMFEMUtilities::CreateRigidFacesFromAllElements)
        ;

    class_<BenchmarkUtils, BenchmarkUtils::Pointer>(m, "BenchmarkUtils")
        .def(init<>())
        .def("ComputeHydrodynamicForces", &BenchmarkUtils::ComputeHydrodynamicForces)
        ;

    class_<ReorderConsecutiveFromGivenIdsModelPartIO, ReorderConsecutiveFromGivenIdsModelPartIO::Pointer, ReorderConsecutiveModelPartIO>(m, "ReorderConsecutiveFromGivenIdsModelPartIO")
        .def(init<std::string const& >())
        .def(init<std::string const&, const int, const int, const int>())
        ;

    class_<AuxiliaryUtilities, AuxiliaryUtilities::Pointer>(m, "AuxiliaryUtilities")
        .def(init<>())
        .def("GetIthSubModelPartIsForceIntegrationGroup", &AuxiliaryUtilities::GetIthSubModelPartIsForceIntegrationGroup)
        .def("GetIthSubModelPartName", &AuxiliaryUtilities::GetIthSubModelPartName)
        .def("GetIthSubModelPartIdentifier", &AuxiliaryUtilities::GetIthSubModelPartIdentifier)
        .def("GetIthSubModelPartData", &AuxiliaryUtilities::GetIthSubModelPartData<double>)
        .def("GetIthSubModelPartData", &AuxiliaryUtilities::GetIthSubModelPartData<int>)
        .def("GetIthSubModelPartData", &AuxiliaryUtilities::GetIthSubModelPartData<array_1d<double,3> >)
        .def("GetIthSubModelPartData", &AuxiliaryUtilities::GetIthSubModelPartData<std::string>)
        .def("GetIthSubModelPartNodes", &AuxiliaryUtilities::GetIthSubModelPartNodes)
        ;

    class_<PropertiesProxiesManager, PropertiesProxiesManager::Pointer>(m, "PropertiesProxiesManager")
        .def(init<>())
        .def("CreatePropertiesProxies", CreatePropertiesProxies1)
        .def("CreatePropertiesProxies", CreatePropertiesProxies2)
        ;

    class_<ExcavatorUtility, ExcavatorUtility::Pointer >(m, "ExcavatorUtility")
        .def(init<ModelPart&, const double, const double, const double, const double, const double, const double, const double, const double, const double, const double, const double, const double, const double>())
        .def("ExecuteInitializeSolutionStep", &ExcavatorUtility::ExecuteInitializeSolutionStep)
        ;

    class_<AnalyticWatcher, AnalyticWatcher::Pointer>(m, "AnalyticWatcher")
        .def(init<>())
        ;

    class_<ParticlesHistoryWatcher, ParticlesHistoryWatcher::Pointer, AnalyticWatcher>(m, "ParticlesHistoryWatcher")
        .def(init<>())
        .def("GetNewParticlesData", &ParticlesHistoryWatcher::GetNewParticlesData)
        ;
    }



/*ModelPart::NodesContainerType::Pointer ModelPartGetNodes1(ModelPart& rModelPart)
{
	return rModelPart.pNodes();
}*/

}  // namespace Python

} // Namespace Kratos
