//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

// External includes 
#include <boost/python.hpp>

// Project includes

#include "includes/model_part.h"
#include "custom_python/add_custom_utilities_to_python.h"
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

#include "boost/python/list.hpp"
#include "boost/python/extract.hpp"


namespace Kratos{

namespace Python{

typedef ModelPart::NodesContainerType::iterator      PointIterator;
typedef std::vector<array_1d<double, 3 > >           ComponentVectorType;
typedef std::vector<array_1d<double, 3 > >::iterator ComponentIteratorType;
typedef SpatialSearch::NodesContainerType            NodesArrayType;

boost::python::list Aux_MeasureTopHeight(PreUtilities& ThisPreUtils, ModelPart& rModelPart)
{
  double subtotal = 0.0;
  double weight = 0.0;
  ThisPreUtils.MeasureTopHeight(rModelPart,subtotal,weight);
  
  // Copy output to a Python list
  boost::python::list Out;
  boost::python::object py_subtotal(subtotal);
  boost::python::object py_weight(weight);
  Out.append( py_subtotal );
  Out.append( py_weight );
  return Out;
}

boost::python::list Aux_MeasureBotHeight(PreUtilities& ThisPreUtils, ModelPart& rModelPart)
{
  double subtotal = 0.0;
  double weight = 0.0;
  ThisPreUtils.MeasureBotHeight(rModelPart,subtotal,weight);
  
  // Copy output to a Python list
  boost::python::list Out;
  boost::python::object py_subtotal(subtotal);
  boost::python::object py_weight(weight);
  Out.append( py_subtotal );
  Out.append( py_weight );
  return Out;
}


void  AddCustomUtilitiesToPython() {
    
    using namespace boost::python;
    
      class_<ParticleCreatorDestructor, boost::noncopyable >
        ("ParticleCreatorDestructor", init<>())
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
        ;
      
      class_<DEM_Inlet, boost::noncopyable >
        ("DEM_Inlet", init<ModelPart&>())
        .def("CreateElementsFromInletMesh", &DEM_Inlet::CreateElementsFromInletMesh)        
        .def("InitializeDEM_Inlet", &DEM_Inlet::InitializeDEM_Inlet) 
        ;

     class_<SphericElementGlobalPhysicsCalculator, boost::noncopyable >
        ("SphericElementGlobalPhysicsCalculator", init<ModelPart&>())
        .def("CalculateTotalVolume", &SphericElementGlobalPhysicsCalculator::CalculateTotalVolume)
        .def("CalculateTotalMass", &SphericElementGlobalPhysicsCalculator::CalculateTotalMass)
        .def("CalculateMaxNodalVariable", &SphericElementGlobalPhysicsCalculator::CalculateMaxNodalVariable)
        .def("CalculateMinNodalVariable", &SphericElementGlobalPhysicsCalculator::CalculateMinNodalVariable)
        .def("CalculateD50", &SphericElementGlobalPhysicsCalculator::CalculateD50)
        .def("CalculateCenterOfMass", &SphericElementGlobalPhysicsCalculator::CalculateCenterOfMass)
        .def("GetInitialCenterOfMass", &SphericElementGlobalPhysicsCalculator::GetInitialCenterOfMass)
        .def("CalculateKinematicEnergy", &SphericElementGlobalPhysicsCalculator::CalculateKinematicEnergy)
        .def("CalculateGravitationalPotentialEnergy", &SphericElementGlobalPhysicsCalculator::CalculateGravitationalPotentialEnergy)
        .def("CalculateTotalMomentum", &SphericElementGlobalPhysicsCalculator::CalculateTotalMomentum)
        .def("CalulateTotalAngularMomentum", &SphericElementGlobalPhysicsCalculator::CalulateTotalAngularMomentum)
        ;

     void (DemSearchUtilities::*SearchNodeNeigboursDistancesMM)(ModelPart&,ModelPart&,const double&,const Variable<double>&) = &DemSearchUtilities::SearchNodeNeigboursDistances<Variable<double> >;
     void (DemSearchUtilities::*SearchNodeNeigboursDistancesML)(NodesArrayType&,ModelPart&,const double&,const Variable<double>&) = &DemSearchUtilities::SearchNodeNeigboursDistances<Variable<double> >;
     void (DemSearchUtilities::*SearchNodeNeigboursDistancesLM)(ModelPart&,NodesArrayType&,const double&,const Variable<double>&) = &DemSearchUtilities::SearchNodeNeigboursDistances<Variable<double> >;
     void (DemSearchUtilities::*SearchNodeNeigboursDistancesLL)(NodesArrayType&,NodesArrayType&,const double&,const Variable<double>&) = &DemSearchUtilities::SearchNodeNeigboursDistances<Variable<double> >;

     class_<DemSearchUtilities, boost::noncopyable >
         ("DemSearchUtilities", init<SpatialSearch::Pointer>())
         .def("SearchNodeNeighboursDistances", SearchNodeNeigboursDistancesMM)
         .def("SearchNodeNeighboursDistances", SearchNodeNeigboursDistancesML)
         .def("SearchNodeNeighboursDistances", SearchNodeNeigboursDistancesLM)
         .def("SearchNodeNeighboursDistances", SearchNodeNeigboursDistancesLL)
         ;

     class_<DEM_FEM_Search, boost::noncopyable >
        ("DEM_FEM_Search", init<>())
        .def("GetBBHighPoint", &DEM_FEM_Search::GetBBHighPoint)
        .def("GetBBLowPoint", &DEM_FEM_Search::GetBBLowPoint)
        ;
        
      class_<PreUtilities, boost::noncopyable >
        ("PreUtilities", init<ModelPart&>())
        .def("MeasureTopHeigh", Aux_MeasureTopHeight)        
        .def("MeasureBotHeigh", Aux_MeasureBotHeight) 
        ;

     class_<PostUtilities, boost::noncopyable >
        ("PostUtilities", init<>())
        .def("VelocityTrap", &PostUtilities::VelocityTrap)
        .def("AddModelPartToModelPart", &PostUtilities::AddModelPartToModelPart)
        .def("QuasiStaticAdimensionalNumber", &PostUtilities::QuasiStaticAdimensionalNumber)
        .def("IntegrationOfForces", &PostUtilities::IntegrationOfForces)
        ;
     
     class_<DEMFEMUtilities, boost::noncopyable >
        ("DEMFEMUtilities", init<>())
        .def("MoveAllMeshes", &DEMFEMUtilities::MoveAllMeshes)
        .def("CreateRigidFacesFromAllElements", &DEMFEMUtilities::CreateRigidFacesFromAllElements)     
        ;


     class_<BenchmarkUtils, boost::noncopyable >
        ("BenchmarkUtils", init<>())
        .def("ComputeHydrodynamicForces", &BenchmarkUtils::ComputeHydrodynamicForces)
        ;
        

    }

}  // namespace Python.

} // Namespace Kratos
