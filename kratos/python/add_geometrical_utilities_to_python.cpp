//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "processes/process.h"
#include "python/add_geometrical_utilities_to_python.h"

//Geometrical utilities
#include "utilities/normal_calculation_utils.h"
#include "utilities/body_normal_calculation_utils.h"
#include "utilities/body_distance_calculation_utils.h"
#include "utilities/signed_distance_calculation_utils.h"
#include "utilities/brute_force_point_locator.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/binbased_fast_point_locator_conditions.h"
#include "utilities/binbased_nodes_in_element_locator.h"
#include "utilities/embedded_skin_utility.h"
#include "utilities/geometry_tester.h"
#include "utilities/cutting_utility.h"
#include "utilities/geometrical_transformation_utilities.h"
#include "utilities/iso_printer.h"
#include "utilities/interval_utility.h"
#include "utilities/convect_particles_utilities.h"
#include "utilities/mls_shape_functions_utility.h"
#include "utilities/tessellation_utilities/delaunator_utilities.h"

namespace Kratos {
namespace Python {

// Embedded skin utility auxiliar functions
template<std::size_t TDim>
void InterpolateMeshVariableToSkinDouble(
    EmbeddedSkinUtility<TDim> &rEmbeddedSkinUtility,
    const Variable<double> &rVariable,
    const Variable<double> &rEmbeddedVariable)
{
    rEmbeddedSkinUtility.InterpolateMeshVariableToSkin(rVariable, rEmbeddedVariable);
}

template<std::size_t TDim>
void InterpolateMeshVariableToSkinArray(
    EmbeddedSkinUtility<TDim> &rEmbeddedSkinUtility,
    const Variable<array_1d<double,3>> &rVariable,
    const Variable<array_1d<double,3>> &rEmbeddedVariable)
{
    rEmbeddedSkinUtility.InterpolateMeshVariableToSkin(rVariable, rEmbeddedVariable);
}

template<std::size_t TDim>
void InterpolateDiscontinuousMeshVariableToSkinDouble(
    EmbeddedSkinUtility<TDim> &rEmbeddedSkinUtility,
    const Variable<double> &rVariable,
    const Variable<double> &rEmbeddedVariable,
    const std::string &rInterfaceSide)
{
    rEmbeddedSkinUtility.InterpolateDiscontinuousMeshVariableToSkin(rVariable, rEmbeddedVariable, rInterfaceSide);
}

template<std::size_t TDim>
void InterpolateDiscontinuousMeshVariableToSkinArray(
    EmbeddedSkinUtility<TDim> &rEmbeddedSkinUtility,
    const Variable<array_1d<double,3>> &rVariable,
    const Variable<array_1d<double,3>> &rEmbeddedVariable,
    const std::string &rInterfaceSide)
{
    rEmbeddedSkinUtility.InterpolateDiscontinuousMeshVariableToSkin(rVariable, rEmbeddedVariable, rInterfaceSide);
}

void AddGeometricalUtilitiesToPython(pybind11::module &m)
{
    namespace py = pybind11;

    py::class_<NormalCalculationUtils > (m,"NormalCalculationUtils")
        .def(py::init<>())
        .def("CalculateNormalsInConditions", &NormalCalculationUtils::CalculateNormalsInContainer<ModelPart::ConditionsContainerType>)
        .def("CalculateNormalsInElements", &NormalCalculationUtils::CalculateNormalsInContainer<ModelPart::ElementsContainerType>)
        .def("CalculateNormals", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart, const bool EnforceGenericAlgorithm, const NormalCalculationUtils::NormalVariableType& rNormalVariable){
            rNormalCalculationUtils.CalculateNormals<ModelPart::ConditionsContainerType,true>(rModelPart, EnforceGenericAlgorithm, rNormalVariable);})
        .def("CalculateNormals", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart, const bool EnforceGenericAlgorithm){
            rNormalCalculationUtils.CalculateNormals<ModelPart::ConditionsContainerType,true>(rModelPart, EnforceGenericAlgorithm);})
        .def("CalculateNormals", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart){
            rNormalCalculationUtils.CalculateNormals<ModelPart::ConditionsContainerType,true>(rModelPart);})
        .def("CalculateNormalsNonHistorical", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart, const bool EnforceGenericAlgorithm, const NormalCalculationUtils::NormalVariableType& rNormalVariable){
            rNormalCalculationUtils.CalculateNormals<ModelPart::ConditionsContainerType,false>(rModelPart, EnforceGenericAlgorithm, rNormalVariable);})
        .def("CalculateNormalsNonHistorical", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart, const bool EnforceGenericAlgorithm){
            rNormalCalculationUtils.CalculateNormals<ModelPart::ConditionsContainerType,false>(rModelPart, EnforceGenericAlgorithm);})
        .def("CalculateNormalsNonHistorical", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart){
            rNormalCalculationUtils.CalculateNormals<ModelPart::ConditionsContainerType,false>(rModelPart);})
        .def("CalculateUnitNormals", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart, const bool EnforceGenericAlgorithm, const NormalCalculationUtils::NormalVariableType& rNormalVariable){
            rNormalCalculationUtils.CalculateUnitNormals<ModelPart::ConditionsContainerType,true>(rModelPart, EnforceGenericAlgorithm, rNormalVariable);})
        .def("CalculateUnitNormals", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart, const bool EnforceGenericAlgorithm){
            rNormalCalculationUtils.CalculateUnitNormals<ModelPart::ConditionsContainerType,true>(rModelPart, EnforceGenericAlgorithm);})
        .def("CalculateUnitNormals", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart){
            rNormalCalculationUtils.CalculateUnitNormals<ModelPart::ConditionsContainerType,true>(rModelPart);})
        .def("CalculateUnitNormalsNonHistorical", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart, const bool EnforceGenericAlgorithm, const NormalCalculationUtils::NormalVariableType& rNormalVariable){
            rNormalCalculationUtils.CalculateUnitNormals<ModelPart::ConditionsContainerType,false>(rModelPart, EnforceGenericAlgorithm, rNormalVariable);})
        .def("CalculateUnitNormalsNonHistorical", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart, const bool EnforceGenericAlgorithm){
            rNormalCalculationUtils.CalculateUnitNormals<ModelPart::ConditionsContainerType,false>(rModelPart, EnforceGenericAlgorithm);})
        .def("CalculateUnitNormalsNonHistorical", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart){
            rNormalCalculationUtils.CalculateUnitNormals<ModelPart::ConditionsContainerType,false>(rModelPart);})
        .def("CalculateOnSimplex", [](NormalCalculationUtils& rNormalCalculationUtils, NormalCalculationUtils::ConditionsArrayType& rConditions,const std::size_t Dimension){
            rNormalCalculationUtils.CalculateOnSimplex(rConditions, Dimension);})
        .def("CalculateOnSimplex", [](NormalCalculationUtils& rNormalCalculationUtils, NormalCalculationUtils::ConditionsArrayType& rConditions,const std::size_t Dimension, NormalCalculationUtils::NormalVariableType& rNormalVariable){
            rNormalCalculationUtils.CalculateOnSimplex(rConditions, Dimension, rNormalVariable);})
        .def("CalculateOnSimplex", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart,const std::size_t Dimension){
            rNormalCalculationUtils.CalculateOnSimplex(rModelPart, Dimension);})
        .def("CalculateOnSimplex", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart,const std::size_t Dimension, NormalCalculationUtils::NormalVariableType& rNormalVariable){
            rNormalCalculationUtils.CalculateOnSimplex(rModelPart, Dimension, rNormalVariable);})
        .def("CalculateOnSimplex", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart){
            rNormalCalculationUtils.CalculateOnSimplex(rModelPart);})
        .def("CalculateOnSimplex", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart, const NormalCalculationUtils::NormalVariableType& rNormalVariable){
            rNormalCalculationUtils.CalculateOnSimplex(rModelPart, rNormalVariable);})
        .def("CalculateOnSimplexNonHistorical", [](NormalCalculationUtils& rNormalCalculationUtils, NormalCalculationUtils::ConditionsArrayType& rConditions,const std::size_t Dimension){
            rNormalCalculationUtils.CalculateOnSimplexNonHistorical(rConditions, Dimension);})
        .def("CalculateOnSimplexNonHistorical", [](NormalCalculationUtils& rNormalCalculationUtils, NormalCalculationUtils::ConditionsArrayType& rConditions,const std::size_t Dimension, NormalCalculationUtils::NormalVariableType& rNormalVariable){
            rNormalCalculationUtils.CalculateOnSimplexNonHistorical(rConditions, Dimension, rNormalVariable);})
        .def("CalculateOnSimplexNonHistorical", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart,const std::size_t Dimension){
            rNormalCalculationUtils.CalculateOnSimplexNonHistorical(rModelPart, Dimension);})
        .def("CalculateOnSimplexNonHistorical", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart,const std::size_t Dimension, NormalCalculationUtils::NormalVariableType& rNormalVariable){
            rNormalCalculationUtils.CalculateOnSimplexNonHistorical(rModelPart, Dimension, rNormalVariable);})
        .def("CalculateOnSimplexNonHistorical", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart){
            rNormalCalculationUtils.CalculateOnSimplexNonHistorical(rModelPart);})
        .def("CalculateOnSimplexNonHistorical", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart, const NormalCalculationUtils::NormalVariableType& rNormalVariable){
            rNormalCalculationUtils.CalculateOnSimplexNonHistorical(rModelPart, rNormalVariable);})
        .def("CalculateOnSimplex", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart,const std::size_t Dimension,Variable<double>& rVariable){rNormalCalculationUtils.CalculateOnSimplex(rModelPart, Dimension, rVariable);})
        .def("CalculateOnSimplex", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart,const std::size_t Dimension,Variable<int>& rVariable){rNormalCalculationUtils.CalculateOnSimplex(rModelPart, Dimension, rVariable);})
        .def("CalculateOnSimplex", [](NormalCalculationUtils& rNormalCalculationUtils, ModelPart& rModelPart,const std::size_t Dimension,Variable<double>& rVariable,const double Zero,const double Alpha){rNormalCalculationUtils.CalculateOnSimplex(rModelPart, Dimension, rVariable, Zero, Alpha);})
        .def("CalculateNormalShapeDerivativesOnSimplex", &NormalCalculationUtils::CalculateNormalShapeDerivativesOnSimplex)
        .def("SwapNormals", &NormalCalculationUtils::SwapNormals)
        ;

    py::class_<BodyNormalCalculationUtils > (m,"BodyNormalCalculationUtils")
        .def(py::init<>())
        .def("CalculateBodyNormals", &BodyNormalCalculationUtils::CalculateBodyNormals)
        ;

    py::class_<BodyDistanceCalculationUtils > (m,"BodyDistanceCalculationUtils")
        .def(py::init<>())
        .def("CalculateDistances2D", &BodyDistanceCalculationUtils::CalculateDistances < 2 >)
        .def("CalculateDistances3D", &BodyDistanceCalculationUtils::CalculateDistances < 3 >)
        ;

    py::class_<SignedDistanceCalculationUtils < 2 > >(m,"SignedDistanceCalculationUtils2D")
        .def(py::init<>())
        .def("CalculateDistances", &SignedDistanceCalculationUtils < 2 > ::CalculateDistances)
        .def("FindMaximumEdgeSize", &SignedDistanceCalculationUtils < 2 > ::FindMaximumEdgeSize)
        ;

    py::class_<SignedDistanceCalculationUtils < 3 > >(m,"SignedDistanceCalculationUtils3D")
        .def(py::init<>())
        .def("CalculateDistances", &SignedDistanceCalculationUtils < 3 > ::CalculateDistances)
        .def("FindMaximumEdgeSize", &SignedDistanceCalculationUtils < 3 > ::FindMaximumEdgeSize)
        ;

    py::enum_<Globals::Configuration>( m, "Configuration" )
        .value( "Initial", Globals::Configuration::Initial )
        .value( "Current", Globals::Configuration::Current )
        ;

    //brute force point locator
    py::class_<BruteForcePointLocator> (m, "BruteForcePointLocator")
        .def(py::init<ModelPart& >())
        .def("FindNode", &BruteForcePointLocator::FindNode)
        .def("FindElement", &BruteForcePointLocator::FindElement)
        .def("FindCondition", &BruteForcePointLocator::FindCondition)
        ;

    //isoprinter
    py::class_<IsosurfacePrinterApplication >(m,"IsosurfacePrinterApplication")
        .def(py::init<ModelPart& >() )
        .def("AddScalarVarIsosurface", &IsosurfacePrinterApplication::AddScalarVarIsosurface)
        .def("AddScalarVarIsosurfaceAndLower", &IsosurfacePrinterApplication::AddScalarVarIsosurfaceAndLower)
        .def("AddScalarVarIsosurfaceAndHigher", &IsosurfacePrinterApplication::AddScalarVarIsosurfaceAndHigher)
        .def("ClearData", &IsosurfacePrinterApplication::ClearData)
        .def("AddSkinConditions", &IsosurfacePrinterApplication::AddSkinConditions)
        .def("CreateNodesArray", &IsosurfacePrinterApplication::CreateNodesArray)
        ;

    //binbased locators
    py::class_< BinBasedFastPointLocator < 2 >, BinBasedFastPointLocator < 2 >::Pointer >(m,"BinBasedFastPointLocator2D")
        .def(py::init<ModelPart& >())
        .def("UpdateSearchDatabase", &BinBasedFastPointLocator < 2 > ::UpdateSearchDatabase)
        .def("UpdateSearchDatabaseAssignedSize", &BinBasedFastPointLocator < 2 > ::UpdateSearchDatabaseAssignedSize)
        .def("FindPointOnMesh", &BinBasedFastPointLocator < 2 > ::FindPointOnMeshSimplified)
        ;

    py::class_< BinBasedFastPointLocator < 3 >, BinBasedFastPointLocator < 3 >::Pointer >(m,"BinBasedFastPointLocator3D")
        .def(py::init<ModelPart&  >())
        .def("UpdateSearchDatabase", &BinBasedFastPointLocator < 3 > ::UpdateSearchDatabase)
        .def("FindPointOnMesh", &BinBasedFastPointLocator < 3 > ::FindPointOnMeshSimplified)
        .def("UpdateSearchDatabaseAssignedSize", &BinBasedFastPointLocator < 3 > ::UpdateSearchDatabaseAssignedSize)
        ;

    py::class_< BinBasedFastPointLocatorConditions < 2 > >(m,"BinBasedFastPointLocatorConditions2D")
        .def(py::init<ModelPart& >())
        .def("UpdateSearchDatabase", &BinBasedFastPointLocatorConditions < 2 > ::UpdateSearchDatabase)
        .def("UpdateSearchDatabaseAssignedSize", &BinBasedFastPointLocatorConditions < 2 > ::UpdateSearchDatabaseAssignedSize)
        .def("FindPointOnMesh", &BinBasedFastPointLocatorConditions < 2 > ::FindPointOnMeshSimplified)
        ;

    py::class_< BinBasedFastPointLocatorConditions < 3 > >(m,"BinBasedFastPointLocatorConditions3D")
        .def(py::init<ModelPart&  >())
        .def("UpdateSearchDatabase", &BinBasedFastPointLocatorConditions < 3 > ::UpdateSearchDatabase)
        .def("FindPointOnMesh", &BinBasedFastPointLocatorConditions < 3 > ::FindPointOnMeshSimplified)
        .def("UpdateSearchDatabaseAssignedSize", &BinBasedFastPointLocatorConditions < 3 > ::UpdateSearchDatabaseAssignedSize)
        ;

    py::class_< BinBasedNodesInElementLocator < 2 > >(m,"BinBasedNodesInElementLocator2D")
        .def(py::init<ModelPart& >())
        .def("UpdateSearchDatabase", &BinBasedNodesInElementLocator < 2 > ::UpdateSearchDatabase)
        .def("FindNodesInElement", &BinBasedNodesInElementLocator < 2 > ::FindNodesInElement)
        .def("UpdateSearchDatabaseAssignedSize", &BinBasedNodesInElementLocator < 2 > ::UpdateSearchDatabaseAssignedSize)
        ;

    py::class_< BinBasedNodesInElementLocator < 3 > >(m,"BinBasedNodesInElementLocator3D")
        .def(py::init<ModelPart&  >())
        .def("UpdateSearchDatabase", &BinBasedNodesInElementLocator < 3 > ::UpdateSearchDatabase)
        .def("FindNodesInElement", &BinBasedNodesInElementLocator < 3 > ::FindNodesInElement)
        .def("UpdateSearchDatabaseAssignedSize", &BinBasedNodesInElementLocator < 3 > ::UpdateSearchDatabaseAssignedSize)
        ;

    //embeded skin utilities
    py::class_< EmbeddedSkinUtility < 2 > >(m,"EmbeddedSkinUtility2D")
        .def(py::init< ModelPart&, ModelPart&, const std::string >())
        .def("GenerateSkin", &EmbeddedSkinUtility < 2 > ::GenerateSkin)
        .def("InterpolateMeshVariableToSkin", InterpolateMeshVariableToSkinArray< 2 > )
        .def("InterpolateMeshVariableToSkin", InterpolateMeshVariableToSkinDouble< 2 > )
        .def("InterpolateDiscontinuousMeshVariableToSkin", InterpolateDiscontinuousMeshVariableToSkinArray< 2 > )
        .def("InterpolateDiscontinuousMeshVariableToSkin", InterpolateDiscontinuousMeshVariableToSkinDouble< 2 > )
        ;

    py::class_< EmbeddedSkinUtility <3 > >(m,"EmbeddedSkinUtility3D")
        .def(py::init< ModelPart&, ModelPart&, const std::string >())
        .def("GenerateSkin", &EmbeddedSkinUtility < 3 > ::GenerateSkin)
        .def("InterpolateMeshVariableToSkin", InterpolateMeshVariableToSkinArray< 3 > )
        .def("InterpolateMeshVariableToSkin", InterpolateMeshVariableToSkinDouble< 3 > )
        .def("InterpolateDiscontinuousMeshVariableToSkin", InterpolateDiscontinuousMeshVariableToSkinArray< 3 > )
        .def("InterpolateDiscontinuousMeshVariableToSkin", InterpolateDiscontinuousMeshVariableToSkinDouble< 3 > )
        ;

    //Geometry tester
    py::class_< GeometryTesterUtility>(m,"GeometryTesterUtility")
        .def(py::init< >())
        .def("RunTest", &GeometryTesterUtility::RunTest)
        .def("TestTriangle2D3N", &GeometryTesterUtility::TestTriangle2D3N)
        .def("TestTriangle2D6N", &GeometryTesterUtility::TestTriangle2D6N)
        .def("TestTetrahedra3D4N", &GeometryTesterUtility::TestTetrahedra3D4N)
        .def("TestTetrahedra3D10N", &GeometryTesterUtility::TestTetrahedra3D10N)
        .def("TestHexahedra3D8N", &GeometryTesterUtility::TestHexahedra3D8N)
        .def("TestHexahedra3D27N", &GeometryTesterUtility::TestHexahedra3D27N)
        .def("TestHexahedra3D20N", &GeometryTesterUtility::TestHexahedra3D20N)
        .def("TestQuadrilateralInterface2D4N", &GeometryTesterUtility::TestQuadrilateralInterface2D4N)
        .def("TestPrismInterface3D6N", &GeometryTesterUtility::TestPrismInterface3D6N)
        .def("TestHexahedraInterface3D8N", &GeometryTesterUtility::TestHexahedraInterface3D8N)
        ;

    //Cutting utility
    py::class_<CuttingUtility >(m,"CuttingUtility")
        .def(py::init< >())
        .def("GenerateCut", &CuttingUtility::GenerateCut)
        .def("UpdateCutData", &CuttingUtility ::UpdateCutData)
        .def("AddSkinConditions", &CuttingUtility ::AddSkinConditions)
        .def("AddVariablesToCutModelPart", &CuttingUtility::AddVariablesToCutModelPart )
        .def("FindSmallestEdge", &CuttingUtility ::FindSmallestEdge)
        ;

    // Interval utility
    #define KRATOS_DEFINE_INTERVAL_UTILITY_PYTHON_BINDINGS(CLASS_NAME) \
        py::class_<CLASS_NAME>(m, #CLASS_NAME)                         \
            .def(py::init<Parameters>())                               \
            .def("GetIntervalBegin", &CLASS_NAME ::GetIntervalBegin)   \
            .def("GetIntervalEnd", &CLASS_NAME ::GetIntervalEnd)       \
            .def("IsInInterval", &CLASS_NAME ::IsInInterval)           \
            ;

    KRATOS_DEFINE_INTERVAL_UTILITY_PYTHON_BINDINGS(IntervalUtility)

    KRATOS_DEFINE_INTERVAL_UTILITY_PYTHON_BINDINGS(DiscreteIntervalUtility)

    #undef KRATOS_DEFINE_INTERVAL_UTILITY_PYTHON_BINDINGS

    //particle convect utility
    py::class_<ParticleConvectUtily<2> >(m,"ParticleConvectUtily2D")
        .def(py::init< BinBasedFastPointLocator < 2 >::Pointer >())
        .def("MoveParticles_Substepping", &ParticleConvectUtily<2>::MoveParticles_Substepping)
        .def("MoveParticles_RK4", &ParticleConvectUtily<2>::MoveParticles_RK4)
        ;

    py::class_<ParticleConvectUtily<3> >(m,"ParticleConvectUtily3D")
        .def(py::init< BinBasedFastPointLocator < 3 >::Pointer >())
        .def("MoveParticles_Substepping", &ParticleConvectUtily<3>::MoveParticles_Substepping)
        .def("MoveParticles_RK4", &ParticleConvectUtily<3>::MoveParticles_RK4)
        ;

    // Delaunator utilities
    auto mod_delaunator = m.def_submodule("CreateTriangleMeshFromNodes");
    mod_delaunator.def("CreateTriangleMeshFromNodes",&DelaunatorUtilities::CreateTriangleMeshFromNodes);

    // GeometricalTransformationUtilities
    auto mod_geom_trans_utils = m.def_submodule("GeometricalTransformationUtilities");
    mod_geom_trans_utils.def("CalculateTranslationMatrix", &GeometricalTransformationUtilities::CalculateTranslationMatrix );
    mod_geom_trans_utils.def("CalculateRotationMatrix", &GeometricalTransformationUtilities::CalculateRotationMatrix );

    // MLS shape functions utility
    py::class_<MLSShapeFunctionsUtility>(m,"MLSShapeFunctionsUtility")
        .def_static("CalculateShapeFunctions2D", &MLSShapeFunctionsUtility::CalculateShapeFunctions<2>)
        .def_static("CalculateShapeFunctions3D", &MLSShapeFunctionsUtility::CalculateShapeFunctions<3>)
        .def_static("CalculateShapeFunctionsAndGradients2D", &MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<2>)
        .def_static("CalculateShapeFunctionsAndGradients3D", &MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<3>)
        ;
}

} // namespace Python.
} // Namespace Kratos
