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
#include <pybind11/numpy.h>
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
#include "utilities/rbf_shape_functions_utility.h"
#include "utilities/assign_master_slave_constraints_to_neighbours_utility.h"

namespace Kratos::Python {

template<std::size_t TDim>
    std::tuple<pybind11::array_t<int>, pybind11::array_t<int>, pybind11::array_t<double> >
    VectorizedFindHelper(BinBasedFastPointLocator < TDim >& rSelf, const Matrix& rCoords )
        {
            pybind11::array_t<int> element_ids_container(rCoords.size1());
            auto element_ids = element_ids_container.mutable_unchecked<1>();

            unsigned int nnodes=TDim+1; //only for triangles and tets
            unsigned int ncoords=rCoords.size1();
            KRATOS_ERROR_IF(rCoords.size2() != TDim) << "the expected size of coordinates is not correct ";

            pybind11::array_t<int> node_ids_container({ncoords,nnodes}); //this is designed to work with Tets
            pybind11::array_t<double> shape_function_values_container({ncoords,nnodes});
            auto node_ids = node_ids_container.mutable_unchecked<2>();
            auto shape_function_values = shape_function_values_container.mutable_unchecked<2>();

            Vector N(nnodes);
            IndexPartition(ncoords).for_each(N,[&](unsigned int i, Vector& N)
            {
                Element::Pointer pelem;
                const bool is_found = rSelf.FindPointOnMeshSimplified(row(rCoords,i),N,pelem);
                const auto& r_geom = pelem->GetGeometry();
                if(is_found){
                    for(unsigned int k = 0; k<nnodes; ++k){
                        node_ids(i,k) = r_geom[k].Id();
                        shape_function_values(i,k) = N[k];
                        element_ids[i] = pelem->Id();
                    }
                }
                else
                {
                    for(unsigned int k = 0; k<nnodes; ++k){
                        node_ids(i,k) = 1; //deliberately set to 1
                        shape_function_values(i,k) = 0.0;
                        element_ids[i] = -1;
                    }
                }
            });
            return std::tuple{element_ids_container, node_ids_container, shape_function_values_container};
        }

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

// MLS shape functions utility
void AuxiliaryCalculateMLSShapeFunctions(
    const std::size_t Dim,
    const std::size_t Order,
    const Matrix& rPoints,
    const array_1d<double,3>& rX,
    const double h,
    Vector& rN)
{
    switch (Dim)
    {
    case 2:
        switch (Order)
        {
        case 1:
            MLSShapeFunctionsUtility::CalculateShapeFunctions<2,1>(rPoints, rX, h, rN);
            break;
        case 2:
            MLSShapeFunctionsUtility::CalculateShapeFunctions<2,2>(rPoints, rX, h, rN);
            break;
        default:
            KRATOS_ERROR << "Wrong order: " << Order << ". Only \'1\' and \'2\' are supported." << std::endl;
            break;
        }
    case 3:
        switch (Order)
        {
        case 1:
            MLSShapeFunctionsUtility::CalculateShapeFunctions<3,1>(rPoints, rX, h, rN);
            break;
        case 2:
            MLSShapeFunctionsUtility::CalculateShapeFunctions<3,2>(rPoints, rX, h, rN);
            break;
        default:
            KRATOS_ERROR << "Wrong order: " << Order << ". Only \'1\' and \'2\' are supported." << std::endl;
            break;
        }
    default:
        KRATOS_ERROR << "Wrong dimension: " << Dim << std::endl;
        break;
    }
}

void AuxiliaryCalculateMLSShapeFunctionsAndGradients(
    const std::size_t Dim,
    const std::size_t Order,
    const Matrix& rPoints,
    const array_1d<double,3>& rX,
    const double h,
    Vector& rN,
    Matrix& rDNDX)
{
    if (Dim == 2) {
        switch (Order)
        {
        case 1:
            MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<2,1>(rPoints, rX, h, rN, rDNDX);
            break;
        case 2:
            MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<2,2>(rPoints, rX, h, rN, rDNDX);
            break;
        default:
            KRATOS_ERROR << "Wrong order: " << Order << ". Only \'1\' and \'2\' are supported." << std::endl;
            break;
        }
    } else if (Dim == 3) {
        switch (Order)
        {
        case 1:
            MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<3,1>(rPoints, rX, h, rN, rDNDX);
            break;
        case 2:
            MLSShapeFunctionsUtility::CalculateShapeFunctionsAndGradients<3,2>(rPoints, rX, h, rN, rDNDX);
            break;
        default:
            KRATOS_ERROR << "Wrong order: " << Order << ". Only \'1\' and \'2\' are supported." << std::endl;
            break;
        }
    } else {
        KRATOS_ERROR << "Wrong dimension: " << Dim << std::endl;
    }
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

    // IsosurfacePrinterApplication
    py::class_<IsosurfacePrinterApplication >(m,"IsosurfacePrinterApplication")
        .def(py::init<ModelPart& >() )
        .def("AddScalarVarIsosurface", &IsosurfacePrinterApplication::AddScalarVarIsosurface)
        .def("AddScalarVarIsosurfaceAndLower", &IsosurfacePrinterApplication::AddScalarVarIsosurfaceAndLower)
        .def("AddScalarVarIsosurfaceAndHigher", &IsosurfacePrinterApplication::AddScalarVarIsosurfaceAndHigher)
        .def("ClearData", &IsosurfacePrinterApplication::ClearData)
        .def("AddSkinConditions", &IsosurfacePrinterApplication::AddSkinConditions)
        .def("CreateNodesArray", &IsosurfacePrinterApplication::CreateNodesArray)
        ;

    // BinBasedFastPointLocator
    py::class_< BinBasedFastPointLocator < 2 >, BinBasedFastPointLocator < 2 >::Pointer >(m,"BinBasedFastPointLocator2D")
        .def(py::init<ModelPart& >())
        .def("UpdateSearchDatabase", &BinBasedFastPointLocator < 2 > ::UpdateSearchDatabase)
        .def("UpdateSearchDatabaseAssignedSize", &BinBasedFastPointLocator < 2 > ::UpdateSearchDatabaseAssignedSize)
        .def("FindPointOnMesh", [](BinBasedFastPointLocator < 2 >& rSelf, const array_1d<double,3>& rCoords ){
            Element::Pointer p_elem = nullptr;
            Vector N;
            const bool is_found = rSelf.FindPointOnMeshSimplified(rCoords,N,p_elem);
            return std::tuple<bool,Vector,Element::Pointer>{is_found,N,p_elem};
        })
        .def("VectorizedFind", &VectorizedFindHelper<2>)
        ;

    py::class_< BinBasedFastPointLocator < 3 >, BinBasedFastPointLocator < 3 >::Pointer >(m,"BinBasedFastPointLocator3D")
        .def(py::init<ModelPart&  >())
        .def("UpdateSearchDatabase", &BinBasedFastPointLocator < 3 > ::UpdateSearchDatabase)
        .def("UpdateSearchDatabaseAssignedSize", &BinBasedFastPointLocator < 3 > ::UpdateSearchDatabaseAssignedSize)
        .def("FindPointOnMesh", [](BinBasedFastPointLocator < 3 >& rSelf, const array_1d<double,3>& rCoords ){
            Element::Pointer p_elem = nullptr;
            Vector N;
            const bool is_found = rSelf.FindPointOnMeshSimplified(rCoords,N,p_elem);
            return std::tuple<bool,Vector,Element::Pointer>{is_found,N,p_elem};
        })
        .def("VectorizedFind", &VectorizedFindHelper<3>)
        ;

    py::class_< BinBasedFastPointLocatorConditions < 2 > >(m,"BinBasedFastPointLocatorConditions2D")
        .def(py::init<ModelPart& >())
        .def("UpdateSearchDatabase", &BinBasedFastPointLocatorConditions < 2 > ::UpdateSearchDatabase)
        .def("UpdateSearchDatabaseAssignedSize", &BinBasedFastPointLocatorConditions < 2 > ::UpdateSearchDatabaseAssignedSize)
        .def("FindPointOnMesh", [](BinBasedFastPointLocatorConditions < 2 >& rSelf, const array_1d<double,3>& rCoords ){
            Condition::Pointer p_cond = nullptr;
            Vector N;
            const bool is_found = rSelf.FindPointOnMeshSimplified(rCoords,N,p_cond);
            return std::tuple<bool,Vector,Condition::Pointer>{is_found,N,p_cond};
        })
        ;

    py::class_< BinBasedFastPointLocatorConditions < 3 > >(m,"BinBasedFastPointLocatorConditions3D")
        .def(py::init<ModelPart&  >())
        .def("UpdateSearchDatabase", &BinBasedFastPointLocatorConditions < 3 > ::UpdateSearchDatabase)
        .def("UpdateSearchDatabaseAssignedSize", &BinBasedFastPointLocatorConditions < 3 > ::UpdateSearchDatabaseAssignedSize)
        .def("FindPointOnMesh", [](BinBasedFastPointLocatorConditions < 3 >& rSelf, const array_1d<double,3>& rCoords ){
            Condition::Pointer p_cond = nullptr;
            Vector N;
            const bool is_found = rSelf.FindPointOnMeshSimplified(rCoords,N,p_cond);
            return std::tuple<bool,Vector,Condition::Pointer>{is_found,N,p_cond};
        })
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




    //embedded skin utilities
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

    //interval utility
    py::class_<IntervalUtility >(m,"IntervalUtility")
        .def(py::init<Parameters >())
        .def("GetIntervalBegin", &IntervalUtility::GetIntervalBegin)
        .def("GetIntervalEnd", &IntervalUtility::GetIntervalEnd)
        .def("IsInInterval", &IntervalUtility ::IsInInterval)
        .def("__str__", PrintObject<IntervalUtility>)
        ;

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
        .def_static("CalculateShapeFunctions", &AuxiliaryCalculateMLSShapeFunctions)
        .def_static("CalculateShapeFunctionsAndGradients", &AuxiliaryCalculateMLSShapeFunctionsAndGradients)
        ;

    // Radial Basis FUnctions utility
    using DenseQRPointerType = typename RBFShapeFunctionsUtility::DenseQRPointerType;
    py::class_<RBFShapeFunctionsUtility>(m,"RBFShapeFunctionsUtility")
        .def_static("CalculateShapeFunctions", [](const Matrix& rPoints, const array_1d<double,3>& rX, Vector& rN){
            return RBFShapeFunctionsUtility::CalculateShapeFunctions(rPoints, rX, rN);})
        .def_static("CalculateShapeFunctions", [](const Matrix& rPoints, const array_1d<double,3>& rX, Vector& rN, DenseQRPointerType pDenseQR){
            return RBFShapeFunctionsUtility::CalculateShapeFunctions(rPoints, rX, rN, pDenseQR);})
        .def_static("CalculateShapeFunctions", [](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN){
            return RBFShapeFunctionsUtility::CalculateShapeFunctions(rPoints, rX, h, rN);})
        .def_static("CalculateShapeFunctions", [](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN, DenseQRPointerType pDenseQR){
            return RBFShapeFunctionsUtility::CalculateShapeFunctions(rPoints, rX, h, rN, pDenseQR);})
        ;

    // Radial Node Search and MasterSlaveConstraints Assignation Utility
    using NodesContainerType = typename AssignMasterSlaveConstraintsToNeighboursUtility::NodesContainerType;
    py::class_<AssignMasterSlaveConstraintsToNeighboursUtility>(m, "AssignMasterSlaveConstraintsToNeighboursUtility")
        .def(py::init<ModelPart::NodesContainerType&>())
        .def("AssignMasterSlaveConstraintsToNodes", [](AssignMasterSlaveConstraintsToNeighboursUtility& rAssignMasterSlaveConstraintsToNeighboursUtility, NodesContainerType pNodes, double const Radius, ModelPart& rComputingModelPart, const std::vector<std::reference_wrapper<const Kratos::Variable<double>>>& rVariableList, double const MinNumOfNeighNodes){
            return rAssignMasterSlaveConstraintsToNeighboursUtility.AssignMasterSlaveConstraintsToNodes(pNodes, Radius, rComputingModelPart, rVariableList, MinNumOfNeighNodes);})
        ;
}

} // namespace Kratos::Python.
