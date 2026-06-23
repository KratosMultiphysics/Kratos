// // KRATOS  ___|  |                   |                   |
// //       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
// //             | |   |    |   | (    |   |   | |   (   | |
// //       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
// //
// //  License:         BSD License
// //                   license: StructuralMechanicsApplication/license.txt
// //
// //  Main authors:    Andi Makarim Katili
// //

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "mpm_application_variables.h"
#include "containers/model.h"


// Application includes
#include "custom_utilities/material_point_search_utility.h"
#include "utilities/quadrature_points_utility.h"
#include "custom_utilities/mpm_volume_sum_utility.h"
#include "custom_processes/material_point_erase_process.h"
namespace Kratos
{
namespace Testing
{

void ResetGridVolume(ModelPart& rBackgroundModelPart)
{
    for (auto& r_element : rBackgroundModelPart.Elements()) {
        r_element.GetGeometry().SetValue(TOTAL_MP_VOLUME, 0.0);
    }
}

void GetVolumeVector(ModelPart& rBackgroundModelPart, Vector& rVolumeVector)
{
    if (rVolumeVector.size() != rBackgroundModelPart.NumberOfElements()) rVolumeVector.resize(rBackgroundModelPart.NumberOfElements(),false);

    IndexType i = 0;
    for (auto& r_element : rBackgroundModelPart.Elements()) {
        rVolumeVector[i] = r_element.GetValue(TOTAL_MP_VOLUME);
        i++;
    }
}

void GenerateGridnodes(
    ModelPart& rBackgroundModelPart,
    const double dx,
    const double dy,
    const double dz,
    const size_t NumberOfRows,
    const size_t NumberOfColumns,
    const size_t NumberOfLayers = 1)
{
    // Nodes
    IndexType point_index = 1;
    for (size_t layer = 0; layer < NumberOfLayers; ++layer){
        for (size_t row = 0; row < NumberOfRows; ++row){
            for (size_t col = 0; col < NumberOfColumns; ++col){
                rBackgroundModelPart.CreateNewNode(point_index, double(col) * dx, double(row) * dy, double(layer) * dz);
                point_index += 1;
            }
        }
    }
}

void PrepareBackgroundModelPart(ModelPart& rBackgroundModelPart, const int caseIndex )
{
    rBackgroundModelPart.AddNodalSolutionStepVariable(TOTAL_MP_VOLUME);

    // General grid scheme:
    //  13--14--15--16
    //  |  /|  /|  /|
    //  | / | / | / |
    //  |/  |/  |/  |
    //  9---10- 11--12 (if skew then x += 0.1 in this row)
    //  |  /|  /|  /|
    //  | / | / | / |
    //  |/  |/  |/  |
    //  5---6---7---8 (if skew then x += 0.1 in this row)
    //  |  /|  /|  /|
    //  | / | / | / |
    //  |/  |/  |/  |
    //  1---2---3---4

    std::string element;
    switch (caseIndex)
    {
    case 0:
        element = "Element2D4N";
        break;
    case 10:
        element = "Element2D3N";
        break;
    case 20:
        element = "Element3D8N";
        break;
    case 30:
        element = "Element3D4N";
        break;
    default:
        break;
    }

    if (element == "Element2D4N")
    {
        // Nodes
        GenerateGridnodes(rBackgroundModelPart, 1.0, 1.0, 1.0, 4, 4, 1);
        // Grid elements
        rBackgroundModelPart.CreateNewElement( element, 1, { 1, 2, 6, 5 }, nullptr);
        rBackgroundModelPart.CreateNewElement( element, 2, { 2, 3, 7, 6 }, nullptr);
        rBackgroundModelPart.CreateNewElement( element, 3, { 3, 4, 8, 7 }, nullptr);

        rBackgroundModelPart.CreateNewElement( element, 4, { 5, 6, 10, 9 }, nullptr);
        rBackgroundModelPart.CreateNewElement( element, 5, { 6, 7, 11, 10 }, nullptr);
        rBackgroundModelPart.CreateNewElement( element, 6, { 7, 8, 12, 11 }, nullptr);

        rBackgroundModelPart.CreateNewElement( element, 7, { 9, 10, 14, 13 }, nullptr);
        rBackgroundModelPart.CreateNewElement( element, 8, { 10, 11, 15, 14 }, nullptr);
        rBackgroundModelPart.CreateNewElement( element, 9, { 11, 12, 16, 15 }, nullptr);
    }
    else if (element == "Element2D3N")
    {
        // Nodes
        GenerateGridnodes(rBackgroundModelPart, 1.0, 1.0, 1.0, 4, 4, 1);

        // Grid elements
        rBackgroundModelPart.CreateNewElement(element, 1, { 1, 2, 6 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 2, { 2, 3, 7 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 3, { 3, 4, 8 }, nullptr);

        rBackgroundModelPart.CreateNewElement(element, 4, { 1, 6, 5 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 5, { 2, 7, 6 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 6, { 3, 8, 7 }, nullptr);

        rBackgroundModelPart.CreateNewElement(element, 7, { 5, 6, 10 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 8, { 6, 7, 11 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 9, { 7, 8, 12 }, nullptr);

        rBackgroundModelPart.CreateNewElement(element, 10, { 5, 10, 9 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 11, { 6, 11, 10 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 12, { 7, 12, 11 }, nullptr);

        rBackgroundModelPart.CreateNewElement(element, 13, { 9, 10, 14 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 14, { 10, 11, 15 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 15, { 11, 12, 16 }, nullptr);

        rBackgroundModelPart.CreateNewElement(element, 16, { 9, 14, 13 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 17, { 10, 15, 14 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 18, { 11, 16, 15 }, nullptr);
    }
    else if (element == "Element3D8N")
    {
        // Nodes
        GenerateGridnodes(rBackgroundModelPart, 1.0, 1.0, 1.0, 3, 3, 3);
        KRATOS_WATCH(rBackgroundModelPart.Nodes())

        rBackgroundModelPart.CreateNewElement(element, 1, { 1, 2, 5, 4, 10, 11, 14, 13}, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 2, { 2, 3, 6, 5, 11, 12, 15, 14}, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 3, { 5, 6, 9, 8, 14, 15, 18, 17}, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 4, { 4, 5, 8, 7, 13, 14, 17, 16}, nullptr);
        
        rBackgroundModelPart.CreateNewElement(element, 5, { 10, 11, 14, 13, 19, 20, 23, 22}, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 6, { 11, 12, 15, 14, 20, 21, 24, 23}, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 7, { 14, 15, 18, 17, 23, 24, 27, 26}, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 8, { 13, 14, 17, 16, 22, 23, 26, 25}, nullptr);
    }
    else if (element == "Element3D4N")
    {
        // Nodes
        rBackgroundModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
        rBackgroundModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
        rBackgroundModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
        rBackgroundModelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
        rBackgroundModelPart.CreateNewNode(5, 1.0, 1.0, 0.0);
        rBackgroundModelPart.CreateNewNode(6, 1.0, 0.0, 1.0);
        rBackgroundModelPart.CreateNewNode(7, 0.0, 1.0, 1.0);
        rBackgroundModelPart.CreateNewNode(8, 1.0, 1.0, 1.0);
        rBackgroundModelPart.CreateNewNode(9, 2.0, 0.0, 0.0);
        rBackgroundModelPart.CreateNewNode(10, 0.0, 2.0, 0.0);
        rBackgroundModelPart.CreateNewNode(11, 0.0, 0.0, 2.0);
        rBackgroundModelPart.CreateNewNode(12, 2.0, 1.0, 0.0);
        rBackgroundModelPart.CreateNewNode(13, 1.0, 2.0, 0.0);
        rBackgroundModelPart.CreateNewNode(14, 2.0, 0.0, 1.0);
        rBackgroundModelPart.CreateNewNode(15, 1.0, 0.0, 2.0);
        rBackgroundModelPart.CreateNewNode(16, 0.0, 2.0, 1.0);
        rBackgroundModelPart.CreateNewNode(17, 0.0, 1.0, 2.0);
        rBackgroundModelPart.CreateNewNode(18, 2.0, 1.0, 1.0);
        rBackgroundModelPart.CreateNewNode(19, 1.0, 2.0, 1.0);
        rBackgroundModelPart.CreateNewNode(20, 1.0, 1.0, 2.0);
        rBackgroundModelPart.CreateNewNode(21, 2.0, 2.0, 0.0);
        rBackgroundModelPart.CreateNewNode(22, 2.0, 0.0, 2.0);
        rBackgroundModelPart.CreateNewNode(23, 0.0, 2.0, 2.0);
        rBackgroundModelPart.CreateNewNode(24, 2.0, 2.0, 1.0);
        rBackgroundModelPart.CreateNewNode(25, 2.0, 1.0, 2.0);
        rBackgroundModelPart.CreateNewNode(26, 1.0, 2.0, 2.0);
        rBackgroundModelPart.CreateNewNode(27, 2.0, 2.0, 2.0);

        // Grid elements
        rBackgroundModelPart.CreateNewElement(element,  1, { 1,  2,  5, 4 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element,  2, { 1,  2,  5, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element,  3, { 4,  7,  8, 1 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element,  4, { 1,  3,  7, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element,  5, { 1,  4,  6, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element,  6, { 3,  5,  8, 1 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element,  7, {12, 18,  8, 9 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element,  8, { 9, 12,  5, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element,  9, {14,  6,  8, 9 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 10, { 9,  2,  6, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 11, { 9, 14, 18, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 12, { 2,  5,  8, 9 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 13, { 7, 17, 11, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 14, { 8,  7,  4, 1 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 15, {20, 15, 11, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 16, { 8,  6, 15, 1 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 17, { 8, 20, 17, 1 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 18, { 6,  4, 11, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 19, { 6, 15, 22, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 20, { 8,  6, 14, 2 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 21, {20, 25, 22, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 22, { 8, 18, 25, 2 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 23, { 8, 20, 15, 2 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 24, {18, 14, 22, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 25, { 3,  7,  8, 0 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 26, {10,  3,  5, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 27, {16, 19,  8, 0 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 28, {10, 13, 19, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 29, {10, 16,  7, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 30, {13,  5,  8, 0 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 31, {13, 19,  8, 1 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 32, {21, 13,  5, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 33, {24, 18,  8, 1 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 34, {21, 12, 18, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 35, {21, 24, 19, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 36, {12,  5,  8, 1 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 37, {19, 26, 23, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 38, { 8, 19, 16, 3 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 39, {20, 17, 23, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 40, { 8,  7, 17, 3 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 41, { 8, 20, 26, 3 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 42, { 7, 16, 23, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 43, {18, 25, 27, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 44, { 8, 18, 24, 7 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 45, {20, 26, 27, 8 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 46, { 8, 19, 26, 7 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 47, { 8, 20, 25, 7 }, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 48, {19, 24, 27, 8 }, nullptr);
    }
    ResetGridVolume(rBackgroundModelPart);
}

void PrepareMPMModelPart(ModelPart& rMPMModelPart, Geometry<Node>::Pointer pInitialGeometry,GeometryData::IntegrationMethod rIntegrationMethod)
{
    const Element& new_element = KratosComponents<Element>::Get("MPMUpdatedLagrangian");
    SizeType new_element_id = rMPMModelPart.NumberOfElements();
    const Geometry<Node>::IntegrationPointsArrayType& rThisIntegrationPoints = pInitialGeometry->IntegrationPoints(rIntegrationMethod);
    
    std::vector<array_1d<double, 3>> coordinate = { ZeroVector(3) };
    Matrix shape_functions_values = pInitialGeometry->ShapeFunctionsValues(rIntegrationMethod);

    std::vector<double> mp_volume(1);
    mp_volume[0] = 1.0; 
    auto p_elem_prop = rMPMModelPart.pGetProperties(0);
    
    for ( IndexType gp_number = 0; gp_number < rThisIntegrationPoints.size(); ++gp_number )
    {
        coordinate[0].clear();
        
        // coordinate calculation
        for (unsigned int dimension = 0; dimension < pInitialGeometry->WorkingSpaceDimension(); dimension++)
        {
            for (unsigned int j = 0; j < pInitialGeometry->size(); j++)
            {
                coordinate[0][dimension] = coordinate[0][dimension] + shape_functions_values(gp_number, j) * (*pInitialGeometry)[j].Coordinates()[dimension];            
            }
        }
        auto p_quadrature_point_geometry = CreateQuadraturePointsUtility<Node>::CreateFromCoordinates( pInitialGeometry, coordinate[0], mp_volume[0]);
        
        Element::Pointer p_element = new_element.Create(new_element_id, p_quadrature_point_geometry, p_elem_prop);

        p_element->SetValuesOnIntegrationPoints(MP_COORD, coordinate, rMPMModelPart.GetProcessInfo());
        p_element->SetValuesOnIntegrationPoints(MP_VOLUME, mp_volume, rMPMModelPart.GetProcessInfo());
        rMPMModelPart.AddElement(p_element);
        new_element_id++;
        mp_volume[0] += 0.5; 
        
    }
    // throw std::runtime_error("\n\n"
    // "         |\\__/,|   (`\\\n"
    // "       _.|o o  |_   ) )\n"
    // "     -(((---(((--------\n"
    // "\n"
    // "Runtime error: Cat refuses to continue!\n"
    // );
}


void MoveMaterialPoints(ModelPart& rMPMModelPart, const Vector& rDisplacement)
{
    for (auto& r_element : rMPMModelPart.Elements())
    {
        std::vector<array_1d<double, 3>> mp_coordinate;
        r_element.CalculateOnIntegrationPoints(MP_COORD, mp_coordinate, rMPMModelPart.GetProcessInfo());
        mp_coordinate[0] = mp_coordinate[0] + rDisplacement;
        r_element.SetValuesOnIntegrationPoints(MP_COORD, mp_coordinate, rMPMModelPart.GetProcessInfo());
    }
}

template<std::size_t TDimension>
void SearchAndCalculateVolume(ModelPart& rBackgroundModelPart, ModelPart& rMPMModelPart)
{
    // Search MP
    int max_number_of_search_results = 1000;
    double searching_tolerance = 1e-6;
    MPMSearchElementUtility::SearchElement<TDimension>(rBackgroundModelPart, rMPMModelPart, max_number_of_search_results, searching_tolerance);

    // Calculate mp volume at grid 
    ResetGridVolume(rBackgroundModelPart);
    MPMVolumeSumUtility::AddModelPartMPMVolumeIntoGrid(rMPMModelPart);
}

KRATOS_TEST_CASE_IN_SUITE(TestVolumeSumUtility2D4N, KratosMPMFastSuite)
{
    // Model Preparation
    Model current_model;
    ModelPart& r_mpm_model_part = current_model.CreateModelPart("MPMModelPart");
    ModelPart& r_background_model_part = current_model.CreateModelPart("BackgroundModelPart");
    const IndexType grid_case = 0; //2D4N
    PrepareBackgroundModelPart(r_background_model_part, grid_case);
    
    auto initial_geometry = r_background_model_part.GetElement(1).pGetGeometry();
    auto rThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_4; // 16 MP

    r_mpm_model_part.CreateNewProperties(0);
    PrepareMPMModelPart(r_mpm_model_part, initial_geometry, rThisIntegrationMethod);
    
    // Initial location
    Vector grid_volume_vector;
    SearchAndCalculateVolume<2>(r_background_model_part, r_mpm_model_part);
    GetVolumeVector(r_background_model_part, grid_volume_vector);
    
    Vector ref_grid_volume_vector = ZeroVector(9);
    ref_grid_volume_vector[0] = 76;
    KRATOS_EXPECT_VECTOR_EQ(grid_volume_vector, ref_grid_volume_vector);
    

    // Move material points and recalculate mp volume at grid
    Vector mp_displacement(3);
    mp_displacement[0] = 0.1;
    mp_displacement[1] = 0.1;
    mp_displacement[2] = 0.0;
    MoveMaterialPoints(r_mpm_model_part, mp_displacement);

    SearchAndCalculateVolume<2>(r_background_model_part, r_mpm_model_part);
    GetVolumeVector(r_background_model_part, grid_volume_vector);
    ref_grid_volume_vector = ZeroVector(9);
    ref_grid_volume_vector[0] = 31.5;
    ref_grid_volume_vector[1] = 22.5;
    ref_grid_volume_vector[3] = 13.5;
    ref_grid_volume_vector[4] = 8.5;

    KRATOS_EXPECT_VECTOR_EQ(grid_volume_vector, ref_grid_volume_vector);


    // Move material points and recalculate mp volume at grid
    mp_displacement[0] = 1.4;
    mp_displacement[1] = 1.4;
    mp_displacement[2] = 0.0;
    MoveMaterialPoints(r_mpm_model_part, mp_displacement);

    SearchAndCalculateVolume<2>(r_background_model_part, r_mpm_model_part);
    GetVolumeVector(r_background_model_part, grid_volume_vector);
    ref_grid_volume_vector = ZeroVector(9);
    ref_grid_volume_vector[4] = 9.0;
    ref_grid_volume_vector[5] = 25.0;
    ref_grid_volume_vector[7] = 13.0;
    ref_grid_volume_vector[8] = 29.0;

    KRATOS_EXPECT_VECTOR_EQ(grid_volume_vector, ref_grid_volume_vector);

    // Delete, search material points and recalculate mp volume at grid
    auto& mp_10 = r_mpm_model_part.GetElement(10);
    mp_10.GetGeometry().clear();
    mp_10.Reset(ACTIVE);
    mp_10.Set(TO_ERASE);
    MaterialPointEraseProcess(r_mpm_model_part).Execute();

    SearchAndCalculateVolume<2>(r_background_model_part, r_mpm_model_part);
    GetVolumeVector(r_background_model_part, grid_volume_vector);
    ref_grid_volume_vector[8] = 23.0;
    KRATOS_EXPECT_VECTOR_EQ(grid_volume_vector, ref_grid_volume_vector);
}

KRATOS_TEST_CASE_IN_SUITE(TestVolumeSumUtility2D3N, KratosMPMFastSuite)
{
    Model current_model;
    ModelPart& r_mpm_model_part = current_model.CreateModelPart("MPMModelPart");
    ModelPart& r_background_model_part = current_model.CreateModelPart("BackgroundModelPart");
    const IndexType grid_case = 10; //2D3N
    PrepareBackgroundModelPart(r_background_model_part, grid_case);
    
    auto initial_geometry1 = r_background_model_part.GetElement(1).pGetGeometry();
    auto initial_geometry2 = r_background_model_part.GetElement(4).pGetGeometry();
    auto rThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;

    r_mpm_model_part.CreateNewProperties(0);
    PrepareMPMModelPart(r_mpm_model_part, initial_geometry1, rThisIntegrationMethod);
    PrepareMPMModelPart(r_mpm_model_part, initial_geometry2, rThisIntegrationMethod);

    // Initial location
    Vector grid_volume_vector;
    SearchAndCalculateVolume<2>(r_background_model_part, r_mpm_model_part);
    GetVolumeVector(r_background_model_part, grid_volume_vector);

    Vector ref_grid_volume_vector = ZeroVector(18);
    ref_grid_volume_vector[0] = 4.5;
    ref_grid_volume_vector[3] = 4.5;
    KRATOS_EXPECT_VECTOR_EQ(grid_volume_vector, ref_grid_volume_vector);
    

    // Move material points and recalculate mp volume at grid
    Vector mp_displacement(3);
    mp_displacement[0] = 0.2;
    mp_displacement[1] = 0.2;
    mp_displacement[2] = 0.0;
    MoveMaterialPoints(r_mpm_model_part, mp_displacement);

    SearchAndCalculateVolume<2>(r_background_model_part, r_mpm_model_part);
    GetVolumeVector(r_background_model_part, grid_volume_vector);
    ref_grid_volume_vector = ZeroVector(18);
    ref_grid_volume_vector[0] = 1.0;
    ref_grid_volume_vector[3] = 1.0;
    ref_grid_volume_vector[4] = 3.5;
    ref_grid_volume_vector[6] = 3.5;

    KRATOS_EXPECT_VECTOR_EQ(grid_volume_vector, ref_grid_volume_vector);


    // Move material points and recalculate mp volume at grid
    mp_displacement[0] = 0.2;
    mp_displacement[1] = 0.2;
    mp_displacement[2] = 0.0;
    MoveMaterialPoints(r_mpm_model_part, mp_displacement);

    SearchAndCalculateVolume<2>(r_background_model_part, r_mpm_model_part);
    GetVolumeVector(r_background_model_part, grid_volume_vector);
    ref_grid_volume_vector = ZeroVector(18);
    ref_grid_volume_vector[0] = 1.0;
    ref_grid_volume_vector[3] = 1.0;
    ref_grid_volume_vector[4] = 1.5;
    ref_grid_volume_vector[6] = 2.0;
    ref_grid_volume_vector[7] = 2.0;
    ref_grid_volume_vector[10] = 1.5;

    KRATOS_EXPECT_VECTOR_EQ(grid_volume_vector, ref_grid_volume_vector);

    // Delete, search material points and recalculate mp volume at grid
    auto& mp_2 = r_mpm_model_part.GetElement(2);
    mp_2.GetGeometry().clear();
    mp_2.Reset(ACTIVE);
    mp_2.Set(TO_ERASE);
    MaterialPointEraseProcess(r_mpm_model_part).Execute();

    SearchAndCalculateVolume<2>(r_background_model_part, r_mpm_model_part);
    GetVolumeVector(r_background_model_part, grid_volume_vector);
    ref_grid_volume_vector[7] = 0.0;
    KRATOS_EXPECT_VECTOR_EQ(grid_volume_vector, ref_grid_volume_vector);
}

KRATOS_TEST_CASE_IN_SUITE(TestVolumeSumUtility3D8N, KratosMPMFastSuite)
{
    Model current_model;
    ModelPart& r_mpm_model_part = current_model.CreateModelPart("MPMModelPart");
    ModelPart& r_background_model_part = current_model.CreateModelPart("BackgroundModelPart");
    const IndexType grid_case = 20; //3D8N
    PrepareBackgroundModelPart(r_background_model_part, grid_case);
    
    auto initial_geometry1 = r_background_model_part.GetElement(1).pGetGeometry();
    auto rThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;

    r_mpm_model_part.CreateNewProperties(0);
    PrepareMPMModelPart(r_mpm_model_part, initial_geometry1, rThisIntegrationMethod);

    // Initial location
    Vector grid_volume_vector;
    SearchAndCalculateVolume<3>(r_background_model_part, r_mpm_model_part);
    GetVolumeVector(r_background_model_part, grid_volume_vector);

    Vector ref_grid_volume_vector = ZeroVector(8);
    ref_grid_volume_vector[0] = 22.0;
    KRATOS_EXPECT_VECTOR_EQ(grid_volume_vector, ref_grid_volume_vector);
    
    // Move material points and recalculate mp volume at grid
    Vector mp_displacement(3);
    mp_displacement[0] = 0.0;
    mp_displacement[1] = 0.0;
    mp_displacement[2] = 1.0;
    MoveMaterialPoints(r_mpm_model_part, mp_displacement);

    SearchAndCalculateVolume<3>(r_background_model_part, r_mpm_model_part);
    GetVolumeVector(r_background_model_part, grid_volume_vector);

    ref_grid_volume_vector = ZeroVector(8);
    ref_grid_volume_vector[4] = 22.0;

    KRATOS_EXPECT_VECTOR_EQ(grid_volume_vector, ref_grid_volume_vector);

    // Move material points and recalculate mp volume at grid
    mp_displacement[0] = 0.5;
    mp_displacement[1] = 0.5;
    mp_displacement[2] = 0.0;
    MoveMaterialPoints(r_mpm_model_part, mp_displacement);

    std::vector<array_1d<double, 3 >> mp_coordinate;
    std::vector<double> mp_volume;
    for (auto& r_element : r_mpm_model_part.Elements())
    {
        r_element.CalculateOnIntegrationPoints(MP_COORD, mp_coordinate, r_mpm_model_part.GetProcessInfo());
        r_element.CalculateOnIntegrationPoints(MP_VOLUME, mp_volume, r_mpm_model_part.GetProcessInfo());
        KRATOS_WATCH(r_element.Id())
        KRATOS_WATCH(mp_coordinate[0])
        KRATOS_WATCH(mp_volume)
    }
    KRATOS_WATCH(r_background_model_part.Elements())


    SearchAndCalculateVolume<3>(r_background_model_part, r_mpm_model_part);
    GetVolumeVector(r_background_model_part, grid_volume_vector);
    ref_grid_volume_vector = ZeroVector(8);
    ref_grid_volume_vector[4] = 4.0;
    ref_grid_volume_vector[5] = 5.0;
    ref_grid_volume_vector[6] = 6.0;
    ref_grid_volume_vector[7] = 7.0;

    KRATOS_EXPECT_VECTOR_EQ(grid_volume_vector, ref_grid_volume_vector);

    // Delete, search material points and recalculate mp volume at grid
    auto& mp_2 = r_mpm_model_part.GetElement(7);
    mp_2.GetGeometry().clear();
    mp_2.Reset(ACTIVE);
    mp_2.Set(TO_ERASE);
    MaterialPointEraseProcess(r_mpm_model_part).Execute();

    SearchAndCalculateVolume<3>(r_background_model_part, r_mpm_model_part);
    GetVolumeVector(r_background_model_part, grid_volume_vector);
    ref_grid_volume_vector = ZeroVector(8);
    ref_grid_volume_vector[4] = 4.0;
    ref_grid_volume_vector[5] = 5.0;
    ref_grid_volume_vector[6] = 6.0;
    ref_grid_volume_vector[7] = 7.0 - 4.5;

    KRATOS_EXPECT_VECTOR_EQ(grid_volume_vector, ref_grid_volume_vector);
}
// TODO:
// KRATOS_TEST_CASE_IN_SUITE(TestVolumeSumUtility3D4N, KratosMPMFastSuite)
// KRATOS_TEST_CASE_IN_SUITE(TestVolumeSumUtility3D8N, KratosMPMFastSuite)



// throw std::runtime_error("\n\n"
// "         |\\__/,|   (`\\\n"
// "       _.|o o  |_   ) )\n"
// "     -(((---(((--------\n"
// "\n"
// "Runtime error: Cat refuses to continue!\n"
// );

// for (auto& r_element : r_mpm_model_part.Elements()) {
//     r_element.CalculateOnIntegrationPoints(MP_COORD,mp_coordinate, r_mpm_model_part.GetProcessInfo());
//     r_element.CalculateOnIntegrationPoints(MP_VOLUME,current_mp_volume, r_mpm_model_part.GetProcessInfo());
//     KRATOS_WATCH(r_element.Id())
//     KRATOS_WATCH(mp_coordinate[0])
//     KRATOS_WATCH(current_mp_volume)
//     KRATOS_WATCH(r_element.GetGeometry().GetGeometryParent(0).Id())

// }
} // namespace Testing
} // namespace Kratos
