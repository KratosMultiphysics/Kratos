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

    bool is_skew = false;
    std::string element;
    IndexType number_of_nodes = 16;
    switch (caseIndex)
    {
    case 0:
        element = "Element2D4N";
        break;
    case 1:
        is_skew = true;
        element = "Element2D4N";
        break;
    case 10:
        element = "Element2D3N";
        break;
    case 11:
        is_skew = true;
        element = "Element2D3N";
        break;
    case 20:
        element = "Element3D8N";
        number_of_nodes *= 2;
        break;
    case 21:
        is_skew = true;
        element = "Element3D8N";
        number_of_nodes *= 2;
        break;
    default:
        break;
    }


    // Nodes
    std::vector<Node::Pointer> point_vector(number_of_nodes);
    const double dx = 1.0;
    const double dy = 1.0;
    IndexType point_index = 1;
    for (size_t row = 0; row < 4; ++row) {
        for (size_t col = 0; col < 4; ++col) {
            double myskew = (col > 0 && row == 1 && is_skew) ? 0.1 : 0.0;
            point_vector[point_index-1] = rBackgroundModelPart.CreateNewNode(point_index,
                double(col)*dx + myskew, double(row)*dy, 0.0);
            point_index += 1;
        }
    }

    if (element == "Element2D4N")
    {
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
        IndexType point_index = number_of_nodes/2 + 1;
        for (size_t row = 0; row < 4; ++row) {
            for (size_t col = 0; col < 4; ++col) {
                double myskew = (col > 0 && row == 1 && is_skew) ? 0.1 : 0.0;
                point_vector[point_index - 1] = rBackgroundModelPart.CreateNewNode(point_index,
                    double(col) * dx + myskew, double(row) * dy, 1.0);
                point_index += 1;
            }
        }
        rBackgroundModelPart.CreateNewElement(element, 1, { 1, 2, 6, 5, 17,18,22,21}, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 2, { 2, 3, 7, 6 , 18,19,23,22}, nullptr);

        rBackgroundModelPart.CreateNewElement(element, 4, { 5, 6, 10, 9 , 21,22,26,25}, nullptr);
        rBackgroundModelPart.CreateNewElement(element, 5, { 6, 7, 11, 10 ,22,23,27,26}, nullptr);
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



    IndexType grid_id = 0;
    for (auto& r_element : r_background_model_part.Elements()) {
        r_element.GetGeometry().SetId(grid_id);
        grid_id++;
    }

    std::vector<array_1d<double, 3 >> mp_coordinate;
    std::vector<double> current_mp_volume;
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
