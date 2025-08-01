//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes
#include "spaces/ublas_space.h"

// Project includes
#include "utilities/geometry_tester.h"
#include "includes/element.h"

#include "geometries/geometry_data.h"
#include "utilities/geometry_utilities.h"
#include "geometries/geometry.h"

#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"

#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_9.h"
#include "geometries/quadrilateral_interface_2d_4.h"

#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"
#include "geometries/hexahedra_interface_3d_8.h"

#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"
#include "geometries/prism_interface_3d_6.h"

namespace Kratos
{

bool GeometryTesterUtility::RunTest(Model& rModel)
{
    //create a cloud of 27 nodes, to be used in testing the geometries, so that 1 10 19 are on the same vertical
    //side has a lenght 0f 2.0/3.0
    //  25  26  27
    // 22  23  24
    //19--20--21
    //|  16--17--18
    //| 13  14  15
    //10--11--12
    //| 7---8---9
    //|4   5   6
    //1---2---3
    ModelPart& model_part = rModel.CreateModelPart("aux_model_part");
    GenerateNodes(model_part);

    bool successful = true;

    std::stringstream error_message;

    if(StreamTestTriangle2D3N(model_part, error_message) == false) successful=false;
    if(StreamTestTriangle2D6N(model_part, error_message) == false) successful=false;

    if(StreamTestQuadrilateral2D4N(model_part, error_message) == false) successful=false;
    if(StreamTestQuadrilateral2D9N(model_part, error_message) == false) successful=false;
    if(StreamTestQuadrilateralInterface2D4N(model_part, error_message) == false) successful=false;

    if(StreamTestTetrahedra3D4N(model_part, error_message) == false) successful=false;
    if(StreamTestTetrahedra3D10N(model_part, error_message) == false) successful=false;

    if(StreamTestHexahedra3D8N(model_part, error_message) == false) successful=false;
    if(StreamTestHexahedra3D20N(model_part, error_message) == false) successful=false;
    if(StreamTestHexahedra3D27N(model_part, error_message) == false) successful=false;
    if(StreamTestHexahedraInterface3D8N(model_part, error_message) == false) successful=false;

    if(StreamTestPrism3D6N(model_part, error_message) == false) successful=false;
    if(StreamTestPrism3D15N(model_part, error_message) == false) successful=false;
    if(StreamTestPrismInterface3D6N(model_part, error_message) == false) successful=false;

    KRATOS_WARNING_IF("GeometryTesterUtility", !successful) << "Some errors were detected in the GeometryTester Utility\n" << error_message.str() << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::StreamTestTetrahedra3D4N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    Tetrahedra3D4<NodeType> geometry( rModelPart.pGetNode(4), rModelPart.pGetNode(3), rModelPart.pGetNode(17), rModelPart.pGetNode(19) );

    bool successful = true;

    //this fast function only exists for simplices. Do not use it in other tests
    BoundedMatrix<double, 4,3 > DN_DX;
    array_1d<double, 4 > N;
    double Area;
    GeometryUtils::CalculateGeometryData(geometry, DN_DX, N, Area);

    //compute area by the method area
    const double expected_area = Area;

    if(std::abs(geometry.Area() - expected_area) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geometry) << " --> " << " error: area returned by the function geometry.Area() does not deliver the correct result " << std::endl;

    //now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_area, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    array_1d<double,3> point_in(3,1.0/3.0);
    if( !VerifyShapeFunctionsSecondDerivativesValues(geometry,point_in,rErrorMessage) ) successful = false;
    rErrorMessage << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::StreamTestTetrahedra3D10N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    Tetrahedra3D10<NodeType> geometry( rModelPart.pGetNode(1), rModelPart.pGetNode(3), rModelPart.pGetNode(7), rModelPart.pGetNode(19),
                                    rModelPart.pGetNode(2), rModelPart.pGetNode(5), rModelPart.pGetNode(4), rModelPart.pGetNode(10),
                                    rModelPart.pGetNode(11), rModelPart.pGetNode(13)
                                );

    bool successful = true;

    //compute area by the method area
    const double area_base = 0.5*std::pow(2.0/3.0,2);
    const double expected_area = area_base*(2.0/3.0) /3.0;

    if(std::abs(geometry.Area() - expected_area) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geometry) << " --> " << " error: area returned by the function geometry.Area() does not deliver the correct result " << std::endl;

    //now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_area, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    rErrorMessage << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::StreamTestTriangle2D3N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    Triangle2D3<NodeType> geometry( rModelPart.pGetNode(4), rModelPart.pGetNode(3), rModelPart.pGetNode(8) );

    bool successful = true;

    //this fast function only exists for simplices. Do not use it in other tests
    BoundedMatrix<double, 3, 2 > DN_DX;
    array_1d<double, 3 > N;
    double Area;
    GeometryUtils::CalculateGeometryData(geometry, DN_DX, N, Area);

    // Compute area by the method area
    const double expected_area = Area;

    if(std::abs(geometry.Area() - expected_area) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geometry) << " --> " << " error: area returned by the function geometry.Area() does not deliver the correct result " << std::endl;

    //now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_area, rErrorMessage) ) successful=false;
//         if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_area, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
//         VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    array_1d<double,3> point_in(3,1.0/3.0);
    if( !VerifyShapeFunctionsSecondDerivativesValues(geometry,point_in,rErrorMessage) ) successful = false;

    rErrorMessage << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::StreamTestTriangle2D6N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    Triangle2D6<NodeType> geometry( rModelPart.pGetNode(1), rModelPart.pGetNode(3), rModelPart.pGetNode(7),
                                rModelPart.pGetNode(2), rModelPart.pGetNode(5), rModelPart.pGetNode(4) );

    bool successful = true;

    // Compute area by the method area
    const double expected_area = 0.5*2.0/3.0*2.0/3.0;

    if(std::abs(geometry.Area() - expected_area) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geometry) << " --> " << " error: area returned by the function geometry.Area() does not deliver the correct result " << std::endl;

    //now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_area, rErrorMessage) ) successful=false;
//         if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_area, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
//         VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    array_1d<double,3> point_in(3,1.0/3.0);
    if( !VerifyShapeFunctionsSecondDerivativesValues(geometry,point_in,rErrorMessage) ) successful = false;

    rErrorMessage << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::StreamTestQuadrilateral2D4N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    Quadrilateral2D4<NodeType> geometry( rModelPart.pGetNode(2), rModelPart.pGetNode(6), rModelPart.pGetNode(7), rModelPart.pGetNode(4));

    bool successful = true;

    // Compute area by the method area
    const double expected_area = 2.0/3.0*2.0/3.0 - 0.5*1.0/3.0*1.0/3.0 - 0.5* 1.0/3.0*1.0/3.0 - 0.5*2.0/3.0*1.0/3.0;

    if(std::abs(geometry.Area() - expected_area) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geometry) << " --> " << " error: area returned by the function geometry.Area() does not deliver the correct result " << std::endl;

    //now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_area, rErrorMessage) ) successful=false;
//         if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_area, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
//         VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    array_1d<double,3> point_in(3,1.0/3.0);
    if( !VerifyShapeFunctionsSecondDerivativesValues(geometry,point_in,rErrorMessage) ) successful = false;

    rErrorMessage << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::StreamTestQuadrilateral2D9N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    Quadrilateral2D9<NodeType> geometry( rModelPart.pGetNode(1), rModelPart.pGetNode(3), rModelPart.pGetNode(9), rModelPart.pGetNode(7),
                                        rModelPart.pGetNode(2), rModelPart.pGetNode(6), rModelPart.pGetNode(8), rModelPart.pGetNode(4),
                                        rModelPart.pGetNode(9));

    bool successful = true;

    // Compute area by the method area
    const double expected_area = 2.0/3.0*2.0/3.0;

    if(std::abs(geometry.Area() - expected_area) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geometry) << " --> " << " error: area returned by the function geometry.Area() does not deliver the correct result " << std::endl;

    //now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_area, rErrorMessage) ) successful=false;
//         if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_area, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
//         VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    array_1d<double,3> point_in(3,1.0/3.0);
    if( !VerifyShapeFunctionsSecondDerivativesValues(geometry,point_in,rErrorMessage) ) successful = false;
    Quadrilateral2D9<NodeType> exact_geometry( rModelPart.pGetNode(1), rModelPart.pGetNode(3), rModelPart.pGetNode(9), rModelPart.pGetNode(7),
                                        rModelPart.pGetNode(2), rModelPart.pGetNode(6), rModelPart.pGetNode(8), rModelPart.pGetNode(4),
                                        rModelPart.pGetNode(5));
    if( !VerifyShapeFunctionsSecondDerivativesInterpolation(exact_geometry,rErrorMessage)) successful = false;

    rErrorMessage << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::StreamTestQuadrilateralInterface2D4N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    QuadrilateralInterface2D4<NodeType> geometry( rModelPart.pGetNode(1), rModelPart.pGetNode(3), rModelPart.pGetNode(6), rModelPart.pGetNode(4));

    bool successful = true;

    // Compute area (length in interface geometries)
    const double expected_area = 2.0/3.0;

    if(std::abs(geometry.Area() - expected_area) > 1e-14) {
        rErrorMessage << "Geometry Type = " << "Kratos_QuadrilateralInterface3D4" << " --> "
                    << " error: area returned by the function geometry.Area() does not deliver the correct result " << std::endl;
        successful=false;
    }

    array_1d<double,3> point_in(3,1.0/3.0);
    if( !VerifyShapeFunctionsSecondDerivativesValues(geometry,point_in,rErrorMessage) ) successful = false;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::StreamTestHexahedra3D8N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    Hexahedra3D8<NodeType> geometry( rModelPart.pGetNode(2), rModelPart.pGetNode(6), rModelPart.pGetNode(7), rModelPart.pGetNode(4),
                                    rModelPart.pGetNode(11), rModelPart.pGetNode(15), rModelPart.pGetNode(16), rModelPart.pGetNode(13));

    bool successful = true;

    // Compute analytical volume
    const double base_area = 2.0/3.0*2.0/3.0 - 0.5*1.0/3.0*1.0/3.0 - 0.5* 1.0/3.0*1.0/3.0 - 0.5*2.0/3.0*1.0/3.0;
    const double expected_vol = base_area*1.0/3.0;

    if(std::abs(geometry.Volume() - expected_vol) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geometry) << " --> " << " error: area returned by the function geometry.Area() does not deliver the correct result " << std::endl;

    //now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_vol, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    array_1d<double,3> point_in(3,1.0/3.0);
    array_1d<double,3> point_out(3,5.0);
    if( !VerifyIsInside( geometry, point_in, true, rErrorMessage) ) successful=false;
    if( !VerifyIsInside( geometry, point_out, false, rErrorMessage) ) successful=false;
    if( !VerfiyShapeFunctionsValues(geometry,point_in,rErrorMessage) ) successful = false;
    if( !VerifyShapeFunctionsSecondDerivativesValues(geometry,point_in,rErrorMessage) ) successful = false;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::StreamTestHexahedra3D20N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    Hexahedra3D20<NodeType> geometry( rModelPart.pGetNode(1), rModelPart.pGetNode(3), rModelPart.pGetNode(9), rModelPart.pGetNode(7),
                                    rModelPart.pGetNode(19), rModelPart.pGetNode(21), rModelPart.pGetNode(27), rModelPart.pGetNode(25),
                                    rModelPart.pGetNode(2), rModelPart.pGetNode(6), rModelPart.pGetNode(8), rModelPart.pGetNode(4),
                                    rModelPart.pGetNode(10), rModelPart.pGetNode(12), rModelPart.pGetNode(18), rModelPart.pGetNode(16),
                                    rModelPart.pGetNode(20), rModelPart.pGetNode(24), rModelPart.pGetNode(26), rModelPart.pGetNode(22)
                                );

    bool successful = true;

    // Compute analytical volume
    const double expected_vol = std::pow(2.0/3.0,3);

    if(std::abs(geometry.Volume() - expected_vol) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geometry) << " --> " << " error: area returned by the function geometry.Area() does not deliver the correct result " << std::endl;

    //now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_vol, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    rErrorMessage << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::StreamTestHexahedra3D27N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    Hexahedra3D27<NodeType> geometry( rModelPart.pGetNode(1), rModelPart.pGetNode(3), rModelPart.pGetNode(9), rModelPart.pGetNode(7),
                                    rModelPart.pGetNode(19), rModelPart.pGetNode(21), rModelPart.pGetNode(27), rModelPart.pGetNode(25),
                                    rModelPart.pGetNode(2), rModelPart.pGetNode(6), rModelPart.pGetNode(8), rModelPart.pGetNode(4),
                                    rModelPart.pGetNode(10), rModelPart.pGetNode(12), rModelPart.pGetNode(18), rModelPart.pGetNode(16),
                                    rModelPart.pGetNode(20), rModelPart.pGetNode(24), rModelPart.pGetNode(26), rModelPart.pGetNode(22),
                                    rModelPart.pGetNode(5), rModelPart.pGetNode(11), rModelPart.pGetNode(15), rModelPart.pGetNode(17),
                                    rModelPart.pGetNode(13), rModelPart.pGetNode(23), rModelPart.pGetNode(14)
                                );

    bool successful = true;

    // Compute analytical volume
    const double expected_vol = std::pow(2.0/3.0,3);

    if(std::abs(geometry.Volume() - expected_vol) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geometry) << " --> " << " error: area returned by the function geometry.Area() does not deliver the correct result " << std::endl;

    // Now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_vol, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    array_1d<double,3> point_in(3,1.0/3.0);
    if( !VerifyShapeFunctionsSecondDerivativesValues(geometry,point_in,rErrorMessage) ) successful = false;
    if( !VerifyShapeFunctionsSecondDerivativesInterpolation(geometry,rErrorMessage)) successful = false;

    rErrorMessage << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::StreamTestHexahedraInterface3D8N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    HexahedraInterface3D8<NodeType> geometry( rModelPart.pGetNode(1), rModelPart.pGetNode(19), rModelPart.pGetNode(21), rModelPart.pGetNode(3),
                                        rModelPart.pGetNode(4), rModelPart.pGetNode(22), rModelPart.pGetNode(24), rModelPart.pGetNode(6) );

    bool successful = true;

    //compute volume (area in interface geometries)
    const double expected_vol = 2.0/3.0*2.0/3.0;

    if(std::abs(geometry.Volume() - expected_vol) > 1e-14) {
        rErrorMessage << "Geometry Type = " << "Kratos_HexahedraInterface3D8" << " --> "
                    << " error: volume returned by the function geometry.Volume() does not deliver the correct result " << std::endl;
        successful=false;
    }

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::StreamTestPrism3D6N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    Prism3D6<NodeType> geometry( rModelPart.pGetNode(1), rModelPart.pGetNode(2), rModelPart.pGetNode(4),
                                rModelPart.pGetNode(10),rModelPart.pGetNode(11), rModelPart.pGetNode(13)
                                );

    bool successful = true;

    // Compute analytical volume

    const double expected_vol = 1.0/54.0;

    if(std::abs(geometry.Volume() - expected_vol) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geometry) << " --> " << " error: area returned by the function geometry.Area() does not deliver the correct result " << std::endl;

    //now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_vol, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    rErrorMessage << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::StreamTestPrism3D15N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    Prism3D15<NodeType> geometry( rModelPart.pGetNode(1),  rModelPart.pGetNode(3),  rModelPart.pGetNode(7),
                            rModelPart.pGetNode(19),  rModelPart.pGetNode(21),  rModelPart.pGetNode(25),
                            rModelPart.pGetNode(2), rModelPart.pGetNode(5), rModelPart.pGetNode(4),
                            rModelPart.pGetNode(10), rModelPart.pGetNode(12), rModelPart.pGetNode(16),
                            rModelPart.pGetNode(20), rModelPart.pGetNode(23), rModelPart.pGetNode(22)
                            );

    bool successful = true;

    // Compute analytical volume
    const double expected_vol = std::pow(2.0/3.0,3)/2.0;

    if(std::abs(geometry.Volume() - expected_vol) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geometry) << " --> " << " error: area returned by the function geometry.Area() does not deliver the correct result " << std::endl;

    //now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_vol, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
    VerifyStrainExactness( geometry, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    rErrorMessage << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::StreamTestPrismInterface3D6N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    PrismInterface3D6<NodeType> geometry( rModelPart.pGetNode(1), rModelPart.pGetNode(19), rModelPart.pGetNode(3),
                                        rModelPart.pGetNode(4), rModelPart.pGetNode(22), rModelPart.pGetNode(6) );

    bool successful = true;

    //compute volume (area in interface geometries)
    const double expected_vol = 0.5*2.0/3.0*2.0/3.0;

    if(std::abs(geometry.Volume() - expected_vol) > 1e-14) {
        rErrorMessage << "Geometry Type = " << "Kratos_PrismInterface3D6" << " --> "
                    << " error: volume returned by the function geometry.Volume() does not deliver the correct result " << std::endl;
        successful=false;
    }

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

void GeometryTesterUtility::GenerateNodes(ModelPart& rModelPart)
{
    const double dx = 1.0/3.0;
    const double dy = 1.0/3.0;
    const double dz = 1.0/3.0;
    std::size_t counter = 1;
    for(std::size_t k=0; k<3; k++) {
        for(std::size_t j=0; j<3; j++) {
            for(std::size_t i=0; i<3; i++) {
                rModelPart.CreateNewNode(counter++, i*dx, j*dy,k*dz);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::VerifyAreaByIntegration(
    GeometryType& rGeometry,
    GeometryType::IntegrationMethod ThisMethod,
    const double ReferenceArea,
    std::stringstream& rErrorMessage
    )
{
    KRATOS_ERROR_IF(rGeometry.WorkingSpaceDimension() != rGeometry.LocalSpaceDimension()) << "VerifyStrainExactness can not be used if LocalSpaceDimension and WorkingSpaceDimension do not coincide --> geometry is " << GetGeometryName(rGeometry) << std::endl;

    double area = 0.0;
    const Element::GeometryType::IntegrationPointsArrayType& integration_points = rGeometry.IntegrationPoints( ThisMethod );

    if ( integration_points.size() == 0 ) {
        rErrorMessage << "Geometry Type = " << GetGeometryName(rGeometry) << " - IntegrationMethod = " << GetIntegrationName(rGeometry,ThisMethod) << " -- the integration method is not supported " << std::endl;
        return false;
    }

    //resizing jacobian inverses containers
    Matrix InvJ0(rGeometry.WorkingSpaceDimension(), rGeometry.WorkingSpaceDimension());

    Element::GeometryType::JacobiansType J0;
    J0= rGeometry.Jacobian( J0, ThisMethod );

    Vector determinants;
    rGeometry.DeterminantOfJacobian(determinants, ThisMethod);

    for ( std::size_t point_number = 0; point_number < integration_points.size(); point_number++ ) {
        const double integration_weight = integration_points[point_number].Weight();

        // Calculating and storing inverse of the jacobian and the parameters needed
        const double DetJ0 = MathUtils<double>::Det( J0[point_number] );

        if( std::abs(determinants[point_number] - DetJ0)/std::abs(DetJ0) > 1e-13) {
            rErrorMessage << "Geometry Type = " << GetGeometryName(rGeometry) << " - IntegrationMethod = " << GetIntegrationName(rGeometry,ThisMethod) << " --> " << " determinant as computed from DeterminantOfJacobian does not match the value computed by taking the determinant of J "  << std::endl;
            return true;
        }

        // Calculating the total area
        area += DetJ0 * integration_weight;
    }

    if( std::abs(area - ReferenceArea)/ReferenceArea < 1e-13) {
        rErrorMessage << "Geometry Type = " << GetGeometryName(rGeometry) << " - IntegrationMethod = " << GetIntegrationName(rGeometry,ThisMethod) << " --> " << " Area Calculation Test: OK "  << std::endl;
        return true;
    } else {
        rErrorMessage << "Geometry Type = " << GetGeometryName(rGeometry) << " - IntegrationMethod = " << GetIntegrationName(rGeometry,ThisMethod) << " --> " << " error: the area value " << std::endl;
        rErrorMessage << "                            " << area << " was obtained by integration, while the reference data was "  << ReferenceArea << std::endl;
        return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void GeometryTesterUtility::VerifyStrainExactness(
    GeometryType& rGeometry,
    GeometryType::IntegrationMethod ThisMethod,
    std::stringstream& rErrorMessage
    )
{
    const Element::GeometryType::IntegrationPointsArrayType& r_integration_points = rGeometry.IntegrationPoints( ThisMethod );
    const std::size_t number_of_nodes = rGeometry.PointsNumber();
    const std::size_t dim = rGeometry.WorkingSpaceDimension();

    KRATOS_ERROR_IF(dim != rGeometry.LocalSpaceDimension()) << "VerifyStrainExactness can not be used if LocalSpaceDimension and WorkingSpaceDimension do not coincide " << GetGeometryName(rGeometry) << std::endl;

    if ( r_integration_points.size() == 0 ) {
        rErrorMessage << "Geometry Type = " << GetGeometryName(rGeometry) << " - IntegrationMethod = " << GetIntegrationName(rGeometry,ThisMethod) << " -- the integration method is not supported " << std::endl;
    } else {
        std::size_t strain_size;
        if(dim == 2) strain_size = 3;
        else strain_size = 6;

        // Definition of the expected strain
        Matrix MatrixA(dim,dim);
        Vector VectorB(dim);
        for(std::size_t i=0; i<dim; i++) {
            VectorB[i]=i*i+0.567; //arbitrary values
            for(std::size_t j=0; j<dim; j++)
                MatrixA(i,j)=i*j + 0.12345; //initialization fo the values of this matrix is arbitrary
        }


        Vector expected_strain(strain_size);
        if(dim == 2) {
            expected_strain[0] = MatrixA(0,0);
            expected_strain[1] = MatrixA(1,1);
            expected_strain[2] = MatrixA(0,1)+MatrixA(1,0);
        } else {
            expected_strain[0] = MatrixA(0,0);
            expected_strain[1] = MatrixA(1,1);
            expected_strain[2] = MatrixA(2,2);
            expected_strain[3] = MatrixA(0,1)+MatrixA(1,0);
            expected_strain[4] = MatrixA(1,2)+MatrixA(2,1);
            expected_strain[5] = MatrixA(0,2)+MatrixA(2,0);
        }

        // Resizing jacobian inverses containers
        Matrix InvJ0(dim,dim);
        double DetJ0;
        Matrix B;
        Matrix DN_DX;
        Vector displacements(dim*number_of_nodes);

        const Element::GeometryType::ShapeFunctionsGradientsType& DN_De = rGeometry.ShapeFunctionsLocalGradients( ThisMethod );
        const Matrix& Ncontainer = rGeometry.ShapeFunctionsValues( ThisMethod );

        Element::GeometryType::JacobiansType J0;
        rGeometry.Jacobian( J0, ThisMethod );

        //untested functions to be tested
        Element::GeometryType::ShapeFunctionsGradientsType DN_DX_geometry;
        rGeometry.ShapeFunctionsIntegrationPointsGradients( DN_DX_geometry, ThisMethod );

        Element::GeometryType::JacobiansType Jinv;
        rGeometry.InverseOfJacobian(Jinv, ThisMethod);

        bool successful = true;
        for ( std::size_t point_number = 0; point_number < r_integration_points.size(); point_number++ ) {
            //check that shape functions sum to 1
            double sum = 0.0;
            for(std::size_t k = 0; k<number_of_nodes; k++) {
                sum += Ncontainer(point_number,k);
            }
            if(std::abs(sum-1.0)>1e-14)
                rErrorMessage << "Geometry Type = " << GetGeometryName(rGeometry) << " - IntegrationMethod = " << GetIntegrationName(rGeometry,ThisMethod) << " --> " << " error: shape functions do not sum to 1 on gauss point" << std::endl;

            //calculating and storing inverse of the jacobian and the parameters needed
            MathUtils<double>::InvertMatrix( J0[point_number], InvJ0, DetJ0 );
            DN_DX  = prod( DN_De[point_number], InvJ0 );

            //check that the shape function gradients as obtained from the rGeometryety match what is obtained here starting from the local_gradients
            if(norm_frobenius(DN_DX_geometry[point_number] - DN_DX)/norm_frobenius(DN_DX) > 1e-13) {
                rErrorMessage << "Geometry Type = " << GetGeometryName(rGeometry) << " - IntegrationMethod = " << GetIntegrationName(rGeometry,ThisMethod) << " -->  " << std::endl;
                    rErrorMessage << "     error: shape function gradients are wrongly calculated in function ShapeFunctionsIntegrationPointsGradients: DN_DX_geometry " << DN_DX_geometry[point_number] << " vs " << DN_DX << std::endl;
                rErrorMessage << " norm_frobenius(DN_DX_geometry[point_number] - DN_DX)/norm_frobenius(DN_DX) = " << norm_frobenius(DN_DX_geometry[point_number] - DN_DX)/norm_frobenius(DN_DX) <<std::endl;
            }
            if(norm_frobenius(Jinv[point_number] - InvJ0)/norm_frobenius(InvJ0) > 1e-13) {
                rErrorMessage << "Geometry Type = " << GetGeometryName(rGeometry) << " - IntegrationMethod = " << GetIntegrationName(rGeometry,ThisMethod) << " --> " << std::endl;
                    rErrorMessage << "     error: shape function gradients are wrongly calculated in function ShapeFunctionsIntegrationPointsGradients: DN_DX_geometry " << DN_DX_geometry[point_number] << " vs " << DN_DX << std::endl;
                rErrorMessage << " norm_frobenius(Jinv[point_number] - InvJ0)/norm_frobenius(InvJ0) = " << norm_frobenius(Jinv[point_number] - InvJ0)/norm_frobenius(InvJ0) <<std::endl;
            }

            CalculateB(B, DN_DX, number_of_nodes, dim);

            //calculate a displacement_field which varies linearly in the space
            for(std::size_t i=0; i<number_of_nodes; i++) {
                const array_1d<double,3>& coords = rGeometry[i].Coordinates();
                Vector disp(dim);
                for(std::size_t k=0; k<dim; k++) {
                    disp[k] = VectorB[k];
                    for(std::size_t l=0; l<dim; l++) {
                        disp[k] += MatrixA(k,l)*coords[l] ;
                    }
                }
//                     Vector disp = prod(MatrixA,) + VectorB;
                for(std::size_t k=0; k<dim; k++) {
                    displacements[i*dim+k] = disp[k];
                }
            }

            Vector strain = prod(B,displacements);

            Vector strain_err = strain-expected_strain;

            if( norm_2(strain_err)/norm_2(expected_strain) < 1e-14) {
                //do nothing
            } else {
                successful = false;
                rErrorMessage << "Geometry Type = " << GetGeometryName(rGeometry) << " - IntegrationMethod = " << GetIntegrationName(rGeometry,ThisMethod) << " --> " << " error: expected strain found was not correctly recovered on gauss point. recovered strain = " << strain << " expected value "  << expected_strain << std::endl;
            }
        }

        if(successful == true)
            rErrorMessage << "Geometry Type = " << GetGeometryName(rGeometry) << " - IntegrationMethod = " << GetIntegrationName(rGeometry,ThisMethod) << " --> " << " Strain Calculation Test: OK "  << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void GeometryTesterUtility::CalculateB(
    Matrix& rB,
    Matrix& rDN_DX,
    const std::size_t NumberOfNodes,
    const std::size_t Dimension
)
{
    KRATOS_TRY

    if ( Dimension == 2 )
        rB.resize(3, 2*NumberOfNodes, false);
    else
        rB.resize(6, 3*NumberOfNodes, false);

    for ( std::size_t i = 0; i < NumberOfNodes; i++ ) {
        std::size_t index = Dimension * i;

        if ( Dimension == 2 ) {
            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 0, index + 1 ) = 0.0;
            rB( 1, index + 0 ) = 0.0;
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rDN_DX( i, 1 ) ;
            rB( 2, index + 1 ) = rDN_DX( i, 0 );
        } else {
            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 0, index + 1 ) = 0.0;
            rB( 0, index + 2 ) = 0.0;
            rB( 1, index + 0 ) = 0.0;
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 1, index + 2 ) = 0.0;
            rB( 2, index + 0 ) = 0.0;
            rB( 2, index + 1 ) = 0.0;
            rB( 2, index + 2 ) = rDN_DX( i, 2 );
            rB( 3, index + 0 ) = rDN_DX( i, 1 );
            rB( 3, index + 1 ) = rDN_DX( i, 0 );
            rB( 3, index + 2 ) = 0.0;
            rB( 4, index + 0 ) = 0.0;
            rB( 4, index + 1 ) = rDN_DX( i, 2 );
            rB( 4, index + 2 ) = rDN_DX( i, 1 );
            rB( 5, index + 0 ) = rDN_DX( i, 2 );
            rB( 5, index + 1 ) = 0.0;
            rB( 5, index + 2 ) = rDN_DX( i, 0 );
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::VerifyIsInside(
    GeometryType& rGeometry,
    GeometryType::CoordinatesArrayType& rGlobalCoordinates,
    bool ExpectedResult,
    std::stringstream& rErrorMessage)
{
    GeometryType::CoordinatesArrayType local_coordinates;
    if( rGeometry.IsInside(rGlobalCoordinates,local_coordinates) == ExpectedResult ) {
        return true;
    } else {
        rErrorMessage << "Geometry Type = " << GetGeometryName(rGeometry) << " and point = " << rGlobalCoordinates << std::endl;
        rErrorMessage << "Failed VerifyIsInside test. Expected result was: ";
        rErrorMessage << ( (ExpectedResult) ? "inside" : "outside" );
        return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::VerfiyShapeFunctionsValues(
    GeometryType& rGeometry,
    GeometryType::CoordinatesArrayType& rGlobalCoordinates,
    std::stringstream& rErrorMessage
    )
{
    GeometryType::CoordinatesArrayType local_coordinates;
    rGeometry.PointLocalCoordinates( local_coordinates, rGlobalCoordinates );

    Vector shape_functions = ZeroVector(rGeometry.size());
    rGeometry.ShapeFunctionsValues(shape_functions,local_coordinates);

    array_1d<double,3> residual = rGlobalCoordinates;
    for(std::size_t i=0; i<rGeometry.size(); i++) {
        residual -= shape_functions[i]*rGeometry[i].Coordinates();
    }

    if( norm_2(residual) < 1e-15 ) {
        return true;
    } else {
        rErrorMessage << "Geometry Type = " << GetGeometryName(rGeometry) << " and point = " << rGlobalCoordinates << std::endl;
        rErrorMessage << "Failed VerfiyShapeFunctionsValues test." << std::endl;
        rErrorMessage << "The difference between exact and interpolated coordinates was : " << residual << std::endl;
        return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::VerifyShapeFunctionsSecondDerivativesValues(
    GeometryType& rGeometry,
    GeometryType::CoordinatesArrayType& rGlobalCoordinates,
    std::stringstream& rErrorMessage
    ) const
{
    GeometryType::CoordinatesArrayType local_coordinates;
    rGeometry.PointLocalCoordinates( local_coordinates, rGlobalCoordinates );

    Vector f_1 = ZeroVector(rGeometry.size());
    Vector f_2 = ZeroVector(rGeometry.size());
    Vector f_3 = ZeroVector(rGeometry.size());
    Vector f_4 = ZeroVector(rGeometry.size());

    DenseVector<Matrix> DDN_DX, H;
    Vector ei,ej,f;

    double delta = 1e-1;
    unsigned int dim = rGeometry.WorkingSpaceDimension();
    rGeometry.ShapeFunctionsSecondDerivatives(DDN_DX,local_coordinates);

    if ( H.size() != rGeometry.size() )
    {
        // KLUDGE: While there is a bug in
        // ublas vector resize, I have to put this beside resizing!!
        DenseVector<Matrix> temp( rGeometry.size() );
        H.swap( temp );
    }


    for ( unsigned int i = 0; i < rGeometry.size(); i++ ) H[i].resize(dim, dim, false);

    for (unsigned int i = 0; i<dim;i++){
        for (unsigned int j = 0; j<dim;j++){
            ei = ZeroVector(3);
            ei[i] = 1.0;
            ej = ZeroVector(3);
            ej[j] = 1.0;
            std::transform(ei.begin(), ei.end(), ei.begin(), [delta](double c){ return c*delta; });
            std::transform(ej.begin(), ej.end(), ej.begin(), [delta](double c){ return c*delta; });
            rGeometry.ShapeFunctionsValues(f_1,local_coordinates + ei + ej);
            rGeometry.ShapeFunctionsValues(f_2,local_coordinates + ei - ej);
            rGeometry.ShapeFunctionsValues(f_3,local_coordinates - ei + ej);
            rGeometry.ShapeFunctionsValues(f_4,local_coordinates - ei - ej);
            f = f_1-f_2-f_3+f_4;
            std::transform(f.begin(), f.end(), f.begin(), [delta](double c){ return c/(4.0*std::pow(delta,2)); });
            for (unsigned int k = 0; k<rGeometry.size();k++){
                H[k](i,j) = f[k];
            }
        }
    }

    for (unsigned int i = 0; i<rGeometry.size();i++){
        if(norm_frobenius(DDN_DX[i] - H[i]) > 1e-13) {
            rErrorMessage << "     error: shape function second derivatives are wrongly calculated in function ShapeFunctionsSecondDerivatives: DDN_DX[point_number] " << DDN_DX[i] << " vs " << H[i] << std::endl;
            rErrorMessage << " norm_frobenius(DDN_DX[point_number] - H[point_number])/norm_frobenius(H[point_number]) = " << norm_frobenius(DDN_DX[i] - H[i])/norm_frobenius(H[i]) <<std::endl;
            return false;
        }
    }

    return true;
}

bool GeometryTesterUtility::VerifyShapeFunctionsSecondDerivativesInterpolation(
    GeometryType& rGeometry,
    std::stringstream& rErrorMessage
    ) const
{
    DenseVector<DenseVector<Matrix>> DDN_DDX;
    DenseVector<Matrix> DN_DX;

    const unsigned int Dim = rGeometry.WorkingSpaceDimension();
    const unsigned int NumNodes = rGeometry.PointsNumber();

    const Geometry<Node>::IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_4;

    Matrix exact_hess = ZeroMatrix(Dim, Dim);

    GeometryUtils::ShapeFunctionsSecondDerivativesTransformOnAllIntegrationPoints( DDN_DDX, rGeometry, integration_method );

    Matrix NContainer = rGeometry.ShapeFunctionsValues(integration_method);

    const Geometry<Node>::IntegrationPointsArrayType integration_points = rGeometry.IntegrationPoints(integration_method);
    const unsigned int number_of_integration_points = integration_points.size();

    Matrix gauss_point_coordinates = ZeroMatrix(number_of_integration_points,Dim);

    for(unsigned int g = 0; g < number_of_integration_points; g++){
            Matrix hess = ZeroMatrix(Dim, Dim);
            for (unsigned int i = 0; i < NumNodes; ++i){
                array_1d<double, 3>& r_coordinates = rGeometry[i].Coordinates();
                const double x1 = r_coordinates[0];
                const double x2 = r_coordinates[1];
                double x3 = 0.0;
                if (Dim == 3) x3 = r_coordinates[2];

                const auto r_function = std::pow(x1,2)*std::pow(x2,2) + std::pow(x1,2)*std::pow(x3,2) + std::pow(x2,2)*std::pow(x3,2) + std::pow(x1,2)*x2 + std::pow(x1,2)*x3 + std::pow(x2,2)*x1 + std::pow(x2,2)*x3 + std::pow(x3,2)*x1 + std::pow(x3,2)*x2 + std::pow(x1,2) + std::pow(x2,2) + std::pow(x3,2) + x1*x2*x3 + x1*x2 + x2*x3 + x1*x3 + x1 + x2 + x3 + 1.0;

                for (unsigned int d = 0; d < Dim; ++d){
                    gauss_point_coordinates(g,d) += NContainer(g,i) * r_coordinates[d];
                    for (unsigned int e = 0; e < Dim; ++e){
                        hess(d,e) += DDN_DDX[g][i](d,e) * r_function;
                        }
                    }
                }

            const double x1 = gauss_point_coordinates(g,0);
            const double x2 = gauss_point_coordinates(g,1);
            double x3 = 0.0;
            if (Dim == 3) x3 = gauss_point_coordinates(g,2);

            if (Dim == 2){
                exact_hess(0,0) = 2.0 * std::pow(x2,2) + 2.0 * std::pow(x3,2) + 2.0 * x2 + 2.0 * x3 + 2.0;
                exact_hess(0,1) = 4.0 * x1 * x2 + 2.0 * x1 + 2.0 * x2 + x3 + 1.0;
                exact_hess(1,0) = 4.0 * x1 * x2 + 2.0 * x1 + 2.0 * x2 + x3 + 1.0;
                exact_hess(1,1) = 2.0 * std::pow(x1,2) + 2.0 * std::pow(x3,2) + 2.0 * x1 + 2.0 * x3 + 2.0;
            }
            else if (Dim == 3){
                exact_hess(0,0) = 2.0 * std::pow(x2,2) + 2.0 * std::pow(x3,2) + 2.0 * x2 + 2.0 * x3 + 2.0;
                exact_hess(0,1) = 4.0 * x1 * x2 + 2.0 * x1 + 2.0 * x2 + x3 + 1.0;
                exact_hess(0,2) = 4.0 * x1 * x3 + 2.0 * x1 + 2.0 * x3 + x2 + 1.0;
                exact_hess(1,0) = 4.0 * x1 * x2 + 2.0 * x1 + 2.0 * x2 + x3 + 1.0;
                exact_hess(1,1) = 2.0 * std::pow(x1,2) + 2.0 * std::pow(x3,2) + 2.0 * x1 + 2.0 * x3 + 2.0;
                exact_hess(1,2) = 4.0 * x2 * x3 + 2.0 * x2 + 2.0 * x3 + x1 + 1.0;
                exact_hess(2,0) = 4.0 * x1 * x3 + 2.0 * x1 + 2.0 * x3 + x2 + 1.0;
                exact_hess(2,1) = 4.0 * x3 * x2 + 2.0 * x3 + 2.0 * x2 + x1 + 1.0;
                exact_hess(2,2) = 2.0 * std::pow(x1,2) + 2.0 * std::pow(x2,2) + 2.0 * x1 + 2.0 * x2 + 2.0;
            }

            for (unsigned int d = 0; d < Dim; ++d){
                for (unsigned int e = 0; e < Dim; ++e){
                    if ((hess(d,e) - exact_hess(d,e)) >= 1e-10){
                        rErrorMessage << " error: shape function second derivatives are wrongly calculated in function ShapeFunctionsIntegrationPointsSecondDerivatives: hess(d,e) " << hess(d,e) << " vs " << exact_hess(d,e) << std::endl;
                        return false;
                    }
                }
            }

        }

    return true;
}

/***********************************************************************************/
/***********************************************************************************/

std::string GeometryTesterUtility::GetIntegrationName(
    GeometryType& rGeometry,
    GeometryType::IntegrationMethod ThisMethod
    )
{
    switch(ThisMethod)
    {
    case GeometryData::IntegrationMethod::GI_GAUSS_1 :
        return std::string("GI_GAUSS_1");
    case GeometryData::IntegrationMethod::GI_GAUSS_2 :
        return std::string("GI_GAUSS_2");
    case GeometryData::IntegrationMethod::GI_GAUSS_3 :
        return std::string("GI_GAUSS_3");
    case GeometryData::IntegrationMethod::GI_GAUSS_4 :
        return std::string("GI_GAUSS_4");
    case GeometryData::IntegrationMethod::GI_GAUSS_5 :
        return std::string("GI_GAUSS_5");
    case GeometryData::IntegrationMethod::GI_LOBATTO_1 :
        return std::string("GI_LOBATTO_1");
    case GeometryData::IntegrationMethod::NumberOfIntegrationMethods :
        return std::string("NumberOfIntegrationMethods");
    };

    return std::string("UnknownIntegrationMethod");
}

/***********************************************************************************/
/***********************************************************************************/

std::string GeometryTesterUtility::GetGeometryName(GeometryType& rGeometry)
{
    return GeometryUtils::GetGeometryName(rGeometry.GetGeometryType());
}

} // namespace Kratos.
