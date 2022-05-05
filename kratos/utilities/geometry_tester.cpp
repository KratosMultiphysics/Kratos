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

    if(TestTriangle2D3N(model_part, error_message) == false) successful=false;
    if(TestTriangle2D6N(model_part, error_message) == false) successful=false;
    if(TestQuadrilateral2D4N(model_part, error_message) == false) successful=false;
    if(TestQuadrilateral2D9N(model_part, error_message) == false) successful=false;
    if(TestQuadrilateralInterface2D4N(model_part, error_message) == false) successful=false;

    if(TestTetrahedra3D4N(model_part, error_message) == false) successful=false;
    if(TestTetrahedra3D10N(model_part, error_message) == false) successful=false;
    if(TestHexahedra3D8N(model_part, error_message) == false) successful=false;
    if(TestHexahedra3D20N(model_part, error_message) == false) successful=false;
    if(TestHexahedra3D27N(model_part, error_message) == false) successful=false;
    if(TestHexahedraInterface3D8N(model_part, error_message) == false) successful=false;

    if(TestPrism3D6N(model_part, error_message) == false) successful=false;
    //if(TestPrism3D15N(model_part, error_message) == false) successful=false; 
    if(TestPrismInterface3D6N(model_part, error_message) == false) successful=false;

    KRATOS_WARNING_IF("GeometryTesterUtility", !successful) << "Some errors were detected in the GeometryTester Utility\n" << error_message.str() << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::TestTetrahedra3D4N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    Tetrahedra3D4<Node<3> > geom( rModelPart.pGetNode(4), rModelPart.pGetNode(3), rModelPart.pGetNode(17), rModelPart.pGetNode(19) );

    bool successful = true;

    //this fast function only exists for simplices. Do not use it in other tests
    BoundedMatrix<double, 4,3 > DN_DX;
    array_1d<double, 4 > N;
    double Area;
    GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);

    //compute area by the method area
    const double expected_area = Area;

    if(std::abs(geom.Area() - expected_area) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

    //now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_area, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    rErrorMessage << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::TestTetrahedra3D10N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    Tetrahedra3D10<Node<3> > geom( rModelPart.pGetNode(1), rModelPart.pGetNode(3), rModelPart.pGetNode(7), rModelPart.pGetNode(19),
                                    rModelPart.pGetNode(2), rModelPart.pGetNode(5), rModelPart.pGetNode(4), rModelPart.pGetNode(10),
                                    rModelPart.pGetNode(11), rModelPart.pGetNode(13)
                                );

    bool successful = true;

    //compute area by the method area
    const double area_base = 0.5*std::pow(2.0/3.0,2);
    const double expected_area = area_base*(2.0/3.0) /3.0;

    if(std::abs(geom.Area() - expected_area) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

    //now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_area, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    rErrorMessage << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::TestTriangle2D3N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    Triangle2D3<Node<3> > geom( rModelPart.pGetNode(4), rModelPart.pGetNode(3), rModelPart.pGetNode(8) );

    bool successful = true;

    //this fast function only exists for simplices. Do not use it in other tests
    BoundedMatrix<double, 3, 2 > DN_DX;
    array_1d<double, 3 > N;
    double Area;
    GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);

    //compute area by the method area
    const double expected_area = Area;

    if(std::abs(geom.Area() - expected_area) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

    //now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_area, rErrorMessage) ) successful=false;
//         if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_area, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
//         VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    rErrorMessage << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::TestTriangle2D6N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    Triangle2D6<Node<3> > geom( rModelPart.pGetNode(1), rModelPart.pGetNode(3), rModelPart.pGetNode(7),
                                rModelPart.pGetNode(2), rModelPart.pGetNode(5), rModelPart.pGetNode(4) );

    bool successful = true;

    //compute area by the method area
    const double expected_area = 0.5*2.0/3.0*2.0/3.0;

    if(std::abs(geom.Area() - expected_area) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

    //now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_area, rErrorMessage) ) successful=false;
//         if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_area, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
//         VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    rErrorMessage << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::TestQuadrilateral2D4N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    Quadrilateral2D4<Node<3> > geom( rModelPart.pGetNode(2), rModelPart.pGetNode(6), rModelPart.pGetNode(7), rModelPart.pGetNode(4));

    bool successful = true;

    //compute area by the method area
    const double expected_area = 2.0/3.0*2.0/3.0 - 0.5*1.0/3.0*1.0/3.0 - 0.5* 1.0/3.0*1.0/3.0 - 0.5*2.0/3.0*1.0/3.0;

    if(std::abs(geom.Area() - expected_area) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

    //now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_area, rErrorMessage) ) successful=false;
//         if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_area, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
//         VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    rErrorMessage << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::TestQuadrilateral2D9N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    Quadrilateral2D9<Node<3> > geom( rModelPart.pGetNode(1), rModelPart.pGetNode(3), rModelPart.pGetNode(9), rModelPart.pGetNode(7),
                                        rModelPart.pGetNode(2), rModelPart.pGetNode(6), rModelPart.pGetNode(8), rModelPart.pGetNode(4),
                                        rModelPart.pGetNode(9));

    bool successful = true;

    //compute area by the method area
    const double expected_area = 2.0/3.0*2.0/3.0;

    if(std::abs(geom.Area() - expected_area) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

    //now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_area, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_area, rErrorMessage) ) successful=false;
//         if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_area, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
//         VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    rErrorMessage << std::endl;

    return successful;
}

bool GeometryTesterUtility::TestQuadrilateralInterface2D4N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    QuadrilateralInterface2D4<Node<3> > geom( rModelPart.pGetNode(1), rModelPart.pGetNode(3), rModelPart.pGetNode(6), rModelPart.pGetNode(4));

    bool successful = true;

    //compute area (length in interface geometries)
    const double expected_area = 2.0/3.0;

    if(std::abs(geom.Area() - expected_area) > 1e-14) {
        rErrorMessage << "Geometry Type = " << "Kratos_QuadrilateralInterface3D4" << " --> "
                    << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;
        successful=false;
    }

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::TestHexahedra3D8N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    Hexahedra3D8<Node<3> > geom( rModelPart.pGetNode(2), rModelPart.pGetNode(6), rModelPart.pGetNode(7), rModelPart.pGetNode(4),
                                    rModelPart.pGetNode(11), rModelPart.pGetNode(15), rModelPart.pGetNode(16), rModelPart.pGetNode(13));

    bool successful = true;

    //compute analytical volume
    const double base_area = 2.0/3.0*2.0/3.0 - 0.5*1.0/3.0*1.0/3.0 - 0.5* 1.0/3.0*1.0/3.0 - 0.5*2.0/3.0*1.0/3.0;
    const double expected_vol = base_area*1.0/3.0;

    if(std::abs(geom.Volume() - expected_vol) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

    //now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_vol, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    array_1d<double,3> point_in(3,1.0/3.0);
    array_1d<double,3> point_out(3,5.0);
    if( !VerifyIsInside( geom, point_in, true, rErrorMessage) ) successful=false;
    if( !VerifyIsInside( geom, point_out, false, rErrorMessage) ) successful=false;
    if( !VerfiyShapeFunctionsValues(geom,point_in,rErrorMessage) ) successful = false;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::TestHexahedra3D20N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    Hexahedra3D20<Node<3> > geom( rModelPart.pGetNode(1), rModelPart.pGetNode(3), rModelPart.pGetNode(9), rModelPart.pGetNode(7),
                                    rModelPart.pGetNode(19), rModelPart.pGetNode(21), rModelPart.pGetNode(27), rModelPart.pGetNode(25),
                                    rModelPart.pGetNode(2), rModelPart.pGetNode(6), rModelPart.pGetNode(8), rModelPart.pGetNode(4),
                                    rModelPart.pGetNode(10), rModelPart.pGetNode(12), rModelPart.pGetNode(18), rModelPart.pGetNode(16),
                                    rModelPart.pGetNode(20), rModelPart.pGetNode(24), rModelPart.pGetNode(26), rModelPart.pGetNode(22)
                                );

    bool successful = true;

    //compute analytical volume

    const double expected_vol = std::pow(2.0/3.0,3);

    if(std::abs(geom.Volume() - expected_vol) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

    //now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_vol, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    rErrorMessage << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::TestHexahedra3D27N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    Hexahedra3D27<Node<3> > geom( rModelPart.pGetNode(1), rModelPart.pGetNode(3), rModelPart.pGetNode(9), rModelPart.pGetNode(7),
                                    rModelPart.pGetNode(19), rModelPart.pGetNode(21), rModelPart.pGetNode(27), rModelPart.pGetNode(25),
                                    rModelPart.pGetNode(2), rModelPart.pGetNode(6), rModelPart.pGetNode(8), rModelPart.pGetNode(4),
                                    rModelPart.pGetNode(10), rModelPart.pGetNode(12), rModelPart.pGetNode(18), rModelPart.pGetNode(16),
                                    rModelPart.pGetNode(20), rModelPart.pGetNode(24), rModelPart.pGetNode(26), rModelPart.pGetNode(22),
                                    rModelPart.pGetNode(5), rModelPart.pGetNode(11), rModelPart.pGetNode(15), rModelPart.pGetNode(17),
                                    rModelPart.pGetNode(13), rModelPart.pGetNode(23), rModelPart.pGetNode(14)
                                );

    bool successful = true;

    //compute analytical volume

    const double expected_vol = std::pow(2.0/3.0,3);

    if(std::abs(geom.Volume() - expected_vol) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

    // Now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_vol, rErrorMessage) ) successful=false;
//         if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_1, expected_vol, rErrorMessage) ) successful=false;
//         if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_2, expected_vol, rErrorMessage) ) successful=false;
//         if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_3, expected_vol, rErrorMessage) ) successful=false;
//         if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_4, expected_vol, rErrorMessage) ) successful=false;
//         if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_5, expected_vol, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);
//         VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_1, rErrorMessage);
//         VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_2, rErrorMessage);
//         VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_3, rErrorMessage);
//         VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_4, rErrorMessage);
//         VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_5, rErrorMessage);

    rErrorMessage << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::TestHexahedraInterface3D8N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    HexahedraInterface3D8<Node<3> > geom( rModelPart.pGetNode(1), rModelPart.pGetNode(19), rModelPart.pGetNode(21), rModelPart.pGetNode(3),
                                        rModelPart.pGetNode(4), rModelPart.pGetNode(22), rModelPart.pGetNode(24), rModelPart.pGetNode(6) );

    bool successful = true;

    //compute volume (area in interface geometries)
    const double expected_vol = 2.0/3.0*2.0/3.0;

    if(std::abs(geom.Volume() - expected_vol) > 1e-14) {
        rErrorMessage << "Geometry Type = " << "Kratos_HexahedraInterface3D8" << " --> "
                    << " error: volume returned by the function geom.Volume() does not deliver the correct result " << std::endl;
        successful=false;
    }

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::TestPrism3D6N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    Prism3D6<Node<3> > geom( rModelPart.pGetNode(1), rModelPart.pGetNode(2), rModelPart.pGetNode(4),
                                rModelPart.pGetNode(10),rModelPart.pGetNode(11), rModelPart.pGetNode(13)
                                );

    bool successful = true;

    //compute analytical volume

    const double expected_vol = 1.0/54.0;

    if(std::abs(geom.Volume() - expected_vol) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

    //now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_vol, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    rErrorMessage << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::TestPrism3D15N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    Prism3D15<Node<3> > geom( rModelPart.pGetNode(1),  rModelPart.pGetNode(2),  rModelPart.pGetNode(3),
                            rModelPart.pGetNode(5),  rModelPart.pGetNode(7),  rModelPart.pGetNode(4),
                            rModelPart.pGetNode(10), rModelPart.pGetNode(12), rModelPart.pGetNode(16),
                            rModelPart.pGetNode(19), rModelPart.pGetNode(20), rModelPart.pGetNode(21),
                            rModelPart.pGetNode(23), rModelPart.pGetNode(25), rModelPart.pGetNode(22)
                            );

    bool successful = true;

    // Compute analytical volume
    const double expected_vol = std::pow(2.0/3.0,3)/2.0;

    if(std::abs(geom.Volume() - expected_vol) > 1e-14)
        rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

    //now let's verify that all integration methods give the same
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, expected_vol, rErrorMessage) ) successful=false;
    if( !VerifyAreaByIntegration( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, expected_vol, rErrorMessage) ) successful=false;

    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_1, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_2, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_3, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_4, rErrorMessage);
    VerifyStrainExactness( geom, GeometryData::IntegrationMethod::GI_GAUSS_5, rErrorMessage);

    rErrorMessage << std::endl;

    return successful;
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::TestPrismInterface3D6N(
    ModelPart& rModelPart,
    std::stringstream& rErrorMessage
    )
{
    GenerateNodes(rModelPart);

    PrismInterface3D6<Node<3> > geom( rModelPart.pGetNode(1), rModelPart.pGetNode(19), rModelPart.pGetNode(3),
                                        rModelPart.pGetNode(4), rModelPart.pGetNode(22), rModelPart.pGetNode(6) );

    bool successful = true;

    //compute volume (area in interface geometries)
    const double expected_vol = 0.5*2.0/3.0*2.0/3.0;

    if(std::abs(geom.Volume() - expected_vol) > 1e-14) {
        rErrorMessage << "Geometry Type = " << "Kratos_PrismInterface3D6" << " --> "
                    << " error: volume returned by the function geom.Volume() does not deliver the correct result " << std::endl;
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

bool GeometryTesterUtility::VerifyAreaByIntegration( Geometry<Node<3> >& geom, Geometry<Node<3> >::IntegrationMethod ThisMethod, const double reference_area, std::stringstream& rErrorMessage)
{
    KRATOS_ERROR_IF(geom.WorkingSpaceDimension() != geom.LocalSpaceDimension()) << "VerifyStrainExactness can not be used if LocalSpaceDimension and WorkingSpaceDimension do not coincide --> geometry is " << GetGeometryName(geom) << std::endl;

    double area = 0.0;
    const Element::GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints( ThisMethod );

    if ( integration_points.size() == 0 ) {
        rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " -- the integration method is not supported " << std::endl;
        return false;
    }

    //resizing jacobian inverses containers
    Matrix InvJ0(geom.WorkingSpaceDimension(), geom.WorkingSpaceDimension());

    Element::GeometryType::JacobiansType J0;
    J0= geom.Jacobian( J0, ThisMethod );

    Vector determinants;
    geom.DeterminantOfJacobian(determinants, ThisMethod);

    for ( std::size_t PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ ) {
        const double IntegrationWeight = integration_points[PointNumber].Weight();

        // Calculating and storing inverse of the jacobian and the parameters needed
        const double DetJ0 = MathUtils<double>::Det( J0[PointNumber] );

        if( std::abs(determinants[PointNumber] - DetJ0)/std::abs(DetJ0) > 1e-13) {
            rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " --> " << " determinant as computed from DeterminantOfJacobian does not match the value computed by taking the determinant of J "  << std::endl;
            return true;
        }

        //calculating the total area
        area += DetJ0 * IntegrationWeight;
    }

    if( std::abs(area - reference_area)/reference_area < 1e-13) {
        rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " --> " << " Area Calculation Test: OK "  << std::endl;
        return true;
    } else {
        rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " --> " << " error: the area value " << std::endl;
        rErrorMessage << "                            " << area << " was obtained by integration, while the reference data was "  << reference_area << std::endl;
        return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void GeometryTesterUtility::VerifyStrainExactness( Geometry<Node<3> >& geom,  Geometry<Node<3> >::IntegrationMethod ThisMethod, std::stringstream& rErrorMessage)
{
    const Element::GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints( ThisMethod );
    const std::size_t number_of_nodes = geom.PointsNumber();
    const std::size_t dim = geom.WorkingSpaceDimension();

    KRATOS_ERROR_IF(dim != geom.LocalSpaceDimension()) << "VerifyStrainExactness can not be used if LocalSpaceDimension and WorkingSpaceDimension do not coincide " << GetGeometryName(geom) << std::endl;

    if ( integration_points.size() == 0 ) {
        rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " -- the integration method is not supported " << std::endl;
    } else {
        std::size_t strain_size;
        if(dim == 2) strain_size = 3;
        else strain_size = 6;

        //definition of the expected strain
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

        const Element::GeometryType::ShapeFunctionsGradientsType& DN_De = geom.ShapeFunctionsLocalGradients( ThisMethod );
        const Matrix& Ncontainer = geom.ShapeFunctionsValues( ThisMethod );

        Element::GeometryType::JacobiansType J0;
        geom.Jacobian( J0, ThisMethod );

        //untested functions to be tested
        Element::GeometryType::ShapeFunctionsGradientsType DN_DX_geom;
        geom.ShapeFunctionsIntegrationPointsGradients( DN_DX_geom, ThisMethod );

        Element::GeometryType::JacobiansType Jinv;
        geom.InverseOfJacobian(Jinv, ThisMethod);

        bool successful = true;
        for ( std::size_t PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ ) {
            //check that shape functions sum to 1
            double sum = 0.0;
            for(std::size_t k = 0; k<number_of_nodes; k++) {
                sum += Ncontainer(PointNumber,k);
            }
            if(std::abs(sum-1.0)>1e-14)
                rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " --> " << " error: shape functions do not sum to 1 on gauss point" << std::endl;

            //calculating and storing inverse of the jacobian and the parameters needed
            MathUtils<double>::InvertMatrix( J0[PointNumber], InvJ0, DetJ0 );
            DN_DX  = prod( DN_De[PointNumber], InvJ0 );

            //check that the shape function gradients as obtained from the geomety match what is obtained here starting from the local_gradients
            if(norm_frobenius(DN_DX_geom[PointNumber] - DN_DX)/norm_frobenius(DN_DX) > 1e-13) {
                rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " -->  " << std::endl;
                    rErrorMessage << "     error: shape function gradients are wrongly calculated in function ShapeFunctionsIntegrationPointsGradients: DN_DX_geom " << DN_DX_geom[PointNumber] << " vs " << DN_DX << std::endl;
                rErrorMessage << " norm_frobenius(DN_DX_geom[PointNumber] - DN_DX)/norm_frobenius(DN_DX) = " << norm_frobenius(DN_DX_geom[PointNumber] - DN_DX)/norm_frobenius(DN_DX) <<std::endl;
            }
            if(norm_frobenius(Jinv[PointNumber] - InvJ0)/norm_frobenius(InvJ0) > 1e-13) {
                rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " --> " << std::endl;
                    rErrorMessage << "     error: shape function gradients are wrongly calculated in function ShapeFunctionsIntegrationPointsGradients: DN_DX_geom " << DN_DX_geom[PointNumber] << " vs " << DN_DX << std::endl;
                rErrorMessage << " norm_frobenius(Jinv[PointNumber] - InvJ0)/norm_frobenius(InvJ0) = " << norm_frobenius(Jinv[PointNumber] - InvJ0)/norm_frobenius(InvJ0) <<std::endl;
            }

            CalculateB(B, DN_DX, number_of_nodes, dim);

            //calculate a displacement_field which varies linearly in the space
            for(std::size_t i=0; i<number_of_nodes; i++) {
                const array_1d<double,3>& coords = geom[i].Coordinates();
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
                rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " --> " << " error: expected strain found was not correctly recovered on gauss point. recovered strain = " << strain << " expected value "  << expected_strain << std::endl;
            }
        }

        if(successful == true)
            rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " --> " << " Strain Calculation Test: OK "  << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void GeometryTesterUtility::CalculateB(
    Matrix& B,
    Matrix& DN_DX,
    const std::size_t number_of_nodes,
    const std::size_t dimension
)
{
    KRATOS_TRY

    if ( dimension == 2 )
        B.resize(3, 2*number_of_nodes, false);
    else
        B.resize(6, 3*number_of_nodes, false);

    for ( std::size_t i = 0; i < number_of_nodes; i++ ) {
        std::size_t index = dimension * i;

        if ( dimension == 2 ) {
            B( 0, index + 0 ) = DN_DX( i, 0 );
            B( 0, index + 1 ) = 0.0;
            B( 1, index + 0 ) = 0.0;
            B( 1, index + 1 ) = DN_DX( i, 1 );
            B( 2, index + 0 ) = DN_DX( i, 1 ) ;
            B( 2, index + 1 ) = DN_DX( i, 0 );
        } else {
            B( 0, index + 0 ) = DN_DX( i, 0 );
            B( 0, index + 1 ) = 0.0;
            B( 0, index + 2 ) = 0.0;
            B( 1, index + 0 ) = 0.0;
            B( 1, index + 1 ) = DN_DX( i, 1 );
            B( 1, index + 2 ) = 0.0;
            B( 2, index + 0 ) = 0.0;
            B( 2, index + 1 ) = 0.0;
            B( 2, index + 2 ) = DN_DX( i, 2 );
            B( 3, index + 0 ) = DN_DX( i, 1 );
            B( 3, index + 1 ) = DN_DX( i, 0 );
            B( 3, index + 2 ) = 0.0;
            B( 4, index + 0 ) = 0.0;
            B( 4, index + 1 ) = DN_DX( i, 2 );
            B( 4, index + 2 ) = DN_DX( i, 1 );
            B( 5, index + 0 ) = DN_DX( i, 2 );
            B( 5, index + 1 ) = 0.0;
            B( 5, index + 2 ) = DN_DX( i, 0 );
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::VerifyIsInside(
    Geometry< Node<3> >& geom,
    Geometry< Node<3> >::CoordinatesArrayType& global_coordinates,
    bool expected_result,
    std::stringstream& rErrorMessage)
{
    Geometry< Node<3> >::CoordinatesArrayType local_coordinates;
    if( geom.IsInside(global_coordinates,local_coordinates) == expected_result ) {
        return true;
    } else {
        rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " and point = " << global_coordinates << std::endl;
        rErrorMessage << "Failed VerifyIsInside test. Expected result was: ";
        rErrorMessage << ( (expected_result) ? "inside" : "outside" );
        return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool GeometryTesterUtility::VerfiyShapeFunctionsValues(
    Geometry< Node<3> >& geom,
    Geometry< Node<3> >::CoordinatesArrayType& global_coordinates,
    std::stringstream& rErrorMessage)
{
    Geometry< Node<3> >::CoordinatesArrayType local_coordinates;
    geom.PointLocalCoordinates( local_coordinates, global_coordinates );

    Vector shape_functions = ZeroVector(geom.size());
    geom.ShapeFunctionsValues(shape_functions,local_coordinates);

    array_1d<double,3> residual = global_coordinates;
    for(std::size_t i=0; i<geom.size(); i++) {
        residual -= shape_functions[i]*geom[i].Coordinates();
    }

    if( norm_2(residual) < 1e-15 ) {
    return true;
    } else {
    rErrorMessage << "Geometry Type = " << GetGeometryName(geom) << " and point = " << global_coordinates << std::endl;
    rErrorMessage << "Failed VerfiyShapeFunctionsValues test." << std::endl;
    rErrorMessage << "The difference between exact and interpolated coordinates was : " << residual << std::endl;
    return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

std::string GeometryTesterUtility::GetIntegrationName(Geometry< Node<3> >& geom, Geometry<Node<3> >::IntegrationMethod ThisMethod)
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
    case GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_1 :
        return std::string("GI_EXTENDED_GAUSS_1");
    case GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_2 :
        return std::string("GI_EXTENDED_GAUSS_2");
    case GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_3 :
        return std::string("GI_EXTENDED_GAUSS_3");
    case GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_4 :
        return std::string("GI_EXTENDED_GAUSS_4");
    case GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_5 :
        return std::string("GI_EXTENDED_GAUSS_5");
    case GeometryData::IntegrationMethod::NumberOfIntegrationMethods :
        return std::string("NumberOfIntegrationMethods");
    };

    return std::string("UnknownIntegrationMethod");
}

/***********************************************************************************/
/***********************************************************************************/

std::string GeometryTesterUtility::GetGeometryName(Geometry< Node<3> >& geom)
{
    GeometryData::KratosGeometryType geom_type = geom.GetGeometryType();
    switch(geom_type)
    {
    case GeometryData::KratosGeometryType::Kratos_generic_type :
        return std::string("Kratos_generic_type");
    case GeometryData::KratosGeometryType::Kratos_Hexahedra3D20 :
        return std::string("Kratos_Hexahedra3D20");
    case GeometryData::KratosGeometryType::Kratos_Hexahedra3D27 :
        return std::string("Kratos_Hexahedra3D27");
    case GeometryData::KratosGeometryType::Kratos_Hexahedra3D8 :
        return std::string("Kratos_Hexahedra3D8");
    case GeometryData::KratosGeometryType::Kratos_Prism3D15 :
        return std::string("Kratos_Prism3D15");
    case GeometryData::KratosGeometryType::Kratos_Prism3D6 :
        return std::string("Kratos_Prism3D6");
    case GeometryData::KratosGeometryType::Kratos_Pyramid3D13 :
        return std::string("Kratos_Pyramid3D13");
    case GeometryData::KratosGeometryType::Kratos_Pyramid3D5 :
        return std::string("Kratos_Pyramid3D5");
    case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4 :
        return std::string("Kratos_Quadrilateral2D4");
    case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8 :
        return std::string("Kratos_Quadrilateral2D8");
    case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9 :
        return std::string("Kratos_Quadrilateral2D9");
    case GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4 :
        return std::string("Kratos_Quadrilateral3D4");
    case GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8 :
        return std::string("Kratos_Quadrilateral3D8");
    case GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9 :
        return std::string("Kratos_Quadrilateral3D9");
    case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10 :
        return std::string("Kratos_Tetrahedra3D10");
    case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4 :
        return std::string("Kratos_Tetrahedra3D4");
    case GeometryData::KratosGeometryType::Kratos_Triangle2D3 :
        return std::string("Kratos_Triangle2D3");
    case GeometryData::KratosGeometryType::Kratos_Triangle2D6 :
        return std::string("Kratos_Triangle2D6");
    case GeometryData::KratosGeometryType::Kratos_Triangle3D3 :
        return std::string("Kratos_Triangle3D3");
    case GeometryData::KratosGeometryType::Kratos_Triangle3D6 :
        return std::string("Kratos_Triangle3D6");
    case GeometryData::KratosGeometryType::Kratos_Line2D2 :
        return std::string("Kratos_Line2D2");
    case GeometryData::KratosGeometryType::Kratos_Line2D3 :
        return std::string("Kratos_Line2D3");
    case GeometryData::KratosGeometryType::Kratos_Line3D2 :
        return std::string("Kratos_Line3D2");
    case GeometryData::KratosGeometryType::Kratos_Line3D3 :
        return std::string("Kratos_Line3D3");
    case GeometryData::KratosGeometryType::Kratos_Point2D :
        return std::string("Kratos_Point2D");
    case GeometryData::KratosGeometryType::Kratos_Point3D :
        return std::string("Kratos_Point3D");
    case GeometryData::KratosGeometryType::Kratos_Sphere3D1 :
        return std::string("Kratos_Sphere3D1");
    case GeometryData::KratosGeometryType::Kratos_Nurbs_Curve:
        return std::string("Kratos_Nurbs_Curve");
    case GeometryData::KratosGeometryType::Kratos_Nurbs_Surface:
        return std::string("Kratos_Nurbs_Surface");
    case GeometryData::KratosGeometryType::Kratos_Nurbs_Volume:
        return std::string("Kratos_Nurbs_Volume");
    case GeometryData::KratosGeometryType::Kratos_Nurbs_Curve_On_Surface:
        return std::string("Kratos_Nurbs_Curve_On_Surface");
    case GeometryData::KratosGeometryType::Kratos_Surface_In_Nurbs_Volume:
        return std::string("Kratos_Surface_In_Nurbs_Volume");
    case GeometryData::KratosGeometryType::Kratos_Brep_Curve:
        return std::string("Kratos_Brep_Curve");
    case GeometryData::KratosGeometryType::Kratos_Brep_Surface:
        return std::string("Kratos_Brep_Surface");
    case GeometryData::KratosGeometryType::Kratos_Brep_Curve_On_Surface:
        return std::string("Kratos_Brep_Curve_On_Surface");
    case GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Geometry:
        return std::string("Kratos_Quadrature_Point_Geometry");
    case GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Curve_On_Surface_Geometry:
        return std::string("Kratos_Quadrature_Point_Curve_On_Surface_Geometry");
    case GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Surface_In_Volume_Geometry:
        return std::string("Kratos_Quadrature_Point_Surface_In_Volume_Geometry");
    case GeometryData::KratosGeometryType::Kratos_Coupling_Geometry:
        return std::string("Kratos_Quadrature_Point_Surface_In_Volume_Geometry");
    case GeometryData::KratosGeometryType::NumberOfGeometryTypes:
        return std::string("UnknownGeometry");
    };

    return std::string("UnknownGeometry");
}

} // namespace Kratos.
