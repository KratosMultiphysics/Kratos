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
//

#if !defined(KRATOS_GEOMETRY_TEST_H_INCLUDED )
#define  KRATOS_GEOMETRY_TEST_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>

// External includes


// Project includes
#include "includes/define.h"
#include "containers/model.h"
#include "includes/element.h"
#include "geometries/geometry_data.h"
#include "utilities/geometry_utilities.h"
#include "geometries/geometry.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_9.h"

#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"

#include "geometries/prism_3d_6.h"
//#include "geometries/prism_3d_15.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class GeometryTesterUtility
 * @ingroup KratosCore
 * @brief This utility tests the geometries
 * @todo Migrate to the new cpp test interface
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
*/
class GeometryTesterUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometryTesterUtility
    KRATOS_CLASS_POINTER_DEFINITION(GeometryTesterUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    GeometryTesterUtility()
    {
    }

    /// Default constructor.
    bool RunTest(Model& rModel) //std::string& out_error_msg)
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

        bool succesful = true;

        if(TestTriangle2D3N(model_part) == false) succesful=false;
        if(TestTriangle2D6N(model_part) == false) succesful=false;
        if(TestQuadrilateral2D4N(model_part) == false) succesful=false;
        if(TestQuadrilateral2D9N(model_part) == false) succesful=false;

        if(TestTetrahedra3D4N(model_part) == false) succesful=false;
        if(TestTetrahedra3D10N(model_part) == false) succesful=false;
        if(TestHexahedra3D8N(model_part) == false) succesful=false;
        if(TestHexahedra3D20N(model_part) == false) succesful=false;
        if(TestHexahedra3D27N(model_part) == false) succesful=false;

        if(TestPrism3D6N(model_part) == false) succesful=false;
//        if(TestPrism3D15N( error_msg ) == false) succesful=false;

        if(succesful == false)
            std::cout << "*** some errors were detected in the GeometryTester Utility ***" << std::endl;

        return succesful;
    }

    /// Destructor.
    virtual ~GeometryTesterUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    bool TestTetrahedra3D4N(ModelPart& rModelPart)
    {
        GenerateNodes(rModelPart);

        std::stringstream error_msg;
        Tetrahedra3D4<Node<3> > geom( rModelPart.pGetNode(4), rModelPart.pGetNode(3), rModelPart.pGetNode(17), rModelPart.pGetNode(19) );

        bool succesful = true;

        //this fast function only exists for simplices. Do not use it in other tests
        BoundedMatrix<double, 4,3 > DN_DX;
        array_1d<double, 4 > N;
        double Area;
        GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);

        //compute area by the method area
        const double expected_area = Area;

        if(std::abs(geom.Area() - expected_area) > 1e-14)
            error_msg << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

        //now let's verify that all integration methods give the same
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_1, expected_area, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_2, expected_area, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_3, expected_area, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_4, expected_area, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_5, expected_area, error_msg) ) succesful=false;

        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_1, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_2, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_3, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_4, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_5, error_msg);

        if(succesful == false)
            KRATOS_THROW_ERROR(std::logic_error,"", error_msg.str());

        error_msg << std::endl;

        return succesful;

    }

    bool TestTetrahedra3D10N(ModelPart& rModelPart)
    {
        GenerateNodes(rModelPart);

        std::stringstream error_msg;

        Tetrahedra3D10<Node<3> > geom( rModelPart.pGetNode(1), rModelPart.pGetNode(3), rModelPart.pGetNode(7), rModelPart.pGetNode(19),
                                        rModelPart.pGetNode(2), rModelPart.pGetNode(5), rModelPart.pGetNode(4), rModelPart.pGetNode(10),
                                        rModelPart.pGetNode(11), rModelPart.pGetNode(13)
                                    );

        bool succesful = true;

        //compute area by the method area
        const double area_base = 0.5*pow(2.0/3.0,2);
        const double expected_area = area_base*(2.0/3.0) /3.0;

        if(std::abs(geom.Area() - expected_area) > 1e-14)
            error_msg << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

        //now let's verify that all integration methods give the same
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_1, expected_area, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_2, expected_area, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_3, expected_area, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_4, expected_area, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_5, expected_area, error_msg) ) succesful=false;

        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_1, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_2, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_3, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_4, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_5, error_msg);

        error_msg << std::endl;

        if(succesful == false)
            KRATOS_THROW_ERROR(std::logic_error,"", error_msg.str());

        return succesful;

    }

    bool TestTriangle2D3N(ModelPart& rModelPart)
    {
        GenerateNodes(rModelPart);

        std::stringstream error_msg;
        Triangle2D3<Node<3> > geom( rModelPart.pGetNode(4), rModelPart.pGetNode(3), rModelPart.pGetNode(8) );

        bool succesful = true;

        //this fast function only exists for simplices. Do not use it in other tests
        BoundedMatrix<double, 3, 2 > DN_DX;
        array_1d<double, 3 > N;
        double Area;
        GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);

        //compute area by the method area
        const double expected_area = Area;

        if(std::abs(geom.Area() - expected_area) > 1e-14)
            error_msg << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

        //now let's verify that all integration methods give the same
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_1, expected_area, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_2, expected_area, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_3, expected_area, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_4, expected_area, error_msg) ) succesful=false;
//         if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_5, expected_area, error_msg) ) succesful=false;

        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_1, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_2, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_3, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_4, error_msg);
//         VerifyStrainExactness( geom, GeometryData::GI_GAUSS_5, error_msg);

        error_msg << std::endl;

        if(succesful == false)
            KRATOS_THROW_ERROR(std::logic_error,"", error_msg.str());

        return succesful;

    }

    bool TestTriangle2D6N(ModelPart& rModelPart)
    {
        GenerateNodes(rModelPart);

        std::stringstream error_msg;
        Triangle2D6<Node<3> > geom( rModelPart.pGetNode(1), rModelPart.pGetNode(3), rModelPart.pGetNode(7),
                                    rModelPart.pGetNode(2), rModelPart.pGetNode(5), rModelPart.pGetNode(4) );

        bool succesful = true;

        //compute area by the method area
        const double expected_area = 0.5*2.0/3.0*2.0/3.0;

        if(std::abs(geom.Area() - expected_area) > 1e-14)
            error_msg << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

        //now let's verify that all integration methods give the same
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_1, expected_area, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_2, expected_area, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_3, expected_area, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_4, expected_area, error_msg) ) succesful=false;
//         if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_5, expected_area, error_msg) ) succesful=false;

        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_1, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_2, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_3, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_4, error_msg);
//         VerifyStrainExactness( geom, GeometryData::GI_GAUSS_5, error_msg);

        error_msg << std::endl;

        if(succesful == false)
            KRATOS_THROW_ERROR(std::logic_error,"", error_msg.str());

        return succesful;

    }

    bool TestQuadrilateral2D4N(ModelPart& rModelPart)
    {
        GenerateNodes(rModelPart);

        std::stringstream error_msg;
        Quadrilateral2D4<Node<3> > geom( rModelPart.pGetNode(2), rModelPart.pGetNode(6), rModelPart.pGetNode(7), rModelPart.pGetNode(4));

        bool succesful = true;

        //compute area by the method area
        const double expected_area = 2.0/3.0*2.0/3.0 - 0.5*1.0/3.0*1.0/3.0 - 0.5* 1.0/3.0*1.0/3.0 - 0.5*2.0/3.0*1.0/3.0;

        if(std::abs(geom.Area() - expected_area) > 1e-14)
            error_msg << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

        //now let's verify that all integration methods give the same
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_1, expected_area, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_2, expected_area, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_3, expected_area, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_4, expected_area, error_msg) ) succesful=false;
//         if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_5, expected_area, error_msg) ) succesful=false;

        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_1, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_2, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_3, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_4, error_msg);
//         VerifyStrainExactness( geom, GeometryData::GI_GAUSS_5, error_msg);

        error_msg << std::endl;

        if(succesful == false)
            KRATOS_THROW_ERROR(std::logic_error,"", error_msg.str());

        return succesful;

    }

    bool TestQuadrilateral2D9N(ModelPart& rModelPart)
    {
        GenerateNodes(rModelPart);

        std::stringstream error_msg;
        Quadrilateral2D9<Node<3> > geom( rModelPart.pGetNode(1), rModelPart.pGetNode(3), rModelPart.pGetNode(9), rModelPart.pGetNode(7),
                                         rModelPart.pGetNode(2), rModelPart.pGetNode(6), rModelPart.pGetNode(8), rModelPart.pGetNode(4),
                                         rModelPart.pGetNode(9));

        bool succesful = true;

        //compute area by the method area
        const double expected_area = 2.0/3.0*2.0/3.0;

        if(std::abs(geom.Area() - expected_area) > 1e-14)
            error_msg << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

        //now let's verify that all integration methods give the same
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_1, expected_area, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_2, expected_area, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_3, expected_area, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_4, expected_area, error_msg) ) succesful=false;
//         if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_5, expected_area, error_msg) ) succesful=false;

        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_1, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_2, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_3, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_4, error_msg);
//         VerifyStrainExactness( geom, GeometryData::GI_GAUSS_5, error_msg);

        error_msg << std::endl;

        if(succesful == false)
            KRATOS_THROW_ERROR(std::logic_error,"", error_msg.str());

        return succesful;

    }

    bool TestHexahedra3D8N(ModelPart& rModelPart)
    {
        GenerateNodes(rModelPart);

        std::stringstream error_msg;
        Hexahedra3D8<Node<3> > geom( rModelPart.pGetNode(2), rModelPart.pGetNode(6), rModelPart.pGetNode(7), rModelPart.pGetNode(4),
                                     rModelPart.pGetNode(11), rModelPart.pGetNode(15), rModelPart.pGetNode(16), rModelPart.pGetNode(13));

        bool succesful = true;

        //compute analytical volume
        const double base_area = 2.0/3.0*2.0/3.0 - 0.5*1.0/3.0*1.0/3.0 - 0.5* 1.0/3.0*1.0/3.0 - 0.5*2.0/3.0*1.0/3.0;
        const double expected_vol = base_area*1.0/3.0;

        if(std::abs(geom.Volume() - expected_vol) > 1e-14)
            error_msg << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

        //now let's verify that all integration methods give the same
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_1, expected_vol, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_2, expected_vol, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_3, expected_vol, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_4, expected_vol, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_5, expected_vol, error_msg) ) succesful=false;

        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_1, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_2, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_3, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_4, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_5, error_msg);

        array_1d<double,3> point_in(3,1.0/3.0);
        array_1d<double,3> point_out(3,5.0);
        if( !VerifyIsInside( geom, point_in, true, error_msg) ) succesful=false;
        if( !VerifyIsInside( geom, point_out, false, error_msg) ) succesful=false;
        if( !VerfiyShapeFunctionsValues(geom,point_in,error_msg) ) succesful = false;

        if(succesful == false)
            KRATOS_THROW_ERROR(std::logic_error,"", error_msg.str());

        return succesful;

    }

    bool TestHexahedra3D20N(ModelPart& rModelPart)
    {
        GenerateNodes(rModelPart);

        std::stringstream error_msg;
        Hexahedra3D20<Node<3> > geom( rModelPart.pGetNode(1), rModelPart.pGetNode(3), rModelPart.pGetNode(9), rModelPart.pGetNode(7),
                                      rModelPart.pGetNode(19), rModelPart.pGetNode(21), rModelPart.pGetNode(27), rModelPart.pGetNode(25),
                                      rModelPart.pGetNode(2), rModelPart.pGetNode(6), rModelPart.pGetNode(8), rModelPart.pGetNode(4),
                                      rModelPart.pGetNode(10), rModelPart.pGetNode(12), rModelPart.pGetNode(18), rModelPart.pGetNode(16),
                                      rModelPart.pGetNode(20), rModelPart.pGetNode(24), rModelPart.pGetNode(26), rModelPart.pGetNode(22)
                                    );

        bool succesful = true;

        //compute analytical volume

        const double expected_vol = pow(2.0/3.0,3);

        if(std::abs(geom.Volume() - expected_vol) > 1e-14)
            error_msg << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

        //now let's verify that all integration methods give the same
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_1, expected_vol, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_2, expected_vol, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_3, expected_vol, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_4, expected_vol, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_5, expected_vol, error_msg) ) succesful=false;

        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_1, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_2, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_3, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_4, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_5, error_msg);

        error_msg << std::endl;

        if(succesful == false)
            KRATOS_THROW_ERROR(std::logic_error,"", error_msg.str());

        return succesful;

    }


    bool TestHexahedra3D27N(ModelPart& rModelPart)
    {
        GenerateNodes(rModelPart);

        std::stringstream error_msg;
        Hexahedra3D27<Node<3> > geom( rModelPart.pGetNode(1), rModelPart.pGetNode(3), rModelPart.pGetNode(9), rModelPart.pGetNode(7),
                                      rModelPart.pGetNode(19), rModelPart.pGetNode(21), rModelPart.pGetNode(27), rModelPart.pGetNode(25),
                                      rModelPart.pGetNode(2), rModelPart.pGetNode(6), rModelPart.pGetNode(8), rModelPart.pGetNode(4),
                                      rModelPart.pGetNode(10), rModelPart.pGetNode(12), rModelPart.pGetNode(18), rModelPart.pGetNode(16),
                                      rModelPart.pGetNode(20), rModelPart.pGetNode(24), rModelPart.pGetNode(26), rModelPart.pGetNode(22),
                                      rModelPart.pGetNode(5), rModelPart.pGetNode(11), rModelPart.pGetNode(15), rModelPart.pGetNode(17),
                                      rModelPart.pGetNode(13), rModelPart.pGetNode(23), rModelPart.pGetNode(14)
                                    );

        bool succesful = true;

        //compute analytical volume

        const double expected_vol = pow(2.0/3.0,3);

        if(std::abs(geom.Volume() - expected_vol) > 1e-14)
            error_msg << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

        // Now let's verify that all integration methods give the same
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_1, expected_vol, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_2, expected_vol, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_3, expected_vol, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_4, expected_vol, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_5, expected_vol, error_msg) ) succesful=false;
//         if( !VerifyAreaByIntegration( geom, GeometryData::GI_EXTENDED_GAUSS_1, expected_vol, error_msg) ) succesful=false;
//         if( !VerifyAreaByIntegration( geom, GeometryData::GI_EXTENDED_GAUSS_2, expected_vol, error_msg) ) succesful=false;
//         if( !VerifyAreaByIntegration( geom, GeometryData::GI_EXTENDED_GAUSS_3, expected_vol, error_msg) ) succesful=false;
//         if( !VerifyAreaByIntegration( geom, GeometryData::GI_EXTENDED_GAUSS_4, expected_vol, error_msg) ) succesful=false;
//         if( !VerifyAreaByIntegration( geom, GeometryData::GI_EXTENDED_GAUSS_5, expected_vol, error_msg) ) succesful=false;

        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_1, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_2, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_3, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_4, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_5, error_msg);
//         VerifyStrainExactness( geom, GeometryData::GI_EXTENDED_GAUSS_1, error_msg);
//         VerifyStrainExactness( geom, GeometryData::GI_EXTENDED_GAUSS_2, error_msg);
//         VerifyStrainExactness( geom, GeometryData::GI_EXTENDED_GAUSS_3, error_msg);
//         VerifyStrainExactness( geom, GeometryData::GI_EXTENDED_GAUSS_4, error_msg);
//         VerifyStrainExactness( geom, GeometryData::GI_EXTENDED_GAUSS_5, error_msg);

        error_msg << std::endl;

        if(succesful == false)
        {
            KRATOS_THROW_ERROR(std::logic_error,"", error_msg.str());
        }

        return succesful;

    }

    bool TestPrism3D6N(ModelPart& rModelPart)
    {
        GenerateNodes(rModelPart);

        std::stringstream error_msg;
        Prism3D6<Node<3> > geom( rModelPart.pGetNode(1), rModelPart.pGetNode(2), rModelPart.pGetNode(4),
                                 rModelPart.pGetNode(10),rModelPart.pGetNode(11), rModelPart.pGetNode(13)
                                    );

        bool succesful = true;

        //compute analytical volume

        const double expected_vol = 1.0/54.0;

        if(std::abs(geom.Volume() - expected_vol) > 1e-14)
            error_msg << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

        //now let's verify that all integration methods give the same
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_1, expected_vol, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_2, expected_vol, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_3, expected_vol, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_4, expected_vol, error_msg) ) succesful=false;
        if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_5, expected_vol, error_msg) ) succesful=false;

        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_1, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_2, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_3, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_4, error_msg);
        VerifyStrainExactness( geom, GeometryData::GI_GAUSS_5, error_msg);

        error_msg << std::endl;

        if(succesful == false)
            KRATOS_THROW_ERROR(std::logic_error,"", error_msg.str());

        return succesful;

    }

//    bool TestPrism3D15N(ModelPart& rModelPart)
//     {
//         std::stringstream error_msg;
//          Prism3D15<Node<3> > geom( rModelPart.pGetNode(1),  rModelPart.pGetNode(2),  rModelPart.pGetNode(3),
//                                    rModelPart.pGetNode(5),  rModelPart.pGetNode(7),  rModelPart.pGetNode(4),
//                                    rModelPart.pGetNode(10), rModelPart.pGetNode(12), rModelPart.pGetNode(16),
//                                    rModelPart.pGetNode(19), rModelPart.pGetNode(20), rModelPart.pGetNode(21),
//                                    rModelPart.pGetNode(23), rModelPart.pGetNode(25), rModelPart.pGetNode(22)
//                                    );

//          bool succesful = true;

//          //compute analytical volume

//          const double expected_vol = pow(2.0/3.0,3)/2.0;

//          if(std::abs(geom.Volume() - expected_vol) > 1e-14)
//              error_msg << "Geometry Type = " << GetGeometryName(geom) << " --> " << " error: area returned by the function geom.Area() does not deliver the correct result " << std::endl;

//          //now let's verify that all integration methods give the same
//          if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_1, expected_vol, error_msg) ) succesful=false;
//          if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_2, expected_vol, error_msg) ) succesful=false;
//          if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_3, expected_vol, error_msg) ) succesful=false;
//          if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_4, expected_vol, error_msg) ) succesful=false;
//          if( !VerifyAreaByIntegration( geom, GeometryData::GI_GAUSS_5, expected_vol, error_msg) ) succesful=false;

//          VerifyStrainExactness( geom, GeometryData::GI_GAUSS_1, error_msg);
//          VerifyStrainExactness( geom, GeometryData::GI_GAUSS_2, error_msg);
//          VerifyStrainExactness( geom, GeometryData::GI_GAUSS_3, error_msg);
//          VerifyStrainExactness( geom, GeometryData::GI_GAUSS_4, error_msg);
//          VerifyStrainExactness( geom, GeometryData::GI_GAUSS_5, error_msg);

//          error_msg << std::endl;

//          return succesful;
//    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "GeometryTesterUtility" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "GeometryTesterUtility";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


    ///@}
    ///@name Friends
    ///@{
    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    void GenerateNodes(ModelPart& rModelPart)
    {
        const double dx = 0.333333333333333333333;
        const double dy = 0.333333333333333333333;
        const double dz = 0.333333333333333333333;
        std::size_t counter = 1;
        for(unsigned int k=0; k<3; k++)
        {
            for(unsigned int j=0; j<3; j++)
            {
                for(unsigned int i=0; i<3; i++)
                {
                    rModelPart.CreateNewNode(counter++, i*dx, j*dy,k*dz);
                }
            }
        }
    }

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{
    bool VerifyAreaByIntegration( Geometry<Node<3> >& geom, Geometry<Node<3> >::IntegrationMethod ThisMethod, const double reference_area, std::stringstream& error_msg)
    {
        if(geom.WorkingSpaceDimension() != geom.LocalSpaceDimension())
            KRATOS_THROW_ERROR(std::logic_error,"VerifyStrainExactness can not be used if LocalSpaceDimension and WorkingSpaceDimension do not coincide --> geometry is ",GetGeometryName(geom) );

        double area = 0.0;
        const Element::GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints( ThisMethod );

        if ( integration_points.size() == 0 )
        {
            error_msg << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " -- the integration method is not supported " << std::endl;
            return false;
        }

        //resizing jacobian inverses containers
        Matrix InvJ0(geom.WorkingSpaceDimension(), geom.WorkingSpaceDimension());

        Element::GeometryType::JacobiansType J0;
        J0= geom.Jacobian( J0, ThisMethod );

        Vector determinants;
        geom.DeterminantOfJacobian(determinants, ThisMethod);

        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            const double IntegrationWeight = integration_points[PointNumber].Weight();

            //calculating and storing inverse of the jacobian and the parameters needed
            double DetJ0 = MathUtils<double>::Det( J0[PointNumber] );

            if( std::abs(determinants[PointNumber] - DetJ0)/std::abs(DetJ0) > 1e-13)
            {
                error_msg << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " --> " << " determinant as computed from DeterminantOfJacobian does not match the value computed by taking the determinant of J "  << std::endl;
                return true;
            }


            //calculating the total area
            area += DetJ0 * IntegrationWeight;
        }

        if( std::abs(area - reference_area)/reference_area < 1e-13)
        {
            error_msg << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " --> " << " Area Calculation Test: OK "  << std::endl;
            return true;
        }
        else
        {
            error_msg << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " --> " << " error: the area value " << std::endl;
            error_msg << "                            " << area << " was obtained by integration, while the reference data was "  << reference_area << std::endl;
            return false;
        }

    }


    //here we verify that a  "displacement field" which varies linearly in space, produces the expected strain distribution.
    //this shall be considered a test for shape function derivatives
    void VerifyStrainExactness( Geometry<Node<3> >& geom,  Geometry<Node<3> >::IntegrationMethod ThisMethod, std::stringstream& error_msg)
    {

        const Element::GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints( ThisMethod );
        const unsigned int number_of_nodes = geom.PointsNumber();
        const unsigned int dim = geom.WorkingSpaceDimension();

        if(dim != geom.LocalSpaceDimension())
            KRATOS_THROW_ERROR(std::logic_error,"VerifyStrainExactness can not be used if LocalSpaceDimension and WorkingSpaceDimension do not coincide ",GetGeometryName(geom) );

        if ( integration_points.size() == 0 )
        {
            error_msg << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " -- the integration method is not supported " << std::endl;
        }
        else
        {


            unsigned int strain_size;
            if(dim == 2) strain_size = 3;
            else strain_size = 6;

            //definition of the expected strain
            Matrix MatrixA(dim,dim);
            Vector VectorB(dim);
            for(unsigned int i=0; i<dim; i++)
            {
                VectorB[i]=i*i+0.567; //arbitrary values
                for(unsigned int j=0; j<dim; j++)
                    MatrixA(i,j)=i*j + 0.12345; //initialization fo the values of this matrix is arbitrary
            }


            Vector expected_strain(strain_size);
            if(dim == 2)
            {
                expected_strain[0] = MatrixA(0,0);
                expected_strain[1] = MatrixA(1,1);
                expected_strain[2] = MatrixA(0,1)+MatrixA(1,0);
            }
            else
            {
                expected_strain[0] = MatrixA(0,0);
                expected_strain[1] = MatrixA(1,1);
                expected_strain[2] = MatrixA(2,2);
                expected_strain[3] = MatrixA(0,1)+MatrixA(1,0);
                expected_strain[4] = MatrixA(1,2)+MatrixA(2,1);
                expected_strain[5] = MatrixA(0,2)+MatrixA(2,0);
            }

            //resizing jacobian inverses containers
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

            bool succesful = true;
            for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
            {
                //check that shape functions sum to 1
                double sum = 0.0;
                for(unsigned int k = 0; k<number_of_nodes; k++)
                {
                    sum += Ncontainer(PointNumber,k);
                }
                if(std::abs(sum-1.0)>1e-14)
                    error_msg << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " --> " << " error: shape functions do not sum to 1 on gauss point" << std::endl;

                //calculating and storing inverse of the jacobian and the parameters needed
                MathUtils<double>::InvertMatrix( J0[PointNumber], InvJ0, DetJ0 );
                DN_DX  = prod( DN_De[PointNumber], InvJ0 );

                //check that the shape function gradients as obtained from the geomety match what is obtained here starting from the local_gradients
                if(norm_frobenius(DN_DX_geom[PointNumber] - DN_DX)/norm_frobenius(DN_DX) > 1e-13)
                {
                    error_msg << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " -->  " << std::endl;
                     error_msg << "     error: shape function gradients are wrongly calculated in function ShapeFunctionsIntegrationPointsGradients: DN_DX_geom " << DN_DX_geom[PointNumber] << " vs " << DN_DX << std::endl;
                    error_msg << " norm_frobenius(DN_DX_geom[PointNumber] - DN_DX)/norm_frobenius(DN_DX) = " << norm_frobenius(DN_DX_geom[PointNumber] - DN_DX)/norm_frobenius(DN_DX) <<std::endl;
                }
                if(norm_frobenius(Jinv[PointNumber] - InvJ0)/norm_frobenius(InvJ0) > 1e-13)
                {
                    error_msg << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " --> " << std::endl;
                     error_msg << "     error: shape function gradients are wrongly calculated in function ShapeFunctionsIntegrationPointsGradients: DN_DX_geom " << DN_DX_geom[PointNumber] << " vs " << DN_DX << std::endl;
                    error_msg << " norm_frobenius(Jinv[PointNumber] - InvJ0)/norm_frobenius(InvJ0) = " << norm_frobenius(Jinv[PointNumber] - InvJ0)/norm_frobenius(InvJ0) <<std::endl;
                }

                CalculateB(B, DN_DX, number_of_nodes, dim);

                //calculate a displacement_field which varies linearly in the space
                for(unsigned int i=0; i<number_of_nodes; i++)
                {
                    const array_1d<double,3>& coords = geom[i].Coordinates();
                    Vector disp(dim);
                    for(unsigned int k=0; k<dim; k++)
                    {
                        disp[k] = VectorB[k];
                        for(unsigned int l=0; l<dim; l++)
                        {
                            disp[k] += MatrixA(k,l)*coords[l] ;
                        }
                    }
//                     Vector disp = prod(MatrixA,) + VectorB;
                    for(unsigned int k=0; k<dim; k++)
                    {
                        displacements[i*dim+k] = disp[k];
                    }
                }

                Vector strain = prod(B,displacements);

                Vector strain_err = strain-expected_strain;

                if( norm_2(strain_err)/norm_2(expected_strain) < 1e-14)
                {
                    //do nothing
                }
                else
                {
                    succesful = false;
                    error_msg << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " --> " << " error: expected strain found was not correctly recovered on gauss point. recovered strain = " << strain << " expected value "  << expected_strain << std::endl;
                }



            }

            if(succesful == true)
                error_msg << "Geometry Type = " << GetGeometryName(geom) << " - IntegrationMethod = " << GetIntegrationName(geom,ThisMethod) << " --> " << " Strain Calculation Test: OK "  << std::endl;

        }
    }




    //this function computes the linear strain matrix - useful to verify that a constant strain can be correctly reproduced
    void CalculateB(
        Matrix& B,
        Matrix& DN_DX,
        const unsigned int number_of_nodes,
        const unsigned int dimension
    )
    {
        KRATOS_TRY

        if ( dimension == 2 )
            B.resize(3, 2*number_of_nodes, false);
        else
            B.resize(6, 3*number_of_nodes, false);

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = dimension * i;

            if ( dimension == 2 )
            {
                B( 0, index + 0 ) = DN_DX( i, 0 );
                B( 0, index + 1 ) = 0.0;
                B( 1, index + 0 ) = 0.0;
                B( 1, index + 1 ) = DN_DX( i, 1 );
                B( 2, index + 0 ) = DN_DX( i, 1 ) ;
                B( 2, index + 1 ) = DN_DX( i, 0 );
            }
            else
            {
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

    bool VerifyIsInside(
        Geometry< Node<3> >& geom,
        Geometry< Node<3> >::CoordinatesArrayType& global_coordinates,
        bool expected_result,
        std::stringstream& error_msg)
    {
      Geometry< Node<3> >::CoordinatesArrayType local_coordinates;
      if( geom.IsInside(global_coordinates,local_coordinates) == expected_result )
        return true;
      else
      {
        error_msg << "Geometry Type = " << GetGeometryName(geom) << " and point = " << global_coordinates << std::endl;
        error_msg << "Failed VerifyIsInside test. Expected result was: ";
        error_msg << ( (expected_result) ? "inside" : "outside" );
        return false;
      }
    }

    bool VerfiyShapeFunctionsValues(
        Geometry< Node<3> >& geom,
        Geometry< Node<3> >::CoordinatesArrayType& global_coordinates,
        std::stringstream& error_msg)
    {
      Geometry< Node<3> >::CoordinatesArrayType local_coordinates;
      geom.PointLocalCoordinates( local_coordinates, global_coordinates );

      Vector shape_functions = ZeroVector(geom.size());
      geom.ShapeFunctionsValues(shape_functions,local_coordinates);

      array_1d<double,3> residual = global_coordinates;
      for(unsigned int i=0; i<geom.size(); i++)
      {
           residual -= shape_functions[i]*geom[i].Coordinates();
      }

      if( norm_2(residual) < 1e-15 )
      {
        return true;
      }
      else
      {
        error_msg << "Geometry Type = " << GetGeometryName(geom) << " and point = " << global_coordinates << std::endl;
        error_msg << "Failed VerfiyShapeFunctionsValues test." << std::endl;
        error_msg << "The difference between exact and interpolated coordinates was : " << residual << std::endl;
        return false;
      }
    }

    std::string GetIntegrationName(Geometry< Node<3> >& geom, Geometry<Node<3> >::IntegrationMethod ThisMethod)
    {
        switch(ThisMethod)
        {
        case GeometryData::GI_GAUSS_1 :
            return std::string("GI_GAUSS_1");
        case GeometryData::GI_GAUSS_2 :
            return std::string("GI_GAUSS_2");
        case GeometryData::GI_GAUSS_3 :
            return std::string("GI_GAUSS_3");
        case GeometryData::GI_GAUSS_4 :
            return std::string("GI_GAUSS_4");
        case GeometryData::GI_GAUSS_5 :
            return std::string("GI_GAUSS_5");
        case GeometryData::GI_EXTENDED_GAUSS_1 :
            return std::string("GI_EXTENDED_GAUSS_1");
        case GeometryData::GI_EXTENDED_GAUSS_2 :
            return std::string("GI_EXTENDED_GAUSS_2");
        case GeometryData::GI_EXTENDED_GAUSS_3 :
            return std::string("GI_EXTENDED_GAUSS_3");
        case GeometryData::GI_EXTENDED_GAUSS_4 :
            return std::string("GI_EXTENDED_GAUSS_4");
        case GeometryData::GI_EXTENDED_GAUSS_5 :
            return std::string("GI_EXTENDED_GAUSS_5");
        case GeometryData::NumberOfIntegrationMethods :
            return std::string("NumberOfIntegrationMethods");
        };

        return std::string("UnknownIntegrationMethod");
    }
    std::string GetGeometryName(Geometry< Node<3> >& geom)
    {
        GeometryData::KratosGeometryType geom_type = geom.GetGeometryType();
        switch(geom_type)
        {
        case GeometryData::Kratos_generic_type :
            return std::string("Kratos_generic_type");
        case GeometryData::Kratos_Hexahedra3D20 :
            return std::string("Kratos_Hexahedra3D20");
        case GeometryData::Kratos_Hexahedra3D27 :
            return std::string("Kratos_Hexahedra3D27");
        case GeometryData::Kratos_Hexahedra3D8 :
            return std::string("Kratos_Hexahedra3D8");
        case GeometryData::Kratos_Prism3D15 :
            return std::string("Kratos_Prism3D15");
        case GeometryData::Kratos_Prism3D6 :
            return std::string("Kratos_Prism3D6");
        case GeometryData::Kratos_Quadrilateral2D4 :
            return std::string("Kratos_Quadrilateral2D4");
        case GeometryData::Kratos_Quadrilateral2D8 :
            return std::string("Kratos_Quadrilateral2D8");
        case GeometryData::Kratos_Quadrilateral2D9 :
            return std::string("Kratos_Quadrilateral2D9");
        case GeometryData::Kratos_Quadrilateral3D4 :
            return std::string("Kratos_Quadrilateral3D4");
        case GeometryData::Kratos_Quadrilateral3D8 :
            return std::string("Kratos_Quadrilateral3D8");
        case GeometryData::Kratos_Quadrilateral3D9 :
            return std::string("Kratos_Quadrilateral3D9");
        case GeometryData::Kratos_Tetrahedra3D10 :
            return std::string("Kratos_Tetrahedra3D10");
        case GeometryData::Kratos_Tetrahedra3D4 :
            return std::string("Kratos_Tetrahedra3D4");
        case GeometryData::Kratos_Triangle2D3 :
            return std::string("Kratos_Triangle2D3");
        case GeometryData::Kratos_Triangle2D6 :
            return std::string("Kratos_Triangle2D6");
        case GeometryData::Kratos_Triangle3D3 :
            return std::string("Kratos_Triangle3D3");
        case GeometryData::Kratos_Triangle3D6 :
            return std::string("Kratos_Triangle3D6");
        case GeometryData::Kratos_Line2D2 :
            return std::string("Kratos_Line2D2");
        case GeometryData::Kratos_Line2D3 :
            return std::string("Kratos_Line2D3");
        case GeometryData::Kratos_Line3D2 :
            return std::string("Kratos_Line3D2");
        case GeometryData::Kratos_Line3D3 :
            return std::string("Kratos_Line3D3");
        case GeometryData::Kratos_Point2D :
            return std::string("Kratos_Point2D");
        case GeometryData::Kratos_Point3D :
            return std::string("Kratos_Point3D");
        case GeometryData::Kratos_Sphere3D1 :
            return std::string("Kratos_Sphere3D1");
        };

        return std::string("UnknownGeometry");
    }


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    GeometryTesterUtility& operator=(GeometryTesterUtility const& rOther) {return *this;}

    /// Copy constructor.
    GeometryTesterUtility(GeometryTesterUtility const& rOther) {}


    ///@}

}; // Class GeometryTesterUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  GeometryTesterUtility& rThis) {
  return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GeometryTesterUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_TEST_H_INCLUDED  defined
