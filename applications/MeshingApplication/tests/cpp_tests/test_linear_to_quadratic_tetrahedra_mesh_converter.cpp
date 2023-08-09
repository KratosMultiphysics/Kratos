// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortés
//
//System inlcudes

// Project includes
#include "includes/checks.h"
#include "testing/testing.h"
#include "containers/model.h"
#include "custom_utilities/linear_to_quadratic_tetrahedra_mesh_converter_utility.h"


namespace Kratos {
namespace Testing {

namespace {
    //Type definitions 
    typedef Node NodeType;
    typedef Node::Pointer NodePtrType;
    typedef Geometry<NodeType> GeometryType;
    typedef GeometryType::Pointer GeometryPtrType;
    typedef GeometryType::GeometriesArrayType GeometryArrayType;
    typedef GeometryType::PointsArrayType PointsArrayType;

    //Auxiliar functions
    /**
     * @brief Returns the distance between two 3D points.
     * @param rPoint0 reference to the first point
     * @param rPoint1 reference an array of 3 coordinates representing the second point
     * @return Distance 
     */  
    static double Distance(const NodeType& Point0, const array_1d<double,3>& Point1) {
        const double lx = Point0.X() - Point1[0];
        const double ly = Point0.Y() - Point1[1];
        const double lz = Point0.Z() - Point1[2];

        const double length = lx * lx + ly * ly + lz * lz;

        return std::sqrt( length );
    }
} //unnamed namespace

    KRATOS_TEST_CASE_IN_SUITE(LinearToQuadraticTetrahedraMeshConverter, KratosMeshingApplicationFastSuite)
    {
        Model MyModel;
        ModelPart& modelpart = MyModel.CreateModelPart("Tetrahedras"); 
        ModelPart& subMp0 = modelpart.CreateSubModelPart("level1");
        ModelPart& subMp1 = subMp0.CreateSubModelPart("level2");
        ModelPart& subMp2 = subMp1.CreateSubModelPart("level3");
        modelpart.CreateNewNode(1, 2, 0, 0); 
        modelpart.CreateNewNode(2, 0, 1, 0);
        modelpart.CreateNewNode(3, 0, 0, 1);
        modelpart.CreateNewNode(4, 0, 0, 2);
        modelpart.CreateNewNode(5, -40, 0, 0);
        //subMp to remain unmodified:
        ModelPart& subMp_unmod = modelpart.CreateSubModelPart("unmod");
        modelpart.CreateNewNode(6, 10, 0, 0); 
        modelpart.CreateNewNode(7, 11, 0, 0);
        modelpart.CreateNewNode(8, 10, 1, 0);
        modelpart.CreateNewNode(9, 10, 0, 1);

        Properties::Pointer p_properties_1(new Properties(0)); 
        Element::Pointer tetra1 = subMp1.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties_1);
        Element::Pointer tetra2 = subMp2.CreateNewElement("Element3D4N", 2, {2, 3, 4, 5}, p_properties_1);
        subMp_unmod.CreateNewElement("Element3D4N", 3, {6, 7, 8, 9}, p_properties_1);
        subMp_unmod.CreateNewElement("Element3D3N", 4, {6, 7, 8}, p_properties_1);

        subMp2.AddNode(modelpart.pGetNode(5));
        subMp2.AddNode(modelpart.pGetNode(2));
        subMp2.AddNode(modelpart.pGetNode(3));
        subMp2.AddNode(modelpart.pGetNode(4));
        subMp1.AddNode(modelpart.pGetNode(1));
        subMp_unmod.AddNode(modelpart.pGetNode(6));
        subMp_unmod.AddNode(modelpart.pGetNode(7));
        subMp_unmod.AddNode(modelpart.pGetNode(8));
        subMp_unmod.AddNode(modelpart.pGetNode(9));

    
        std::vector<double> volumes(4); 
        GeometryPtrType geom1 = tetra1->pGetGeometry();
        volumes[tetra1->Id()] = geom1->Volume();
        GeometryPtrType geom2 = tetra2->pGetGeometry();
        volumes[tetra2->Id()] = geom2->Volume(); 

        Condition::Pointer cond1 = subMp2.CreateNewCondition("SurfaceCondition3D3N", 3, {2, 3, 5}, p_properties_1);

        GeometryPtrType geom3 = cond1->pGetGeometry();
        volumes[cond1->Id()] = geom3->Area(); 

        for(auto elem : subMp0.Elements()) {
            elem.SetValue(SPLIT_ELEMENT,true);
        }
        for(auto cond : subMp0.Conditions()) {
            cond.SetValue(SPLIT_ELEMENT,true);
        }

        LinearToQuadraticTetrahedraMeshConverter refineTetra(modelpart); 
        refineTetra.LocalConvertLinearToQuadraticTetrahedraMesh(false,false);

        KRATOS_CHECK_EQUAL(modelpart.Nodes().size(),18); //There are 4 + 14 nodes (10 for each tetra but 6 are shared)
        KRATOS_CHECK_EQUAL(subMp0.Nodes().size(),14); //Also in the first level submodelpart
        KRATOS_CHECK_EQUAL(subMp1.Nodes().size(),14); //Also in the second level submodelpart
        KRATOS_CHECK_EQUAL(subMp2.Nodes().size(),10); //In the third level submodelpart only 10 nodes
        KRATOS_CHECK_EQUAL(modelpart.Elements().size(),4); //No new elements are added
        KRATOS_CHECK_EQUAL(subMp1.Elements().size(),2); //No new elements are added
        KRATOS_CHECK_EQUAL(subMp2.Elements().size(),1); //No new elements are added
        KRATOS_CHECK_EQUAL(subMp2.Conditions().size(),1); //No new conditions are added
        
        for(auto elem : modelpart.Elements()) {
            GeometryPtrType geom = elem.pGetGeometry();
            const auto geometryType = geom->GetGeometryType();
            if(elem.Id()<3){ //replaced elems
                KRATOS_CHECK_EQUAL(geometryType, GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10);
                const double vol = geom->Volume();
                KRATOS_CHECK_NEAR(volumes[elem.Id()], vol ,0.001);

                auto points = geom->Points();
                KRATOS_CHECK_EQUAL(Distance(points[0],points[1]), Distance(points[0],points[4]) + Distance(points[4],points[1]));
                KRATOS_CHECK_EQUAL(Distance(points[1],points[2]), Distance(points[1],points[5]) + Distance(points[5],points[2]));
                KRATOS_CHECK_EQUAL(Distance(points[2],points[3]), Distance(points[2],points[9]) + Distance(points[9],points[3]));
                KRATOS_CHECK_EQUAL(Distance(points[3],points[0]), Distance(points[3],points[7]) + Distance(points[7],points[0]));
                KRATOS_CHECK_EQUAL(Distance(points[0],points[2]), Distance(points[0],points[6]) + Distance(points[6],points[0]));
                KRATOS_CHECK_EQUAL(Distance(points[1],points[3]), Distance(points[1],points[8]) + Distance(points[8],points[1]));
            }else if(elem.Id()==3) { //original linear tet
                KRATOS_CHECK_EQUAL(geometryType, GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4);
            } else {
                KRATOS_CHECK_EQUAL(geometryType, GeometryData::KratosGeometryType::Kratos_Triangle3D3);
            }
        }

        for(auto cond : modelpart.Conditions()) {
            GeometryPtrType geom = cond.pGetGeometry();
            const auto geometryType = geom->GetGeometryType();
            KRATOS_CHECK_EQUAL(geometryType, GeometryData::KratosGeometryType::Kratos_Triangle3D6);
            const double vol = geom->Area();
            KRATOS_CHECK_NEAR(volumes[3], vol ,0.001);

            auto points = geom->Points();
            KRATOS_CHECK_EQUAL(Distance(points[0],points[1]), Distance(points[0],points[3]) + Distance(points[3],points[1]));
            KRATOS_CHECK_EQUAL(Distance(points[1],points[2]), Distance(points[1],points[4]) + Distance(points[4],points[2]));
            KRATOS_CHECK_EQUAL(Distance(points[2],points[0]), Distance(points[2],points[5]) + Distance(points[5],points[0]));

        }

    }

}  // namespace Testing.
}  // namespace Kratos.