//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortés
//
//

// Project includes
#include "custom_utilities/tet10_refinement_utility.h"


namespace Kratos {
namespace Testing {

    //Type definitions 
    typedef Node<3> NodeType;
    typedef Node<3>::Pointer NodePtrType;
    typedef Geometry<NodeType> GeometryType;
    typedef GeometryType::Pointer GeometryPtrType;
    typedef GeometryType::GeometriesArrayType GeometryArrayType;
    typedef GeometryType::PointsArrayType PointsArrayType;

    KRATOS_TEST_CASE_IN_SUITE(Tet10RefinementUtility, KratosMeshingApplicationFastSuite)
    {
        Model MyModel;
        ModelPart& modelpart = MyModel.CreateModelPart("Tetrahedras"); 
        modelpart.CreateNewNode(1, 1, 0, 0); 
        modelpart.CreateNewNode(2, 0, 1, 0);
        modelpart.CreateNewNode(3, 0, 0, 1);
        modelpart.CreateNewNode(4, 0, 0, 2);
        modelpart.CreateNewNode(5, 0, 2, 0);
        Properties::Pointer p_properties_1(new Properties(0)); 
        Element::Pointer tetra1 = modelpart.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties_1);
        GeometryPtrType ptetra1 = tetra1->pGetGeometry();
        Element::Pointer tetra2 = modelpart.CreateNewElement("Element3D4N", 2, {2, 3, 4, 5}, p_properties_1);
        GeometryPtrType ptetra2 = tetra2->pGetGeometry();
        
        //tetra1->SetValue(SPLIT_ELEMENT,true);
        //tetra2->SetValue(SPLIT_ELEMENT,true);

        Tet10RefinementUtility refineTetra(modelpart); 
        refineTetra.LocalRefineTet10Mesh(true);

        KRATOS_CHECK_EQUAL(modelpart.Nodes().size(),14); //There are 14 nodes (10 for each tetra but 6 are shared) 
        KRATOS_CHECK_EQUAL(modelpart.Elements().size(),2); //No new elements are added
        
        for(auto elem : modelpart.Elements()) {
            KRATOS_CHECK_EQUAL(typeid(elem.GetGeometry()), typeid(Tetrahedra3D10<NodeType>));
            
        }

    }

}  // namespace Testing.
}  // namespace Kratos.