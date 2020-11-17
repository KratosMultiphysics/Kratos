//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/pointer_vector.h"
#include "geometries/nurbs_volume_geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_volume_grid_modeler.h"

#include "tests/cpp_tests/geometries/test_geometry.h"
#include "geometries/point.h"
#include "containers/model.h"

namespace Kratos {
namespace Testing {

typedef Node<3> NodeType;
typedef NurbsVolumeGeometry<PointerVector<NodeType>> NurbsVolumeGeometryType;
typedef typename NurbsVolumeGeometry<PointerVector<NodeType>>::Pointer NurbsVolumeGeometryPointerType;

KRATOS_TEST_CASE_IN_SUITE(NurbsVolumeGridModelerCreateGrid, KratosCoreNurbsGeometriesFastSuite) {
    Point A(0.0, 0.0, 0.0);
    Point B(1.0, 1.0, 1.0);
    SizeType order_u = 8;
    SizeType order_v = 1;
    SizeType order_w = 1;
    SizeType num_elements_u = 5;
    SizeType num_elements_v = 1;
    SizeType num_elements_w = 1;
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Grid",2);
    auto geom = NurbsVolumeGridModeler::CreateGrid(model_part, A, B, order_u, order_v, order_w, num_elements_u, num_elements_v, num_elements_w);
    for( int i = 0; i < geom->size(); ++i){
        std::cout << "CP: " << (*geom)[i] << std::endl;
    }
}

} // End namespace testing
} // End namespace kratos