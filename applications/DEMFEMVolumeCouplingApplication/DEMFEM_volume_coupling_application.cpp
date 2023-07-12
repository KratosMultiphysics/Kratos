/*
 * Author: Miguel Angel Celigueta
 *
 *  maceli@cimne.upc.edu
 */

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "DEMFEM_volume_coupling_application.h"

#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/pyramid_3d_5.h"
#include "geometries/pyramid_3d_13.h"
#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"

#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/line_2d_4.h"
#include "geometries/line_2d_5.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/triangle_2d_10.h"
#include "geometries/triangle_2d_15.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"

namespace Kratos {

      mVolumeCouplingElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mVolumeCouplingElement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mVolumeCouplingElement2D6N(0, Element::GeometryType::Pointer(new Triangle2D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mVolumeCouplingElement2D8N(0, Element::GeometryType::Pointer(new Quadrilateral2D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mVolumeCouplingElement2D9N(0, Element::GeometryType::Pointer(new Quadrilateral2D9<NodeType >(Element::GeometryType::PointsArrayType(9)))),
      mVolumeCouplingElement2D10N(0, Element::GeometryType::Pointer(new Triangle2D10<NodeType >(Element::GeometryType::PointsArrayType(10)))),
      mVolumeCouplingElement2D15N(0, Element::GeometryType::Pointer(new Triangle2D15<NodeType >(Element::GeometryType::PointsArrayType(15)))),
      mVolumeCouplingElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mVolumeCouplingElement3D5N(0, Element::GeometryType::Pointer(new Pyramid3D5<NodeType >(Element::GeometryType::PointsArrayType(5)))),
      mVolumeCouplingElement3D6N(0, Element::GeometryType::Pointer(new Prism3D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mVolumeCouplingElement3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mVolumeCouplingElement3D10N(0, Element::GeometryType::Pointer(new Tetrahedra3D10<NodeType >(Element::GeometryType::PointsArrayType(10)))),
      mVolumeCouplingElement3D13N(0, Element::GeometryType::Pointer(new Pyramid3D13<NodeType >(Element::GeometryType::PointsArrayType(13)))),
      mVolumeCouplingElement3D15N(0, Element::GeometryType::Pointer(new Prism3D15<NodeType >(Element::GeometryType::PointsArrayType(15)))),
      mVolumeCouplingElement3D20N(0, Element::GeometryType::Pointer(new Hexahedra3D20<NodeType >(Element::GeometryType::PointsArrayType(20)))),
      mVolumeCouplingElement3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27<NodeType >(Element::GeometryType::PointsArrayType(27)))),

    // // We define the node type
    // typedef Node NodeType;

    // // STRUCTURAL COUPLING
    // KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DEM_SURFACE_LOAD)
    // KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(BACKUP_LAST_STRUCTURAL_VELOCITY)
    // KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(BACKUP_LAST_STRUCTURAL_DISPLACEMENT)
    // KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SMOOTHED_STRUCTURAL_VELOCITY)

    // KratosDemStructuresCouplingApplication::KratosDemStructuresCouplingApplication()
    //     : KratosApplication("DemStructuresCouplingApplication"),

    //     // Adding line load conditions
    //     mLineLoadFromDEMCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),

    //     // Adding surface load conditions
    //     mSurfaceLoadFromDEMCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))) {}

    // void KratosDemStructuresCouplingApplication::Register() {
    //     // STRUCTURAL COUPLING
    //     KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DEM_SURFACE_LOAD)
    //     KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(BACKUP_LAST_STRUCTURAL_VELOCITY)
    //     KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(BACKUP_LAST_STRUCTURAL_DISPLACEMENT)
    //     KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SMOOTHED_STRUCTURAL_VELOCITY)

    //     // Line loads
    //     KRATOS_REGISTER_CONDITION("LineLoadFromDEMCondition2D2N", mLineLoadFromDEMCondition2D2N)

    //     // Surface loads
    //     KRATOS_REGISTER_CONDITION("SurfaceLoadFromDEMCondition3D3N", mSurfaceLoadFromDEMCondition3D3N)

    //     KRATOS_INFO("Dem-Struct") << std::endl;
    //     KRATOS_INFO("Dem-Struct") << "     KRATOS DEM STRUCTURES COUPLING APPLICATION " << std::endl;
    //     KRATOS_INFO("Dem-Struct") << std::endl;
    //     KRATOS_INFO("Dem-Struct") << "Importing DemStructuresCouplingApplication... ";
    //     KRATOS_INFO("") << " done." << std::endl;
    // }
}  // namespace Kratos
