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

     void FEMDEMVolumeCouplingApplication::Register() {

    KRATOS_REGISTER_ELEMENT("VolumeCouplingElementElement2D3N", mVolumeCouplingElement2D3N)
    KRATOS_REGISTER_ELEMENT("VolumeCouplingElementElement2D4N", mVolumeCouplingElement2D4N)
    KRATOS_REGISTER_ELEMENT("VolumeCouplingElementElement2D6N", mVolumeCouplingElement2D6N)
    KRATOS_REGISTER_ELEMENT("VolumeCouplingElementElement2D8N", mVolumeCouplingElement2D8N)
    KRATOS_REGISTER_ELEMENT("VolumeCouplingElementElement2D9N", mVolumeCouplingElement2D9N)
    KRATOS_REGISTER_ELEMENT("VolumeCouplingElementElement2D10N", mVolumeCouplingElement2D10N)
    KRATOS_REGISTER_ELEMENT("VolumeCouplingElementElement2D15N", mVolumeCouplingElement2D15N)
    KRATOS_REGISTER_ELEMENT("VolumeCouplingElementElement3D4N", mVolumeCouplingElement3D4N)
    KRATOS_REGISTER_ELEMENT("VolumeCouplingElementElement3D5N", mVolumeCouplingElement3D5N)
    KRATOS_REGISTER_ELEMENT("VolumeCouplingElementElement3D6N", mVolumeCouplingElement3D6N)
    KRATOS_REGISTER_ELEMENT("VolumeCouplingElementElement3D8N", mVolumeCouplingElement3D8N)
    KRATOS_REGISTER_ELEMENT("VolumeCouplingElementElement3D10N", mVolumeCouplingElement3D10N)
    KRATOS_REGISTER_ELEMENT("VolumeCouplingElementElement3D13N", mVolumeCouplingElement3D13N)
    KRATOS_REGISTER_ELEMENT("VolumeCouplingElementElement3D15N", mVolumeCouplingElement3D15N)
    KRATOS_REGISTER_ELEMENT("VolumeCouplingElementElement3D20N", mVolumeCouplingElement3D20N)
    KRATOS_REGISTER_ELEMENT("VolumeCouplingElementElement3D27N", mVolumeCouplingElement3D27N)

     }
}  // namespace Kratos
