//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes

// External includes

// Project includes
/* Kratos base classes */
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_flags.h"
#include "includes/kernel.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/constitutive_law.h"
#include "includes/geometrical_object.h"
#include "includes/master_slave_constraint.h"

/* Geometries definition */
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "geometries/point.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "geometries/sphere_3d_1.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"
#include "geometries/quadrature_point_geometry.h"

/* Factories */
#include "factories/standard_linear_solver_factory.h"
#include "factories/standard_preconditioner_factory.h"

namespace Kratos {
typedef Node<3> NodeType;
typedef Geometry<NodeType> GeometryType;

KratosApplication::KratosApplication(const std::string ApplicationName)
    : mApplicationName(ApplicationName),
      // Point conditions
      mPointCondition2D1N( 0, GeometryType::Pointer(new Point2D<NodeType >(GeometryType::PointsArrayType(1)))),
      mPointCondition3D1N( 0, GeometryType::Pointer(new Point3D<NodeType >(GeometryType::PointsArrayType(1)))),
      // Line conditions
      mLineCondition2D2N( 0, GeometryType::Pointer(new Line2D2<NodeType >(GeometryType::PointsArrayType(2)))),
      mLineCondition2D3N( 0, GeometryType::Pointer(new Line2D3<NodeType >(GeometryType::PointsArrayType(3)))),
      mLineCondition3D2N( 0, GeometryType::Pointer(new Line3D2<NodeType >(GeometryType::PointsArrayType(2)))),
      mLineCondition3D3N( 0, GeometryType::Pointer(new Line3D3<NodeType >(GeometryType::PointsArrayType(3)))),
      // Surface conditions
      mSurfaceCondition3D3N( 0, GeometryType::Pointer(new Triangle3D3<NodeType >(GeometryType::PointsArrayType(3)))),
      mSurfaceCondition3D6N( 0, GeometryType::Pointer(new Triangle3D6<NodeType >(GeometryType::PointsArrayType(6)))),
      mSurfaceCondition3D4N( 0, GeometryType::Pointer(new Quadrilateral3D4<NodeType >(GeometryType::PointsArrayType(4)))),
      mSurfaceCondition3D8N( 0, GeometryType::Pointer(new Quadrilateral3D8<NodeType >(GeometryType::PointsArrayType(8)))),
      mSurfaceCondition3D9N( 0, GeometryType::Pointer(new Quadrilateral3D9<NodeType >(GeometryType::PointsArrayType(9)))),

      // Master-Slave Constraint
      mMasterSlaveConstraint(),
      mLinearMasterSlaveConstraint(),

      // Deprecated conditions start
      mCondition2D( 0, GeometryType::Pointer(new Geometry<NodeType >(GeometryType::PointsArrayType(2)))),
      mCondition2D2N( 0, GeometryType::Pointer(new Line2D2<NodeType >(GeometryType::PointsArrayType(2)))),
      mCondition2D3N( 0, GeometryType::Pointer(new Line2D3<NodeType >(GeometryType::PointsArrayType(3)))),
      mCondition3D( 0, GeometryType::Pointer(new Triangle3D3<NodeType >(GeometryType::PointsArrayType(3)))),  // Note: Could be interesting to change the name to mCondition3D3N (conflict with quadratic line)
      mCondition3D2N( 0, GeometryType::Pointer(new Line3D2<NodeType >(GeometryType::PointsArrayType(2)))),
      mCondition3D3N( 0, GeometryType::Pointer(new Line3D3<NodeType >(GeometryType::PointsArrayType(3)))),
      mCondition3D6N( 0, GeometryType::Pointer(new Triangle3D6<NodeType >(GeometryType::PointsArrayType(6)))),
      mCondition3D4N( 0, GeometryType::Pointer(new Quadrilateral3D4<NodeType >(GeometryType::PointsArrayType(4)))),
      mCondition3D8N( 0, GeometryType::Pointer(new Quadrilateral3D8<NodeType >(GeometryType::PointsArrayType(8)))),
      mCondition3D9N( 0, GeometryType::Pointer(new Quadrilateral3D9<NodeType >(GeometryType::PointsArrayType(9)))),
      // Deprecated conditions end

      // Periodic conditions
      mPeriodicCondition( 0, GeometryType::Pointer(new Line2D2<NodeType >(GeometryType::PointsArrayType(2)))),
      mPeriodicConditionEdge( 0, GeometryType::Pointer(new Quadrilateral3D4<NodeType >(GeometryType::PointsArrayType(4)))),
      mPeriodicConditionCorner( 0, GeometryType::Pointer(new Hexahedra3D8<NodeType >(GeometryType::PointsArrayType(8)))),

      // Elements
      mElement2D2N( 0, GeometryType::Pointer(new Line2D2<NodeType >(GeometryType::PointsArrayType(2)))),
      mElement2D3N( 0, GeometryType::Pointer(new Triangle2D3<NodeType >(GeometryType::PointsArrayType(3)))),
      mElement2D4N( 0, GeometryType::Pointer(new Quadrilateral2D4<NodeType >(GeometryType::PointsArrayType(4)))),
      mElement3D2N( 0, GeometryType::Pointer(new Line3D2<NodeType >(GeometryType::PointsArrayType(2)))),
      mElement3D3N( 0, GeometryType::Pointer(new Triangle3D3<NodeType >(GeometryType::PointsArrayType(3)))),
      mElement3D4N( 0, GeometryType::Pointer(new Tetrahedra3D4<NodeType >(GeometryType::PointsArrayType(4)))),
      mElement3D6N( 0, GeometryType::Pointer(new Prism3D6<NodeType >(GeometryType::PointsArrayType(6)))),
      mElement3D8N( 0, GeometryType::Pointer(new Hexahedra3D8<NodeType >(GeometryType::PointsArrayType(8)))),
      mElement3D10N( 0, GeometryType::Pointer(new Tetrahedra3D10<NodeType >(GeometryType::PointsArrayType(10)))),
      mDistanceCalculationElementSimplex2D3N( 0, GeometryType::Pointer(new Triangle2D3<NodeType >(GeometryType::PointsArrayType(3)))),
      mDistanceCalculationElementSimplex3D4N( 0, GeometryType::Pointer(new Tetrahedra3D4<NodeType >(GeometryType::PointsArrayType(4)))),
      mLevelSetConvectionElementSimplex2D3N( 0, GeometryType::Pointer(new Triangle2D3<NodeType >(GeometryType::PointsArrayType(3)))),
      mLevelSetConvectionElementSimplex3D4N( 0, GeometryType::Pointer(new Tetrahedra3D4<NodeType >(GeometryType::PointsArrayType(4)))),

      // Components
      mpVariableData(KratosComponents<VariableData>::pGetComponents()),
      mpIntVariables(KratosComponents<Variable<int> >::pGetComponents()),
      mpUnsignedIntVariables(KratosComponents<Variable<unsigned int> >::pGetComponents()),
      mpDoubleVariables(KratosComponents<Variable<double> >::pGetComponents()),
      mpArray1DVariables(KratosComponents<Variable<array_1d<double, 3> > >::pGetComponents()),
      mpArray1D4Variables(KratosComponents<Variable<array_1d<double, 4> > >::pGetComponents()),
      mpArray1D6Variables(KratosComponents<Variable<array_1d<double, 6> > >::pGetComponents()),
      mpArray1D9Variables(KratosComponents<Variable<array_1d<double, 9> > >::pGetComponents()),
      mpQuaternionVariables(KratosComponents<Variable<Quaternion<double> > >::pGetComponents()),
      mpVectorVariables(KratosComponents<Variable<Vector> >::pGetComponents()),
      mpMatrixVariables(KratosComponents<Variable<Matrix> >::pGetComponents()),
      mpArray1DVariableComponents(KratosComponents<VariableComponent<VectorComponentAdaptor< array_1d<double, 3> > > >::pGetComponents()),
      mpArray1D4VariableComponents(KratosComponents<VariableComponent<VectorComponentAdaptor< array_1d<double, 4> > > >::pGetComponents()),
      mpArray1D6VariableComponents(KratosComponents<VariableComponent<VectorComponentAdaptor< array_1d<double, 6> > > >::pGetComponents()),
      mpArray1D9VariableComponents(KratosComponents<VariableComponent<VectorComponentAdaptor< array_1d<double, 9> > > >::pGetComponents()),
      mpElements(KratosComponents<Element>::pGetComponents()),
      mpConditions(KratosComponents<Condition>::pGetComponents()),
      mpRegisteredObjects(&(Serializer::GetRegisteredObjects())),
      mpRegisteredObjectsName(&(Serializer::GetRegisteredObjectsName())) {}

void KratosApplication::RegisterKratosCore() {

    // Registering all the variables
    KratosApplication::RegisterVariables();

    // Register linear solvers and preconditioners
    RegisterLinearSolvers();
    RegisterPreconditioners();

    //Register objects with general definition
    Serializer::Register("Node", NodeType());
    Serializer::Register("Dof", Dof<double>());
    Serializer::Register("Element", Element());
    Serializer::Register("Condition", Condition());
    Serializer::Register("Properties", Properties());
    Serializer::Register("GeometricalObject", GeometricalObject());

    //Register objects with specific definition ( non essential, must be deleted in future )
    Serializer::Register("Node3D", NodeType());
    Serializer::Register("DofDouble", Dof<double>());

    Serializer::Register("MasterSlaveConstraint", MasterSlaveConstraint());

    //Register specific conditions ( must be completed : conditions defined in kratos_application.h)

    //point conditions
    KRATOS_REGISTER_CONDITION("PointCondition2D1N", mPointCondition2D1N);
    KRATOS_REGISTER_CONDITION("PointCondition3D1N", mPointCondition3D1N);
    //line conditions
    KRATOS_REGISTER_CONDITION("LineCondition2D2N", mLineCondition2D2N);
    KRATOS_REGISTER_CONDITION("LineCondition2D3N", mLineCondition2D3N);
    KRATOS_REGISTER_CONDITION("LineCondition3D2N", mLineCondition3D2N);
    KRATOS_REGISTER_CONDITION("LineCondition3D3N", mLineCondition3D3N);
    //surface conditions
    KRATOS_REGISTER_CONDITION("SurfaceCondition3D3N", mSurfaceCondition3D3N);
    KRATOS_REGISTER_CONDITION("SurfaceCondition3D6N", mSurfaceCondition3D6N);
    KRATOS_REGISTER_CONDITION("SurfaceCondition3D4N", mSurfaceCondition3D4N);
    KRATOS_REGISTER_CONDITION("SurfaceCondition3D8N", mSurfaceCondition3D8N);
    KRATOS_REGISTER_CONDITION("SurfaceCondition3D9N", mSurfaceCondition3D9N);

    //master-slave constraints
    KRATOS_REGISTER_CONSTRAINT("MasterSlaveConstraint",mMasterSlaveConstraint);
    KRATOS_REGISTER_CONSTRAINT("LinearMasterSlaveConstraint",mLinearMasterSlaveConstraint);

    //deprecated conditions start
    KRATOS_REGISTER_CONDITION("Condition2D", mCondition2D);
    KRATOS_REGISTER_CONDITION("Condition2D2N", mCondition2D2N);
    KRATOS_REGISTER_CONDITION("Condition2D3N", mCondition2D3N);
    KRATOS_REGISTER_CONDITION("Condition3D",
        mCondition3D);  // Note: The name could be changed to Condition3D3N (conflict with the quadratic line)
    KRATOS_REGISTER_CONDITION("Condition3D2N", mCondition3D2N);
    KRATOS_REGISTER_CONDITION("Condition3D3N", mCondition3D3N);
    KRATOS_REGISTER_CONDITION("Condition3D6N", mCondition3D6N);
    KRATOS_REGISTER_CONDITION("Condition3D4N", mCondition3D4N);
    KRATOS_REGISTER_CONDITION("Condition3D8N", mCondition3D8N);
    KRATOS_REGISTER_CONDITION("Condition3D9N", mCondition3D9N);
    //deprecated conditions start

    KRATOS_REGISTER_CONDITION("PeriodicCondition", mPeriodicCondition)
    KRATOS_REGISTER_CONDITION("PeriodicConditionEdge", mPeriodicConditionEdge)
    KRATOS_REGISTER_CONDITION("PeriodicConditionCorner", mPeriodicConditionCorner)

    //Register specific elements ( must be completed : elements defined in kratos_appliction.h)
    KRATOS_REGISTER_ELEMENT("Element2D2N", mElement2D2N)
    KRATOS_REGISTER_ELEMENT("Element2D3N", mElement2D3N)
    KRATOS_REGISTER_ELEMENT("Element2D4N", mElement2D4N)
    KRATOS_REGISTER_ELEMENT("Element3D2N", mElement3D2N)
    KRATOS_REGISTER_ELEMENT("Element3D3N", mElement3D3N)
    KRATOS_REGISTER_ELEMENT("Element3D4N", mElement3D4N)
    KRATOS_REGISTER_ELEMENT("Element3D6N", mElement3D6N)
    KRATOS_REGISTER_ELEMENT("Element3D8N", mElement3D8N)
    KRATOS_REGISTER_ELEMENT("Element3D10N", mElement3D10N)

    KRATOS_REGISTER_ELEMENT("DistanceCalculationElementSimplex2D3N", mDistanceCalculationElementSimplex2D3N)
    KRATOS_REGISTER_ELEMENT("DistanceCalculationElementSimplex3D4N", mDistanceCalculationElementSimplex3D4N)
    KRATOS_REGISTER_ELEMENT("LevelSetConvectionElementSimplex2D3N", mLevelSetConvectionElementSimplex2D3N)
    KRATOS_REGISTER_ELEMENT("LevelSetConvectionElementSimplex3D4N", mLevelSetConvectionElementSimplex3D4N)

    //Register general geometries:

    //Points:
    Serializer::Register("Point", Point());

    Point2D<NodeType > Point2DPrototype(GeometryType::PointsArrayType(1));
    Serializer::Register("Point2D", Point2DPrototype);

    Point3D<NodeType > Point3DPrototype(GeometryType::PointsArrayType(1));
    Serializer::Register("Point3D", Point3DPrototype);

    //Sphere
    Sphere3D1<NodeType > Sphere3D1Prototype(GeometryType::PointsArrayType(1));
    Serializer::Register("Sphere3D1", Sphere3D1Prototype);

    //Lines:
    Line2D2<NodeType > Line2D2Prototype(GeometryType::PointsArrayType(2));
    Serializer::Register("Line2D2", Line2D2Prototype);

    Line2D3<NodeType > Line2D3Prototype(GeometryType::PointsArrayType(3));
    Serializer::Register("Line2D3", Line2D3Prototype);

    Line3D2<NodeType > Line3D2Prototype(GeometryType::PointsArrayType(2));
    Serializer::Register("Line3D2", Line3D2Prototype);

    Line3D3<NodeType > Line3D3Prototype(GeometryType::PointsArrayType(3));
    Serializer::Register("Line3D3", Line3D3Prototype);

    //Triangles:
    Triangle2D3<NodeType > Triangle2D3Prototype(GeometryType::PointsArrayType(3));
    Serializer::Register("Triangle2D3", Triangle2D3Prototype);

    Triangle2D6<NodeType > Triangle2D6Prototype(GeometryType::PointsArrayType(6));
    Serializer::Register("Triangle2D6", Triangle2D6Prototype);

    Triangle3D3<NodeType > Triangle3D3Prototype(GeometryType::PointsArrayType(3));
    Serializer::Register("Triangle3D3", Triangle3D3Prototype);

    Triangle3D6<NodeType > Triangle3D6Prototype( GeometryType::PointsArrayType(6));
    Serializer::Register("Triangle3D6", Triangle3D6Prototype);

    //Quadrilaterals:
    Quadrilateral2D4<NodeType > Quadrilateral2D4Prototype( GeometryType::PointsArrayType(4));
    Serializer::Register("Quadrilateral2D4", Quadrilateral2D4Prototype);

    Quadrilateral2D8<NodeType > Quadrilateral2D8Prototype( GeometryType::PointsArrayType(8));
    Serializer::Register("Quadrilateral2D8", Quadrilateral2D8Prototype);

    Quadrilateral2D9<NodeType > Quadrilateral2D9Prototype( GeometryType::PointsArrayType(9));
    Serializer::Register("Quadrilateral2D9", Quadrilateral2D9Prototype);

    Quadrilateral3D4<NodeType > Quadrilateral3D4Prototype( GeometryType::PointsArrayType(4));
    Serializer::Register("Quadrilateral3D4", Quadrilateral3D4Prototype);

    Quadrilateral3D8<NodeType > Quadrilateral3D8Prototype( GeometryType::PointsArrayType(8));
    Serializer::Register("Quadrilateral3D8", Quadrilateral3D8Prototype);

    Quadrilateral3D9<NodeType > Quadrilateral3D9Prototype( GeometryType::PointsArrayType(9));
    Serializer::Register("Quadrilateral3D9", Quadrilateral3D9Prototype);

    //Tetrahedra:
    Tetrahedra3D4<NodeType > Tetrahedra3D4Prototype( GeometryType::PointsArrayType(4));
    Serializer::Register("Tetrahedra3D4", Tetrahedra3D4Prototype);

    Tetrahedra3D10<NodeType > Tetrahedra3D10Prototype( GeometryType::PointsArrayType(10));
    Serializer::Register("Tetrahedra3D10", Tetrahedra3D10Prototype);

    //Prisms:
    Prism3D6<NodeType > Prism3D6Prototype( GeometryType::PointsArrayType(6));
    Serializer::Register("Prism3D6", Prism3D6Prototype);

    Prism3D15<NodeType > Prism3D15Prototype( GeometryType::PointsArrayType(15));
    Serializer::Register("Prism3D15", Prism3D15Prototype);

    //Hexahedra:
    Hexahedra3D8<NodeType > Hexahedra3D8Prototype( GeometryType::PointsArrayType(8));
    Serializer::Register("Hexahedra3D8", Hexahedra3D8Prototype);

    Hexahedra3D20<NodeType > Hexahedra3D20Prototype( GeometryType::PointsArrayType(20));
    Serializer::Register("Hexahedra3D20", Hexahedra3D20Prototype);

    Hexahedra3D27<NodeType > Hexahedra3D27Prototype( GeometryType::PointsArrayType(27));
    Serializer::Register("Hexahedra3D27", Hexahedra3D27Prototype);

    //QuadraturePointGeometry:
    QuadraturePointGeometry<NodeType, 3> QuadraturePointGeometryVolume3dPrototype(GeometryType::PointsArrayType(1));
    Serializer::Register("QuadraturePointGeometryVolume3d", QuadraturePointGeometryVolume3dPrototype);

    QuadraturePointGeometry< NodeType, 3, 2 > QuadraturePointGeometrySurface3dPrototype(GeometryType::PointsArrayType(1));
    Serializer::Register("QuadraturePointGeometrySurface3d", QuadraturePointGeometrySurface3dPrototype);

    QuadraturePointGeometry< NodeType, 2 > QuadraturePointGeometrySurface2dPrototype(GeometryType::PointsArrayType(1));
    Serializer::Register("QuadraturePointGeometrySurface2d", QuadraturePointGeometrySurface2dPrototype);

    QuadraturePointGeometry< NodeType, 3, 1 > QuadraturePointGeometryCurve3dPrototype(GeometryType::PointsArrayType(1));
    Serializer::Register("QuadraturePointGeometryCurve3d", QuadraturePointGeometryCurve3dPrototype);

    QuadraturePointGeometry< NodeType, 2, 1 > QuadraturePointGeometryCurve2dPrototype(GeometryType::PointsArrayType(1));
    Serializer::Register("QuadraturePointGeometryCurve2d", QuadraturePointGeometryCurve2dPrototype);

    // Register flags:
    KRATOS_REGISTER_FLAG(STRUCTURE);
    KRATOS_REGISTER_FLAG(FLUID);
    KRATOS_REGISTER_FLAG(THERMAL);
    KRATOS_REGISTER_FLAG(VISITED);
    KRATOS_REGISTER_FLAG(SELECTED);
    KRATOS_REGISTER_FLAG(BOUNDARY);
    KRATOS_REGISTER_FLAG(INLET);
    KRATOS_REGISTER_FLAG(OUTLET);
    KRATOS_REGISTER_FLAG(SLIP);
    KRATOS_REGISTER_FLAG(INTERFACE);
    KRATOS_REGISTER_FLAG(CONTACT);
    KRATOS_REGISTER_FLAG(TO_SPLIT);
    KRATOS_REGISTER_FLAG(TO_ERASE);
    KRATOS_REGISTER_FLAG(TO_REFINE);
    KRATOS_REGISTER_FLAG(NEW_ENTITY);
    KRATOS_REGISTER_FLAG(OLD_ENTITY);
    KRATOS_REGISTER_FLAG(ACTIVE);
    KRATOS_REGISTER_FLAG(MODIFIED);
    KRATOS_REGISTER_FLAG(RIGID);
    KRATOS_REGISTER_FLAG(SOLID);
    KRATOS_REGISTER_FLAG(MPI_BOUNDARY);
    KRATOS_REGISTER_FLAG(INTERACTION);
    KRATOS_REGISTER_FLAG(ISOLATED);
    KRATOS_REGISTER_FLAG(MASTER);
    KRATOS_REGISTER_FLAG(SLAVE);
    KRATOS_REGISTER_FLAG(INSIDE);
    KRATOS_REGISTER_FLAG(FREE_SURFACE);
    KRATOS_REGISTER_FLAG(BLOCKED);
    KRATOS_REGISTER_FLAG(MARKER);
    KRATOS_REGISTER_FLAG(PERIODIC);

    // Note: using internal macro for these two because they do not have a NOT_ version
    KRATOS_ADD_FLAG_TO_KRATOS_COMPONENTS(ALL_DEFINED);
    KRATOS_ADD_FLAG_TO_KRATOS_COMPONENTS(ALL_TRUE);

    // Register ConstitutiveLaw BaseClass
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ConstitutiveLaw", mConstitutiveLaw);
}
}  // namespace Kratos.
