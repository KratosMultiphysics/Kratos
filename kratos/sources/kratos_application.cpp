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

/* Factories */
#include "factories/standard_linear_solver_factory.h"
#include "factories/standard_preconditioner_factory.h"

namespace Kratos {

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

      // Periodic conditions
      mPeriodicCondition( 0, GeometryType::Pointer(new Line2D2<NodeType >(GeometryType::PointsArrayType(2)))),
      mPeriodicConditionEdge( 0, GeometryType::Pointer(new Quadrilateral3D4<NodeType >(GeometryType::PointsArrayType(4)))),
      mPeriodicConditionCorner( 0, GeometryType::Pointer(new Hexahedra3D8<NodeType >(GeometryType::PointsArrayType(8)))),

      // Elements
      mElement2D2N( 0, GeometryType::Pointer(new Line2D2<NodeType >(GeometryType::PointsArrayType(2)))),
      mElement2D3N( 0, GeometryType::Pointer(new Triangle2D3<NodeType >(GeometryType::PointsArrayType(3)))),
      mElement2D6N( 0, GeometryType::Pointer(new Triangle2D6<NodeType >(GeometryType::PointsArrayType(6)))),
      mElement2D4N( 0, GeometryType::Pointer(new Quadrilateral2D4<NodeType >(GeometryType::PointsArrayType(4)))),
      mElement2D8N( 0, GeometryType::Pointer(new Quadrilateral2D8<NodeType >(GeometryType::PointsArrayType(8)))),
      mElement2D9N( 0, GeometryType::Pointer(new Quadrilateral2D9<NodeType >(GeometryType::PointsArrayType(9)))),
      mElement3D2N( 0, GeometryType::Pointer(new Line3D2<NodeType >(GeometryType::PointsArrayType(2)))),
      mElement3D3N( 0, GeometryType::Pointer(new Triangle3D3<NodeType >(GeometryType::PointsArrayType(3)))),
      mElement3D4N( 0, GeometryType::Pointer(new Tetrahedra3D4<NodeType >(GeometryType::PointsArrayType(4)))),
      mElement3D6N( 0, GeometryType::Pointer(new Prism3D6<NodeType >(GeometryType::PointsArrayType(6)))),
      mElement3D8N( 0, GeometryType::Pointer(new Hexahedra3D8<NodeType >(GeometryType::PointsArrayType(8)))),
      mElement3D10N( 0, GeometryType::Pointer(new Tetrahedra3D10<NodeType >(GeometryType::PointsArrayType(10)))),
      mElement3D15N( 0, GeometryType::Pointer(new Prism3D15<NodeType >(GeometryType::PointsArrayType(15)))),
      mElement3D20N( 0, GeometryType::Pointer(new Hexahedra3D20<NodeType >(GeometryType::PointsArrayType(20)))),
      mElement3D27N( 0, GeometryType::Pointer(new Hexahedra3D27<NodeType >(GeometryType::PointsArrayType(27)))),
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
      mpGeometries(KratosComponents<GeometryType>::pGetComponents()),
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
    Serializer::Register("Geometry", GeometryType());
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

    KRATOS_REGISTER_CONDITION("PeriodicCondition", mPeriodicCondition)
    KRATOS_REGISTER_CONDITION("PeriodicConditionEdge", mPeriodicConditionEdge)
    KRATOS_REGISTER_CONDITION("PeriodicConditionCorner", mPeriodicConditionCorner)

    //Register specific elements ( must be completed : elements defined in kratos_appliction.h)
    KRATOS_REGISTER_ELEMENT("Element2D2N", mElement2D2N)
    KRATOS_REGISTER_ELEMENT("Element2D3N", mElement2D3N)
    KRATOS_REGISTER_ELEMENT("Element2D6N", mElement2D6N)
    KRATOS_REGISTER_ELEMENT("Element2D4N", mElement2D4N)
    KRATOS_REGISTER_ELEMENT("Element2D8N", mElement2D8N)
    KRATOS_REGISTER_ELEMENT("Element2D9N", mElement2D9N)
    KRATOS_REGISTER_ELEMENT("Element3D2N", mElement3D2N)
    KRATOS_REGISTER_ELEMENT("Element3D3N", mElement3D3N)
    KRATOS_REGISTER_ELEMENT("Element3D4N", mElement3D4N)
    KRATOS_REGISTER_ELEMENT("Element3D6N", mElement3D6N)
    KRATOS_REGISTER_ELEMENT("Element3D8N", mElement3D8N)
    KRATOS_REGISTER_ELEMENT("Element3D10N", mElement3D10N)
    KRATOS_REGISTER_ELEMENT("Element3D15N", mElement3D15N)
    KRATOS_REGISTER_ELEMENT("Element3D20N", mElement3D20N)
    KRATOS_REGISTER_ELEMENT("Element3D27N", mElement3D27N)

    KRATOS_REGISTER_ELEMENT("DistanceCalculationElementSimplex2D3N", mDistanceCalculationElementSimplex2D3N)
    KRATOS_REGISTER_ELEMENT("DistanceCalculationElementSimplex3D4N", mDistanceCalculationElementSimplex3D4N)
    KRATOS_REGISTER_ELEMENT("LevelSetConvectionElementSimplex2D3N", mLevelSetConvectionElementSimplex2D3N)
    KRATOS_REGISTER_ELEMENT("LevelSetConvectionElementSimplex3D4N", mLevelSetConvectionElementSimplex3D4N)

    //Register general geometries:
    // Point register:
    Serializer::Register("Point", mPointPrototype);
    
    // Register + KratosComponents
    KRATOS_REGISTER_GEOMETRY("Point2D", mPoint2DPrototype);
    KRATOS_REGISTER_GEOMETRY("Point3D", mPoint3DPrototype);
    KRATOS_REGISTER_GEOMETRY("Sphere3D1", mSphere3D1Prototype);
    KRATOS_REGISTER_GEOMETRY("Line2D2", mLine2D2Prototype);
    KRATOS_REGISTER_GEOMETRY("Line2D3", mLine2D3Prototype);
    KRATOS_REGISTER_GEOMETRY("Line3D2", mLine3D2Prototype);
    KRATOS_REGISTER_GEOMETRY("Line3D3", mLine3D3Prototype);
    KRATOS_REGISTER_GEOMETRY("Triangle2D3", mTriangle2D3Prototype);
    KRATOS_REGISTER_GEOMETRY("Triangle2D6", mTriangle2D6Prototype);
    KRATOS_REGISTER_GEOMETRY("Triangle3D3", mTriangle3D3Prototype);
    KRATOS_REGISTER_GEOMETRY("Triangle3D6", mTriangle3D6Prototype);
    KRATOS_REGISTER_GEOMETRY("Quadrilateral2D4", mQuadrilateral2D4Prototype);
    KRATOS_REGISTER_GEOMETRY("Quadrilateral2D8", mQuadrilateral2D8Prototype);
    KRATOS_REGISTER_GEOMETRY("Quadrilateral2D9", mQuadrilateral2D9Prototype);
    KRATOS_REGISTER_GEOMETRY("Quadrilateral3D4", mQuadrilateral3D4Prototype);
    KRATOS_REGISTER_GEOMETRY("Quadrilateral3D8", mQuadrilateral3D8Prototype);
    KRATOS_REGISTER_GEOMETRY("Quadrilateral3D9", mQuadrilateral3D9Prototype);
    KRATOS_REGISTER_GEOMETRY("Tetrahedra3D4", mTetrahedra3D4Prototype);
    KRATOS_REGISTER_GEOMETRY("Tetrahedra3D10", mTetrahedra3D10Prototype);
    KRATOS_REGISTER_GEOMETRY("Prism3D6", mPrism3D6Prototype);
    KRATOS_REGISTER_GEOMETRY("Prism3D15", mPrism3D15Prototype);
    KRATOS_REGISTER_GEOMETRY("Hexahedra3D8", mHexahedra3D8Prototype);
    KRATOS_REGISTER_GEOMETRY("Hexahedra3D20", mHexahedra3D20Prototype);
    KRATOS_REGISTER_GEOMETRY("Hexahedra3D27", mHexahedra3D27Prototype);

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
