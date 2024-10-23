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

#include "modeler/coloring/voxel_mesher_coloring_factory.h"
#include "modeler/key_plane_generation/key_plane_generation_factory.h"
#include "modeler/entity_generation/voxel_mesher_entity_generation_factory.h"
#include "modeler/operation/voxel_mesher_operation_factory.h"

namespace Kratos {

KratosApplication::KratosApplication(const std::string& ApplicationName)
    : mApplicationName(ApplicationName),
      // Generic condition
      mGenericCondition( 0, GeometryType::Pointer(new Geometry<NodeType>())),
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
      //prisms
      mPrismCondition2D4N( 0, GeometryType::Pointer(new Quadrilateral2D4<NodeType >(GeometryType::PointsArrayType(4)))),
      mPrismCondition3D6N( 0, GeometryType::Pointer(new Prism3D6<NodeType >(GeometryType::PointsArrayType(6)))),

      // Master-Slave Constraint
      mMasterSlaveConstraint(),
      mLinearMasterSlaveConstraint(),

      // Periodic conditions
      mPeriodicCondition( 0, GeometryType::Pointer(new Line2D2<NodeType >(GeometryType::PointsArrayType(2)))),
      mPeriodicConditionEdge( 0, GeometryType::Pointer(new Quadrilateral3D4<NodeType >(GeometryType::PointsArrayType(4)))),
      mPeriodicConditionCorner( 0, GeometryType::Pointer(new Hexahedra3D8<NodeType >(GeometryType::PointsArrayType(8)))),

      // Elements
      mGenericElement( 0, GeometryType::Pointer(new Geometry<NodeType>())),

      mElement2D1N( 0, GeometryType::Pointer(new Point2D<NodeType >(GeometryType::PointsArrayType(1)))),
      mElement2D2N( 0, GeometryType::Pointer(new Line2D2<NodeType >(GeometryType::PointsArrayType(2)))),
      mElement2D3N( 0, GeometryType::Pointer(new Triangle2D3<NodeType >(GeometryType::PointsArrayType(3)))),
      mElement2D6N( 0, GeometryType::Pointer(new Triangle2D6<NodeType >(GeometryType::PointsArrayType(6)))),
      mElement2D4N( 0, GeometryType::Pointer(new Quadrilateral2D4<NodeType >(GeometryType::PointsArrayType(4)))),
      mElement2D8N( 0, GeometryType::Pointer(new Quadrilateral2D8<NodeType >(GeometryType::PointsArrayType(8)))),
      mElement2D9N( 0, GeometryType::Pointer(new Quadrilateral2D9<NodeType >(GeometryType::PointsArrayType(9)))),

      mElement3D1N( 0, GeometryType::Pointer(new Point3D<NodeType >(GeometryType::PointsArrayType(1)))),
      mElement3D2N( 0, GeometryType::Pointer(new Line3D2<NodeType >(GeometryType::PointsArrayType(2)))),
      mElement3D3N( 0, GeometryType::Pointer(new Triangle3D3<NodeType >(GeometryType::PointsArrayType(3)))),
      mElement3D4N( 0, GeometryType::Pointer(new Tetrahedra3D4<NodeType >(GeometryType::PointsArrayType(4)))),
      mElement3D5N( 0, GeometryType::Pointer(new Pyramid3D5<NodeType >(GeometryType::PointsArrayType(5)))),
      mElement3D6N( 0, GeometryType::Pointer(new Prism3D6<NodeType >(GeometryType::PointsArrayType(6)))),
      mElement3D8N( 0, GeometryType::Pointer(new Hexahedra3D8<NodeType >(GeometryType::PointsArrayType(8)))),
      mElement3D10N( 0, GeometryType::Pointer(new Tetrahedra3D10<NodeType >(GeometryType::PointsArrayType(10)))),
      mElement3D13N( 0, GeometryType::Pointer(new Pyramid3D13<NodeType >(GeometryType::PointsArrayType(13)))),
      mElement3D15N( 0, GeometryType::Pointer(new Prism3D15<NodeType >(GeometryType::PointsArrayType(15)))),
      mElement3D20N( 0, GeometryType::Pointer(new Hexahedra3D20<NodeType >(GeometryType::PointsArrayType(20)))),
      mElement3D27N( 0, GeometryType::Pointer(new Hexahedra3D27<NodeType >(GeometryType::PointsArrayType(27)))),

      mPointElement2D1N( 0, GeometryType::Pointer(new Point2D<NodeType >(GeometryType::PointsArrayType(1)))),
      mPointElement3D1N( 0, GeometryType::Pointer(new Point3D<NodeType >(GeometryType::PointsArrayType(1)))),

      mLineElement2D2N( 0, GeometryType::Pointer(new Line2D2<NodeType >(GeometryType::PointsArrayType(2)))),
      mLineElement2D3N( 0, GeometryType::Pointer(new Line2D3<NodeType >(GeometryType::PointsArrayType(3)))),
      mLineElement3D2N( 0, GeometryType::Pointer(new Line3D2<NodeType >(GeometryType::PointsArrayType(2)))),
      mLineElement3D3N( 0, GeometryType::Pointer(new Line3D3<NodeType >(GeometryType::PointsArrayType(3)))),

      mSurfaceElement3D3N( 0, GeometryType::Pointer(new Triangle3D3<NodeType >(GeometryType::PointsArrayType(3)))),
      mSurfaceElement3D6N( 0, GeometryType::Pointer(new Triangle3D6<NodeType >(GeometryType::PointsArrayType(6)))),
      mSurfaceElement3D4N( 0, GeometryType::Pointer(new Quadrilateral3D4<NodeType >(GeometryType::PointsArrayType(4)))),
      mSurfaceElement3D8N( 0, GeometryType::Pointer(new Quadrilateral3D8<NodeType >(GeometryType::PointsArrayType(8)))),
      mSurfaceElement3D9N( 0, GeometryType::Pointer(new Quadrilateral3D9<NodeType >(GeometryType::PointsArrayType(9)))),

      mDistanceCalculationElementSimplex2D3N( 0, GeometryType::Pointer(new Triangle2D3<NodeType >(GeometryType::PointsArrayType(3)))),
      mDistanceCalculationElementSimplex3D4N( 0, GeometryType::Pointer(new Tetrahedra3D4<NodeType >(GeometryType::PointsArrayType(4)))),

      mEdgeBasedGradientRecoveryElement2D2N( 0, GeometryType::Pointer(new Line2D2<NodeType>(GeometryType::PointsArrayType(2)))),
      mEdgeBasedGradientRecoveryElement3D2N( 0, GeometryType::Pointer(new Line3D2<NodeType>(GeometryType::PointsArrayType(2)))),

      mLevelSetConvectionElementSimplex2D3N( 0, GeometryType::Pointer(new Triangle2D3<NodeType >(GeometryType::PointsArrayType(3)))),
      mLevelSetConvectionElementSimplex3D4N( 0, GeometryType::Pointer(new Tetrahedra3D4<NodeType >(GeometryType::PointsArrayType(4)))),
      mLevelSetConvectionElementSimplexAlgebraicStabilization2D3N( 0, GeometryType::Pointer(new Triangle2D3<NodeType >(GeometryType::PointsArrayType(3)))),
      mLevelSetConvectionElementSimplexAlgebraicStabilization3D4N( 0, GeometryType::Pointer(new Tetrahedra3D4<NodeType >(GeometryType::PointsArrayType(4)))),

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
      mpGeometries(KratosComponents<GeometryType>::pGetComponents()),
      mpElements(KratosComponents<Element>::pGetComponents()),
      mpConditions(KratosComponents<Condition>::pGetComponents()),
      mpModelers(KratosComponents<Modeler>::pGetComponents()),
      mpRegisteredObjects(&(Serializer::GetRegisteredObjects())),
      mpRegisteredObjectsName(&(Serializer::GetRegisteredObjectsName())) {

        Registry::SetCurrentSource(mApplicationName);

        for (auto component : {"geometries", "elements", "conditions", "constraints", "modelers", "constitutive_laws"}) {
            if (!Registry::HasItem(std::string(component))) {
                Registry::AddItem<RegistryItem>(std::string(component)+"."+mApplicationName);
            }
        }
      }

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
    Serializer::Register("Modeler", Modeler());
    Serializer::Register("Properties", Properties());
    Serializer::Register("GeometricalObject", GeometricalObject());

    //Register objects with specific definition ( non essential, must be deleted in future )
    Serializer::Register("Node3D", NodeType());
    Serializer::Register("DofDouble", Dof<double>());

    Serializer::Register("MasterSlaveConstraint", MasterSlaveConstraint());

    //Register specific conditions ( must be completed : conditions defined in kratos_application.h)
    //generic condition
    KRATOS_REGISTER_CONDITION("GenericCondition", mGenericCondition);
    
    // Point conditions
    KRATOS_REGISTER_CONDITION("PointCondition2D1N", mPointCondition2D1N);
    KRATOS_REGISTER_CONDITION("PointCondition3D1N", mPointCondition3D1N);

    // Line conditions
    KRATOS_REGISTER_CONDITION("LineCondition2D2N", mLineCondition2D2N);
    KRATOS_REGISTER_CONDITION("LineCondition2D3N", mLineCondition2D3N);
    KRATOS_REGISTER_CONDITION("LineCondition3D2N", mLineCondition3D2N);
    KRATOS_REGISTER_CONDITION("LineCondition3D3N", mLineCondition3D3N);

    // Surface conditions
    KRATOS_REGISTER_CONDITION("SurfaceCondition3D3N", mSurfaceCondition3D3N);
    KRATOS_REGISTER_CONDITION("SurfaceCondition3D6N", mSurfaceCondition3D6N);
    KRATOS_REGISTER_CONDITION("SurfaceCondition3D4N", mSurfaceCondition3D4N);
    KRATOS_REGISTER_CONDITION("SurfaceCondition3D8N", mSurfaceCondition3D8N);
    KRATOS_REGISTER_CONDITION("SurfaceCondition3D9N", mSurfaceCondition3D9N);

    //Prism conditions
    KRATOS_REGISTER_CONDITION("PrismCondition2D4N", mPrismCondition2D4N);
    KRATOS_REGISTER_CONDITION("PrismCondition3D6N", mPrismCondition3D6N);

    //Master-slave constraints
    KRATOS_REGISTER_CONSTRAINT("MasterSlaveConstraint",mMasterSlaveConstraint);
    KRATOS_REGISTER_CONSTRAINT("LinearMasterSlaveConstraint",mLinearMasterSlaveConstraint);

    KRATOS_REGISTER_CONDITION("PeriodicCondition", mPeriodicCondition)
    KRATOS_REGISTER_CONDITION("PeriodicConditionEdge", mPeriodicConditionEdge)
    KRATOS_REGISTER_CONDITION("PeriodicConditionCorner", mPeriodicConditionCorner)

    //Register specific elements ( must be completed : elements defined in kratos_application.h)
    KRATOS_REGISTER_ELEMENT("GenericElement", mGenericElement)

    KRATOS_REGISTER_ELEMENT("Element2D1N", mElement2D1N)
    KRATOS_REGISTER_ELEMENT("Element2D2N", mElement2D2N)
    KRATOS_REGISTER_ELEMENT("Element2D3N", mElement2D3N)
    KRATOS_REGISTER_ELEMENT("Element2D6N", mElement2D6N)
    KRATOS_REGISTER_ELEMENT("Element2D4N", mElement2D4N)
    KRATOS_REGISTER_ELEMENT("Element2D8N", mElement2D8N)
    KRATOS_REGISTER_ELEMENT("Element2D9N", mElement2D9N)

    KRATOS_REGISTER_ELEMENT("Element3D1N", mElement3D1N)
    KRATOS_REGISTER_ELEMENT("Element3D2N", mElement3D2N)
    KRATOS_REGISTER_ELEMENT("Element3D3N", mElement3D3N)
    KRATOS_REGISTER_ELEMENT("Element3D4N", mElement3D4N)
    KRATOS_REGISTER_ELEMENT("Element3D5N", mElement3D5N)
    KRATOS_REGISTER_ELEMENT("Element3D6N", mElement3D6N)
    KRATOS_REGISTER_ELEMENT("Element3D8N", mElement3D8N)
    KRATOS_REGISTER_ELEMENT("Element3D10N", mElement3D10N)
    KRATOS_REGISTER_ELEMENT("Element3D13N", mElement3D13N)
    KRATOS_REGISTER_ELEMENT("Element3D15N", mElement3D15N)
    KRATOS_REGISTER_ELEMENT("Element3D20N", mElement3D20N)
    KRATOS_REGISTER_ELEMENT("Element3D27N", mElement3D27N)

    KRATOS_REGISTER_ELEMENT("PointElement2D1N", mPointElement2D1N)
    KRATOS_REGISTER_ELEMENT("PointElement3D1N", mPointElement3D1N)

    KRATOS_REGISTER_ELEMENT("LineElement2D2N", mLineElement2D2N)
    KRATOS_REGISTER_ELEMENT("LineElement2D3N", mLineElement2D3N)
    KRATOS_REGISTER_ELEMENT("LineElement3D2N", mLineElement3D2N)
    KRATOS_REGISTER_ELEMENT("LineElement3D3N", mLineElement3D3N)

    KRATOS_REGISTER_ELEMENT("SurfaceElement3D3N", mSurfaceElement3D3N)
    KRATOS_REGISTER_ELEMENT("SurfaceElement3D6N", mSurfaceElement3D6N)
    KRATOS_REGISTER_ELEMENT("SurfaceElement3D4N", mSurfaceElement3D4N)
    KRATOS_REGISTER_ELEMENT("SurfaceElement3D8N", mSurfaceElement3D8N)
    KRATOS_REGISTER_ELEMENT("SurfaceElement3D9N", mSurfaceElement3D9N)

    KRATOS_REGISTER_ELEMENT("DistanceCalculationElementSimplex2D3N", mDistanceCalculationElementSimplex2D3N)
    KRATOS_REGISTER_ELEMENT("DistanceCalculationElementSimplex3D4N", mDistanceCalculationElementSimplex3D4N)

    KRATOS_REGISTER_ELEMENT("EdgeBasedGradientRecoveryElement2D2N", mEdgeBasedGradientRecoveryElement2D2N)
    KRATOS_REGISTER_ELEMENT("EdgeBasedGradientRecoveryElement3D2N", mEdgeBasedGradientRecoveryElement3D2N)

    KRATOS_REGISTER_ELEMENT("LevelSetConvectionElementSimplex2D3N", mLevelSetConvectionElementSimplex2D3N)
    KRATOS_REGISTER_ELEMENT("LevelSetConvectionElementSimplex3D4N", mLevelSetConvectionElementSimplex3D4N)
    KRATOS_REGISTER_ELEMENT("LevelSetConvectionElementSimplexAlgebraicStabilization2D3N", mLevelSetConvectionElementSimplexAlgebraicStabilization2D3N)
    KRATOS_REGISTER_ELEMENT("LevelSetConvectionElementSimplexAlgebraicStabilization3D4N", mLevelSetConvectionElementSimplexAlgebraicStabilization3D4N)

    KRATOS_REGISTER_MODELER("Modeler", mModeler);
    KRATOS_REGISTER_MODELER("CadIoModeler", mCadIoModeler);
#if USE_TRIANGLE_NONFREE_TPL
    KRATOS_REGISTER_MODELER("CadTessellationModeler", mCadTessellationModeler);
#endif
    KRATOS_REGISTER_MODELER("SerialModelPartCombinatorModeler", mSerialModelPartCombinatorModeler);
    KRATOS_REGISTER_MODELER("CombineModelPartModeler", mCombineModelPartModeler);
    KRATOS_REGISTER_MODELER("ConnectivityPreserveModeler", mConnectivityPreserveModeler);
    KRATOS_REGISTER_MODELER("VoxelMeshGeneratorModeler", mVoxelMeshGeneratorModeler);
    KRATOS_REGISTER_MODELER("CleanUpProblematicTrianglesModeler", mCleanUpProblematicTrianglesModeler);

    // Register general geometries:
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
    KRATOS_REGISTER_GEOMETRY("Pyramid3D5", mPyramid3D5Prototype);
    KRATOS_REGISTER_GEOMETRY("Pyramid3D13", mPyramid3D13Prototype);
    KRATOS_REGISTER_GEOMETRY("Hexahedra3D8", mHexahedra3D8Prototype);
    KRATOS_REGISTER_GEOMETRY("Hexahedra3D20", mHexahedra3D20Prototype);
    KRATOS_REGISTER_GEOMETRY("Hexahedra3D27", mHexahedra3D27Prototype);
    KRATOS_REGISTER_GEOMETRY("QuadraturePointGeometryPoint1D", mQuadraturePointGeometryPoint1D);
    KRATOS_REGISTER_GEOMETRY("QuadraturePointGeometryPoint2D", mQuadraturePointGeometryPoint2D);
    KRATOS_REGISTER_GEOMETRY("QuadraturePointGeometryPoint3D", mQuadraturePointGeometryPoint3D);
    KRATOS_REGISTER_GEOMETRY("QuadraturePointGeometrySurface2D", mQuadraturePointGeometrySurface2D);
    KRATOS_REGISTER_GEOMETRY("QuadraturePointGeometrySurface3D", mQuadraturePointGeometrySurface3D);
    KRATOS_REGISTER_GEOMETRY("QuadraturePointGeometryVolume3D", mQuadraturePointGeometryVolume3D);
    KRATOS_REGISTER_GEOMETRY("CouplingGeometry", mCouplingGeometry);

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
    KRATOS_REGISTER_FLAG(WALL);

    // Note: using internal macro for these two because they do not have a NOT_ version
    KRATOS_ADD_FLAG_TO_KRATOS_COMPONENTS(ALL_DEFINED);
    KRATOS_ADD_FLAG_TO_KRATOS_COMPONENTS(ALL_TRUE);

    // Register ConstitutiveLaw BaseClass
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ConstitutiveLaw", mConstitutiveLaw);

    //Register Voxel Modeler modular components
    RegisterVoxelMesherColoring();
    RegisterVoxelMesherKeyPlaneGeneration();
    RegisterVoxelMesherEntityGeneration();
    RegisterVoxelMesherOperation();
}

template<class TComponentsContainer>
void KratosApplication::DeregisterComponent(std::string const & rComponentName) {
    auto path = std::string(rComponentName)+"."+mApplicationName;

    // Remove only if the application has this type of components registered
    if (Registry::HasItem(path)) {

        // Generate a temporal list with all the keys to avoid invalidating the iterator (Convert this into a transform range when C++20 is available)
        std::vector<std::string> keys;
        std::transform(Registry::GetItem(path).cbegin(), Registry::GetItem(path).cend(), std::back_inserter(keys), [](auto & key){return std::string(key.first);});

        for (auto & key : keys) {
            auto cmpt_key = "components."+key;
            auto type_key = path+"."+key;

            // Remove from KratosComponents
            KratosComponents<TComponentsContainer>::Remove(key);

            // Remove from registry general component list
            if (Registry::HasItem(cmpt_key)) {
                Registry::RemoveItem(cmpt_key);
            } else {
                KRATOS_ERROR << "Trying ro remove: " << cmpt_key << " which was not found in registry" << std::endl;
            }

            // Remove from registry component typed list
            if (Registry::HasItem(type_key)) {
                Registry::RemoveItem(type_key);
            } else {
                KRATOS_ERROR << "Trying ro remove: " << type_key << " which was not found in registry" << std::endl;
            }
        }

        // Finally, remove the entry all together
        Registry::RemoveItem(path);
    }
}

void KratosApplication::DeregisterCommonComponents()
{
    KRATOS_INFO("") << "Deregistering " << mApplicationName << std::endl;

    DeregisterComponent<Geometry<Node>>("geometries");
    DeregisterComponent<Element>("elements");
    DeregisterComponent<Condition>("conditions");
    DeregisterComponent<MasterSlaveConstraint>("constraints");
    DeregisterComponent<Modeler>("modelers");
    DeregisterComponent<ConstitutiveLaw>("constitutive_laws");
}

void KratosApplication::DeregisterApplication() {
    DeregisterMappers();
    // DeregisterLinearSolvers();
    // DeregisterPreconditioners();
}

void KratosApplication::DeregisterMappers() {
    // Unload the mpi branch first to avoid having a special case later
    const std::string mpi_path = "mappers."+mApplicationName+".mpi";
    if (Registry::HasItem(mpi_path)) {
        auto& r_mappers = Registry::GetItem(mpi_path);
        // Iterate over items at path. For each item, remove it from the mappers.all.mpi branch too
        for (auto i_key = r_mappers.KeyConstBegin(); i_key != r_mappers.KeyConstEnd(); ++i_key) {
            Registry::RemoveItem("mappers.all.mpi."+*i_key);
        }
        Registry::RemoveItem(mpi_path);
    }

    const std::string path = "mappers."+mApplicationName;
    if (Registry::HasItem(path)) {
        auto& r_mappers = Registry::GetItem(path);
        // Iterate over items at path. For each item, remove it from the mappers.all branch too
        for (auto i_key = r_mappers.KeyConstBegin(); i_key != r_mappers.KeyConstEnd(); ++i_key) {
            Registry::RemoveItem("mappers.all."+*i_key);
        }
        Registry::RemoveItem(path);
    }
}

}  // namespace Kratos.