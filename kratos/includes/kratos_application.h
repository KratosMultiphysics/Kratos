//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

#pragma once

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/kratos_components.h"
#include "includes/geometrical_object.h"
#include "includes/periodic_condition.h"
#include "includes/master_slave_constraint.h"
#include "input_output/logger.h"
#include "utilities/quaternion.h"
#include "constraints/linear_master_slave_constraint.h"

// Geometries definition
#include "geometries/register_kratos_components_for_geometry.h"
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
#include "geometries/pyramid_3d_5.h"
#include "geometries/pyramid_3d_13.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"
#include "geometries/quadrature_point_geometry.h"
#include "geometries/coupling_geometry.h"

// Elements
#include "elements/mesh_element.h"
#include "elements/distance_calculation_element_simplex.h"
#include "elements/edge_based_gradient_recovery_element.h"
#include "elements/levelset_convection_element_simplex.h"
#include "elements/levelset_convection_element_simplex_algebraic_stabilization.h"

// Conditions
#include "conditions/mesh_condition.h"

// Modelers
#include "modeler/modeler.h"
#include "modeler/cad_io_modeler.h"
#include "modeler/cad_tessellation_modeler.h"
#include "modeler/serial_model_part_combinator_modeler.h"
#include "modeler/combine_model_part_modeler.h"
#include "modeler/connectivity_preserve_modeler.h"
#include "modeler/voxel_mesh_generator_modeler.h"
#include "modeler/clean_up_problematic_triangles_modeler.h"

namespace Kratos {
///@name Kratos Classes
///@{

/**
 * @class KratosApplication
 * @brief This class defines the interface with kernel for all applications in Kratos.
 * @details The application class defines the interface necessary for providing the information needed by Kernel in order to configure the whole system correctly.
 * @ingroup KratosCore
 * @author Pooyan Dadvand
 * @author Riccardo Rossi
 */
class KRATOS_API(KRATOS_CORE) KratosApplication {
   public:
    ///@name Type Definitions
    ///@{

    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;

    /// Pointer definition of KratosApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit KratosApplication(const std::string& ApplicationName);

    KratosApplication() = delete;

    /// Copy constructor.
    KratosApplication(KratosApplication const& rOther)
        : mpVariableData(rOther.mpVariableData),
          mpIntVariables(rOther.mpIntVariables),
          mpUnsignedIntVariables(rOther.mpUnsignedIntVariables),
          mpDoubleVariables(rOther.mpDoubleVariables),
          mpArray1DVariables(rOther.mpArray1DVariables),
          mpArray1D4Variables(rOther.mpArray1D4Variables),
          mpArray1D6Variables(rOther.mpArray1D6Variables),
          mpArray1D9Variables(rOther.mpArray1D9Variables),
          mpVectorVariables(rOther.mpVectorVariables),
          mpMatrixVariables(rOther.mpMatrixVariables),
          mpGeometries(rOther.mpGeometries),
          mpElements(rOther.mpElements),
          mpConditions(rOther.mpConditions),
          mpMasterSlaveConstraints(rOther.mpMasterSlaveConstraints),
          mpModelers(rOther.mpModelers) {}

    /// Destructor.
    virtual ~KratosApplication()
    {
        // This must be commented until tests have been fixed.
        DeregisterCommonComponents();
        DeregisterApplication();
    }

    ///@}
    ///@name Operations
    ///@{

    virtual void Register()
    {
        RegisterKratosCore();
    }

    void RegisterKratosCore();

    template<class TComponentsContainer>
    void DeregisterComponent(std::string const & rComponentName);

    /**
     * @brief This method is used to unregister common components of the application.
     * @details This method is used to unregister common components of the application.
     * The list of unregistered components are the ones exposed in the common KratosComponents interface:
     * - Geometries
     * - Elements
     * - Conditions
     * - MasterSlaveConstraints
     * - Modelers
     * - ConstitutiveLaws
     */
    void DeregisterCommonComponents();

    /**
     * @brief This method is used to unregister specific application components.
     * @details This method is used to unregister specific application components.
     */
    virtual void DeregisterApplication();

    /**
     * @brief This method is used to unregister mappers registered by the application..
     * @details This method is used to unregister mappers registered by the application..
     */
    void DeregisterMappers();

    ///////////////////////////////////////////////////////////////////
    void RegisterVariables();  // This contains the whole list of common variables in the Kratos Core
    void RegisterDeprecatedVariables();           //TODO: remove, this variables should not be there
    void RegisterCFDVariables();                  //TODO: move to application
    void RegisterALEVariables();                  //TODO: move to application
    void RegisterMappingVariables();              //TODO: move to application
    void RegisterDEMVariables();                  //TODO: move to application
    void RegisterFSIVariables();                  //TODO: move to application
    void RegisterMATVariables();                  //TODO: move to application
    void RegisterGlobalPointerVariables();

    const std::string& Name() const { return mApplicationName; }

    ///@}
    ///@name Access
    ///@{

    // I have to see why the above version is not working for multi thread ...
    // Anyway its working with these functions.Pooyan.
    KratosComponents<Variable<int> >::ComponentsContainerType& GetComponents(
        Variable<int> const& rComponentType) {
        return *mpIntVariables;
    }

    KratosComponents<Variable<unsigned int> >::ComponentsContainerType&
    GetComponents(Variable<unsigned int> const& rComponentType) {
        return *mpUnsignedIntVariables;
    }

    KratosComponents<Variable<double> >::ComponentsContainerType& GetComponents(
        Variable<double> const& rComponentType) {
        return *mpDoubleVariables;
    }

    KratosComponents<Variable<array_1d<double, 3> > >::ComponentsContainerType&
    GetComponents(Variable<array_1d<double, 3> > const& rComponentType) {
        return *mpArray1DVariables;
    }

    KratosComponents<Variable<array_1d<double, 4> > >::ComponentsContainerType&
    GetComponents(Variable<array_1d<double, 4> > const& rComponentType) {
        return *mpArray1D4Variables;
    }

    KratosComponents<Variable<array_1d<double, 6> > >::ComponentsContainerType&
    GetComponents(Variable<array_1d<double, 6> > const& rComponentType) {
        return *mpArray1D6Variables;
    }

    KratosComponents<Variable<array_1d<double, 9> > >::ComponentsContainerType&
    GetComponents(Variable<array_1d<double, 9> > const& rComponentType) {
        return *mpArray1D9Variables;
    }

    KratosComponents<Variable<Quaternion<double> > >::ComponentsContainerType&
    GetComponents(Variable<Quaternion<double> > const& rComponentType) {
        return *mpQuaternionVariables;
    }

    KratosComponents<Variable<Vector> >::ComponentsContainerType& GetComponents(
        Variable<Vector> const& rComponentType) {
        return *mpVectorVariables;
    }

    KratosComponents<Variable<Matrix> >::ComponentsContainerType& GetComponents(
        Variable<Matrix> const& rComponentType) {
        return *mpMatrixVariables;
    }

    KratosComponents<VariableData>::ComponentsContainerType& GetVariables() {
        return *mpVariableData;
    }

    KratosComponents<Geometry<Node>>::ComponentsContainerType& GetGeometries() {
        return *mpGeometries;
    }

    KratosComponents<Element>::ComponentsContainerType& GetElements() {
        return *mpElements;
    }

    KratosComponents<Condition>::ComponentsContainerType& GetConditions() {
        return *mpConditions;
    }

    KratosComponents<MasterSlaveConstraint>::ComponentsContainerType& GetMasterSlaveConstraints() {
        return *mpMasterSlaveConstraints;
    }

    KratosComponents<Modeler>::ComponentsContainerType& GetModelers() {
        return *mpModelers;
    }

    void SetComponents(
        KratosComponents<VariableData>::ComponentsContainerType const&
            VariableDataComponents)

    {
        for (auto it = mpVariableData->begin(); it != mpVariableData->end(); it++) {
            std::string const& r_variable_name = it->second->Name();
            auto it_variable = VariableDataComponents.find(r_variable_name);
            KRATOS_ERROR_IF(it_variable == VariableDataComponents.end()) << "This variable is not registered in Kernel : " << *(it_variable->second) << std::endl;
        }
    }

    void SetComponents(KratosComponents<Geometry<Node>>::ComponentsContainerType const& GeometryComponents)
    {
        mpGeometries->insert(GeometryComponents.begin(), GeometryComponents.end());
    }

    void SetComponents(KratosComponents<Element>::ComponentsContainerType const&
            ElementComponents)

    {
        // It's better to make a loop over new components and add them if they are NOT already exist in application. Or make an ERROR for incompatibility between applications.
        mpElements->insert(ElementComponents.begin(), ElementComponents.end());
    }

    void SetComponents(KratosComponents<MasterSlaveConstraint>::ComponentsContainerType const&
            MasterSlaveConstraintComponents)

    {
        mpMasterSlaveConstraints->insert(MasterSlaveConstraintComponents.begin(), MasterSlaveConstraintComponents.end());
    }

    void SetComponents(KratosComponents<Modeler>::ComponentsContainerType const& ModelerComponents)
    {
        mpModelers->insert(ModelerComponents.begin(), ModelerComponents.end());
    }

    void SetComponents(
        KratosComponents<Condition>::ComponentsContainerType const&
            ConditionComponents)

    {
        mpConditions->insert(
            ConditionComponents.begin(), ConditionComponents.end());
    }

    Serializer::RegisteredObjectsContainerType& GetRegisteredObjects() {
        return *mpRegisteredObjects;
    }

    Serializer::RegisteredObjectsNameContainerType& GetRegisteredObjectsName() {
        return *mpRegisteredObjectsName;
    }

    ///@}

    ///@name Inquiry

    ///@{

    ///@}

    ///@name Input and output

    ///@{

    /// Turn back information as a string.

    virtual std::string Info() const

    {
        return "KratosApplication";
    }

    /// Print information about this object.

    virtual void PrintInfo(std::ostream& rOStream) const

    {
        rOStream << Info();
    }

    /// Print object's data.

    virtual void PrintData(std::ostream& rOStream) const

    {
        rOStream << "Variables:" << std::endl;

        KratosComponents<VariableData>().PrintData(rOStream);

        rOStream << std::endl;

        rOStream << "Geometries:" << std::endl;

        KratosComponents<Geometry<Node>>().PrintData(rOStream);

        rOStream << "Elements:" << std::endl;

        KratosComponents<Element>().PrintData(rOStream);

        rOStream << std::endl;

        rOStream << "Conditions:" << std::endl;

        KratosComponents<Condition>().PrintData(rOStream);

        rOStream << std::endl;

        rOStream << "MasterSlaveConstraints:" << std::endl;

        KratosComponents<MasterSlaveConstraint>().PrintData(rOStream);

        rOStream << std::endl;

        rOStream << "Modelers:" << std::endl;

        KratosComponents<Modeler>().PrintData(rOStream);
    }

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

    std::string mApplicationName;

    // General geometries must be defined
    //Points:
    const Point mPointPrototype;
    const Point2D<NodeType> mPoint2DPrototype = Point2D<NodeType>(GeometryType::PointsArrayType(1));
    const Point3D<NodeType> mPoint3DPrototype = Point3D<NodeType>(GeometryType::PointsArrayType(1));
    //Sphere
    const Sphere3D1<NodeType> mSphere3D1Prototype = Sphere3D1<NodeType>(GeometryType::PointsArrayType(1));
    //Lines:
    const Line2D2<NodeType> mLine2D2Prototype = Line2D2<NodeType>(GeometryType::PointsArrayType(2));
    const Line2D3<NodeType> mLine2D3Prototype = Line2D3<NodeType>(GeometryType::PointsArrayType(3));
    const Line3D2<NodeType> mLine3D2Prototype = Line3D2<NodeType>(GeometryType::PointsArrayType(2));
    const Line3D3<NodeType> mLine3D3Prototype = Line3D3<NodeType>(GeometryType::PointsArrayType(3));
    //Triangles:
    const Triangle2D3<NodeType> mTriangle2D3Prototype = Triangle2D3<NodeType>(GeometryType::PointsArrayType(3));
    const Triangle2D6<NodeType> mTriangle2D6Prototype = Triangle2D6<NodeType>(GeometryType::PointsArrayType(6));
    const Triangle3D3<NodeType> mTriangle3D3Prototype = Triangle3D3<NodeType>(GeometryType::PointsArrayType(3));
    const Triangle3D6<NodeType> mTriangle3D6Prototype = Triangle3D6<NodeType>( GeometryType::PointsArrayType(6));
    //Quadrilaterals:
    const Quadrilateral2D4<NodeType> mQuadrilateral2D4Prototype = Quadrilateral2D4<NodeType>( GeometryType::PointsArrayType(4));
    const Quadrilateral2D8<NodeType> mQuadrilateral2D8Prototype = Quadrilateral2D8<NodeType>( GeometryType::PointsArrayType(8));
    const Quadrilateral2D9<NodeType> mQuadrilateral2D9Prototype = Quadrilateral2D9<NodeType>( GeometryType::PointsArrayType(9));
    const Quadrilateral3D4<NodeType> mQuadrilateral3D4Prototype = Quadrilateral3D4<NodeType>( GeometryType::PointsArrayType(4));
    const Quadrilateral3D8<NodeType> mQuadrilateral3D8Prototype = Quadrilateral3D8<NodeType>( GeometryType::PointsArrayType(8));
    const Quadrilateral3D9<NodeType> mQuadrilateral3D9Prototype = Quadrilateral3D9<NodeType>( GeometryType::PointsArrayType(9));
    //Coupling
    const CouplingGeometry<NodeType> mCouplingGeometry = CouplingGeometry<NodeType>();
    //Tetrahedra:
    const Tetrahedra3D4<NodeType> mTetrahedra3D4Prototype = Tetrahedra3D4<NodeType>( GeometryType::PointsArrayType(4));
    const Tetrahedra3D10<NodeType> mTetrahedra3D10Prototype = Tetrahedra3D10<NodeType>( GeometryType::PointsArrayType(10));
    //Prisms:
    const Prism3D6<NodeType> mPrism3D6Prototype = Prism3D6<NodeType>( GeometryType::PointsArrayType(6));
    const Prism3D15<NodeType> mPrism3D15Prototype = Prism3D15<NodeType>( GeometryType::PointsArrayType(15));
    //Pyramids:
    const Pyramid3D5<NodeType> mPyramid3D5Prototype = Pyramid3D5<NodeType>( GeometryType::PointsArrayType(5));
    const Pyramid3D13<NodeType> mPyramid3D13Prototype = Pyramid3D13<NodeType>( GeometryType::PointsArrayType(13));
    //Hexahedra:
    const Hexahedra3D8<NodeType> mHexahedra3D8Prototype = Hexahedra3D8<NodeType>( GeometryType::PointsArrayType(8));
    const Hexahedra3D20<NodeType> mHexahedra3D20Prototype = Hexahedra3D20<NodeType>( GeometryType::PointsArrayType(20));
    const Hexahedra3D27<NodeType> mHexahedra3D27Prototype = Hexahedra3D27<NodeType>( GeometryType::PointsArrayType(27));
    //QuadraturePointGeometries:
    const QuadraturePointGeometry<Node,1> mQuadraturePointGeometryPoint1D = QuadraturePointGeometry<Node,1>(GeometryType::PointsArrayType(),
        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>(GeometryData::IntegrationMethod::GI_GAUSS_1, IntegrationPoint<3>(), Matrix(), Matrix()));
    const QuadraturePointGeometry<Node,2,1> mQuadraturePointGeometryPoint2D = QuadraturePointGeometry<Node,2,1>(GeometryType::PointsArrayType(),
        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>(GeometryData::IntegrationMethod::GI_GAUSS_1, IntegrationPoint<3>(), Matrix(), Matrix()));
    const QuadraturePointGeometry<Node,3,1> mQuadraturePointGeometryPoint3D = QuadraturePointGeometry<Node,3,1>(GeometryType::PointsArrayType(),
        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>(GeometryData::IntegrationMethod::GI_GAUSS_1, IntegrationPoint<3>(), Matrix(), Matrix()));
    const QuadraturePointGeometry<Node,2> mQuadraturePointGeometrySurface2D = QuadraturePointGeometry<Node,2>(GeometryType::PointsArrayType(),
        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>(GeometryData::IntegrationMethod::GI_GAUSS_1, IntegrationPoint<3>(), Matrix(), Matrix()));
    const QuadraturePointGeometry<Node,3,2> mQuadraturePointGeometrySurface3D = QuadraturePointGeometry<Node,3,2>(GeometryType::PointsArrayType(),
        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>(GeometryData::IntegrationMethod::GI_GAUSS_1, IntegrationPoint<3>(), Matrix(), Matrix()));
    const QuadraturePointGeometry<Node,3> mQuadraturePointGeometryVolume3D = QuadraturePointGeometry<Node,3>(GeometryType::PointsArrayType(),
        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>(GeometryData::IntegrationMethod::GI_GAUSS_1, IntegrationPoint<3>(), Matrix(), Matrix()));

    // General conditions must be defined
    // Generic condition
    const MeshCondition mGenericCondition;

    // Point conditions
    const MeshCondition mPointCondition2D1N;
    const MeshCondition mPointCondition3D1N;

    // Line conditions
    const MeshCondition mLineCondition2D2N;
    const MeshCondition mLineCondition2D3N;
    const MeshCondition mLineCondition3D2N;
    const MeshCondition mLineCondition3D3N;

    // Surface conditions
    const MeshCondition mSurfaceCondition3D3N;
    const MeshCondition mSurfaceCondition3D6N;
    const MeshCondition mSurfaceCondition3D4N;
    const MeshCondition mSurfaceCondition3D8N;
    const MeshCondition mSurfaceCondition3D9N;

    // Prisms conditions
    const MeshCondition mPrismCondition2D4N;
    const MeshCondition mPrismCondition3D6N;


    // Master-Slave base constraint
    const MasterSlaveConstraint mMasterSlaveConstraint;
    const LinearMasterSlaveConstraint mLinearMasterSlaveConstraint;

    // Periodic Condition
    const PeriodicCondition mPeriodicCondition;
    const PeriodicCondition mPeriodicConditionEdge;
    const PeriodicCondition mPeriodicConditionCorner;

    // General elements must be defined
    const MeshElement mGenericElement;

    // 2D elements
    const MeshElement mElement2D1N; // NOTE: This would be duplicated mPointElement2D1N
    const MeshElement mElement2D2N; // NOTE: This would be duplicated mLineElement2D2N
    const MeshElement mElement2D3N;
    const MeshElement mElement2D6N;
    const MeshElement mElement2D4N;
    const MeshElement mElement2D8N;
    const MeshElement mElement2D9N;

    // 3D elements
    const MeshElement mElement3D1N; // NOTE: This would be duplicated mPointElement3D1N
    const MeshElement mElement3D2N; // NOTE: This would be duplicated mLineElement3D2N
    const MeshElement mElement3D3N; // NOTE: This would be duplicated mSurfaceElement3D3N
    const MeshElement mElement3D4N;
    const MeshElement mElement3D5N;
    const MeshElement mElement3D6N;
    const MeshElement mElement3D8N;
    const MeshElement mElement3D10N;
    const MeshElement mElement3D13N;
    const MeshElement mElement3D15N;
    const MeshElement mElement3D20N;
    const MeshElement mElement3D27N;

    // For sake of consistency and in order to define missing geometries for elements
    // Point elements
    const MeshElement mPointElement2D1N;
    const MeshElement mPointElement3D1N;

    // Line elements
    const MeshElement mLineElement2D2N;
    const MeshElement mLineElement2D3N;
    const MeshElement mLineElement3D2N;
    const MeshElement mLineElement3D3N;

    // Surface elements
    const MeshElement mSurfaceElement3D3N;
    const MeshElement mSurfaceElement3D6N;
    const MeshElement mSurfaceElement3D4N;
    const MeshElement mSurfaceElement3D8N;
    const MeshElement mSurfaceElement3D9N;

    // Distance calculation elements
    const DistanceCalculationElementSimplex<2> mDistanceCalculationElementSimplex2D3N;
    const DistanceCalculationElementSimplex<3> mDistanceCalculationElementSimplex3D4N;

    // Edge based gradient recovery elements
    const EdgeBasedGradientRecoveryElement<2> mEdgeBasedGradientRecoveryElement2D2N;
    const EdgeBasedGradientRecoveryElement<3> mEdgeBasedGradientRecoveryElement3D2N;

    // Level set convection elements
    const LevelSetConvectionElementSimplex<2,3> mLevelSetConvectionElementSimplex2D3N;
    const LevelSetConvectionElementSimplex<3,4> mLevelSetConvectionElementSimplex3D4N;
    const LevelSetConvectionElementSimplexAlgebraicStabilization<2,3> mLevelSetConvectionElementSimplexAlgebraicStabilization2D3N;
    const LevelSetConvectionElementSimplexAlgebraicStabilization<3,4> mLevelSetConvectionElementSimplexAlgebraicStabilization3D4N;

    // Modeler
    const Modeler mModeler;
    const CadIoModeler mCadIoModeler;
#if USE_TRIANGLE_NONFREE_TPL
    const CadTessellationModeler mCadTessellationModeler;
#endif
    const SerialModelPartCombinatorModeler mSerialModelPartCombinatorModeler;
    const CombineModelPartModeler mCombineModelPartModeler;
    const ConnectivityPreserveModeler mConnectivityPreserveModeler;
    const VoxelMeshGeneratorModeler mVoxelMeshGeneratorModeler;
    const CleanUpProblematicTrianglesModeler mCleanUpProblematicTrianglesModeler;

    // Base constitutive law definition
    const ConstitutiveLaw mConstitutiveLaw;

    // KratosComponents definition
    KratosComponents<VariableData>::ComponentsContainerType* mpVariableData;

    KratosComponents<Variable<int> >::ComponentsContainerType* mpIntVariables;

    KratosComponents<Variable<unsigned int> >::ComponentsContainerType* mpUnsignedIntVariables;

    KratosComponents<Variable<double> >::ComponentsContainerType* mpDoubleVariables;

    KratosComponents<Variable<array_1d<double, 3> > >::ComponentsContainerType* mpArray1DVariables;

    KratosComponents<Variable<array_1d<double, 4> > >::ComponentsContainerType* mpArray1D4Variables;

    KratosComponents<Variable<array_1d<double, 6> > >::ComponentsContainerType* mpArray1D6Variables;

    KratosComponents<Variable<array_1d<double, 9> > >::ComponentsContainerType* mpArray1D9Variables;

    KratosComponents<Variable<Quaternion<double> > >::ComponentsContainerType* mpQuaternionVariables;

    KratosComponents<Variable<Vector> >::ComponentsContainerType* mpVectorVariables;

    KratosComponents<Variable<Matrix> >::ComponentsContainerType* mpMatrixVariables;

    KratosComponents<Geometry<Node>>::ComponentsContainerType* mpGeometries;

    KratosComponents<Element>::ComponentsContainerType* mpElements;

    KratosComponents<Condition>::ComponentsContainerType* mpConditions;

    KratosComponents<MasterSlaveConstraint>::ComponentsContainerType* mpMasterSlaveConstraints;

    KratosComponents<Modeler>::ComponentsContainerType* mpModelers;

    // Serialization
    Serializer::RegisteredObjectsContainerType* mpRegisteredObjects;

    Serializer::RegisteredObjectsNameContainerType* mpRegisteredObjectsName;

    ///@}

    ///@name Protected Operators

    ///@{

    ///@}

    ///@name Protected Operations

    ///@{

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

    KratosApplication& operator=(KratosApplication const& rOther);

    ///@}

};  // Class KratosApplication

///@}

///@name Type Definitions

///@{

///@}

///@name Input and output

///@{

/// input stream function

inline std::istream& operator>>(std::istream& rIStream,

    KratosApplication& rThis);

/// output stream function

inline std::ostream& operator<<(std::ostream& rOStream,

    const KratosApplication& rThis)

{
    rThis.PrintInfo(rOStream);

    rOStream << std::endl;

    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.