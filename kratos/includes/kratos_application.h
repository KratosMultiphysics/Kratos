//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

#if !defined(KRATOS_KRATOS_APPLICATION_H_INCLUDED)
#define KRATOS_KRATOS_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/kratos_components.h"
#include "includes/element.h"
#include "elements/mesh_element.h"
#include "elements/distance_calculation_element_simplex.h"
#include "elements/levelset_convection_element_simplex.h"
#include "includes/condition.h"
#include "conditions/mesh_condition.h"
#include "includes/periodic_condition.h"
#include "utilities/quaternion.h"
#include "includes/master_slave_constraint.h"
#include "includes/linear_master_slave_constraint.h"
#include "includes/geometrical_object.h"

/* Geometries definition */
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
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"

namespace Kratos {
///@name Kratos Classes
///@{
 
/** 
 * @class KratosApplication
 * @brief This class defines the interface with kernel for all applications in Kratos.
 * @details The application class defines the interface necessary for providing the information needed by Kernel in order to configure the whole sistem correctly.
 * @ingroup KratosCore
 * @author Pooyan Dadvand
 * @author Riccardo Rossi
 */
class KRATOS_API(KRATOS_CORE) KratosApplication {
   public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;
       
    /// Pointer definition of KratosApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit KratosApplication(const std::string ApplicationName);

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
          mpArray1DVariableComponents(rOther.mpArray1DVariableComponents),
          mpArray1D4VariableComponents(rOther.mpArray1D4VariableComponents),
          mpArray1D6VariableComponents(rOther.mpArray1D6VariableComponents),
          mpArray1D9VariableComponents(rOther.mpArray1D9VariableComponents),
          mpGeometries(rOther.mpGeometries),
          mpElements(rOther.mpElements),
          mpConditions(rOther.mpConditions),
          mpMasterSlaveConstraints(rOther.mpMasterSlaveConstraints) {}

    /// Destructor.
    virtual ~KratosApplication() {}

    ///@}
    ///@name Operations
    ///@{

    virtual void Register()
    {
        RegisterKratosCore();
    }

    void RegisterKratosCore();

    ///////////////////////////////////////////////////////////////////
    void RegisterVariables();  // This contains the whole list of common variables in the Kratos Core
    void RegisterDeprecatedVariables();           //TODO: remove, this variables should not be there
    void RegisterC2CVariables();                  //TODO: move to application
    void RegisterCFDVariables();                  //TODO: move to application
    void RegisterALEVariables();                  //TODO: move to application
    void RegisterMappingVariables();              //TODO: move to application
    void RegisterDEMVariables();                  //TODO: move to application
    void RegisterFSIVariables();                  //TODO: move to application
    void RegisterMATVariables();                  //TODO: move to application
    void RegisterLegacyStructuralAppVariables();  //TODO: move to application
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

    KratosComponents<VariableComponent<VectorComponentAdaptor< array_1d<double, 3> > > >::ComponentsContainerType& GetComponents(
        VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > const& rComponentType) {
        return *mpArray1DVariableComponents;
    }

    KratosComponents<VariableComponent<VectorComponentAdaptor< array_1d<double, 4> > > >::ComponentsContainerType& GetComponents(
        VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > const& rComponentType) {
        return *mpArray1D4VariableComponents;
    }

    KratosComponents<VariableComponent<VectorComponentAdaptor< array_1d<double, 6> > > >::ComponentsContainerType& GetComponents(
        VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > const& rComponentType) {
        return *mpArray1D6VariableComponents;
    }

    KratosComponents<VariableComponent<VectorComponentAdaptor< array_1d<double, 9> > > >::ComponentsContainerType& GetComponents(
        VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > const& rComponentType) {
        return *mpArray1D9VariableComponents;
    }

    KratosComponents<VariableData>::ComponentsContainerType& GetVariables() {
        return *mpVariableData;
    }

    KratosComponents<Geometry<Node<3>>>::ComponentsContainerType& GetGeometries() {
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

    void SetComponents(KratosComponents<Geometry<Node<3>>>::ComponentsContainerType const& GeometryComponents)
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

        KratosComponents<Geometry<Node<3>>>().PrintData(rOStream);
        
        rOStream << "Elements:" << std::endl;

        KratosComponents<Element>().PrintData(rOStream);

        rOStream << std::endl;

        rOStream << "Conditions:" << std::endl;

        KratosComponents<Condition>().PrintData(rOStream);

        rOStream << std::endl;

        rOStream << "MasterSlaveConstraints:" << std::endl;

        KratosComponents<MasterSlaveConstraint>().PrintData(rOStream);
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
    //Tetrahedra:
    const Tetrahedra3D4<NodeType> mTetrahedra3D4Prototype = Tetrahedra3D4<NodeType>( GeometryType::PointsArrayType(4));
    const Tetrahedra3D10<NodeType> mTetrahedra3D10Prototype = Tetrahedra3D10<NodeType>( GeometryType::PointsArrayType(10));
    //Prisms:
    const Prism3D6<NodeType> mPrism3D6Prototype = Prism3D6<NodeType>( GeometryType::PointsArrayType(6));
    const Prism3D15<NodeType> mPrism3D15Prototype = Prism3D15<NodeType>( GeometryType::PointsArrayType(15));
    //Hexahedra:
    const Hexahedra3D8<NodeType> mHexahedra3D8Prototype = Hexahedra3D8<NodeType>( GeometryType::PointsArrayType(8));
    const Hexahedra3D20<NodeType> mHexahedra3D20Prototype = Hexahedra3D20<NodeType>( GeometryType::PointsArrayType(20));
    const Hexahedra3D27<NodeType> mHexahedra3D27Prototype = Hexahedra3D27<NodeType>( GeometryType::PointsArrayType(27));

    // General conditions must be defined

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

    // Master-Slave base constraint
    const MasterSlaveConstraint mMasterSlaveConstraint;
    const LinearMasterSlaveConstraint mLinearMasterSlaveConstraint;

    // Periodic Condition
    const PeriodicCondition mPeriodicCondition;
    const PeriodicCondition mPeriodicConditionEdge;
    const PeriodicCondition mPeriodicConditionCorner;

    // General elements must be defined
    const MeshElement mElement;
    const MeshElement mElement2D2N;
    const MeshElement mElement2D3N;
    const MeshElement mElement2D6N;
    const MeshElement mElement2D4N;
    const MeshElement mElement2D8N;
    const MeshElement mElement2D9N;

    const MeshElement mElement3D2N;
    const MeshElement mElement3D3N;
    const MeshElement mElement3D4N;
    const MeshElement mElement3D6N;
    const MeshElement mElement3D8N;
    const MeshElement mElement3D10N;
    const MeshElement mElement3D15N;
    const MeshElement mElement3D20N;
    const MeshElement mElement3D27N;

    const DistanceCalculationElementSimplex<2> mDistanceCalculationElementSimplex2D3N;
    const DistanceCalculationElementSimplex<3> mDistanceCalculationElementSimplex3D4N;

    const LevelSetConvectionElementSimplex<2,3> mLevelSetConvectionElementSimplex2D3N;
    const LevelSetConvectionElementSimplex<3,4> mLevelSetConvectionElementSimplex3D4N;

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

    KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >::ComponentsContainerType* mpArray1DVariableComponents;

    KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 4> > > >::ComponentsContainerType* mpArray1D4VariableComponents;

    KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 6> > > >::ComponentsContainerType* mpArray1D6VariableComponents;

    KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 9> > > >::ComponentsContainerType* mpArray1D9VariableComponents;

    KratosComponents<Geometry<Node<3>>>::ComponentsContainerType* mpGeometries;
    
    KratosComponents<Element>::ComponentsContainerType* mpElements;

    KratosComponents<Condition>::ComponentsContainerType* mpConditions;

    KratosComponents<MasterSlaveConstraint>::ComponentsContainerType* mpMasterSlaveConstraints;

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

#endif  // KRATOS_KRATOS_APPLICATION_H_INCLUDED  defined
