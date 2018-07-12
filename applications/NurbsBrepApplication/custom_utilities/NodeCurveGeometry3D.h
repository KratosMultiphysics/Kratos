#if !defined(KRATOS_NODE_CURVE_GEOMETRY_3D_H_INCLUDED)
#define KRATOS_NODE_CURVE_GEOMETRY_3D_H_INCLUDED

// System includes

// External includes
#include <ANurbs/Core>

// Project includes
#include "includes/define.h"
#include "includes/Node.h"
#include "includes/variables.h"

namespace Kratos {

class NodeCurveGeometry3D
: public ANurbs::CurveGeometryBase<double, Kratos::array_1d<double, 3>>
{
protected:
    using NodePointer = typename Node<3>::Pointer;

public:
    using CurveGeometryBaseType = ANurbs::CurveGeometryBase<double,
        Kratos::array_1d<double, 3>>;
    using typename CurveGeometryBaseType::KnotsType;
    using typename CurveGeometryBaseType::ScalarType;
    using typename CurveGeometryBaseType::VectorType;
 
protected:
    std::vector<NodePointer> mNodes;
 
public:
    NodeCurveGeometry3D(
        const int Degree,
        const int NumberOfNodes)
        : CurveGeometryBaseType(Degree, NumberOfNodes)
        , mNodes(NumberOfNodes)
    {
    }
 
    NodePointer
    Node(
        const int Index
    ) const
    {
        return mNodes[Index];
    }
 
    void
    SetNode(
        const int Index,
        NodePointer Value
    )
    {
        mNodes[Index] = Value;
    }
 
    VectorType
    Pole(
        const int Index) const override
    {
        auto& node = *Node(Index);
 
        VectorType pole;
        pole[0] = node[0];
        pole[1] = node[1];
        pole[2] = node[2];
 
        return pole;
    }
 
    void
    SetPole(
        const int Index,
        const VectorType& Value) override
    {
        auto& node = *Node(Index);
 
        node[0] = Value[0];
        node[1] = Value[1];
        node[2] = Value[2];
    }
 
    bool
    IsRational() const override
    {
        return true;
    }
 
    ScalarType
    Weight(
        const int Index) const override
    {
        auto node = Node(Index);
 
        return node->GetValue(Kratos::NURBS_CONTROLPOINT_WEIGHT);
    }

    void
    SetWeight(
        const int Index,
        const ScalarType Value) override
    {
        auto node = Node(Index);
 
        node->SetValue(Kratos::NURBS_CONTROLPOINT_WEIGHT, Value);
    }
    
    template <typename TDataType, typename TVariableType = Variable<TDataType>>
    TDataType
    ValueAt(
        const TVariableType& Variable,
        const double T)
    {
        return EvaluateAt<TDataType>([&](int i) -> TDataType {
            return Node(i)->GetValue(Variable);
        }, T);
    }
    
    template <typename TDataType, typename TVariableType = Variable<TDataType>>
    std::vector<TDataType>
    ValueAt(
        const TVariableType& Variable,
        const double T,
        const int Order)
    {
        return EvaluateAt<TDataType>([&](int i) -> TDataType {
            return Node(i)->GetValue(Variable);
        }, T, Order);
    }
};
 
} // namespace Kratos

#endif // !defined(KRATOS_NODE_CURVE_GEOMETRY_3D_H_INCLUDED)
