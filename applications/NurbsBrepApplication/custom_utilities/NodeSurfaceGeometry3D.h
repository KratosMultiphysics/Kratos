#if !defined(KRATOS_NODE_SURFACE_GEOMETRY_3D_H_INCLUDED)
#define KRATOS_NODE_SURFACE_GEOMETRY_3D_H_INCLUDED

// System includes

// External includes
#include <ANurbs/Core>

// Project includes
#include "includes/define.h"
#include "includes/Node.h"
#include "includes/variables.h"

namespace Kratos {

class NodeSurfaceGeometry3D
    : public ANurbs::SurfaceGeometryBase<double, Kratos::array_1d<double, 3>>
{
protected:
    using NodePointer = typename Node<3>::Pointer;

public:
    using SurfaceGeometryBaseType = SurfaceGeometryBase<double,
        Kratos::array_1d<double, 3>>;
    using typename SurfaceGeometryBaseType::KnotsType;
    using typename SurfaceGeometryBaseType::ScalarType;
    using typename SurfaceGeometryBaseType::VectorType;
 
protected:
    ANurbs::Grid<NodePointer> mNodes;
 
public:
    NodeSurfaceGeometry3D(
        const int DegreeU,
        const int DegreeV,
        const int NumberOfNodesU,
        const int NumberOfNodesV)
        : SurfaceGeometryBaseType(DegreeU, DegreeV, NumberOfNodesU,
            NumberOfNodesV)
        , mNodes(NumberOfNodesU, NumberOfNodesV)
    {
    }

    NodePointer
    Node(
        const int IndexU,
        const int IndexV) const
    {
        return mNodes(IndexU, IndexV);
    }
 
    void
    SetNode(
        const int IndexU,
        const int IndexV,
        NodePointer Value)
    {
        mNodes.SetValue(IndexU, IndexV, Value);
    }
 
    VectorType
    Pole(
        const int IndexU,
        const int IndexV) const override
    {
        auto& node = *Node(IndexU, IndexV);
 
        VectorType pole;
        pole[0] = node[0];
        pole[1] = node[1];
        pole[2] = node[2];
 
        return pole;
    }

    void
    SetPole(
        const int IndexU,
        const int IndexV,
        const VectorType& Value) override
    {
        auto& node = *Node(IndexU, IndexV);
 
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
        const int IndexU,
        const int IndexV) const override
    {
        auto node = Node(IndexU, IndexV);
 
        return node->GetValue(Kratos::NURBS_CONTROLPOINT_WEIGHT);
    }

    void
    SetWeight(
        const int IndexU,
        const int IndexV,
        const ScalarType Value) override
    {
        auto node = Node(IndexU, IndexV);
 
        node->SetValue(Kratos::NURBS_CONTROLPOINT_WEIGHT, Value);
    }
    
    template <typename TDataType, typename TVariableType = Variable<TDataType>>
    TDataType
    ValueAt(
        const TVariableType& Variable,
        const double U,
        const double V)
    {
        return EvaluateAt<TDataType>([&](int i, int j) -> TDataType {
            return Node(i, j)->GetValue(Variable);
        }, U, V);
    }
    
    template <typename TDataType, typename TVariableType = Variable<TDataType>>
    std::vector<TDataType>
    ValueAt(
        const TVariableType& Variable,
        const double U,
        const double V,
        const int Order)
    {
        return EvaluateAt<TDataType>([&](int i, int j) -> TDataType {
            return Node(i, j)->GetValue(Variable);
        }, U, V, Order);
    }
};

}

#endif // !defined(KRATOS_NODE_SURFACE_GEOMETRY_3D_H_INCLUDED)
