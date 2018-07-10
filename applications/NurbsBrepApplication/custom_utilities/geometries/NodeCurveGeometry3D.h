/*
//  KRATOS _____________ 
//        /  _/ ____/   |
//        / // / __/ /| |
//      _/ // /_/ / ___ |
//     /___/\____/_/  |_| StructuralMechanicsApplication
//
//  Author: Thomas Oberbichler
*/

#pragma once

#include <ANurbs/Core>


namespace Kratos {

using NodalCurveControlPoints = std::vector<typename Node<3>::Pointer>;

class NodeCurveGeometry3D
: public ANurbs::CurveGeometryBase<double, ANurbs::Point3D>
{
protected:
    using typename NodePointer = typename Node<3>::Pointer;

public:
    using CurveGeometryBaseType = ANurbs::CurveGeometryBase<double, ANurbs::Point3D>;
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
        auto node = Node(Index);

        VectorType pole;
        pole[0] = node->X();
        pole[1] = node->Y();
        pole[2] = node->Z();

        return pole;
    }

    void
    SetPole(
        const int Index,
        const VectorType& Value) override
    {
        auto node = Node(Index);

        node->X() = Value.X();
        node->Y() = Value.Y();
        node->Z() = Value.Z();
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

        return node->GetValue(Kratos::INTEGRATION_WEIGHT); // FIXME use WEIGHT
    }

    void
    SetWeight(
        const int Index,
        const ScalarType Value) override
    {
        auto node = Node(Index);

        node->SetValue(Kratos::INTEGRATION_WEIGHT, Value); // FIXME use WEIGHT
    }
    
    template <typename TDataType>
    TDataType
    ValueAt(
        const Variable<TDataType>& Variable,
        const double& T)
    {
        return EvaluateAt<TDataType>([&](int i) -> TDataType {
            return Node(i)->GetValue(Variable);
        }, T);
    }
    
    template <typename TDataType>
    std::vector<TDataType>
    ValueAt2( // FIXME use pybind overloading
        const Variable<TDataType>& Variable,
        const double& T,
        const int Order)
    {
        return EvaluateAt<TDataType>([&](int i) -> TDataType {
            return Node(i)->GetValue(Variable);
        }, T, Order);
    }
};

} // namespace Kratos
