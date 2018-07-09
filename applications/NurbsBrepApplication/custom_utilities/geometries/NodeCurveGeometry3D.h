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
    std::vector<NodePointer> m_nodes;

public:
    NodeCurveGeometry3D(
        const int& degree,
        const int& nbNodes)
        : CurveGeometryBaseType(degree, nbNodes)
        , m_nodes(nbNodes)
    {
    }

    NodePointer
    Node(
        const int& index
    ) const
    {
        return m_nodes[index];
    }

    void
    SetNode(
        const int& index,
        NodePointer value
    )
    {
        m_nodes[index] = value;
    }

    VectorType
    Pole(
        const int& index) const override
    {
        auto node = Node(index);

        VectorType pole;
        pole[0] = node->X();
        pole[1] = node->Y();
        pole[2] = node->Z();

        return pole;
    }

    void
    SetPole(
        const int& index,
        const VectorType& value) override
    {
        auto node = Node(index);

        node->X() = value.X();
        node->Y() = value.Y();
        node->Z() = value.Z();
    }

    bool
    IsRational() const override
    {
        return true;
    }

    ScalarType
    Weight(
        const int& index) const override
    {
        auto node = Node(index);

        return node->GetValue(Kratos::INTEGRATION_WEIGHT); // FIXME use WEIGHT
    }

    void
    SetWeight(
        const int& index,
        const ScalarType& value) override
    {
        auto node = Node(index);

        node->SetValue(Kratos::INTEGRATION_WEIGHT, value); // FIXME use WEIGHT
    }
    
    template <typename TDataType>
    TDataType
    ValueAt(
        const Variable<TDataType>& variable,
        const double& t)
    {
        return EvaluateAt<TDataType>([&](int i) -> TDataType {
            return Node(i)->GetValue(variable);
        }, t);
    }
    
    template <typename TDataType>
    std::vector<TDataType>
    ValueAt2( // FIXME use pybind overloading
        const Variable<TDataType>& variable,
        const double& t,
        const int& order)
    {
        return EvaluateAt<TDataType>([&](int i) -> TDataType {
            return Node(i)->GetValue(variable);
        }, t, order);
    }
};

} // namespace Kratos
