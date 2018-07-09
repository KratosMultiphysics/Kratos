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

class NodeSurfaceGeometry3D
    : public ANurbs::SurfaceGeometryBase<double, ANurbs::Point3D>
{
protected:
    using typename NodePointer = typename Node<3>::Pointer;

public:
    using SurfaceGeometryBaseType = SurfaceGeometryBase<double, ANurbs::Point3D>;
    using typename SurfaceGeometryBaseType::KnotsType;
    using typename SurfaceGeometryBaseType::ScalarType;
    using typename SurfaceGeometryBaseType::VectorType;

protected:
    ANurbs::Grid<NodePointer> m_nodes;

public:
    NodeSurfaceGeometry3D(
        const int& degreeU,
        const int& degreeV,
        const int& nbNodesU,
        const int& nbNodesV)
        : SurfaceGeometryBaseType(degreeU, degreeV, nbNodesU, nbNodesV)
        , m_nodes(nbNodesU, nbNodesV)
    {
    }

    NodePointer
    Node(
        const int& indexU,
        const int& indexV) const
    {
        return m_nodes(indexU, indexV);
    }

    void
    SetNode(
        const int& indexU,
        const int& indexV,
        NodePointer value)
    {
        m_nodes(indexU, indexV) = value;
    }

    VectorType
    Pole(
        const int& indexU,
        const int& indexV) const override
    {
        auto node = Node(indexU, indexV);

        VectorType pole;
        pole[0] = node->X();
        pole[1] = node->Y();
        pole[2] = node->Z();

        return pole;
    }

    void
    SetPole(
        const int& indexU,
        const int& indexV,
        const VectorType& value) override
    {
        auto node = Node(indexU, indexV);

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
        const int& indexU,
        const int& indexV) const override
    {
        auto node = Node(indexU, indexV);

        return node->GetValue(Kratos::INTEGRATION_WEIGHT); // FIXME use WEIGHT
    }

    void
    SetWeight(
        const int& indexU,
        const int& indexV,
        const ScalarType& value) override
    {
        auto node = Node(indexU, indexV);

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
