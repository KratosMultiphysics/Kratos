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
    ANurbs::Grid<NodePointer> mNodes;

public:
    NodeSurfaceGeometry3D(
        const int& DegreeU,
        const int& DegreeV,
        const int& NbNodesU,
        const int& NbNodesV)
        : SurfaceGeometryBaseType(DegreeU, DegreeV, NbNodesU, NbNodesV)
        , mNodes(NbNodesU, NbNodesV)
    {
    }

    NodePointer
    Node(
        const int& IndexU,
        const int& IndexV) const
    {
        return mNodes(IndexU, IndexV);
    }

    void
    SetNode(
        const int& IndexU,
        const int& IndexV,
        NodePointer Value)
    {
        mNodes(IndexU, IndexV) = Value;
    }

    VectorType
    Pole(
        const int& IndexU,
        const int& IndexV) const override
    {
        auto node = Node(IndexU, IndexV);

        VectorType pole;
        pole[0] = node->X();
        pole[1] = node->Y();
        pole[2] = node->Z();

        return pole;
    }

    void
    SetPole(
        const int& IndexU,
        const int& IndexV,
        const VectorType& Value) override
    {
        auto node = Node(IndexU, IndexV);

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
        const int& IndexU,
        const int& IndexV) const override
    {
        auto node = Node(IndexU, IndexV);

        return node->GetValue(Kratos::INTEGRATION_WEIGHT); // FIXME use WEIGHT
    }

    void
    SetWeight(
        const int& IndexU,
        const int& IndexV,
        const ScalarType& Value) override
    {
        auto node = Node(IndexU, IndexV);

        node->SetValue(Kratos::INTEGRATION_WEIGHT, Value); // FIXME use WEIGHT
    }
    
    template <typename TDataType>
    TDataType
    ValueAt(
        const Variable<TDataType>& Variable,
        const double& U,
        const double& V)
    {
        return EvaluateAt<TDataType>([&](int i, int j) -> TDataType {
            return Node(i, j)->GetValue(Variable);
        }, U, V);
    }
    
    template <typename TDataType>
    std::vector<TDataType>
    ValueAt2( // FIXME use pybind overloading
        const Variable<TDataType>& Variable,
        const double& U,
        const double& V,
        const int& Order)
    {
        return EvaluateAt<TDataType>([&](int i, int j) -> TDataType {
            return Node(i, j)->GetValue(Variable);
        }, U, V, Order);
    }
};

} // namespace Kratos
