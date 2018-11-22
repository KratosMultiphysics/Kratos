/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
*/

#if !defined(KRATOS_NODE_SURFACE_GEOMETRY_3D_H_INCLUDED)
#define KRATOS_NODE_SURFACE_GEOMETRY_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "anurbs.h"

namespace Kratos {

/// Spatial NURBS surface geometry with Kratos Nodes as control points.
/** Unlike the ANurbs::SurfaceGeometry, this class does not use static Points as
 *  control points but Finite Element Nodes. That means that the Geometry
 *  changes whenever the Nodes are moving.
 */
class NodeSurfaceGeometry3D
    : public ANurbs::SurfaceGeometryBase<Kratos::array_1d<double, 3>>
{
protected:
    using NodePointer = typename Node<3>::Pointer;

public:
    using NodeType = Node<3>;
    using SurfaceGeometryBaseType = SurfaceGeometryBase<
        Kratos::array_1d<double, 3>>;
    using typename SurfaceGeometryBaseType::KnotsType;
    using typename SurfaceGeometryBaseType::ScalarType;
    using typename SurfaceGeometryBaseType::VectorType;

protected:
    ANurbs::Grid<NodePointer> mNodes;

public:
    /** Creates a new NodeSurfaceGeometry3D.
     *
     *  @param DegreeU Degree in u direction
     *  @param DegreeV Degree in v direction
     *  @param NumberOfNodesU Number of nodes in u direction
     *  @param NumberOfNodesV Number of nodes in v direction
     */
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

    /** Gets the Kratos node at a given index.
     * 
     * @param IndexU Index in u direction
     * @param IndexV Index in v direction
     * 
     * @return Kratos node at the given index.
     */
    NodePointer GetNode(
        const int IndexU,
        const int IndexV) const
    {
        KRATOS_DEBUG_ERROR_IF(IndexU < 0 || NbPolesU() <= IndexU) <<
            "IndexU out of range" << std::endl;
        KRATOS_DEBUG_ERROR_IF(IndexV < 0 || NbPolesV() <= IndexV) <<
            "IndexV out of range" << std::endl;

        return mNodes(IndexU, IndexV);
    }

    /** Sets the Kratos node at a given index.
     * 
     * @param IndexU Index in u direction
     * @param IndexV Index in v direction
     */
    void SetNode(
        const int IndexU,
        const int IndexV,
        NodePointer Value)
    {
        KRATOS_DEBUG_ERROR_IF(IndexU < 0 || NbPolesU() <= IndexU) <<
            "IndexU out of range" << std::endl;
        KRATOS_DEBUG_ERROR_IF(IndexV < 0 || NbPolesV() <= IndexV) <<
            "IndexV out of range" << std::endl;

        mNodes.SetValue(IndexU, IndexV, Value);
    }

    /** Gets the location of the Kratos node at a given index.
     * 
     * @param IndexU Index in u direction
     * @param IndexV Index in v direction
     * 
     * @return Location of the Kratos node at the given index.
     */
    VectorType Pole(
        const int IndexU,
        const int IndexV) const override
    {
        const NodeType& node = *GetNode(IndexU, IndexV);

        VectorType pole;
        for (std::size_t i = 0; i < 3; i++) {
            pole[i] = node[i];
        }

        return pole;
    }

    /** Sets the location of the Kratos node at a given index.
     * 
     * @param IndexU Index in u direction
     * @param IndexV Index in v direction
     * @param Value New location of the Kratos node
     */
    void SetPole(
        const int IndexU,
        const int IndexV,
        const VectorType& Value) override
    {
        NodeType& node = *GetNode(IndexU, IndexV);

        for (std::size_t i = 0; i < 3; i++) {
            node[i] = Value[i];
        }
    }

    /** Gets a value indicating whether or not the NURBS surface is rational.
     * 
     * @return True whether the surface is rational, otherwise false.
     */
    bool IsRational() const override
    {
        return true;
    }

    /** Gets the weight of the Kratos node at a given index.
     * 
     * @param IndexU Index in u direction
     * @param IndexV Index in v direction
     * 
     * @return Weight of the Kratos node at the given index.
     */
    ScalarType Weight(
        const int IndexU,
        const int IndexV) const override
    {
        const NodeType& node = *GetNode(IndexU, IndexV);

        if (node.Has(Kratos::NURBS_CONTROL_POINT_WEIGHT)) {
            return node.GetValue(Kratos::NURBS_CONTROL_POINT_WEIGHT);
        } else {
            return 1;
        }
    }

    /** Sets the weight of the Kratos node at a given index.
     * 
     * @param IndexU Index in u direction
     * @param IndexV Index in v direction
     * @param Value New weight of the Kratos node
     */
    void SetWeight(
        const int IndexU,
        const int IndexV,
        const ScalarType Value) override
    {
        NodeType& node = *GetNode(IndexU, IndexV);

        node.SetValue(Kratos::NURBS_CONTROL_POINT_WEIGHT, Value);
    }

    /** Gets the value of a nodal Kratos variable on a point at the surface.
     * 
     * @param Variable Kratos variable
     * @param U Surface parameter in u direction
     * @param V Surface parameter in v direction
     * 
     * @return The value of the variable at the given surface point.
     */
    template <typename TDataType, typename TVariableType = Variable<TDataType>>
    TDataType ValueAt(
        const TVariableType& Variable,
        const double U,
        const double V) const
    {
        return EvaluateAt<TDataType>([&](int i, int j) -> TDataType {
            return GetNode(i, j)->GetValue(Variable);
        }, U, V);
    }

    /** Gets the derivatives of a nodal Kratos variable on a point at the
     * surface.
     * 
     * @param Variable Kratos variable
     * @param U Surface parameter in u direction
     * @param V Surface parameter in v direction
     * @param Order Order of the highest derivative to compute
     * 
     * @return The value and the derivatives of the variable at the given
     * surface point.
     */
    template <typename TDataType, typename TVariableType = Variable<TDataType>>
    std::vector<TDataType> ValueAt(
        const TVariableType& Variable,
        const double U,
        const double V,
        const int Order) const
    {
        return EvaluateAt<TDataType>([&](int i, int j) -> TDataType {
            return GetNode(i, j)->GetValue(Variable);
        }, U, V, Order);
    }
};

}

#endif // !defined(KRATOS_NODE_SURFACE_GEOMETRY_3D_H_INCLUDED)
