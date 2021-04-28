// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//
// ==============================================================================

#ifndef SYMMETRY_REVOLUTION_H
#define SYMMETRY_REVOLUTION_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"
#include "spaces/ublas_space.h"
#include "symmetry_base.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/

class SymmetryRevolution : public SymmetryBase
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(SymmetryRevolution);

    SymmetryRevolution(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, Parameters Settings)
    : SymmetryBase(rOriginModelPart, rDestinationModelPart, Settings)
    {
        mPoint = mSettings["point"].GetVector();
        mAxis = mSettings["axis"].GetVector();
        mAxis /= norm_2(mAxis);

        mPlaneVector1 = ZeroVector(3);
        int index = 0;
        double max_value = 0;
        for (int i=0; i<3; ++i){
            if (std::abs(mAxis[i]) > max_value){
                index = i;
                max_value = std::abs(mAxis[i]);
            }
        }

        const int index2 = (index == 2) ? 0 : index + 1;

        mPlaneVector1[index2] = mAxis[index];
        mPlaneVector1[index] = -mAxis[index2];

        // create transformed node vecs
        mOriginNodes.resize(mrOriginModelPart.Nodes().size());
        mTransformedOriginNodes.resize(mrOriginModelPart.Nodes().size());
        for (auto& r_node_i : mrOriginModelPart.Nodes()){
            const int mapping_id = r_node_i.GetValue(MAPPING_ID);
            mOriginNodes[mapping_id] = &r_node_i;
            mTransformedOriginNodes[mapping_id] = GetTransformedNode(r_node_i);
        }

        mDestinationNodes.resize(mrDestinationModelPart.Nodes().size());
        mTransformedDestinationNodes.resize(mrDestinationModelPart.Nodes().size());
        for (auto& r_node_i : mrDestinationModelPart.Nodes()){
            const int mapping_id = r_node_i.GetValue(MAPPING_ID);
            mDestinationNodes[mapping_id] = &r_node_i;
            mTransformedDestinationNodes[mapping_id] = GetTransformedNode(r_node_i);
        }
    }

    NodeVector& GetOriginSearchNodes() override {
        return mTransformedOriginNodes;
    }

    std::vector<std::pair<array_3d, bool>> GetDestinationSearchNodes(const size_t MappingId) override {
        return {
            std::make_pair(mTransformedDestinationNodes[MappingId]->Coordinates(), true),
        };
    }

    BoundedMatrix<double, 3, 3> TransformationMatrix(const size_t DestinationMappingId, const size_t OriginMappingId) const override
    {
        const NodeType& r_origin_node = *mOriginNodes[OriginMappingId];
        const NodeType& r_destination_node = *mDestinationNodes[DestinationMappingId];

        // find angle between ortho
        const array_3d origin_v1 = r_origin_node.Coordinates() - mPoint;
        const array_3d origin_along_axis = inner_prod(mAxis, origin_v1) * mAxis;
        array_3d origin_ortho = origin_v1 - origin_along_axis;

        const double norm_origin = norm_2(origin_ortho);
        if (norm_origin < std::numeric_limits<double>::epsilon()) {
            BoundedMatrix<double, 3, 3 > r = ZeroMatrix(3,3);
            r(0,0) = mAxis[0];
            r(1,1) = mAxis[1];
            r(2,2) = mAxis[2];
            return r;
        }
        origin_ortho /= norm_origin;

        const array_3d destination_v1 = r_destination_node.Coordinates() - mPoint;
        const array_3d destination_along_axis = inner_prod(mAxis, destination_v1) * mAxis;
        array_3d destination_ortho = destination_v1 - destination_along_axis;

        const double norm_destination = norm_2(destination_ortho);
        if (norm_destination < std::numeric_limits<double>::epsilon()) {
            BoundedMatrix<double, 3, 3 > r = ZeroMatrix(3,3);
            r(0,0) = mAxis[0];
            r(1,1) = mAxis[1];
            r(2,2) = mAxis[2];
            return r;
        }
        destination_ortho /= norm_destination;

        double dot_prod = inner_prod(origin_ortho, destination_ortho);
        // result of inner_prod can be slightly larger (smaller) than 1 (-1). This results in NAN for acos!
        if (dot_prod >= 1.0) dot_prod = 1.0;
        else if (dot_prod <= -1.0) dot_prod = -1.0;

        double angle = acos(dot_prod);
        if (inner_prod(mAxis, MathUtils<double>::CrossProduct(origin_ortho, destination_ortho)) < 0) { // Or > 0
            angle = -angle;
        }

        BoundedMatrix<double, 3, 3 > r;
		const double c=cos(angle);
		const double s=sin(angle);
        const double t = 1-c;
        const array_3d& k = mAxis;

        r(0,0) = t*k[0]*k[0] + c;         r(0,1) = t*k[0]*k[1] - k[2]*s;    r(0,2) = t*k[0]*k[2] + k[1]*s;
        r(1,0) = t*k[0]*k[1] + k[2]*s;    r(1,1) = t*k[1]*k[1] + c;         r(1,2) = t*k[1]*k[2] - k[0]*s;
        r(2,0) = t*k[0]*k[2] - k[1]*s;    r(2,1) = t*k[1]*k[2] + k[0]*s;    r(2,2) = t*k[2]*k[2] + c;

        return r;
    }

    NodeTypePointer GetTransformedNode(NodeType& rNode) {
        NodeTypePointer p_new_node = Kratos::make_intrusive<NodeType>(rNode.Id(), rNode[0], rNode[1], rNode[2]);
        p_new_node->SetValue(MAPPING_ID, rNode.GetValue(MAPPING_ID));

        const array_3d v1 = p_new_node->Coordinates() - mPoint;
        const array_3d along_axis = inner_prod(mAxis, v1) * mAxis;
        const array_3d ortho = v1 - along_axis;

        p_new_node->Coordinates() = mPoint + along_axis + mPlaneVector1 * norm_2(ortho);

        return p_new_node;
    }

    array_3d mPoint;
    array_3d mAxis;
    array_3d mPlaneVector1; // vector in the plane orthogonal to the mAxis

    NodeVector mOriginNodes;
    NodeVector mDestinationNodes;
    NodeVector mTransformedOriginNodes;
    NodeVector mTransformedDestinationNodes;


}; // Class SymmetryRevolution

}  // namespace Kratos.

#endif // SYMMETRY_REVOLUTION_H
