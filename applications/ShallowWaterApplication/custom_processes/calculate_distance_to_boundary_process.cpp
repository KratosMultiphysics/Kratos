//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "includes/model_part.h"
#include "geometries/line_2d_2.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "calculate_distance_to_boundary_process.h"


namespace Kratos
{

void CalculateDistanceToBoundaryProcess::FindApproximatingGeometry(
    GeometryType::Pointer& pEntity,
    const ModelPart& rModelPart)
{
    double low, left = std::numeric_limits<double>::max();
    double top, right = std::numeric_limits<double>::lowest();

    typedef CombinedReduction< MinReduction<double>,
                               MinReduction<double>,
                               MaxReduction<double>,
                               MaxReduction<double> > BoundingBoxReduction;

    std::tie(low, left, top, right) = block_for_each<BoundingBoxReduction>(rModelPart.Nodes(), [&](NodeType& rNode){
        double x = rNode.X();
        double y = rNode.Y();
        return std::make_tuple(y, x, y, x);
    });

    auto point_0 = Kratos::make_shared<Point>(left, low);
    auto point_1 = Kratos::make_shared<Point>(right, low);
    auto point_2 = Kratos::make_shared<Point>(right, top);
    auto point_3 = Kratos::make_shared<Point>(left, top);

    Line2D2<Point> line_a(point_0, point_2);
    Line2D2<Point> line_b(point_1, point_3);

    double distance_a, distance_b;

    distance_a = block_for_each<SumReduction<double>>(rModelPart.Nodes(), [&](NodeType& rNode){
        Point projected;
        return (GeometricalProjectionUtilities::FastProjectOnLine2D(line_a, rNode, projected));
    });

    distance_b = block_for_each<SumReduction<double>>(rModelPart.Nodes(), [&](NodeType& rNode){
        Point projected;
        return (GeometricalProjectionUtilities::FastProjectOnLine2D(line_b, rNode, projected));
    });

    if (distance_a < distance_b) {
        pEntity = Kratos::make_shared<Line2D2<Point>>(line_a);
    } else {
        pEntity = Kratos::make_shared<Line2D2<Point>>(line_b);
    }
}


}  // namespace Kratos.
