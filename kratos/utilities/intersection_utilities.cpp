//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/node.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "utilities/intersection_utilities.h"

namespace Kratos
{

template<>
PointerVector<Point> IntersectionUtilities::ComputeShortestLineBetweenTwoLines<Line2D2<Point>>(
    const Line2D2<Point>& rSegment1,
    const Line2D2<Point>& rSegment2
    )
{
    // Zero tolerance
    const double zero_tolerance = std::numeric_limits<double>::epsilon();

    // Resulting line segment
    auto pa = Kratos::make_shared<Point>(0.0, 0.0, 0.0);
    auto pb = Kratos::make_shared<Point>(0.0, 0.0, 0.0);
    auto resulting_line = PointerVector<Point>();

    // Variable definitions
    array_1d<double, 2> p13,p43,p21;
    double d1343,d4321,d1321,d4343,d2121;
    double mua, mub;
    double numer,denom;

    // Points segments
    const Point& p1 = rSegment1[0];
    const Point& p2 = rSegment1[1];
    const Point& p3 = rSegment2[0];
    const Point& p4 = rSegment2[1];

    p13[0] = p1.X() - p3.X();
    p13[1] = p1.Y() - p3.Y();

    p43[0] = p4.X() - p3.X();
    p43[1] = p4.Y() - p3.Y();
    if (std::abs(p43[0]) < zero_tolerance && std::abs(p43[1]) < zero_tolerance)
        return resulting_line;

    p21[0] = p2.X() - p1.X();
    p21[1] = p2.Y() - p1.Y();
    if (std::abs(p21[0]) < zero_tolerance && std::abs(p21[1]) < zero_tolerance)
        return resulting_line;

    d1343 = p13[0] * p43[0] + p13[1] * p43[1];
    d4321 = p43[0] * p21[0] + p43[1] * p21[1];
    d1321 = p13[0] * p21[0] + p13[1] * p21[1];
    d4343 = p43[0] * p43[0] + p43[1] * p43[1];
    d2121 = p21[0] * p21[0] + p21[1] * p21[1];

    denom = d2121 * d4343 - d4321 * d4321;
    if (std::abs(denom) < zero_tolerance)
        return resulting_line;
    numer = d1343 * d4321 - d1321 * d4343;

    mua = numer / denom;
    mub = (d1343 + d4321 * mua) / d4343;

    pa->X() = p1.X() + mua * p21[0];
    pa->Y() = p1.Y() + mua * p21[1];
    pb->X() = p3.X() + mub * p43[0];
    pb->Y() = p3.Y() + mub * p43[1];

    resulting_line.push_back(pa);
    resulting_line.push_back(pb);
    return resulting_line;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVector<Point> IntersectionUtilities::ComputeShortestLineBetweenTwoLines<Line3D2<Point>>(
    const Line3D2<Point>& rSegment1,
    const Line3D2<Point>& rSegment2
    )
{
    // Zero tolerance
    const double zero_tolerance = std::numeric_limits<double>::epsilon();

    // Resulting line segment
    auto pa = Kratos::make_shared<Point>(0.0, 0.0, 0.0);
    auto pb = Kratos::make_shared<Point>(0.0, 0.0, 0.0);
    auto resulting_line = PointerVector<Point>();

    // Variable definitions
    array_1d<double, 3> p13,p43,p21;
    double d1343,d4321,d1321,d4343,d2121;
    double mua, mub;
    double numer,denom;

    // Points segments
    const Point& p1 = rSegment1[0];
    const Point& p2 = rSegment1[1];
    const Point& p3 = rSegment2[0];
    const Point& p4 = rSegment2[1];

    p13[0] = p1.X() - p3.X();
    p13[1] = p1.Y() - p3.Y();
    p13[2] = p1.Z() - p3.Z();

    p43[0] = p4.X() - p3.X();
    p43[1] = p4.Y() - p3.Y();
    p43[2] = p4.Z() - p3.Z();
    if (std::abs(p43[0]) < zero_tolerance && std::abs(p43[1]) < zero_tolerance && std::abs(p43[2]) < zero_tolerance)
        return resulting_line;

    p21[0] = p2.X() - p1.X();
    p21[1] = p2.Y() - p1.Y();
    p21[2] = p2.Z() - p1.Z();
    if (std::abs(p21[0]) < zero_tolerance && std::abs(p21[1]) < zero_tolerance && std::abs(p21[2]) < zero_tolerance)
        return resulting_line;

    d1343 = p13[0] * p43[0] + p13[1] * p43[1] + p13[2] * p43[2];
    d4321 = p43[0] * p21[0] + p43[1] * p21[1] + p43[2] * p21[2];
    d1321 = p13[0] * p21[0] + p13[1] * p21[1] + p13[2] * p21[2];
    d4343 = p43[0] * p43[0] + p43[1] * p43[1] + p43[2] * p43[2];
    d2121 = p21[0] * p21[0] + p21[1] * p21[1] + p21[2] * p21[2];

    denom = d2121 * d4343 - d4321 * d4321;
    if (std::abs(denom) < zero_tolerance)
        return resulting_line;
    numer = d1343 * d4321 - d1321 * d4343;

    mua = numer / denom;
    mub = (d1343 + d4321 * mua) / d4343;

    pa->X() = p1.X() + mua * p21[0];
    pa->Y() = p1.Y() + mua * p21[1];
    pa->Z() = p1.Z() + mua * p21[2];
    pb->X() = p3.X() + mub * p43[0];
    pb->Y() = p3.Y() + mub * p43[1];
    pb->Z() = p3.Z() + mub * p43[2];

    resulting_line.push_back(pa);
    resulting_line.push_back(pb);
    return resulting_line;
}

/***********************************************************************************/
/*** This duplicates code for Node<3>. Explicit instantation gives compi. errors ***/
/***********************************************************************************/

template<>
PointerVector<Point> IntersectionUtilities::ComputeShortestLineBetweenTwoLines<Line2D2<Node<3>>>(
    const Line2D2<Node<3>>& rSegment1,
    const Line2D2<Node<3>>& rSegment2
    )
{
    // Zero tolerance
    const double zero_tolerance = std::numeric_limits<double>::epsilon();

    // Resulting line segment
    auto pa = Kratos::make_shared<Point>(0.0, 0.0, 0.0);
    auto pb = Kratos::make_shared<Point>(0.0, 0.0, 0.0);
    auto resulting_line = PointerVector<Point>();

    // Variable definitions
    array_1d<double, 2> p13,p43,p21;
    double d1343,d4321,d1321,d4343,d2121;
    double mua, mub;
    double numer,denom;

    // Points segments
    const Point& p1 = rSegment1[0];
    const Point& p2 = rSegment1[1];
    const Point& p3 = rSegment2[0];
    const Point& p4 = rSegment2[1];

    p13[0] = p1.X() - p3.X();
    p13[1] = p1.Y() - p3.Y();

    p43[0] = p4.X() - p3.X();
    p43[1] = p4.Y() - p3.Y();
    if (std::abs(p43[0]) < zero_tolerance && std::abs(p43[1]) < zero_tolerance)
        return resulting_line;

    p21[0] = p2.X() - p1.X();
    p21[1] = p2.Y() - p1.Y();
    if (std::abs(p21[0]) < zero_tolerance && std::abs(p21[1]) < zero_tolerance)
        return resulting_line;

    d1343 = p13[0] * p43[0] + p13[1] * p43[1];
    d4321 = p43[0] * p21[0] + p43[1] * p21[1];
    d1321 = p13[0] * p21[0] + p13[1] * p21[1];
    d4343 = p43[0] * p43[0] + p43[1] * p43[1];
    d2121 = p21[0] * p21[0] + p21[1] * p21[1];

    denom = d2121 * d4343 - d4321 * d4321;
    if (std::abs(denom) < zero_tolerance)
        return resulting_line;
    numer = d1343 * d4321 - d1321 * d4343;

    mua = numer / denom;
    mub = (d1343 + d4321 * mua) / d4343;

    pa->X() = p1.X() + mua * p21[0];
    pa->Y() = p1.Y() + mua * p21[1];
    pb->X() = p3.X() + mub * p43[0];
    pb->Y() = p3.Y() + mub * p43[1];

    resulting_line.push_back(pa);
    resulting_line.push_back(pb);
    return resulting_line;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVector<Point> IntersectionUtilities::ComputeShortestLineBetweenTwoLines<Line3D2<Node<3>>>(
    const Line3D2<Node<3>>& rSegment1,
    const Line3D2<Node<3>>& rSegment2
    )
{
    // Zero tolerance
    const double zero_tolerance = std::numeric_limits<double>::epsilon();

    // Resulting line segment
    auto pa = Kratos::make_shared<Point>(0.0, 0.0, 0.0);
    auto pb = Kratos::make_shared<Point>(0.0, 0.0, 0.0);
    auto resulting_line = PointerVector<Point>();

    // Variable definitions
    array_1d<double, 3> p13,p43,p21;
    double d1343,d4321,d1321,d4343,d2121;
    double mua, mub;
    double numer,denom;

    // Points segments
    const Point& p1 = rSegment1[0];
    const Point& p2 = rSegment1[1];
    const Point& p3 = rSegment2[0];
    const Point& p4 = rSegment2[1];

    p13[0] = p1.X() - p3.X();
    p13[1] = p1.Y() - p3.Y();
    p13[2] = p1.Z() - p3.Z();

    p43[0] = p4.X() - p3.X();
    p43[1] = p4.Y() - p3.Y();
    p43[2] = p4.Z() - p3.Z();
    if (std::abs(p43[0]) < zero_tolerance && std::abs(p43[1]) < zero_tolerance && std::abs(p43[2]) < zero_tolerance)
        return resulting_line;

    p21[0] = p2.X() - p1.X();
    p21[1] = p2.Y() - p1.Y();
    p21[2] = p2.Z() - p1.Z();
    if (std::abs(p21[0]) < zero_tolerance && std::abs(p21[1]) < zero_tolerance && std::abs(p21[2]) < zero_tolerance)
        return resulting_line;

    d1343 = p13[0] * p43[0] + p13[1] * p43[1] + p13[2] * p43[2];
    d4321 = p43[0] * p21[0] + p43[1] * p21[1] + p43[2] * p21[2];
    d1321 = p13[0] * p21[0] + p13[1] * p21[1] + p13[2] * p21[2];
    d4343 = p43[0] * p43[0] + p43[1] * p43[1] + p43[2] * p43[2];
    d2121 = p21[0] * p21[0] + p21[1] * p21[1] + p21[2] * p21[2];

    denom = d2121 * d4343 - d4321 * d4321;
    if (std::abs(denom) < zero_tolerance)
        return resulting_line;
    numer = d1343 * d4321 - d1321 * d4343;

    mua = numer / denom;
    mub = (d1343 + d4321 * mua) / d4343;

    pa->X() = p1.X() + mua * p21[0];
    pa->Y() = p1.Y() + mua * p21[1];
    pa->Z() = p1.Z() + mua * p21[2];
    pb->X() = p3.X() + mub * p43[0];
    pb->Y() = p3.Y() + mub * p43[1];
    pb->Z() = p3.Z() + mub * p43[2];

    resulting_line.push_back(pa);
    resulting_line.push_back(pb);
    return resulting_line;
}

}  /* namespace Kratos.*/