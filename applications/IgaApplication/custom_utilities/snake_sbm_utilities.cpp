//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Andrea Gorgi
//                   Nicoló Antonelli
//

// System includes
#include <algorithm>
#include <cmath>

// Project includes
#include "custom_utilities/snake_sbm_utilities.h"

namespace Kratos
{

//-------------------------------

double SnakeSbmUtilities::Orientation(
    const Node& rPointP,
    const Node& rPointQ,
    const Node& rPointR)
{
    return (rPointQ.X() - rPointP.X()) * (rPointR.Y() - rPointP.Y()) -
           (rPointQ.Y() - rPointP.Y()) * (rPointR.X() - rPointP.X());
}

//-------------------------------

bool SnakeSbmUtilities::OnSegment(
    const Node& rPointP,
    const Node& rPointQ,
    const Node& rPointR)
{
    return std::min(rPointP.X(), rPointR.X()) <= rPointQ.X() && rPointQ.X() <= std::max(rPointP.X(), rPointR.X()) &&
           std::min(rPointP.Y(), rPointR.Y()) <= rPointQ.Y() && rPointQ.Y() <= std::max(rPointP.Y(), rPointR.Y());
}

//-------------------------------

bool SnakeSbmUtilities::SegmentsIntersect(
    const Node& rPointA,
    const Node& rPointB,
    const Node& rPointC,
    const Node& rPointD)
{
    const double orientation_1 = Orientation(rPointA, rPointB, rPointC);
    const double orientation_2 = Orientation(rPointA, rPointB, rPointD);
    const double orientation_3 = Orientation(rPointC, rPointD, rPointA);
    const double orientation_4 = Orientation(rPointC, rPointD, rPointB);

    if (orientation_1 * orientation_2 < 0.0 && orientation_3 * orientation_4 < 0.0) {
        return true;
    }

    if (std::abs(orientation_1) < 1e-14 && OnSegment(rPointA, rPointC, rPointB)) return true;
    if (std::abs(orientation_2) < 1e-14 && OnSegment(rPointA, rPointD, rPointB)) return true;
    if (std::abs(orientation_3) < 1e-14 && OnSegment(rPointC, rPointA, rPointD)) return true;
    if (std::abs(orientation_4) < 1e-14 && OnSegment(rPointC, rPointB, rPointD)) return true;
    return false;
}

//-------------------------------

bool SnakeSbmUtilities::SegmentsIntersect(
    const Node::Pointer& pPointA,
    const Node::Pointer& pPointB,
    const Node::Pointer& pPointC,
    const Node::Pointer& pPointD)
{
    return SegmentsIntersect(*pPointA, *pPointB, *pPointC, *pPointD);
}

//-------------------------------

bool SnakeSbmUtilities::SegmentsIntersect(
    const Node::Pointer& pPointA,
    const Node::Pointer& pPointB,
    const Node::Pointer& pPointC,
    const Node::Pointer& pPointD,
    array_1d<double, 3>& rIntersection)
{
    const double x1 = pPointA->X(), y1 = pPointA->Y();
    const double x2 = pPointB->X(), y2 = pPointB->Y();
    const double x3 = pPointC->X(), y3 = pPointC->Y();
    const double x4 = pPointD->X(), y4 = pPointD->Y();
    const double dx12 = x2 - x1;
    const double dy12 = y2 - y1;
    const double dx34 = x4 - x3;
    const double dy34 = y4 - y3;
    const double denominator = dx12 * dy34 - dy12 * dx34;
    constexpr double epsilon = 1.0e-14;

    if (std::abs(denominator) < epsilon) {
        const double orientation_1 = Orientation(*pPointA, *pPointB, *pPointC);
        const double orientation_2 = Orientation(*pPointA, *pPointB, *pPointD);
        if (std::abs(orientation_1) < epsilon || std::abs(orientation_2) < epsilon) {
            if (OnSegment(*pPointA, *pPointC, *pPointB)) {
                rIntersection[0] = x3;
                rIntersection[1] = y3;
                rIntersection[2] = 0.0;
                return true;
            }
            if (OnSegment(*pPointA, *pPointD, *pPointB)) {
                rIntersection[0] = x4;
                rIntersection[1] = y4;
                rIntersection[2] = 0.0;
                return true;
            }
            if (OnSegment(*pPointC, *pPointA, *pPointD)) {
                rIntersection[0] = x1;
                rIntersection[1] = y1;
                rIntersection[2] = 0.0;
                return true;
            }
            if (OnSegment(*pPointC, *pPointB, *pPointD)) {
                rIntersection[0] = x2;
                rIntersection[1] = y2;
                rIntersection[2] = 0.0;
                return true;
            }
        }
        return false;
    }

    const double t = ((x3 - x1) * dy34 - (y3 - y1) * dx34) / denominator;
    const double u = ((x3 - x1) * dy12 - (y3 - y1) * dx12) / denominator;
    if (t < -epsilon || t > 1.0 + epsilon || u < -epsilon || u > 1.0 + epsilon) {
        return false;
    }

    rIntersection[0] = x1 + t * dx12;
    rIntersection[1] = y1 + t * dy12;
    rIntersection[2] = pPointA->Z() + t * (pPointB->Z() - pPointA->Z());
    return true;
}

} // namespace Kratos
