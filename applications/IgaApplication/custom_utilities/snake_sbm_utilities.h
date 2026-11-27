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
#pragma once

// Project includes
#include "includes/node.h"
#include "includes/define.h"

namespace Kratos
{

class KRATOS_API(IGA_APPLICATION) SnakeSbmUtilities
{
public:
//-------------------------------
    static double Orientation(const Node& rPointP, const Node& rPointQ, const Node& rPointR);
//-------------------------------
    static bool OnSegment(const Node& rPointP, const Node& rPointQ, const Node& rPointR);
//-------------------------------
    static bool SegmentsIntersect(
        const Node& rPointA,
        const Node& rPointB,
        const Node& rPointC,
        const Node& rPointD);
//-------------------------------
    static bool SegmentsIntersect(
        const Node::Pointer& pPointA,
        const Node::Pointer& pPointB,
        const Node::Pointer& pPointC,
        const Node::Pointer& pPointD);
//-------------------------------
    static bool SegmentsIntersect(
        const Node::Pointer& pPointA,
        const Node::Pointer& pPointB,
        const Node::Pointer& pPointC,
        const Node::Pointer& pPointD,
        array_1d<double, 3>& rIntersection);
};

} // namespace Kratos
