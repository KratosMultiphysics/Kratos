// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#pragma once

#include "geometries/point.h"
#include "includes/element.h"

#include <vector>

namespace Kratos::Testing
{

class ElementSetupUtilities
{
public:
    static std::vector<Kratos::Point> CreatePointsFor2D3NElement();
    static std::vector<Kratos::Point> CreatePointsFor2D6NElement();
    static std::vector<Kratos::Point> CreatePointsFor2D10NElement();
    static std::vector<Kratos::Point> CreatePointsFor2D15NElement();
    static std::vector<Kratos::Point> CreatePointsFor3D10NElement();

    static Element::Pointer Create2D3NElement(const PointerVector<Node>& rNodes,
                                              const Properties::Pointer& rProperties);
    static Element::Pointer Create2D3NElement();

    static Element::Pointer Create2D6NElement(const PointerVector<Node>& rNodes,
                                              const Properties::Pointer& rProperties);
    static Element::Pointer Create2D6NElement();

    static Element::Pointer Create2D10NElement(const PointerVector<Node>& rNodes,
                                               const Properties::Pointer& rProperties);
    static Element::Pointer Create2D10NElement();

    static Element::Pointer Create2D15NElement(const PointerVector<Node>& rNodes,
                                               const Properties::Pointer& rProperties);
    static Element::Pointer Create2D15NElement();

    static Element::Pointer Create3D10NElement(const PointerVector<Node>& rNodes,
                                               const Properties::Pointer& rProperties);
    static Element::Pointer Create3D10NElement();
};

} // namespace Kratos::Testing
