//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/element.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "includes/ublas_interface.h"

// Application includes

namespace Kratos
{
///@name Kratos Classes
///@{

class ScalarWallFluxConditionData
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

    ///@}
    ///@name Classes
    ///@{

    struct Parameters
    {
        double mYPlus;
        double mWallTurbulentViscosity;
        double mDensity;
        double mKinematicViscosity;
        double mKappa;
    };

    ///@}
    ///@name Life Cycle
    ///@{

    ScalarWallFluxConditionData(
        const GeometryType& rGeometry,
        const Properties& rConditionProperties,
        const ProcessInfo& rProcessInfo)
        : mrGeometry(rGeometry),
          mrConditionProperties(rConditionProperties),
          mrElementProperties(rGeometry.GetValue(NEIGHBOUR_ELEMENTS)[0].GetProperties())
    {
    }

    ///@}
    ///@name Access
    ///@{

    const GeometryType& GetGeometry() const
    {
        return mrGeometry;
    }

    const Properties& GetElementProperties() const
    {
        return mrElementProperties;
    }

    const Properties& GetConditionProperties() const
    {
        return mrConditionProperties;
    }

    ///@}

private:
    ///@name Private Members
    ///@{

    const GeometryType& mrGeometry;
    const Properties& mrConditionProperties;
    const Properties& mrElementProperties;

    ///@}
};
} // namespace Kratos