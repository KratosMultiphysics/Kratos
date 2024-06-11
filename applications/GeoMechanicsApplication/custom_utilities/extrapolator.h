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

namespace Kratos
{

class Extrapolator
{
public:
    using NodeType     = Node;
    using GeometryType = Geometry<NodeType>;
    using SizeType     = std::size_t;
    using IndexType    = std::size_t;

    struct TLSType {
        Vector vector_J;
        Vector N;
    };

    Matrix CalculateElementExtrapolationMatrix(GeometryType& r_this_geometry,
                                               SizeType      integration_points_number,
                                               GeometryType::IntegrationPointsArrayType& integration_points,
                                               GeometryData::IntegrationMethod this_integration_method,
                                               SizeType number_of_nodes) const;
};

} // namespace Kratos
