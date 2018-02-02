//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//

// System includes

// External includes

// Project includes
#include "move_mesh_utilities.h"

namespace Kratos {
namespace MoveMeshUtilities {
void CheckJacobianDimension(GeometryType::JacobiansType &rInvJ0,
                            VectorType &rDetJ0, GeometryType &rGeometry) {
  const IntegrationMethod this_integration_method =
      rGeometry.GetDefaultIntegrationMethod();
  const GeometryType::IntegrationPointsArrayType &integration_points =
      rGeometry.IntegrationPoints(this_integration_method);

  if (rInvJ0.size() != integration_points.size())
    rInvJ0.resize(integration_points.size());
  if (rDetJ0.size() != integration_points.size())
    rDetJ0.resize(integration_points.size());
}

} // namespace Move Mesh Utilities.

} // namespace Kratos.