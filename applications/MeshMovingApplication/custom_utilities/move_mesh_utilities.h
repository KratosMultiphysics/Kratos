//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//

#if !defined(KRATOS_MESHMOVING_UTILITIES_H_INCLUDED)
#define KRATOS_MESHMOVING_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos {
namespace MoveMeshUtilities {

typedef Element::GeometryType GeometryType;
typedef Element::VectorType VectorType;

void KRATOS_API(MESH_MOVING_APPLICATION) CheckJacobianDimension(GeometryType::JacobiansType &rInvJ0,
                            VectorType &rDetJ0, const GeometryType &rGeometry);

void KRATOS_API(MESH_MOVING_APPLICATION) MoveMesh(ModelPart::NodesContainerType &rNodes);

KRATOS_API(MESH_MOVING_APPLICATION) ModelPart* GenerateMeshPart(ModelPart &rModelPart,
                                    const std::string &rElementName);

void KRATOS_API(MESH_MOVING_APPLICATION) SuperImposeVariables(ModelPart &rModelPart, const Variable< array_1d<double, 3> >& rVariable,
                                                 const Variable< array_1d<double, 3> >& rVariableToSuperImpose);

void KRATOS_API(MESH_MOVING_APPLICATION) SuperImposeMeshDisplacement(ModelPart &rModelPart, const Variable< array_1d<double, 3> >& rVariableToSuperImpose);

void KRATOS_API(MESH_MOVING_APPLICATION) SuperImposeMeshVelocity(ModelPart &rModelPart, const Variable< array_1d<double, 3> >& rVariableToSuperImpose);

} // namespace Move Mesh Utilities.

} // namespace Kratos.

#endif // KRATOS_MESHMOVING_UTILITIES_H_INCLUDED  defined
