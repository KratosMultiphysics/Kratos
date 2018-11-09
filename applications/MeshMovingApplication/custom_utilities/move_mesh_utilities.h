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
#include "mesh_moving_application.h"
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos {
namespace MoveMeshUtilities {

typedef Element BaseType;
typedef Element::GeometryType GeometryType;
typedef GeometryData::IntegrationMethod IntegrationMethod;
typedef BaseType::VectorType VectorType;

void CheckJacobianDimension(GeometryType::JacobiansType &rInvJ0,
                            VectorType &rDetJ0, GeometryType &rGeometry);

void CalculateMeshVelocities(ModelPart &rMeshModelPart,
                             const int TimeOrder, const double DeltaTime);

void CalculateMeshVelocities(ModelPart* pMeshModelPart,
                             const int TimeOrder, const double DeltaTime);

void MoveMesh(const ModelPart::NodesContainerType &rNodes);

KRATOS_DEPRECATED_MESSAGE("This is legacy version, please use the \"Kratos.VariableUtils.UpdateCurrentToInitialConfiguration\"")
void SetMeshToInitialConfiguration(const ModelPart::NodesContainerType &rNodes);

ModelPart* GenerateMeshPart(ModelPart &rModelPart,
                                    const std::string &rElementName);

KRATOS_DEPRECATED_MESSAGE("This is legacy version, please use the \"Kratos.VariableUtils.UpdateInitialToCurrentConfiguration\"")
void UpdateReferenceMesh(ModelPart &rModelPart);

} // namespace Move Mesh Utilities.

} // namespace Kratos.

#endif // KRATOS_MESHMOVING_UTILITIES_H_INCLUDED  defined
