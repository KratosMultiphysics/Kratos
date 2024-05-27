//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
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
#include "affine_transform.h"
#include "parametric_affine_transform.h"

namespace Kratos::MoveMeshUtilities {

typedef Element::GeometryType GeometryType;
typedef Element::VectorType VectorType;

void KRATOS_API(MESH_MOVING_APPLICATION) CheckJacobianDimension(GeometryType::JacobiansType &rInvJ0,
                            VectorType &rDetJ0, const GeometryType &rGeometry);

void KRATOS_API(MESH_MOVING_APPLICATION) MoveMesh(ModelPart::NodesContainerType &rNodes);

/** @brief Impose mesh movement on all nodes of a model part.
 *  @details the movement is a linear transformation (see @ref AffineTransform )
 *           defined by the specified by the arguments.
 *  @param rModelPart model part containing the nodes to set MESH_DISPLACEMENT on.
 *  @param rRotationAxis axis of rotation.
 *  @param rotationAngle angle of rotation (radians).
 *  @param rReferencePoint point on the axis of rotation.
 *  @param rTranslationVector translation vector (applied AFTER the rotation).
 */
void KRATOS_API(MESH_MOVING_APPLICATION) MoveModelPart(
    ModelPart& rModelPart,
    const array_1d<double,3>& rRotationAxis,
    const double rotationAngle,
    const array_1d<double,3>& rReferencePoint,
    const array_1d<double,3>& rTranslationVector);

/** @brief Impose parametric mesh movement on all nodes of a model part, as a function of the current position, time, and initial position.
 *  @details the movement is a linear transformation (see @ref ParametricAffineTransform )
 *           defined by the specified by the arguments.
 *  @param rModelPart model part containing the nodes to set MESH_DISPLACEMENT on.
 *  @param rRotationAxis axis of rotation (vector of size 3).
 *  @param rRotationAngle angle of rotation (scalar - radians).
 *  @param rReferencePoint point on the axis of rotation (vector of size 3).
 *  @param rTranslationVector translation vector (vector of size 3 - applied AFTER the rotation).
 */
void KRATOS_API(MESH_MOVING_APPLICATION) MoveModelPart(
    ModelPart& rModelPart,
    const Parameters rotationAxis,
    const Parameters rotationAngle,
    const Parameters referencePoint,
    const Parameters translationVector);

/// @brief Impose mesh movement on all nodes of a model part.
void KRATOS_API(MESH_MOVING_APPLICATION) MoveModelPart(
    ModelPart& rModelPart,
    const AffineTransform& rTransform);

/// @brief Impose parametric mesh movement on all nodes of a model part.
void KRATOS_API(MESH_MOVING_APPLICATION) MoveModelPart(
    ModelPart& rModelPart,
    ParametricAffineTransform& rTransform);

KRATOS_API(MESH_MOVING_APPLICATION) ModelPart* GenerateMeshPart(ModelPart &rModelPart,
                                    const std::string &rElementName);

KRATOS_API(MESH_MOVING_APPLICATION) void InitializeMeshPartWithElements(
    ModelPart& rDestinationModelPart,
    ModelPart& rOriginModelPart,
    Properties::Pointer pProperties,
    const std::string& rElementName);

void KRATOS_API(MESH_MOVING_APPLICATION) SuperImposeVariables(ModelPart &rModelPart, const Variable< array_1d<double, 3> >& rVariable,
                                                 const Variable< array_1d<double, 3> >& rVariableToSuperImpose);

void KRATOS_API(MESH_MOVING_APPLICATION) SuperImposeMeshDisplacement(ModelPart &rModelPart, const Variable< array_1d<double, 3> >& rVariableToSuperImpose);

void KRATOS_API(MESH_MOVING_APPLICATION) SuperImposeMeshVelocity(ModelPart &rModelPart, const Variable< array_1d<double, 3> >& rVariableToSuperImpose);

} // namespace Kratos::MoveMeshUtilities

#endif // KRATOS_MESHMOVING_UTILITIES_H_INCLUDED  defined
