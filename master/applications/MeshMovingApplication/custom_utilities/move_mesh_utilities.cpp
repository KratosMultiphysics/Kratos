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

// System includes

// External includes

// Project includes
#include "move_mesh_utilities.h"
#include "containers/model.h"
#include "custom_utilities/parametric_affine_transform.h"
#include "includes/mesh_moving_variables.h" // TODO remove after mesh-vel-comp-functions are removed
#include "utilities/parallel_utilities.h"

namespace Kratos {
namespace MoveMeshUtilities {

//******************************************************************************
//******************************************************************************
void CheckJacobianDimension(GeometryType::JacobiansType &rInvJ0,
                            VectorType &rDetJ0, const GeometryType &rGeometry) {
  KRATOS_TRY;

  const auto this_integration_method = rGeometry.GetDefaultIntegrationMethod();
  const auto& r_integration_points = rGeometry.IntegrationPoints(this_integration_method);

  if (rInvJ0.size() != r_integration_points.size()) {
    rInvJ0.resize(r_integration_points.size());
  }
  if (rDetJ0.size() != r_integration_points.size()) {
    rDetJ0.resize(r_integration_points.size());
  }

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
void MoveMesh(ModelPart::NodesContainerType& rNodes) {
    KRATOS_TRY;

    block_for_each(rNodes, [](Node& rNode ){
        noalias(rNode.Coordinates()) = rNode.GetInitialPosition() + rNode.FastGetSolutionStepValue(MESH_DISPLACEMENT);
    });

    KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
void MoveModelPart(
    ModelPart& rModelPart,
    const array_1d<double,3>& rRotationAxis,
    const double rotationAngle,
    const array_1d<double,3>& rReferencePoint,
    const array_1d<double,3>& rTranslationVector)
{
    KRATOS_TRY

    const AffineTransform transform(
        rRotationAxis,
        rotationAngle,
        rReferencePoint,
        rTranslationVector);

    MoveModelPart(rModelPart, transform);

    KRATOS_CATCH("");
}

void MoveModelPart(
    ModelPart& rModelPart,
    const Parameters rotationAxis,
    const Parameters rotationAngle,
    const Parameters referencePoint,
    const Parameters translationVector)
{
    KRATOS_TRY

    ParametricAffineTransform transform(
        rotationAxis,
        rotationAngle,
        referencePoint,
        translationVector);

    MoveModelPart(rModelPart, transform);

    KRATOS_CATCH("");
}

void MoveModelPart(
    ModelPart& rModelPart,
    const AffineTransform& rTransform)
{
    KRATOS_TRY

    block_for_each(
        rModelPart.Nodes(),
        [&rTransform](Node& rNode){
            const array_1d<double,3>& initial_position = rNode.GetInitialPosition();
            noalias(rNode.GetSolutionStepValue(MESH_DISPLACEMENT)) = rTransform.Apply(initial_position) - initial_position;
        });

    KRATOS_CATCH("");
}

void MoveModelPart(
    ModelPart& rModelPart,
    ParametricAffineTransform& rTransform)
{
    KRATOS_TRY

    const double time = rModelPart.GetProcessInfo().GetValue(TIME);

    block_for_each(
        rModelPart.Nodes(),
        rTransform,
        [time](Node& rNode, ParametricAffineTransform& rTLSTransform){
            const array_1d<double,3>& initial_position = rNode.GetInitialPosition();
            noalias(rNode.GetSolutionStepValue(MESH_DISPLACEMENT)) = rTLSTransform.Apply(
                initial_position,
                time,
                rNode.X0(),
                rNode.Y0(),
                rNode.Z0()) - initial_position;
        });

    KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
ModelPart* GenerateMeshPart(ModelPart &rModelPart,
                                    const std::string &rElementName) {
  KRATOS_TRY;

  ModelPart* pmesh_model_part = &(rModelPart.GetModel().CreateModelPart(rModelPart.Name()+"_MeshPart", 1));

  // initializing mesh nodes and variables
  pmesh_model_part->Nodes() = rModelPart.Nodes();

  // creating mesh elements
  ModelPart::ElementsContainerType &rmesh_elements =
      pmesh_model_part->Elements();

  const Element &r_reference_element =
      KratosComponents<Element>::Get(rElementName);


  // create a new property for the mesh motion
  Properties::Pointer p_mesh_motion_property = pmesh_model_part->CreateNewProperties(0);

  for (int i = 0; i < (int)rModelPart.Elements().size(); i++) {
    ModelPart::ElementsContainerType::iterator it =
        rModelPart.ElementsBegin() + i;
    Element::Pointer p_element = r_reference_element.Create(
        it->Id(), it->pGetGeometry(), p_mesh_motion_property);
    rmesh_elements.push_back(p_element);
  }

  return std::move(pmesh_model_part);

  KRATOS_CATCH("");
}

void InitializeMeshPartWithElements(
    ModelPart& rDestinationModelPart,
    ModelPart& rOriginModelPart,
    Properties::Pointer pProperties,
    const std::string& rElementName)
{
    KRATOS_TRY

    // initializing mesh nodes and variables
    rDestinationModelPart.Nodes() = rOriginModelPart.Nodes();

    auto& r_elements = rDestinationModelPart.Elements();

    // clear all existing elements
    r_elements.clear();

    const Element& r_reference_element = KratosComponents<Element>::Get(rElementName);

    KRATOS_ERROR_IF(rOriginModelPart.GetCommunicator().GlobalNumberOfElements() == 0)
        << "No elements are found in " << rOriginModelPart.Name()
        << " to initialize " << rDestinationModelPart.Name() << ".\n";

    for (auto& r_origin_element : rOriginModelPart.Elements()) {
        auto p_destination_element = r_reference_element.Create(
            r_origin_element.Id(), r_origin_element.pGetGeometry(), pProperties);
        r_elements.push_back(p_destination_element);
    }

    KRATOS_CATCH("");
}

void SuperImposeVariables(ModelPart &rModelPart, const Variable< array_1d<double, 3> >& rVariable,
                                                 const Variable< array_1d<double, 3> >& rVariableToSuperImpose)
{
    KRATOS_TRY;

    block_for_each(rModelPart.Nodes(), [&](Node& rNode){
        if (rNode.Has(rVariableToSuperImpose)) {
            rNode.GetSolutionStepValue(rVariable,0) += rNode.GetValue(rVariableToSuperImpose);
        }
    });

  KRATOS_CATCH("");
}

void SuperImposeMeshDisplacement(ModelPart &rModelPart, const Variable< array_1d<double, 3> >& rVariableToSuperImpose)
{
  SuperImposeVariables(rModelPart, MESH_DISPLACEMENT, rVariableToSuperImpose);
}

void SuperImposeMeshVelocity(ModelPart &rModelPart, const Variable< array_1d<double, 3> >& rVariableToSuperImpose)
{
  SuperImposeVariables(rModelPart, MESH_VELOCITY, rVariableToSuperImpose);
}

} // namespace Move Mesh Utilities.

} // namespace Kratos.
