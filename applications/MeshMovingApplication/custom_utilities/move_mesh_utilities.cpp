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

    block_for_each(rNodes, [](Node<3>& rNode ){
        noalias(rNode.Coordinates()) = rNode.GetInitialPosition() + rNode.FastGetSolutionStepValue(MESH_DISPLACEMENT);
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

void SuperImposeVariables(ModelPart &rModelPart, const Variable< array_1d<double, 3> >& rVariable,
                                                 const Variable< array_1d<double, 3> >& rVariableToSuperImpose)
{
    KRATOS_TRY;

    block_for_each(rModelPart.Nodes(), [&](Node<3>& rNode){
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
