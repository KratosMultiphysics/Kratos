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
void MoveMesh(const ModelPart::NodesContainerType& rNodes) {
    KRATOS_TRY;

    const int num_nodes = rNodes.size();
    const auto nodes_begin = rNodes.begin();

    #pragma omp parallel for
    for (int i=0; i<num_nodes; i++) {
        const auto it_node  = nodes_begin + i;
        noalias(it_node->Coordinates()) = it_node->GetInitialPosition()
            + it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT);
    }
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
  auto r_nodes = rModelPart.Nodes();
  const int num_nodes = r_nodes.size();
  const auto nodes_begin = r_nodes.begin();

  #pragma omp parallel for
  for (int i=0; i<num_nodes; i++) {
      const auto it_node  = nodes_begin + i;
      if(it_node->Has(rVariableToSuperImpose))
          it_node->GetSolutionStepValue(rVariable,0) += it_node->GetValue(rVariableToSuperImpose);
  }
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
