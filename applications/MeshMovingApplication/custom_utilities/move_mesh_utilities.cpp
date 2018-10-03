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

namespace Kratos {
namespace MoveMeshUtilities {

//******************************************************************************
//******************************************************************************
void CheckJacobianDimension(GeometryType::JacobiansType &rInvJ0,
                            VectorType &rDetJ0, GeometryType &rGeometry) {
  KRATOS_TRY;

  const IntegrationMethod this_integration_method =
      rGeometry.GetDefaultIntegrationMethod();
  const GeometryType::IntegrationPointsArrayType &integration_points =
      rGeometry.IntegrationPoints(this_integration_method);

  if (rInvJ0.size() != integration_points.size())
    rInvJ0.resize(integration_points.size());
  if (rDetJ0.size() != integration_points.size())
    rDetJ0.resize(integration_points.size());

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************

void CalculateMeshVelocities(ModelPart* pMeshModelPart,
                             const int TimeOrder, const double DeltaTime) {

    CalculateMeshVelocities(*pMeshModelPart, TimeOrder, DeltaTime);
}

void CalculateMeshVelocities(ModelPart &rMeshModelPart,
                             const int TimeOrder, const double DeltaTime) {
  KRATOS_TRY;

  KRATOS_ERROR_IF(DeltaTime <= 0.0) << "Invalid DELTA_TIME." << std::endl;

  const double coeff = 1 / DeltaTime;

  if (TimeOrder == 1) {
    for (ModelPart::NodeIterator i =
             rMeshModelPart.GetCommunicator().LocalMesh().NodesBegin();
         i != rMeshModelPart.GetCommunicator().LocalMesh().NodesEnd(); ++i) {

      array_1d<double, 3> &mesh_v =
          (i)->FastGetSolutionStepValue(MESH_VELOCITY);
      array_1d<double, 3> &disp =
          (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT);
      array_1d<double, 3> &dispold =
          (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
      noalias(mesh_v) = disp - dispold;
      mesh_v *= coeff;
    }
  } else if (TimeOrder == 2) {
    const double c1 = 1.50 * coeff;
    const double c2 = -2.0 * coeff;
    const double c3 = 0.50 * coeff;

    for (ModelPart::NodeIterator i =
             rMeshModelPart.GetCommunicator().LocalMesh().NodesBegin();
         i != rMeshModelPart.GetCommunicator().LocalMesh().NodesEnd(); ++i) {

      array_1d<double, 3> &mesh_v =
          (i)->FastGetSolutionStepValue(MESH_VELOCITY);
      noalias(mesh_v) = c1 * (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT);
      noalias(mesh_v) +=
          c2 * (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
      noalias(mesh_v) +=
          c3 * (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT, 2);
    }
  } else {
    KRATOS_ERROR << "Wrong TimeOrder: Acceptable values are: 1 and 2"
                 << std::endl;
  }

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
void MoveMesh(const ModelPart::NodesContainerType &rNodes) {
  KRATOS_TRY;

  for (auto &rnode : rNodes) {
    noalias(rnode.Coordinates()) = rnode.GetInitialPosition()
                     + rnode.FastGetSolutionStepValue(MESH_DISPLACEMENT);
  }

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
void SetMeshToInitialConfiguration(
    const ModelPart::NodesContainerType &rNodes) {
  KRATOS_TRY;

  for (auto &rnode : rNodes) {
    noalias(rnode.Coordinates()) = rnode.GetInitialPosition();
  }

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
ModelPart* GenerateMeshPart(ModelPart &rModelPart,
                                    const std::string &rElementName) {
  KRATOS_TRY;

  ModelPart* pmesh_model_part = &(rModelPart.GetOwnerModel().CreateModelPart("MeshPart", 1));

  // initializing mesh nodes and variables
  pmesh_model_part->Nodes() = rModelPart.Nodes();

  // creating mesh elements
  ModelPart::ElementsContainerType &rmesh_elements =
      pmesh_model_part->Elements();

  const Element &r_reference_element =
      KratosComponents<Element>::Get(rElementName);

  for (int i = 0; i < (int)rModelPart.Elements().size(); i++) {
    ModelPart::ElementsContainerType::iterator it =
        rModelPart.ElementsBegin() + i;
    Element::Pointer p_element = r_reference_element.Create(
        it->Id(), it->pGetGeometry(), it->pGetProperties());
    rmesh_elements.push_back(p_element);
  }

  return std::move(pmesh_model_part);

  KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
  void UpdateReferenceMesh(ModelPart &rModelPart)
{

  for (ModelPart::NodeIterator i = rModelPart.NodesBegin();
          i != rModelPart.NodesEnd(); ++i){

      (i)->X0() = (i)->X();
      (i)->Y0() = (i)->Y();
      (i)->Z0() = (i)->Z();

  }
}

} // namespace Move Mesh Utilities.

} // namespace Kratos.
