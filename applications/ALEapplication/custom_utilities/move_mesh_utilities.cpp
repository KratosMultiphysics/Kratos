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


    void CalculateMeshVelocities(ModelPart::Pointer pMeshModelPart, const int TimeOrder, const double DeltaTime) {
    KRATOS_TRY;

    KRATOS_ERROR_IF(DeltaTime <= 0.0)<< "Invalid DELTA_TIME." << std::endl;

    const double coeff = 1 / DeltaTime;

    if (TimeOrder == 1) {
      for (ModelPart::NodeIterator i = (*pMeshModelPart).GetCommunicator().LocalMesh().NodesBegin();
           i != (*pMeshModelPart).GetCommunicator().LocalMesh().NodesEnd(); ++i) {

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

      for (ModelPart::NodeIterator i = (*pMeshModelPart).GetCommunicator().LocalMesh().NodesBegin();
           i != (*pMeshModelPart).GetCommunicator().LocalMesh().NodesEnd(); ++i) {

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


  void MoveMesh(const ModelPart::NodesContainerType& rNodes){
    for (auto& r_node : rNodes) {
      (r_node).X() = (r_node).X0() + r_node.GetSolutionStepValue(MESH_DISPLACEMENT_X);
      (r_node).Y() = (r_node).Y0() + r_node.GetSolutionStepValue(MESH_DISPLACEMENT_Y);
      (r_node).Z() = (r_node).Z0() + r_node.GetSolutionStepValue(MESH_DISPLACEMENT_Z);
    }
  }

  void SetMeshToInitialConfiguration(const ModelPart::NodesContainerType& rNodes) {
    for (auto& r_node : rNodes) {
      (r_node).X() = (r_node).X0();
      (r_node).Y() = (r_node).Y0();
      (r_node).Z() = (r_node).Z0();
    }
  }

  void UpdateReferenceMesh(const ModelPart::NodesContainerType& rNodes) {
    for (auto& r_node : rNodes) {
      (r_node).X0() = (r_node).X();
      (r_node).Y0() = (r_node).Y();
      (r_node).Z0() = (r_node).Z();
    }
  }


/*   void GenerateMeshPart(ModelPart::Pointer pMeshModelPart, const ModelPart::ElementsContainerType& rElements) {
    pMeshModelPart = ModelPart::Pointer(new ModelPart("MeshPart", 1));

    // initializing mesh nodes
    //pMeshModelPart->Nodes() = BaseType::GetModelPart().Nodes();

    // creating mesh elements
    ModelPart::ElementsContainerType &MeshElems =
        pMeshModelPart->Elements();
    Element::Pointer pElem;

    for (ModelPart::ElementsContainerType::iterator it =
             rElements.ElementsBegin();
         it != rElements.ElementsEnd(); ++it) {

      pElem = Kratos::make_shared<LaplacianMeshMovingElement>(
          (*it).Id(), (*it).pGetGeometry(), (*it).pGetProperties());
      MeshElems.push_back(pElem);
    }
  } */
} // namespace Move Mesh Utilities.

} // namespace Kratos.