//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Riccardo Rossi
//

// System includes


// External includes 


// Project includes
#include "NurbsBrepModeler.h"
#include "utilities/math_utils.h"
#include "../../kratos/includes/node.h"

//#include "nurbs_utilities.h"

#include "BrepModelGeometryReader.h"


namespace Kratos
{

//TODO: for testing in c++ use the testing mechanism of kratos
  void NurbsBrepModeler::CreateMeshedPoints(ModelPart& model_part)
  {
    //std::cout << "start createMeshedPoints: number of breps" << m_brep_model_vector.size() <<std::endl;
    unsigned int u_resolution = 10;
    unsigned int v_resolution = 10;

    unsigned int cad_node_counter = 0;

    for (unsigned int brep_itr = 0; brep_itr < m_brep_model_vector.size(); brep_itr++)
    {

      for (unsigned int face_itr = 0; face_itr < m_brep_model_vector[brep_itr].GetFaceVector().size(); face_itr++)
      {
        BrepFace& face = m_brep_model_vector[brep_itr].GetFaceVector()[face_itr];
        unsigned int face_id = face.GetId();

        Vector& knot_vector_u = face.GetUKnotVector();
        Vector& knot_vector_v = face.GetVKnotVector();

        double u_min = knot_vector_u[0];
        double u_max = knot_vector_u[knot_vector_u.size()-1];
        double v_min = knot_vector_v[0];
        double v_max = knot_vector_v[knot_vector_v.size()-1];

        double delta_u = (u_max - u_min) / u_resolution;
        double delta_v = (v_max - v_min) / v_resolution;

        // Loop over all u & v according to specified resolution
        for (unsigned int i = 1; i < u_resolution; i++)
        {
          // current u-value
          double u_i = u_min + i*delta_u;

          for (unsigned int j = 1; j < v_resolution; j++)
          {
            // current v-value
            double v_j = v_min + j*delta_v;

            // Check if u_i and v_j represent a point inside the closed boundary loop
            Vector point_of_interest = ZeroVector(2);
            point_of_interest[0] = u_i;
            point_of_interest[1] = v_j;

            bool point_is_inside = face.CheckIfPointIsInside(point_of_interest);

            //if (point_is_inside)
            //{
            //  // compute unique point in CAD-model for given u&v
            //  ++cad_node_counter;
            //  Node<3>::Pointer cad_point = Node<3>::Pointer(new Node<3>(cad_node_counter, 0,0,0));

            //  face.EvaluateSurfacePoint(cad_point, u_i, v_j);

            //  cad_point->SetValue(LOCAL_PARAMETERS, point_of_interest);
            //  cad_point->SetValue(FACE_BREP_ID, face_id);
            //  // Add id to point --> node. Add node to list of CAD nodes
            //  model_part.AddNode(cad_point);// .CreateNewNode(cad_node_counter, cad_point_coordinates[0], cad_point_coordinates[1], cad_point_coordinates[2]);
            //}
          }
        }
      }
    }
    unsigned int nodes = model_part.NumberOfNodes();
  }

  void NurbsBrepModeler::CreateIntegrationDomain(const int& shapefunction_order, ModelPart& model_part)
  {
    int id_itr = 1;
    //model_part.CreateSubModelPart("FACES");
    //model_part.CreateSubModelPart("COUPLING_EDGES");
    //model_part.CreateSubModelPart("EDGES");
    ModelPart::Pointer model_part_faces = model_part.CreateSubModelPart("FACES");
    ModelPart::Pointer model_part_coupling_edges = model_part.CreateSubModelPart("COUPLING_EDGES");
    ModelPart::Pointer model_part_edges = model_part.CreateSubModelPart("EDGES");

    for (unsigned int brep_itr = 0; brep_itr < m_brep_model_vector.size(); brep_itr++)
    {
      for (unsigned int face_itr = 0; face_itr < m_brep_model_vector[brep_itr].GetFaceVector().size(); face_itr++)
      {
        BrepFace& face = m_brep_model_vector[brep_itr].GetFaceVector()[face_itr];

        //model_part_faces.CreateSubModelPart("FACE_" + std::to_string(face.Id()));
        ModelPart::Pointer model_part_face_id = model_part_faces->CreateSubModelPart("FACE_" + std::to_string(face.Id()));

        std::vector<Node<3>::Pointer> NodeVectorElement = face.GetQuadraturePointsTrimmed(shapefunction_order);
        for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
        {
          NodeVectorElement[k]->SetId(id_itr);
          id_itr++;
          model_part_face_id->AddNode(NodeVectorElement[k]);
        }
        //model_part_faces.CreateSubModelPart("FACE_" + std::to_string(face.Id()) + "_EMBEDDED");
        ModelPart::Pointer model_part_face_id_embedded = model_part_faces->CreateSubModelPart("FACE_" + std::to_string(face.Id()) + "_EMBEDDED");
        std::vector<Node<3>::Pointer> NodeVectorEmbeddedElement = face.GetQuadraturePointsEmbedded(shapefunction_order);
        for (unsigned int k = 0; k < NodeVectorEmbeddedElement.size(); k++)
        {
          NodeVectorEmbeddedElement[k]->SetId(id_itr);
          id_itr++;
          model_part_face_id_embedded->AddNode(NodeVectorEmbeddedElement[k]);
        }
      }

      for (unsigned int edge_itr = 0; edge_itr < m_brep_model_vector[brep_itr].GetEdgeVector().size(); edge_itr++)
      {
        BrepEdge& edge = m_brep_model_vector[brep_itr].GetEdgeVector()[edge_itr];

        if (edge.isCouplingEdge())
        {
          //model_part_coupling_edges.CreateSubModelPart("COUPLING_EDGE_" + std::to_string(edge.Id()));
          ModelPart::Pointer model_part_coupling_edge_id = model_part_coupling_edges->CreateSubModelPart("COUPLING_EDGE_" + std::to_string(edge.Id()));
          int face_id_slave;
          int trim_index_slave;
          edge.GetEdgeInformation(1, face_id_slave, trim_index_slave);
          BrepFace& face_slave = GetFace(face_id_slave);
          std::vector<Point> points = face_slave.GetIntersectionPoints(trim_index_slave);

          int face_id_master;
          int trim_index_master;
          edge.GetEdgeInformation(0, face_id_master, trim_index_master);
          BrepFace& face_master = GetFace(face_id_master);

          std::vector<Node<3>::Pointer> NodeVectorElement = face_master.GetQuadraturePointsOfTrimmingCurveWithPoints(
            shapefunction_order, trim_index_master, points);

          //for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
          //{
          //  KRATOS_WATCH(NodeVectorElement[k]->GetValue(SHAPE_FUNCTION_VALUES))
          //  KRATOS_WATCH(NodeVectorElement[k]->GetValue(SHAPE_FUNCTION_SLAVE))
          //}

          face_slave.EnhanceShapeFunctionsSlave(NodeVectorElement, trim_index_slave, 2);
          for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
          {
            //KRATOS_WATCH(NodeVectorElement[k]->GetValue(SHAPE_FUNCTION_VALUES))
            //KRATOS_WATCH(NodeVectorElement[k]->GetValue(SHAPE_FUNCTION_SLAVE))
            NodeVectorElement[k]->SetId(id_itr);
            id_itr++;
            model_part_coupling_edge_id->AddNode(NodeVectorElement[k]);
          }
        }
        else
        {
          //model_part_edges.CreateSubModelPart("EDGE_" + std::to_string(edge.Id()));
          ModelPart::Pointer model_part_edge_id = model_part_edges->CreateSubModelPart("EDGE_" + std::to_string(edge.Id()));
          int face_id;
          int trim_index;
          edge.GetEdgeInformation(0, face_id, trim_index);

          BrepFace& face = GetFace(face_id);

          std::vector<Node<3>::Pointer> NodeVectorElement = face.GetQuadraturePointsOfTrimmingCurve(shapefunction_order, trim_index);
          for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
          {
            NodeVectorElement[k]->SetId(id_itr);
            id_itr++;
            model_part_edge_id->AddNode(NodeVectorElement[k]);
          }
        }
      }
    }
  }

  BrepFace& NurbsBrepModeler::GetFace(const unsigned int face_id)
  {
    for (unsigned int i = 0; i < m_brep_model_vector.size(); i++)
    {
      BrepFacesVector& face_vector = m_brep_model_vector[i].GetFaceVector();
      for (unsigned int j = 0; j < face_vector.size(); j++)
      {
        if (face_vector[j].GetId() == face_id)
        {
          return face_vector[j];
        }
      }
    }
    KRATOS_THROW_ERROR(std::logic_error, "NurbsBrepModeler::GetFace: face with face_id " + std::to_string(face_id) + " does not exist", "");
  }
  //void NurbsBrepModeler::GetClosestPoint(const Node<3>::Pointer& node, Node<3>::Pointer& node_on_geometry)
  void NurbsBrepModeler::MapNode(const Node<3>::Pointer& node, Node<3>::Pointer& node_on_geometry)
  {
    ModelPart local_model_part("MeshingPoints");

    CreateMeshedPoints(local_model_part);

    std::vector<Node<3>::Pointer> list_of_nodes;
    list_of_nodes.reserve(local_model_part.NumberOfNodes());
    //unsigned int number_of_nodes2 = local_model_part.NumberOfNodes();
    for (ModelPart::NodesContainerType::iterator i_node = local_model_part.NodesBegin(); i_node != local_model_part.NodesEnd(); i_node++)
    {
      (list_of_nodes).push_back(*(i_node.base()));
    }
    //KRATOS_WATCH(list_of_nodes.size())
    //  unsigned int number_of_nodes = list_of_nodes.size();
    //boost::shared_ptr<Node<3>>& nodes_local = local_model_part.NodesArray;
    int bucket_size = 20;
    tree search_tree(list_of_nodes.begin(), list_of_nodes.end(), bucket_size);
    node_on_geometry = search_tree.SearchNearestPoint(*node);
    //std::cout << "Point found!" << std::endl;
    Vector local_parameters_of_nearest_point = node_on_geometry->GetValue(LOCAL_PARAMETERS);
    unsigned int face_id_of_nearest_point = node_on_geometry->GetValue(FACE_BREP_ID);
    std::cout << "face_id_of_nearest_point: " << face_id_of_nearest_point << std::endl;
    std::cout << "local_parameters_of_nearest_point: " << local_parameters_of_nearest_point << std::endl;

    BrepFace& face = GetFace(face_id_of_nearest_point);

    face.MapNodeNewtonRaphson(node, node_on_geometry);

    local_parameters_of_nearest_point = node_on_geometry->GetValue(LOCAL_PARAMETERS);
    face_id_of_nearest_point = node_on_geometry->GetValue(FACE_BREP_ID);
    std::cout << "face_id_of_nearest_point: " << face_id_of_nearest_point << std::endl;
    std::cout << "local_parameters_of_nearest_point: " << local_parameters_of_nearest_point << std::endl;

    Vector N = node_on_geometry->GetValue(SHAPE_FUNCTION_VALUES);
    KRATOS_WATCH(N)
  }

NurbsBrepModeler::NurbsBrepModeler(BrepModelGeometryReader& brep_model_geometry_reader, ModelPart& model_part)
  : m_model_part(model_part)
{
  m_brep_model_vector = brep_model_geometry_reader.ReadGeometry(m_model_part);

  //Node<3>::Pointer node = Node<3>::Pointer(new Node<3>(0, 0, 7.3, 0));
  //Node<3>::Pointer closest_point = Node<3>::Pointer(new Node<3>(0, 0.25, 0.1792, 0));

  //MapNode(node, closest_point);
}



NurbsBrepModeler::~NurbsBrepModeler()
{}

}  // namespace Kratos.


