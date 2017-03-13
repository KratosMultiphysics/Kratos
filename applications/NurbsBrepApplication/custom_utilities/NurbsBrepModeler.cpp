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
      //std::cout << "number of faces: " << m_brep_model_vector[brep_itr].GetFaceVector().size() << std::endl;
      //std::cout << "number of edges: " << m_brep_model_vector[brep_itr].GetEdgeVector().size() << std::endl;

      for (unsigned int face_itr = 0; face_itr < m_brep_model_vector[brep_itr].GetFaceVector().size(); face_itr++)
      {
        Face& face = m_brep_model_vector[brep_itr].GetFaceVector()[face_itr];
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

            if (point_is_inside)
            {
              // compute unique point in CAD-model for given u&v
              ++cad_node_counter;
              Point<3> cad_point_coordinates;

              std::cout << "check point inside" << std::endl;

              face.EvaluateSurfacePoint(cad_point_coordinates, u_i, v_j, m_model_part);

              std::cout << "check: " << std::endl;

              // Add id to point --> node. Add node to list of CAD nodes
              Node<3>::Pointer node = model_part.CreateNewNode(cad_node_counter, cad_point_coordinates[0], cad_point_coordinates[1], cad_point_coordinates[2]);
              node->SetValue(LOCAL_PARAMETERS, point_of_interest);
              node->SetValue(FACE_BREP_ID, face_id);
            }
          }
        }
        std::cout << "Ende" << std::endl;
      }
    }
    unsigned int nodes = model_part.NumberOfNodes();
  }

  Face& NurbsBrepModeler::GetFace(const unsigned int face_id)
  {
    for (unsigned int i = 0; i < m_brep_model_vector.size(); i++)
    {
      FacesVector& face_vector = m_brep_model_vector[i].GetFaceVector();
      for (unsigned int j = 0; j < face_vector.size(); j++)
      {
        if (face_vector[j].GetId() == face_id)
        {
          return face_vector[j];
        }
      }
    }
    KRATOS_THROW_ERROR(std::logic_error, "NurbsBrepModeler::GetFace: face_id does not exist", "");
  }

  //Tree< KDTreePartition<BucketType> > NurbsBrepModeler::CreateSearchTree(ModelPart model_part)
  //{
  //  std::cout << "\n> Starting construction of search-tree..." << std::endl;
  //  //boost::timer timer;
  //  typedef Bucket< 3, NodeType, NodeVector, NodeType::Pointer, std::vector<NodeType::Pointer>::iterator, std::vector<double>::iterator > BucketType;
  //  //std::vector<NodeType::Pointer>& node_vector = model_part.NodesArray;
  //  int bucket_size = 20;
  //  Tree< KDTreePartition<BucketType> > search_tree(model_part.NodesArray.NodesStart, model_part.NodesArray.NodesEnd, bucket_size);
  //  return search_tree;
  //  //std::cout << "> Time needed for constructing search-tree: " << timer.elapsed() << " s" << std::endl;
  //}

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
    NodeType::Pointer nearest_point = search_tree.SearchNearestPoint(*node);
    //std::cout << "Point found!" << std::endl;
    Vector local_parameters_of_nearest_point = nearest_point->GetValue(LOCAL_PARAMETERS);
    unsigned int face_id_of_nearest_point = nearest_point->GetValue(FACE_BREP_ID);
    std::cout << "face_id_of_nearest_point: " << face_id_of_nearest_point << std::endl;
    std::cout << "local_parameters_of_nearest_point: " << local_parameters_of_nearest_point << std::endl;

    Face& face = GetFace(face_id_of_nearest_point);

    face.MapNodeNewtonRaphson(nearest_point, node_on_geometry, m_model_part);
  }

NurbsBrepModeler::NurbsBrepModeler(BrepModelGeometryReader& brep_model_geometry_reader, ModelPart& model_part)
  : m_model_part(model_part)
{
  m_brep_model_vector = brep_model_geometry_reader.ReadGeometry(m_model_part);

  Node<3>::Pointer node = Node<3>::Pointer(new Node<3>(0, 1, 0));
  Node<3>::Pointer closest_point = Node<3>::Pointer(new Node<3>(0, 1, 0));

  MapNode(node, closest_point);
}



NurbsBrepModeler::~NurbsBrepModeler()
{}

}  // namespace Kratos.


