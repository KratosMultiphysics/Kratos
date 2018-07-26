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
#include "includes/node.h"

#include "spatial_containers/bins_dynamic_objects.h"
#include "custom_search/bins_iga_configure.h"

//#include "nurbs_utilities.h"

#include "BrepModelGeometryReader.h"


namespace Kratos
{
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

	void NurbsBrepModeler::ApplyGeometryRefinement(Parameters& rRefinementParameters)
	{
		// for (auto refinement = rRefinementParameters.begin(); refinement != rRefinementParameters.end(); ++refinement)
		for (int refinement_i = 0; refinement_i < rRefinementParameters.size(); refinement_i++)
		{
			Parameters refinement_parameters = rRefinementParameters[refinement_i]["parameters"];
			// Parameter description
			BrepFace::GeometryRefinementParameters geometrical_refinement_parameter;
			if (refinement_parameters.Has("knot_insertions_u"))
			{
				unsigned int length_knot_insertions_u = refinement_parameters["knot_insertions_u"].size();
				geometrical_refinement_parameter.knot_insertions_u = ZeroVector(length_knot_insertions_u);
				for (int i = 0; i < length_knot_insertions_u; i++)
					geometrical_refinement_parameter.knot_insertions_u[i] = refinement_parameters["knot_insertions_u"][i].GetDouble();

			}
			if (refinement_parameters.Has("knot_insertions_v"))
			{
				unsigned int length_knot_insertions_v = refinement_parameters["knot_insertions_v"].size();
				geometrical_refinement_parameter.knot_insertions_v = ZeroVector(length_knot_insertions_v);
				for (int i = 0; i < length_knot_insertions_v; i++)
					geometrical_refinement_parameter.knot_insertions_v[i] = refinement_parameters["knot_insertions_v"][i].GetDouble();
			}

			if (refinement_parameters.Has("min_elements_u"))
				geometrical_refinement_parameter.multiply_knots_u = refinement_parameters["min_elements_u"].GetInt();
			if (refinement_parameters.Has("min_elements_v"))
				geometrical_refinement_parameter.multiply_knots_v = refinement_parameters["min_elements_v"].GetInt();

			if (refinement_parameters.Has("multiply_knots_u"))
				geometrical_refinement_parameter.multiply_knots_u = refinement_parameters["multiply_knots_u"].GetDouble();
			if (refinement_parameters.Has("multiply_knots_v"))
				geometrical_refinement_parameter.multiply_knots_v = refinement_parameters["multiply_knots_v"].GetDouble();

			if (refinement_parameters.Has("min_order_p"))
				geometrical_refinement_parameter.min_order_p = refinement_parameters["min_order_p"].GetInt();
			if (refinement_parameters.Has("min_order_q"))
				geometrical_refinement_parameter.min_order_q = refinement_parameters["min_order_q"].GetInt();

			if (refinement_parameters.Has("order_elevation_p"))
				geometrical_refinement_parameter.order_elevation_p = refinement_parameters["order_elevation_p"].GetInt();
			if (refinement_parameters.Has("order_elevation_q"))
				geometrical_refinement_parameter.order_elevation_q = refinement_parameters["order_elevation_q"].GetInt();

			if (rRefinementParameters[refinement_i]["selection"].GetString() == "ALL")
			{
				for (unsigned int brep_itr = 0; brep_itr < m_brep_model_vector.size(); brep_itr++)
				{
					for (unsigned int face_itr = 0; face_itr < m_brep_model_vector[brep_itr].GetFaceVector().size(); face_itr++)
					{
						BrepFace& face = m_brep_model_vector[brep_itr].GetFaceVector()[face_itr];
						face.ApplyGeometryRefinement(geometrical_refinement_parameter);
					}
				}
			}
		}
	}

	void NurbsBrepModeler::CreateIntegrationDomain(Parameters& rIntegrationDomainParameter, ModelPart& rIntegrationDomainModelPart)
	{
		std::cout << "> Start create integration domain" << std::endl;

		const double accuracy = rIntegrationDomainParameter["accuracy"].GetDouble();
		const int shapefunction_order = rIntegrationDomainParameter["shape_function_derivatives_order"].GetInt();
		const int number_projection_iterations = rIntegrationDomainParameter["number_projection_iterations"].GetInt();
		const int polygon_discretization = rIntegrationDomainParameter["polygon_discretization"].GetInt();

		const bool faces = rIntegrationDomainParameter["integration_domains"]["faces"].GetBool();
		const bool faces_embedded = rIntegrationDomainParameter["integration_domains"]["faces_embedded"].GetBool();
		const bool faces_reversed = rIntegrationDomainParameter["integration_domains"]["faces_reversed"].GetBool();
		const bool edges = rIntegrationDomainParameter["integration_domains"]["edges"].GetBool();
		const bool coupling_edges = rIntegrationDomainParameter["integration_domains"]["coupling_edges"].GetBool();

		int id_itr = 1;
		ModelPart::Pointer model_part_faces = rIntegrationDomainModelPart.CreateSubModelPart("FACES");
		ModelPart::Pointer model_part_coupling_edges = rIntegrationDomainModelPart.CreateSubModelPart("COUPLING_EDGES");
		ModelPart::Pointer model_part_edges = rIntegrationDomainModelPart.CreateSubModelPart("EDGES");

		for (unsigned int brep_itr = 0; brep_itr < m_brep_model_vector.size(); brep_itr++)
		{
			for (unsigned int face_itr = 0; face_itr < m_brep_model_vector[brep_itr].GetFaceVector().size(); face_itr++)
			{
				BrepFace& face = m_brep_model_vector[brep_itr].GetFaceVector()[face_itr];
				if (faces)
				{
					ModelPart::Pointer model_part_face_id = model_part_faces->CreateSubModelPart("FACE_" + std::to_string(face.Id()));

					std::vector<Node<3>::Pointer> NodeVectorElement = face.GetIntegrationNodesSurface(shapefunction_order, polygon_discretization);
					for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
					{
						NodeVectorElement[k]->SetId(id_itr);
						id_itr++;
						model_part_face_id->AddNode(NodeVectorElement[k]);
					}
				}
				if (faces_embedded)
				{
					ModelPart::Pointer model_part_face_id_embedded = model_part_faces->CreateSubModelPart("FACE_" + std::to_string(face.Id()) + "_EMBEDDED");
					std::vector<Node<3>::Pointer> NodeVectorEmbeddedElement = face.GetIntegrationNodesEmbedded(shapefunction_order, polygon_discretization);
					for (unsigned int k = 0; k < NodeVectorEmbeddedElement.size(); k++)
					{
						NodeVectorEmbeddedElement[k]->SetId(id_itr);
						id_itr++;
						model_part_face_id_embedded->AddNode(NodeVectorEmbeddedElement[k]);
					}
				}
				if (faces_reversed)
				{
					ModelPart::Pointer model_part_face_id_reversed = model_part_faces->CreateSubModelPart("FACE_" + std::to_string(face.Id()) + "_REVERSE");
					std::vector<Node<3>::Pointer> NodeVectorReversedElement = face.GetIntegrationNodesSurfaceReversed(shapefunction_order, polygon_discretization);
					for (unsigned int k = 0; k < NodeVectorReversedElement.size(); k++)
					{
						NodeVectorReversedElement[k]->SetId(id_itr);
						id_itr++;
						model_part_face_id_reversed->AddNode(NodeVectorReversedElement[k]);
					}
				}
			}

			for (unsigned int edge_itr = 0; edge_itr < m_brep_model_vector[brep_itr].GetEdgeVector().size(); edge_itr++)
			{
				BrepEdge& edge = m_brep_model_vector[brep_itr].GetEdgeVector()[edge_itr];

				if (edge.IsCouplingEdge())
				{
					if (coupling_edges)
					{
						ModelPart::Pointer model_part_coupling_edge_id = model_part_coupling_edges->CreateSubModelPart("COUPLING_EDGE_" + std::to_string(edge.Id()));
						BrepEdge::Topology slave_topology = edge.GetEdgeInformation(1);
						BrepFace& face_slave = GetFace(slave_topology.face_id);
						std::vector<Point> points = face_slave.GetKnotIntersections(slave_topology.trim_index);

						int slave_p_q = face_slave.GetP() + face_slave.GetQ() + 1;

						BrepEdge::Topology master_topology = edge.GetEdgeInformation(0);
						BrepFace& face_master = GetFace(master_topology.face_id);

						std::vector<Node<3>::Pointer> NodeVectorElement = face_master.GetIntegrationNodesTrimmingCurveMaster(
							points, shapefunction_order, master_topology.trim_index, slave_p_q, accuracy, m_model_tolerance, number_projection_iterations);

						face_slave.EvaluateIntegrationNodesTrimmingCurveSlave(NodeVectorElement, shapefunction_order, slave_topology.trim_index, accuracy, m_model_tolerance, number_projection_iterations);
						for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
						{
							NodeVectorElement[k]->SetId(id_itr);
							id_itr++;
							model_part_coupling_edge_id->AddNode(NodeVectorElement[k]);
						}
					}
				}
				else
				{
					if (edges)
					{
						ModelPart::Pointer model_part_edge_id = model_part_edges->CreateSubModelPart("EDGE_" + std::to_string(edge.Id()));
						BrepEdge::Topology edge_topology = edge.GetEdgeInformation(0);

						BrepFace& face = GetFace(edge_topology.face_id);

						std::vector<Node<3>::Pointer> NodeVectorElement = face.GetIntegrationNodesTrimmingCurve(shapefunction_order, edge_topology.trim_index, accuracy);
						for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
						{
							NodeVectorElement[k]->SetId(id_itr);
							id_itr++;
							model_part_edge_id->AddNode(NodeVectorElement[k]);
						}
					}
				}
			}

			for (unsigned int vertex_itr = 0; vertex_itr < m_brep_model_vector[brep_itr].GetVertexVector().size(); vertex_itr++)
			{
				ModelPart::Pointer model_part_points;
				if (!rIntegrationDomainModelPart.HasSubModelPart("POINTS"))
					model_part_points = rIntegrationDomainModelPart.CreateSubModelPart("POINTS");
				else
					model_part_points = rIntegrationDomainModelPart.pGetSubModelPart("POINTS");

				BrepVertex& vertex = m_brep_model_vector[brep_itr].GetVertexVector()[vertex_itr];
				BrepVertex::Topology vertex_topology = vertex.GetVertexInformation(0);

				ModelPart::Pointer model_part_point_id = model_part_points->CreateSubModelPart("POINT_" + std::to_string(vertex.Id()));


				BrepFace& face = GetFace(vertex_topology.face_id);
				Node<3>::Pointer node = face.GetIntegrationNodePoint(vertex_topology.trim_index, shapefunction_order);
				node->SetId(id_itr);
				model_part_point_id->AddNode(node);
			}
		}
		std::cout << "> Finished create integration domain" << std::endl;
	}

	void NurbsBrepModeler::CreateIntegrationDomainElements(Parameters& rIntegrationDomainParameter, ModelPart& rModelPart)
	{
		std::cout << "> Start create integration domain" << std::endl;

		const double accuracy = rIntegrationDomainParameter["accuracy"].GetDouble();
		const int shapefunction_order = rIntegrationDomainParameter["shape_function_derivatives_order"].GetInt();
		const int number_projection_iterations = rIntegrationDomainParameter["number_projection_iterations"].GetInt();
		const int polygon_discretization = rIntegrationDomainParameter["polygon_discretization"].GetInt();

		const bool faces = rIntegrationDomainParameter["integration_domains"]["faces"].GetBool();
		const bool faces_embedded = rIntegrationDomainParameter["integration_domains"]["faces_embedded"].GetBool();
		const bool faces_reversed = rIntegrationDomainParameter["integration_domains"]["faces_reversed"].GetBool();
		const bool edges = rIntegrationDomainParameter["integration_domains"]["edges"].GetBool();
		const bool coupling_edges = rIntegrationDomainParameter["integration_domains"]["coupling_edges"].GetBool();

		int id_itr = 1;
		ModelPart::Pointer model_part_faces = rModelPart.CreateSubModelPart("FACES");
		ModelPart::Pointer model_part_coupling_edges = rModelPart.CreateSubModelPart("COUPLING_EDGES");
		ModelPart::Pointer model_part_edges = rModelPart.CreateSubModelPart("EDGES");

		for (unsigned int brep_itr = 0; brep_itr < m_brep_model_vector.size(); brep_itr++)
		{
			for (unsigned int face_itr = 0; face_itr < m_brep_model_vector[brep_itr].GetFaceVector().size(); face_itr++)
			{
				BrepFace& face = m_brep_model_vector[brep_itr].GetFaceVector()[face_itr];

				if (faces)
				{
					ModelPart::Pointer model_part_face_id = model_part_faces->CreateSubModelPart("FACE_" + std::to_string(face.Id()));

					std::vector<NodeType::Pointer> NodeVectorElement = face.GetIntegrationNodesSurface(shapefunction_order, polygon_discretization);
					for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
					{
						Vector node_ids = NodeVectorElement[k]->GetValue(CONTROL_POINT_IDS);
						std::vector<NodeType::Pointer> node_vector;
						for (int i = 0; i < node_ids.size(); i++)
						{
							node_vector.push_back(rModelPart.pGetNode(node_ids[i]));
							model_part_faces->AddNode(rModelPart.pGetNode(node_ids[i]));
						}
						//model_part_faces->CreateNewElement("Element", 1, {1,2},0);
						//NodeVectorElement[k]->SetId(id_itr);
						//id_itr++;
						//model_part_face_id->AddNode(NodeVectorElement[k]);
					}
				}
				//if (faces_embedded)
				//{
				//	ModelPart::Pointer model_part_face_id_embedded = model_part_faces->CreateSubModelPart("FACE_" + std::to_string(face.Id()) + "_EMBEDDED");
				//	std::vector<Node<3>::Pointer> NodeVectorEmbeddedElement = face.GetIntegrationNodesEmbedded(shapefunction_order, polygon_discretization);
				//	for (unsigned int k = 0; k < NodeVectorEmbeddedElement.size(); k++)
				//	{
				//		NodeVectorEmbeddedElement[k]->SetId(id_itr);
				//		id_itr++;
				//		model_part_face_id_embedded->AddNode(NodeVectorEmbeddedElement[k]);
				//	}
				//}
				//if (faces_reversed)
				//{
				//	ModelPart::Pointer model_part_face_id_reversed = model_part_faces->CreateSubModelPart("FACE_" + std::to_string(face.Id()) + "_REVERSE");
				//	std::vector<Node<3>::Pointer> NodeVectorReversedElement = face.GetIntegrationNodesSurfaceReversed(shapefunction_order, polygon_discretization);
				//	for (unsigned int k = 0; k < NodeVectorReversedElement.size(); k++)
				//	{
				//		NodeVectorReversedElement[k]->SetId(id_itr);
				//		id_itr++;
				//		model_part_face_id_reversed->AddNode(NodeVectorReversedElement[k]);
				//	}
				//}
			}

			//for (unsigned int edge_itr = 0; edge_itr < m_brep_model_vector[brep_itr].GetEdgeVector().size(); edge_itr++)
			//{
			//	BrepEdge& edge = m_brep_model_vector[brep_itr].GetEdgeVector()[edge_itr];

			//	if (edge.IsCouplingEdge())
			//	{
			//		if (coupling_edges)
			//		{
			//			ModelPart::Pointer model_part_coupling_edge_id = model_part_coupling_edges->CreateSubModelPart("COUPLING_EDGE_" + std::to_string(edge.Id()));
			//			BrepEdge::Topology slave_topology = edge.GetEdgeInformation(1);
			//			BrepFace& face_slave = GetFace(slave_topology.face_id);
			//			std::vector<Point> points = face_slave.GetKnotIntersections(slave_topology.trim_index);

			//			int slave_p_q = face_slave.GetP() + face_slave.GetQ() + 1;

			//			BrepEdge::Topology master_topology = edge.GetEdgeInformation(0);
			//			BrepFace& face_master = GetFace(master_topology.face_id);


			//			std::vector<Node<3>::Pointer> NodeVectorElement = face_master.GetIntegrationNodesTrimmingCurveMaster(
			//				points, shapefunction_order, master_topology.trim_index, slave_p_q, accuracy, m_model_tolerance, number_projection_iterations);

			//			face_slave.EvaluateIntegrationNodesTrimmingCurveSlave(NodeVectorElement, shapefunction_order, slave_topology.trim_index, accuracy, m_model_tolerance, number_projection_iterations);
			//			for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
			//			{
			//				NodeVectorElement[k]->SetId(id_itr);
			//				id_itr++;
			//				model_part_coupling_edge_id->AddNode(NodeVectorElement[k]);
			//			}
			//		}
			//	}
			//	else
			//	{
			//		if (edges)
			//		{
			//			ModelPart::Pointer model_part_edge_id = model_part_edges->CreateSubModelPart("EDGE_" + std::to_string(edge.Id()));
			//			BrepEdge::Topology edge_topology = edge.GetEdgeInformation(0);

			//			BrepFace& face = GetFace(edge_topology.face_id);

			//			std::vector<Node<3>::Pointer> NodeVectorElement = face.GetIntegrationNodesTrimmingCurve(shapefunction_order, edge_topology.trim_index, accuracy);
			//			for (unsigned int k = 0; k < NodeVectorElement.size(); k++)
			//			{
			//				NodeVectorElement[k]->SetId(id_itr);
			//				id_itr++;
			//				model_part_edge_id->AddNode(NodeVectorElement[k]);
			//			}
			//		}
			//	}
			//}

			//for (unsigned int vertex_itr = 0; vertex_itr < m_brep_model_vector[brep_itr].GetVertexVector().size(); vertex_itr++)
			//{
			//	ModelPart::Pointer model_part_points;
			//	if (!rIntegrationDomainModelPart.HasSubModelPart("POINTS"))
			//		model_part_points = rIntegrationDomainModelPart.CreateSubModelPart("POINTS");
			//	else
			//		*model_part_points = rIntegrationDomainModelPart.GetSubModelPart("POINTS");

			//	BrepVertex& vertex = m_brep_model_vector[brep_itr].GetVertexVector()[vertex_itr];
			//	BrepVertex::Topology vertex_topology = vertex.GetVertexInformation(0);

			//	ModelPart::Pointer model_part_point_id = model_part_points->CreateSubModelPart("POINT_" + std::to_string(vertex.Id()));


			//	BrepFace& face = GetFace(vertex_topology.face_id);
			//	Node<3>::Pointer node = face.GetIntegrationNodePoint(vertex_topology.trim_index, shapefunction_order);
			//	node->SetId(id_itr);
			//	model_part_point_id->AddNode(node);
			//}
		}
		std::cout << "> Finished create integration domain" << std::endl;
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
		KRATOS_ERROR << "NurbsBrepModeler::GetFace: face with face_id " << std::to_string(face_id) << " does not exist" << std::endl;
	}

	void NurbsBrepModeler::MapNode(const Node<3>::Pointer& node, Node<3>::Pointer& node_on_geometry, ModelPart& rSearchModelPart)
	{
		std::vector<Node<3>::Pointer> list_of_nodes;
		list_of_nodes.reserve(rSearchModelPart.NumberOfNodes());
		for (ModelPart::NodesContainerType::iterator i_node = rSearchModelPart.NodesBegin(); i_node != rSearchModelPart.NodesEnd(); i_node++)
		{
			(list_of_nodes).push_back(*(i_node.base()));
		}
		int bucket_size = 20;
		tree search_tree(list_of_nodes.begin(), list_of_nodes.end(), bucket_size);
		node_on_geometry = search_tree.SearchNearestPoint(*node);
		//std::cout << "Point found!" << std::endl;
		Vector local_parameters_of_nearest_point = node_on_geometry->GetValue(LOCAL_PARAMETERS);
		unsigned int face_id_of_nearest_point = node_on_geometry->GetValue(FACE_BREP_ID);
		std::cout << "face_id_of_nearest_point: " << face_id_of_nearest_point << std::endl;
		std::cout << "local_parameters_of_nearest_point: " << local_parameters_of_nearest_point << std::endl;

		BrepFace& face = GetFace(face_id_of_nearest_point);

		face.GetClosestIntegrationNode(node_on_geometry, node, 2, 1e-7, 30);

		local_parameters_of_nearest_point = node_on_geometry->GetValue(LOCAL_PARAMETERS);
		face_id_of_nearest_point = node_on_geometry->GetValue(FACE_BREP_ID);
		std::cout << "face_id_of_nearest_point: " << face_id_of_nearest_point << std::endl;
		std::cout << "local_parameters_of_nearest_point: " << local_parameters_of_nearest_point << std::endl;

		Vector N = node_on_geometry->GetValue(NURBS_SHAPE_FUNCTIONS);
		KRATOS_WATCH(N)
	}

	void NurbsBrepModeler::GetInterfaceConditions(ModelPart& rParticleModelPart, ModelPart& rConditionModelPart, ModelPart& rSearchModelPart)
	{
		//std::cout << "Number of elements: " << rParticleModelPart.Elements().size() << std::endl;
        tree search_tree(list_of_nodes.begin(), list_of_nodes.end(), bucket_size);

        for iga_elements:
            BinsIGAObjs.push_back(BinsIgaObject(iga_elem_ptr));


        iga_bins_structure = BinsObjectDynamic<BinsIgaConfigure>(BinsIGAObjs->begin(), BinsIGAObjs->end();

        num_interface_obj_bin = BinsIGAObjs.size()

        BinsIgaConfigure::ResultContainerType neighbor_results(num_interface_obj_bin);
        std::vector<double> neighbor_distances(num_interface_obj_bin);

        auto particle_obj = Kratos::make_shared<BinsIgaObject>(array_1d<double, 3>(0.0));

            auto results_itr = neighbor_results.begin();
            auto distance_itr = neighbor_distances.begin();

            particle_obj->UpdateCoordinates(particle.GetGeometry[0].Cooridnates());

            const SizeType number_of_results = iga_bins_structure.SearchObjectsInRadius(
                particle_obj, search_radius, results_itr,
                distance_itr, num_interface_obj_bin);

            for i in number_of_results:
                particle.SetValue(results_itr[i]);


		for (auto element = rParticleModelPart.ElementsBegin(); element != rParticleModelPart.ElementsEnd(); element++)
		{
			//std::cout << "Check here for more information: X: " << element->GetGeometry()[0].X() << ", Y: " << element->GetGeometry()[0].Y() << ", Z: " << element->GetGeometry()[0].Z() << std::endl;
			std::vector<Condition*> new_conditions;
			std::vector<array_1d<double, 3>> new_elastic_forces;
			std::vector<array_1d<double, 3>> new_total_forces;
			bool check = true;
			if (check)
			{
				Vector friction = ZeroVector(3);
				//Vector local_parameters_of_nearest_point;
				//unsigned int face_id_of_nearest_point;
				Node<3>::Pointer node = Kratos::make_shared<Node<3>>(0,0,0,0); // ::Pointer(new Node<3>(0));
				Node<3>::Pointer node_on_geometry = Node<3>::Pointer(new Node<3>(0,0,0,0));

				std::vector<Node<3>::Pointer> list_of_nodes;
				list_of_nodes.reserve(rSearchModelPart.NumberOfNodes());
				for (ModelPart::NodesContainerType::iterator i_node = rSearchModelPart.NodesBegin(); i_node != rSearchModelPart.NodesEnd(); i_node++)
				{
					(list_of_nodes).push_back(*(i_node.base()));
				}
				const int bucket_size = 20;
				*node = element->GetGeometry()[0];
				node_on_geometry = search_tree.SearchNearestPoint(*node);
				//}
				//std::cout << "Closest point: X: " << node_on_geometry->X() << ", Y: " << node_on_geometry->Y() << ", Z: " << node_on_geometry->Z() << std::endl;
				double distance_to_closest_element = sqrt(pow(node->X() - node_on_geometry->X(), 2) + pow(node->Y() - node_on_geometry->Y(), 2) + pow(node->Z() - node_on_geometry->Z(), 2));
				//KRATOS_WATCH(distance_to_closest_element)
				double radius = element->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
				double search_radius = radius*4;
				if (distance_to_closest_element < search_radius)
				{
					unsigned int face_id_of_nearest_point = node_on_geometry->GetValue(FACE_BREP_ID);

					BrepFace& face = GetFace(face_id_of_nearest_point);

					//std::cout << "Closest point: X: " << node_on_geometry->X() << ", Y: " << node_on_geometry->Y() << ", Z: " << node_on_geometry->Z() << std::endl;

					face.GetClosestIntegrationNode(node_on_geometry, node, 2, 1e-7, 30);

					//KRATOS_WATCH("CLOSEST NODE FOUND")
					//KRATOS_WATCH(radius)
					std::cout << "Projected Point: " << node_on_geometry->X() << ", " << node_on_geometry->Y() << ", " << node_on_geometry->Z() << std::endl;
					std::cout << "External Point: " << node->X() << ", " << node->Y() << ", " << node->Z() << std::endl;
					//KRATOS_WATCH(node)
					//KRATOS_ERROR << "Here ends the game" << std::endl;

					double distance_radius = std::sqrt(std::pow(node->X() - node_on_geometry->X(), 2) +
					std::pow(node->Y() - node_on_geometry->Y(), 2) +
					std::pow(node->Z() - node_on_geometry->Z(), 2));

					//KRATOS_WATCH(distance_radius)
					//std::cout << "Projected point: X: " << node_on_geometry->X() << ", Y: " << node_on_geometry->Y() << ", Z: " << node_on_geometry->Z() << std::endl;
					if (distance_radius <= radius + 1e-7)
					{
						//std::cout << "output here!" << std::endl;
						Vector node_ids = node_on_geometry->GetValue(CONTROL_POINT_IDS);
						std::vector<std::size_t> node_ids_int(node_ids.size());
						for (int i = 0; i < node_ids.size(); i++)
						{
							node_ids_int[i] = static_cast<std::size_t>(node_ids[i]);
							ModelPart& mp = rConditionModelPart.GetRootModelPart();
							Node<3>::Pointer this_node = mp.pGetNode(node_ids_int[i]);
							rConditionModelPart.AddNode(this_node);
						}
						std::string condition_name = "MeshlessForceInterfaceCondition";
						int id = 1;
						//std::cout << "number of conditions: " << rConditionModelPart.Conditions().size() << std::endl;
						if (rConditionModelPart.Conditions().size()>0)
							id = rConditionModelPart.GetRootModelPart().Conditions().back().Id() + 1;
						rConditionModelPart.AddProperties(element->pGetProperties());

						//std::cout << "check create new conditions" << std::endl;
						Condition::Pointer cond = rConditionModelPart.CreateNewCondition(condition_name, id, node_ids_int, element->pGetProperties());
						Vector external_force_vector = ZeroVector(3);
						cond->SetValue(LOCAL_PARAMETERS, node_on_geometry->GetValue(LOCAL_PARAMETERS));
						cond->SetValue(FACE_BREP_ID, node_on_geometry->GetValue(FACE_BREP_ID));
						cond->SetValue(SHAPE_FUNCTION_VALUES, node_on_geometry->GetValue(NURBS_SHAPE_FUNCTIONS));
						cond->SetValue(EXTERNAL_FORCES_VECTOR, external_force_vector);

						new_conditions.push_back(&*cond);

						//std::cout << "check change of old conditions" << std::endl;
						int condition_length = element->GetValue(WALL_POINT_CONDITION_POINTERS).size();
						std::vector<Vector> coords;
						ProcessInfo emptyProcessInfo = ProcessInfo();
						bool success = false;
						for (int i = 0; i < condition_length; ++i)
						{
							//std::cout << "coordinates" << std::endl;
							element->GetValue(WALL_POINT_CONDITION_POINTERS)[i]->GetValueOnIntegrationPoints(COORDINATES, coords, emptyProcessInfo);
							//std::cout << "new condition on same place" << std::endl;
							//KRATOS_WATCH(coords[0])
							double distance_radius = std::sqrt(std::pow(coords[0][0] - node_on_geometry->X(), 2) +
								std::pow(coords[0][1] - node_on_geometry->Y(), 2) +
								std::pow(coords[0][2] - node_on_geometry->Z(), 2));
							if (distance_radius < radius / 3)
							{
								//std::cout << "new condition on same place" << std::endl;
								//KRATOS_WATCH(element->GetValue(WALL_POINT_CONDITION_ELASTIC_FORCES)[i])
								//KRATOS_WATCH(element->GetValue(WALL_POINT_CONDITION_TOTAL_FORCES)[i])
								new_elastic_forces.push_back(element->GetValue(WALL_POINT_CONDITION_ELASTIC_FORCES)[i]);
								new_total_forces.push_back(element->GetValue(WALL_POINT_CONDITION_TOTAL_FORCES)[i]);
								success = true;
							}
						}

						if (!success)
						{
							new_elastic_forces.push_back(ZeroVector(3));
							new_total_forces.push_back(ZeroVector(3));
						}
					}
				}
			}

			int condition_length = element->GetValue(WALL_POINT_CONDITION_POINTERS).size();
			std::vector<Vector> coords;
			ProcessInfo emptyProcessInfo = ProcessInfo();
			for (int i = 0; i < condition_length; ++i)
			{
				rConditionModelPart.RemoveConditionFromAllLevels(*(element->GetValue(WALL_POINT_CONDITION_POINTERS)[i]));
			}
			//std::cout << "check 1" << std::endl;
			element->SetValue(WALL_POINT_CONDITION_ELASTIC_FORCES, new_elastic_forces);
			element->SetValue(WALL_POINT_CONDITION_TOTAL_FORCES, new_total_forces);
			element->SetValue(WALL_POINT_CONDITION_POINTERS, new_conditions);

			//std::cout << "check 2" << std::endl;

		}
	}

	void NurbsBrepModeler::ComputeArea(ModelPart& rModelPart)
	{
		ModelPart& faces = rModelPart.GetSubModelPart("FACES");
		std::vector<std::string> sub_model_part_names = faces.GetSubModelPartNames();
		for (int i = 0; i < sub_model_part_names.size(); i++)
		{
			ModelPart& sub_model_part = faces.GetSubModelPart(sub_model_part_names[i]);
			double area = 0.0;
			for (auto node = sub_model_part.NodesBegin(); node != sub_model_part.NodesEnd(); ++node)
			{
				area += node->GetValue(INTEGRATION_WEIGHT);
			}
			std::cout << "> Area of sub model part " << sub_model_part_names[i] << ": " << area << std::endl;
		}
	}

	void NurbsBrepModeler::LoadGeometry(BrepModelGeometryReader& rBrepModelGeometryReader)
	{
		BrepModelVector brep_model_vector = rBrepModelGeometryReader.ReadGeometry(*mp_model_part);
		for (auto brep_model = brep_model_vector.begin(); brep_model != brep_model_vector.end(); brep_model++)
		{
			m_brep_model_vector.push_back(*brep_model);
		}

		m_model_tolerance = rBrepModelGeometryReader.ReadModelTolerance();
	}

	NurbsBrepModeler::NurbsBrepModeler(ModelPart::Pointer rpModelPart)
		: mp_model_part(rpModelPart)
	{
	}



NurbsBrepModeler::~NurbsBrepModeler()
{}

}  // namespace Kratos.


