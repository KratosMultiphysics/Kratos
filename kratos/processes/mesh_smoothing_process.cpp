//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                    
//
	           
// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/exception.h"
#include "includes/kratos_flags.h"
#include "processes/mesh_smoothing_process.h"
#include "processes/measure_mesh_quality_process.h"
#include "processes/find_nodal_neighbours_process.h"

namespace Kratos
{

	KRATOS_CREATE_LOCAL_FLAG(MeshSmoothingProcess, LAPLACIAN_SMOOTHING, 0); 
	KRATOS_CREATE_LOCAL_FLAG(MeshSmoothingProcess, EDGE_LENGTH_SMOOTHING, 1);
	KRATOS_CREATE_LOCAL_FLAG(MeshSmoothingProcess, MOVEMENT_SMOOTHING, 2);
	KRATOS_CREATE_LOCAL_FLAG(MeshSmoothingProcess, MULTI_LEVEL_SMOOTHING, 3);
	KRATOS_CREATE_LOCAL_FLAG(MeshSmoothingProcess, COARSE_MESH_NODE, 10);

	MeshSmoothingProcess::MeshSmoothingProcess(ModelPart& rModelPart, Flags Options, std::size_t IterationsNumber)
		:mrModelPart(rModelPart), mOptions(Options), mMaxIterationsNumber(IterationsNumber)
	{

	}
  
	MeshSmoothingProcess::~MeshSmoothingProcess()
	{

	}

	void MeshSmoothingProcess::Execute()
	{
		KRATOS_TRY

		PerformSmoothing();

		MeasureMeshQualityProcess measure_mesh_quality_process(mrModelPart, 2);
		measure_mesh_quality_process.Execute();

		std::cout << measure_mesh_quality_process << std::endl;

		KRATOS_CATCH("");
	}

	/// Turn back information as a string.
	std::string MeshSmoothingProcess::Info() const
	{
		return "MeshSmoothingProcess";
	}

	/// Print information about this object.
	void MeshSmoothingProcess::PrintInfo(std::ostream& rOStream) const
	{
		rOStream << Info();
	}

	/// Print object's data.
	void MeshSmoothingProcess::PrintData(std::ostream& rOStream) const
	{
		rOStream << Info();
	}

	void MeshSmoothingProcess::PerformSmoothing()
	{
		const double epsilon = 1e-12;
		if (mOptions.Is(MULTI_LEVEL_SMOOTHING))
		{
			ModelPart::MeshType& r_original_mesh = mrModelPart.GetMesh();
			ModelPart::MeshType& r_coars_mesh1 = mrModelPart.GetMesh(1);
			ModelPart::MeshType& r_coars_mesh2 = mrModelPart.GetMesh(2);
			if (r_coars_mesh1.NumberOfNodes() == 0)
			{
				PerformCoarsening(r_original_mesh, r_coars_mesh1);
				ModelPart coarse_model_part;
				coarse_model_part.GetMeshes().clear();
				coarse_model_part.GetMeshes().push_back(ModelPart::MeshType::Pointer(new ModelPart::MeshType(r_coars_mesh1)));
				FindNodalNeighboursProcess find_neighbour_process(coarse_model_part);
				find_neighbour_process.Execute();
				PerformCoarsening(r_coars_mesh1, r_coars_mesh2);
			}

			ModelPart coarse_model_part;
			coarse_model_part.GetMeshes().clear();
			coarse_model_part.GetMeshes().push_back(ModelPart::MeshType::Pointer(new ModelPart::MeshType(r_coars_mesh2)));
			FindNodalNeighboursProcess find_neighbour_process(coarse_model_part);
			find_neighbour_process.Execute();
			MeshSmoothingProcess coarse_smoothing_process(coarse_model_part, LAPLACIAN_SMOOTHING, 10);
			coarse_smoothing_process.Execute();

			KRATOS_WATCH(r_coars_mesh1)
			KRATOS_WATCH(r_coars_mesh2)

			return;
		}
		for (ModelPart::ElementIterator i_element = mrModelPart.ElementsBegin(); i_element != mrModelPart.ElementsEnd(); i_element++)
			if (i_element->GetGeometry().Area() < epsilon)
				for (Element::GeometryType::iterator i_node = i_element->GetGeometry().begin(); i_node != i_element->GetGeometry().end(); i_node++)
					i_node->Set(SELECTED);

		PointsVectorType optimal_positions;
		Vector weights;
		for (std::size_t i = 0; i < mMaxIterationsNumber; i++)
		{
			for (ModelPart::NodeIterator i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++)
			{
				NeighboursVectorType const& r_neighbours = i_node->GetValue(NEIGHBOUR_NODES);
				FindOptimumPositionsAndWeights(*i_node, r_neighbours, optimal_positions, weights);
				MoveNode(*i_node, r_neighbours, optimal_positions, weights);
			}
		}
	}

//		const double epsilon = 1e-12;
//		if (mOptions.Is(MOVEMENT_SMOOTHING))
//		{
//			std::cout << "Performing MOVEMENT_SMOOTHING" << std::endl;
//			for (std::size_t i = 0; i < 100 /*mMaxIterationsNumber*/; i++)
//			{
//				for (ModelPart::NodeIterator i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++)
//				{
//					NeighboursVectorType const& r_neighbours = i_node->GetValue(NEIGHBOUR_NODES);
//					const std::size_t size = r_neighbours.size();
//					if ((size > 0) && (i_node->IsNot(BOUNDARY)))
//					{
//						array_1d<double, 3> velocity = ZeroVector(3);
//						array_1d<double, 3> weight_sum = ZeroVector(3);
//						for (NeighboursVectorType::const_iterator i_neighbour_node = r_neighbours.begin(); i_neighbour_node != r_neighbours.end(); i_neighbour_node++)
//						{
//							array_1d<double, 3> connection = *i_neighbour_node - *i_node;
//							for (int i = 0; i < 3; i++)
//								connection[i] = std::fabs(connection[i]);
//							array_1d<double, 3> neighbour_velocity = i_neighbour_node->GetSolutionStepValue(VELOCITY);
//
//							double l = std::sqrt(inner_prod(connection, connection));
//
//							array_1d<double, 3> weight = connection / l;
//							weight += ScalarVector(3,0.1);
//							//double norm_v = std::sqrt(inner_prod(neighbour_velocity, neighbour_velocity));
//
//							//double weight = fabs(inner_prod(connection, connection));
//
//							//if (weight < epsilon)
//							//	weight = epsilon;
//
//							//if (norm_v > epsilon)
//							//	weight = fabs(inner_prod(connection, neighbour_velocity)) / weight;
//							//else
//							//	weight = 1.0 / weight;
//							for (int i = 0; i < 3; i++)
//							{
//								
//								velocity[i] += weight[i] * neighbour_velocity[i];
//							}
//							weight_sum += weight;
//						}
//						for (int i = 0; i < 3; i++)
//							i_node->GetSolutionStepValue(VELOCITY)[i] = velocity[i] / weight_sum[i];
//					}
//				}
//			}
//		}
//		for (ModelPart::ElementIterator i_element = mrModelPart.ElementsBegin(); i_element != mrModelPart.ElementsEnd(); i_element++)
//			if (i_element->GetGeometry().Area() < epsilon)
//				for (Element::GeometryType::iterator i_node = i_element->GetGeometry().begin(); i_node != i_element->GetGeometry().end(); i_node++)
//					i_node->Set(SELECTED);
//		
////		for (std::size_t i = 0; i < mMaxIterationsNumber; i++)
//		{
//			for (ModelPart::NodeIterator i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++)
//			{
//				NeighboursVectorType const& r_neighbours = i_node->GetValue(NEIGHBOUR_NODES);
//				MoveNode(*i_node, r_neighbours);
//			}
//		}
//	}


	void MeshSmoothingProcess::PerformCoarsening(ModelPart::MeshType& rOriginalMesh, ModelPart::MeshType& rCoarsMesh)
	{
		SelectCoarseMeshNodes(rOriginalMesh);
		CollapseNodes(rOriginalMesh);
		CreateCoarseMeshElements(rOriginalMesh, rCoarsMesh);

		for (ModelPart::NodeIterator i_node = rOriginalMesh.NodesBegin(); i_node != rOriginalMesh.NodesEnd(); i_node++)
		{
			if (i_node->Is(COARSE_MESH_NODE))
			{
				rCoarsMesh.AddNode(*(i_node.base()));
				i_node->GetSolutionStepValue(PRESSURE) = 100;
			}
		}

		//return;
		//// This algorithm can be return with less loop but has been return like this to be MPI compatible
		//// The VISITED flag is used 
		//int loop_counter = 0;
		//for (ModelPart::NodeIterator i_node = rOriginalMesh.NodesBegin(); i_node != rOriginalMesh.NodesEnd(); i_node++)
		//{
		//	i_node->GetValue(FATHER_NODES).clear();
		//	i_node->Set(VISITED, false);
		//	if (i_node->Is(BOUNDARY))
		//	{
		//		i_node->Set(VISITED);
		//		rCoarsMesh.AddNode(*(i_node.base()));
		//		i_node->GetValue(FATHER_NODES).push_back(*(i_node.base()));
		//		GlobalPointersVector< Node<3> >& r_neighbours = i_node->GetValue(NEIGHBOUR_NODES);
		//		for (GlobalPointersVector<Node<3> >::iterator i_neighbour_node = r_neighbours.begin(); i_neighbour_node != r_neighbours.end(); i_neighbour_node++)
		//		{
		//			loop_counter++;
		//			if ((i_neighbour_node->IsNot(VISITED)) && (i_neighbour_node->IsNot(BOUNDARY)))
		//			{
		//				i_neighbour_node->Set(VISITED);
		//				if (i_neighbour_node->GetValue(FATHER_NODES).empty())
		//					i_neighbour_node->GetValue(FATHER_NODES).push_back(*(i_node.base()));
		//				//i_neighbour_node->GetSolutionStepValue(PRESSURE) = cluster_number;
		//			}
		//		}
		//	}
		//}

		//int max_number_of_iterations = rOriginalMesh.NumberOfNodes();
		//bool is_finished = false;
		//int layer_category = 0;
		//int cluster_number = 0;

		//for (int i = 0; i < max_number_of_iterations; i++)
		//{
		//	KRATOS_WATCH(i)
		//		if (is_finished)
		//			break;
		//	is_finished = true;
		//	for (ModelPart::NodeIterator i_node = rOriginalMesh.NodesBegin(); i_node != rOriginalMesh.NodesEnd(); i_node++)
		//		if (i_node->Is(VISITED))
		//		{
		//			i_node->GetSolutionStepValue(PRESSURE) = 10.00; // cluster_number;
		//			GlobalPointersVector< Node<3> >& r_neighbours = i_node->GetValue(NEIGHBOUR_NODES);
		//			for (GlobalPointersVector<Node<3> >::iterator i_neighbour_node = r_neighbours.begin(); i_neighbour_node != r_neighbours.end(); i_neighbour_node++)
		//			{
		//				loop_counter++;
		//				if ((i_neighbour_node->IsNot(VISITED)))
		//				{
		//					is_finished = false;
		//					break;
		//				}
		//			}
		//		}
		//		else 
		//		{
		//			bool is_seed = true; // as a candidate
		//			GlobalPointersVector< Node<3> >& r_neighbours = i_node->GetValue(NEIGHBOUR_NODES);
		//			for (GlobalPointersVector<Node<3> >::iterator i_neighbour_node = r_neighbours.begin(); i_neighbour_node != r_neighbours.end(); i_neighbour_node++)
		//			{
		//				loop_counter++;
		//				if ((i_neighbour_node->Is(VISITED)))
		//				{
		//					is_seed = false;
		//					is_finished = false;
		//					break;
		//				}
		//			}
		//			if (is_seed)
		//			{
		//				i_node->GetSolutionStepValue(PRESSURE) = 100.00; // cluster_number;
		//				i_node->Set(VISITED);
		//				i_node->GetValue(FATHER_NODES).push_back(*(i_node.base()));
		//				rCoarsMesh.AddNode(*(i_node.base()));
		//				for (GlobalPointersVector<Node<3> >::iterator i_neighbour_node = r_neighbours.begin(); i_neighbour_node != r_neighbours.end(); i_neighbour_node++)
		//				{
		//					loop_counter++;
		//						i_neighbour_node->Set(VISITED);
		//						is_finished = false;
		//						i_neighbour_node->GetValue(FATHER_NODES).push_back(*(i_node.base()));
		//				}
		//			}

		//		}
		//}

		//KRATOS_WATCH(loop_counter);


	}

	void MeshSmoothingProcess::SelectCoarseMeshNodes(ModelPart::MeshType& rOriginalMesh)
	{
		// initializing and selecting all boundary nodes as seed
		for (ModelPart::NodeIterator i_node = rOriginalMesh.NodesBegin(); i_node != rOriginalMesh.NodesEnd(); i_node++)
		{
			i_node->GetValue(FATHER_NODES).clear();
			i_node->Set(COARSE_MESH_NODE, false);
			if (i_node->Is(BOUNDARY)) // We want to keep all the boundary nodes 
				i_node->Set(COARSE_MESH_NODE);
		}

		// Now selecting the seeds
		// A node without any seed in its neighbour is a seed
		for (ModelPart::NodeIterator i_node = rOriginalMesh.NodesBegin(); i_node != rOriginalMesh.NodesEnd(); i_node++)
		{
			if (i_node->IsNot(COARSE_MESH_NODE))
			{
				bool is_seed = true; // as a candidate
				GlobalPointersVector< Node<3> >& r_neighbours = i_node->GetValue(NEIGHBOUR_NODES);
				for (GlobalPointersVector<Node<3> >::iterator i_neighbour_node = r_neighbours.begin(); i_neighbour_node != r_neighbours.end(); i_neighbour_node++)
				{
					if ((i_neighbour_node->Is(COARSE_MESH_NODE)))
					{
						is_seed = false;
						break;
					}
				}
				if (is_seed)
				{
					i_node->Set(COARSE_MESH_NODE);
					i_node->GetValue(FATHER_NODES).push_back(*(i_node.base()));
				}
			}
		}
	}

	void MeshSmoothingProcess::CollapseNodes(ModelPart::MeshType& rOriginalMesh)
	{
		for (ModelPart::NodeIterator i_node = rOriginalMesh.NodesBegin(); i_node != rOriginalMesh.NodesEnd(); i_node++)
		{
			if (i_node->IsNot(COARSE_MESH_NODE))
			{
				CollapseToNearestCoarseNode(*i_node);
			}
		}
	}


	void  MeshSmoothingProcess::CollapseToNearestCoarseNode(NodeType& rThisNode)
	{
		double min_distance = std::numeric_limits<double>::max();
		
		GlobalPointersVector< Node<3> >& r_neighbours = rThisNode.GetValue(NEIGHBOUR_NODES);
		GlobalPointersVector<Node<3> >::iterator i_coarse_node = r_neighbours.end();

		for (GlobalPointersVector<Node<3> >::iterator i_neighbour_node = r_neighbours.begin(); i_neighbour_node != r_neighbours.end(); i_neighbour_node++)
		{
			if ((i_neighbour_node->Is(COARSE_MESH_NODE)))
			{
				Point d = *i_neighbour_node - rThisNode;
				double distance2 = inner_prod(d,d);
				if (distance2 < min_distance)
				{
					min_distance = distance2;
					i_coarse_node = i_neighbour_node;
				}
			}
		}

		if (i_coarse_node == r_neighbours.end())
		{
			KRATOS_ERROR << "The node " << rThisNode.Id() << "cannot be collapsed" << std::endl;
		}
		else
		{
			rThisNode.GetValue(FATHER_NODES).push_back(*(i_coarse_node.base()));
		}
		
	}

	void MeshSmoothingProcess::CreateCoarseMeshElements(ModelPart::MeshType& rOriginalMesh, ModelPart::MeshType& rCoarsMesh)
	{
		int id = 0;
		for (ModelPart::ElementIterator i_element = rOriginalMesh.ElementsBegin(); i_element != rOriginalMesh.ElementsEnd(); i_element++)
		{

			if (IsNotCollapsed(*i_element))
			{
				Element::GeometryType::ContainerType element_nodes;
				Element::GeometryType& r_geometry = i_element->GetGeometry();
				for (unsigned int i = 0; i < r_geometry.size(); i++)
				{
					if (r_geometry[i].GetValue(FATHER_NODES).empty())
					{
						id = r_geometry[i].Id();
						rCoarsMesh.AddNode(r_geometry.pGetPoint(i));
					}
					else
						id = r_geometry[i].GetValue(FATHER_NODES).front().Id();
					element_nodes.push_back(rOriginalMesh.pGetNode(id));
				}

				Element::Pointer p_new_element = i_element->Clone(i_element->Id(), element_nodes);

				if (p_new_element->GetGeometry().DomainSize() < 0.1 * i_element->GetGeometry().DomainSize())
					ChangeElmentCoarseNodes(*i_element); // Find another coarse nodes to collapse in order to get better result.

				rCoarsMesh.AddElement(p_new_element);
			} 
		}


	}

	void MeshSmoothingProcess::ChangeElmentCoarseNodes(Element& ThisElement)
	{
		Element::GeometryType::ContainerType element_nodes;
		Element::GeometryType& r_geometry = ThisElement.GetGeometry();
		//for (unsigned int i = 0; i < r_geometry.size(); i++)
		//{
		//	if (r_geometry[i].GetValue(FATHER_NODES).empty())
		//	{
		//		id = r_geometry[i].Id();
		//		rCoarsMesh.AddNode(r_geometry.pGetPoint(i));
		//	}
		//	else
		//		id = r_geometry[i].GetValue(FATHER_NODES).front().Id();
		//	element_nodes.push_back(rOriginalMesh.pGetNode(id));
		//}

	}

	bool MeshSmoothingProcess::IsNotCollapsed(Element& ThisElement)
	{
		Element::GeometryType& r_geometry = ThisElement.GetGeometry();
		//for (unsigned int i = 0; i < r_geometry.size(); i++)
		//{
		//	if (r_geometry[i].Is(SELECTED))
		//		return false; // All element around a seed will be collapsed
		//	if (r_geometry[i].GetValue(FATHER_NODES).empty())
		//		KRATOS_ERROR << "A non collapsed node has been found: Node#" << r_geometry[i].Id() << std::endl;
		//}
		// I can do this in the previous loop but with checking the father of r_geometry[j].
		int id1, id2;
		for (unsigned int i = 0; i < r_geometry.size(); i++)
		{
			if (r_geometry[i].GetValue(FATHER_NODES).empty())
				id1 = r_geometry[i].Id();
			else
				id1 = r_geometry[i].GetValue(FATHER_NODES).front().Id();

			for (unsigned int j = i + 1; j < r_geometry.size(); j++)
			{
				if (r_geometry[j].GetValue(FATHER_NODES).empty())
					id2 = r_geometry[j].Id();
				else
					id2 = r_geometry[j].GetValue(FATHER_NODES).front().Id();

				if (id1 == id2)
					return false;
			}
		}

		return true;
	}

	void MeshSmoothingProcess::MoveNode(NodeType& rNode, NeighboursVectorType const& rNeighbours, PointsVectorType const& rOptimumPoints, Vector const& rWeights)
	{
		if (rNode.IsNot(BOUNDARY))
		{
			std::size_t size = rOptimumPoints.size();
			if (size > 0)
			{
				Point optimal_position = ZeroVector(3);
				double weight_sum = 0.00;

				for (std::size_t i = 0; i < size; i++)
				{
					optimal_position += rOptimumPoints[i] * rWeights[i];
					weight_sum += rWeights[i];
				}

				rNode.Coordinates() = optimal_position / weight_sum;
			}
		}
	}

	void MeshSmoothingProcess::FindOptimumPositionsAndWeights(NodeType& rNode, NeighboursVectorType const& rNeighbours, PointsVectorType& rOptimumPoints, Vector& rWeights)
	{
		const std::size_t size = rNeighbours.size();
		rOptimumPoints.resize(size);
		rWeights.resize(size);
		if (mOptions.Is(LAPLACIAN_SMOOTHING))
			LaplacianSmoothingPositionsAndWeights(rNode, rNeighbours, rOptimumPoints, rWeights);
		else if (mOptions.Is(EDGE_LENGTH_SMOOTHING))
			EdgeLengthSmoothingPositionsAndWeights(rNode, rNeighbours, rOptimumPoints, rWeights);
		else if (mOptions.Is(MOVEMENT_SMOOTHING))
			MovementSmoothingPositionsAndWeights(rNode, rNeighbours, rOptimumPoints, rWeights);
		else
			KRATOS_ERROR << "Invalid type of smoothing in options." << std::endl
			<< "The valid ones are:" << std::endl
			<< "    LAPLACIAN_SMOOTHING" << std::endl
			<< "    EDGE_LENGTH_SMOOTHING" << std::endl
			<< "    MOVEMENT_SMOOTHING" << std::endl;

	}

	void MeshSmoothingProcess::LaplacianSmoothingPositionsAndWeights(NodeType& rNode, NeighboursVectorType const& rNeighbours, PointsVectorType& rOptimumPoints, Vector& rWeights)
	{
		const std::size_t size = rNeighbours.size();
		for (std::size_t i = 0; i < size; i++)
		{
			rOptimumPoints[i] = rNeighbours[i];
			rWeights[i] = 1.00;
		}
	}

	void MeshSmoothingProcess::EdgeLengthSmoothingPositionsAndWeights(NodeType& rNode, NeighboursVectorType const& rNeighbours, PointsVectorType& rOptimumPoints, Vector& rWeights)
	{
		const double epsilon = 1e-12;
		const std::size_t size = rNeighbours.size();
		double average_h = 0.00;
		for (std::size_t i = 0; i < size; i++)
		{
			rWeights[i] = 1.00;
			if (rNode.Is(SELECTED))
				if (rNeighbours[i].IsNot(SELECTED))
					rOptimumPoints[i] = rNode - rNeighbours[i];
				else
				{
					rWeights[i] = 2.00;
					rOptimumPoints[i] = rNeighbours[i] - rNode;
				}
			else
				rOptimumPoints[i] = rNode - rNeighbours[i];
			double l = norm_2(rOptimumPoints[i]);
			if (l < epsilon)
			{
				rOptimumPoints[i] = ZeroVector(3);
				rWeights[i] = 0.00;
			}
			else
			{
				rOptimumPoints[i] /= l;// *= rWeights[i];
			}
			average_h += l;
		}

		average_h /= size;
		for (int i = 0; i < size; i++)
		{
			rOptimumPoints[i] *= average_h;
			rOptimumPoints[i] += rNeighbours[i];
		}
	}

	void MeshSmoothingProcess::MovementSmoothingPositionsAndWeights(NodeType& rNode, NeighboursVectorType const& rNeighbours, PointsVectorType& rOptimumPoints, Vector& rWeights)
	{
	}


}  // namespace Kratos.


