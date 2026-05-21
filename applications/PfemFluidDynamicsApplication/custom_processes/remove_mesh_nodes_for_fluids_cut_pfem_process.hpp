//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:                   AFranci $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                October 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined(KRATOS_REMOVE_MESH_NODES_FOR_FLUIDS_CUT_PFEM_PROCESS_H_INCLUDED)
#define KRATOS_REMOVE_MESH_NODES_FOR_FLUIDS_CUT_PFEM_PROCESS_H_INCLUDED

// External includes

// System includes

// Project includes
#include "custom_processes/remove_mesh_nodes_for_fluids_process.hpp"

/// VARIABLES used:
// Data:     NORMAL, MASTER_NODES, NEIGHBOUR_NODES, NEIGBOUR_ELEMENTS
// Flags:    (checked) TO_ERASE, BOUNDARY, STRUCTURE, TO_SPLIT, CONTACT, NEW_ENTITY, BLOCKED
//           (set)     TO_ERASE(conditions,nodes)(set), NEW_ENTITY(conditions,nodes)(set), BLOCKED(nodes)->locally, VISITED(nodes)(set)
//           (modified)
//           (reset)   BLOCKED->locally
//(set):=(set in this process)

namespace Kratos
{

	///@name Kratos Classes
	///@{

	/// Remove Mesh Nodes Process for 2D and 3D cases
	/** The process labels the nodes to be erased (TO_ERASE)
	if they are too close (mRemoveOnDistance == true)
	if the error of the patch they belong is very small (REMOVE_NODES_ON_ERROR)
	In the interior of the domain or in the boundary (REMOVE_BOUNDARY_NODES) ...

	Additional treatment of the nonconvex boundaries is also going to erase nodes.

	At the end of the execution nodes are cleaned (only in the current mesh)
*/

	class RemoveMeshNodesForFluidsCutPfemProcess
		: public RemoveMeshNodesForFluidsProcess
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of Process
		KRATOS_CLASS_POINTER_DEFINITION(RemoveMeshNodesForFluidsCutPfemProcess);

		typedef ModelPart::ConditionType ConditionType;
		typedef ModelPart::PropertiesType PropertiesType;
		typedef ConditionType::GeometryType GeometryType;
		typedef Bucket<3, Node, std::vector<Node::Pointer>, Node::Pointer, std::vector<Node::Pointer>::iterator, std::vector<double>::iterator> BucketType;
		typedef Tree<KDTreePartition<BucketType>> KdtreeType; // Kdtree
		typedef ModelPart::MeshType::GeometryType::PointsArrayType PointsArrayType;
		typedef GlobalPointersVector<Node> NodeWeakPtrVectorType;
		typedef GlobalPointersVector<Element> ElementWeakPtrVectorType;
		typedef std::size_t SizeType;

		typedef RemoveMeshNodesForFluidsProcess BaseType;

		///@}
		///@name Life Cycle
		///@{

		/// Default constructor.
		RemoveMeshNodesForFluidsCutPfemProcess(ModelPart &rModelPart,
											   MesherUtilities::MeshingParameters &rRemeshingParameters,
											   int EchoLevel)
			: BaseType(rModelPart, rRemeshingParameters, EchoLevel)
		{
			KRATOS_INFO("RemoveMeshNodesForFluidsCutPfemProcess") << " activated " << std::endl;
		}

		/// Destructor.
		virtual ~RemoveMeshNodesForFluidsCutPfemProcess() {}

		///@}
		///@name Operators
		///@{

		/// This operator is provided to call the process as a function and simply calls the Execute method.
		void operator()()
		{
			Execute();
		}

		///@}
		///@name Operations
		///@{

		/// Execute method is used to execute the Process algorithms.
		void Execute() override
		{

			KRATOS_TRY
			if (mEchoLevel > 1)
			{
				std::cout << " [ REMOVE CLOSE NODES: " << std::endl;
			}

			SizeType inside_nodes_removed = 0;
			SizeType boundary_nodes_removed = 0;
			SizeType nodes_removed_inlet_zone = 0;

			// if the remove_node switch is activated, we check if the nodes got too close
			if (mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_NODES))
			{
				bool some_node_is_removed = false;

				if (mEchoLevel > 1)
					std::cout << " REMOVE_NODES is TRUE " << std::endl;

				if (mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_NODES_ON_DISTANCE))
				{
					if (mEchoLevel > 1)
						std::cout << " REMOVE_NODES_ON_DISTANCE is TRUE " << std::endl;

					some_node_is_removed = RemoveNodesOnDistance(inside_nodes_removed, boundary_nodes_removed, nodes_removed_inlet_zone);
				}

				if (some_node_is_removed || mrRemesh.UseBoundingBox == true)
					this->CleanRemovedNodes(mrModelPart);
			}

			// number of removed nodes:
			mrRemesh.Info->RemovedNodes += inside_nodes_removed + boundary_nodes_removed - nodes_removed_inlet_zone;
			int distance_remove = inside_nodes_removed + boundary_nodes_removed;

			if (mEchoLevel > 1)
			{
				std::cout << "   [ NODES      ( removed : " << mrRemesh.Info->RemovedNodes << " ) ]" << std::endl;
				std::cout << "   [ Distance(removed: " << distance_remove << "; inside: " << inside_nodes_removed << "; boundary: " << boundary_nodes_removed << ") ]" << std::endl;
				std::cout << "   REMOVE CLOSE NODES ]; " << std::endl;
			}

			KRATOS_CATCH(" ")
		}

		///@}
		///@name Access
		///@{

		///@}
		///@name Inquiry
		///@{

		///@}
		///@name Input and output
		///@{

		/// Turn back information as a string.
		std::string Info() const override
		{
			return "RemoveMeshNodesForFluidsCutPfemProcess";
		}

		/// Print information about this object.
		void PrintInfo(std::ostream &rOStream) const override
		{
			rOStream << "RemoveMeshNodesForFluidsCutPfemProcess";
		}

		/// Print object's data.
		void PrintData(std::ostream &rOStream) const override
		{
		}

		///@}
		///@name Friends
		///@{

		///@}

	private:
		///@name Static Member Variables
		///@{

		///@}
		///@name Static Member Variables
		///@{
		//**************************************************************************
		//**************************************************************************

		//**************************************************************************
		//**************************************************************************

		bool RemoveNodesOnDistance(SizeType &inside_nodes_removed,
								   SizeType &boundary_nodes_removed,
								   SizeType &nodes_removed_inlet_zone)
		{
			KRATOS_TRY

			const SizeType dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

			bool some_node_is_removed = false;

			// bucket size definition:
			SizeType bucket_size = 20;

			// create the list of the nodes to be check during the search
			std::vector<Node::Pointer> list_of_nodes;
			list_of_nodes.reserve(mrModelPart.NumberOfNodes());
			for (ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++)
			{
				(list_of_nodes).push_back(*(i_node.base()));
			}

			KdtreeType nodes_tree(list_of_nodes.begin(), list_of_nodes.end(), bucket_size);
			////////////////////////////////////////////////////////////
			// all of the nodes in this list will be preserved
			SizeType num_neighbours = 100;
			std::vector<Node::Pointer> neighbours(num_neighbours);
			std::vector<double> neighbour_distances(num_neighbours);

			// radius means the distance, if the distance between two nodes is closer to radius -> mark for removing
			double radius = 0;
			Node work_point(0, 0.0, 0.0, 0.0);
			SizeType n_points_in_radius;

			const ProcessInfo &rCurrentProcessInfo = mrModelPart.GetProcessInfo();
			const double currentTime = rCurrentProcessInfo[TIME];
			double meshSize = mrRemesh.Refine->CriticalRadius;

			bool refiningBox = false;
			for (SizeType index = 0; index < mrRemesh.UseRefiningBox.size(); index++)
			{
				if (mrRemesh.UseRefiningBox[index] == true && currentTime > mrRemesh.RefiningBoxInitialTime[index] && currentTime < mrRemesh.RefiningBoxFinalTime[index])
				{
					refiningBox = true;
				}
			}

			const SizeType principalModelPartId = rCurrentProcessInfo[MAIN_MATERIAL_PROPERTY];
			for (ModelPart::ElementsContainerType::const_iterator ie = mrModelPart.ElementsBegin();
				 ie != mrModelPart.ElementsEnd(); ie++)
			{
				SizeType rigidNodes = 0;
				// coordinates
				for (SizeType i = 0; i < ie->GetGeometry().size(); i++)
				{
					if ((ie->GetGeometry()[i].Is(RIGID) && ie->GetGeometry()[i].IsNot(PFEMFlags::LAGRANGIAN_INLET)) ||
						ie->GetGeometry()[i].Is(SOLID) ||
						ie->GetGeometry()[i].Is(PFEMFlags::EULERIAN_INLET) ||
						ie->GetGeometry()[i].GetSolutionStepValue(DISTANCE) < 0)
					{
						rigidNodes++;
					}
				}
				if (dimension == 2)
				{
					if (rigidNodes > 0)
						this->EraseCriticalNodes2D(ie->GetGeometry(), inside_nodes_removed, nodes_removed_inlet_zone);
				}
				else if (dimension == 3)
				{
					if (rigidNodes > 1)
						EraseCriticalNodes3D(ie->GetGeometry(), inside_nodes_removed, nodes_removed_inlet_zone, rigidNodes);
				}
			}

			for (ModelPart::NodesContainerType::const_iterator in = mrModelPart.NodesBegin(); in != mrModelPart.NodesEnd(); in++)
			{

				const SizeType propertyIdNode = in->FastGetSolutionStepValue(PROPERTY_ID);
				meshSize = mrRemesh.Refine->CriticalRadius;
				double rigidNodeLocalMeshSize = 0;
				double rigidNodeMeshCounter = 0;
				NodeWeakPtrVectorType &rN = in->GetValue(NEIGHBOUR_NODES);

				for (SizeType i = 0; i < rN.size(); i++)
				{
					if (rN[i].Is(RIGID))
					{
						rigidNodeLocalMeshSize += rN[i].FastGetSolutionStepValue(NODAL_H_WALL);
						rigidNodeMeshCounter += 1.0;
					}
				}

				if (refiningBox == true)
				{
					array_1d<double, 3> NodeCoordinates = in->Coordinates();
					if (dimension == 2)
					{
						bool insideTransitionZone = false;
						mMesherUtilities.DefineMeshSizeInTransitionZones2D(mrRemesh, currentTime, NodeCoordinates, meshSize, insideTransitionZone);
					}
					else if (dimension == 3)
					{
						bool insideTransitionZone = false;
						mMesherUtilities.DefineMeshSizeInTransitionZones3D(mrRemesh, currentTime, NodeCoordinates, meshSize, insideTransitionZone);
					}
				}

				const double distance_tolerance = 0.05 * meshSize;

				this->SetMeshSizeNearBoundaries(meshSize, rigidNodeLocalMeshSize, rigidNodeMeshCounter);

				const double size_for_distance_boundary = 0.6 * meshSize;
				const double size_for_wall_tip_contact_side = 0.15 * mrRemesh.Refine->CriticalSide;

				if (in->Is(TO_ERASE))
				{
					some_node_is_removed = true;
				}
				bool on_contact_tip = false;

				if (in->Is(TO_SPLIT) || in->Is(CONTACT))
					on_contact_tip = true;

				if (in->IsNot(NEW_ENTITY) && in->IsNot(INLET) && in->IsNot(ISOLATED) && in->IsNot(RIGID) && in->IsNot(SOLID) && (dimension == 2 || (dimension == 3 && in->GetSolutionStepValue(DISTANCE) > 0)))
				{
					SizeType neighErasedNodes = 0;
					SizeType freeSurfaceNeighNodes = 0;
					bool inletElement = false;
					bool interfaceElement = false;

					radius = 0.6 * meshSize;
					work_point[0] = in->X();
					work_point[1] = in->Y();
					work_point[2] = in->Z();

					if (in->Is(FREE_SURFACE))
					{
						// it must be more difficult to erase a free_surface node, otherwise, lot of volume is lost
						// this value has a strong effect on volume variation due to remeshing
						radius = 0.475 * meshSize; // compared with element radius
						NodeWeakPtrVectorType &neighb_nodes = in->GetValue(NEIGHBOUR_NODES);
						SizeType countRigid = 0;

						for (NodeWeakPtrVectorType::iterator nn = neighb_nodes.begin(); nn != neighb_nodes.end(); nn++)
						{
							if ((nn)->Is(RIGID) || (nn)->Is(SOLID))
							{
								countRigid++;
							}
							if ((nn)->Is(TO_ERASE))
							{
								neighErasedNodes++;
							}
						}
						if (countRigid == neighb_nodes.size())
						{
							radius = 0.15 * meshSize;
						}
					}
					else
					{
						NodeWeakPtrVectorType &neighb_nodes = in->GetValue(NEIGHBOUR_NODES);
						for (NodeWeakPtrVectorType::iterator nn = neighb_nodes.begin(); nn != neighb_nodes.end(); nn++)
						{
							const SizeType propertyIdSecondNode = (nn)->FastGetSolutionStepValue(PROPERTY_ID);
							if ((nn)->Is(FREE_SURFACE))
							{
								freeSurfaceNeighNodes++;
							}
							if ((nn)->Is(TO_ERASE))
							{
								neighErasedNodes++;
							}
							if ((nn)->Is(PFEMFlags::EULERIAN_INLET))
							{
								inletElement = true;
							}
							if (propertyIdNode != propertyIdSecondNode && (nn)->IsNot(RIGID))
							{
								interfaceElement = true;
							}
						}
					}

					if (freeSurfaceNeighNodes > 1)
					{
						radius = 0.5 * meshSize;
					}
					else if (interfaceElement == true)
					{
						if (dimension == 2)
							radius = 0.54 * meshSize; // 10% less than normal nodes
						if (dimension == 3)
							radius = 0.48 * meshSize; // 20% less than normal nodes
					}
					if (in->GetSolutionStepValue(DISTANCE) < distance_tolerance)
					{
						if (dimension == 2)
							radius = 0.2 * meshSize; // 1/3 of normal nodes
						if (dimension == 3)
							radius = 0.15 * meshSize; // 1/3 of normal nodes
					}
					else if (in->GetSolutionStepValue(DISTANCE) < distance_tolerance * 2.0)
					{
						if (dimension == 2)
							radius = 0.4 * meshSize; // 1/3 of normal nodes
						if (dimension == 3)
							radius = 0.3 * meshSize; // 1/3 of normal nodes
					}
					n_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, neighbours.begin(), neighbour_distances.begin(), num_neighbours);

					if (n_points_in_radius > 1 && neighErasedNodes == 0 && in->IsNot(INLET))
					{

						if (in->IsNot(RIGID) && in->IsNot(SOLID) && in->IsNot(ISOLATED))
						{
							if (mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_NODES_ON_DISTANCE))
							{

								if (in->IsNot(FREE_SURFACE) && in->IsNot(RIGID) && freeSurfaceNeighNodes == dimension)
								{
									this->SetVariablesToFreeSurfaceElements(in, dimension);
								}
								else
								{
									in->Set(TO_ERASE);
									some_node_is_removed = true;
									inside_nodes_removed++;
									if (inletElement)
									{
										nodes_removed_inlet_zone++;
									}
									if (propertyIdNode != principalModelPartId) // this is to conserve the number of nodes of the smaller domain in case of a two-fluid analysis
									{
										mrRemesh.Info->BalancePrincipalSecondaryPartsNodes += -1;
									}
								}
							}
						}
						else
						{
							SizeType k = 0;
							SizeType counter = 0;
							for (std::vector<Node::Pointer>::iterator nn = neighbours.begin(); nn != neighbours.begin() + n_points_in_radius; nn++)
							{
								bool nn_on_contact_tip = false;

								if ((*nn)->Is(TO_SPLIT) || (*nn)->Is(CONTACT))
									nn_on_contact_tip = true;

								if ((*nn)->Is(BOUNDARY) && !nn_on_contact_tip && neighbour_distances[k] < size_for_distance_boundary && neighbour_distances[k] > 0.0)
								{
									if ((*nn)->IsNot(TO_ERASE))
									{
										counter += 1;
									}
								}

								if ((*nn)->Is(BOUNDARY) && nn_on_contact_tip && neighbour_distances[k] < size_for_wall_tip_contact_side)
								{
									if ((*nn)->IsNot(TO_ERASE))
									{
										counter += 1;
									}
								}
								k++;
							}

							if (counter > 1 && in->IsNot(RIGID) && in->IsNot(SOLID) && in->IsNot(NEW_ENTITY) && !on_contact_tip && (dimension == 2 || (dimension == 3 && in->GetSolutionStepValue(DISTANCE) > distance_tolerance)))
							{ // Can be inserted in the boundary refine
								in->Set(TO_ERASE);
								if (mEchoLevel > 1)
									std::cout << "     Removed Boundary Node [" << in->Id() << "] on Distance " << std::endl;
								some_node_is_removed = true;
								boundary_nodes_removed++;

								if (propertyIdNode != principalModelPartId) // this is to conserve the number of nodes of the smaller domain in case of a two-fluid analysis
								{
									mrRemesh.Info->BalancePrincipalSecondaryPartsNodes += -1;
								}
							}
						}
					}
				}

				if (in->Is(ISOLATED) && in->GetSolutionStepValue(DISTANCE) < distance_tolerance && in->IsNot(TO_ERASE))
				{
					in->Set(TO_ERASE);
					some_node_is_removed = true;
					inside_nodes_removed++;
					if (propertyIdNode != principalModelPartId) // this is to conserve the number of nodes of the smaller domain in case of a two-fluid analysis
					{
						mrRemesh.Info->BalancePrincipalSecondaryPartsNodes += -1;
					}
				}
			}

			if (boundary_nodes_removed > 0 || inside_nodes_removed > 0)
			{
				some_node_is_removed = true;
			}
			// Build boundary after removing boundary nodes due distance criterion
			if (mEchoLevel > 1)
			{
				std::cout << "boundary_nodes_removed " << boundary_nodes_removed << std::endl;
				std::cout << "inside_nodes_removed " << inside_nodes_removed << std::endl;
			}
			return some_node_is_removed;

			KRATOS_CATCH(" ")
		}

		void EraseCriticalNodes3D(Element::GeometryType &eElement, SizeType &inside_nodes_removed, SizeType &nodes_removed_inlet_zone, SizeType rigidNodes)
		{

			KRATOS_TRY
			double safetyCoefficient3D = 0.6;

			const ProcessInfo &rCurrentProcessInfo = mrModelPart.GetProcessInfo();
			SizeType principalModelPartId = rCurrentProcessInfo[MAIN_MATERIAL_PROPERTY];

			SizeType freeSurfaceNodes = 0;
			SizeType numNodes = eElement.size();
			double elementVolume = eElement.Volume();
			double criticalVolume = 0.1 * mrRemesh.Refine->MeanVolume;

			std::vector<array_1d<double, 3>> rigidNodesCoordinates;
			std::vector<array_1d<double, 3>> rigidNodesNormals;
			array_1d<double, 3> notRigidNodesCoordinates(3, 0.0);
			SizeType notRigidNodeId = 0;
			rigidNodesCoordinates.resize(3);
			rigidNodesNormals.resize(3);

			this->CriticalVolumeForRefinedBoxes(eElement, criticalVolume, safetyCoefficient3D);

			bool inletElement = false;
			if (eElement[0].Is(PFEMFlags::EULERIAN_INLET) || eElement[1].Is(PFEMFlags::EULERIAN_INLET) || eElement[2].Is(PFEMFlags::EULERIAN_INLET) || eElement[3].Is(PFEMFlags::EULERIAN_INLET))
			{
				inletElement = true;
			}

			SizeType rigidNode = 0;
			array_1d<double, 3> WallBaricenter = ZeroVector(3);

			for (SizeType i = 0; i < numNodes; i++)
			{
				if (eElement[i].Is(FREE_SURFACE))
				{
					freeSurfaceNodes++;
				}
				if (rigidNodes == 3)
				{
					if (eElement[i].Is(RIGID))
					{
						WallBaricenter += eElement[i].Coordinates() / 3.0;
						rigidNodesCoordinates[rigidNode] = eElement[i].Coordinates();
						rigidNodesNormals[rigidNode] = eElement[i].FastGetSolutionStepValue(NORMAL);
						rigidNode++;
					}
					else
					{
						notRigidNodesCoordinates = eElement[i].Coordinates();
						notRigidNodeId = i;
					}
				}
				if (elementVolume < criticalVolume)
				{

					if (eElement[i].IsNot(RIGID) && eElement[i].IsNot(SOLID) && eElement[i].IsNot(TO_ERASE) && eElement[i].GetSolutionStepValue(DISTANCE) > 0)
					{
						eElement[i].Set(TO_ERASE);
						if (mEchoLevel > 1)
							std::cout << "erase this layer node because it may be potentially dangerous and pass through the solid contour" << std::endl;
						inside_nodes_removed++;

						if (inletElement)
						{
							nodes_removed_inlet_zone++;
						}

						const SizeType propertyIdNode = eElement[i].FastGetSolutionStepValue(PROPERTY_ID);
						if (propertyIdNode != principalModelPartId) // this is to conserve the number of nodes of the smaller domain in case of a two-fluid analysis
						{
							mrRemesh.Info->BalancePrincipalSecondaryPartsNodes += -1;
						}
					}
				}
			}
			double wallNodeDistance = -1;
			if (rigidNode == 3 && eElement[notRigidNodeId].IsNot(TO_ERASE))
			{
				this->CheckSliverElements(eElement,
										  inside_nodes_removed,
										  nodes_removed_inlet_zone,
										  inletElement,
										  notRigidNodeId,
										  principalModelPartId,
										  wallNodeDistance,
										  WallBaricenter,
										  rigidNodesCoordinates,
										  rigidNodesNormals,
										  notRigidNodesCoordinates);
			}

			array_1d<double, 6> Edges(6, 0.0);
			array_1d<SizeType, 6> FirstEdgeNode(6, 0);
			array_1d<SizeType, 6> SecondEdgeNode(6, 0);
			int removableNode = -1;
			double wallLength = 0;
			double minimumLength = 0;
			this->DetectRemovableNodes(eElement, Edges, FirstEdgeNode, SecondEdgeNode, removableNode, wallLength, minimumLength, safetyCoefficient3D);

			if ((minimumLength < (0.35 * wallLength) || (wallNodeDistance < (0.25 * wallLength) && wallNodeDistance > 0)) && removableNode > -1)
			{
				if (eElement[removableNode].IsNot(RIGID) && eElement[removableNode].IsNot(TO_ERASE) && eElement[removableNode].IsNot(SOLID) && eElement[removableNode].IsNot(ISOLATED) && eElement[removableNode].IsNot(FREE_SURFACE))
				{
					eElement[removableNode].Set(TO_ERASE);
					inside_nodes_removed++;

					if (inletElement)
					{
						nodes_removed_inlet_zone++;
					}
					const SizeType propertyIdNode = eElement[removableNode].FastGetSolutionStepValue(PROPERTY_ID);
					if (propertyIdNode != principalModelPartId) // this is to conserve the number of nodes of the smaller domain in case of a two-fluid analysis
					{
						mrRemesh.Info->BalancePrincipalSecondaryPartsNodes += -1;
					}
				}
			}

			// ////////  to compare the non-wall length to wall edge length /////////
			this->CompareEdgeLengthsNearWalls(eElement, Edges, FirstEdgeNode, SecondEdgeNode, freeSurfaceNodes, safetyCoefficient3D, wallLength, inletElement, inside_nodes_removed, nodes_removed_inlet_zone);

			KRATOS_CATCH("")
		}

		/// Assignment operator.
		RemoveMeshNodesForFluidsCutPfemProcess &operator=(RemoveMeshNodesForFluidsCutPfemProcess const &rOther);

		/// this function is a private function

		/// Copy constructor.
		// Process(Process const& rOther);

		///@}

	}; // Class Process

	///@}

	///@name Type Definitions
	///@{

	///@}
	///@name Input and output
	///@{

	/// input stream function
	inline std::istream &operator>>(std::istream &rIStream,
									RemoveMeshNodesForFluidsCutPfemProcess &rThis);

	/// output stream function
	inline std::ostream &operator<<(std::ostream &rOStream,
									const RemoveMeshNodesForFluidsCutPfemProcess &rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@}

} // namespace Kratos.

#endif // KRATOS_REMOVE_MESH_NODES_FOR_FLUIDS_CUT_PFEM_PROCESS_H_INCLUDED  defined
