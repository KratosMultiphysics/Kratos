//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:                   AFranci $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                October 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined(KRATOS_REMOVE_MESH_NODES_FOR_FLUIDS_PROCESS_H_INCLUDED)
#define KRATOS_REMOVE_MESH_NODES_FOR_FLUIDS_PROCESS_H_INCLUDED

// External includes

// System includes

// Project includes
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "includes/model_part.h"
#include "custom_utilities/mesh_error_calculation_utilities.hpp"
#include "custom_utilities/mesher_utilities.hpp"

#include "pfem_fluid_dynamics_application_variables.h"
#include "custom_processes/mesher_process.hpp"

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

	class RemoveMeshNodesForFluidsProcess
		: public MesherProcess
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of Process
		KRATOS_CLASS_POINTER_DEFINITION(RemoveMeshNodesForFluidsProcess);

		typedef ModelPart::ConditionType ConditionType;
		typedef ModelPart::PropertiesType PropertiesType;
		typedef ConditionType::GeometryType GeometryType;
		typedef Bucket<3, Node, std::vector<Node::Pointer>, Node::Pointer, std::vector<Node::Pointer>::iterator, std::vector<double>::iterator> BucketType;
		typedef Tree<KDTreePartition<BucketType>> KdtreeType; // Kdtree
		typedef ModelPart::MeshType::GeometryType::PointsArrayType PointsArrayType;
		typedef GlobalPointersVector<Node> NodeWeakPtrVectorType;
		typedef GlobalPointersVector<Element> ElementWeakPtrVectorType;
		typedef std::size_t SizeType;
		///@}
		///@name Life Cycle
		///@{

		/// Default constructor.
		RemoveMeshNodesForFluidsProcess(ModelPart &rModelPart,
										MesherUtilities::MeshingParameters &rRemeshingParameters,
										int EchoLevel)
			: mrModelPart(rModelPart),
			  mrRemesh(rRemeshingParameters)
		{
			KRATOS_INFO("RemoveMeshNodesForFluidsProcess") << " activated " << std::endl;

			mEchoLevel = EchoLevel;
		}

		/// Destructor.
		virtual ~RemoveMeshNodesForFluidsProcess() {}

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
			return "RemoveMeshNodesForFluidsProcess";
		}

		/// Print information about this object.
		void PrintInfo(std::ostream &rOStream) const override
		{
			rOStream << "RemoveMeshNodesForFluidsProcess";
		}

		/// Print object's data.
		void PrintData(std::ostream &rOStream) const override
		{
		}

		///@}
		///@name Friends
		///@{

		///@}

	protected:
		///@name Static Member Variables
		///@{

		///@}
		///@name Static Member Variables
		///@{
		ModelPart &mrModelPart;

		MesherUtilities::MeshingParameters &mrRemesh;

		MesherUtilities mMesherUtilities;

		int mEchoLevel;

		///@}
		///@name Un accessible methods
		///@{

		//**************************************************************************
		//**************************************************************************

		void CleanRemovedNodes(ModelPart &rModelPart)
		{
			KRATOS_TRY

			// MESH 0 total domain mesh
			ModelPart::NodesContainerType temporal_nodes;
			temporal_nodes.reserve(rModelPart.Nodes().size());

			temporal_nodes.swap(rModelPart.Nodes());

			for (ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin(); i_node != temporal_nodes.end(); i_node++)
			{
				if (i_node->IsNot(TO_ERASE))
				{

					/////////////////////////////////////////// here for BOUNDING BOX ///////////////////////////////////////////
					bool boundingBox = mrRemesh.UseBoundingBox;
					if (boundingBox == true && i_node->IsNot(RIGID))
					{
						const ProcessInfo &rCurrentProcessInfo = mrModelPart.GetProcessInfo();
						double currentTime = rCurrentProcessInfo[TIME];
						double initialTime = mrRemesh.BoundingBoxInitialTime;
						double finalTime = mrRemesh.BoundingBoxFinalTime;
						if (currentTime > initialTime && currentTime < finalTime)
						{
							array_1d<double, 3> BoundingBoxLowerPoint = mrRemesh.BoundingBoxLowerPoint;
							array_1d<double, 3> BoundingBoxUpperPoint = mrRemesh.BoundingBoxUpperPoint;
							if (i_node->X() < BoundingBoxLowerPoint[0] || i_node->Y() < BoundingBoxLowerPoint[1] || i_node->Z() < BoundingBoxLowerPoint[2] ||
								i_node->X() > BoundingBoxUpperPoint[0] || i_node->Y() > BoundingBoxUpperPoint[1] || i_node->Z() > BoundingBoxUpperPoint[2])
							{
								i_node->Set(TO_ERASE);
							}
							else
							{
								(rModelPart.Nodes()).push_back(*(i_node.base()));
							}
						}
						else
						{
							(rModelPart.Nodes()).push_back(*(i_node.base()));
						}
					}
					else
					{
						(rModelPart.Nodes()).push_back(*(i_node.base()));
					}
				}
			}

			rModelPart.Nodes().Sort();

			KRATOS_CATCH("")
		}

		void EraseCriticalNodes2D(Element::GeometryType &eElement, SizeType &inside_nodes_removed, SizeType &nodes_removed_inlet_zone)
		{

			KRATOS_TRY

			double safetyCoefficient2D = 0.5;
			double elementVolume = eElement.Area();

			const ProcessInfo &rCurrentProcessInfo = mrModelPart.GetProcessInfo();
			const SizeType principalModelPartId = rCurrentProcessInfo[MAIN_MATERIAL_PROPERTY];

			array_1d<double, 3> Edges(3, 0.0);
			array_1d<SizeType, 3> FirstEdgeNode(3, 0);
			array_1d<SizeType, 3> SecondEdgeNode(3, 0);
			double wallLength = 0;
			array_1d<double, 3> CoorDifference = eElement[1].Coordinates() - eElement[0].Coordinates();
			double SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1];
			Edges[0] = std::sqrt(SquaredLength);
			FirstEdgeNode[0] = 0;
			SecondEdgeNode[0] = 1;
			if (eElement[0].Is(RIGID) && eElement[1].Is(RIGID))
			{
				wallLength = Edges[0];
			}
			bool inletElement = false;
			if (eElement[0].Is(PFEMFlags::EULERIAN_INLET) || eElement[1].Is(PFEMFlags::EULERIAN_INLET) || eElement[2].Is(PFEMFlags::EULERIAN_INLET))
			{
				inletElement = true;
			}
			SizeType counter = 0;
			for (SizeType i = 2; i < eElement.size(); i++)
			{
				for (SizeType j = 0; j < i; j++)
				{
					noalias(CoorDifference) = eElement[i].Coordinates() - eElement[j].Coordinates();
					SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1];
					counter += 1;
					Edges[counter] = std::sqrt(SquaredLength);
					FirstEdgeNode[counter] = j;
					SecondEdgeNode[counter] = i;
					if (eElement[i].Is(RIGID) && eElement[j].Is(RIGID) && Edges[counter] > wallLength)
					{
						wallLength = Edges[counter];
					}
				}
			}

			////////  to compare the triangle height to wall edge length /////////
			for (SizeType i = 0; i < eElement.size(); i++)
			{
				if (eElement[i].IsNot(RIGID) && eElement[i].IsNot(TO_ERASE) && eElement[i].IsNot(SOLID) && eElement[i].IsNot(ISOLATED))
				{
					double height = elementVolume * 2.0 / wallLength;

					//////it is evident when a freesurface particle in touch with wall is erased --> reduce the safety coeff
					if (eElement[i].Is(FREE_SURFACE))
					{
						NodeWeakPtrVectorType &neighb_nodes = eElement[i].GetValue(NEIGHBOUR_NODES);
						SizeType countRigid = 0;
						SizeType countFreeSurface = 0;
						for (NodeWeakPtrVectorType::iterator nn = neighb_nodes.begin(); nn != neighb_nodes.end(); nn++)
						{
							if ((nn)->Is(RIGID) || (nn)->Is(SOLID))
							{
								countRigid++;
							}
							if ((nn)->Is(FREE_SURFACE) && (nn)->IsNot(RIGID) && (nn)->IsNot(SOLID))
							{
								countFreeSurface++;
							}
						}
						if ((countRigid + countFreeSurface) == neighb_nodes.size() && countRigid > 0)
						{
							safetyCoefficient2D = 0.25;
						}
					}

					////// if the node is very close to the wall is erased in any case
					if (height < (0.5 * safetyCoefficient2D * wallLength))
					{
						eElement[i].Set(TO_ERASE);
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

			bool longDamBreak = false; // to attivate in case of long dam breaks to avoid separeted elements in the water front
			if (longDamBreak)
				ControlForLongDamBreakProblems2D(eElement, inletElement, inside_nodes_removed, nodes_removed_inlet_zone);

			KRATOS_CATCH("")
		}

		void ControlForLongDamBreakProblems2D(Element::GeometryType &eElement,
											  bool inletElement,
											  SizeType &inside_nodes_removed,
											  SizeType &nodes_removed_inlet_zone)
		{
			KRATOS_TRY

			const SizeType principalModelPartId = mrModelPart.GetProcessInfo()[MAIN_MATERIAL_PROPERTY];

			for (SizeType i = 0; i < eElement.size(); i++)
			{
				if (eElement[i].Is(FREE_SURFACE) && eElement[i].IsNot(RIGID))
				{

					GlobalPointersVector<Element> &neighb_elems = eElement[i].GetValue(NEIGHBOUR_ELEMENTS);
					GlobalPointersVector<Node> &neighb_nodes = eElement[i].GetValue(NEIGHBOUR_NODES);

					if (neighb_elems.size() < 2)
					{
						eElement[i].Set(TO_ERASE);
						std::cout << "erased an isolated element node" << std::endl;
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
					else
					{
						if (neighb_nodes.size() < 4)
						{
							for (SizeType j = 0; j < neighb_nodes.size(); j++)
							{
								if (neighb_nodes[j].IsNot(FREE_SURFACE) && neighb_nodes[j].IsNot(RIGID))
								{
									break;
								}
								if (j == (neighb_nodes.size() - 1))
								{
									eElement[i].Set(TO_ERASE);
									std::cout << "_________________________          erased an isolated element node" << std::endl;
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
					}
				}
			}

			KRATOS_CATCH("")
		}

		void ControlForLongDamBreakProblems3D(Element::GeometryType &eElement,
											  bool inletElement,
											  SizeType &inside_nodes_removed,
											  SizeType &nodes_removed_inlet_zone)
		{
			KRATOS_TRY

			const SizeType principalModelPartId = mrModelPart.GetProcessInfo()[MAIN_MATERIAL_PROPERTY];
			for (SizeType i = 0; i < eElement.size(); i++)
			{
				if (eElement[i].Is(FREE_SURFACE) && eElement[i].IsNot(RIGID))
				{

					GlobalPointersVector<Element> &neighb_elems = eElement[i].GetValue(NEIGHBOUR_ELEMENTS);
					GlobalPointersVector<Node> &neighb_nodes = eElement[i].GetValue(NEIGHBOUR_NODES);

					if (neighb_elems.size() < 2)
					{
						eElement[i].Set(TO_ERASE);
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
					else
					{
						if (neighb_nodes.size() < 10)
						{
							SizeType freeSurfaceNodesNeigh = 0;
							for (SizeType j = 0; j < neighb_nodes.size(); j++)
							{
								if (neighb_nodes[j].Is(FREE_SURFACE) && neighb_nodes[j].IsNot(RIGID))
								{
									freeSurfaceNodesNeigh++;
								}
								if (neighb_nodes[j].IsNot(FREE_SURFACE) && neighb_nodes[j].IsNot(RIGID) && neighb_nodes[j].IsNot(TO_ERASE))
								{
									break;
								}
								if (j == (neighb_nodes.size() - 1) && freeSurfaceNodesNeigh < 2)
								{
									eElement[i].Set(TO_ERASE);
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
						else
						{
							for (SizeType j = 0; j < neighb_nodes.size(); j++)
							{
								if (neighb_nodes[j].IsNot(RIGID))
								{
									break;
								}
								if (j == (neighb_nodes.size() - 1))
								{
									eElement[i].Set(TO_ERASE);
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
					}
				}
			}

			KRATOS_CATCH("")
		}

		void CompareEdgeLengthsNearWalls(Element::GeometryType &eElement,
										 array_1d<double, 6> &Edges,
										 array_1d<SizeType, 6> &FirstEdgeNode,
										 array_1d<SizeType, 6> &SecondEdgeNode,
										 SizeType freeSurfaceNodes,
										 double &safetyCoefficient3D,
										 double wallLength,
										 bool inletElement,
										 SizeType &inside_nodes_removed,
										 SizeType &nodes_removed_inlet_zone)
		{
			KRATOS_TRY
			const SizeType principalModelPartId = mrModelPart.GetProcessInfo()[MAIN_MATERIAL_PROPERTY];

			// ////////  to compare the non-wall length to wall edge length /////////
			for (SizeType i = 0; i < Edges.size(); i++)
			{

				if (mrRemesh.ExecutionOptions.Is(MesherUtilities::REFINE_WALL_CORNER))
				{
					const ProcessInfo &rCurrentProcessInfo = mrModelPart.GetProcessInfo();
					double currentTime = rCurrentProcessInfo[TIME];
					double deltaTime = rCurrentProcessInfo[DELTA_TIME];
					if (freeSurfaceNodes == 0 && currentTime > 2.0 * deltaTime)
					{
						safetyCoefficient3D *= 0.7;
					}
					if (currentTime < 2.0 * deltaTime)
					{
						safetyCoefficient3D *= 1.05;
					}
				}
				if (((eElement[FirstEdgeNode[i]].Is(RIGID) && eElement[SecondEdgeNode[i]].IsNot(RIGID)) ||
					 (eElement[SecondEdgeNode[i]].Is(RIGID) && eElement[FirstEdgeNode[i]].IsNot(RIGID))) &&
					eElement[FirstEdgeNode[i]].IsNot(TO_ERASE) &&
					eElement[SecondEdgeNode[i]].IsNot(TO_ERASE) &&
					Edges[i] < safetyCoefficient3D * wallLength)
				{
					if (eElement[FirstEdgeNode[i]].IsNot(RIGID) && eElement[FirstEdgeNode[i]].IsNot(SOLID) && eElement[FirstEdgeNode[i]].IsNot(TO_ERASE) && eElement[FirstEdgeNode[i]].IsNot(ISOLATED))
					{

						eElement[FirstEdgeNode[i]].Set(TO_ERASE);
						inside_nodes_removed++;

						if (inletElement)
						{
							nodes_removed_inlet_zone++;
						}

						const SizeType propertyIdNode = eElement[FirstEdgeNode[i]].FastGetSolutionStepValue(PROPERTY_ID);
						if (propertyIdNode != principalModelPartId) // this is to conserve the number of nodes of the smaller domain in case of a two-fluid analysis
						{
							mrRemesh.Info->BalancePrincipalSecondaryPartsNodes += -1;
						}
					}
					else if (eElement[SecondEdgeNode[i]].IsNot(RIGID) && eElement[SecondEdgeNode[i]].IsNot(SOLID) && eElement[SecondEdgeNode[i]].IsNot(TO_ERASE) && eElement[SecondEdgeNode[i]].IsNot(ISOLATED))
					{
						eElement[SecondEdgeNode[i]].Set(TO_ERASE);
						inside_nodes_removed++;

						if (inletElement)
						{
							nodes_removed_inlet_zone++;
						}

						const SizeType propertyIdNode = eElement[SecondEdgeNode[i]].FastGetSolutionStepValue(PROPERTY_ID);
						if (propertyIdNode != principalModelPartId) // this is to conserve the number of nodes of the smaller domain in case of a two-fluid analysis
						{
							mrRemesh.Info->BalancePrincipalSecondaryPartsNodes += -1;
						}
					}
				}
			}

			KRATOS_CATCH("")
		}

		bool CheckForMovingLayerNodes(Node &CheckedNode, const double wallLength)
		{
			KRATOS_TRY
			const SizeType dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

			NodeWeakPtrVectorType &neighb_nodes = CheckedNode.GetValue(NEIGHBOUR_NODES);
			bool eraseNode = true;
			double maxSquaredDistance = 0;
			NodeWeakPtrVectorType::iterator j = neighb_nodes.begin();
			for (NodeWeakPtrVectorType::iterator nn = neighb_nodes.begin(); nn != neighb_nodes.end(); nn++)
			{
				if ((nn)->IsNot(RIGID) && (nn)->IsNot(SOLID))
				{
					// std::cout<<"neigh coordinates: "<<(nn)->X()<<" "<<(nn)->Y()<<std::endl;
					array_1d<double, 3> CoorNeighDifference = CheckedNode.Coordinates() - (nn)->Coordinates();
					double squaredDistance = CoorNeighDifference[0] * CoorNeighDifference[0] + CoorNeighDifference[1] * CoorNeighDifference[1];
					if (dimension == 3)
					{
						squaredDistance += CoorNeighDifference[2] * CoorNeighDifference[2];
					}
					if (squaredDistance > maxSquaredDistance)
					{
						// std::cout<<"(distances:  "<<squaredDistance<<" vs "<<maxSquaredDistance<<")"<<std::endl;
						maxSquaredDistance = squaredDistance;
						j = nn;
					}
				}
			}
			// I have looked for the biggest edge for moving there the layer node
			double maxNeighDistance = std::sqrt(maxSquaredDistance);
			if (maxNeighDistance > wallLength && wallLength > 0)
			{
				for (NodeWeakPtrVectorType::iterator nn = neighb_nodes.begin(); nn != neighb_nodes.end(); nn++)
				{
					if (nn == j)
					{

						SizeType idMaster = CheckedNode.GetId();
						SizeType idSlave = (j)->GetId();
						InterpolateFromTwoNodes(idMaster, idMaster, idSlave);
						std::vector<double> NewCoordinates(3);
						NewCoordinates[0] = (CheckedNode.X() + (j)->X()) * 0.5;
						NewCoordinates[1] = (CheckedNode.Y() + (j)->Y()) * 0.5;
						CheckedNode.X() = NewCoordinates[0];
						CheckedNode.Y() = NewCoordinates[1];
						CheckedNode.X0() = NewCoordinates[0];
						CheckedNode.Y0() = NewCoordinates[1];
						CheckedNode.FastGetSolutionStepValue(DISPLACEMENT_X, 0) = 0;
						CheckedNode.FastGetSolutionStepValue(DISPLACEMENT_Y, 0) = 0;
						CheckedNode.FastGetSolutionStepValue(DISPLACEMENT_X, 1) = 0;
						CheckedNode.FastGetSolutionStepValue(DISPLACEMENT_Y, 1) = 0;
						if (dimension == 3)
						{
							NewCoordinates[2] = (CheckedNode.Z() + (j)->Z()) * 0.5;
							CheckedNode.Z() = NewCoordinates[2];
							CheckedNode.Z0() = NewCoordinates[2];
							CheckedNode.FastGetSolutionStepValue(DISPLACEMENT_Z, 0) = 0;
							CheckedNode.FastGetSolutionStepValue(DISPLACEMENT_Z, 1) = 0;
						}
						// std::cout<<"new coordinates: "<<CheckedNode.X()<<" "<<CheckedNode.Y()<<std::endl;
						// std::cout<<"(distances:  "<<maxNeighDistance<<" vs "<<wallLength<<")"<<std::endl;
					}
					eraseNode = false;
				}
			}
			return eraseNode;

			KRATOS_CATCH("")
		}

		void InterpolateFromTwoNodes(SizeType idMaster, SizeType idSlave1, SizeType idSlave2)
		{

			KRATOS_TRY

			Node::Pointer MasterNode = mrModelPart.pGetNode(idMaster);
			Node::Pointer SlaveNode1 = mrModelPart.pGetNode(idSlave1);
			Node::Pointer SlaveNode2 = mrModelPart.pGetNode(idSlave2);

			VariablesList &rVariablesList = mrModelPart.GetNodalSolutionStepVariablesList();

			SizeType buffer_size = MasterNode->GetBufferSize();

			for (VariablesList::const_iterator i_variable = rVariablesList.begin(); i_variable != rVariablesList.end(); i_variable++)
			{
				// std::cout<<" name "<<i_variable->Name()<<std::endl;
				// std::cout<<" type "<<typeid(*i_variable).name()<<std::endl;
				std::string variable_name = i_variable->Name();
				if (KratosComponents<Variable<double>>::Has(variable_name))
				{
					// std::cout<<"double"<<std::endl;
					const Variable<double> &variable = KratosComponents<Variable<double>>::Get(variable_name);
					for (SizeType step = 0; step < buffer_size; step++)
					{
						// getting the data of the solution step
						double &node_data = MasterNode->FastGetSolutionStepValue(variable, step);

						double node0_data = SlaveNode1->FastGetSolutionStepValue(variable, step);
						double node1_data = SlaveNode2->FastGetSolutionStepValue(variable, step);

						node_data = (0.5 * node0_data + 0.5 * node1_data);
					}
				}
				else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name))
				{
					// std::cout<<"array1d"<<std::endl;
					const Variable<array_1d<double, 3>> &variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(variable_name);
					for (SizeType step = 0; step < buffer_size; step++)
					{
						// getting the data of the solution step
						array_1d<double, 3> &node_data = MasterNode->FastGetSolutionStepValue(variable, step);
						const array_1d<double, 3> &node0_data = SlaveNode1->FastGetSolutionStepValue(variable, step);
						const array_1d<double, 3> &node1_data = SlaveNode2->FastGetSolutionStepValue(variable, step);
						noalias(node_data) = (0.5 * node0_data + 0.5 * node1_data);
					}
				}
				else if (KratosComponents<Variable<int>>::Has(variable_name))
				{
					// std::cout<<"int"<<std::endl;
					// NO INTERPOLATION
				}
				else if (KratosComponents<Variable<bool>>::Has(variable_name))
				{
					// std::cout<<"bool"<<std::endl;
					// NO INTERPOLATION
				}
				else if (KratosComponents<Variable<Matrix>>::Has(variable_name))
				{
					// std::cout<<"Matrix"<<std::endl;
					const Variable<Matrix> &variable = KratosComponents<Variable<Matrix>>::Get(variable_name);
					for (SizeType step = 0; step < buffer_size; step++)
					{
						// getting the data of the solution step
						Matrix &node_data = MasterNode->FastGetSolutionStepValue(variable, step);

						Matrix &node0_data = SlaveNode1->FastGetSolutionStepValue(variable, step);
						Matrix &node1_data = SlaveNode2->FastGetSolutionStepValue(variable, step);

						if (node_data.size1() > 0 && node_data.size2())
						{
							if (node_data.size1() == node0_data.size1() && node_data.size2() == node0_data.size2() &&
								node_data.size1() == node1_data.size1() && node_data.size2() == node1_data.size2())
							{
								noalias(node_data) = (0.5 * node0_data + 0.5 * node1_data);
							}
						}
					}
				}
				else if (KratosComponents<Variable<Vector>>::Has(variable_name))
				{
					// std::cout<<"Vector"<<std::endl;
					const Variable<Vector> &variable = KratosComponents<Variable<Vector>>::Get(variable_name);
					for (SizeType step = 0; step < buffer_size; step++)
					{
						// getting the data of the solution step
						Vector &node_data = MasterNode->FastGetSolutionStepValue(variable, step);

						Vector &node0_data = SlaveNode1->FastGetSolutionStepValue(variable, step);
						Vector &node1_data = SlaveNode2->FastGetSolutionStepValue(variable, step);

						if (node_data.size() > 0)
						{
							if (node_data.size() == node0_data.size() &&
								node_data.size() == node1_data.size())
							{
								noalias(node_data) = (0.5 * node0_data + 0.5 * node1_data);
								// node_data = (0.5*node0_data + 0.5*node1_data);
							}
						}
					}
				}
			}

			KRATOS_CATCH("")
		}

		void SetVariablesToFreeSurfaceElements(ModelPart::NodesContainerType::const_iterator &in, const SizeType dimension)
		{
			KRATOS_TRY

			NodeWeakPtrVectorType &neighb_nodes = in->GetValue(NEIGHBOUR_NODES);
			array_1d<double, 3> sumOfCoordinates = in->Coordinates();
			array_1d<double, 3> sumOfCurrentVelocities = in->FastGetSolutionStepValue(VELOCITY, 0);
			array_1d<double, 3> sumOfPreviousVelocities = in->FastGetSolutionStepValue(VELOCITY, 1);
			double sumOfPressures = in->FastGetSolutionStepValue(PRESSURE, 0);
			double counter = 1.0;

			for (NodeWeakPtrVectorType::iterator nn = neighb_nodes.begin(); nn != neighb_nodes.end(); nn++)
			{
				counter += 1.0;
				noalias(sumOfCoordinates) += (nn)->Coordinates();
				noalias(sumOfCurrentVelocities) += (nn)->FastGetSolutionStepValue(VELOCITY, 0);
				noalias(sumOfPreviousVelocities) += (nn)->FastGetSolutionStepValue(VELOCITY, 1);
				sumOfPressures += (nn)->FastGetSolutionStepValue(PRESSURE, 0);
			}
			in->X() = sumOfCoordinates[0] / counter;
			in->Y() = sumOfCoordinates[1] / counter;
			in->X0() = sumOfCoordinates[0] / counter;
			in->Y0() = sumOfCoordinates[1] / counter;
			in->FastGetSolutionStepValue(DISPLACEMENT_X, 0) = 0;
			in->FastGetSolutionStepValue(DISPLACEMENT_Y, 0) = 0;
			in->FastGetSolutionStepValue(DISPLACEMENT_X, 1) = 0;
			in->FastGetSolutionStepValue(DISPLACEMENT_Y, 1) = 0;
			in->FastGetSolutionStepValue(VELOCITY_X, 0) = sumOfCurrentVelocities[0] / counter;
			in->FastGetSolutionStepValue(VELOCITY_Y, 0) = sumOfCurrentVelocities[1] / counter;
			in->FastGetSolutionStepValue(VELOCITY_X, 1) = sumOfPreviousVelocities[0] / counter;
			in->FastGetSolutionStepValue(VELOCITY_Y, 1) = sumOfPreviousVelocities[1] / counter;
			in->FastGetSolutionStepValue(PRESSURE, 0) = sumOfPressures / counter;
			if (dimension == 3)
			{
				in->Z() = sumOfCoordinates[2] / counter;
				in->Z0() = sumOfCoordinates[2] / counter;
				in->FastGetSolutionStepValue(DISPLACEMENT_Z, 0) = 0;
				in->FastGetSolutionStepValue(DISPLACEMENT_Z, 1) = 0;
				in->FastGetSolutionStepValue(VELOCITY_Z, 0) = sumOfCurrentVelocities[2] / counter;
				in->FastGetSolutionStepValue(VELOCITY_Z, 1) = sumOfPreviousVelocities[2] / counter;
			}

			KRATOS_CATCH("")
		}

		void SetMeshSizeNearBoundaries(double &meshSize,
									   double rigidNodeLocalMeshSize,
									   double rigidNodeMeshCounter)
		{
			KRATOS_TRY

			if (rigidNodeMeshCounter > 0)
			{
				const double rigidWallMeshSize = rigidNodeLocalMeshSize / rigidNodeMeshCounter;
				const double ratio = rigidWallMeshSize / meshSize;
				const double tolerance = 2.0;
				if (ratio > tolerance)
				{
					meshSize *= 0.5;
					meshSize += 0.5 * rigidWallMeshSize;
				}
			}

			KRATOS_CATCH("")
		}

		void CriticalVolumeForRefinedBoxes(Element::GeometryType eElement,
										   double &criticalVolume,
										   double &safetyCoefficient3D)
		{

			const SizeType numberOfRefiningBoxes = mrRemesh.UseRefiningBox.size();
			for (SizeType index = 0; index < numberOfRefiningBoxes; index++)
			{
				double currentTime = mrModelPart.GetProcessInfo()[TIME];
				double initialTime = mrRemesh.RefiningBoxInitialTime[index];
				double finalTime = mrRemesh.RefiningBoxFinalTime[index];
				bool refiningBox = mrRemesh.UseRefiningBox[index];

				if (refiningBox == true && currentTime > initialTime && currentTime < finalTime)
				{
					double baricenterX = 0.25 * (eElement[0].X() + eElement[1].X() + eElement[2].X() + eElement[3].X());
					double baricenterY = 0.25 * (eElement[0].Y() + eElement[1].Y() + eElement[2].Y() + eElement[3].Y());
					double baricenterZ = 0.25 * (eElement[0].Z() + eElement[1].Z() + eElement[2].Z() + eElement[3].Z());

					double meanMeshSize = mrRemesh.RefiningBoxMeshSize[index] * 0.5 + mrRemesh.Refine->CriticalRadius * 0.5;
					array_1d<double, 3> minExternalPoint = mrRemesh.RefiningBoxShiftedMinimumPoint[index];
					array_1d<double, 3> maxExternalPoint = mrRemesh.RefiningBoxShiftedMaximumPoint[index];
					array_1d<double, 3> RefiningBoxMinimumPoint = mrRemesh.RefiningBoxMinimumPoint[index];
					array_1d<double, 3> RefiningBoxMaximumPoint = mrRemesh.RefiningBoxMaximumPoint[index];
					if (baricenterX > RefiningBoxMinimumPoint[0] && baricenterX < RefiningBoxMaximumPoint[0] &&
						baricenterY > RefiningBoxMinimumPoint[1] && baricenterY < RefiningBoxMaximumPoint[1] &&
						baricenterZ > RefiningBoxMinimumPoint[2] && baricenterZ < RefiningBoxMaximumPoint[2])
					{
						criticalVolume = 0.01 * (std::pow(mrRemesh.RefiningBoxMeshSize[index], 3) / (6.0 * std::sqrt(2))); // mean Volume of a regular tetrahedral per node with 0.01 of penalization
					}
					else if ((baricenterX < RefiningBoxMinimumPoint[0] && baricenterX > minExternalPoint[0]) || (baricenterX > RefiningBoxMaximumPoint[0] && baricenterX < maxExternalPoint[0]) ||
							 (baricenterY < RefiningBoxMinimumPoint[1] && baricenterY > minExternalPoint[1]) || (baricenterY > RefiningBoxMaximumPoint[1] && baricenterY < maxExternalPoint[1]) ||
							 (baricenterZ < RefiningBoxMinimumPoint[2] && baricenterZ > minExternalPoint[2]) || (baricenterZ > RefiningBoxMaximumPoint[2] && baricenterZ < maxExternalPoint[2])) // transition zone
					{
						criticalVolume = 0.005 * (std::pow(meanMeshSize, 3) / (6.0 * std::sqrt(2)));
						safetyCoefficient3D *= 0.5;
					}
					else
					{
						criticalVolume = 0.01 * (std::pow(meanMeshSize, 3) / (6.0 * std::sqrt(2)));
					}
				}
			}
		}

		void CheckSliverElements(Element::GeometryType eElement,
								 SizeType &inside_nodes_removed,
								 SizeType &nodes_removed_inlet_zone,
								 bool inletElement,
								 SizeType notRigidNodeId,
								 const SizeType principalModelPartId,
								 double &wallNodeDistance,
								 array_1d<double, 3> &WallBaricenter,
								 std::vector<array_1d<double, 3>> &rigidNode,
								 std::vector<array_1d<double, 3>> &rigidNodesNormals,
								 array_1d<double, 3> &notRigidNode)
		{
			KRATOS_TRY

			double a1 = 0; // slope x,y,z for the plane composed by rigid nodes only
			double b1 = 0;
			double c1 = 0;
			a1 = (rigidNode[1][1] - rigidNode[0][1]) * (rigidNode[2][2] - rigidNode[0][2]) - (rigidNode[2][1] - rigidNode[0][1]) * (rigidNode[1][2] - rigidNode[0][2]);
			b1 = (rigidNode[1][2] - rigidNode[0][2]) * (rigidNode[2][0] - rigidNode[0][0]) - (rigidNode[2][2] - rigidNode[0][2]) * (rigidNode[1][0] - rigidNode[0][0]);
			c1 = (rigidNode[1][0] - rigidNode[0][0]) * (rigidNode[2][1] - rigidNode[0][1]) - (rigidNode[2][0] - rigidNode[0][0]) * (rigidNode[1][1] - rigidNode[0][1]);
			double a2 = 0; // slope x,y,z between 2 rigid nodes and a not rigid node
			double b2 = 0;
			double c2 = 0;
			a2 = (rigidNode[1][1] - rigidNode[0][1]) * (notRigidNode[2] - rigidNode[0][2]) - (notRigidNode[1] - rigidNode[0][1]) * (rigidNode[1][2] - rigidNode[0][2]);
			b2 = (rigidNode[1][2] - rigidNode[0][2]) * (notRigidNode[0] - rigidNode[0][0]) - (notRigidNode[2] - rigidNode[0][2]) * (rigidNode[1][0] - rigidNode[0][0]);
			c2 = (rigidNode[1][0] - rigidNode[0][0]) * (notRigidNode[1] - rigidNode[0][1]) - (notRigidNode[0] - rigidNode[0][0]) * (rigidNode[1][1] - rigidNode[0][1]);
			double a3 = 0; // slope x,y,z between 2 rigid nodes and a not rigid node
			double b3 = 0;
			double c3 = 0;
			a3 = (rigidNode[1][1] - rigidNode[2][1]) * (notRigidNode[2] - rigidNode[2][2]) - (notRigidNode[1] - rigidNode[2][1]) * (rigidNode[1][2] - rigidNode[2][2]);
			b3 = (rigidNode[1][2] - rigidNode[2][2]) * (notRigidNode[0] - rigidNode[2][0]) - (notRigidNode[2] - rigidNode[2][2]) * (rigidNode[1][0] - rigidNode[2][0]);
			c3 = (rigidNode[1][0] - rigidNode[2][0]) * (notRigidNode[1] - rigidNode[2][1]) - (notRigidNode[0] - rigidNode[2][0]) * (rigidNode[1][1] - rigidNode[2][1]);
			double a4 = 0; // slope x,y,z between 2 rigid nodes and a not rigid node
			double b4 = 0;
			double c4 = 0;
			a4 = (rigidNode[0][1] - rigidNode[2][1]) * (notRigidNode[2] - rigidNode[2][2]) - (notRigidNode[1] - rigidNode[2][1]) * (rigidNode[0][2] - rigidNode[2][2]);
			b4 = (rigidNode[0][2] - rigidNode[2][2]) * (notRigidNode[0] - rigidNode[2][0]) - (notRigidNode[2] - rigidNode[2][2]) * (rigidNode[0][0] - rigidNode[2][0]);
			c4 = (rigidNode[0][0] - rigidNode[2][0]) * (notRigidNode[1] - rigidNode[2][1]) - (notRigidNode[0] - rigidNode[2][0]) * (rigidNode[0][1] - rigidNode[2][1]);

			// angle between the plane composed by rigid nodes only and the other plans. If the angle is small, the particle can pass through the wall
			double cosAngle12 = (a1 * a2 + b1 * b2 + c1 * c2) / (std::sqrt(std::pow(a1, 2) + std::pow(b1, 2) + std::pow(c1, 2)) * std::sqrt(std::pow(a2, 2) + std::pow(b2, 2) + std::pow(c2, 2)));
			double cosAngle13 = (a1 * a3 + b1 * b3 + c1 * c3) / (std::sqrt(std::pow(a1, 2) + std::pow(b1, 2) + std::pow(c1, 2)) * std::sqrt(std::pow(a3, 2) + std::pow(b3, 2) + std::pow(c3, 2)));
			double cosAngle14 = (a1 * a4 + b1 * b4 + c1 * c4) / (std::sqrt(std::pow(a1, 2) + std::pow(b1, 2) + std::pow(c1, 2)) * std::sqrt(std::pow(a4, 2) + std::pow(b4, 2) + std::pow(c4, 2)));

			// angle between the normals of the rigid nodes. I want to avoid rigid elements at the corner
			double cosAngleBetweenNormals01 = (rigidNodesNormals[0][0] * rigidNodesNormals[1][0] + rigidNodesNormals[0][1] * rigidNodesNormals[1][1]) /
											  (std::sqrt(std::pow(rigidNodesNormals[0][0], 2) + std::pow(rigidNodesNormals[0][1], 2)) *
											   std::sqrt(std::pow(rigidNodesNormals[1][0], 2) + std::pow(rigidNodesNormals[1][1], 2)));
			double cosAngleBetweenNormals02 = (rigidNodesNormals[0][0] * rigidNodesNormals[2][0] + rigidNodesNormals[0][1] * rigidNodesNormals[2][1]) /
											  (std::sqrt(std::pow(rigidNodesNormals[0][0], 2) + std::pow(rigidNodesNormals[0][1], 2)) *
											   std::sqrt(std::pow(rigidNodesNormals[2][0], 2) + std::pow(rigidNodesNormals[2][1], 2)));
			double cosAngleBetweenNormals12 = (rigidNodesNormals[1][0] * rigidNodesNormals[2][0] + rigidNodesNormals[1][1] * rigidNodesNormals[2][1]) /
											  (std::sqrt(std::pow(rigidNodesNormals[1][0], 2) + std::pow(rigidNodesNormals[1][1], 2)) *
											   std::sqrt(std::pow(rigidNodesNormals[2][0], 2) + std::pow(rigidNodesNormals[2][1], 2)));

			if ((std::abs(cosAngle12) > 0.995 || std::abs(cosAngle13) > 0.995 || std::abs(cosAngle14) > 0.995) && (cosAngleBetweenNormals01 > 0.99 && cosAngleBetweenNormals02 > 0.99 && cosAngleBetweenNormals12 > 0.99))
			{
				eElement[notRigidNodeId].Set(TO_ERASE);
				inside_nodes_removed++;

				if (inletElement)
				{
					nodes_removed_inlet_zone++;
				}

				const SizeType propertyIdNode = eElement[notRigidNodeId].FastGetSolutionStepValue(PROPERTY_ID);
				if (propertyIdNode != principalModelPartId) // this is to conserve the number of nodes of the smaller domain in case of a two-fluid analysis
				{
					mrRemesh.Info->BalancePrincipalSecondaryPartsNodes += -1;
				}
			}

			double pwdDistance = 0.0f;
			for (std::size_t i = 0; i < 3; i++)
			{
				pwdDistance += std::pow(eElement[notRigidNodeId].Coordinates()[i] - WallBaricenter[i], 2);
			}
			wallNodeDistance = std::sqrt(pwdDistance);

			KRATOS_CATCH("")
		}

		void DetectRemovableNodes(Element::GeometryType eElement,
								  array_1d<double, 6> &Edges,
								  array_1d<SizeType, 6> &FirstEdgeNode,
								  array_1d<SizeType, 6> &SecondEdgeNode,
								  int &removableNode,
								  double &wallLength,
								  double &minimumLength,
								  double &safetyCoefficient3D)
		{
			KRATOS_TRY

			// ////////  to compute the length of the wall edge /////////
			array_1d<double, 3> CoorDifference = eElement[1].Coordinates() - eElement[0].Coordinates();

			double SquaredLength = CoorDifference[0] * CoorDifference[0] +
								   CoorDifference[1] * CoorDifference[1] +
								   CoorDifference[2] * CoorDifference[2];
			Edges[0] = std::sqrt(SquaredLength);
			minimumLength = Edges[0] * 10.0;
			FirstEdgeNode[0] = 0;
			SecondEdgeNode[0] = 1;
			if (eElement[0].Is(RIGID) && eElement[1].Is(RIGID))
			{
				wallLength = Edges[0];
			}
			if ((eElement[0].Is(RIGID) && eElement[1].IsNot(RIGID)) ||
				(eElement[1].Is(RIGID) && eElement[0].IsNot(RIGID)))
			{
				minimumLength = Edges[0];
				if (eElement[0].IsNot(RIGID))
				{
					removableNode = 0;
				}
				else
				{
					removableNode = 1;
				}
			}
			SizeType counter = 0;
			for (SizeType i = 2; i < eElement.size(); i++)
			{
				for (SizeType j = 0; j < i; j++)
				{
					noalias(CoorDifference) = eElement[i].Coordinates() - eElement[j].Coordinates();
					SquaredLength = CoorDifference[0] * CoorDifference[0] +
									CoorDifference[1] * CoorDifference[1] +
									CoorDifference[2] * CoorDifference[2];
					counter += 1;
					Edges[counter] = std::sqrt(SquaredLength);
					FirstEdgeNode[counter] = j;
					SecondEdgeNode[counter] = i;
					if (eElement[i].Is(RIGID) && eElement[j].Is(RIGID) && (wallLength == 0 || Edges[counter] > wallLength))
					{
						wallLength = Edges[counter];
					}
					if (((eElement[i].Is(RIGID) && eElement[j].IsNot(RIGID)) ||
						 (eElement[j].Is(RIGID) && eElement[i].IsNot(RIGID))) &&
						eElement[i].IsNot(TO_ERASE) && eElement[j].IsNot(TO_ERASE) &&
						(Edges[counter] < minimumLength || minimumLength == 0))
					{
						minimumLength = Edges[counter];
						if (eElement[i].IsNot(RIGID))
						{
							removableNode = i;
						}
						else
						{
							removableNode = j;
						}
					}
				}
			}

			////////  to avoid the elimination of isolated free-surface-rigid elements /////////
			for (SizeType i = 0; i < eElement.size(); i++)
			{
				if (eElement[i].IsNot(RIGID) && eElement[i].IsNot(TO_ERASE) && eElement[i].IsNot(SOLID) && eElement[i].IsNot(ISOLATED))
				{
					//////it is evident when a freesurface particle in touch with wall is erased --> reduce the safety coeff
					if (eElement[i].Is(FREE_SURFACE))
					{
						NodeWeakPtrVectorType &neighb_nodes = eElement[i].GetValue(NEIGHBOUR_NODES);
						SizeType countRigid = 0;
						SizeType countFreeSurface = 0;
						for (NodeWeakPtrVectorType::iterator nn = neighb_nodes.begin(); nn != neighb_nodes.end(); nn++)
						{
							if ((nn)->Is(RIGID) || (nn)->Is(SOLID))
							{
								countRigid++;
							}
							if ((nn)->Is(FREE_SURFACE) && (nn)->IsNot(RIGID) && (nn)->IsNot(SOLID))
							{
								countFreeSurface++;
							}
						}
						if ((countRigid + countFreeSurface) == neighb_nodes.size() && countRigid > 0)
						{
							safetyCoefficient3D = 0.25;
						}
					}
				}
			}
			KRATOS_CATCH("")
		}

		//**************************************************************************
	private:
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
					if ((ie->GetGeometry()[i].Is(RIGID) && ie->GetGeometry()[i].IsNot(PFEMFlags::LAGRANGIAN_INLET)) || ie->GetGeometry()[i].Is(SOLID) || ie->GetGeometry()[i].Is(PFEMFlags::EULERIAN_INLET))
					{
						rigidNodes++;
					}
				}
				if (dimension == 2)
				{
					if (rigidNodes > 0)
						EraseCriticalNodes2D(ie->GetGeometry(), inside_nodes_removed, nodes_removed_inlet_zone);
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

				SetMeshSizeNearBoundaries(meshSize, rigidNodeLocalMeshSize, rigidNodeMeshCounter);

				const double size_for_distance_boundary = 0.6 * meshSize;
				const double size_for_wall_tip_contact_side = 0.15 * mrRemesh.Refine->CriticalSide;

				if (in->Is(TO_ERASE))
				{
					some_node_is_removed = true;
				}
				bool on_contact_tip = false;

				if (in->Is(TO_SPLIT) || in->Is(CONTACT))
					on_contact_tip = true;

				if (in->IsNot(NEW_ENTITY) && in->IsNot(INLET) && in->IsNot(ISOLATED) && in->IsNot(RIGID) && in->IsNot(SOLID))
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

					n_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, neighbours.begin(), neighbour_distances.begin(), num_neighbours);

					if (n_points_in_radius > 1 && neighErasedNodes == 0 && in->IsNot(INLET))
					{

						if (in->IsNot(RIGID) && in->IsNot(SOLID) && in->IsNot(ISOLATED))
						{
							if (mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_NODES_ON_DISTANCE))
							{

								if (in->IsNot(FREE_SURFACE) && in->IsNot(RIGID) && freeSurfaceNeighNodes == dimension)
								{
									SetVariablesToFreeSurfaceElements(in, dimension);
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

							if (counter > 1 && in->IsNot(RIGID) && in->IsNot(SOLID) && in->IsNot(NEW_ENTITY) && !on_contact_tip)
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

			CriticalVolumeForRefinedBoxes(eElement, criticalVolume, safetyCoefficient3D);

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

					if (eElement[i].IsNot(RIGID) && eElement[i].IsNot(SOLID) && eElement[i].IsNot(TO_ERASE))
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
				CheckSliverElements(eElement,
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

			bool longDamBreak = false; // to attivate in case of long dam breaks to avoid separated elements in the water front
			if (longDamBreak == true && freeSurfaceNodes > 2 && rigidNodes > 1)
			{
				ControlForLongDamBreakProblems3D(eElement, inletElement, inside_nodes_removed, nodes_removed_inlet_zone);
			}

			array_1d<double, 6> Edges(6, 0.0);
			array_1d<SizeType, 6> FirstEdgeNode(6, 0);
			array_1d<SizeType, 6> SecondEdgeNode(6, 0);
			int removableNode = -1;
			double wallLength = 0;
			double minimumLength = 0;
			DetectRemovableNodes(eElement, Edges, FirstEdgeNode, SecondEdgeNode, removableNode, wallLength, minimumLength, safetyCoefficient3D);

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
			CompareEdgeLengthsNearWalls(eElement, Edges, FirstEdgeNode, SecondEdgeNode, freeSurfaceNodes, safetyCoefficient3D, wallLength, inletElement, inside_nodes_removed, nodes_removed_inlet_zone);

			KRATOS_CATCH("")
		}

		/// Assignment operator.
		RemoveMeshNodesForFluidsProcess &operator=(RemoveMeshNodesForFluidsProcess const &rOther);

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
									RemoveMeshNodesForFluidsProcess &rThis);

	/// output stream function
	inline std::ostream &operator<<(std::ostream &rOStream,
									const RemoveMeshNodesForFluidsProcess &rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@}

} // namespace Kratos.

#endif // KRATOS_REMOVE_MESH_NODES_FOR_FLUIDS_PROCESS_H_INCLUDED  defined
