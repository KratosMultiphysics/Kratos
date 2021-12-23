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

///VARIABLES used:
//Data:     NORMAL, MASTER_NODES, NEIGHBOUR_NODES, NEIGBOUR_ELEMENTS
//Flags:    (checked) TO_ERASE, BOUNDARY, STRUCTURE, TO_SPLIT, CONTACT, NEW_ENTITY, BLOCKED
//          (set)     TO_ERASE(conditions,nodes)(set), NEW_ENTITY(conditions,nodes)(set), BLOCKED(nodes)->locally, VISITED(nodes)(set)
//          (modified)
//          (reset)   BLOCKED->locally
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
		typedef Bucket<3, Node<3>, std::vector<Node<3>::Pointer>, Node<3>::Pointer, std::vector<Node<3>::Pointer>::iterator, std::vector<double>::iterator> BucketType;
		typedef Tree<KDTreePartition<BucketType>> KdtreeType; //Kdtree
		typedef ModelPart::MeshType::GeometryType::PointsArrayType PointsArrayType;

		typedef GlobalPointersVector<Node<3>> NodeWeakPtrVectorType;
		typedef GlobalPointersVector<Element> ElementWeakPtrVectorType;
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

			// double NumberOfNodes = mrModelPart.NumberOfNodes();

			bool any_node_removed = false;

			int error_nodes_removed = 0;
			int inside_nodes_removed = 0;
			int boundary_nodes_removed = 0;

			//if the remove_node switch is activated, we check if the nodes got too close
			if (mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_NODES))
			{
				if (mEchoLevel > 1)
					std::cout << " REMOVE_NODES is TRUE " << std::endl;
				// bool any_node_removed_on_error = false;
				// ////////////////////////////////////////////////////////////
				// if (mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_NODES_ON_ERROR))
				// {
				// 	if (mEchoLevel > 1)
				// 		std::cout << " REMOVE_NODES_ON_ERROR is TRUE " << std::endl;

				// 	any_node_removed_on_error = RemoveNodesOnError(error_nodes_removed); //2D and 3D
				// }
				////////////////////////////////////////////////////////////
				if (mEchoLevel > 1)
					std::cout << "error_nodes_removed :" << error_nodes_removed << std::endl;

				bool any_node_removed_on_distance = false;
				////////////////////////////////////////////////////////////
				if (mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_NODES_ON_DISTANCE))
				{
					if (mEchoLevel > 1)
						std::cout << " REMOVE_NODES_ON_DISTANCE is TRUE " << std::endl;
					// double  MeanRadius=0;
					any_node_removed_on_distance = RemoveNodesOnDistance(inside_nodes_removed, boundary_nodes_removed);
				}
				// REMOVE ON DISTANCE
				////////////////////////////////////////////////////////////

				if (any_node_removed_on_distance)
					any_node_removed = true;

				if (any_node_removed || mrRemesh.UseBoundingBox == true)
					this->CleanRemovedNodes(mrModelPart);
			}

			// number of removed nodes:
			// mrRemesh.Info->RemovedNodes = NumberOfNodes - mrModelPart.NumberOfNodes();
			mrRemesh.Info->RemovedNodes += error_nodes_removed + inside_nodes_removed + boundary_nodes_removed;
			int distance_remove = inside_nodes_removed + boundary_nodes_removed;

			if (mEchoLevel > 1)
			{
				std::cout << "   [ NODES      ( removed : " << mrRemesh.Info->RemovedNodes << " ) ]" << std::endl;
				std::cout << "   [ Error(removed: " << error_nodes_removed << "); Distance(removed: " << distance_remove << "; inside: " << inside_nodes_removed << "; boundary: " << boundary_nodes_removed << ") ]" << std::endl;

				//std::cout<<"   Nodes after  erasing : "<<mrModelPart.Nodes().size()<<std::endl;
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

	private:
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

			//MESH 0 total domain mesh
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

					/////////////////////////////////////////// here for BOUNDING BOX ///////////////////////////////////////////
				}
			}

			rModelPart.Nodes().Sort();

			KRATOS_CATCH("")
		}

		//**************************************************************************
		//**************************************************************************

		// bool RemoveNodesOnError(int &error_removed_nodes)
		// {
		// 	KRATOS_TRY

		// 	//***SIZES :::: parameters do define the tolerance in mesh size:
		// 	double size_for_criterion_error = 2.0 * mrRemesh.Refine->CriticalRadius; //compared with mean node radius

		// 	bool any_node_removed = false;

		// 	MeshErrorCalculationUtilities MeshErrorDistribution;
		// 	MeshErrorDistribution.SetEchoLevel(mEchoLevel);

		// 	std::vector<double> NodalError;
		// 	std::vector<int> nodes_ids;

		// 	MeshErrorDistribution.NodalErrorCalculation(mrModelPart, NodalError, nodes_ids, mrRemesh.Refine->GetErrorVariable());

		// 	for (ModelPart::NodesContainerType::const_iterator in = mrModelPart.NodesBegin(); in != mrModelPart.NodesEnd(); in++)
		// 	{

		// 		NodeWeakPtrVectorType &rN = in->GetValue(NEIGHBOUR_NODES);
		// 		int erased_nodes = 0;
		// 		for (unsigned int i = 0; i < rN.size(); i++)
		// 		{
		// 			if (rN[i].Is(TO_ERASE))
		// 				erased_nodes += 1;
		// 		}

		// 		if (in->IsNot(BOUNDARY) && in->IsNot(STRUCTURE) && erased_nodes < 1)
		// 		{
		// 			double &MeanError = in->FastGetSolutionStepValue(MEAN_ERROR);
		// 			MeanError = NodalError[nodes_ids[in->Id()]];

		// 			ElementWeakPtrVectorType &neighb_elems = in->GetValue(NEIGHBOUR_ELEMENTS);
		// 			double mean_node_radius = 0;
		// 			for (ElementWeakPtrVectorType::iterator ne = neighb_elems.begin(); ne != neighb_elems.end(); ne++)
		// 			{
		// 				mean_node_radius += mMesherUtilities.CalculateElementRadius((ne)->GetGeometry()); //Triangle 2D, Tetrahedron 3D
		// 																								  //mean_node_radius+= mMesherUtilities.CalculateTriangleRadius((ne)->GetGeometry());
		// 																								  //mean_node_radius+= mMesherUtilities.CalculateTetrahedronRadius((ne)->GetGeometry());
		// 			}

		// 			mean_node_radius /= double(neighb_elems.size());

		// 			if (NodalError[nodes_ids[in->Id()]] < mrRemesh.Refine->ReferenceError && mean_node_radius < size_for_criterion_error)
		// 			{
		// 				in->Set(TO_ERASE);
		// 				any_node_removed = true;
		// 				error_removed_nodes++;
		// 			}
		// 		}
		// 	}

		// 	return any_node_removed;

		// 	KRATOS_CATCH(" ")
		// }

		bool RemoveNodesOnDistance(int &inside_nodes_removed, int &boundary_nodes_removed)
		{
			KRATOS_TRY

			const unsigned int dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

			//***SIZES :::: parameters do define the tolerance in mesh size:

			// if(dimension==3){
			//   size_for_distance_inside       = 0.5 * initialMeanRadius;//compared to element radius
			//   size_for_distance_boundary     = 0.5 * initialMeanRadius; //compared to element radius
			//   size_for_wall_tip_contact_side = 0.15 * mrRemesh.Refine->CriticalSide;
			// }

			bool derefine_wall_tip_contact = false;

			bool any_node_removed = false;

			//bucket size definition:
			unsigned int bucket_size = 20;

			//create the list of the nodes to be check during the search
			std::vector<Node<3>::Pointer> list_of_nodes;
			list_of_nodes.reserve(mrModelPart.NumberOfNodes());
			for (ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++)
			{
				(list_of_nodes).push_back(*(i_node.base()));
			}

			KdtreeType nodes_tree(list_of_nodes.begin(), list_of_nodes.end(), bucket_size);

			////////////////////////////////////////////////////////////

			//all of the nodes in this list will be preserved
			unsigned int num_neighbours = 100;

			std::vector<Node<3>::Pointer> neighbours(num_neighbours);
			std::vector<double> neighbour_distances(num_neighbours);

			//radius means the distance, if the distance between two nodes is closer to radius -> mark for removing
			double radius = 0;
			Node<3> work_point(0, 0.0, 0.0, 0.0);
			unsigned int n_points_in_radius;

			const ProcessInfo &rCurrentProcessInfo = mrModelPart.GetProcessInfo();
			double currentTime = rCurrentProcessInfo[TIME];

			double initialMeanRadius = mrRemesh.Refine->CriticalRadius;

			double initialTimeForRefinement = mrRemesh.RefiningBoxInitialTime;
			double finalTimeForRefinement = mrRemesh.RefiningBoxFinalTime;
			bool refiningBox = mrRemesh.UseRefiningBox;

			if (!(refiningBox == true && currentTime > initialTimeForRefinement && currentTime < finalTimeForRefinement))
			{
				refiningBox = false;
			}

			unsigned int erased_nodes = 0;
			for (ModelPart::ElementsContainerType::const_iterator ie = mrModelPart.ElementsBegin();
				 ie != mrModelPart.ElementsEnd(); ie++)
			{
				unsigned int rigidNodes = 0;
				//coordinates
				for (unsigned int i = 0; i < ie->GetGeometry().size(); i++)
				{
					if ((ie->GetGeometry()[i].Is(RIGID) && ie->GetGeometry()[i].IsNot(INLET)) || ie->GetGeometry()[i].Is(SOLID))
					{
						rigidNodes++;
					}
				}

				if (dimension == 2)
				{
					if (rigidNodes > 0)
						EraseCriticalNodes2D(ie->GetGeometry(), erased_nodes, inside_nodes_removed);
				}
				else if (dimension == 3)
				{
					if (rigidNodes > 1)
						EraseCriticalNodes3D(ie->GetGeometry(), erased_nodes, inside_nodes_removed, rigidNodes);
				}
			}

			for (ModelPart::NodesContainerType::const_iterator in = mrModelPart.NodesBegin(); in != mrModelPart.NodesEnd(); in++)
			{

				if (refiningBox == true)
				{
					array_1d<double, 3> NodeCoordinates = in->Coordinates();
					if (dimension == 2)
					{
						initialMeanRadius = SetMeshSizeInMeshRefinementArea(NodeCoordinates);
					}
					else if (dimension == 3)
					{
						initialMeanRadius = SetMeshSizeInMeshRefinementVolume(NodeCoordinates);
					}
				}

				double size_for_distance_boundary = 0.6 * initialMeanRadius;
				double size_for_wall_tip_contact_side = 0.15 * mrRemesh.Refine->CriticalSide;

				if (in->Is(TO_ERASE))
				{
					any_node_removed = true;
				}
				bool on_contact_tip = false;

				if (in->Is(TO_SPLIT) || in->Is(CONTACT))
					on_contact_tip = true;

				if (in->IsNot(NEW_ENTITY) && in->IsNot(INLET) && in->IsNot(ISOLATED))
				// if( in->IsNot(NEW_ENTITY) )
				{
					unsigned int neighErasedNodes = 0;
					radius = 0.6 * initialMeanRadius;

					work_point[0] = in->X();
					work_point[1] = in->Y();
					work_point[2] = in->Z();
					unsigned int freeSurfaceNeighNodes = 0;
					// unsigned int rigidNeighNodes=0;

					if (in->Is(FREE_SURFACE))
					{
						// it must be more difficult to erase a free_surface node, otherwise, lot of volume is lost
						// this value has a strong effect on volume variation due to remeshing
						radius = 0.475 * initialMeanRadius; //compared with element radius
						// radius = 0.4  * initialMeanRadius;//compared with element radius
						NodeWeakPtrVectorType &neighb_nodes = in->GetValue(NEIGHBOUR_NODES);
						unsigned int countRigid = 0;
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
							radius = 0.15 * initialMeanRadius;
						}
					}
					else
					{
						NodeWeakPtrVectorType &neighb_nodes = in->GetValue(NEIGHBOUR_NODES);
						for (NodeWeakPtrVectorType::iterator nn = neighb_nodes.begin(); nn != neighb_nodes.end(); nn++)
						{
							if ((nn)->Is(FREE_SURFACE))
							{
								freeSurfaceNeighNodes++;
							}
							if ((nn)->Is(TO_ERASE))
							{
								neighErasedNodes++;
							}
							// if((nn)->Is(RIGID)){
							//   rigidNeighNodes++;
							// }
						}
					}

					if (freeSurfaceNeighNodes > 1)
					{
						radius = 0.5 * initialMeanRadius;
					}

					if (in->Is(INLET))
					{
						radius = 0.3 * initialMeanRadius; //compared with element radius
					}
					n_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, neighbours.begin(), neighbour_distances.begin(), num_neighbours);

					if (n_points_in_radius > 1 && neighErasedNodes == 0)
					{

						if (in->IsNot(INLET) && in->IsNot(RIGID) && in->IsNot(SOLID) && in->IsNot(ISOLATED))
						{
							if (mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_NODES_ON_DISTANCE))
							{

								// if (in->IsNot(FREE_SURFACE) && in->IsNot(RIGID) && (freeSurfaceNeighNodes==dimension || rigidNeighNodes==dimension)){
								if (in->IsNot(FREE_SURFACE) && in->IsNot(RIGID) && freeSurfaceNeighNodes == dimension)
								{
									NodeWeakPtrVectorType &neighb_nodes = in->GetValue(NEIGHBOUR_NODES);
									array_1d<double, 3> sumOfCoordinates = in->Coordinates();
									array_1d<double, 3> sumOfCurrentVelocities = in->FastGetSolutionStepValue(VELOCITY, 0);
									array_1d<double, 3> sumOfPreviousVelocities = in->FastGetSolutionStepValue(VELOCITY, 1);
									double sumOfPressures = in->FastGetSolutionStepValue(PRESSURE, 0);
									double counter = 1.0;
									// std::cout<<"I WAS GOING TO ERASE THIS NODE: "<<std::endl;
									// std::cout<<"sumOfCoordinates: "<< sumOfCoordinates<<std::endl;
									// std::cout<<"sumOfCurrentVelocities: "<< sumOfCurrentVelocities<<std::endl;
									// std::cout<<"sumOfPressures: "<< sumOfPressures<<std::endl;
									// std::cout<<"freeSurfaceNeighNodes "<<freeSurfaceNeighNodes<<"   rigidNeighNodes "<<rigidNeighNodes<<std::endl;
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
									// std::cout<<"NOW ITS VARIABLES ARE: "<<std::endl;
									// std::cout<<"Coordinates: "<< in->X()<<" " << in->Y()<<" " << in->Z()<<std::endl;
									// std::cout<<"CurrentVelocities: "<<in->FastGetSolutionStepValue(VELOCITY_X,0)<<" "<<in->FastGetSolutionStepValue(VELOCITY_Y,0)<<" "<<in->FastGetSolutionStepValue(VELOCITY_Z,0)  <<std::endl;
									// std::cout<<"Pressures: "<<in->FastGetSolutionStepValue(PRESSURE,0) <<std::endl;
								}
								else
								{
									//look if we are already erasing any of the other nodes
									unsigned int contact_nodes = 0;
									for (std::vector<Node<3>::Pointer>::iterator nn = neighbours.begin(); nn != neighbours.begin() + n_points_in_radius; nn++)
									{
										if ((*nn)->Is(BOUNDARY) && (*nn)->Is(CONTACT))
											contact_nodes += 1;
									}

									if (contact_nodes < 1)
									{ //we release the node if no other nodes neighbours are being erased
										in->Set(TO_ERASE);
										any_node_removed = true;
										inside_nodes_removed++;
										//distance_remove++;
									}
								}
							}
						}
						else if (in->IsNot(INLET))
						{
							// else 	 {

							// std::cout<<"  Remove close boundary nodes: Candidate ["<<in->Id()<<"]"<<std::endl;
							//here we loop over the neighbouring nodes and if there are nodes
							//with BOUNDARY flag and closer than 0.2*nodal_h from our node, we remove the node we are considering
							unsigned int k = 0;
							unsigned int counter = 0;
							for (std::vector<Node<3>::Pointer>::iterator nn = neighbours.begin(); nn != neighbours.begin() + n_points_in_radius; nn++)
							{
								bool nn_on_contact_tip = false;

								if ((*nn)->Is(TO_SPLIT) || (*nn)->Is(CONTACT))
									nn_on_contact_tip = true;

								if ((*nn)->Is(BOUNDARY) && !nn_on_contact_tip && neighbour_distances[k] < size_for_distance_boundary && neighbour_distances[k] > 0.0)
								{
									//KRATOS_WATCH( neighbours_distances[k] )
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
							{ //Can be inserted in the boundary refine
								in->Set(TO_ERASE);
								if (mEchoLevel > 1)
									std::cout << "     Removed Boundary Node [" << in->Id() << "] on Distance " << std::endl;
								any_node_removed = true;
								boundary_nodes_removed++;
								//distance_remove ++;
							}
							else if (counter > 2 && in->IsNot(RIGID) && in->IsNot(SOLID) && in->IsNot(NEW_ENTITY) && on_contact_tip && derefine_wall_tip_contact)
							{
								in->Set(TO_ERASE);
								if (mEchoLevel > 1)
									std::cout << "     Removing a TIP POINT due to that criterion [" << in->Id() << "]" << std::endl;
								any_node_removed = true;
								boundary_nodes_removed++;
							}
						}
					}
					// else {
					// if (in->IsNot(FREE_SURFACE) && in->IsNot(RIGID) && freeSurfaceNeighNodes>0 && freeSurfaceNeighNodes<=dimension){
					// 	NodeWeakPtrVectorType& neighb_nodes = in->GetValue(NEIGHBOUR_NODES);
					// 	array_1d<double,3> sumOfCoordinates=in->Coordinates();
					// 	array_1d<double,3> sumOfCurrentVelocities=in->FastGetSolutionStepValue(VELOCITY,0);
					// 	array_1d<double,3> sumOfPreviousVelocities=in->FastGetSolutionStepValue(VELOCITY,1);
					// 	double sumOfPressures=in->FastGetSolutionStepValue(PRESSURE,0);
					// 	double counter=1.0;
					// 	// std::cout<<"I WAS GOING TO ERASE THIS NODE: "<<std::endl;
					// 	// std::cout<<"sumOfCoordinates: "<< sumOfCoordinates<<std::endl;
					// 	// std::cout<<"sumOfCurrentVelocities: "<< sumOfCurrentVelocities<<std::endl;
					// 	// std::cout<<"sumOfPressures: "<< sumOfPressures<<std::endl;
					// 	// std::cout<<"freeSurfaceNeighNodes "<<freeSurfaceNeighNodes<<"   rigidNeighNodes "<<rigidNeighNodes<<std::endl;
					// 	for (NodeWeakPtrVectorType::iterator nn = neighb_nodes.begin();nn != neighb_nodes.end(); nn++)
					// 	  {
					// 	    counter+=1.0;
					// 	    noalias(sumOfCoordinates)+=(*nn)->Coordinates();
					// 	    noalias(sumOfCurrentVelocities)+=(*nn)->FastGetSolutionStepValue(VELOCITY,0);
					// 	    noalias(sumOfPreviousVelocities)+=(*nn)->FastGetSolutionStepValue(VELOCITY,1);
					// 	    sumOfPressures+=(*nn)->FastGetSolutionStepValue(PRESSURE,0);
					// 	  }
					// 	in->X() =sumOfCoordinates[0]/counter;
					// 	in->Y() =sumOfCoordinates[1]/counter;
					// 	in->X0() =sumOfCoordinates[0]/counter;
					// 	in->Y0() =sumOfCoordinates[1]/counter;
					// 	in->FastGetSolutionStepValue(DISPLACEMENT_X,0)=0;
					// 	in->FastGetSolutionStepValue(DISPLACEMENT_Y,0)=0;
					// 	in->FastGetSolutionStepValue(DISPLACEMENT_X,1)=0;
					// 	in->FastGetSolutionStepValue(DISPLACEMENT_Y,1)=0;
					// 	in->FastGetSolutionStepValue(VELOCITY_X,0)=sumOfCurrentVelocities[0]/counter;
					// 	in->FastGetSolutionStepValue(VELOCITY_Y,0)=sumOfCurrentVelocities[1]/counter;
					// 	in->FastGetSolutionStepValue(VELOCITY_X,1)=sumOfPreviousVelocities[0]/counter;
					// 	in->FastGetSolutionStepValue(VELOCITY_Y,1)=sumOfPreviousVelocities[1]/counter;
					// 	in->FastGetSolutionStepValue(PRESSURE,0)=sumOfPressures/counter;
					// 	if(dimension==3){
					// 	  in->Z() =sumOfCoordinates[2]/counter;
					// 	  in->Z0() =sumOfCoordinates[2]/counter;
					// 	  in->FastGetSolutionStepValue(DISPLACEMENT_Z,0)=0;
					// 	  in->FastGetSolutionStepValue(DISPLACEMENT_Z,1)=0;
					// 	  in->FastGetSolutionStepValue(VELOCITY_Z,0)=sumOfCurrentVelocities[2]/counter;
					// 	  in->FastGetSolutionStepValue(VELOCITY_Z,1)=sumOfPreviousVelocities[2]/counter;
					// 	}
					// 	// std::cout<<"NOW ITS VARIABLES ARE: "<<std::endl;
					// 	// std::cout<<"Coordinates: "<< in->X()<<" " << in->Y()<<" " << in->Z()<<std::endl;
					// 	// std::cout<<"CurrentVelocities: "<<in->FastGetSolutionStepValue(VELOCITY_X,0)<<" "<<in->FastGetSolutionStepValue(VELOCITY_Y,0)<<" "<<in->FastGetSolutionStepValue(VELOCITY_Z,0)  <<std::endl;
					// 	// std::cout<<"Pressures: "<<in->FastGetSolutionStepValue(PRESSURE,0) <<std::endl;
					// }
					// }
				}
			}

			if (erased_nodes > 0)
			{
				if (mEchoLevel > 1)
					std::cout << "layer_nodes_removed " << erased_nodes << std::endl;
				any_node_removed = true;
			}
			//Build boundary after removing boundary nodes due distance criterion
			if (mEchoLevel > 1)
			{
				std::cout << "boundary_nodes_removed " << boundary_nodes_removed << std::endl;
				std::cout << "inside_nodes_removed " << inside_nodes_removed << std::endl;
			}
			return any_node_removed;

			KRATOS_CATCH(" ")
		}

		void EraseCriticalNodes2D(Element::GeometryType &eElement, unsigned int &erased_nodes, int &inside_nodes_removed)
		{

			KRATOS_TRY

			// std::cout<<"erased_nodes "<<erased_nodes<<std::endl;
			double safetyCoefficient2D = 0.5;
			double elementVolume = eElement.Area();

			unsigned int numNodes = eElement.size();
			// ////////  it erases nodes in very small elements /////////
			// double criticalVolume=0.1*mrRemesh.Refine->MeanVolume;
			// criticalVolume=0;
			// if(elementVolume<criticalVolume){
			//   for(unsigned int i=0; i<eElement.size(); i++)
			// 	{
			// 	  if(eElement[i].IsNot(RIGID) && eElement[i].IsNot(SOLID) && eElement[i].IsNot(TO_ERASE)){
			// 	    eElement[i].Set(TO_ERASE);
			// 	    if( mEchoLevel > 1 )
			// 	      std::cout<<"erase this layer node because it may be potentially dangerous and pass through the solid contour"<<std::endl;
			// 	    erased_nodes += 1;
			// 	    inside_nodes_removed++;
			// 	    break;
			// 	  }
			// 	}

			// }

			array_1d<double, 3> Edges(3, 0.0);
			array_1d<unsigned int, 3> FirstEdgeNode(3, 0);
			array_1d<unsigned int, 3> SecondEdgeNode(3, 0);
			double wallLength = 0;
			// array_1d<double,3> CoorDifference(3,0.0);

			// ////////  to compute the length of the wall edge /////////
			// noalias(CoorDifference) = eElement[1].Coordinates() - eElement[0].Coordinates();
			array_1d<double, 3> CoorDifference = eElement[1].Coordinates() - eElement[0].Coordinates();
			double SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1];
			Edges[0] = sqrt(SquaredLength);
			FirstEdgeNode[0] = 0;
			SecondEdgeNode[0] = 1;
			if (eElement[0].Is(RIGID) && eElement[1].Is(RIGID))
			{
				wallLength = Edges[0];
			}
			unsigned int counter = 0;
			for (unsigned int i = 2; i < eElement.size(); i++)
			{
				for (unsigned int j = 0; j < i; j++)
				{
					noalias(CoorDifference) = eElement[i].Coordinates() - eElement[j].Coordinates();
					SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1];
					counter += 1;
					Edges[counter] = sqrt(SquaredLength);
					FirstEdgeNode[counter] = j;
					SecondEdgeNode[counter] = i;
					if (eElement[i].Is(RIGID) && eElement[j].Is(RIGID) && Edges[counter] > wallLength)
					{
						wallLength = Edges[counter];
					}
				}
			}

			////////  to compare the triangle height to wall edge length /////////
			for (unsigned int i = 0; i < eElement.size(); i++)
			{
				if (eElement[i].IsNot(RIGID) && eElement[i].IsNot(TO_ERASE) && eElement[i].IsNot(SOLID) && eElement[i].IsNot(ISOLATED))
				{
					double height = elementVolume * 2.0 / wallLength;

					//////it is evident when a freesurface particle in touch with wall is erased --> reduce the safety coeff
					if (eElement[i].Is(FREE_SURFACE))
					{
						NodeWeakPtrVectorType &neighb_nodes = eElement[i].GetValue(NEIGHBOUR_NODES);
						unsigned int countRigid = 0;
						unsigned int countFreeSurface = 0;
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
						erased_nodes += 1;
						inside_nodes_removed++;
					}

					// // if the node is near to the wall but not too close, if possible, it is not erased but just moved in the middle of its largest edge (not shared with a wall node)
					// else if(height<safetyCoefficient2D*wallLength){
					//   bool eraseNode=true;
					//   eraseNode=CheckForMovingLayerNodes(Element[i],wallLength);

					//     if(eraseNode==true){
					//       // std::cout<<"I will erase this node because too close to neighbour nodes "<<std::endl;
					//       std::cout<<"(distances:  "<<height<<" vs "<<wallLength<<")"<<std::endl;
					//     Element[i].Set(TO_ERASE);
					//     erased_nodes += 1;
					//     inside_nodes_removed++;
					//   }
					// }
				}
			}

			bool longDamBreak = false; //to attivate in case of long dam breaks to avoid separeted elelements in the water front
			if (longDamBreak == true)
			{
				for (unsigned int i = 0; i < numNodes; i++)
				{
					if (eElement[i].Is(FREE_SURFACE) && eElement[i].IsNot(RIGID))
					{

						GlobalPointersVector<Element> &neighb_elems = eElement[i].GetValue(NEIGHBOUR_ELEMENTS);
						GlobalPointersVector<Node<3>> &neighb_nodes = eElement[i].GetValue(NEIGHBOUR_NODES);

						if (neighb_elems.size() < 2)
						{
							eElement[i].Set(TO_ERASE);
							std::cout << "erased an isolated element node" << std::endl;
							erased_nodes += 1;
							inside_nodes_removed++;
						}
						else
						{
							if (neighb_nodes.size() < 4)
							{
								for (unsigned int j = 0; j < neighb_nodes.size(); j++)
								{
									if (neighb_nodes[j].IsNot(FREE_SURFACE) && neighb_nodes[j].IsNot(RIGID))
									{
										break;
									}
									if (j == (neighb_nodes.size() - 1))
									{
										eElement[i].Set(TO_ERASE);
										std::cout << "_________________________          erased an isolated element node" << std::endl;
										erased_nodes += 1;
										inside_nodes_removed++;
									}
								}
							}
						}
					}
				}
			}

			// ////////  to compare the non-wall length to wall edge length /////////
			// for (unsigned int i = 0; i < 3; i++){
			//   if(((Element[FirstEdgeNode[i]].Is(RIGID) && Element[SecondEdgeNode[i]].IsNot(RIGID)) ||
			// 	  (Element[SecondEdgeNode[i]].Is(RIGID) && Element[FirstEdgeNode[i]].IsNot(RIGID))) &&
			// 	 Element[FirstEdgeNode[i]].IsNot(TO_ERASE) &&
			// 	 Element[SecondEdgeNode[i]].IsNot(TO_ERASE)&&
			// 	 Edges[i]<safetyCoefficient2D*wallLength){
			// 	if(Element[FirstEdgeNode[i]].IsNot(RIGID) && Element[FirstEdgeNode[i]].IsNot(SOLID) && Element[FirstEdgeNode[i]].IsNot(TO_ERASE)){
			// 	  Element[FirstEdgeNode[i]].Set(TO_ERASE);
			// 	  erased_nodes += 1;
			// 	  inside_nodes_removed++;
			// 	}else if(Element[SecondEdgeNode[i]].IsNot(RIGID) && Element[SecondEdgeNode[i]].IsNot(SOLID) && Element[SecondEdgeNode[i]].IsNot(TO_ERASE)){
			// 	  Element[SecondEdgeNode[i]].Set(TO_ERASE);
			// 	  erased_nodes += 1;
			// 	  inside_nodes_removed++;
			// 	}

			//   }

			// }
			KRATOS_CATCH("")
		}

		double SetMeshSizeInMeshRefinementArea(array_1d<double, 3> NodeCoordinates)
		{
			KRATOS_TRY

			array_1d<double, 3> minExternalPoint = mrRemesh.RefiningBoxMinExternalPoint;
			array_1d<double, 3> minInternalPoint = mrRemesh.RefiningBoxMinInternalPoint;
			array_1d<double, 3> maxExternalPoint = mrRemesh.RefiningBoxMaxExternalPoint;
			array_1d<double, 3> maxInternalPoint = mrRemesh.RefiningBoxMaxInternalPoint;
			array_1d<double, 3> RefiningBoxMinimumPoint = mrRemesh.RefiningBoxMinimumPoint;
			array_1d<double, 3> RefiningBoxMaximumPoint = mrRemesh.RefiningBoxMaximumPoint;

			double meshSize = mrRemesh.Refine->CriticalRadius;
			double distance = 2 * meshSize;
			double seperation = 0;
			double coefficient = 0;
			if (meshSize > mrRemesh.RefiningBoxMeshSize)
			{
				if (NodeCoordinates[0] > RefiningBoxMinimumPoint[0] && NodeCoordinates[0] < RefiningBoxMaximumPoint[0] &&
					NodeCoordinates[1] > RefiningBoxMinimumPoint[1] && NodeCoordinates[1] < RefiningBoxMaximumPoint[1])
				{
					meshSize = mrRemesh.RefiningBoxMeshSize; // in the internal domain the size is the one given by the user
				}
				else if ((NodeCoordinates[0] < RefiningBoxMinimumPoint[0] && NodeCoordinates[0] > (minExternalPoint[0] - distance) && NodeCoordinates[1] > minExternalPoint[1] && NodeCoordinates[1] < maxExternalPoint[1]))
				{
					seperation = NodeCoordinates[0] - RefiningBoxMinimumPoint[0];
					coefficient = fabs(seperation) / (distance + meshSize);
					meshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
				}
				else if ((NodeCoordinates[1] < RefiningBoxMinimumPoint[1] && NodeCoordinates[1] > (minExternalPoint[1] - distance) && NodeCoordinates[0] > minExternalPoint[0] && NodeCoordinates[0] < maxExternalPoint[0]))
				{
					seperation = NodeCoordinates[1] - RefiningBoxMinimumPoint[1];
					coefficient = fabs(seperation) / (distance + meshSize);
					meshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
				}
				else if ((NodeCoordinates[0] > RefiningBoxMaximumPoint[0] && NodeCoordinates[0] < (maxExternalPoint[0] + distance) && NodeCoordinates[1] > minExternalPoint[1] && NodeCoordinates[1] < maxExternalPoint[1]))
				{
					seperation = NodeCoordinates[0] - RefiningBoxMaximumPoint[0];
					coefficient = fabs(seperation) / (distance + meshSize);
					meshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
				}
				else if ((NodeCoordinates[1] > RefiningBoxMaximumPoint[1] && NodeCoordinates[1] < (maxExternalPoint[1] + distance) && NodeCoordinates[0] > minExternalPoint[0] && NodeCoordinates[0] < maxExternalPoint[0]))
				{
					seperation = NodeCoordinates[1] - RefiningBoxMaximumPoint[1];
					coefficient = fabs(seperation) / (distance + meshSize);
					meshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
				}
			}
			else
			{
				distance = 2.0 * mrRemesh.RefiningBoxMeshSize;

				if (NodeCoordinates[0] > (minInternalPoint[0] + distance) && NodeCoordinates[0] < (maxInternalPoint[0] - distance) &&
					NodeCoordinates[1] > (minInternalPoint[1] + distance) && NodeCoordinates[1] < (maxInternalPoint[1] - distance))
				{
					meshSize = mrRemesh.RefiningBoxMeshSize; // in the internal domain the size is the one given by the user
				}
				else if ((NodeCoordinates[0] > RefiningBoxMinimumPoint[0] && NodeCoordinates[0] < (minInternalPoint[0] + distance) && NodeCoordinates[1] > minExternalPoint[1] && NodeCoordinates[1] < maxExternalPoint[1]))
				{
					seperation = (minInternalPoint[0] + distance) - NodeCoordinates[0];
					coefficient = fabs(seperation) / (distance + mrRemesh.RefiningBoxMeshSize);
					meshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
				}
				else if ((NodeCoordinates[1] > RefiningBoxMinimumPoint[1] && NodeCoordinates[1] < (minInternalPoint[1] + distance) && NodeCoordinates[0] > minExternalPoint[0] && NodeCoordinates[0] < maxExternalPoint[0]))
				{
					seperation = (minInternalPoint[1] + distance) - NodeCoordinates[1];
					coefficient = fabs(seperation) / (distance + mrRemesh.RefiningBoxMeshSize);
					meshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
				}
				else if ((NodeCoordinates[0] < RefiningBoxMaximumPoint[0] && NodeCoordinates[0] > (maxInternalPoint[0] - distance) && NodeCoordinates[1] > minExternalPoint[1] && NodeCoordinates[1] < maxExternalPoint[1]))
				{
					seperation = (maxInternalPoint[0] - distance) - NodeCoordinates[0];
					coefficient = fabs(seperation) / (distance + mrRemesh.RefiningBoxMeshSize);
					meshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
				}
				else if ((NodeCoordinates[1] < RefiningBoxMaximumPoint[1] && NodeCoordinates[1] > (maxInternalPoint[1] - distance) && NodeCoordinates[0] > minExternalPoint[0] && NodeCoordinates[0] < maxExternalPoint[0]))
				{
					seperation = (maxInternalPoint[1] - distance) - NodeCoordinates[1];
					coefficient = fabs(seperation) / (distance + mrRemesh.RefiningBoxMeshSize);
					meshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
				}
			}

			return meshSize;

			KRATOS_CATCH("")
		}

		double SetMeshSizeInMeshRefinementVolume(array_1d<double, 3> NodeCoordinates)
		{
			KRATOS_TRY

			array_1d<double, 3> minExternalPoint = mrRemesh.RefiningBoxMinExternalPoint;
			array_1d<double, 3> minInternalPoint = mrRemesh.RefiningBoxMinInternalPoint;
			array_1d<double, 3> maxExternalPoint = mrRemesh.RefiningBoxMaxExternalPoint;
			array_1d<double, 3> maxInternalPoint = mrRemesh.RefiningBoxMaxInternalPoint;
			array_1d<double, 3> RefiningBoxMinimumPoint = mrRemesh.RefiningBoxMinimumPoint;
			array_1d<double, 3> RefiningBoxMaximumPoint = mrRemesh.RefiningBoxMaximumPoint;

			double meshSize = mrRemesh.Refine->CriticalRadius;
			double distance = 2.0 * meshSize;
			double seperation = 0;
			double coefficient = 0;
			if (meshSize > mrRemesh.RefiningBoxMeshSize)
			{
				if (NodeCoordinates[0] > RefiningBoxMinimumPoint[0] && NodeCoordinates[0] < RefiningBoxMaximumPoint[0] &&
					NodeCoordinates[1] > RefiningBoxMinimumPoint[1] && NodeCoordinates[1] < RefiningBoxMaximumPoint[1] &&
					NodeCoordinates[2] > RefiningBoxMinimumPoint[2] && NodeCoordinates[2] < RefiningBoxMaximumPoint[2])
				{
					meshSize = mrRemesh.RefiningBoxMeshSize; // in the internal domain the size is the one given by the user
				}
				else if ((NodeCoordinates[0] < RefiningBoxMinimumPoint[0] && NodeCoordinates[0] > (minExternalPoint[0] - distance) && NodeCoordinates[1] > minExternalPoint[1] && NodeCoordinates[1] < maxExternalPoint[1] && NodeCoordinates[2] > minExternalPoint[2] && NodeCoordinates[2] < maxExternalPoint[2]))
				{
					seperation = NodeCoordinates[0] - RefiningBoxMinimumPoint[0];
					coefficient = fabs(seperation) / (distance + meshSize);
					meshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
				}
				else if ((NodeCoordinates[1] < RefiningBoxMinimumPoint[1] && NodeCoordinates[1] > (minExternalPoint[1] - distance) && NodeCoordinates[0] > minExternalPoint[0] && NodeCoordinates[0] < maxExternalPoint[0] && NodeCoordinates[2] > minExternalPoint[2] && NodeCoordinates[2] < maxExternalPoint[2]))
				{
					seperation = NodeCoordinates[1] - RefiningBoxMinimumPoint[1];
					coefficient = fabs(seperation) / (distance + meshSize);
					meshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
				}
				else if ((NodeCoordinates[2] < RefiningBoxMinimumPoint[2] && NodeCoordinates[2] > (minExternalPoint[2] - distance) && NodeCoordinates[0] > minExternalPoint[0] && NodeCoordinates[0] < maxExternalPoint[0] && NodeCoordinates[1] > minExternalPoint[1] && NodeCoordinates[1] < maxExternalPoint[1]))
				{
					seperation = NodeCoordinates[2] - RefiningBoxMinimumPoint[2];
					coefficient = fabs(seperation) / (distance + meshSize);
					meshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
				}
				else if ((NodeCoordinates[0] > RefiningBoxMaximumPoint[0] && NodeCoordinates[0] < (maxExternalPoint[0] + distance) && NodeCoordinates[1] > minExternalPoint[1] && NodeCoordinates[1] < maxExternalPoint[1] && NodeCoordinates[2] > minExternalPoint[2] && NodeCoordinates[2] < maxExternalPoint[2]))
				{
					seperation = NodeCoordinates[0] - RefiningBoxMaximumPoint[0];
					coefficient = fabs(seperation) / (distance + meshSize);
					meshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
				}
				else if ((NodeCoordinates[1] > RefiningBoxMaximumPoint[1] && NodeCoordinates[1] < (maxExternalPoint[1] + distance) && NodeCoordinates[0] > minExternalPoint[0] && NodeCoordinates[0] < maxExternalPoint[0] && NodeCoordinates[2] > minExternalPoint[2] && NodeCoordinates[2] < maxExternalPoint[2]))
				{
					seperation = NodeCoordinates[1] - RefiningBoxMaximumPoint[1];
					coefficient = fabs(seperation) / (distance + meshSize);
					meshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
				}
				else if ((NodeCoordinates[2] > RefiningBoxMaximumPoint[2] && NodeCoordinates[2] < (maxExternalPoint[2] + distance) && NodeCoordinates[0] > minExternalPoint[0] && NodeCoordinates[0] < maxExternalPoint[0] && NodeCoordinates[1] > minExternalPoint[1] && NodeCoordinates[1] < maxExternalPoint[1]))
				{
					seperation = NodeCoordinates[2] - RefiningBoxMaximumPoint[2];
					coefficient = fabs(seperation) / (distance + meshSize);
					meshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
				}
			}
			else
			{

				distance = 2.0 * mrRemesh.RefiningBoxMeshSize;

				if (NodeCoordinates[0] > (minInternalPoint[0] + distance) && NodeCoordinates[0] < (maxInternalPoint[0] - distance) &&
					NodeCoordinates[1] > (minInternalPoint[1] + distance) && NodeCoordinates[1] < (maxInternalPoint[1] - distance) &&
					NodeCoordinates[2] > (minInternalPoint[2] + distance) && NodeCoordinates[2] < (maxInternalPoint[2] - distance))
				{
					meshSize = mrRemesh.RefiningBoxMeshSize; // in the internal domain the size is the one given by the user
				}
				else if ((NodeCoordinates[0] > RefiningBoxMinimumPoint[0] && NodeCoordinates[0] < (minInternalPoint[0] + distance) && NodeCoordinates[1] > minExternalPoint[1] && NodeCoordinates[1] < maxExternalPoint[1] && NodeCoordinates[2] > minExternalPoint[2] && NodeCoordinates[2] < maxExternalPoint[2]))
				{
					seperation = (minInternalPoint[0] + distance) - NodeCoordinates[0];
					coefficient = fabs(seperation) / (distance + mrRemesh.RefiningBoxMeshSize);
					meshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
				}
				else if ((NodeCoordinates[1] > RefiningBoxMinimumPoint[1] && NodeCoordinates[1] < (minInternalPoint[1] + distance) && NodeCoordinates[0] > minExternalPoint[0] && NodeCoordinates[0] < maxExternalPoint[0] && NodeCoordinates[2] > minExternalPoint[2] && NodeCoordinates[2] < maxExternalPoint[2]))
				{
					seperation = (minInternalPoint[1] + distance) - NodeCoordinates[1];
					coefficient = fabs(seperation) / (distance + mrRemesh.RefiningBoxMeshSize);
					meshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
				}
				else if ((NodeCoordinates[2] > RefiningBoxMinimumPoint[2] && NodeCoordinates[2] < (minInternalPoint[2] + distance) && NodeCoordinates[0] > minExternalPoint[0] && NodeCoordinates[0] < maxExternalPoint[0] && NodeCoordinates[1] > minExternalPoint[1] && NodeCoordinates[1] < maxExternalPoint[1]))
				{
					seperation = (minInternalPoint[2] + distance) - NodeCoordinates[2];
					coefficient = fabs(seperation) / (distance + mrRemesh.RefiningBoxMeshSize);
					meshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
				}
				else if ((NodeCoordinates[0] < RefiningBoxMaximumPoint[0] && NodeCoordinates[0] > (maxInternalPoint[0] - distance) && NodeCoordinates[1] > minExternalPoint[1] && NodeCoordinates[1] < maxExternalPoint[1] && NodeCoordinates[2] > minExternalPoint[2] && NodeCoordinates[2] < maxExternalPoint[2]))
				{
					seperation = (maxInternalPoint[0] - distance) - NodeCoordinates[0];
					coefficient = fabs(seperation) / (distance + mrRemesh.RefiningBoxMeshSize);
					meshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
				}
				else if ((NodeCoordinates[1] < RefiningBoxMaximumPoint[1] && NodeCoordinates[1] > (maxInternalPoint[1] - distance) && NodeCoordinates[0] > minExternalPoint[0] && NodeCoordinates[0] < maxExternalPoint[0] && NodeCoordinates[2] > minExternalPoint[2] && NodeCoordinates[2] < maxExternalPoint[2]))
				{
					seperation = (maxInternalPoint[1] - distance) - NodeCoordinates[1];
					coefficient = fabs(seperation) / (distance + mrRemesh.RefiningBoxMeshSize);
					meshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
				}
				else if ((NodeCoordinates[2] < RefiningBoxMaximumPoint[2] && NodeCoordinates[2] > (maxInternalPoint[2] - distance) && NodeCoordinates[0] > minExternalPoint[0] && NodeCoordinates[0] < maxExternalPoint[0] && NodeCoordinates[1] > minExternalPoint[1] && NodeCoordinates[1] < maxExternalPoint[1]))
				{
					seperation = (maxInternalPoint[2] - distance) - NodeCoordinates[2];
					coefficient = fabs(seperation) / (distance + mrRemesh.RefiningBoxMeshSize);
					meshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
				}
			}
			return meshSize;

			KRATOS_CATCH("")
		}

		void EraseCriticalNodes3D(Element::GeometryType &eElement, unsigned int &erased_nodes, int &inside_nodes_removed, unsigned int rigidNodes)
		{

			KRATOS_TRY
			double safetyCoefficient3D = 0.6;
			// double safetyCoefficient3D=0.7;

			unsigned int freeSurfaceNodes = 0;
			unsigned int numNodes = eElement.size();
			double elementVolume = eElement.Volume();
			double criticalVolume = 0.1 * mrRemesh.Refine->MeanVolume;

			std::vector<array_1d<double, 3>> rigidNodesCoordinates;
			std::vector<array_1d<double, 3>> rigidNodesNormals;
			array_1d<double, 3> notRigidNodeCoordinates(3, 0.0);
			unsigned int notRigidNodeId = 0;
			rigidNodesCoordinates.resize(3);
			rigidNodesNormals.resize(3);
			double baricenterX = 0.25 * (eElement[0].X() + eElement[1].X() + eElement[2].X() + eElement[3].X());
			double baricenterY = 0.25 * (eElement[0].Y() + eElement[1].Y() + eElement[2].Y() + eElement[3].Y());
			double baricenterZ = 0.25 * (eElement[0].Z() + eElement[1].Z() + eElement[2].Z() + eElement[3].Z());

			if (mrRemesh.UseRefiningBox == true)
			{
				array_1d<double, 3> RefiningBoxMinimumPoint = mrRemesh.RefiningBoxMinimumPoint;
				array_1d<double, 3> RefiningBoxMaximumPoint = mrRemesh.RefiningBoxMaximumPoint;
				if (baricenterX > RefiningBoxMinimumPoint[0] && baricenterX < RefiningBoxMaximumPoint[0] &&
					baricenterY > RefiningBoxMinimumPoint[1] && baricenterY < RefiningBoxMaximumPoint[1] &&
					baricenterZ > RefiningBoxMinimumPoint[2] && baricenterZ < RefiningBoxMaximumPoint[2])
				{
					criticalVolume = 0.01 * (pow(mrRemesh.RefiningBoxMeshSize, 3) / (6.0 * sqrt(2))); //mean Volume of a regular tetrahedral per node with 0.01 of penalization
				}
				else
				{
					criticalVolume = 0.01 * (pow(mrRemesh.Refine->CriticalRadius, 3) / (6.0 * sqrt(2)));
				}
				double meanSize = 0.5 * (mrRemesh.RefiningBoxMeshSize + mrRemesh.Refine->CriticalRadius);
				if ((baricenterX > (RefiningBoxMinimumPoint[0] - meanSize) && baricenterX < (RefiningBoxMaximumPoint[0] + meanSize)) ||
					(baricenterY > (RefiningBoxMinimumPoint[1] - meanSize) && baricenterY < (RefiningBoxMaximumPoint[1] + meanSize)) ||
					(baricenterZ > (RefiningBoxMinimumPoint[2] - meanSize) && baricenterZ < (RefiningBoxMaximumPoint[2] + meanSize))) //transition zone
				{
					safetyCoefficient3D *= 0.8;
				}
			}

			unsigned int rigidNode = 0;
			for (unsigned int i = 0; i < numNodes; i++)
			{
				if (eElement[i].Is(FREE_SURFACE))
				{
					freeSurfaceNodes++;
				}
				if (rigidNodes == 3)
				{
					if (eElement[i].Is(RIGID))
					{
						rigidNodesCoordinates[rigidNode] = eElement[i].Coordinates();
						rigidNodesNormals[rigidNode] = eElement[i].FastGetSolutionStepValue(NORMAL);
						rigidNode++;
					}
					else
					{
						notRigidNodeCoordinates = eElement[i].Coordinates();
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
						erased_nodes += 1;
						inside_nodes_removed++;
					}
				}
			}

			if (rigidNode == 3 && eElement[notRigidNodeId].IsNot(TO_ERASE))
			{
				double a1 = 0; //slope x,y,z for the plane composed by rigid nodes only
				double b1 = 0;
				double c1 = 0;
				a1 = (rigidNodesCoordinates[1][1] - rigidNodesCoordinates[0][1]) * (rigidNodesCoordinates[2][2] - rigidNodesCoordinates[0][2]) - (rigidNodesCoordinates[2][1] - rigidNodesCoordinates[0][1]) * (rigidNodesCoordinates[1][2] - rigidNodesCoordinates[0][2]);
				b1 = (rigidNodesCoordinates[1][2] - rigidNodesCoordinates[0][2]) * (rigidNodesCoordinates[2][0] - rigidNodesCoordinates[0][0]) - (rigidNodesCoordinates[2][2] - rigidNodesCoordinates[0][2]) * (rigidNodesCoordinates[1][0] - rigidNodesCoordinates[0][0]);
				c1 = (rigidNodesCoordinates[1][0] - rigidNodesCoordinates[0][0]) * (rigidNodesCoordinates[2][1] - rigidNodesCoordinates[0][1]) - (rigidNodesCoordinates[2][0] - rigidNodesCoordinates[0][0]) * (rigidNodesCoordinates[1][1] - rigidNodesCoordinates[0][1]);
				double a2 = 0; //slope x,y,z between 2 rigid nodes and a not rigid node
				double b2 = 0;
				double c2 = 0;
				a2 = (rigidNodesCoordinates[1][1] - rigidNodesCoordinates[0][1]) * (notRigidNodeCoordinates[2] - rigidNodesCoordinates[0][2]) - (notRigidNodeCoordinates[1] - rigidNodesCoordinates[0][1]) * (rigidNodesCoordinates[1][2] - rigidNodesCoordinates[0][2]);
				b2 = (rigidNodesCoordinates[1][2] - rigidNodesCoordinates[0][2]) * (notRigidNodeCoordinates[0] - rigidNodesCoordinates[0][0]) - (notRigidNodeCoordinates[2] - rigidNodesCoordinates[0][2]) * (rigidNodesCoordinates[1][0] - rigidNodesCoordinates[0][0]);
				c2 = (rigidNodesCoordinates[1][0] - rigidNodesCoordinates[0][0]) * (notRigidNodeCoordinates[1] - rigidNodesCoordinates[0][1]) - (notRigidNodeCoordinates[0] - rigidNodesCoordinates[0][0]) * (rigidNodesCoordinates[1][1] - rigidNodesCoordinates[0][1]);
				double a3 = 0; //slope x,y,z between 2 rigid nodes and a not rigid node
				double b3 = 0;
				double c3 = 0;
				a3 = (rigidNodesCoordinates[1][1] - rigidNodesCoordinates[2][1]) * (notRigidNodeCoordinates[2] - rigidNodesCoordinates[2][2]) - (notRigidNodeCoordinates[1] - rigidNodesCoordinates[2][1]) * (rigidNodesCoordinates[1][2] - rigidNodesCoordinates[2][2]);
				b3 = (rigidNodesCoordinates[1][2] - rigidNodesCoordinates[2][2]) * (notRigidNodeCoordinates[0] - rigidNodesCoordinates[2][0]) - (notRigidNodeCoordinates[2] - rigidNodesCoordinates[2][2]) * (rigidNodesCoordinates[1][0] - rigidNodesCoordinates[2][0]);
				c3 = (rigidNodesCoordinates[1][0] - rigidNodesCoordinates[2][0]) * (notRigidNodeCoordinates[1] - rigidNodesCoordinates[2][1]) - (notRigidNodeCoordinates[0] - rigidNodesCoordinates[2][0]) * (rigidNodesCoordinates[1][1] - rigidNodesCoordinates[2][1]);
				double a4 = 0; //slope x,y,z between 2 rigid nodes and a not rigid node
				double b4 = 0;
				double c4 = 0;
				a4 = (rigidNodesCoordinates[0][1] - rigidNodesCoordinates[2][1]) * (notRigidNodeCoordinates[2] - rigidNodesCoordinates[2][2]) - (notRigidNodeCoordinates[1] - rigidNodesCoordinates[2][1]) * (rigidNodesCoordinates[0][2] - rigidNodesCoordinates[2][2]);
				b4 = (rigidNodesCoordinates[0][2] - rigidNodesCoordinates[2][2]) * (notRigidNodeCoordinates[0] - rigidNodesCoordinates[2][0]) - (notRigidNodeCoordinates[2] - rigidNodesCoordinates[2][2]) * (rigidNodesCoordinates[0][0] - rigidNodesCoordinates[2][0]);
				c4 = (rigidNodesCoordinates[0][0] - rigidNodesCoordinates[2][0]) * (notRigidNodeCoordinates[1] - rigidNodesCoordinates[2][1]) - (notRigidNodeCoordinates[0] - rigidNodesCoordinates[2][0]) * (rigidNodesCoordinates[0][1] - rigidNodesCoordinates[2][1]);

				//angle between the plane composed by rigid nodes only and the other plans. If the angle is small, the particle can pass through the wall
				double cosAngle12 = (a1 * a2 + b1 * b2 + c1 * c2) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a2, 2) + pow(b2, 2) + pow(c2, 2)));
				double cosAngle13 = (a1 * a3 + b1 * b3 + c1 * c3) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a3, 2) + pow(b3, 2) + pow(c3, 2)));
				double cosAngle14 = (a1 * a4 + b1 * b4 + c1 * c4) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a4, 2) + pow(b4, 2) + pow(c4, 2)));

				//angle between the normals of the rigid nodes. I want to avoid rigid elements at the corner
				double cosAngleBetweenNormals01 = (rigidNodesNormals[0][0] * rigidNodesNormals[1][0] + rigidNodesNormals[0][1] * rigidNodesNormals[1][1]) /
												  (sqrt(pow(rigidNodesNormals[0][0], 2) + pow(rigidNodesNormals[0][1], 2)) *
												   sqrt(pow(rigidNodesNormals[1][0], 2) + pow(rigidNodesNormals[1][1], 2)));
				double cosAngleBetweenNormals02 = (rigidNodesNormals[0][0] * rigidNodesNormals[2][0] + rigidNodesNormals[0][1] * rigidNodesNormals[2][1]) /
												  (sqrt(pow(rigidNodesNormals[0][0], 2) + pow(rigidNodesNormals[0][1], 2)) *
												   sqrt(pow(rigidNodesNormals[2][0], 2) + pow(rigidNodesNormals[2][1], 2)));
				double cosAngleBetweenNormals12 = (rigidNodesNormals[1][0] * rigidNodesNormals[2][0] + rigidNodesNormals[1][1] * rigidNodesNormals[2][1]) /
												  (sqrt(pow(rigidNodesNormals[1][0], 2) + pow(rigidNodesNormals[1][1], 2)) *
												   sqrt(pow(rigidNodesNormals[2][0], 2) + pow(rigidNodesNormals[2][1], 2)));

				if ((fabs(cosAngle12) > 0.995 || fabs(cosAngle13) > 0.995 || fabs(cosAngle14) > 0.995) && (cosAngleBetweenNormals01 > 0.99 && cosAngleBetweenNormals02 > 0.99 && cosAngleBetweenNormals12 > 0.99))
				{
					eElement[notRigidNodeId].Set(TO_ERASE);
					//std::cout << eElement[notRigidNodeId].Id() << " nodeId is erased because it may pass through the solid contour. Coordinates are: " << notRigidNodeCoordinates << std::endl;
					erased_nodes += 1;
					inside_nodes_removed++;
				}
			}

			bool longDamBreak = false; //to attivate in case of long dam breaks to avoid separated elements in the water front
			if (longDamBreak == true && freeSurfaceNodes > 2 && rigidNodes > 1)
			{
				for (unsigned int i = 0; i < numNodes; i++)
				{
					if (eElement[i].Is(FREE_SURFACE) && eElement[i].IsNot(RIGID))
					{

						GlobalPointersVector<Element> &neighb_elems = eElement[i].GetValue(NEIGHBOUR_ELEMENTS);
						GlobalPointersVector<Node<3>> &neighb_nodes = eElement[i].GetValue(NEIGHBOUR_NODES);

						if (neighb_elems.size() < 2)
						{
							eElement[i].Set(TO_ERASE);
							// std::cout<<"erased an isolated element node"<<std::endl;
							erased_nodes += 1;
							inside_nodes_removed++;
						}
						else
						{
							if (neighb_nodes.size() < 10)
							{
								unsigned int freeSurfaceNodesNeigh = 0;
								for (unsigned int j = 0; j < neighb_nodes.size(); j++)
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
										// std::cout<<"_________________________          erased an isolated element node"<<std::endl;
										erased_nodes += 1;
										inside_nodes_removed++;
									}
								}
							}
							else
							{
								for (unsigned int j = 0; j < neighb_nodes.size(); j++)
								{
									if (neighb_nodes[j].IsNot(RIGID))
									{
										break;
									}
									if (j == (neighb_nodes.size() - 1))
									{
										eElement[i].Set(TO_ERASE);
										// std::cout<<"_________________________          erased an isolated wall element node"<<std::endl;
										erased_nodes += 1;
										inside_nodes_removed++;
									}
								}
							}
						}
					}
				}
			}

			array_1d<double, 6> Edges(6, 0.0);
			array_1d<unsigned int, 6> FirstEdgeNode(6, 0);
			array_1d<unsigned int, 6> SecondEdgeNode(6, 0);
			double wallLength = 0;
			double minimumLength = 0;
			// array_1d<double,3> CoorDifference(3,0.0);

			// ////////  to compute the length of the wall edge /////////
			// CoorDifference = Element[1].Coordinates() - Element[0].Coordinates();
			array_1d<double, 3> CoorDifference = eElement[1].Coordinates() - eElement[0].Coordinates();

			double SquaredLength = CoorDifference[0] * CoorDifference[0] +
								   CoorDifference[1] * CoorDifference[1] +
								   CoorDifference[2] * CoorDifference[2];
			Edges[0] = sqrt(SquaredLength);
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
			}
			unsigned int counter = 0;
			for (unsigned int i = 2; i < eElement.size(); i++)
			{
				for (unsigned int j = 0; j < i; j++)
				{
					noalias(CoorDifference) = eElement[i].Coordinates() - eElement[j].Coordinates();
					// CoorDifference = Element[i].Coordinates() - Element[j].Coordinates();
					SquaredLength = CoorDifference[0] * CoorDifference[0] +
									CoorDifference[1] * CoorDifference[1] +
									CoorDifference[2] * CoorDifference[2];
					counter += 1;
					Edges[counter] = sqrt(SquaredLength);
					FirstEdgeNode[counter] = j;
					SecondEdgeNode[counter] = i;
					if (eElement[i].Is(RIGID) && eElement[j].Is(RIGID) && wallLength == 0)
					{
						wallLength = Edges[counter];
					}
					if (((eElement[i].Is(RIGID) && eElement[j].IsNot(RIGID)) ||
						 (eElement[j].Is(RIGID) && eElement[i].IsNot(RIGID))) &&
						(Edges[counter] < minimumLength || minimumLength == 0))
					{
						minimumLength = Edges[counter];
					}
				}
			}

			////////  to avoid the elimination of isolated free-surface-rigid elements /////////
			for (unsigned int i = 0; i < eElement.size(); i++)
			{
				if (eElement[i].IsNot(RIGID) && eElement[i].IsNot(TO_ERASE) && eElement[i].IsNot(SOLID) && eElement[i].IsNot(ISOLATED))
				{
					//////it is evident when a freesurface particle in touch with wall is erased --> reduce the safety coeff
					if (eElement[i].Is(FREE_SURFACE))
					{
						NodeWeakPtrVectorType &neighb_nodes = eElement[i].GetValue(NEIGHBOUR_NODES);
						unsigned int countRigid = 0;
						unsigned int countFreeSurface = 0;
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

			// ////////// ////////// ////////// ////////// ////////// ////////
			// if(minimumLength<(0.5*safetyCoefficient3D*wallLength)){
			//   std::cout<<"(1. minimumLength:  "<<minimumLength<<" vs "<<wallLength<<")"<<std::endl;

			//   for (unsigned int i = 0; i < Element.size(); i++){
			// 	if(Element[i].IsNot(RIGID) && Element[i].IsNot(TO_ERASE) && Element[i].IsNot(SOLID) && Element[i].IsNot(ISOLATED)){

			// 	  Element[i].Set(TO_ERASE);
			// 	  inside_nodes_removed++;
			// 	  erased_nodes += 1;
			// 	}
			//   }
			// }
			// else if(minimumLength<safetyCoefficient3D*wallLength){

			//   std::cout<<"(2. minimumLength:  "<<minimumLength<<" vs "<<wallLength<<")"<<std::endl;

			//   for (unsigned int i = 0; i < Element.size(); i++){
			// 	if(Element[i].IsNot(RIGID) && Element[i].IsNot(TO_ERASE) && Element[i].IsNot(SOLID) && Element[i].IsNot(ISOLATED)){
			// 	  bool eraseNode=true;
			// 	  eraseNode=CheckForMovingLayerNodes(Element[i],wallLength);

			// 	  if(eraseNode==true){
			// 	    std::cout<<"I will erase this node because too close to neighbour nodes "<<std::endl;
			// 	    Element[i].Set(TO_ERASE);
			// 	    erased_nodes += 1;
			// 	    inside_nodes_removed++;
			// 	  }
			// 	}
			//   }
			// }
			// ////////// ////////// ////////// ////////// ////////// ////////

			// ////////  to compare the non-wall length to wall edge length /////////
			for (unsigned int i = 0; i < Edges.size(); i++)
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
						erased_nodes += 1;
					}
					else if (eElement[SecondEdgeNode[i]].IsNot(RIGID) && eElement[SecondEdgeNode[i]].IsNot(SOLID) && eElement[SecondEdgeNode[i]].IsNot(TO_ERASE) && eElement[SecondEdgeNode[i]].IsNot(ISOLATED))
					{
						eElement[SecondEdgeNode[i]].Set(TO_ERASE);
						inside_nodes_removed++;
						erased_nodes += 1;
					}
				}
			}

			KRATOS_CATCH("")
		}

		bool CheckForMovingLayerNodes(Node<3> &CheckedNode, const double wallLength)
		{
			KRATOS_TRY
			const unsigned int dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

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
			//I have looked for the biggest edge for moving there the layer node
			double maxNeighDistance = sqrt(maxSquaredDistance);
			if (maxNeighDistance > wallLength && wallLength > 0)
			{
				for (NodeWeakPtrVectorType::iterator nn = neighb_nodes.begin(); nn != neighb_nodes.end(); nn++)
				{
					if (nn == j)
					{

						unsigned int idMaster = CheckedNode.GetId();
						unsigned int idSlave = (j)->GetId();
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

		void InterpolateFromTwoNodes(unsigned int idMaster, unsigned int idSlave1, unsigned int idSlave2)
		{

			KRATOS_TRY

			Node<3>::Pointer MasterNode = mrModelPart.pGetNode(idMaster);
			Node<3>::Pointer SlaveNode1 = mrModelPart.pGetNode(idSlave1);
			Node<3>::Pointer SlaveNode2 = mrModelPart.pGetNode(idSlave2);

			VariablesList &rVariablesList = mrModelPart.GetNodalSolutionStepVariablesList();

			unsigned int buffer_size = MasterNode->GetBufferSize();

			for (VariablesList::const_iterator i_variable = rVariablesList.begin(); i_variable != rVariablesList.end(); i_variable++)
			{
				//std::cout<<" name "<<i_variable->Name()<<std::endl;
				//std::cout<<" type "<<typeid(*i_variable).name()<<std::endl;
				std::string variable_name = i_variable->Name();
				if (KratosComponents<Variable<double>>::Has(variable_name))
				{
					//std::cout<<"double"<<std::endl;
					const Variable<double> &variable = KratosComponents<Variable<double>>::Get(variable_name);
					for (unsigned int step = 0; step < buffer_size; step++)
					{
						//getting the data of the solution step
						double &node_data = MasterNode->FastGetSolutionStepValue(variable, step);

						double node0_data = SlaveNode1->FastGetSolutionStepValue(variable, step);
						double node1_data = SlaveNode2->FastGetSolutionStepValue(variable, step);

						node_data = (0.5 * node0_data + 0.5 * node1_data);
					}
				}
				else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name))
				{
					//std::cout<<"array1d"<<std::endl;
					const Variable<array_1d<double, 3>> &variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(variable_name);
					for (unsigned int step = 0; step < buffer_size; step++)
					{
						//getting the data of the solution step
						array_1d<double, 3> &node_data = MasterNode->FastGetSolutionStepValue(variable, step);

						const array_1d<double, 3> &node0_data = SlaveNode1->FastGetSolutionStepValue(variable, step);
						const array_1d<double, 3> &node1_data = SlaveNode2->FastGetSolutionStepValue(variable, step);
						noalias(node_data) = (0.5 * node0_data + 0.5 * node1_data);
						// node_data = (0.5*node0_data + 0.5*node1_data);
					}
				}
				else if (KratosComponents<Variable<int>>::Has(variable_name))
				{
					//std::cout<<"int"<<std::endl;
					//NO INTERPOLATION
				}
				else if (KratosComponents<Variable<bool>>::Has(variable_name))
				{
					//std::cout<<"bool"<<std::endl;
					//NO INTERPOLATION
				}
				else if (KratosComponents<Variable<Matrix>>::Has(variable_name))
				{
					//std::cout<<"Matrix"<<std::endl;
					const Variable<Matrix> &variable = KratosComponents<Variable<Matrix>>::Get(variable_name);
					for (unsigned int step = 0; step < buffer_size; step++)
					{
						//getting the data of the solution step
						Matrix &node_data = MasterNode->FastGetSolutionStepValue(variable, step);

						Matrix &node0_data = SlaveNode1->FastGetSolutionStepValue(variable, step);
						Matrix &node1_data = SlaveNode2->FastGetSolutionStepValue(variable, step);

						if (node_data.size1() > 0 && node_data.size2())
						{
							if (node_data.size1() == node0_data.size1() && node_data.size2() == node0_data.size2() &&
								node_data.size1() == node1_data.size1() && node_data.size2() == node1_data.size2())
							{
								noalias(node_data) = (0.5 * node0_data + 0.5 * node1_data);
								// node_data = (0.5*node0_data + 0.5*node1_data);
							}
						}
					}
				}
				else if (KratosComponents<Variable<Vector>>::Has(variable_name))
				{
					//std::cout<<"Vector"<<std::endl;
					const Variable<Vector> &variable = KratosComponents<Variable<Vector>>::Get(variable_name);
					for (unsigned int step = 0; step < buffer_size; step++)
					{
						//getting the data of the solution step
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

		/// Assignment operator.
		RemoveMeshNodesForFluidsProcess &operator=(RemoveMeshNodesForFluidsProcess const &rOther);

		/// this function is a private function

		/// Copy constructor.
		//Process(Process const& rOther);

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
