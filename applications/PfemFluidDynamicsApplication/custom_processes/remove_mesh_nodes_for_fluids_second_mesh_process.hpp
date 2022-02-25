//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:                   AFranci $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               November 2021 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined(KRATOS_REMOVE_MESH_NODES_FOR_FLUIDS_SECOND_MESH_PROCESS_H_INCLUDED)
#define KRATOS_REMOVE_MESH_NODES_FOR_FLUIDS_SECOND_MESH_PROCESS_H_INCLUDED

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

	class RemoveMeshNodesForFluidsSecondMeshProcess
		: public MesherProcess
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of Process
		KRATOS_CLASS_POINTER_DEFINITION(RemoveMeshNodesForFluidsSecondMeshProcess);

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
		RemoveMeshNodesForFluidsSecondMeshProcess(ModelPart &rModelPart,
												  MesherUtilities::MeshingParameters &rRemeshingParameters,
												  int EchoLevel)
			: mrModelPart(rModelPart),
			  mrRemesh(rRemeshingParameters)
		{
			KRATOS_INFO("RemoveMeshNodesForFluidsSecondMeshProcess") << " activated " << std::endl;

			mEchoLevel = EchoLevel;
		}

		/// Destructor.
		virtual ~RemoveMeshNodesForFluidsSecondMeshProcess() {}

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

			bool any_node_removed = false;
			bool any_node_removed_on_distance = false;

			int critical_nodes_removed = 0;
			any_node_removed_on_distance = RemoveCriticalNodes(critical_nodes_removed);

			if (any_node_removed_on_distance)
				any_node_removed = true;

			if (any_node_removed || mrRemesh.UseBoundingBox == true)
				this->CleanRemovedNodes(mrModelPart);

			// number of removed nodes:
			mrRemesh.Info->RemovedNodes += critical_nodes_removed;
			int distance_remove = critical_nodes_removed;

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
			return "RemoveMeshNodesForFluidsSecondMeshProcess";
		}

		/// Print information about this object.
		void PrintInfo(std::ostream &rOStream) const override
		{
			rOStream << "RemoveMeshNodesForFluidsSecondMeshProcess";
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


		bool RemoveCriticalNodes(int &inside_nodes_removed)
		{
			KRATOS_TRY

			const unsigned int dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
			bool any_node_removed = false;
			unsigned int erased_nodes = 0;

			if (dimension == 3)
			{
				for (ModelPart::ElementsContainerType::const_iterator ie = mrModelPart.ElementsBegin(); ie != mrModelPart.ElementsEnd(); ie++)
				{
					EraseSliverNodes(ie->GetGeometry(), erased_nodes, inside_nodes_removed);
				}
			}

			if (erased_nodes > 0)
			{
				any_node_removed = true;
				if (mEchoLevel > 1)
					std::cout << "layer_nodes_removed " << erased_nodes << std::endl;
			}
			//Build boundary after removing boundary nodes due distance criterion
			if (mEchoLevel > 1)
			{
				std::cout << "inside_nodes_removed " << inside_nodes_removed << std::endl;
			}
			return any_node_removed;

			KRATOS_CATCH(" ")
		}

		void EraseSliverNodes(Element::GeometryType &eElement, unsigned int &erased_nodes, int &inside_nodes_removed)
		{
			double a1 = 0; //slope x for plane on the first triangular face of the tetrahedra (nodes A,B,C)
			double b1 = 0; //slope y for plane on the first triangular face of the tetrahedra (nodes A,B,C)
			double c1 = 0; //slope z for plane on the first triangular face of the tetrahedra (nodes A,B,C)
			a1 = (eElement[1].Y() - eElement[0].Y()) * (eElement[2].Z() - eElement[0].Z()) - (eElement[2].Y() - eElement[0].Y()) * (eElement[1].Z() - eElement[0].Z());
			b1 = (eElement[1].Z() - eElement[0].Z()) * (eElement[2].X() - eElement[0].X()) - (eElement[2].Z() - eElement[0].Z()) * (eElement[1].X() - eElement[0].X());
			c1 = (eElement[1].X() - eElement[0].X()) * (eElement[2].Y() - eElement[0].Y()) - (eElement[2].X() - eElement[0].X()) * (eElement[1].Y() - eElement[0].Y());
			double a2 = 0; //slope x for plane on the second triangular face of the tetrahedra (nodes A,B,D)
			double b2 = 0; //slope y for plane on the second triangular face of the tetrahedra (nodes A,B,D)
			double c2 = 0; //slope z for plane on the second triangular face of the tetrahedra (nodes A,B,D)
			a2 = (eElement[1].Y() - eElement[0].Y()) * (eElement[3].Z() - eElement[0].Z()) - (eElement[3].Y() - eElement[0].Y()) * (eElement[1].Z() - eElement[0].Z());
			b2 = (eElement[1].Z() - eElement[0].Z()) * (eElement[3].X() - eElement[0].X()) - (eElement[3].Z() - eElement[0].Z()) * (eElement[1].X() - eElement[0].X());
			c2 = (eElement[1].X() - eElement[0].X()) * (eElement[3].Y() - eElement[0].Y()) - (eElement[3].X() - eElement[0].X()) * (eElement[1].Y() - eElement[0].Y());
			double a3 = 0; //slope x for plane on the third triangular face of the tetrahedra (nodes B,C,D)
			double b3 = 0; //slope y for plane on the third triangular face of the tetrahedra (nodes B,C,D)
			double c3 = 0; //slope z for plane on the third triangular face of the tetrahedra (nodes B,C,D)
			a3 = (eElement[1].Y() - eElement[2].Y()) * (eElement[3].Z() - eElement[2].Z()) - (eElement[3].Y() - eElement[2].Y()) * (eElement[1].Z() - eElement[2].Z());
			b3 = (eElement[1].Z() - eElement[2].Z()) * (eElement[3].X() - eElement[2].X()) - (eElement[3].Z() - eElement[2].Z()) * (eElement[1].X() - eElement[2].X());
			c3 = (eElement[1].X() - eElement[2].X()) * (eElement[3].Y() - eElement[2].Y()) - (eElement[3].X() - eElement[2].X()) * (eElement[1].Y() - eElement[2].Y());
			double a4 = 0; //slope x for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
			double b4 = 0; //slope y for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
			double c4 = 0; //slope z for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
			a4 = (eElement[0].Y() - eElement[2].Y()) * (eElement[3].Z() - eElement[2].Z()) - (eElement[3].Y() - eElement[2].Y()) * (eElement[0].Z() - eElement[2].Z());
			b4 = (eElement[0].Z() - eElement[2].Z()) * (eElement[3].X() - eElement[2].X()) - (eElement[3].Z() - eElement[2].Z()) * (eElement[0].X() - eElement[2].X());
			c4 = (eElement[0].X() - eElement[2].X()) * (eElement[3].Y() - eElement[2].Y()) - (eElement[3].X() - eElement[2].X()) * (eElement[0].Y() - eElement[2].Y());

			double cosAngle12 = (a1 * a2 + b1 * b2 + c1 * c2) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a2, 2) + pow(b2, 2) + pow(c2, 2)));
			double cosAngle13 = (a1 * a3 + b1 * b3 + c1 * c3) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a3, 2) + pow(b3, 2) + pow(c3, 2)));
			double cosAngle14 = (a1 * a4 + b1 * b4 + c1 * c4) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a4, 2) + pow(b4, 2) + pow(c4, 2)));
			double cosAngle23 = (a3 * a2 + b3 * b2 + c3 * c2) / (sqrt(pow(a3, 2) + pow(b3, 2) + pow(c3, 2)) * sqrt(pow(a2, 2) + pow(b2, 2) + pow(c2, 2)));
			double cosAngle24 = (a4 * a2 + b4 * b2 + c4 * c2) / (sqrt(pow(a4, 2) + pow(b4, 2) + pow(c4, 2)) * sqrt(pow(a2, 2) + pow(b2, 2) + pow(c2, 2)));
			double cosAngle34 = (a4 * a3 + b4 * b3 + c4 * c3) / (sqrt(pow(a4, 2) + pow(b4, 2) + pow(c4, 2)) * sqrt(pow(a3, 2) + pow(b3, 2) + pow(c3, 2)));

			double tolerance = 0.999;
			if (fabs(cosAngle12) > tolerance || fabs(cosAngle13) > tolerance || fabs(cosAngle14) > tolerance || fabs(cosAngle23) > tolerance || fabs(cosAngle24) > tolerance || fabs(cosAngle34) > tolerance) // if two faces are coplanar, I will erase the element (which is probably a sliver)
			{
				if (eElement[0].IsNot(FREE_SURFACE) && eElement[0].IsNot(RIGID) && eElement[0].IsNot(TO_ERASE))
				{
					eElement[0].Set(TO_ERASE);
					erased_nodes += 1;
					inside_nodes_removed++;
					std::cout << "node erased for sliver criterion" << std::endl;
				}
				else if (eElement[1].IsNot(FREE_SURFACE) && eElement[1].IsNot(RIGID) && eElement[1].IsNot(TO_ERASE))
				{
					eElement[1].Set(TO_ERASE);
					erased_nodes += 1;
					inside_nodes_removed++;
					std::cout << "node erased for sliver criterion" << std::endl;
				}
				else if (eElement[2].IsNot(FREE_SURFACE) && eElement[2].IsNot(RIGID) && eElement[2].IsNot(TO_ERASE))
				{
					eElement[2].Set(TO_ERASE);
					erased_nodes += 1;
					inside_nodes_removed++;
					std::cout << "node erased for sliver criterion" << std::endl;
				}
				else if (eElement[3].IsNot(FREE_SURFACE) && eElement[3].IsNot(RIGID) && eElement[3].IsNot(TO_ERASE))
				{
					eElement[3].Set(TO_ERASE);
					erased_nodes += 1;
					inside_nodes_removed++;
					std::cout << "node erased for sliver criterion" << std::endl;
				}
				else
				{
					std::cout << "all nodes are FS or RIGID" << std::endl;
				}
			}
		}

		/// Assignment operator.
		RemoveMeshNodesForFluidsSecondMeshProcess &operator=(RemoveMeshNodesForFluidsSecondMeshProcess const &rOther);

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
									RemoveMeshNodesForFluidsSecondMeshProcess &rThis);

	/// output stream function
	inline std::ostream &operator<<(std::ostream &rOStream,
									const RemoveMeshNodesForFluidsSecondMeshProcess &rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@}

} // namespace Kratos.

#endif // KRATOS_REMOVE_MESH_NODES_FOR_FLUIDS_SECOND_MESH_PROCESS_H_INCLUDED  defined
