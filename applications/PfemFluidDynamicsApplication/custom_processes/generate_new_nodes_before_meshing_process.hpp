//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:                   AFranci $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                October 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined(KRATOS_GENERATE_NEW_NODES_BEFORE_MESHING_PROCESS_H_INCLUDED)
#define KRATOS_GENERATE_NEW_NODES_BEFORE_MESHING_PROCESS_H_INCLUDED

// External includes

// System includes

// Project includes
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "includes/model_part.h"
#include "custom_utilities/mesh_error_calculation_utilities.hpp"
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_processes/mesher_process.hpp"

/// VARIABLES used:
// Data:
// Flags:    (checked)
//           (set)
//           (modified)
//           (reset)
//(set):=(set in this process)

namespace Kratos
{

	///@name Kratos Classes
	///@{

	/// Refine Mesh Elements Process 2D and 3D
	/** The process labels the nodes to be refined (TO_REFINE)
	if the ThresholdVariable  is larger than a ReferenceThreshold
*/

	class GenerateNewNodesBeforeMeshingProcess
		: public MesherProcess
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of Process
		KRATOS_CLASS_POINTER_DEFINITION(GenerateNewNodesBeforeMeshingProcess);

		typedef ModelPart::NodeType NodeType;
		typedef ModelPart::ConditionType ConditionType;
		typedef ModelPart::PropertiesType PropertiesType;
		typedef ConditionType::GeometryType GeometryType;
		typedef std::size_t SizeType;

		///@}
		///@name Life Cycle
		///@{

		/// Default constructor.
		GenerateNewNodesBeforeMeshingProcess(ModelPart &rModelPart,
											 MesherUtilities::MeshingParameters &rRemeshingParameters,
											 int EchoLevel)
			: mrModelPart(rModelPart),
			  mrRemesh(rRemeshingParameters)
		{
			mEchoLevel = EchoLevel;
		}

		/// Destructor.
		virtual ~GenerateNewNodesBeforeMeshingProcess() {}

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
				std::cout << " [ GENERATE NEW NODES for homomgeneous mesh: " << std::endl;

			if (mrModelPart.Name() != mrRemesh.SubModelPartName)
				std::cout << " ModelPart Supplied do not corresponds to the Meshing Domain: (" << mrModelPart.Name() << " != " << mrRemesh.SubModelPartName << ")" << std::endl;

			const ProcessInfo &rCurrentProcessInfo = mrModelPart.GetProcessInfo();
			double currentTime = rCurrentProcessInfo[TIME];
			double timeInterval = rCurrentProcessInfo[DELTA_TIME];
			const SizeType dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

			bool refiningBox = false;
			for (SizeType index = 0; index < mrRemesh.UseRefiningBox.size(); index++)
			{
				if (mrRemesh.UseRefiningBox[index] == true && currentTime > mrRemesh.RefiningBoxInitialTime[index] && currentTime < mrRemesh.RefiningBoxFinalTime[index])
				{
					refiningBox = true;
				}
			}

			if (currentTime < 2 * timeInterval)
			{
				mrRemesh.Info->RemovedNodes = 0;
				mrRemesh.Info->BalancePrincipalSecondaryPartsNodes = 0;

				if (mEchoLevel > 1)
					std::cout << " First meshes: I repare the mesh without adding new nodes" << std::endl;
				mrRemesh.Info->InitialNumberOfNodes = mrRemesh.Info->NumberOfNodes;
			}

			int ElementsToRefine = mrRemesh.Info->RemovedNodes;
			SizeType eulerianInletNodes = mrRemesh.Info->NumberOfEulerianInletNodes;

			int initialNumberOfNodes = mrRemesh.Info->InitialNumberOfNodes;
			int numberOfNodes = mrRemesh.Info->NumberOfNodes;
			int extraNodes = initialNumberOfNodes - numberOfNodes;
			int toleredExtraNodes = int(0.05 * mrRemesh.Info->InitialNumberOfNodes);

			if (mrRemesh.ExecutionOptions.Is(MesherUtilities::REFINE_WALL_CORNER))
			{
				if ((ElementsToRefine - extraNodes) > toleredExtraNodes && refiningBox == false)
				{
					ElementsToRefine = toleredExtraNodes + extraNodes;
					if (ElementsToRefine < 0)
					{
						ElementsToRefine = 0;
					}
				}
			}

			if (ElementsToRefine > 0 && mEchoLevel > 1)
				std::cout << " I will look for " << ElementsToRefine << " new nodes" << std::endl;

			if (refiningBox == false)
			{

				if (ElementsToRefine > 0 || eulerianInletNodes > 0)
				{
					std::vector<array_1d<double, 3>> new_positions;
					std::vector<double> biggest_volumes;
					std::vector<array_1d<SizeType, 4>> nodes_id_to_interpolate;

					int CountNodes = 0;

					new_positions.resize(ElementsToRefine);
					biggest_volumes.resize(ElementsToRefine, false);
					nodes_id_to_interpolate.resize(ElementsToRefine);

					for (int nn = 0; nn < ElementsToRefine; nn++)
					{
						biggest_volumes[nn] = -1.0;
					}

					// std::vector<array_1d<double, 3>> CornerWallNewPositions;
					// std::vector<array_1d<SizeType, 4>> CornerWallNodesIDToInterpolate;
					// std::vector<Node::DofsContainerType> CornerWallNewDofs;
					// int cornerWallNewNodes = 0;
					// int maxOfNewWallNodes = toleredExtraNodes;
					// if (mrRemesh.ExecutionOptions.Is(MesherUtilities::REFINE_WALL_CORNER))
					// {
					// 	CornerWallNewPositions.resize(maxOfNewWallNodes, false);
					// 	CornerWallNodesIDToInterpolate.resize(maxOfNewWallNodes, false);
					// 	CornerWallNewDofs.resize(maxOfNewWallNodes, false);
					// }
					SizeType addedNodesAtEulerianInlet = 0;
					ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin();
					for (ModelPart::ElementsContainerType::const_iterator ie = element_begin; ie != mrModelPart.ElementsEnd(); ie++)
					{

						//////// choose the right (big and safe) elements to refine and compute the new node position and variables ////////
						if (dimension == 2)
						{
							SelectEdgeToRefine2D(ie->GetGeometry(), new_positions, biggest_volumes, nodes_id_to_interpolate, CountNodes, ElementsToRefine, addedNodesAtEulerianInlet);
						}
						else if (dimension == 3)
						{
							SelectEdgeToRefine3D(ie->GetGeometry(), new_positions, biggest_volumes, nodes_id_to_interpolate, CountNodes, ElementsToRefine, addedNodesAtEulerianInlet);
						}

					} // elements loop
					ElementsToRefine += addedNodesAtEulerianInlet;

					mrRemesh.Info->RemovedNodes -= ElementsToRefine;
					if (CountNodes < ElementsToRefine)
					{
						mrRemesh.Info->RemovedNodes += ElementsToRefine - CountNodes;
						new_positions.resize(CountNodes);
						biggest_volumes.resize(CountNodes, false);
						nodes_id_to_interpolate.resize(CountNodes);
					}
					SizeType maxId = 0;
					CreateAndAddNewNodes(new_positions, nodes_id_to_interpolate, ElementsToRefine, maxId);

					// if (mrRemesh.ExecutionOptions.Is(MesherUtilities::REFINE_WALL_CORNER))
					// {
					// 	if (cornerWallNewNodes < maxOfNewWallNodes)
					// 	{
					// 		CornerWallNewPositions.resize(cornerWallNewNodes, false);
					// 		CornerWallNewDofs.resize(cornerWallNewNodes, false);
					// 		CornerWallNodesIDToInterpolate.resize(cornerWallNewNodes, false);
					// 	}
					// 	CreateAndAddNewNodesInCornerWall(CornerWallNewPositions, CornerWallNodesIDToInterpolate, CornerWallNewDofs, cornerWallNewNodes, maxId);
					// }
				}
			}
			else
			{

				std::vector<array_1d<double, 3>> new_positions;
				std::vector<array_1d<SizeType, 4>> nodes_id_to_interpolate;

				int CountNodes = 0;

				new_positions.resize(0);
				nodes_id_to_interpolate.resize(0);

				ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin();
				int nodesInTransitionZone = 0;
				for (ModelPart::ElementsContainerType::const_iterator ie = element_begin; ie != mrModelPart.ElementsEnd(); ie++)
				{
					const SizeType dimension = ie->GetGeometry().WorkingSpaceDimension();
					//////// choose the right (big and safe) elements to refine and compute the new node position and variables ////////
					if (dimension == 2)
					{
						SelectEdgeToRefine2DWithRefinement(ie->GetGeometry(), new_positions, nodes_id_to_interpolate, CountNodes, ElementsToRefine, nodesInTransitionZone);
					}
					else if (dimension == 3)
					{
						SelectEdgeToRefine3DWithRefinement(ie->GetGeometry(), new_positions, nodes_id_to_interpolate, CountNodes, ElementsToRefine, nodesInTransitionZone);
					}

				} // elements loop

				mrRemesh.Info->RemovedNodes -= CountNodes;
				if (CountNodes < ElementsToRefine)
				{
					mrRemesh.Info->RemovedNodes += ElementsToRefine - CountNodes;
					new_positions.resize(CountNodes);
					nodes_id_to_interpolate.resize(CountNodes);
				}
				SizeType maxId = 0;
				CreateAndAddNewNodes(new_positions, nodes_id_to_interpolate, CountNodes, maxId);
			}

			mrRemesh.InputInitializedFlag = false;

			if (mEchoLevel > 1)
				std::cout << "   GENERATE NEW NODES ]; " << std::endl;

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
			return "GenerateNewNodesBeforeMeshingProcess";
		}

		/// Print information about this object.
		void PrintInfo(std::ostream &rOStream) const override
		{
			rOStream << "GenerateNewNodesBeforeMeshingProcess";
		}

		/// Print object's data.s
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
		///@name Private Operators
		///@{

		///@}
		///@name Private Operations
		///@{

		void CopyDofs(Node::DofsContainerType const &From, Node::DofsContainerType &To)
		{
			for (auto &p_dof : From)
			{
				To.push_back(Kratos::unique_ptr<Dof<double>>(new Dof<double>(*p_dof)));
			}
		}

		void InsertNodeInCornerElement2D(Element::GeometryType &Element,
										 std::vector<array_1d<double, 3>> &new_positions,
										 std::vector<array_1d<SizeType, 4>> &nodes_id_to_interpolate,
										 std::vector<Node::DofsContainerType> &NewDofs,
										 int &CountNodes)
		{
			KRATOS_TRY

			const SizeType nds = Element.size();

			SizeType rigidNodes = 0;
			SizeType freesurfaceNodes = 0;

			for (SizeType pn = 0; pn < nds; ++pn)
			{
				if (Element[pn].Is(RIGID))
				{
					rigidNodes++;
				}
				if (Element[pn].Is(FREE_SURFACE))
				{
					freesurfaceNodes++;
				}
			}
			double cosTolerance = 0.01;

			if (rigidNodes == 2 && freesurfaceNodes == 0)
			{
				array_1d<double, 2> NormalA(2, 0.0);
				array_1d<double, 2> NormalB(2, 0.0);
				double cosAngle = 1;

				if ((Element[0].Is(RIGID) && Element[1].Is(RIGID)) || (Element[0].Is(INLET) && Element[1].Is(INLET)))
				{
					NormalA = Element[0].FastGetSolutionStepValue(NORMAL);
					NormalB = Element[1].FastGetSolutionStepValue(NORMAL);
					cosAngle = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1];
					if (cosAngle < cosTolerance && cosAngle > -cosTolerance)
					{
						array_1d<double, 3> new_position = (Element[0].Coordinates() + Element[1].Coordinates()) * 0.5;
						nodes_id_to_interpolate[CountNodes][0] = Element[0].GetId();
						nodes_id_to_interpolate[CountNodes][1] = Element[1].GetId();
						if (Element[2].IsNot(TO_ERASE))
						{
							nodes_id_to_interpolate[CountNodes][2] = Element[2].GetId();
						}
						else
						{
							nodes_id_to_interpolate[CountNodes][2] = Element[0].GetId();
						}
						CopyDofs(Element[2].GetDofs(), NewDofs[CountNodes]);
						new_positions[CountNodes] = new_position;
						CountNodes++;
					}
				}
				else if (Element[0].Is(RIGID) && Element[2].Is(RIGID))
				{
					NormalA = Element[0].FastGetSolutionStepValue(NORMAL);
					NormalB = Element[2].FastGetSolutionStepValue(NORMAL);
					cosAngle = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1];
					if (cosAngle < cosTolerance && cosAngle > -cosTolerance)
					{
						array_1d<double, 3> new_position = (Element[0].Coordinates() + Element[1].Coordinates()) * 0.5;
						nodes_id_to_interpolate[CountNodes][0] = Element[0].GetId();
						nodes_id_to_interpolate[CountNodes][1] = Element[2].GetId();
						if (Element[1].IsNot(TO_ERASE))
						{
							nodes_id_to_interpolate[CountNodes][2] = Element[1].GetId();
						}
						else
						{
							nodes_id_to_interpolate[CountNodes][2] = Element[0].GetId();
						}
						CopyDofs(Element[1].GetDofs(), NewDofs[CountNodes]);
						new_positions[CountNodes] = new_position;
						CountNodes++;
					}
				}
				else if (Element[1].Is(RIGID) && Element[2].Is(RIGID))
				{
					NormalA = Element[1].FastGetSolutionStepValue(NORMAL);
					NormalB = Element[2].FastGetSolutionStepValue(NORMAL);
					cosAngle = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1];
					if (cosAngle < cosTolerance && cosAngle > -cosTolerance)
					{
						array_1d<double, 3> new_position = (Element[2].Coordinates() + Element[1].Coordinates()) * 0.5;
						nodes_id_to_interpolate[CountNodes][0] = Element[2].GetId();
						nodes_id_to_interpolate[CountNodes][1] = Element[1].GetId();
						if (Element[0].IsNot(TO_ERASE))
						{
							nodes_id_to_interpolate[CountNodes][2] = Element[0].GetId();
						}
						else
						{
							nodes_id_to_interpolate[CountNodes][2] = Element[2].GetId();
						}
						CopyDofs(Element[0].GetDofs(), NewDofs[CountNodes]);
						new_positions[CountNodes] = new_position;
						CountNodes++;
					}
				}
			}

			KRATOS_CATCH("")
		}

		void InsertNodeInCornerElement3D(Element::GeometryType &Element,
										 std::vector<array_1d<double, 3>> &new_positions,
										 std::vector<array_1d<SizeType, 4>> &nodes_id_to_interpolate,
										 std::vector<Node::DofsContainerType> &NewDofs,
										 int &CountNodes)
		{
			KRATOS_TRY

			const SizeType nds = Element.size();

			SizeType rigidNodes = 0;
			SizeType freesurfaceNodes = 0;
			SizeType toEraseNodes = 0;

			for (SizeType pn = 0; pn < nds; ++pn)
			{
				if (Element[pn].Is(RIGID))
				{
					rigidNodes++;
				}
				if (Element[pn].Is(FREE_SURFACE))
				{
					freesurfaceNodes++;
				}
				if (Element[pn].Is(TO_ERASE))
				{
					toEraseNodes++;
				}
			}

			if (rigidNodes == 3 && freesurfaceNodes == 0 && toEraseNodes == 0)
			{
				array_1d<double, 3> NormalA(3, 0.0);
				array_1d<double, 3> NormalB(3, 0.0);
				double normNormalA = 0;
				double normNormalB = 0;
				double cos = 1.0;
				double minCos = 1.0;
				array_1d<SizeType, 2> idsWallNodes(2, 0);
				SizeType idFreeNode = 0;
				double cosTolerance = 0.1;
				if (Element[0].IsNot(RIGID))
				{
					NormalA = Element[1].FastGetSolutionStepValue(NORMAL);
					NormalB = Element[2].FastGetSolutionStepValue(NORMAL);
					normNormalA = NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
					normNormalB = NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
					cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
					if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA > 0.99 && normNormalA < 1.01) && (normNormalB > 0.99 && normNormalB < 1.01))
					{
						minCos = cos;
						idsWallNodes[0] = 1;
						idsWallNodes[1] = 2;
						idFreeNode = 0;
					}

					NormalA = Element[1].FastGetSolutionStepValue(NORMAL);
					NormalB = Element[3].FastGetSolutionStepValue(NORMAL);
					normNormalA = NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
					normNormalB = NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
					cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
					if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA > 0.99 && normNormalA < 1.01) && (normNormalB > 0.99 && normNormalB < 1.01))
					{
						minCos = cos;
						idsWallNodes[0] = 1;
						idsWallNodes[1] = 3;
						idFreeNode = 0;
					}

					NormalA = Element[2].FastGetSolutionStepValue(NORMAL);
					NormalB = Element[3].FastGetSolutionStepValue(NORMAL);
					normNormalA = NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
					normNormalB = NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
					cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
					if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA > 0.99 && normNormalA < 1.01) && (normNormalB > 0.99 && normNormalB < 1.01))
					{
						minCos = cos;
						idsWallNodes[0] = 2;
						idsWallNodes[1] = 3;
						idFreeNode = 0;
					}
				}
				else if (Element[1].IsNot(RIGID))
				{
					NormalA = Element[0].FastGetSolutionStepValue(NORMAL);
					NormalB = Element[2].FastGetSolutionStepValue(NORMAL);
					normNormalA = NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
					normNormalB = NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
					cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
					if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA > 0.99 && normNormalA < 1.01) && (normNormalB > 0.99 && normNormalB < 1.01))
					{
						minCos = cos;
						idsWallNodes[0] = 0;
						idsWallNodes[1] = 2;
						idFreeNode = 1;
					}

					NormalA = Element[0].FastGetSolutionStepValue(NORMAL);
					NormalB = Element[3].FastGetSolutionStepValue(NORMAL);
					normNormalA = NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
					normNormalB = NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
					cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
					if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA > 0.99 && normNormalA < 1.01) && (normNormalB > 0.99 && normNormalB < 1.01))
					{
						minCos = cos;
						idsWallNodes[0] = 0;
						idsWallNodes[1] = 3;
						idFreeNode = 1;
					}

					NormalA = Element[2].FastGetSolutionStepValue(NORMAL);
					NormalB = Element[3].FastGetSolutionStepValue(NORMAL);
					normNormalA = NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
					normNormalB = NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
					cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
					if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA > 0.99 && normNormalA < 1.01) && (normNormalB > 0.99 && normNormalB < 1.01))
					{
						minCos = cos;
						idsWallNodes[0] = 2;
						idsWallNodes[1] = 3;
						idFreeNode = 1;
					}
				}
				else if (Element[2].IsNot(RIGID))
				{

					NormalA = Element[0].FastGetSolutionStepValue(NORMAL);
					NormalB = Element[1].FastGetSolutionStepValue(NORMAL);
					normNormalA = NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
					normNormalB = NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
					cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
					if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA > 0.99 && normNormalA < 1.01) && (normNormalB > 0.99 && normNormalB < 1.01))
					{
						minCos = cos;
						idsWallNodes[0] = 0;
						idsWallNodes[1] = 1;
						idFreeNode = 2;
					}

					NormalA = Element[0].FastGetSolutionStepValue(NORMAL);
					NormalB = Element[3].FastGetSolutionStepValue(NORMAL);
					normNormalA = NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
					normNormalB = NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
					cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
					if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA > 0.99 && normNormalA < 1.01) && (normNormalB > 0.99 && normNormalB < 1.01))
					{
						minCos = cos;
						idsWallNodes[0] = 0;
						idsWallNodes[1] = 3;
						idFreeNode = 2;
					}

					NormalA = Element[1].FastGetSolutionStepValue(NORMAL);
					NormalB = Element[3].FastGetSolutionStepValue(NORMAL);
					normNormalA = NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
					normNormalB = NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
					cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
					if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA > 0.99 && normNormalA < 1.01) && (normNormalB > 0.99 && normNormalB < 1.01))
					{
						minCos = cos;
						idsWallNodes[0] = 1;
						idsWallNodes[1] = 3;
						idFreeNode = 2;
					}
				}
				else if (Element[3].IsNot(RIGID))
				{

					NormalA = Element[0].FastGetSolutionStepValue(NORMAL);
					NormalB = Element[1].FastGetSolutionStepValue(NORMAL);
					normNormalA = NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
					normNormalB = NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
					cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
					if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA > 0.99 && normNormalA < 1.01) && (normNormalB > 0.99 && normNormalB < 1.01))
					{
						minCos = cos;
						idsWallNodes[0] = 0;
						idsWallNodes[1] = 1;
						idFreeNode = 3;
					}

					NormalA = Element[0].FastGetSolutionStepValue(NORMAL);
					NormalB = Element[2].FastGetSolutionStepValue(NORMAL);
					normNormalA = NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
					normNormalB = NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
					cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
					if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA > 0.99 && normNormalA < 1.01) && (normNormalB > 0.99 && normNormalB < 1.01))
					{
						minCos = cos;
						idsWallNodes[0] = 0;
						idsWallNodes[1] = 2;
						idFreeNode = 3;
					}

					NormalA = Element[1].FastGetSolutionStepValue(NORMAL);
					NormalB = Element[2].FastGetSolutionStepValue(NORMAL);
					normNormalA = NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
					normNormalB = NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
					cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
					if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA > 0.99 && normNormalA < 1.01) && (normNormalB > 0.99 && normNormalB < 1.01))
					{
						minCos = cos;
						idsWallNodes[0] = 1;
						idsWallNodes[1] = 2;
						idFreeNode = 3;
					}
				}

				if (minCos < cosTolerance && minCos > -cosTolerance)
				{

					bool alreadyAddedNode = false;
					SizeType idA = Element[idsWallNodes[0]].GetId();
					SizeType idB = Element[idsWallNodes[1]].GetId();
					double minimumDistanceToInstert = 1.3 * mrRemesh.Refine->CriticalRadius;
					array_1d<double, 3> CoorDifference = Element[idsWallNodes[0]].Coordinates() - Element[idsWallNodes[1]].Coordinates();
					double SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1];
					double separation = std::sqrt(SquaredLength);
					SizeType idC = Element[idFreeNode].GetId();
					if (separation > minimumDistanceToInstert)
					{

						for (SizeType i = 0; i < unsigned(CountNodes); i++)
						{
							if (idA == nodes_id_to_interpolate[i][0] || idA == nodes_id_to_interpolate[i][1] || idB == nodes_id_to_interpolate[i][0] || idB == nodes_id_to_interpolate[i][1])
							{
								alreadyAddedNode = true;
								break;
							}
						}
						if (alreadyAddedNode == false)
						{
							array_1d<double, 3> new_position = (Element[idsWallNodes[0]].Coordinates() + Element[idsWallNodes[1]].Coordinates()) * 0.5;
							nodes_id_to_interpolate[CountNodes][0] = idA;
							nodes_id_to_interpolate[CountNodes][1] = idB;
							if (Element[idFreeNode].IsNot(TO_ERASE))
							{
								nodes_id_to_interpolate[CountNodes][2] = idC;
							}
							else
							{
								nodes_id_to_interpolate[CountNodes][2] = idA;
							}
							CopyDofs(Element[idFreeNode].GetDofs(), NewDofs[CountNodes]);
							new_positions[CountNodes] = new_position;
							CountNodes++;
						}
					}
				}
			}

			KRATOS_CATCH("")
		}

		void SelectEdgeToRefine2D(Element::GeometryType &Element,
								  std::vector<array_1d<double, 3>> &new_positions,
								  std::vector<double> &biggest_volumes,
								  std::vector<array_1d<SizeType, 4>> &nodes_id_to_interpolate,
								  int &CountNodes,
								  int &ElementsToRefine,
								  SizeType &addedNodesAtEulerianInlet)
		{
			KRATOS_TRY

			const SizeType nds = Element.size();
			double meanMeshSize = mrRemesh.Refine->CriticalRadius;

			SizeType rigidNodes = 0;
			SizeType freesurfaceNodes = 0;
			SizeType lagrangianInletNodes = 0;
			SizeType eulerianInletNodes = 0;
			bool toEraseNodeFound = false;
			double rigidNodeLocalMeshSize = 0;
			double rigidNodeMeshCounter = 0;

			for (SizeType pn = 0; pn < nds; ++pn)
			{
				if (Element[pn].Is(RIGID))
				{
					rigidNodes++;
					rigidNodeLocalMeshSize += Element[pn].FastGetSolutionStepValue(NODAL_H_WALL);
					rigidNodeMeshCounter += 1.0;
				}
				if (Element[pn].Is(TO_ERASE))
				{
					toEraseNodeFound = true;
				}
				if (Element[pn].Is(FREE_SURFACE))
				{
					freesurfaceNodes++;
				}
				if (Element[pn].Is(PFEMFlags::LAGRANGIAN_INLET))
				{
					lagrangianInletNodes++;
				}
				if (Element[pn].Is(PFEMFlags::EULERIAN_INLET))
				{
					eulerianInletNodes++;
				}
			}

			if (rigidNodeMeshCounter > 0)
			{
				const double rigidWallMeshSize = rigidNodeLocalMeshSize / rigidNodeMeshCounter;
				const double ratio = rigidWallMeshSize / meanMeshSize;
				const double tolerance = 1.8;
				if (ratio > tolerance)
				{
					meanMeshSize *= 0.5;
					meanMeshSize += 0.5 * rigidWallMeshSize;
				}
			}

			double penalization = 1.0; // to penalize adding node, penalization here should be smaller than 1
			if (rigidNodes > 1)
			{
				penalization = 0.8;
				if (lagrangianInletNodes > 0)
				{
					penalization = 0.9;
				}
			}
			else if (rigidNodes > 0 && freesurfaceNodes > 0 && eulerianInletNodes == 0)
			{
				penalization = 0;
			}
			else if (freesurfaceNodes > 0)
			{
				penalization = 0.875;
			}

			array_1d<double, 3> Edges(3, 0.0);
			array_1d<SizeType, 3> FirstEdgeNode(3, 0);
			array_1d<SizeType, 3> SecondEdgeNode(3, 0);
			double WallCharacteristicDistance = 0;
			ComputeWallCharacteristicDistance2D(Element, WallCharacteristicDistance, Edges, FirstEdgeNode, SecondEdgeNode);

			bool dangerousElement = false;
			DetectDangerousElements2D(Element, WallCharacteristicDistance, Edges, FirstEdgeNode, SecondEdgeNode, rigidNodes, penalization, dangerousElement);

			if (dangerousElement == false && toEraseNodeFound == false)
			{
				ManageDangerousElements2D(Element,
										  new_positions,
										  biggest_volumes,
										  nodes_id_to_interpolate,
										  CountNodes,
										  ElementsToRefine,
										  meanMeshSize,
										  addedNodesAtEulerianInlet,
										  eulerianInletNodes,
										  rigidNodes,
										  freesurfaceNodes,
										  Edges,
										  FirstEdgeNode,
										  SecondEdgeNode,
										  penalization);
			}

			KRATOS_CATCH("")
		}

		void ComputeWallCharacteristicDistance2D(Element::GeometryType &Element,
												 double &WallCharacteristicDistance,
												 array_1d<double, 3> &Edges,
												 array_1d<SizeType, 3> &FirstEdgeNode,
												 array_1d<SizeType, 3> &SecondEdgeNode)
		{
			KRATOS_TRY
			const SizeType nds = Element.size();
			array_1d<double, 3> CoorDifference = Element[1].Coordinates() - Element[0].Coordinates();
			double SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1];
			Edges[0] = std::sqrt(SquaredLength);
			FirstEdgeNode[0] = 0;
			SecondEdgeNode[0] = 1;
			if ((Element[0].Is(RIGID) && Element[1].Is(RIGID)) || (Element[0].Is(INLET) && Element[1].Is(INLET)))
			{
				WallCharacteristicDistance = Edges[0];
			}
			SizeType Counter = 0;
			for (SizeType i = 2; i < nds; i++)
			{
				for (SizeType j = 0; j < i; j++)
				{
					noalias(CoorDifference) = Element[i].Coordinates() - Element[j].Coordinates();
					SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1];
					Counter += 1;
					Edges[Counter] = std::sqrt(SquaredLength);
					FirstEdgeNode[Counter] = j;
					SecondEdgeNode[Counter] = i;
					if (((Element[i].Is(RIGID) && Element[j].Is(RIGID)) || (Element[i].Is(INLET) && Element[j].Is(INLET))) && Edges[Counter] > WallCharacteristicDistance)
					{
						WallCharacteristicDistance = Edges[Counter];
					}
				}
			}

			KRATOS_CATCH("")
		}

		void ComputeWallCharacteristicDistance2DWithRefinement(Element::GeometryType &Element,
															   double &WallCharacteristicDistance,
															   array_1d<double, 3> &Edges,
															   array_1d<SizeType, 3> &FirstEdgeNode,
															   array_1d<SizeType, 3> &SecondEdgeNode)
		{
			KRATOS_TRY
			const SizeType nds = Element.size();
			array_1d<double, 3> CoorDifference = Element[1].Coordinates() - Element[0].Coordinates();
			double SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1];
			Edges[0] = std::sqrt(SquaredLength);
			FirstEdgeNode[0] = 0;
			SecondEdgeNode[0] = 1;
			if ((Element[0].Is(RIGID) && Element[1].Is(RIGID)) || (Element[0].Is(INLET) && Element[1].Is(INLET)))
			{
				WallCharacteristicDistance = Edges[0];
			}
			SizeType Counter = 0;
			for (SizeType i = 2; i < nds; i++)
			{
				for (SizeType j = 0; j < i; j++)
				{
					noalias(CoorDifference) = Element[i].Coordinates() - Element[j].Coordinates();
					SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1];
					Counter += 1;
					Edges[Counter] = std::sqrt(SquaredLength);
					FirstEdgeNode[Counter] = j;
					SecondEdgeNode[Counter] = i;
					if (Element[i].Is(RIGID) && Element[j].Is(RIGID) && Edges[Counter] > WallCharacteristicDistance)
					{
						WallCharacteristicDistance = Edges[Counter];
					}
				}
			}

			KRATOS_CATCH("")
		}

		void DetectDangerousElements2D(Element::GeometryType &Element,
									   double &WallCharacteristicDistance,
									   array_1d<double, 3> &Edges,
									   array_1d<SizeType, 3> &FirstEdgeNode,
									   array_1d<SizeType, 3> &SecondEdgeNode,
									   SizeType rigidNodes,
									   double &penalization,
									   bool &dangerousElement)
		{
			KRATOS_TRY
			const double safetyCoefficient2D = 1.5;

			for (SizeType i = 0; i < 3; i++)
			{
				if (rigidNodes > 1)
				{
					if ((Edges[i] < WallCharacteristicDistance * safetyCoefficient2D && (Element[FirstEdgeNode[i]].Is(RIGID) || Element[SecondEdgeNode[i]].Is(RIGID))) ||
						((Element[FirstEdgeNode[i]].Is(RIGID) && Element[SecondEdgeNode[i]].Is(RIGID)) || (Element[FirstEdgeNode[i]].Is(INLET) && Element[SecondEdgeNode[i]].Is(INLET))))
					{
						Edges[i] = 0;
					}
					if ((Element[FirstEdgeNode[i]].Is(FREE_SURFACE) || Element[FirstEdgeNode[i]].Is(RIGID)) &&
						(Element[SecondEdgeNode[i]].Is(FREE_SURFACE) || Element[SecondEdgeNode[i]].Is(RIGID)) &&
						(Element[FirstEdgeNode[i]].IsNot(PFEMFlags::EULERIAN_INLET) && Element[SecondEdgeNode[i]].IsNot(PFEMFlags::EULERIAN_INLET)))
					{
						Edges[i] = 0;
					}
				}
				else if (rigidNodes == 0)
				{
					SizeType propertyIdFirstNode = Element[FirstEdgeNode[i]].FastGetSolutionStepValue(PROPERTY_ID);
					SizeType propertyIdSecondNode = Element[SecondEdgeNode[i]].FastGetSolutionStepValue(PROPERTY_ID);
					if (propertyIdFirstNode != propertyIdSecondNode)
					{
						penalization = 0.9; // 10% less than normal nodes
					}
				}
			}
			if ((Edges[0] == 0 && Edges[1] == 0 && Edges[2] == 0) || rigidNodes == 3)
			{
				dangerousElement = true;
			}

			KRATOS_CATCH("")
		}

		void DetectDangerousElements2DWithRefinement(Element::GeometryType &Element,
													 double &WallCharacteristicDistance,
													 array_1d<double, 3> &Edges,
													 array_1d<SizeType, 3> &FirstEdgeNode,
													 array_1d<SizeType, 3> &SecondEdgeNode,
													 SizeType rigidNodes,
													 double &penalization,
													 bool &dangerousElement)
		{
			KRATOS_TRY
			const double safetyCoefficient2D = 1.5;

			for (SizeType i = 0; i < 3; i++)
			{
				if (rigidNodes > 1)
				{
					if ((Edges[i] < WallCharacteristicDistance * safetyCoefficient2D && (Element[FirstEdgeNode[i]].Is(RIGID) || Element[SecondEdgeNode[i]].Is(RIGID))) ||
						(Element[FirstEdgeNode[i]].Is(RIGID) && Element[SecondEdgeNode[i]].Is(RIGID)))
					{
						Edges[i] = 0;
					}
					if ((Element[FirstEdgeNode[i]].Is(FREE_SURFACE) || Element[FirstEdgeNode[i]].Is(RIGID)) &&
						(Element[SecondEdgeNode[i]].Is(FREE_SURFACE) || Element[SecondEdgeNode[i]].Is(RIGID)))
					{
						Edges[i] = 0;
					}
				}
				else if (rigidNodes == 0)
				{
					SizeType propertyIdFirstNode = Element[FirstEdgeNode[i]].FastGetSolutionStepValue(PROPERTY_ID);
					SizeType propertyIdSecondNode = Element[SecondEdgeNode[i]].FastGetSolutionStepValue(PROPERTY_ID);
					if (propertyIdFirstNode != propertyIdSecondNode)
					{
						penalization = 1.1; // 10% more than normal nodes
					}
				}
			}
			if ((Edges[0] == 0 && Edges[1] == 0 && Edges[2] == 0) || rigidNodes == 3)
			{
				dangerousElement = true;
			}

			KRATOS_CATCH("")
		}

		void ManageDangerousElements2D(Element::GeometryType &Element,
									   std::vector<array_1d<double, 3>> &new_positions,
									   std::vector<double> &biggest_volumes,
									   std::vector<array_1d<SizeType, 4>> &nodes_id_to_interpolate,
									   int &CountNodes,
									   int &ElementsToRefine,
									   double &meanMeshSize,
									   SizeType &addedNodesAtEulerianInlet,
									   SizeType eulerianInletNodes,
									   SizeType rigidNodes,
									   SizeType freesurfaceNodes,
									   array_1d<double, 3> &Edges,
									   array_1d<SizeType, 3> &FirstEdgeNode,
									   array_1d<SizeType, 3> &SecondEdgeNode,
									   double &penalization)
		{
			KRATOS_TRY

			SizeType maxCount = 3;
			double LargestEdge = 0;
			const double limitEdgeLength = 1.4 * meanMeshSize;
			double ElementalVolume = Element.Area();
			bool suitableElementForSecondAdd = true;

			for (SizeType i = 0; i < 3; i++)
			{
				if (Edges[i] > LargestEdge)
				{
					maxCount = i;
					LargestEdge = Edges[i];
				}
			}

			if (CountNodes < ElementsToRefine && LargestEdge > limitEdgeLength && eulerianInletNodes == 0)
			{

				array_1d<double, 3> new_position = (Element[FirstEdgeNode[maxCount]].Coordinates() + Element[SecondEdgeNode[maxCount]].Coordinates()) * 0.5;

				bool suitableElement = true;
				for (int j = 0; j < CountNodes; j++)
				{
					const double diffX = std::abs(new_positions[j][0] - new_position[0]) - meanMeshSize * 0.5;
					const double diffY = std::abs(new_positions[j][1] - new_position[1]) - meanMeshSize * 0.5;
					if (diffX < 0 && diffY < 0) //  the node is in the same zone of a previously inserted node
					{
						suitableElement = false;
					}
				}

				if (suitableElement)
				{
					nodes_id_to_interpolate[CountNodes][0] = Element[FirstEdgeNode[maxCount]].GetId();
					nodes_id_to_interpolate[CountNodes][1] = Element[SecondEdgeNode[maxCount]].GetId();
					biggest_volumes[CountNodes] = ElementalVolume;
					new_positions[CountNodes] = new_position;
					CountNodes++;
				}
			}
			else if (freesurfaceNodes < 3 && rigidNodes < 3 && eulerianInletNodes == 0)
			{
				ElementalVolume *= penalization;
				for (int nn = 0; nn < ElementsToRefine; nn++)
				{
					if (ElementalVolume > biggest_volumes[nn])
					{

						if (maxCount < 3 && LargestEdge > limitEdgeLength)
						{
							array_1d<double, 3> new_position = (Element[FirstEdgeNode[maxCount]].Coordinates() + Element[SecondEdgeNode[maxCount]].Coordinates()) * 0.5;
							if (ElementsToRefine > 1 && CountNodes > 0)
							{
								for (int j = 0; j < ElementsToRefine; j++)
								{
									const double diffX = std::abs(new_positions[j][0] - new_position[0]) - meanMeshSize * 0.5;
									const double diffY = std::abs(new_positions[j][1] - new_position[1]) - meanMeshSize * 0.5;
									if (diffX < 0 && diffY < 0) // the node is in the same zone of a previously inserted node
									{
										suitableElementForSecondAdd = false;
									}
								}
							}

							if (suitableElementForSecondAdd)
							{
								nodes_id_to_interpolate[nn][0] = Element[FirstEdgeNode[maxCount]].GetId();
								nodes_id_to_interpolate[nn][1] = Element[SecondEdgeNode[maxCount]].GetId();
								biggest_volumes[nn] = ElementalVolume;
								new_positions[nn] = new_position;
							}
						}

						break;
					}
				}
			}

			if (eulerianInletNodes > 0 && LargestEdge > (2.0 * meanMeshSize))
			{

				array_1d<double, 3> new_position = (Element[FirstEdgeNode[maxCount]].Coordinates() + Element[SecondEdgeNode[maxCount]].Coordinates()) * 0.5;

				bool suitableElement = true;
				for (int j = 0; j < CountNodes; j++)
				{
					const double diffX = std::abs(new_positions[j][0] - new_position[0]) - meanMeshSize * 0.5;
					const double diffY = std::abs(new_positions[j][1] - new_position[1]) - meanMeshSize * 0.5;
					// if (diffX < 0 && diffY < 0)																															//  the node is in the same zone of a previously inserted node
					if ((diffX < 0 && diffY < 0) || (Element[FirstEdgeNode[maxCount]].IsNot(INLET) && Element[SecondEdgeNode[maxCount]].IsNot(INLET))) //  the node is in the same zone of a previously inserted node

					{
						suitableElement = false;
					}
				}

				if (suitableElement)
				{
					if (CountNodes >= ElementsToRefine)
					{
						new_positions.resize(CountNodes + 1);
						biggest_volumes.resize(CountNodes + 1, false);
						nodes_id_to_interpolate.resize(CountNodes + 1);
						addedNodesAtEulerianInlet++;
					}
					nodes_id_to_interpolate[CountNodes][0] = Element[FirstEdgeNode[maxCount]].GetId();
					nodes_id_to_interpolate[CountNodes][1] = Element[SecondEdgeNode[maxCount]].GetId();
					biggest_volumes[CountNodes] = ElementalVolume;
					new_positions[CountNodes] = new_position;
					CountNodes++;
				}
			}

			KRATOS_CATCH("")
		}

		void SelectEdgeToRefine3D(Element::GeometryType &Element,
								  std::vector<array_1d<double, 3>> &new_positions,
								  std::vector<double> &biggest_volumes,
								  std::vector<array_1d<SizeType, 4>> &nodes_id_to_interpolate,
								  int &CountNodes,
								  int &ElementsToRefine,
								  SizeType &addedNodesAtEulerianInlet)
		{
			KRATOS_TRY

			const SizeType nds = Element.size();
			double meanMeshSize = mrRemesh.Refine->CriticalRadius;

			SizeType rigidNodes = 0;
			SizeType freesurfaceNodes = 0;
			SizeType lagrangianInletNodes = 0;
			SizeType eulerianInletNodes = 0;
			bool toEraseNodeFound = false;
			double rigidNodeLocalMeshSize = 0;
			double rigidNodeMeshCounter = 0;

			for (SizeType pn = 0; pn < nds; ++pn)
			{
				if (Element[pn].Is(RIGID))
				{
					rigidNodes++;
					rigidNodeLocalMeshSize += Element[pn].FastGetSolutionStepValue(NODAL_H_WALL);
					rigidNodeMeshCounter += 1.0;
				}
				if (Element[pn].Is(TO_ERASE))
				{
					toEraseNodeFound = true;
				}
				if (Element[pn].Is(FREE_SURFACE))
				{
					freesurfaceNodes++;
				}
				if (Element[pn].Is(PFEMFlags::LAGRANGIAN_INLET))
				{
					lagrangianInletNodes++;
				}
				if (Element[pn].Is(PFEMFlags::EULERIAN_INLET))
				{
					eulerianInletNodes++;
				}
			}

			if (rigidNodeMeshCounter > 0)
			{
				const double rigidWallMeshSize = rigidNodeLocalMeshSize / rigidNodeMeshCounter;
				const double ratio = rigidWallMeshSize / meanMeshSize;
				const double tolerance = 1.8;
				if (ratio > tolerance)
				{
					meanMeshSize *= 0.5;
					meanMeshSize += 0.5 * rigidWallMeshSize;
				}
			}

			double penalization = 1.0; // penalization here should be smaller than 1
			if (rigidNodes > 2)
			{
				penalization = 0.7;
				if (lagrangianInletNodes > 0)
				{
					penalization = 0.9;
				}
			}
			else if (rigidNodes > 0 && freesurfaceNodes > 0 && eulerianInletNodes == 0)
			{
				penalization = 0;
			}
			else if (freesurfaceNodes > 0)
			{
				penalization = 0.95;
			}

			array_1d<double, 6> Edges(6, 0.0);
			array_1d<SizeType, 6> FirstEdgeNode(6, 0);
			array_1d<SizeType, 6> SecondEdgeNode(6, 0);
			double WallCharacteristicDistance = 0;
			ComputeWallCharacteristicDistance3D(Element, WallCharacteristicDistance, Edges, FirstEdgeNode, SecondEdgeNode);

			// Edges connectivity: Edges[0]=d01, Edges[1]=d20, Edges[2]=d21, Edges[3]=d30, Edges[4]=d31, Edges[5]=d32
			bool dangerousElement = false;
			DetectDangerousElements3D(Element, WallCharacteristicDistance, Edges, FirstEdgeNode, SecondEdgeNode, rigidNodes, penalization, dangerousElement);

			// just to fill the vector
			if (dangerousElement == false && toEraseNodeFound == false)
			{
				ManageDangerousElements3D(Element,
										  new_positions,
										  biggest_volumes,
										  nodes_id_to_interpolate,
										  CountNodes,
										  ElementsToRefine,
										  meanMeshSize,
										  addedNodesAtEulerianInlet,
										  eulerianInletNodes,
										  rigidNodes,
										  freesurfaceNodes,
										  Edges,
										  FirstEdgeNode,
										  SecondEdgeNode,
										  penalization);
			}

			KRATOS_CATCH("")
		}

		void ComputeWallCharacteristicDistance3D(Element::GeometryType &Element,
												 double &WallCharacteristicDistance,
												 array_1d<double, 6> &Edges,
												 array_1d<SizeType, 6> &FirstEdgeNode,
												 array_1d<SizeType, 6> &SecondEdgeNode)
		{
			KRATOS_TRY
			const SizeType nds = Element.size();
			array_1d<double, 3> CoorDifference = Element[1].Coordinates() - Element[0].Coordinates();
			double SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1] + CoorDifference[2] * CoorDifference[2];
			Edges[0] = std::sqrt(SquaredLength);
			FirstEdgeNode[0] = 0;
			SecondEdgeNode[0] = 1;
			if ((Element[0].Is(RIGID) && Element[1].Is(RIGID)) || (Element[0].Is(INLET) && Element[1].Is(INLET)))
			{
				WallCharacteristicDistance = Edges[0];
			}
			SizeType Counter = 0;
			for (SizeType i = 2; i < nds; i++)
			{
				for (SizeType j = 0; j < i; j++)
				{
					noalias(CoorDifference) = Element[i].Coordinates() - Element[j].Coordinates();
					SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1] + CoorDifference[2] * CoorDifference[2];
					Counter += 1;
					Edges[Counter] = std::sqrt(SquaredLength);
					FirstEdgeNode[Counter] = j;
					SecondEdgeNode[Counter] = i;
					if (((Element[i].Is(RIGID) && Element[j].Is(RIGID)) || (Element[i].Is(INLET) && Element[j].Is(INLET))) && Edges[Counter] > WallCharacteristicDistance)
					{
						WallCharacteristicDistance = Edges[Counter];
					}
				}
			}

			KRATOS_CATCH("")
		}

		void ComputeWallCharacteristicDistance3DWithRefinement(Element::GeometryType &Element,
															   double &WallCharacteristicDistance,
															   array_1d<double, 6> &Edges,
															   array_1d<SizeType, 6> &FirstEdgeNode,
															   array_1d<SizeType, 6> &SecondEdgeNode)
		{
			KRATOS_TRY
			const SizeType nds = Element.size();
			array_1d<double, 3> CoorDifference = Element[1].Coordinates() - Element[0].Coordinates();
			double SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1] + CoorDifference[2] * CoorDifference[2];
			Edges[0] = std::sqrt(SquaredLength);
			FirstEdgeNode[0] = 0;
			SecondEdgeNode[0] = 1;
			if ((Element[0].Is(RIGID) && Element[1].Is(RIGID)) || (Element[0].Is(INLET) && Element[1].Is(INLET)))
			{
				WallCharacteristicDistance = Edges[0];
			}
			SizeType Counter = 0;
			for (SizeType i = 2; i < nds; i++)
			{
				for (SizeType j = 0; j < i; j++)
				{
					noalias(CoorDifference) = Element[i].Coordinates() - Element[j].Coordinates();
					SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1] + CoorDifference[2] * CoorDifference[2];
					Counter += 1;
					Edges[Counter] = std::sqrt(SquaredLength);
					FirstEdgeNode[Counter] = j;
					SecondEdgeNode[Counter] = i;
					if (Element[i].Is(RIGID) && Element[j].Is(RIGID) && Edges[Counter] > WallCharacteristicDistance)
					{
						WallCharacteristicDistance = Edges[Counter];
					}
				}
			}

			KRATOS_CATCH("")
		}

		void DetectDangerousElements3D(Element::GeometryType &Element,
									   double &WallCharacteristicDistance,
									   array_1d<double, 6> &Edges,
									   array_1d<SizeType, 6> &FirstEdgeNode,
									   array_1d<SizeType, 6> &SecondEdgeNode,
									   SizeType rigidNodes,
									   double &penalization,
									   bool &dangerousElement)
		{
			KRATOS_TRY
			const double safetyCoefficient3D = 1.6;

			for (SizeType i = 0; i < 6; i++)
			{
				if (rigidNodes > 1)
				{
					if ((Edges[i] < WallCharacteristicDistance * safetyCoefficient3D && (Element[FirstEdgeNode[i]].Is(RIGID) || Element[SecondEdgeNode[i]].Is(RIGID))) ||
						((Element[FirstEdgeNode[i]].Is(RIGID) && Element[SecondEdgeNode[i]].Is(RIGID)) || (Element[FirstEdgeNode[i]].Is(INLET) && Element[SecondEdgeNode[i]].Is(INLET))))
					{
						Edges[i] = 0;
					}
					if ((Element[FirstEdgeNode[i]].Is(FREE_SURFACE) || Element[FirstEdgeNode[i]].Is(RIGID)) &&
						(Element[SecondEdgeNode[i]].Is(FREE_SURFACE) || Element[SecondEdgeNode[i]].Is(RIGID)) &&
						(Element[FirstEdgeNode[i]].IsNot(PFEMFlags::EULERIAN_INLET) && Element[SecondEdgeNode[i]].IsNot(PFEMFlags::EULERIAN_INLET)))
					{
						Edges[i] = 0;
					}
				}
				else if (rigidNodes == 0)
				{
					SizeType propertyIdFirstNode = Element[FirstEdgeNode[i]].FastGetSolutionStepValue(PROPERTY_ID);
					SizeType propertyIdSecondNode = Element[SecondEdgeNode[i]].FastGetSolutionStepValue(PROPERTY_ID);
					if (propertyIdFirstNode != propertyIdSecondNode)
					{
						penalization = 0.8; // 20% less than normal nodes
					}
				}
			}
			if (rigidNodes == 1)
			{
				if (Element[0].Is(RIGID))
				{
					Edges[0] = 0;
					Edges[1] = 0;
					Edges[3] = 0;
				}
				if (Element[1].Is(RIGID))
				{
					Edges[0] = 0;
					Edges[2] = 0;
					Edges[4] = 0;
				}
				if (Element[2].Is(RIGID))
				{
					Edges[1] = 0;
					Edges[2] = 0;
					Edges[5] = 0;
				}
				if (Element[3].Is(RIGID))
				{
					Edges[3] = 0;
					Edges[4] = 0;
					Edges[5] = 0;
				}
			}

			if ((Edges[0] == 0 && Edges[1] == 0 && Edges[2] == 0 && Edges[3] == 0 && Edges[4] == 0 && Edges[5] == 0) || rigidNodes > 2)
			{
				dangerousElement = true;
			}

			KRATOS_CATCH("")
		}

		void DetectDangerousElements3DWithRefinement(Element::GeometryType &Element,
													 double &WallCharacteristicDistance,
													 array_1d<double, 6> &Edges,
													 array_1d<SizeType, 6> &FirstEdgeNode,
													 array_1d<SizeType, 6> &SecondEdgeNode,
													 SizeType rigidNodes,
													 double &penalization,
													 bool &dangerousElement)
		{
			KRATOS_TRY
			const double safetyCoefficient3D = 1.6;

						for (SizeType i = 0; i < 6; i++)
			{
				if (rigidNodes > 1)
				{
					if ((Edges[i] < WallCharacteristicDistance * safetyCoefficient3D && (Element[FirstEdgeNode[i]].Is(RIGID) || Element[SecondEdgeNode[i]].Is(RIGID))) ||
						(Element[FirstEdgeNode[i]].Is(RIGID) && Element[SecondEdgeNode[i]].Is(RIGID)))
					{
						Edges[i] = 0;
					}
					if ((Element[FirstEdgeNode[i]].Is(FREE_SURFACE) || Element[FirstEdgeNode[i]].Is(RIGID)) &&
						(Element[SecondEdgeNode[i]].Is(FREE_SURFACE) || Element[SecondEdgeNode[i]].Is(RIGID)))
					{
						Edges[i] = 0;
					}
				}
				else if (rigidNodes == 0)
				{
					SizeType propertyIdFirstNode = Element[FirstEdgeNode[i]].FastGetSolutionStepValue(PROPERTY_ID);
					SizeType propertyIdSecondNode = Element[SecondEdgeNode[i]].FastGetSolutionStepValue(PROPERTY_ID);
					if (propertyIdFirstNode != propertyIdSecondNode)
					{
						penalization = 1.2; // 20% less than normal nodes
					}
				}
			}
			if (rigidNodes == 1)
			{
				if (Element[0].Is(RIGID))
				{
					Edges[0] = 0;
					Edges[1] = 0;
					Edges[3] = 0;
				}
				if (Element[1].Is(RIGID))
				{
					Edges[0] = 0;
					Edges[2] = 0;
					Edges[4] = 0;
				}
				if (Element[2].Is(RIGID))
				{
					Edges[1] = 0;
					Edges[2] = 0;
					Edges[5] = 0;
				}
				if (Element[3].Is(RIGID))
				{
					Edges[3] = 0;
					Edges[4] = 0;
					Edges[5] = 0;
				}
			}

			if ((Edges[0] == 0 && Edges[1] == 0 && Edges[2] == 0 && Edges[3] == 0 && Edges[4] == 0 && Edges[5] == 0) || rigidNodes > 2)
			{
				dangerousElement = true;
			}

			KRATOS_CATCH("")
		}

		void ManageDangerousElements3D(Element::GeometryType &Element,
									   std::vector<array_1d<double, 3>> &new_positions,
									   std::vector<double> &biggest_volumes,
									   std::vector<array_1d<SizeType, 4>> &nodes_id_to_interpolate,
									   int &CountNodes,
									   int &ElementsToRefine,
									   double &meanMeshSize,
									   SizeType &addedNodesAtEulerianInlet,
									   SizeType eulerianInletNodes,
									   SizeType rigidNodes,
									   SizeType freesurfaceNodes,
									   array_1d<double, 6> &Edges,
									   array_1d<SizeType, 6> &FirstEdgeNode,
									   array_1d<SizeType, 6> &SecondEdgeNode,
									   double &penalization)
		{
			KRATOS_TRY

			SizeType maxCount = 6;
			double LargestEdge = 0;
			const double limitEdgeLength = 1.25 * meanMeshSize;
			double ElementalVolume = Element.Volume();
			bool suitableElementForSecondAdd = true;

			for (SizeType i = 0; i < 6; i++)
			{
				if (Edges[i] > LargestEdge)
				{
					maxCount = i;
					LargestEdge = Edges[i];
				}
			}
			if ((CountNodes < ElementsToRefine && LargestEdge > limitEdgeLength && eulerianInletNodes == 0))
			{

				array_1d<double, 3> new_position = (Element[FirstEdgeNode[maxCount]].Coordinates() + Element[SecondEdgeNode[maxCount]].Coordinates()) * 0.5;

				bool suitableElement = true;
				for (int j = 0; j < CountNodes; j++)
				{
					const double diffX = std::abs(new_positions[j][0] - new_position[0]) - meanMeshSize * 0.5;
					const double diffY = std::abs(new_positions[j][1] - new_position[1]) - meanMeshSize * 0.5;
					const double diffZ = std::abs(new_positions[j][2] - new_position[2]) - meanMeshSize * 0.5;
					if (diffX < 0 && diffY < 0 && diffZ < 0) //  the node is in the same zone of a previously inserted node
					{
						suitableElement = false;
					}
				}

				if (suitableElement)
				{
					nodes_id_to_interpolate[CountNodes][0] = Element[FirstEdgeNode[maxCount]].GetId();
					nodes_id_to_interpolate[CountNodes][1] = Element[SecondEdgeNode[maxCount]].GetId();
					biggest_volumes[CountNodes] = ElementalVolume;
					new_positions[CountNodes] = new_position;
					CountNodes++;
				}
			}
			else if (freesurfaceNodes < 4 && rigidNodes < 4 && eulerianInletNodes == 0)
			{

				ElementalVolume *= penalization;
				for (int nn = 0; nn < ElementsToRefine; nn++)
				{
					if (ElementalVolume > biggest_volumes[nn])
					{
						if (maxCount < 6 && LargestEdge > limitEdgeLength)
						{
							array_1d<double, 3> new_position = (Element[FirstEdgeNode[maxCount]].Coordinates() + Element[SecondEdgeNode[maxCount]].Coordinates()) * 0.5;

							if (ElementsToRefine > 1 && CountNodes > 0)
							{
								for (int j = 0; j < ElementsToRefine; j++)
								{
									const double diffX = std::abs(new_positions[j][0] - new_position[0]) - meanMeshSize * 0.5;
									const double diffY = std::abs(new_positions[j][1] - new_position[1]) - meanMeshSize * 0.5;
									const double diffZ = std::abs(new_positions[j][2] - new_position[2]) - meanMeshSize * 0.5;
									if (diffX < 0 && diffY < 0 && diffZ < 0) // the node is in the same zone of a previously inserted node
									{
										suitableElementForSecondAdd = false;
									}
								}
							}

							if (suitableElementForSecondAdd)
							{
								nodes_id_to_interpolate[nn][0] = Element[FirstEdgeNode[maxCount]].GetId();
								nodes_id_to_interpolate[nn][1] = Element[SecondEdgeNode[maxCount]].GetId();
								biggest_volumes[nn] = ElementalVolume;
								new_positions[nn] = new_position;
							}
						}

						break;
					}
				}
			}

			if (eulerianInletNodes > 0 && LargestEdge > (1.7 * meanMeshSize))
			{

				array_1d<double, 3> new_position = (Element[FirstEdgeNode[maxCount]].Coordinates() + Element[SecondEdgeNode[maxCount]].Coordinates()) * 0.5;

				bool suitableElement = true;
				for (int j = 0; j < CountNodes; j++)
				{
					const double diffX = std::abs(new_positions[j][0] - new_position[0]) - meanMeshSize * 0.5;
					const double diffY = std::abs(new_positions[j][1] - new_position[1]) - meanMeshSize * 0.5;
					const double diffZ = std::abs(new_positions[j][2] - new_position[2]) - meanMeshSize * 0.5;
					if ((diffX < 0 && diffY < 0 && diffZ < 0) || (Element[FirstEdgeNode[maxCount]].IsNot(INLET) && Element[SecondEdgeNode[maxCount]].IsNot(INLET))) //  the node is in the same zone of a previously inserted node
					{
						suitableElement = false;
					}
				}

				if (suitableElement)
				{
					if (CountNodes >= ElementsToRefine)
					{
						new_positions.resize(CountNodes + 1);
						biggest_volumes.resize(CountNodes + 1, false);
						nodes_id_to_interpolate.resize(CountNodes + 1);
						addedNodesAtEulerianInlet++;
					}
					nodes_id_to_interpolate[CountNodes][0] = Element[FirstEdgeNode[maxCount]].GetId();
					nodes_id_to_interpolate[CountNodes][1] = Element[SecondEdgeNode[maxCount]].GetId();
					biggest_volumes[CountNodes] = ElementalVolume;
					new_positions[CountNodes] = new_position;
					CountNodes++;
				}
			}

			KRATOS_CATCH("")
		}

		void SelectEdgeToRefine2DWithRefinement(Element::GeometryType &Element,
												std::vector<array_1d<double, 3>> &new_positions,
												std::vector<array_1d<SizeType, 4>> &nodes_id_to_interpolate,
												int &CountNodes,
												const int ElementsToRefine,
												int &nodesInTransitionZone)
		{
			KRATOS_TRY

			const SizeType nds = Element.size();

			SizeType rigidNodes = 0;
			SizeType freesurfaceNodes = 0;
			bool toEraseNodeFound = false;
			double rigidNodeLocalMeshSize = 0;
			double rigidNodeMeshCounter = 0;
			double meanMeshSize = mrRemesh.Refine->CriticalRadius;
			const ProcessInfo &rCurrentProcessInfo = mrModelPart.GetProcessInfo();
			const double currentTime = rCurrentProcessInfo[TIME];
			bool insideTransitionZone = false;
			for (SizeType pn = 0; pn < nds; ++pn)
			{
				mMesherUtilities.DefineMeshSizeInTransitionZones2D(mrRemesh, currentTime, Element[pn].Coordinates(), meanMeshSize, insideTransitionZone);
			}

			for (SizeType pn = 0; pn < nds; ++pn)
			{
				if (Element[pn].Is(RIGID))
				{
					rigidNodes++;
					rigidNodeLocalMeshSize += Element[pn].FastGetSolutionStepValue(NODAL_H_WALL);
					rigidNodeMeshCounter += 1.0;
				}
				if (Element[pn].Is(TO_ERASE))
				{
					toEraseNodeFound = true;
				}
				if (Element[pn].Is(FREE_SURFACE))
				{
					freesurfaceNodes++;
				}
			}

			if (rigidNodeMeshCounter > 0)
			{
				const double rigidWallMeshSize = rigidNodeLocalMeshSize / rigidNodeMeshCounter;
				const double ratio = rigidWallMeshSize / meanMeshSize;
				const double tolerance = 1.8;
				if (ratio > tolerance)
				{
					meanMeshSize *= 0.5;
					meanMeshSize += 0.5 * rigidWallMeshSize;
				}
			}
			double penalization = 1.0; // penalization here should be greater than 1

			if (freesurfaceNodes > 0)
			{
				penalization = 1.2; // to avoid to gain too much volume during remeshing step
			}
			array_1d<double, 3> Edges(3, 0.0);
			array_1d<SizeType, 3> FirstEdgeNode(3, 0);
			array_1d<SizeType, 3> SecondEdgeNode(3, 0);
			double WallCharacteristicDistance = 0;
			ComputeWallCharacteristicDistance2DWithRefinement(Element, WallCharacteristicDistance, Edges, FirstEdgeNode, SecondEdgeNode);

			bool dangerousElement = false;
			DetectDangerousElements2DWithRefinement(Element, WallCharacteristicDistance, Edges, FirstEdgeNode, SecondEdgeNode, rigidNodes, penalization, dangerousElement);

			const double limitEdgeLength = 1.9 * meanMeshSize * penalization;
			const double extraLimitEdgeLength = 2.5 * meanMeshSize * penalization;

			if (dangerousElement == false && toEraseNodeFound == false)
			{
				SizeType maxCount = 3;
				double LargestEdge = 0;

				for (SizeType i = 0; i < 3; i++)
				{
					if (Edges[i] > LargestEdge)
					{
						maxCount = i;
						LargestEdge = Edges[i];
					}
				}

				if (((CountNodes < (ElementsToRefine + nodesInTransitionZone) || insideTransitionZone == true) && LargestEdge > limitEdgeLength) || LargestEdge > extraLimitEdgeLength)
				{
					bool newNode = true;
					for (SizeType i = 0; i < unsigned(CountNodes); i++)
					{
						if ((nodes_id_to_interpolate[i][0] == Element[FirstEdgeNode[maxCount]].GetId() && nodes_id_to_interpolate[i][1] == Element[SecondEdgeNode[maxCount]].GetId()) ||
							(nodes_id_to_interpolate[i][1] == Element[FirstEdgeNode[maxCount]].GetId() && nodes_id_to_interpolate[i][0] == Element[SecondEdgeNode[maxCount]].GetId()))
						{
							newNode = false;
						}
					}
					if (newNode == true)
					{
						new_positions.resize(CountNodes + 1);
						nodes_id_to_interpolate.resize(CountNodes + 1);
						array_1d<double, 3> new_position = (Element[FirstEdgeNode[maxCount]].Coordinates() + Element[SecondEdgeNode[maxCount]].Coordinates()) * 0.5;
						nodes_id_to_interpolate[CountNodes][0] = Element[FirstEdgeNode[maxCount]].GetId();
						nodes_id_to_interpolate[CountNodes][1] = Element[SecondEdgeNode[maxCount]].GetId();
						new_positions[CountNodes] = new_position;
						CountNodes++;
						if (insideTransitionZone)
						{
							nodesInTransitionZone++;
						}
					}
				}
			}

			KRATOS_CATCH("")
		}

		void SelectEdgeToRefine3DWithRefinement(Element::GeometryType &Element,
												std::vector<array_1d<double, 3>> &new_positions,
												std::vector<array_1d<SizeType, 4>> &nodes_id_to_interpolate,
												int &CountNodes,
												const int ElementsToRefine,
												int &nodesInTransitionZone)
		{
			KRATOS_TRY

			const SizeType nds = Element.size();

			SizeType rigidNodes = 0;
			SizeType freesurfaceNodes = 0;
			bool toEraseNodeFound = false;

			double meanMeshSize = mrRemesh.Refine->CriticalRadius;
			const ProcessInfo &rCurrentProcessInfo = mrModelPart.GetProcessInfo();
			double currentTime = rCurrentProcessInfo[TIME];
			bool insideTransitionZone = false;
			double rigidNodeLocalMeshSize = 0;
			double rigidNodeMeshCounter = 0;
			for (SizeType pn = 0; pn < nds; ++pn)
			{
				mMesherUtilities.DefineMeshSizeInTransitionZones3D(mrRemesh, currentTime, Element[pn].Coordinates(), meanMeshSize, insideTransitionZone);
			}

			for (SizeType pn = 0; pn < nds; ++pn)
			{
				if (Element[pn].Is(RIGID))
				{
					rigidNodes++;
					rigidNodeLocalMeshSize += Element[pn].FastGetSolutionStepValue(NODAL_H_WALL);
					rigidNodeMeshCounter += 1.0;
				}
				if (Element[pn].Is(TO_ERASE))
				{
					toEraseNodeFound = true;
				}
				if (Element[pn].Is(FREE_SURFACE))
				{
					freesurfaceNodes++;
				}
			}

			if (rigidNodeMeshCounter > 0)
			{
				const double rigidWallMeshSize = rigidNodeLocalMeshSize / rigidNodeMeshCounter;
				const double ratio = rigidWallMeshSize / meanMeshSize;
				const double tolerance = 1.8;
				if (ratio > tolerance)
				{
					meanMeshSize *= 0.5;
					meanMeshSize += 0.5 * rigidWallMeshSize;
				}
			}

			double penalization = 1.0; // penalization here should be greater than 1

			if (freesurfaceNodes > 0)
			{
				penalization = 1.2; // to avoid to gain too much volume during remeshing step
			}

			array_1d<double, 6> Edges(6, 0.0);
			array_1d<SizeType, 6> FirstEdgeNode(6, 0);
			array_1d<SizeType, 6> SecondEdgeNode(6, 0);
			double WallCharacteristicDistance = 0;
			ComputeWallCharacteristicDistance3DWithRefinement(Element, WallCharacteristicDistance, Edges, FirstEdgeNode, SecondEdgeNode);

			// Edges connectivity: Edges[0]=d01, Edges[1]=d20, Edges[2]=d21, Edges[3]=d30, Edges[4]=d31, Edges[5]=d32
			bool dangerousElement = false;
			DetectDangerousElements3DWithRefinement(Element, WallCharacteristicDistance, Edges, FirstEdgeNode, SecondEdgeNode, rigidNodes, penalization, dangerousElement);

			const double limitEdgeLength = 1.9 * meanMeshSize * penalization;
			const double extraLimitEdgeLength = 2.5 * meanMeshSize * penalization;

			// just to fill the vector
			if (dangerousElement == false && toEraseNodeFound == false)
			{
				SizeType maxCount = 6;
				double LargestEdge = 0;
				for (SizeType i = 0; i < 6; i++)
				{
					if (Edges[i] > LargestEdge)
					{
						maxCount = i;
						LargestEdge = Edges[i];
					}
				}

				if (((CountNodes < (ElementsToRefine + nodesInTransitionZone) || insideTransitionZone == true) && LargestEdge > limitEdgeLength) || LargestEdge > extraLimitEdgeLength)
				{
					bool newNode = true;
					for (SizeType i = 0; i < unsigned(CountNodes); i++)
					{
						if ((nodes_id_to_interpolate[i][0] == Element[FirstEdgeNode[maxCount]].GetId() && nodes_id_to_interpolate[i][1] == Element[SecondEdgeNode[maxCount]].GetId()) ||
							(nodes_id_to_interpolate[i][1] == Element[FirstEdgeNode[maxCount]].GetId() && nodes_id_to_interpolate[i][0] == Element[SecondEdgeNode[maxCount]].GetId()))
						{
							newNode = false;
						}
					}
					if (newNode == true)
					{
						new_positions.resize(CountNodes + 1);
						nodes_id_to_interpolate.resize(CountNodes + 1);
						array_1d<double, 3> new_position = (Element[FirstEdgeNode[maxCount]].Coordinates() + Element[SecondEdgeNode[maxCount]].Coordinates()) * 0.5;
						nodes_id_to_interpolate[CountNodes][0] = Element[FirstEdgeNode[maxCount]].GetId();
						nodes_id_to_interpolate[CountNodes][1] = Element[SecondEdgeNode[maxCount]].GetId();
						new_positions[CountNodes] = new_position;
						CountNodes++;
						if (insideTransitionZone)
						{
							nodesInTransitionZone++;
						}
					}
				}
			}

			KRATOS_CATCH("")
		}

		void CreateAndAddNewNodesInCornerWall(std::vector<array_1d<double, 3>> &new_positions,
											  std::vector<array_1d<SizeType, 4>> &nodes_id_to_interpolate,
											  std::vector<Node::DofsContainerType> &NewDofs,
											  int ElementsToRefine,
											  SizeType &maxId)
		{
			KRATOS_TRY

			const SizeType dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

			std::vector<Node::Pointer> list_of_new_nodes;

			// assign data to dofs
			VariablesList &VariablesList = mrModelPart.GetNodalSolutionStepVariablesList();

			for (SizeType nn = 0; nn < new_positions.size(); nn++)
			{

				SizeType id = maxId + 1 + nn;

				double x = new_positions[nn][0];
				double y = new_positions[nn][1];
				double z = 0;
				if (dimension == 3)
					z = new_positions[nn][2];

				Node::Pointer pnode = mrModelPart.CreateNewNode(id, x, y, z);
				pnode->Set(NEW_ENTITY); // not boundary
				list_of_new_nodes.push_back(pnode);
				if (mrRemesh.InputInitializedFlag)
				{
					mrRemesh.NodalPreIds.push_back(pnode->Id());
					pnode->SetId(id);
				}

				// //giving model part variables list to the node
				pnode->SetSolutionStepVariablesList(&VariablesList);

				// //set buffer size
				pnode->SetBufferSize(mrModelPart.GetBufferSize());

				Node::DofsContainerType &reference_dofs = NewDofs[nn];

				for (Node::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
				{
					Node::DofType &rDof = **iii;
					pnode->pAddDof(rDof);
				}

				auto slave_node_1 = mrModelPart.pGetNode(nodes_id_to_interpolate[nn][0]);
				auto slave_node_2 = mrModelPart.pGetNode(nodes_id_to_interpolate[nn][1]);
				auto slave_node_3 = mrModelPart.pGetNode(nodes_id_to_interpolate[nn][2]);

				InterpolateFromTwoNodes(pnode, slave_node_1, slave_node_2, VariablesList);

				TakeMaterialPropertiesFromNotRigidNode(pnode, slave_node_3);
			}

			// set the coordinates to the original value
			const array_1d<double, 3> ZeroNormal(3, 0.0);
			for (std::vector<Node::Pointer>::iterator it = list_of_new_nodes.begin(); it != list_of_new_nodes.end(); it++)
			{
				const array_1d<double, 3> &displacement = (*it)->FastGetSolutionStepValue(DISPLACEMENT);
				(*it)->X0() = (*it)->X() - displacement[0];
				(*it)->Y0() = (*it)->Y() - displacement[1];
				(*it)->Z0() = (*it)->Z() - displacement[2];

				(*it)->Set(FLUID);
				(*it)->Set(ACTIVE);
				(*it)->Reset(TO_ERASE);
			}

			KRATOS_CATCH("")
		}

		void CreateAndAddNewNodes(std::vector<array_1d<double, 3>> &new_positions,
								  std::vector<array_1d<SizeType, 4>> &nodes_id_to_interpolate,
								  int ElementsToRefine,
								  SizeType &maxId)
		{
			KRATOS_TRY

			const SizeType dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

			std::vector<Node::Pointer> list_of_new_nodes;
			double NodeIdParent = MesherUtilities::GetMaxNodeId(mrModelPart.GetParentModelPart());
			double NodeId = MesherUtilities::GetMaxNodeId(mrModelPart);

			SizeType initial_node_size = NodeIdParent + 1 + ElementsToRefine; // total model part node size

			NodeType::Pointer pnode;
			NodeType::DofsContainerType &ReferenceDofs = mrModelPart.Nodes().front().GetDofs();

			if (NodeId > NodeIdParent)
			{
				initial_node_size = NodeId + 1 + ElementsToRefine;
				std::cout << "initial_node_size  " << initial_node_size << std::endl;
			}

			// assign data to dofs
			VariablesList &VariablesList = mrModelPart.GetNodalSolutionStepVariablesList();

			const ProcessInfo &rCurrentProcessInfo = mrModelPart.GetProcessInfo();
			SizeType principalModelPartId = rCurrentProcessInfo[MAIN_MATERIAL_PROPERTY];

			for (SizeType nn = 0; nn < new_positions.size(); nn++)
			{

				SizeType id = initial_node_size + nn;
				maxId = id;
				double x = new_positions[nn][0];
				double y = new_positions[nn][1];
				double z = 0;
				if (dimension == 3)
					z = new_positions[nn][2];

				// Node::Pointer pnode = mrModelPart.CreateNewNode(id, x, y, z);
				pnode = Kratos::make_intrusive<Node>(id, x, y, z);

				pnode->Set(NEW_ENTITY); // not boundary

				if (mrRemesh.InputInitializedFlag)
				{
					mrRemesh.NodalPreIds.push_back(pnode->Id());
					pnode->SetId(id);
				}

				// //giving model part variables list to the node
				pnode->SetSolutionStepVariablesList(&VariablesList);

				// //set buffer size
				pnode->SetBufferSize(mrModelPart.GetBufferSize());

				// Node::DofsContainerType &reference_dofs = NewDofs[nn];
				// for (Node::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
				// {
				// 	Node::DofType &rDof = **iii;
				// 	pnode->pAddDof(rDof);
				// }
				// generating the dofs
				for (Node::DofsContainerType::iterator i_dof = ReferenceDofs.begin(); i_dof != ReferenceDofs.end(); ++i_dof)
				{
					NodeType::DofType &rDof = **i_dof;
					NodeType::DofType::Pointer pNewDof = pnode->pAddDof(rDof);

					(pNewDof)->FreeDof();
				}

				list_of_new_nodes.push_back(pnode);

				auto slave_node_1 = mrModelPart.pGetNode(nodes_id_to_interpolate[nn][0]);
				auto slave_node_2 = mrModelPart.pGetNode(nodes_id_to_interpolate[nn][1]);
				InterpolateFromTwoNodes(pnode, slave_node_1, slave_node_2, VariablesList);
				if (slave_node_1->Is(RIGID) || slave_node_1->Is(SOLID))
				{
					TakeMaterialPropertiesFromNotRigidNode(pnode, slave_node_2);
				}
				else
				{

					SizeType propertyIdNodeSlave1 = slave_node_1->FastGetSolutionStepValue(PROPERTY_ID);
					if ((mrRemesh.Info->BalancePrincipalSecondaryPartsNodes > 0 && propertyIdNodeSlave1 == principalModelPartId) ||
						(mrRemesh.Info->BalancePrincipalSecondaryPartsNodes < 0 && propertyIdNodeSlave1 != principalModelPartId) ||
						(slave_node_2->Is(RIGID) || slave_node_2->Is(SOLID)))
					{
						TakeMaterialPropertiesFromNotRigidNode(pnode, slave_node_1);
					}
					else
					{
						TakeMaterialPropertiesFromNotRigidNode(pnode, slave_node_2);
					}
				}

				SizeType propertyIdNode = pnode->FastGetSolutionStepValue(PROPERTY_ID);
				if (propertyIdNode != principalModelPartId)
				{
					mrRemesh.Info->BalancePrincipalSecondaryPartsNodes += 1;
				}
			}

			// set the coordinates to the original value
			const array_1d<double, 3> ZeroNormal(3, 0.0);
			for (std::vector<Node::Pointer>::iterator it = list_of_new_nodes.begin(); it != list_of_new_nodes.end(); it++)
			{
				const array_1d<double, 3> &displacement = (*it)->FastGetSolutionStepValue(DISPLACEMENT);
				(*it)->X0() = (*it)->X() - displacement[0];
				(*it)->Y0() = (*it)->Y() - displacement[1];
				(*it)->Z0() = (*it)->Z() - displacement[2];

				(*it)->Set(FLUID);
				(*it)->Set(ACTIVE);
				(*it)->Reset(TO_ERASE);

				mrModelPart.Nodes().push_back(*(it));
			}

			KRATOS_CATCH("")
		}

		void InterpolateFromTwoNodes(Node::Pointer master_node, Node::Pointer slave_node_1, Node::Pointer slave_node_2, VariablesList &rVariablesList)
		{

			KRATOS_TRY

			SizeType buffer_size = master_node->GetBufferSize();

			for (VariablesList::const_iterator i_variable = rVariablesList.begin(); i_variable != rVariablesList.end(); i_variable++)
			{
				std::string variable_name = i_variable->Name();
				if (KratosComponents<Variable<double>>::Has(variable_name))
				{
					const Variable<double> &variable = KratosComponents<Variable<double>>::Get(variable_name);
					for (SizeType step = 0; step < buffer_size; step++)
					{
						// getting the data of the solution step
						double &node_data = master_node->FastGetSolutionStepValue(variable, step);

						double node0_data = slave_node_1->FastGetSolutionStepValue(variable, step);
						double node1_data = slave_node_2->FastGetSolutionStepValue(variable, step);

						node_data = (0.5 * node0_data + 0.5 * node1_data);
					}
				}
				else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name))
				{
					const Variable<array_1d<double, 3>> &variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(variable_name);
					for (SizeType step = 0; step < buffer_size; step++)
					{
						// getting the data of the solution step
						array_1d<double, 3> &node_data = master_node->FastGetSolutionStepValue(variable, step);

						const array_1d<double, 3> &node0_data = slave_node_1->FastGetSolutionStepValue(variable, step);
						const array_1d<double, 3> &node1_data = slave_node_2->FastGetSolutionStepValue(variable, step);

						noalias(node_data) = (0.5 * node0_data + 0.5 * node1_data);
						// node_data = (0.5*node0_data + 0.5*node1_data);
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
						Matrix &node_data = master_node->FastGetSolutionStepValue(variable, step);

						Matrix &node0_data = slave_node_1->FastGetSolutionStepValue(variable, step);
						Matrix &node1_data = slave_node_2->FastGetSolutionStepValue(variable, step);

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
					// std::cout<<"Vector"<<std::endl;
					const Variable<Vector> &variable = KratosComponents<Variable<Vector>>::Get(variable_name);
					for (SizeType step = 0; step < buffer_size; step++)
					{
						// getting the data of the solution step
						Vector &node_data = master_node->FastGetSolutionStepValue(variable, step);

						Vector &node0_data = slave_node_1->FastGetSolutionStepValue(variable, step);
						Vector &node1_data = slave_node_2->FastGetSolutionStepValue(variable, step);

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

		void TakeMaterialPropertiesFromNotRigidNode(Node::Pointer master_node, Node::Pointer SlaveNode)
		{

			KRATOS_TRY
			master_node->FastGetSolutionStepValue(PROPERTY_ID) = SlaveNode->FastGetSolutionStepValue(PROPERTY_ID);
			if (master_node->SolutionStepsDataHas(BULK_MODULUS) && SlaveNode->SolutionStepsDataHas(BULK_MODULUS))
				master_node->FastGetSolutionStepValue(BULK_MODULUS) = SlaveNode->FastGetSolutionStepValue(BULK_MODULUS);
			if (master_node->SolutionStepsDataHas(DENSITY) && SlaveNode->SolutionStepsDataHas(DENSITY))
				master_node->FastGetSolutionStepValue(DENSITY) = SlaveNode->FastGetSolutionStepValue(DENSITY);
			if (master_node->SolutionStepsDataHas(DYNAMIC_VISCOSITY) && SlaveNode->SolutionStepsDataHas(DYNAMIC_VISCOSITY))
				master_node->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = SlaveNode->FastGetSolutionStepValue(DYNAMIC_VISCOSITY);

			if (master_node->SolutionStepsDataHas(YIELD_SHEAR) && SlaveNode->SolutionStepsDataHas(YIELD_SHEAR))
				master_node->FastGetSolutionStepValue(YIELD_SHEAR) = SlaveNode->FastGetSolutionStepValue(YIELD_SHEAR);

			if (master_node->SolutionStepsDataHas(FLOW_INDEX) && SlaveNode->SolutionStepsDataHas(FLOW_INDEX))
				master_node->FastGetSolutionStepValue(FLOW_INDEX) = SlaveNode->FastGetSolutionStepValue(FLOW_INDEX);
			if (master_node->SolutionStepsDataHas(ADAPTIVE_EXPONENT) && SlaveNode->SolutionStepsDataHas(ADAPTIVE_EXPONENT))
				master_node->FastGetSolutionStepValue(ADAPTIVE_EXPONENT) = SlaveNode->FastGetSolutionStepValue(ADAPTIVE_EXPONENT);
			if (master_node->SolutionStepsDataHas(STATIC_FRICTION) && SlaveNode->SolutionStepsDataHas(STATIC_FRICTION))
				master_node->FastGetSolutionStepValue(STATIC_FRICTION) = SlaveNode->FastGetSolutionStepValue(STATIC_FRICTION);
			if (master_node->SolutionStepsDataHas(DYNAMIC_FRICTION) && SlaveNode->SolutionStepsDataHas(DYNAMIC_FRICTION))
				master_node->FastGetSolutionStepValue(DYNAMIC_FRICTION) = SlaveNode->FastGetSolutionStepValue(DYNAMIC_FRICTION);
			if (master_node->SolutionStepsDataHas(INERTIAL_NUMBER_ZERO) && SlaveNode->SolutionStepsDataHas(INERTIAL_NUMBER_ZERO))
				master_node->FastGetSolutionStepValue(INERTIAL_NUMBER_ZERO) = SlaveNode->FastGetSolutionStepValue(INERTIAL_NUMBER_ZERO);
			if (master_node->SolutionStepsDataHas(GRAIN_DIAMETER) && SlaveNode->SolutionStepsDataHas(GRAIN_DIAMETER))
				master_node->FastGetSolutionStepValue(GRAIN_DIAMETER) = SlaveNode->FastGetSolutionStepValue(GRAIN_DIAMETER);
			if (master_node->SolutionStepsDataHas(GRAIN_DENSITY) && SlaveNode->SolutionStepsDataHas(GRAIN_DENSITY))
				master_node->FastGetSolutionStepValue(GRAIN_DENSITY) = SlaveNode->FastGetSolutionStepValue(GRAIN_DENSITY);
			if (master_node->SolutionStepsDataHas(REGULARIZATION_COEFFICIENT) && SlaveNode->SolutionStepsDataHas(REGULARIZATION_COEFFICIENT))
				master_node->FastGetSolutionStepValue(REGULARIZATION_COEFFICIENT) = SlaveNode->FastGetSolutionStepValue(REGULARIZATION_COEFFICIENT);

			if (master_node->SolutionStepsDataHas(INTERNAL_FRICTION_ANGLE) && SlaveNode->SolutionStepsDataHas(INTERNAL_FRICTION_ANGLE))
				master_node->FastGetSolutionStepValue(INTERNAL_FRICTION_ANGLE) = SlaveNode->FastGetSolutionStepValue(INTERNAL_FRICTION_ANGLE);
			if (master_node->SolutionStepsDataHas(COHESION) && SlaveNode->SolutionStepsDataHas(COHESION))
				master_node->FastGetSolutionStepValue(COHESION) = SlaveNode->FastGetSolutionStepValue(COHESION);

			if (master_node->SolutionStepsDataHas(DEVIATORIC_COEFFICIENT) && SlaveNode->SolutionStepsDataHas(DEVIATORIC_COEFFICIENT))
			{
				master_node->FastGetSolutionStepValue(DEVIATORIC_COEFFICIENT) = SlaveNode->FastGetSolutionStepValue(DEVIATORIC_COEFFICIENT);
				master_node->FastGetSolutionStepValue(VOLUMETRIC_COEFFICIENT) = SlaveNode->FastGetSolutionStepValue(VOLUMETRIC_COEFFICIENT);
			}

			if (master_node->SolutionStepsDataHas(YOUNG_MODULUS) && SlaveNode->SolutionStepsDataHas(YOUNG_MODULUS))
			{
				master_node->FastGetSolutionStepValue(YOUNG_MODULUS) = 0;
				master_node->FastGetSolutionStepValue(POISSON_RATIO) = 0;
			}
			if (master_node->SolutionStepsDataHas(SOLID_DENSITY) && SlaveNode->SolutionStepsDataHas(SOLID_DENSITY))
			{
				master_node->FastGetSolutionStepValue(SOLID_DENSITY) = 0;
				master_node->Reset(SOLID);
				master_node->FastGetSolutionStepValue(INTERFACE_NODE) = false;
			}

			KRATOS_CATCH("")
		}

		///@}
		///@name Private  Access
		///@{

		///@}
		///@name Private Inquiry
		///@{

		///@}
		///@name Un accessible methods
		///@{

		/// Assignment operator.
		GenerateNewNodesBeforeMeshingProcess &operator=(GenerateNewNodesBeforeMeshingProcess const &rOther);

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
									GenerateNewNodesBeforeMeshingProcess &rThis);

	/// output stream function
	inline std::ostream &operator<<(std::ostream &rOStream,
									const GenerateNewNodesBeforeMeshingProcess &rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@}

} // namespace Kratos.

#endif // KRATOS_GENERATE_NEW_NODES_BEFORE_MESHING_PROCESS_H_INCLUDED  defined
