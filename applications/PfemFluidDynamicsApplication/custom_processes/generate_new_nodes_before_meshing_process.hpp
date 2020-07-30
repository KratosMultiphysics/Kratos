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

///VARIABLES used:
//Data:
//StepData: CONTACT_FORCE, DISPLACEMENT
//Flags:    (checked)
//          (set)
//          (modified)
//          (reset)
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

		double initialTimeRefiningBox = mrRemesh.RefiningBoxInitialTime;
		double finalTimeRefiningBox = mrRemesh.RefiningBoxFinalTime;
		bool refiningBox = mrRemesh.UseRefiningBox;

		const unsigned int dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

		if (!(refiningBox == true && currentTime > initialTimeRefiningBox && currentTime < finalTimeRefiningBox))
		{
			refiningBox = false;
		}

		if (currentTime < 2 * timeInterval)
		{
			mrRemesh.Info->RemovedNodes = 0;
			if (mEchoLevel > 1)
				std::cout << " First meshes: I repare the mesh without adding new nodes" << std::endl;
			mrRemesh.Info->InitialNumberOfNodes = mrRemesh.Info->NumberOfNodes;
		}

		int ElementsToRefine = mrRemesh.Info->RemovedNodes;

		int initialNumberOfNodes = mrRemesh.Info->InitialNumberOfNodes;
		int numberOfNodes = mrRemesh.Info->NumberOfNodes;
		int extraNodes = numberOfNodes - initialNumberOfNodes;
		int toleredExtraNodes = int(0.05 * mrRemesh.Info->InitialNumberOfNodes);

		if (mrRemesh.ExecutionOptions.Is(MesherUtilities::REFINE_WALL_CORNER))
		{
			if ((extraNodes + ElementsToRefine) > toleredExtraNodes && refiningBox == false)
			{
				ElementsToRefine = toleredExtraNodes - extraNodes;
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

			if (ElementsToRefine > 0)
			{
				std::vector<array_1d<double, 3>> NewPositions;
				std::vector<double> BiggestVolumes;
				std::vector<array_1d<unsigned int, 4>> NodesIDToInterpolate;
				std::vector<Node<3>::DofsContainerType> NewDofs;

				int CountNodes = 0;

				NewPositions.resize(ElementsToRefine);
				BiggestVolumes.resize(ElementsToRefine);
				NodesIDToInterpolate.resize(ElementsToRefine);
				NewDofs.resize(ElementsToRefine);

				for (int nn = 0; nn < ElementsToRefine; nn++)
				{
					BiggestVolumes[nn] = -1.0;
				}

				std::vector<array_1d<double, 3>> CornerWallNewPositions;
				std::vector<array_1d<unsigned int, 4>> CornerWallNodesIDToInterpolate;
				std::vector<Node<3>::DofsContainerType> CornerWallNewDofs;
				int cornerWallNewNodes = 0;
				int maxOfNewWallNodes = toleredExtraNodes;
				if (mrRemesh.ExecutionOptions.Is(MesherUtilities::REFINE_WALL_CORNER))
				{
					CornerWallNewPositions.resize(maxOfNewWallNodes);
					CornerWallNodesIDToInterpolate.resize(maxOfNewWallNodes);
					CornerWallNewDofs.resize(maxOfNewWallNodes);
				}

				ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin();
				// const unsigned int nds = element_begin->GetGeometry().size();
				for (ModelPart::ElementsContainerType::const_iterator ie = element_begin; ie != mrModelPart.ElementsEnd(); ie++)
				{

					//////// choose the right (big and safe) elements to refine and compute the new node position and variables ////////
					if (dimension == 2)
					{
						SelectEdgeToRefine2D(ie->GetGeometry(), NewPositions, BiggestVolumes, NodesIDToInterpolate, NewDofs, CountNodes, ElementsToRefine);

						if (mrRemesh.ExecutionOptions.Is(MesherUtilities::REFINE_WALL_CORNER) && cornerWallNewNodes < maxOfNewWallNodes)
						{
							InsertNodeInCornerElement2D(ie->GetGeometry(), CornerWallNewPositions, CornerWallNodesIDToInterpolate, CornerWallNewDofs, cornerWallNewNodes);
						}
					}
					else if (dimension == 3)
					{
						SelectEdgeToRefine3D(ie->GetGeometry(), NewPositions, BiggestVolumes, NodesIDToInterpolate, NewDofs, CountNodes, ElementsToRefine);

						if (mrRemesh.ExecutionOptions.Is(MesherUtilities::REFINE_WALL_CORNER) && cornerWallNewNodes < maxOfNewWallNodes)
						{
							InsertNodeInCornerElement3D(ie->GetGeometry(), CornerWallNewPositions, CornerWallNodesIDToInterpolate, CornerWallNewDofs, cornerWallNewNodes);
						}
					}

				} // elements loop

				mrRemesh.Info->RemovedNodes -= ElementsToRefine;
				if (CountNodes < ElementsToRefine)
				{
					mrRemesh.Info->RemovedNodes += ElementsToRefine - CountNodes;
					NewPositions.resize(CountNodes);
					BiggestVolumes.resize(CountNodes);
					NodesIDToInterpolate.resize(CountNodes);
					NewDofs.resize(CountNodes);
				}
				unsigned int maxId = 0;
				CreateAndAddNewNodes(NewPositions, NodesIDToInterpolate, NewDofs, ElementsToRefine, maxId);

				if (mrRemesh.ExecutionOptions.Is(MesherUtilities::REFINE_WALL_CORNER))
				{
					if (cornerWallNewNodes < maxOfNewWallNodes)
					{
						CornerWallNewPositions.resize(cornerWallNewNodes);
						CornerWallNewDofs.resize(cornerWallNewNodes);
						CornerWallNodesIDToInterpolate.resize(cornerWallNewNodes);
					}
					CreateAndAddNewNodesInCornerWall(CornerWallNewPositions, CornerWallNodesIDToInterpolate, CornerWallNewDofs, cornerWallNewNodes, maxId);
				}
			}
		}
		else
		{

			if (ElementsToRefine < 10)
			{
				ElementsToRefine = 10;
			}
			else
			{
				ElementsToRefine *= 3;
			}
			std::vector<array_1d<double, 3>> NewPositions;
			std::vector<double> BiggestVolumes;
			std::vector<array_1d<unsigned int, 4>> NodesIDToInterpolate;
			std::vector<Node<3>::DofsContainerType> NewDofs;

			int CountNodes = 0;

			NewPositions.resize(ElementsToRefine);
			BiggestVolumes.resize(ElementsToRefine);
			NodesIDToInterpolate.resize(ElementsToRefine);
			NewDofs.resize(ElementsToRefine);

			for (int nn = 0; nn < ElementsToRefine; nn++)
			{
				BiggestVolumes[nn] = -1.0;
			}

			ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin();
			// const unsigned int nds = element_begin->GetGeometry().size();
			for (ModelPart::ElementsContainerType::const_iterator ie = element_begin; ie != mrModelPart.ElementsEnd(); ie++)
			{

				const unsigned int dimension = ie->GetGeometry().WorkingSpaceDimension();

				//////// choose the right (big and safe) elements to refine and compute the new node position and variables ////////
				if (dimension == 2)
				{
					SelectEdgeToRefine2DWithRefinement(ie->GetGeometry(), NewPositions, BiggestVolumes, NodesIDToInterpolate, NewDofs, CountNodes, ElementsToRefine);
				}
				else if (dimension == 3)
				{
					SelectEdgeToRefine3DWithRefinement(ie->GetGeometry(), NewPositions, BiggestVolumes, NodesIDToInterpolate, NewDofs, CountNodes, ElementsToRefine);
				}

			} // elements loop

			mrRemesh.Info->RemovedNodes -= ElementsToRefine;
			if (CountNodes < ElementsToRefine)
			{
				mrRemesh.Info->RemovedNodes += ElementsToRefine - CountNodes;
				NewPositions.resize(CountNodes);
				BiggestVolumes.resize(CountNodes);
				NodesIDToInterpolate.resize(CountNodes);
				NewDofs.resize(CountNodes);
			}
			unsigned int maxId = 0;
			CreateAndAddNewNodes(NewPositions, NodesIDToInterpolate, NewDofs, ElementsToRefine, maxId);
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
	///@name Private Operators
	///@{

	///@}
	///@name Private Operations
	///@{

	void CopyDofs(Node<3>::DofsContainerType const &From, Node<3>::DofsContainerType &To)
	{
		for (auto &p_dof : From)
		{
			To.push_back(Kratos::unique_ptr<Dof<double>>(new Dof<double>(*p_dof)));
		}
	}

	void InsertNodeInCornerElement2D(Element::GeometryType &Element,
									 std::vector<array_1d<double, 3>> &NewPositions,
									 std::vector<array_1d<unsigned int, 4>> &NodesIDToInterpolate,
									 std::vector<Node<3>::DofsContainerType> &NewDofs,
									 int &CountNodes)
	{
		KRATOS_TRY

		const unsigned int nds = Element.size();

		unsigned int rigidNodes = 0;
		unsigned int freesurfaceNodes = 0;

		for (unsigned int pn = 0; pn < nds; pn++)
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

			if (Element[0].Is(RIGID) && Element[1].Is(RIGID))
			{
				NormalA = Element[0].FastGetSolutionStepValue(NORMAL);
				NormalB = Element[1].FastGetSolutionStepValue(NORMAL);
				cosAngle = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1];
				if (cosAngle < cosTolerance && cosAngle > -cosTolerance)
				{
					array_1d<double, 3> NewPosition = (Element[0].Coordinates() + Element[1].Coordinates()) * 0.5;
					NodesIDToInterpolate[CountNodes][0] = Element[0].GetId();
					NodesIDToInterpolate[CountNodes][1] = Element[1].GetId();
					if (Element[2].IsNot(TO_ERASE))
					{
						NodesIDToInterpolate[CountNodes][2] = Element[2].GetId();
					}
					else
					{
						NodesIDToInterpolate[CountNodes][2] = Element[0].GetId();
					}
					CopyDofs(Element[2].GetDofs(), NewDofs[CountNodes]);
					CopyDofs(Element[2].GetDofs(), NewDofs[CountNodes]);
					NewPositions[CountNodes] = NewPosition;
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
					array_1d<double, 3> NewPosition = (Element[0].Coordinates() + Element[1].Coordinates()) * 0.5;
					NodesIDToInterpolate[CountNodes][0] = Element[0].GetId();
					NodesIDToInterpolate[CountNodes][1] = Element[2].GetId();
					if (Element[1].IsNot(TO_ERASE))
					{
						NodesIDToInterpolate[CountNodes][2] = Element[1].GetId();
					}
					else
					{
						NodesIDToInterpolate[CountNodes][2] = Element[0].GetId();
					}
					CopyDofs(Element[1].GetDofs(), NewDofs[CountNodes]);
					CopyDofs(Element[1].GetDofs(), NewDofs[CountNodes]);
					NewPositions[CountNodes] = NewPosition;
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
					array_1d<double, 3> NewPosition = (Element[2].Coordinates() + Element[1].Coordinates()) * 0.5;
					NodesIDToInterpolate[CountNodes][0] = Element[2].GetId();
					NodesIDToInterpolate[CountNodes][1] = Element[1].GetId();
					if (Element[0].IsNot(TO_ERASE))
					{
						NodesIDToInterpolate[CountNodes][2] = Element[0].GetId();
					}
					else
					{
						NodesIDToInterpolate[CountNodes][2] = Element[2].GetId();
					}
					CopyDofs(Element[0].GetDofs(), NewDofs[CountNodes]);
					CopyDofs(Element[0].GetDofs(), NewDofs[CountNodes]);
					NewPositions[CountNodes] = NewPosition;
					CountNodes++;
				}
			}
		}

		KRATOS_CATCH("")
	}

	void InsertNodeInCornerElement3D(Element::GeometryType &Element,
									 std::vector<array_1d<double, 3>> &NewPositions,
									 std::vector<array_1d<unsigned int, 4>> &NodesIDToInterpolate,
									 std::vector<Node<3>::DofsContainerType> &NewDofs,
									 int &CountNodes)
	{
		KRATOS_TRY

		const unsigned int nds = Element.size();

		unsigned int rigidNodes = 0;
		unsigned int freesurfaceNodes = 0;
		unsigned int toEraseNodes = 0;

		for (unsigned int pn = 0; pn < nds; pn++)
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
			double normNormalA=0;
			double normNormalB=0;
			double cos = 1.0;
			double minCos = 1.0;
			array_1d<unsigned int, 2> idsWallNodes(2, 0);
			unsigned int idFreeNode = 0;
			double cosTolerance = 0.1;
			if (Element[0].IsNot(RIGID))
			{
				NormalA = Element[1].FastGetSolutionStepValue(NORMAL);
				NormalB = Element[2].FastGetSolutionStepValue(NORMAL);
				normNormalA=NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
				normNormalB=NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
				cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
				if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA>0.99 && normNormalA<1.01) && (normNormalB>0.99 && normNormalB<1.01))
				{
					minCos = cos;
					idsWallNodes[0] = 1;
					idsWallNodes[1] = 2;
					idFreeNode = 0;
				}

				NormalA = Element[1].FastGetSolutionStepValue(NORMAL);
				NormalB = Element[3].FastGetSolutionStepValue(NORMAL);
				normNormalA=NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
				normNormalB=NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
				cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
				if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA>0.99 && normNormalA<1.01) && (normNormalB>0.99 && normNormalB<1.01))
				{
					minCos = cos;
					idsWallNodes[0] = 1;
					idsWallNodes[1] = 3;
					idFreeNode = 0;
				}

				NormalA = Element[2].FastGetSolutionStepValue(NORMAL);
				NormalB = Element[3].FastGetSolutionStepValue(NORMAL);
				normNormalA=NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
				normNormalB=NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
				cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
				if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA>0.99 && normNormalA<1.01) && (normNormalB>0.99 && normNormalB<1.01))
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
				normNormalA=NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
				normNormalB=NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
				cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
				if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA>0.99 && normNormalA<1.01) && (normNormalB>0.99 && normNormalB<1.01))
				{
					minCos = cos;
					idsWallNodes[0] = 0;
					idsWallNodes[1] = 2;
					idFreeNode = 1;
				}

				NormalA = Element[0].FastGetSolutionStepValue(NORMAL);
				NormalB = Element[3].FastGetSolutionStepValue(NORMAL);
				normNormalA=NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
				normNormalB=NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
				cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
				if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA>0.99 && normNormalA<1.01) && (normNormalB>0.99 && normNormalB<1.01))
				{
					minCos = cos;
					idsWallNodes[0] = 0;
					idsWallNodes[1] = 3;
					idFreeNode = 1;
				}

				NormalA = Element[2].FastGetSolutionStepValue(NORMAL);
				NormalB = Element[3].FastGetSolutionStepValue(NORMAL);
				normNormalA=NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
				normNormalB=NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
				cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
				if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA>0.99 && normNormalA<1.01) && (normNormalB>0.99 && normNormalB<1.01))
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
				normNormalA=NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
				normNormalB=NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
				cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
				if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA>0.99 && normNormalA<1.01) && (normNormalB>0.99 && normNormalB<1.01))
				{
					minCos = cos;
					idsWallNodes[0] = 0;
					idsWallNodes[1] = 1;
					idFreeNode = 2;
				}

				NormalA = Element[0].FastGetSolutionStepValue(NORMAL);
				NormalB = Element[3].FastGetSolutionStepValue(NORMAL);
				normNormalA=NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
				normNormalB=NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
				cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
				if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA>0.99 && normNormalA<1.01) && (normNormalB>0.99 && normNormalB<1.01))
				{
					minCos = cos;
					idsWallNodes[0] = 0;
					idsWallNodes[1] = 3;
					idFreeNode = 2;
				}

				NormalA = Element[1].FastGetSolutionStepValue(NORMAL);
				NormalB = Element[3].FastGetSolutionStepValue(NORMAL);
				normNormalA=NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
				normNormalB=NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
				cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
				if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA>0.99 && normNormalA<1.01) && (normNormalB>0.99 && normNormalB<1.01))
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
				normNormalA=NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
				normNormalB=NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
				cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
				if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA>0.99 && normNormalA<1.01) && (normNormalB>0.99 && normNormalB<1.01))
				{
					minCos = cos;
					idsWallNodes[0] = 0;
					idsWallNodes[1] = 1;
					idFreeNode = 3;
				}

				NormalA = Element[0].FastGetSolutionStepValue(NORMAL);
				NormalB = Element[2].FastGetSolutionStepValue(NORMAL);
				normNormalA=NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
				normNormalB=NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
				cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
				if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA>0.99 && normNormalA<1.01) && (normNormalB>0.99 && normNormalB<1.01))
				{
					minCos = cos;
					idsWallNodes[0] = 0;
					idsWallNodes[1] = 2;
					idFreeNode = 3;
				}

				NormalA = Element[1].FastGetSolutionStepValue(NORMAL);
				NormalB = Element[2].FastGetSolutionStepValue(NORMAL);
				normNormalA=NormalA[0] * NormalA[0] + NormalA[1] * NormalA[1] + NormalA[2] * NormalA[2];
				normNormalB=NormalB[0] * NormalB[0] + NormalB[1] * NormalB[1] + NormalB[2] * NormalB[2];
				cos = NormalA[0] * NormalB[0] + NormalA[1] * NormalB[1] + NormalA[2] * NormalB[2];
				if (cos < minCos && (cos < cosTolerance && cos > -cosTolerance) && (normNormalA>0.99 && normNormalA<1.01) && (normNormalB>0.99 && normNormalB<1.01))
				{
					minCos = cos;
					idsWallNodes[0] = 1;
					idsWallNodes[1] = 2;
					idFreeNode = 3;
				}
			}

			if (minCos < cosTolerance && minCos> -cosTolerance)
			{

				bool alreadyAddedNode = false;
				unsigned int idA = Element[idsWallNodes[0]].GetId();
				unsigned int idB = Element[idsWallNodes[1]].GetId();
				double minimumDistanceToInstert = 1.3 * mrRemesh.Refine->CriticalRadius;
				array_1d<double, 3> CoorDifference = Element[idsWallNodes[0]].Coordinates() - Element[idsWallNodes[1]].Coordinates();
				double SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1];
				double separation = sqrt(SquaredLength);
				unsigned int idC = Element[idFreeNode].GetId();
				if (separation > minimumDistanceToInstert)
				{

					for (unsigned int i = 0; i < unsigned(CountNodes); i++)
					{
						if (idA == NodesIDToInterpolate[i][0] || idA == NodesIDToInterpolate[i][1] || idB == NodesIDToInterpolate[i][0] || idB == NodesIDToInterpolate[i][1])
						{
							alreadyAddedNode = true;
							break;
						}
					}
					if (alreadyAddedNode == false)
					{
						array_1d<double, 3> NewPosition = (Element[idsWallNodes[0]].Coordinates() + Element[idsWallNodes[1]].Coordinates()) * 0.5;
						NodesIDToInterpolate[CountNodes][0] = idA;
						NodesIDToInterpolate[CountNodes][1] = idB;
						if (Element[idFreeNode].IsNot(TO_ERASE))
						{
							NodesIDToInterpolate[CountNodes][2] = idC;
						}
						else
						{
							NodesIDToInterpolate[CountNodes][2] = idA;
						}
						CopyDofs(Element[idFreeNode].GetDofs(), NewDofs[CountNodes]);
						NewPositions[CountNodes] = NewPosition;
						CountNodes++;
						// std::cout << "  NewPosition  NewPosition NewPosition " << NewPosition << std::endl;
						// std::cout <<idsWallNodes[0] <<" idA " <<idA << std::endl;
						// std::cout <<idsWallNodes[1] <<" idB " <<idB << std::endl;
						// std::cout <<idFreeNode<< " idC " <<idC << std::endl;
					}
				}
			}
		}

		KRATOS_CATCH("")
	}

	void SelectEdgeToRefine2D(Element::GeometryType &Element,
							  std::vector<array_1d<double, 3>> &NewPositions,
							  std::vector<double> &BiggestVolumes,
							  std::vector<array_1d<unsigned int, 4>> &NodesIDToInterpolate,
							  std::vector<Node<3>::DofsContainerType> &NewDofs,
							  int &CountNodes,
							  int ElementsToRefine)
	{
		KRATOS_TRY

		const unsigned int nds = Element.size();

		unsigned int rigidNodes = 0;
		unsigned int boundaryNodes = 0;
		unsigned int freesurfaceNodes = 0;
		unsigned int inletNodes = 0;
		bool toEraseNodeFound = false;

		for (unsigned int pn = 0; pn < nds; pn++)
		{
			if (Element[pn].Is(RIGID))
			{
				rigidNodes++;
			}
			if (Element[pn].Is(BOUNDARY))
			{
				boundaryNodes++;
			}
			if (Element[pn].Is(TO_ERASE))
			{
				toEraseNodeFound = true;
			}
			if (Element[pn].Is(FREE_SURFACE))
			{
				freesurfaceNodes++;
			}
			if (Element[pn].Is(INLET))
			{
				inletNodes++;
			}
		}

		double limitEdgeLength = 1.4 * mrRemesh.Refine->CriticalRadius;
		double safetyCoefficient2D = 1.5;
		double penalization = 1.0;
		if (rigidNodes > 1)
		{
			// penalization=0.7;
			penalization = 0.8;
			if (inletNodes > 0)
			{
				penalization = 0.9;
			}
		}
		else if (rigidNodes > 0 && freesurfaceNodes > 0)
		{
			penalization = 0;
		}
		else if (freesurfaceNodes > 0)
		{
			penalization = 0.85;
		}

		double ElementalVolume = Element.Area();

		array_1d<double, 3> Edges(3, 0.0);
		array_1d<unsigned int, 3> FirstEdgeNode(3, 0);
		array_1d<unsigned int, 3> SecondEdgeNode(3, 0);
		double WallCharacteristicDistance = 0;
		array_1d<double, 3> CoorDifference = Element[1].Coordinates() - Element[0].Coordinates();
		// array_1d<double,3> CoorDifference(3,0.0);
		// noalias(CoorDifference) = Element[1].Coordinates() - Element[0].Coordinates();
		// CoorDifference = Element[1].Coordinates() - Element[0].Coordinates();
		double SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1];
		Edges[0] = sqrt(SquaredLength);
		FirstEdgeNode[0] = 0;
		SecondEdgeNode[0] = 1;
		if (Element[0].Is(RIGID) && Element[1].Is(RIGID))
		{
			WallCharacteristicDistance = Edges[0];
		}
		unsigned int Counter = 0;
		for (unsigned int i = 2; i < nds; i++)
		{
			for (unsigned int j = 0; j < i; j++)
			{
				noalias(CoorDifference) = Element[i].Coordinates() - Element[j].Coordinates();
				// CoorDifference = Element[i].Coordinates() - Element[j].Coordinates();
				SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1];
				Counter += 1;
				Edges[Counter] = sqrt(SquaredLength);
				FirstEdgeNode[Counter] = j;
				SecondEdgeNode[Counter] = i;
				if (Element[i].Is(RIGID) && Element[j].Is(RIGID) && Edges[Counter] > WallCharacteristicDistance)
				{
					WallCharacteristicDistance = Edges[Counter];
				}
			}
		}

		bool dangerousElement = false;
		if (rigidNodes > 1)
		{
			for (unsigned int i = 0; i < 3; i++)
			{
				if ((Edges[i] < WallCharacteristicDistance * safetyCoefficient2D && (Element[FirstEdgeNode[i]].Is(RIGID) || Element[SecondEdgeNode[i]].Is(RIGID))) ||
					(Element[FirstEdgeNode[i]].Is(RIGID) && Element[SecondEdgeNode[i]].Is(RIGID)))
				{
					Edges[i] = 0;
				}
				// if(Element[FirstEdgeNode[i]].Is(FREE_SURFACE) && Element[SecondEdgeNode[i]].Is(FREE_SURFACE)){
				//   Edges[i]=0;
				//   // Edges[i]*=penalizationFreeSurface;
				// }
				if ((Element[FirstEdgeNode[i]].Is(FREE_SURFACE) || Element[FirstEdgeNode[i]].Is(RIGID)) &&
					(Element[SecondEdgeNode[i]].Is(FREE_SURFACE) || Element[SecondEdgeNode[i]].Is(RIGID)))
				{
					Edges[i] = 0;
				}
			}
		}
		if ((Edges[0] == 0 && Edges[1] == 0 && Edges[2] == 0) || rigidNodes == 3)
		{
			dangerousElement = true;
		}

		if (dangerousElement == false && toEraseNodeFound == false)
		{

			// array_1d<double,3> NewPosition(3,0.0);
			unsigned int maxCount = 3;
			double LargestEdge = 0;

			for (unsigned int i = 0; i < 3; i++)
			{
				if (Edges[i] > LargestEdge)
				{
					maxCount = i;
					LargestEdge = Edges[i];
				}
			}

			if (CountNodes < ElementsToRefine && LargestEdge > limitEdgeLength)
			{
				array_1d<double, 3> NewPosition = (Element[FirstEdgeNode[maxCount]].Coordinates() + Element[SecondEdgeNode[maxCount]].Coordinates()) * 0.5;
				// noalias(NewPosition)=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
				// NewPosition=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
				NodesIDToInterpolate[CountNodes][0] = Element[FirstEdgeNode[maxCount]].GetId();
				NodesIDToInterpolate[CountNodes][1] = Element[SecondEdgeNode[maxCount]].GetId();
				if (Element[SecondEdgeNode[maxCount]].IsNot(RIGID))
				{
					CopyDofs(Element[SecondEdgeNode[maxCount]].GetDofs(), NewDofs[CountNodes]);
				}
				else if (Element[FirstEdgeNode[maxCount]].IsNot(RIGID))
				{
					CopyDofs(Element[FirstEdgeNode[maxCount]].GetDofs(), NewDofs[CountNodes]);
				}
				else
				{
					std::cout << "CAUTION! THIS IS A WALL EDGE" << std::endl;
				}
				BiggestVolumes[CountNodes] = ElementalVolume;
				NewPositions[CountNodes] = NewPosition;
				CountNodes++;
			}
			else if (freesurfaceNodes < 3 && rigidNodes < 3)
			{

				ElementalVolume *= penalization;
				for (int nn = 0; nn < ElementsToRefine; nn++)
				{
					if (ElementalVolume > BiggestVolumes[nn])
					{

						bool suitableElement = true;
						if (maxCount < 3 && LargestEdge > limitEdgeLength)
						{
							array_1d<double, 3> NewPosition = (Element[FirstEdgeNode[maxCount]].Coordinates() + Element[SecondEdgeNode[maxCount]].Coordinates()) * 0.5;
							// noalias(NewPosition)=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
							// NewPosition=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
							for (int j = 0; j < ElementsToRefine; j++)
							{
								if (NewPositions[j][0] == NewPosition[0] && NewPositions[j][1] == NewPosition[1])
								{
									suitableElement = false;
								}
							}
							if (suitableElement == true)
							{
								NodesIDToInterpolate[nn][0] = Element[FirstEdgeNode[maxCount]].GetId();
								NodesIDToInterpolate[nn][1] = Element[SecondEdgeNode[maxCount]].GetId();
								if (Element[SecondEdgeNode[maxCount]].IsNot(RIGID))
								{
									CopyDofs(Element[SecondEdgeNode[maxCount]].GetDofs(), NewDofs[nn]);
								}
								else if (Element[FirstEdgeNode[maxCount]].IsNot(RIGID))
								{
									CopyDofs(Element[FirstEdgeNode[maxCount]].GetDofs(), NewDofs[nn]);
								}
								else
								{
									std::cout << "CAUTION! THIS IS A WALL EDGE" << std::endl;
								}
								BiggestVolumes[nn] = ElementalVolume;
								NewPositions[nn] = NewPosition;
							}
						}

						break;
					}
				}
			}
		}

		KRATOS_CATCH("")
	}

	void SelectEdgeToRefine3D(Element::GeometryType &Element,
							  std::vector<array_1d<double, 3>> &NewPositions,
							  std::vector<double> &BiggestVolumes,
							  std::vector<array_1d<unsigned int, 4>> &NodesIDToInterpolate,
							  std::vector<Node<3>::DofsContainerType> &NewDofs,
							  int &CountNodes,
							  int ElementsToRefine)
	{
		KRATOS_TRY

		const unsigned int nds = Element.size();

		unsigned int rigidNodes = 0;
		unsigned int freesurfaceNodes = 0;
		unsigned int inletNodes = 0;
		bool toEraseNodeFound = false;

		for (unsigned int pn = 0; pn < nds; pn++)
		{
			if (Element[pn].Is(RIGID))
			{
				rigidNodes++;
			}
			if (Element[pn].Is(TO_ERASE))
			{
				toEraseNodeFound = true;
			}
			if (Element[pn].Is(FREE_SURFACE))
			{
				freesurfaceNodes++;
			}
			if (Element[pn].Is(INLET))
			{
				inletNodes++;
			}
		}

		double limitEdgeLength = 1.25 * mrRemesh.Refine->CriticalRadius;
		double safetyCoefficient3D = 1.6;
		double penalization = 1.0;
		if (rigidNodes > 2)
		{
			penalization = 0.7;
			if (inletNodes > 0)
			{
				penalization = 0.9;
			}
		}
		else if (rigidNodes > 0 && freesurfaceNodes > 0)
		{
			penalization = 0;
		}
		else if (freesurfaceNodes > 0)
		{
			penalization = 0.95;
		}

		// if(freesurfaceNodes>2){
		//   penalization=0.6;
		// }

		double ElementalVolume = Element.Volume();

		array_1d<double, 6> Edges(6, 0.0);
		array_1d<unsigned int, 6> FirstEdgeNode(6, 0);
		array_1d<unsigned int, 6> SecondEdgeNode(6, 0);
		double WallCharacteristicDistance = 0;
		array_1d<double, 3> CoorDifference = Element[1].Coordinates() - Element[0].Coordinates();
		// array_1d<double,3> CoorDifference(3,0.0);
		// noalias(CoorDifference) = Element[1].Coordinates() - Element[0].Coordinates();
		// CoorDifference = Element[1].Coordinates() - Element[0].Coordinates();
		double SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1] + CoorDifference[2] * CoorDifference[2];
		Edges[0] = sqrt(SquaredLength);
		FirstEdgeNode[0] = 0;
		SecondEdgeNode[0] = 1;
		if (Element[0].Is(RIGID) && Element[1].Is(RIGID))
		{
			WallCharacteristicDistance = Edges[0];
		}
		unsigned int Counter = 0;
		for (unsigned int i = 2; i < nds; i++)
		{
			for (unsigned int j = 0; j < i; j++)
			{
				noalias(CoorDifference) = Element[i].Coordinates() - Element[j].Coordinates();
				// CoorDifference = Element[i].Coordinates() - Element[j].Coordinates();
				SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1] + CoorDifference[2] * CoorDifference[2];
				Counter += 1;
				Edges[Counter] = sqrt(SquaredLength);
				FirstEdgeNode[Counter] = j;
				SecondEdgeNode[Counter] = i;
				if (Element[i].Is(RIGID) && Element[j].Is(RIGID) && Edges[Counter] > WallCharacteristicDistance)
				{
					WallCharacteristicDistance = Edges[Counter];
				}
			}
		}
		//Edges connectivity: Edges[0]=d01, Edges[1]=d20, Edges[2]=d21, Edges[3]=d30, Edges[4]=d31, Edges[5]=d32
		bool dangerousElement = false;
		if (rigidNodes > 1)
		{
			for (unsigned int i = 0; i < 6; i++)
			{
				if ((Edges[i] < WallCharacteristicDistance * safetyCoefficient3D && (Element[FirstEdgeNode[i]].Is(RIGID) || Element[SecondEdgeNode[i]].Is(RIGID))) ||
					(Element[FirstEdgeNode[i]].Is(RIGID) && Element[SecondEdgeNode[i]].Is(RIGID)))
				{
					Edges[i] = 0;
				}
				// if(Element[FirstEdgeNode[i]].Is(FREE_SURFACE) && Element[SecondEdgeNode[i]].Is(FREE_SURFACE)){
				//   Edges[i]=0;
				// }
				if ((Element[FirstEdgeNode[i]].Is(FREE_SURFACE) || Element[FirstEdgeNode[i]].Is(RIGID)) &&
					(Element[SecondEdgeNode[i]].Is(FREE_SURFACE) || Element[SecondEdgeNode[i]].Is(RIGID)))
				{
					Edges[i] = 0;
				}
			}
		}
		else if (rigidNodes == 1)
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

		//just to fill the vector
		if (dangerousElement == false && toEraseNodeFound == false)
		{

			// array_1d<double,3> NewPosition(3,0.0);
			unsigned int maxCount = 6;
			double LargestEdge = 0;

			for (unsigned int i = 0; i < 6; i++)
			{
				if (Edges[i] > LargestEdge)
				{
					maxCount = i;
					LargestEdge = Edges[i];
				}
			}

			if (CountNodes < ElementsToRefine && LargestEdge > limitEdgeLength)
			{
				array_1d<double, 3> NewPosition = (Element[FirstEdgeNode[maxCount]].Coordinates() + Element[SecondEdgeNode[maxCount]].Coordinates()) * 0.5;
				// noalias(NewPosition)=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
				// NewPosition=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
				NodesIDToInterpolate[CountNodes][0] = Element[FirstEdgeNode[maxCount]].GetId();
				NodesIDToInterpolate[CountNodes][1] = Element[SecondEdgeNode[maxCount]].GetId();
				if (Element[SecondEdgeNode[maxCount]].IsNot(RIGID))
				{
					CopyDofs(Element[SecondEdgeNode[maxCount]].GetDofs(), NewDofs[CountNodes]);
				}
				else if (Element[FirstEdgeNode[maxCount]].IsNot(RIGID))
				{
					CopyDofs(Element[FirstEdgeNode[maxCount]].GetDofs(), NewDofs[CountNodes]);
				}
				else
				{
					std::cout << "CAUTION! THIS IS A WALL EDGE" << std::endl;
				}
				BiggestVolumes[CountNodes] = ElementalVolume;
				NewPositions[CountNodes] = NewPosition;
				CountNodes++;
			}
			else if (freesurfaceNodes < 4 && rigidNodes < 4)
			{

				ElementalVolume *= penalization;
				for (int nn = 0; nn < ElementsToRefine; nn++)
				{
					if (ElementalVolume > BiggestVolumes[nn])
					{

						bool suitableElement = true;

						if (maxCount < 6 && LargestEdge > limitEdgeLength)
						{
							array_1d<double, 3> NewPosition = (Element[FirstEdgeNode[maxCount]].Coordinates() + Element[SecondEdgeNode[maxCount]].Coordinates()) * 0.5;
							// noalias(NewPosition)=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
							// NewPosition=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
							for (int j = 0; j < ElementsToRefine; j++)
							{
								if (NewPositions[j][0] == NewPosition[0] && NewPositions[j][1] == NewPosition[1] && NewPositions[j][2] == NewPosition[2])
								{
									suitableElement = false; //this is a repeated node, I have already choose this from another element
								}
							}
							if (suitableElement == true)
							{
								NodesIDToInterpolate[nn][0] = Element[FirstEdgeNode[maxCount]].GetId();
								NodesIDToInterpolate[nn][1] = Element[SecondEdgeNode[maxCount]].GetId();
								if (Element[SecondEdgeNode[maxCount]].IsNot(RIGID))
								{
									CopyDofs(Element[SecondEdgeNode[maxCount]].GetDofs(), NewDofs[nn]);
								}
								else if (Element[FirstEdgeNode[maxCount]].IsNot(RIGID))
								{
									CopyDofs(Element[FirstEdgeNode[maxCount]].GetDofs(), NewDofs[nn]);
								}
								else
								{
									std::cout << "CAUTION! THIS IS A WALL EDGE" << std::endl;
								}
								BiggestVolumes[nn] = ElementalVolume;
								NewPositions[nn] = NewPosition;
							}
						}

						break;
					}
				}
			}
		}

		KRATOS_CATCH("")
	}

	void SelectEdgeToRefine2DWithRefinement(Element::GeometryType &Element,
											std::vector<array_1d<double, 3>> &NewPositions,
											std::vector<double> &BiggestVolumes,
											std::vector<array_1d<unsigned int, 4>> &NodesIDToInterpolate,
											std::vector<Node<3>::DofsContainerType> &NewDofs,
											int &CountNodes,
											int ElementsToRefine)
	{
		KRATOS_TRY

		const unsigned int nds = Element.size();

		unsigned int rigidNodes = 0;
		unsigned int freesurfaceNodes = 0;
		unsigned int inletNodes = 0;
		bool toEraseNodeFound = false;

		double meanMeshSize = mrRemesh.Refine->CriticalRadius;
		const ProcessInfo &rCurrentProcessInfo = mrModelPart.GetProcessInfo();
		double currentTime = rCurrentProcessInfo[TIME];
		double initialTime = mrRemesh.RefiningBoxInitialTime;
		double finalTime = mrRemesh.RefiningBoxFinalTime;
		bool refiningBox = mrRemesh.UseRefiningBox;
		double distance = 2.0 * meanMeshSize;
		bool penalizationRigid = false;
		double seperation = 0;
		double coefficient = 0;
		if (!(refiningBox == true && currentTime > initialTime && currentTime < finalTime))
		{
			refiningBox = false;
		}

		for (unsigned int pn = 0; pn < nds; pn++)
		{
			if (Element[pn].Is(RIGID))
			{
				rigidNodes++;
			}
			if (Element[pn].Is(TO_ERASE))
			{
				toEraseNodeFound = true;
			}
			if (Element[pn].Is(FREE_SURFACE))
			{
				freesurfaceNodes++;
			}
			if (Element[pn].Is(INLET))
			{
				inletNodes++;
			}

			if (refiningBox == true)
			{

				array_1d<double, 3> RefiningBoxMinimumPoint = mrRemesh.RefiningBoxMinimumPoint;
				array_1d<double, 3> RefiningBoxMaximumPoint = mrRemesh.RefiningBoxMaximumPoint;
				array_1d<double, 3> minExternalPoint = mrRemesh.RefiningBoxMinExternalPoint;
				array_1d<double, 3> minInternalPoint = mrRemesh.RefiningBoxMinInternalPoint;
				array_1d<double, 3> maxExternalPoint = mrRemesh.RefiningBoxMaxExternalPoint;
				array_1d<double, 3> maxInternalPoint = mrRemesh.RefiningBoxMaxInternalPoint;
				if (mrRemesh.Refine->CriticalRadius > mrRemesh.RefiningBoxMeshSize)
				{
					if (Element[pn].X() > RefiningBoxMinimumPoint[0] && Element[pn].Y() > RefiningBoxMinimumPoint[1] &&
						Element[pn].X() < RefiningBoxMaximumPoint[0] && Element[pn].Y() < RefiningBoxMaximumPoint[1])
					{
						meanMeshSize = mrRemesh.RefiningBoxMeshSize;
					}
					else if ((Element[pn].X() < RefiningBoxMinimumPoint[0] && Element[pn].X() > (minExternalPoint[0] - distance) && Element[pn].Y() > minExternalPoint[1] && Element[pn].Y() < maxExternalPoint[1]))
					{
						seperation = Element[pn].X() - RefiningBoxMinimumPoint[0];
						coefficient = fabs(seperation) / (distance + meanMeshSize);
						meanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
						penalizationRigid = true;
					}
					else if ((Element[pn].Y() < RefiningBoxMinimumPoint[1] && Element[pn].Y() > (minExternalPoint[1] - distance) && Element[pn].X() > minExternalPoint[0] && Element[pn].X() < maxExternalPoint[0]))
					{
						seperation = Element[pn].Y() - RefiningBoxMinimumPoint[1];
						coefficient = fabs(seperation) / (distance + meanMeshSize);
						meanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
						penalizationRigid = true;
					}
					else if ((Element[pn].X() > RefiningBoxMaximumPoint[0] && Element[pn].X() < (maxExternalPoint[0] + distance) && Element[pn].Y() > minExternalPoint[1] && Element[pn].Y() < maxExternalPoint[1]))
					{
						seperation = Element[pn].X() - RefiningBoxMaximumPoint[0];
						coefficient = fabs(seperation) / (distance + meanMeshSize);
						meanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
						penalizationRigid = true;
					}
					else if ((Element[pn].Y() > RefiningBoxMaximumPoint[1] && Element[pn].Y() < (maxExternalPoint[1] + distance) && Element[pn].X() > minExternalPoint[0] && Element[pn].X() < maxExternalPoint[0]))
					{
						seperation = Element[pn].Y() - RefiningBoxMaximumPoint[1];
						coefficient = fabs(seperation) / (distance + meanMeshSize);
						meanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
						penalizationRigid = true;
					}
				}
				else
				{
					distance = 2.0 * mrRemesh.RefiningBoxMeshSize;
					if (Element[pn].X() > (minInternalPoint[0] + distance) && Element[pn].X() < (maxInternalPoint[0] - distance) &&
						Element[pn].Y() > (minInternalPoint[1] + distance) && Element[pn].Y() < (maxInternalPoint[1] - distance))
					{
						meanMeshSize = mrRemesh.RefiningBoxMeshSize; // in the internal domain the size is the one given by the user
					}
					else if ((Element[pn].X() > RefiningBoxMinimumPoint[0] && Element[pn].X() < (minInternalPoint[0] + distance) && Element[pn].Y() > minExternalPoint[1] && Element[pn].Y() < maxExternalPoint[1]))
					{
						seperation = (minInternalPoint[0] + distance) - Element[pn].X();
						coefficient = fabs(seperation) / (distance + mrRemesh.RefiningBoxMeshSize);
						meanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
						if (meanMeshSize < 0)
						{
							meanMeshSize = mrRemesh.Refine->CriticalRadius;
						}
						penalizationRigid = true;
					}
					else if ((Element[pn].Y() > RefiningBoxMinimumPoint[1] && Element[pn].Y() < (minInternalPoint[1] + distance) && Element[pn].X() > minExternalPoint[0] && Element[pn].X() < maxExternalPoint[0]))
					{
						seperation = (minInternalPoint[1] + distance) - Element[pn].Y();
						coefficient = fabs(seperation) / (distance + mrRemesh.RefiningBoxMeshSize);
						meanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
						if (meanMeshSize < 0)
						{
							meanMeshSize = mrRemesh.Refine->CriticalRadius;
						}
						penalizationRigid = true;
					}
					else if ((Element[pn].X() < RefiningBoxMaximumPoint[0] && Element[pn].X() > (maxInternalPoint[0] - distance) && Element[pn].Y() > minExternalPoint[1] && Element[pn].Y() < maxExternalPoint[1]))
					{
						seperation = (maxInternalPoint[0] - distance) - Element[pn].X();
						coefficient = fabs(seperation) / (distance + mrRemesh.RefiningBoxMeshSize);
						meanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
						if (meanMeshSize < 0)
						{
							meanMeshSize = mrRemesh.Refine->CriticalRadius;
						}
						penalizationRigid = true;
					}
					else if ((Element[pn].Y() < RefiningBoxMaximumPoint[1] && Element[pn].Y() > (maxInternalPoint[1] - distance) && Element[pn].X() > minExternalPoint[0] && Element[pn].X() < maxExternalPoint[0]))
					{
						seperation = (maxInternalPoint[1] - distance) - Element[pn].Y();
						coefficient = fabs(seperation) / (distance + mrRemesh.RefiningBoxMeshSize);
						meanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
						if (meanMeshSize < 0)
						{
							meanMeshSize = mrRemesh.Refine->CriticalRadius;
						}
						penalizationRigid = true;
					}
				}
			}
		}

		double penalization = 1.0;
		if (refiningBox == true)
		{
			if (freesurfaceNodes > 0)
			{
				penalization = 1.2; //to avoid to gain too much volume during remeshing step
			}

			if (rigidNodes > 0 && penalizationRigid == true)
			{
				penalization = 1.2;
			}
		}

		double limitEdgeLength = 1.9 * meanMeshSize * penalization;
		double safetyCoefficient2D = 1.5;

		double ElementalVolume = Element.Area();

		array_1d<double, 3> Edges(3, 0.0);
		array_1d<unsigned int, 3> FirstEdgeNode(3, 0);
		array_1d<unsigned int, 3> SecondEdgeNode(3, 0);
		double WallCharacteristicDistance = 0;
		array_1d<double, 3> CoorDifference = Element[1].Coordinates() - Element[0].Coordinates();
		// array_1d<double,3> CoorDifference(3,0.0);
		// noalias(CoorDifference) = Element[1].Coordinates() - Element[0].Coordinates();
		// CoorDifference = Element[1].Coordinates() - Element[0].Coordinates();
		double SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1];
		Edges[0] = sqrt(SquaredLength);
		FirstEdgeNode[0] = 0;
		SecondEdgeNode[0] = 1;
		if (Element[0].Is(RIGID) && Element[1].Is(RIGID))
		{
			WallCharacteristicDistance = Edges[0];
		}
		unsigned int Counter = 0;
		for (unsigned int i = 2; i < nds; i++)
		{
			for (unsigned int j = 0; j < i; j++)
			{
				noalias(CoorDifference) = Element[i].Coordinates() - Element[j].Coordinates();
				// CoorDifference = Element[i].Coordinates() - Element[j].Coordinates();
				SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1];
				Counter += 1;
				Edges[Counter] = sqrt(SquaredLength);
				FirstEdgeNode[Counter] = j;
				SecondEdgeNode[Counter] = i;
				if (Element[i].Is(RIGID) && Element[j].Is(RIGID) && Edges[Counter] > WallCharacteristicDistance)
				{
					WallCharacteristicDistance = Edges[Counter];
				}
			}
		}

		bool dangerousElement = false;
		if (rigidNodes > 1)
		{
			for (unsigned int i = 0; i < 3; i++)
			{
				if ((Edges[i] < WallCharacteristicDistance * safetyCoefficient2D && (Element[FirstEdgeNode[i]].Is(RIGID) || Element[SecondEdgeNode[i]].Is(RIGID))) ||
					(Element[FirstEdgeNode[i]].Is(RIGID) && Element[SecondEdgeNode[i]].Is(RIGID)))
				{
					Edges[i] = 0;
				}
				// if(Element[FirstEdgeNode[i]].Is(FREE_SURFACE) && Element[SecondEdgeNode[i]].Is(FREE_SURFACE)){
				//   Edges[i]=0;
				//   // Edges[i]*=penalizationFreeSurface;
				// }
				if ((Element[FirstEdgeNode[i]].Is(FREE_SURFACE) || Element[FirstEdgeNode[i]].Is(RIGID)) &&
					(Element[SecondEdgeNode[i]].Is(FREE_SURFACE) || Element[SecondEdgeNode[i]].Is(RIGID)))
				{
					Edges[i] = 0;
				}
			}
		}
		if ((Edges[0] == 0 && Edges[1] == 0 && Edges[2] == 0) || rigidNodes == 3)
		{
			dangerousElement = true;
		}

		if (dangerousElement == false && toEraseNodeFound == false)
		{

			// array_1d<double,3> NewPosition(3,0.0);
			unsigned int maxCount = 3;
			double LargestEdge = 0;

			for (unsigned int i = 0; i < 3; i++)
			{
				if (Edges[i] > LargestEdge)
				{
					maxCount = i;
					LargestEdge = Edges[i];
				}
			}

			if (CountNodes < ElementsToRefine && LargestEdge > limitEdgeLength)
			{

				bool newNode = true;
				for (unsigned int i = 0; i < unsigned(CountNodes); i++)
				{
					if ((NodesIDToInterpolate[i][0] == Element[FirstEdgeNode[maxCount]].GetId() && NodesIDToInterpolate[i][1] == Element[SecondEdgeNode[maxCount]].GetId()) ||
						(NodesIDToInterpolate[i][1] == Element[FirstEdgeNode[maxCount]].GetId() && NodesIDToInterpolate[i][0] == Element[SecondEdgeNode[maxCount]].GetId()))
					{
						newNode = false;
					}
				}
				if (newNode == true)
				{
					array_1d<double, 3> NewPosition = (Element[FirstEdgeNode[maxCount]].Coordinates() + Element[SecondEdgeNode[maxCount]].Coordinates()) * 0.5;
					// noalias(NewPosition)=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
					// NewPosition=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
					NodesIDToInterpolate[CountNodes][0] = Element[FirstEdgeNode[maxCount]].GetId();
					NodesIDToInterpolate[CountNodes][1] = Element[SecondEdgeNode[maxCount]].GetId();
					if (Element[SecondEdgeNode[maxCount]].IsNot(RIGID))
					{
						CopyDofs(Element[SecondEdgeNode[maxCount]].GetDofs(), NewDofs[CountNodes]);
					}
					else if (Element[FirstEdgeNode[maxCount]].IsNot(RIGID))
					{
						CopyDofs(Element[FirstEdgeNode[maxCount]].GetDofs(), NewDofs[CountNodes]);
					}
					else
					{
						std::cout << "CAUTION! THIS IS A WALL EDGE" << std::endl;
					}

					BiggestVolumes[CountNodes] = ElementalVolume;
					NewPositions[CountNodes] = NewPosition;
					CountNodes++;
				}
			}
		}

		KRATOS_CATCH("")
	}

	void SelectEdgeToRefine3DWithRefinement(Element::GeometryType &Element,
											std::vector<array_1d<double, 3>> &NewPositions,
											std::vector<double> &BiggestVolumes,
											std::vector<array_1d<unsigned int, 4>> &NodesIDToInterpolate,
											std::vector<Node<3>::DofsContainerType> &NewDofs,
											int &CountNodes,
											int ElementsToRefine)
	{
		KRATOS_TRY

		const unsigned int nds = Element.size();

		unsigned int rigidNodes = 0;
		unsigned int freesurfaceNodes = 0;
		unsigned int inletNodes = 0;
		bool toEraseNodeFound = false;

		double meanMeshSize = mrRemesh.Refine->CriticalRadius;
		const ProcessInfo &rCurrentProcessInfo = mrModelPart.GetProcessInfo();
		double currentTime = rCurrentProcessInfo[TIME];
		double initialTime = mrRemesh.RefiningBoxInitialTime;
		double finalTime = mrRemesh.RefiningBoxFinalTime;
		bool refiningBox = mrRemesh.UseRefiningBox;
		double distance = 2.0 * meanMeshSize;
		bool penalizationRigid = false;
		double seperation = 0;
		double coefficient = 0;
		if (!(refiningBox == true && currentTime > initialTime && currentTime < finalTime))
		{
			refiningBox = false;
		}

		for (unsigned int pn = 0; pn < nds; pn++)
		{
			if (Element[pn].Is(RIGID))
			{
				rigidNodes++;
			}
			if (Element[pn].Is(TO_ERASE))
			{
				toEraseNodeFound = true;
			}
			if (Element[pn].Is(FREE_SURFACE))
			{
				freesurfaceNodes++;
			}
			if (Element[pn].Is(INLET))
			{
				inletNodes++;
			}

			if (refiningBox == true)
			{

				array_1d<double, 3> RefiningBoxMinimumPoint = mrRemesh.RefiningBoxMinimumPoint;
				array_1d<double, 3> RefiningBoxMaximumPoint = mrRemesh.RefiningBoxMaximumPoint;
				array_1d<double, 3> minExternalPoint = mrRemesh.RefiningBoxMinExternalPoint;
				array_1d<double, 3> minInternalPoint = mrRemesh.RefiningBoxMinInternalPoint;
				array_1d<double, 3> maxExternalPoint = mrRemesh.RefiningBoxMaxExternalPoint;
				array_1d<double, 3> maxInternalPoint = mrRemesh.RefiningBoxMaxInternalPoint;
				if (mrRemesh.Refine->CriticalRadius > mrRemesh.RefiningBoxMeshSize)
				{
					if (Element[pn].X() > RefiningBoxMinimumPoint[0] && Element[pn].X() < RefiningBoxMaximumPoint[0] &&
						Element[pn].Y() > RefiningBoxMinimumPoint[1] && Element[pn].Y() < RefiningBoxMaximumPoint[1] &&
						Element[pn].Z() < RefiningBoxMinimumPoint[2] && Element[pn].Z() < RefiningBoxMaximumPoint[2])
					{
						meanMeshSize = mrRemesh.RefiningBoxMeshSize;
					}
					else if ((Element[pn].X() < RefiningBoxMinimumPoint[0] && Element[pn].X() > (minExternalPoint[0] - distance) && Element[pn].Y() > minExternalPoint[1] && Element[pn].Y() < maxExternalPoint[1] && Element[pn].Z() > minExternalPoint[2] && Element[pn].Z() < maxExternalPoint[2]))
					{
						seperation = Element[pn].X() - RefiningBoxMinimumPoint[0];
						coefficient = fabs(seperation) / (distance + meanMeshSize);
						meanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
						penalizationRigid = true;
					}
					else if ((Element[pn].Y() < RefiningBoxMinimumPoint[1] && Element[pn].Y() > (minExternalPoint[1] - distance) && Element[pn].X() > minExternalPoint[0] && Element[pn].X() < maxExternalPoint[0] && Element[pn].Z() > minExternalPoint[2] && Element[pn].Z() < maxExternalPoint[2]))
					{
						seperation = Element[pn].Y() - RefiningBoxMinimumPoint[1];
						coefficient = fabs(seperation) / (distance + meanMeshSize);
						meanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
						penalizationRigid = true;
					}
					else if ((Element[pn].Z() < RefiningBoxMinimumPoint[2] && Element[pn].Z() > (minExternalPoint[2] - distance) && Element[pn].X() > minExternalPoint[0] && Element[pn].X() < maxExternalPoint[0] && Element[pn].Y() > minExternalPoint[1] && Element[pn].Y() < maxExternalPoint[1]))
					{
						seperation = Element[pn].Z() - RefiningBoxMinimumPoint[2];
						coefficient = fabs(seperation) / (distance + meanMeshSize);
						meanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
						penalizationRigid = true;
					}
					else if ((Element[pn].X() > RefiningBoxMaximumPoint[0] && Element[pn].X() < (maxExternalPoint[0] + distance) && Element[pn].Y() > minExternalPoint[1] && Element[pn].Y() < maxExternalPoint[1] && Element[pn].Z() > minExternalPoint[2] && Element[pn].Z() < maxExternalPoint[2]))
					{
						seperation = Element[pn].X() - RefiningBoxMaximumPoint[0];
						coefficient = fabs(seperation) / (distance + meanMeshSize);
						meanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
						penalizationRigid = true;
					}
					else if ((Element[pn].Y() > RefiningBoxMaximumPoint[1] && Element[pn].Y() < (maxExternalPoint[1] + distance) && Element[pn].X() > minExternalPoint[0] && Element[pn].X() < maxExternalPoint[0] && Element[pn].Z() > minExternalPoint[2] && Element[pn].Z() < maxExternalPoint[2]))
					{
						seperation = Element[pn].Y() - RefiningBoxMaximumPoint[1];
						coefficient = fabs(seperation) / (distance + meanMeshSize);
						meanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
						penalizationRigid = true;
					}
					else if ((Element[pn].Z() > RefiningBoxMaximumPoint[2] && Element[pn].Z() < (maxExternalPoint[2] + distance) && Element[pn].X() > minExternalPoint[0] && Element[pn].X() < maxExternalPoint[0] && Element[pn].Y() > minExternalPoint[1] && Element[pn].Y() < maxExternalPoint[1]))
					{
						seperation = Element[pn].Z() - RefiningBoxMaximumPoint[2];
						coefficient = fabs(seperation) / (distance + meanMeshSize);
						meanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
						penalizationRigid = true;
					}
				}
				else
				{
					distance = 2.0 * mrRemesh.RefiningBoxMeshSize;
					if (Element[pn].X() > (minInternalPoint[0] + distance) && Element[pn].X() < (maxInternalPoint[0] - distance) &&
						Element[pn].Y() > (minInternalPoint[1] + distance) && Element[pn].Y() < (maxInternalPoint[1] - distance) &&
						Element[pn].Z() > (minInternalPoint[2] + distance) && Element[pn].Z() < (maxInternalPoint[2] - distance))
					{
						meanMeshSize = mrRemesh.RefiningBoxMeshSize; // in the internal domain the size is the one given by the user
					}
					else if ((Element[pn].X() > RefiningBoxMinimumPoint[0] && Element[pn].X() < (minInternalPoint[0] + distance) && Element[pn].Y() > minExternalPoint[1] && Element[pn].Y() < maxExternalPoint[1] && Element[pn].Z() > minExternalPoint[2] && Element[pn].Z() < maxExternalPoint[2]))
					{
						seperation = (minInternalPoint[0] + distance) - Element[pn].X();
						coefficient = fabs(seperation) / (distance + mrRemesh.RefiningBoxMeshSize);
						meanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
						if (meanMeshSize < 0)
						{
							meanMeshSize = mrRemesh.Refine->CriticalRadius;
						}
						penalizationRigid = true;
					}
					else if ((Element[pn].Y() > RefiningBoxMinimumPoint[1] && Element[pn].Y() < (minInternalPoint[1] + distance) && Element[pn].X() > minExternalPoint[0] && Element[pn].X() < maxExternalPoint[0] && Element[pn].Z() > minExternalPoint[2] && Element[pn].Z() < maxExternalPoint[2]))
					{
						seperation = (minInternalPoint[1] + distance) - Element[pn].Y();
						coefficient = fabs(seperation) / (distance + mrRemesh.RefiningBoxMeshSize);
						meanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
						if (meanMeshSize < 0)
						{
							meanMeshSize = mrRemesh.Refine->CriticalRadius;
						}
						penalizationRigid = true;
					}
					else if ((Element[pn].Z() > RefiningBoxMinimumPoint[2] && Element[pn].Z() < (minInternalPoint[2] + distance) && Element[pn].X() > minExternalPoint[0] && Element[pn].X() < maxExternalPoint[0] && Element[pn].Y() > minExternalPoint[1] && Element[pn].Y() < maxExternalPoint[1]))
					{
						seperation = (minInternalPoint[2] + distance) - Element[pn].Z();
						coefficient = fabs(seperation) / (distance + mrRemesh.RefiningBoxMeshSize);
						meanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
						if (meanMeshSize < 0)
						{
							meanMeshSize = mrRemesh.Refine->CriticalRadius;
						}
						penalizationRigid = true;
					}
					else if ((Element[pn].X() < RefiningBoxMaximumPoint[0] && Element[pn].X() > (maxInternalPoint[0] - distance) && Element[pn].Y() > minExternalPoint[1] && Element[pn].Y() < maxExternalPoint[1] && Element[pn].Z() > minExternalPoint[2] && Element[pn].Z() < maxExternalPoint[2]))
					{
						seperation = (maxInternalPoint[0] - distance) - Element[pn].X();
						coefficient = fabs(seperation) / (distance + mrRemesh.RefiningBoxMeshSize);
						meanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
						if (meanMeshSize < 0)
						{
							meanMeshSize = mrRemesh.Refine->CriticalRadius;
						}
						penalizationRigid = true;
					}
					else if ((Element[pn].Y() < RefiningBoxMaximumPoint[1] && Element[pn].Y() > (maxInternalPoint[1] - distance) && Element[pn].X() > minExternalPoint[0] && Element[pn].X() < maxExternalPoint[0] && Element[pn].Z() > minExternalPoint[2] && Element[pn].Z() < maxExternalPoint[2]))
					{
						seperation = (maxInternalPoint[1] - distance) - Element[pn].Y();
						coefficient = fabs(seperation) / (distance + mrRemesh.RefiningBoxMeshSize);
						meanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
						if (meanMeshSize < 0)
						{
							meanMeshSize = mrRemesh.Refine->CriticalRadius;
						}
						penalizationRigid = true;
					}
					else if ((Element[pn].Z() < RefiningBoxMaximumPoint[2] && Element[pn].Z() > (maxInternalPoint[2] - distance) && Element[pn].X() > minExternalPoint[0] && Element[pn].X() < maxExternalPoint[0] && Element[pn].Y() > minExternalPoint[1] && Element[pn].Y() < maxExternalPoint[1]))
					{
						seperation = (maxInternalPoint[2] - distance) - Element[pn].Z();
						coefficient = fabs(seperation) / (distance + mrRemesh.RefiningBoxMeshSize);
						meanMeshSize = (1 - coefficient) * mrRemesh.RefiningBoxMeshSize + coefficient * mrRemesh.Refine->CriticalRadius;
						if (meanMeshSize < 0)
						{
							meanMeshSize = mrRemesh.Refine->CriticalRadius;
						}
						penalizationRigid = true;
					}
				}
			}
		}
		double penalization = 1.0;
		if (refiningBox == true)
		{
			if (freesurfaceNodes > 0)
			{
				penalization = 1.2; //to avoid to gain too much volume during remeshing step
			}

			if (rigidNodes > 0 && penalizationRigid == true)
			{
				penalization = 1.15;
			}
		}

		double limitEdgeLength = 1.6 * meanMeshSize * penalization;
		double safetyCoefficient3D = 1.6;

		double ElementalVolume = Element.Volume();

		array_1d<double, 6> Edges(6, 0.0);
		array_1d<unsigned int, 6> FirstEdgeNode(6, 0);
		array_1d<unsigned int, 6> SecondEdgeNode(6, 0);
		double WallCharacteristicDistance = 0;
		array_1d<double, 3> CoorDifference = Element[1].Coordinates() - Element[0].Coordinates();
		// array_1d<double,3> CoorDifference(3,0.0);
		// noalias(CoorDifference) = Element[1].Coordinates() - Element[0].Coordinates();
		// CoorDifference = Element[1].Coordinates() - Element[0].Coordinates();
		double SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1] + CoorDifference[2] * CoorDifference[2];
		Edges[0] = sqrt(SquaredLength);
		FirstEdgeNode[0] = 0;
		SecondEdgeNode[0] = 1;
		if (Element[0].Is(RIGID) && Element[1].Is(RIGID))
		{
			WallCharacteristicDistance = Edges[0];
		}
		unsigned int Counter = 0;
		for (unsigned int i = 2; i < nds; i++)
		{
			for (unsigned int j = 0; j < i; j++)
			{
				noalias(CoorDifference) = Element[i].Coordinates() - Element[j].Coordinates();
				// CoorDifference = Element[i].Coordinates() - Element[j].Coordinates();
				SquaredLength = CoorDifference[0] * CoorDifference[0] + CoorDifference[1] * CoorDifference[1] + CoorDifference[2] * CoorDifference[2];
				Counter += 1;
				Edges[Counter] = sqrt(SquaredLength);
				FirstEdgeNode[Counter] = j;
				SecondEdgeNode[Counter] = i;
				if (Element[i].Is(RIGID) && Element[j].Is(RIGID) && Edges[Counter] > WallCharacteristicDistance)
				{
					WallCharacteristicDistance = Edges[Counter];
				}
			}
		}
		//Edges connectivity: Edges[0]=d01, Edges[1]=d20, Edges[2]=d21, Edges[3]=d30, Edges[4]=d31, Edges[5]=d32
		bool dangerousElement = false;
		if (rigidNodes > 1)
		{
			for (unsigned int i = 0; i < 6; i++)
			{
				if ((Edges[i] < WallCharacteristicDistance * safetyCoefficient3D && (Element[FirstEdgeNode[i]].Is(RIGID) || Element[SecondEdgeNode[i]].Is(RIGID))) ||
					(Element[FirstEdgeNode[i]].Is(RIGID) && Element[SecondEdgeNode[i]].Is(RIGID)))
				{
					Edges[i] = 0;
				}
				// if(Element[FirstEdgeNode[i]].Is(FREE_SURFACE) && Element[SecondEdgeNode[i]].Is(FREE_SURFACE)){
				//   Edges[i]=0;
				// }
				if ((Element[FirstEdgeNode[i]].Is(FREE_SURFACE) || Element[FirstEdgeNode[i]].Is(RIGID)) &&
					(Element[SecondEdgeNode[i]].Is(FREE_SURFACE) || Element[SecondEdgeNode[i]].Is(RIGID)))
				{
					Edges[i] = 0;
				}
			}
		}
		else if (rigidNodes == 1)
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

		//just to fill the vector
		if (dangerousElement == false && toEraseNodeFound == false)
		{

			// array_1d<double,3> NewPosition(3,0.0);
			unsigned int maxCount = 6;
			double LargestEdge = 0;

			for (unsigned int i = 0; i < 6; i++)
			{
				if (Edges[i] > LargestEdge)
				{
					maxCount = i;
					LargestEdge = Edges[i];
				}
			}

			if (CountNodes < ElementsToRefine && LargestEdge > limitEdgeLength)
			{

				bool newNode = true;
				for (unsigned int i = 0; i < unsigned(CountNodes); i++)
				{
					if ((NodesIDToInterpolate[i][0] == Element[FirstEdgeNode[maxCount]].GetId() && NodesIDToInterpolate[i][1] == Element[SecondEdgeNode[maxCount]].GetId()) ||
						(NodesIDToInterpolate[i][1] == Element[FirstEdgeNode[maxCount]].GetId() && NodesIDToInterpolate[i][0] == Element[SecondEdgeNode[maxCount]].GetId()))
					{
						newNode = false;
					}
				}
				if (newNode == true)
				{
					array_1d<double, 3> NewPosition = (Element[FirstEdgeNode[maxCount]].Coordinates() + Element[SecondEdgeNode[maxCount]].Coordinates()) * 0.5;
					// noalias(NewPosition)=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
					// NewPosition=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
					NodesIDToInterpolate[CountNodes][0] = Element[FirstEdgeNode[maxCount]].GetId();
					NodesIDToInterpolate[CountNodes][1] = Element[SecondEdgeNode[maxCount]].GetId();
					if (Element[SecondEdgeNode[maxCount]].IsNot(RIGID))
					{
						CopyDofs(Element[SecondEdgeNode[maxCount]].GetDofs(), NewDofs[CountNodes]);
					}
					else if (Element[FirstEdgeNode[maxCount]].IsNot(RIGID))
					{
						CopyDofs(Element[FirstEdgeNode[maxCount]].GetDofs(), NewDofs[CountNodes]);
					}
					else
					{
						std::cout << "CAUTION! THIS IS A WALL EDGE" << std::endl;
					}

					BiggestVolumes[CountNodes] = ElementalVolume;
					NewPositions[CountNodes] = NewPosition;
					CountNodes++;
				}
			}
		}

		KRATOS_CATCH("")
	}

	void CreateAndAddNewNodesInCornerWall(std::vector<array_1d<double, 3>> &NewPositions,
										  std::vector<array_1d<unsigned int, 4>> &NodesIDToInterpolate,
										  std::vector<Node<3>::DofsContainerType> &NewDofs,
										  int ElementsToRefine,
										  unsigned int &maxId)
	{
		KRATOS_TRY

		const unsigned int dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

		std::vector<Node<3>::Pointer> list_of_new_nodes;

		//assign data to dofs
		VariablesList &VariablesList = mrModelPart.GetNodalSolutionStepVariablesList();

		for (unsigned int nn = 0; nn < NewPositions.size(); nn++)
		{

			unsigned int id = maxId + 1 + nn;

			double x = NewPositions[nn][0];
			double y = NewPositions[nn][1];
			double z = 0;
			if (dimension == 3)
				z = NewPositions[nn][2];

			Node<3>::Pointer pnode = mrModelPart.CreateNewNode(id, x, y, z);
			pnode->Set(NEW_ENTITY); //not boundary
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

			Node<3>::DofsContainerType &reference_dofs = NewDofs[nn];

			for (Node<3>::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
			{
				Node<3>::DofType &rDof = **iii;
				pnode->pAddDof(rDof);
			}

			Node<3>::Pointer SlaveNode1 = mrModelPart.pGetNode(NodesIDToInterpolate[nn][0]);
			Node<3>::Pointer SlaveNode2 = mrModelPart.pGetNode(NodesIDToInterpolate[nn][1]);
			Node<3>::Pointer SlaveNode3 = mrModelPart.pGetNode(NodesIDToInterpolate[nn][2]);

			InterpolateFromTwoNodes(pnode, SlaveNode1, SlaveNode2, VariablesList);

			TakeMaterialPropertiesFromNotRigidNode(pnode, SlaveNode3);
		}

		//set the coordinates to the original value
		const array_1d<double, 3> ZeroNormal(3, 0.0);
		for (std::vector<Node<3>::Pointer>::iterator it = list_of_new_nodes.begin(); it != list_of_new_nodes.end(); it++)
		{
			const array_1d<double, 3> &displacement = (*it)->FastGetSolutionStepValue(DISPLACEMENT);
			(*it)->X0() = (*it)->X() - displacement[0];
			(*it)->Y0() = (*it)->Y() - displacement[1];
			(*it)->Z0() = (*it)->Z() - displacement[2];

			(*it)->Set(FLUID);
			(*it)->Set(ACTIVE);
			(*it)->Reset(TO_ERASE);
			//correct contact_normal interpolation
			if ((*it)->SolutionStepsDataHas(CONTACT_FORCE))
				noalias((*it)->GetSolutionStepValue(CONTACT_FORCE)) = ZeroNormal;
		}

		KRATOS_CATCH("")
	}

	void CreateAndAddNewNodes(std::vector<array_1d<double, 3>> &NewPositions,
							  std::vector<array_1d<unsigned int, 4>> &NodesIDToInterpolate,
							  std::vector<Node<3>::DofsContainerType> &NewDofs,
							  int ElementsToRefine,
							  unsigned int &maxId)
	{
		KRATOS_TRY

		const unsigned int dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

		std::vector<Node<3>::Pointer> list_of_new_nodes;
		double NodeIdParent = MesherUtilities::GetMaxNodeId(*(mrModelPart.GetParentModelPart()));
		double NodeId = MesherUtilities::GetMaxNodeId(mrModelPart);

		unsigned int initial_node_size = NodeIdParent + 1 + ElementsToRefine; //total model part node size

		if (NodeId > NodeIdParent)
		{
			initial_node_size = NodeId + 1 + ElementsToRefine;
			std::cout << "initial_node_size  " << initial_node_size << std::endl;
		}

		//assign data to dofs
		VariablesList &VariablesList = mrModelPart.GetNodalSolutionStepVariablesList();

		for (unsigned int nn = 0; nn < NewPositions.size(); nn++)
		{

			unsigned int id = initial_node_size + nn;
			maxId = id;
			double x = NewPositions[nn][0];
			double y = NewPositions[nn][1];
			double z = 0;
			if (dimension == 3)
				z = NewPositions[nn][2];

			Node<3>::Pointer pnode = mrModelPart.CreateNewNode(id, x, y, z);
			pnode->Set(NEW_ENTITY); //not boundary
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

			// Node<3>::DofsContainerType& reference_dofs = (mrModelPart.NodesBegin())->GetDofs();
			Node<3>::DofsContainerType &reference_dofs = NewDofs[nn];

			for (Node<3>::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
			{
				Node<3>::DofType &rDof = **iii;
				// Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof(rDof);
				// (p_new_dof)->FreeDof();
				pnode->pAddDof(rDof);
			}

			Node<3>::Pointer SlaveNode1 = mrModelPart.pGetNode(NodesIDToInterpolate[nn][0]);
			Node<3>::Pointer SlaveNode2 = mrModelPart.pGetNode(NodesIDToInterpolate[nn][1]);
			InterpolateFromTwoNodes(pnode, SlaveNode1, SlaveNode2, VariablesList);
			if (SlaveNode1->Is(RIGID) || SlaveNode1->Is(SOLID))
			{
				TakeMaterialPropertiesFromNotRigidNode(pnode, SlaveNode2);
			}
			else {
                // Master node's properties are set using the maximum PROPERTY_ID value between SlaveNode1 and SlaveNode2
            	if (SlaveNode1->FastGetSolutionStepValue(PROPERTY_ID) >= SlaveNode2->FastGetSolutionStepValue(PROPERTY_ID)) {
                    TakeMaterialPropertiesFromNotRigidNode(pnode, SlaveNode1);
                } else {
                    TakeMaterialPropertiesFromNotRigidNode(pnode, SlaveNode2);
                }
			}
			// if(SlaveNode2->Is(RIGID) || SlaveNode2->Is(SOLID)){
			//   TakeMaterialPropertiesFromNotRigidNode(pnode,SlaveNode1);
			// }
        }

		//set the coordinates to the original value
		const array_1d<double, 3> ZeroNormal(3, 0.0);
		for (std::vector<Node<3>::Pointer>::iterator it = list_of_new_nodes.begin(); it != list_of_new_nodes.end(); it++)
		{
			const array_1d<double, 3> &displacement = (*it)->FastGetSolutionStepValue(DISPLACEMENT);
			(*it)->X0() = (*it)->X() - displacement[0];
			(*it)->Y0() = (*it)->Y() - displacement[1];
			(*it)->Z0() = (*it)->Z() - displacement[2];

			(*it)->Set(FLUID);
			(*it)->Set(ACTIVE);
			(*it)->Reset(TO_ERASE);
			// std::cout<<"velocity_x "<<(*it)->FastGetSolutionStepValue(VELOCITY_X,0)<<std::endl;
			// std::cout<<"velocity_x "<<(*it)->FastGetSolutionStepValue(VELOCITY_X,1)<<std::endl;
			// std::cout<<"velocity_x "<<(*it)->FastGetSolutionStepValue(VELOCITY_X,2)<<std::endl;
			// std::cout<<"pressure "<<(*it)->FastGetSolutionStepValue(PRESSURE,0)<<std::endl;
			// std::cout<<"pressure "<<(*it)->FastGetSolutionStepValue(PRESSURE,1)<<std::endl;
			// std::cout<<"pressure "<<(*it)->FastGetSolutionStepValue(PRESSURE,2)<<std::endl;
			// std::cout<<"acc "<<(*it)->FastGetSolutionStepValue(ACCELERATION_X,0)<<std::endl;
			// std::cout<<"acc "<<(*it)->FastGetSolutionStepValue(ACCELERATION_X,1)<<std::endl;
			// std::cout<<"acc "<<(*it)->FastGetSolutionStepValue(ACCELERATION_X,2)<<std::endl;
			// std::cout<<"bulkModulus "<<(*it)->FastGetSolutionStepValue(BULK_MODULUS)<<std::endl;
			// std::cout<<"density "<<(*it)->FastGetSolutionStepValue(DENSITY)<<std::endl;
			// std::cout<<"viscosity "<<(*it)->FastGetSolutionStepValue(DYNAMIC_VISCOSITY)<<std::endl;
			//correct contact_normal interpolation
			if ((*it)->SolutionStepsDataHas(CONTACT_FORCE))
				noalias((*it)->GetSolutionStepValue(CONTACT_FORCE)) = ZeroNormal;
		}

		KRATOS_CATCH("")
	}

	void InterpolateFromTwoNodes(Node<3>::Pointer MasterNode, Node<3>::Pointer SlaveNode1, Node<3>::Pointer SlaveNode2, VariablesList &rVariablesList)
	{

		KRATOS_TRY

		unsigned int buffer_size = MasterNode->GetBufferSize();

		for (VariablesList::const_iterator i_variable = rVariablesList.begin(); i_variable != rVariablesList.end(); i_variable++)
		{
			std::string variable_name = i_variable->Name();
			if (KratosComponents<Variable<double>>::Has(variable_name))
			{
				const Variable<double> & variable = KratosComponents<Variable<double>>::Get(variable_name);
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
				const Variable<array_1d<double, 3>> & variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(variable_name);
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
				const Variable<Matrix> & variable = KratosComponents<Variable<Matrix>>::Get(variable_name);
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
				const Variable<Vector> & variable = KratosComponents<Variable<Vector>>::Get(variable_name);
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

	void TakeMaterialPropertiesFromNotRigidNode(Node<3>::Pointer MasterNode, Node<3>::Pointer SlaveNode)
	{

		KRATOS_TRY
        MasterNode->FastGetSolutionStepValue(PROPERTY_ID) = SlaveNode->FastGetSolutionStepValue(PROPERTY_ID);
		MasterNode->FastGetSolutionStepValue(BULK_MODULUS) = SlaveNode->FastGetSolutionStepValue(BULK_MODULUS);
		MasterNode->FastGetSolutionStepValue(DENSITY) = SlaveNode->FastGetSolutionStepValue(DENSITY);
		MasterNode->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = SlaveNode->FastGetSolutionStepValue(DYNAMIC_VISCOSITY);

		MasterNode->FastGetSolutionStepValue(YIELD_SHEAR) = SlaveNode->FastGetSolutionStepValue(YIELD_SHEAR);

		MasterNode->FastGetSolutionStepValue(FLOW_INDEX) = SlaveNode->FastGetSolutionStepValue(FLOW_INDEX);
		MasterNode->FastGetSolutionStepValue(ADAPTIVE_EXPONENT) = SlaveNode->FastGetSolutionStepValue(ADAPTIVE_EXPONENT);
		MasterNode->FastGetSolutionStepValue(STATIC_FRICTION) = SlaveNode->FastGetSolutionStepValue(STATIC_FRICTION);
		MasterNode->FastGetSolutionStepValue(DYNAMIC_FRICTION) = SlaveNode->FastGetSolutionStepValue(DYNAMIC_FRICTION);
		MasterNode->FastGetSolutionStepValue(INERTIAL_NUMBER_ZERO) = SlaveNode->FastGetSolutionStepValue(INERTIAL_NUMBER_ZERO);
		MasterNode->FastGetSolutionStepValue(GRAIN_DIAMETER) = SlaveNode->FastGetSolutionStepValue(GRAIN_DIAMETER);
		MasterNode->FastGetSolutionStepValue(GRAIN_DENSITY) = SlaveNode->FastGetSolutionStepValue(GRAIN_DENSITY);
		MasterNode->FastGetSolutionStepValue(REGULARIZATION_COEFFICIENT) = SlaveNode->FastGetSolutionStepValue(REGULARIZATION_COEFFICIENT);

		if (MasterNode->SolutionStepsDataHas(DEVIATORIC_COEFFICIENT) && SlaveNode->SolutionStepsDataHas(DEVIATORIC_COEFFICIENT))
		{
			MasterNode->FastGetSolutionStepValue(DEVIATORIC_COEFFICIENT) = SlaveNode->FastGetSolutionStepValue(DEVIATORIC_COEFFICIENT);
			MasterNode->FastGetSolutionStepValue(VOLUMETRIC_COEFFICIENT) = SlaveNode->FastGetSolutionStepValue(VOLUMETRIC_COEFFICIENT);
		}

		if (MasterNode->SolutionStepsDataHas(YOUNG_MODULUS) && SlaveNode->SolutionStepsDataHas(YOUNG_MODULUS))
		{
			MasterNode->FastGetSolutionStepValue(YOUNG_MODULUS) = 0;
			MasterNode->FastGetSolutionStepValue(POISSON_RATIO) = 0;
		}
		if (MasterNode->SolutionStepsDataHas(SOLID_DENSITY) && SlaveNode->SolutionStepsDataHas(SOLID_DENSITY))
		{
			MasterNode->FastGetSolutionStepValue(SOLID_DENSITY) = 0;
			MasterNode->Reset(SOLID);
			MasterNode->FastGetSolutionStepValue(INTERFACE_NODE) = false;
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
