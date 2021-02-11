//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:                 AFranci $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:             eptember 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_BUILD_MODEL_PART_BOUNDARY_FOR_FLUIDS_PROCESS_H_INCLUDED)
#define KRATOS_BUILD_MODEL_PART_BOUNDARY_FOR_FLUIDS_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"

#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/line_3d_2.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"

#include "custom_conditions/composite_condition.hpp"
#include "custom_utilities/boundary_normals_calculation_utilities.hpp"
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_processes/mesher_process.hpp"

#include "delaunay_meshing_application_variables.h"

///VARIABLES used:
//Data:     MASTER_ELEMENTS(set), MASTER_NODES(set), NEIGHBOUR_ELEMENTS
//Flags:    (checked) CONTACT
//          (set)     BOUNDARY(set)
//          (modified)
//          (reset)
//(set):=(set in this process)

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
typedef ModelPart::NodesContainerType NodesContainerType;
typedef ModelPart::ElementsContainerType ElementsContainerType;
typedef ModelPart::ConditionsContainerType ConditionsContainerType;

typedef GlobalPointersVector<Node<3>> NodeWeakPtrVectorType;
typedef GlobalPointersVector<Element> ElementWeakPtrVectorType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
   */
class BuildModelPartBoundaryForFluidsProcess
	: public MesherProcess
{
public:
	///@name Type Definitions
	///@{

	/// Pointer definition of BuildModelPartBoundaryForFluidsProcess
	KRATOS_CLASS_POINTER_DEFINITION(BuildModelPartBoundaryForFluidsProcess);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	BuildModelPartBoundaryForFluidsProcess(ModelPart &rModelPart,
										   std::string const rModelPartName,
										   int EchoLevel = 0)
		: mrModelPart(rModelPart)
	{
		mModelPartName = rModelPartName;
		mEchoLevel = EchoLevel;
	}

	/// Destructor.
	virtual ~BuildModelPartBoundaryForFluidsProcess()
	{
	}

	///@}
	///@name Operators
	///@{

	void operator()()
	{
		Execute();
	}

	///@}
	///@name Operations
	///@{

	void Execute() override
	{

		KRATOS_TRY

		bool success = false;

		double begin_time = OpenMPUtils::GetCurrentTime();

		unsigned int NumberOfSubModelParts = mrModelPart.NumberOfSubModelParts();

		this->ResetNodesBoundaryFlag(mrModelPart);

		if (mModelPartName == mrModelPart.Name())
		{

			for (ModelPart::SubModelPartIterator i_mp = mrModelPart.SubModelPartsBegin(); i_mp != mrModelPart.SubModelPartsEnd(); ++i_mp)
			{

				if (mEchoLevel >= 1)
					std::cout << " [ Construct Boundary on ModelPart [" << i_mp->Name() << "] ]" << std::endl;

				success = UniqueSkinSearch(*i_mp);

				if (!success)
				{
					std::cout << "  ERROR: BOUNDARY CONSTRUCTION FAILED ModelPart : [" << i_mp->Name() << "] " << std::endl;
				}
				else
				{
					if (mEchoLevel >= 1)
					{
						double end_time = OpenMPUtils::GetCurrentTime();
						std::cout << " [ Performed in Time = " << end_time - begin_time << " ]" << std::endl;
					}
					//PrintSkin(*i_mp);
				}
			}
		}
		else
		{

			if (mEchoLevel >= 1)
				std::cout << " [ Construct Boundary on ModelPart[" << mModelPartName << "] ]" << std::endl;

			ModelPart &rModelPart = mrModelPart.GetSubModelPart(mModelPartName);
			success = UniqueSkinSearch(rModelPart);

			if (!success)
			{
				std::cout << "  ERROR: BOUNDARY CONSTRUCTION FAILED on ModelPart : [" << rModelPart.Name() << "] " << std::endl;
			}
			else
			{
				if (mEchoLevel >= 1)
				{
					double end_time = OpenMPUtils::GetCurrentTime();
					std::cout << " [ Performed in Time = " << end_time - begin_time << " ]" << std::endl;
				}
				//PrintSkin(rModelPart);
			}
		}

		if (NumberOfSubModelParts > 1)
		{
			SetMainModelPartConditions();
			SetComputingModelPart();
		}

		//ComputeBoundaryNormals BoundUtils;
		BoundaryNormalsCalculationUtilities BoundaryComputation;
		if (mModelPartName == mrModelPart.Name())
		{
			BoundaryComputation.CalculateWeightedBoundaryNormals(mrModelPart, mEchoLevel);
		}
		else
		{
			ModelPart &rModelPart = mrModelPart.GetSubModelPart(mModelPartName);
			BoundaryComputation.CalculateWeightedBoundaryNormals(rModelPart, mEchoLevel);
		}

		if (mEchoLevel >= 1)
			std::cout << "  Boundary Normals Computed and Assigned ] " << std::endl;

		KRATOS_CATCH("")
	}

	//**************************************************************************
	//**************************************************************************

	bool SearchConditionMasters()
	{

		KRATOS_TRY

		int composite_conditions = 0;
		int total_conditions = 0;
		int counter = 0;

		bool found = false;

		for (ModelPart::ConditionsContainerType::iterator i_cond = mrModelPart.ConditionsBegin(); i_cond != mrModelPart.ConditionsEnd(); ++i_cond)
		{

			if (i_cond->Is(BOUNDARY)) //composite condition
				composite_conditions++;

			//std::cout<<" BeforeSearch::Condition ("<<i_cond->Id()<<") ME="<<i_cond->GetValue(MASTER_ELEMENTS)[0]->Id()<<", MN= "<<i_cond->GetValue(MASTER_NODES)[0]->Id()<<std::endl;

			//********************************************************************

			DenseMatrix<unsigned int> lpofa; //connectivities of points defining faces
			DenseVector<unsigned int> lnofa; //number of points defining faces

			Geometry<Node<3>> &rConditionGeometry = i_cond->GetGeometry();
			unsigned int size = rConditionGeometry.size();

			bool perform_search = true;
			for (unsigned int i = 0; i < size; ++i)
			{
				if (rConditionGeometry[i].Is(RIGID)) //if is a rigid wall do not search else do search
					perform_search = false;
			}

			if (i_cond->Is(CONTACT))
				perform_search = false;

			//********************************************************************
			found = false;

			if (perform_search)
			{

				if (size == 2)
				{

					ElementWeakPtrVectorType &rE1 = rConditionGeometry[0].GetValue(NEIGHBOUR_ELEMENTS);
					ElementWeakPtrVectorType &rE2 = rConditionGeometry[1].GetValue(NEIGHBOUR_ELEMENTS);

					if (rE1.size() == 0 || rE2.size() == 0)
						std::cout << " NO SIZE in NEIGHBOUR_ELEMENTS " << std::endl;

					for (ElementWeakPtrVectorType::iterator ie = rE1.begin(); ie != rE1.end(); ++ie)
					{
						for (ElementWeakPtrVectorType::iterator ne = rE2.begin(); ne != rE2.end(); ++ne)
						{

							if ((ne)->Id() == (ie)->Id() && !found)
							{
								ElementWeakPtrVectorType MasterElements;
								MasterElements.push_back(*ie.base());
								if (mEchoLevel >= 1)
								{
									//if(i_cond->GetValue(MASTER_ELEMENTS)[0]->Id() != MasterElements[0]->Id())
									//std::cout<<"Condition "<<i_cond->Id()<<" WARNING: master elements ("<<i_cond->GetValue(MASTER_ELEMENTS)[0]->Id()<<" != "<<MasterElements[0]->Id()<<")"<<std::endl;
								}
								i_cond->SetValue(MASTER_ELEMENTS, MasterElements);

								Geometry<Node<3>> &rElementGeometry = (ie)->GetGeometry();

								//get matrix nodes in faces
								rElementGeometry.NodesInFaces(lpofa);
								rElementGeometry.NumberNodesInFaces(lnofa);

								int node = 0;
								for (unsigned int iface = 0; iface < rElementGeometry.size(); ++iface)
								{
									MesherUtilities MesherUtils;
									found = MesherUtils.FindCondition(rConditionGeometry, rElementGeometry, lpofa, lnofa, iface);

									if (found)
									{
										node = iface;
										break;
									}
								}

								if (found)
								{
									NodeWeakPtrVectorType MasterNodes;
									MasterNodes.push_back(rElementGeometry(lpofa(0, node)));
									if (mEchoLevel >= 1)
									{
										if (i_cond->GetValue(MASTER_NODES)[0].Id() != MasterNodes[0].Id())
											std::cout << "Condition " << i_cond->Id() << " WARNING: master nodes (" << i_cond->GetValue(MASTER_NODES)[0].Id() << " != " << MasterNodes[0].Id() << ")" << std::endl;
										i_cond->SetValue(MASTER_NODES, MasterNodes);
									}
								}
								else
								{
									std::cout << " MASTER_NODE not FOUND : something is wrong " << std::endl;
								}
							}
						}
					}
				}
				if (size == 3)
				{

					ElementWeakPtrVectorType &rE1 = rConditionGeometry[0].GetValue(NEIGHBOUR_ELEMENTS);
					ElementWeakPtrVectorType &rE2 = rConditionGeometry[1].GetValue(NEIGHBOUR_ELEMENTS);
					ElementWeakPtrVectorType &rE3 = rConditionGeometry[2].GetValue(NEIGHBOUR_ELEMENTS);

					if (rE1.size() == 0 || rE2.size() == 0 || rE3.size() == 0)
						std::cout << " NO SIZE in NEIGHBOUR_ELEMENTS " << std::endl;

					for (ElementWeakPtrVectorType::iterator ie = rE1.begin(); ie != rE1.end(); ++ie)
					{
						for (ElementWeakPtrVectorType::iterator je = rE2.begin(); je != rE2.end(); ++je)
						{

							if ((je)->Id() == (ie)->Id() && !found)
							{

								for (ElementWeakPtrVectorType::iterator ke = rE3.begin(); ke != rE3.end(); ++ke)
								{

									if ((ke)->Id() == (ie)->Id() && !found)
									{

										ElementWeakPtrVectorType MasterElements;
										MasterElements.push_back(*ie.base());
										if (mEchoLevel >= 1)
										{
											if (i_cond->GetValue(MASTER_ELEMENTS)[0].Id() != MasterElements[0].Id())
												std::cout << "Condition " << i_cond->Id() << " WARNING: master elements (" << i_cond->GetValue(MASTER_ELEMENTS)[0].Id() << " != " << MasterElements[0].Id() << ")" << std::endl;
										}
										i_cond->SetValue(MASTER_ELEMENTS, MasterElements);

										Geometry<Node<3>> &rElementGeometry = (ie)->GetGeometry();

										//get matrix nodes in faces
										rElementGeometry.NodesInFaces(lpofa);
										rElementGeometry.NumberNodesInFaces(lnofa);

										int node = 0;
										for (unsigned int iface = 0; iface < rElementGeometry.size(); ++iface)
										{
											MesherUtilities MesherUtils;
											found = MesherUtils.FindCondition(rConditionGeometry, rElementGeometry, lpofa, lnofa, iface);

											if (found)
											{
												node = iface;
												break;
											}
										}

										if (found)
										{
											NodeWeakPtrVectorType MasterNodes;
											MasterNodes.push_back(rElementGeometry(lpofa(0, node)));
											if (mEchoLevel >= 1)
											{
												if (i_cond->GetValue(MASTER_NODES)[0].Id() != MasterNodes[0].Id())
													std::cout << "Condition " << i_cond->Id() << " WARNING: master nodes (" << i_cond->GetValue(MASTER_NODES)[0].Id() << " != " << MasterNodes[0].Id() << ")" << std::endl;
											}
											i_cond->SetValue(MASTER_NODES, MasterNodes);
										}
										else
										{
											std::cout << " MASTER_NODE not FOUND : something is wrong " << std::endl;
										}
									}
								}
							}
						}
					}
				}

				total_conditions++;
			}

			//********************************************************************

			//std::cout<<" AfterSearch::Condition ("<<i_cond->Id()<<") : ME="<<i_cond->GetValue(MASTER_ELEMENTS)[0].Id()<<", MN= "<<i_cond->GetValue(MASTER_NODES)[0].Id()<<std::endl;

			if (found)
				counter++;
		}

		if (counter == total_conditions)
		{
			if (mEchoLevel >= 1)
				std::cout << "   Condition Masters (ModelPart " << mrModelPart.Name() << "): LOCATED [" << counter << "]" << std::endl;
			found = true;
		}
		else
		{
			if (mEchoLevel >= 1)
				std::cout << "   Condition Masters (ModelPart " << mrModelPart.Name() << "): not LOCATED [" << counter - total_conditions << "]" << std::endl;
			found = false;
		}

		if (counter != composite_conditions)
			if (mEchoLevel >= 1)
				std::cout << "   Condition Masters (ModelPart " << mrModelPart.Name() << "): LOCATED [" << counter << "] COMPOSITE [" << composite_conditions << "] NO MATCH" << std::endl;

		return found;

		std::cout << " Condition Masters Found " << std::endl;

		KRATOS_CATCH("")
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
		return "BuildModelPartBoundaryForFluidsProcess";
	}

	/// Print information about this object.
	void PrintInfo(std::ostream &rOStream) const override
	{
		rOStream << "BuildModelPartBoundaryForFluidsProcess";
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
	///@name Protected static Member Variables
	///@{

	///@}
	///@name Protected member Variables
	///@{

	ModelPart &mrModelPart;

	std::string mModelPartName;

	int mEchoLevel;

	///@}
	///@name Protected Operators
	///@{

	///@}
	///@name Protected Operations
	///@{

	//**************************************************************************
	//**************************************************************************

	bool ClearMasterEntities(ModelPart &rModelPart, ModelPart::ConditionsContainerType &rTemporaryConditions)
	{
		KRATOS_TRY

		for (ModelPart::ConditionsContainerType::iterator ic = rTemporaryConditions.begin(); ic != rTemporaryConditions.end(); ++ic)
		{
			ElementWeakPtrVectorType &MasterElements = ic->GetValue(MASTER_ELEMENTS);
			MasterElements.erase(MasterElements.begin(), MasterElements.end());

			NodeWeakPtrVectorType &MasterNodes = ic->GetValue(MASTER_NODES);
			MasterNodes.erase(MasterNodes.begin(), MasterNodes.end());
		}

		return true;

		KRATOS_CATCH("")
	}

	//**************************************************************************
	//**************************************************************************

	bool UniqueSkinSearch(ModelPart &rModelPart)
	{

		KRATOS_TRY

		if (mEchoLevel > 0)
		{
			std::cout << " [ Initial Conditions : " << rModelPart.Conditions().size() << std::endl;
		}

		if (!rModelPart.Elements().size() || (rModelPart.Is(ACTIVE)))
		{
			if (mEchoLevel > 0)
			{
				std::cout << " [ Final Conditions   : " << rModelPart.Conditions().size() << std::endl;
			}
			return true;
		}

		//check if a remesh process has been performed and there is any node to erase
		bool any_node_to_erase = false;
		for (ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(); in != rModelPart.NodesEnd(); ++in)
		{
			if (any_node_to_erase == false)
				if (in->Is(TO_ERASE))
					any_node_to_erase = true;
		}

		this->SetBoundaryAndFreeSurface(rModelPart);
		// //swap conditions for a temporary use
		// unsigned int ConditionId=1;
		// ModelPart::ConditionsContainerType TemporaryConditions;

		// //if there are no conditions check main modelpart mesh conditions
		// if( !rModelPart.Conditions().size() ){

		// 	for(ModelPart::ConditionsContainerType::iterator i_cond = rModelPart.GetParentModelPart().ConditionsBegin(); i_cond!= rModelPart.GetParentModelPart().ConditionsEnd(); ++i_cond)
		// 	  {
		// 	    TemporaryConditions.push_back(*(i_cond.base()));
		// 	    i_cond->SetId(ConditionId);
		// 	    ConditionId++;
		// 	  }

		// }
		// else{

		// 	TemporaryConditions.reserve(rModelPart.Conditions().size());
		// 	TemporaryConditions.swap(rModelPart.Conditions());

		// 	//set consecutive ids in the mesh conditions
		// 	if( any_node_to_erase ){
		// 	  for(ModelPart::ConditionsContainerType::iterator i_cond = TemporaryConditions.begin(); i_cond!= TemporaryConditions.end(); ++i_cond)
		// 	    {
		// 	      Geometry< Node<3> >& rConditionGeometry = i_cond->GetGeometry();
		// 	      for( unsigned int i=0; i<rConditionGeometry.size(); i++ )
		// 		{
		// 		  if( rConditionGeometry[i].Is(TO_ERASE)){
		// 		    i_cond->Set(TO_ERASE);
		// 		    break;
		// 		  }
		// 		}

		// 	      i_cond->SetId(ConditionId);
		// 	      ConditionId++;
		// 	    }
		// 	}
		// 	else{
		// 	  for(ModelPart::ConditionsContainerType::iterator i_cond = TemporaryConditions.begin(); i_cond!= TemporaryConditions.end(); ++i_cond)
		// 	    {

		// 	      i_cond->SetId(ConditionId);
		// 	      ConditionId++;
		// 	    }
		// 	}

		// }

		// //control the previous mesh conditions
		// std::vector<int> PreservedConditions( TemporaryConditions.size() + 1 );
		// std::fill( PreservedConditions.begin(), PreservedConditions.end(), 0 );

		// //build new skin for the Modelpart
		// this->BuildCompositeConditions(rModelPart, TemporaryConditions, PreservedConditions, ConditionId);

		// //add other conditions out of the skin space dimension
		// this->AddOtherConditions(rModelPart, TemporaryConditions, PreservedConditions, ConditionId);

		return true;

		KRATOS_CATCH("")
	}

	//**************************************************************************
	//**************************************************************************

	virtual bool SetBoundaryAndFreeSurface(ModelPart &rModelPart)
	{

		KRATOS_TRY

		//properties to be used in the generation
		int number_properties = rModelPart.GetParentModelPart().NumberOfProperties();
		Properties::Pointer properties = rModelPart.GetParentModelPart().pGetProperties(number_properties - 1);

		ModelPart::ElementsContainerType::iterator elements_begin = rModelPart.ElementsBegin();
		ModelPart::ElementsContainerType::iterator elements_end = rModelPart.ElementsEnd();

		for (ModelPart::ElementsContainerType::iterator ie = elements_begin; ie != elements_end; ie++)
		{
			Geometry<Node<3>> &rElementGeometry = ie->GetGeometry();

			if (rElementGeometry.FacesNumber() >= 3)
			{ //3 or 4

				//********************************************************************
				/*each face is opposite to the corresponding node number so in 2D triangle
	      0 ----- 1 2
	      1 ----- 2 0
	      2 ----- 0 1
	    */

				/*each face is opposite to the corresponding node number so in 3D tetrahedron
	      0 ----- 1 2 3
	      1 ----- 2 0 3
	      2 ----- 0 1 3
	      3 ----- 0 2 1
	    */
				//********************************************************************

				//finding boundaries and creating the "skin"
				boost::numeric::ublas::matrix<unsigned int> lpofa; //connectivities of points defining faces
				boost::numeric::ublas::vector<unsigned int> lnofa; //number of points defining faces

				ElementWeakPtrVectorType &rE = ie->GetValue(NEIGHBOUR_ELEMENTS);

				//get matrix nodes in faces
				rElementGeometry.NodesInFaces(lpofa);
				rElementGeometry.NumberNodesInFaces(lnofa);

				//loop on neighbour elements of an element
				unsigned int iface = 0;
				for (ElementWeakPtrVectorType::iterator ne = rE.begin(); ne != rE.end(); ne++)
				{
					unsigned int NumberNodesInFace = lnofa[iface];

					if ((ne)->Id() == ie->Id())
					{
						//if no neighbour is present => the face is free surface
						bool freeSurfaceFace = false;
						for (unsigned int j = 1; j <= NumberNodesInFace; j++)
						{
							rElementGeometry[lpofa(j, iface)].Set(BOUNDARY);
							if (rElementGeometry[lpofa(j, iface)].IsNot(RIGID))
							{
								freeSurfaceFace = true;
							}
						}
						if (freeSurfaceFace == true)
						{
							for (unsigned int j = 1; j <= NumberNodesInFace; j++)
							{
								rElementGeometry[lpofa(j, iface)].Set(FREE_SURFACE);
							}
						}

					} //end face condition

					iface += 1;
				} //end loop neighbours
			}
		}

		return true;

		KRATOS_CATCH("")
	}

	virtual bool BuildCompositeConditions(ModelPart &rModelPart, ModelPart::ConditionsContainerType &rTemporaryConditions, std::vector<int> &rPreservedConditions, unsigned int &rConditionId)
	{

		KRATOS_TRY

		//master conditions must be deleted and set them again in the build
		this->ClearMasterEntities(rModelPart, rTemporaryConditions);

		//properties to be used in the generation
		int number_properties = rModelPart.GetParentModelPart().NumberOfProperties();
		Properties::Pointer properties = rModelPart.GetParentModelPart().pGetProperties(number_properties - 1);

		ModelPart::ElementsContainerType::iterator elements_begin = rModelPart.ElementsBegin();
		ModelPart::ElementsContainerType::iterator elements_end = rModelPart.ElementsEnd();

		//clear nodal boundary flag
		for (ModelPart::ElementsContainerType::iterator ie = elements_begin; ie != elements_end; ++ie)
		{
			Geometry<Node<3>> &rElementGeometry = ie->GetGeometry();

			for (unsigned int j = 0; j < rElementGeometry.size(); ++j)
			{
				rElementGeometry[j].Reset(BOUNDARY);
			}
		}

		rConditionId = 0;
		for (ModelPart::ElementsContainerType::iterator ie = elements_begin; ie != elements_end; ++ie)
		{
			Geometry<Node<3>> &rElementGeometry = ie->GetGeometry();

			const unsigned int dimension = rElementGeometry.WorkingSpaceDimension();

			if (rElementGeometry.FacesNumber() >= (dimension + 1))
			{ //3 or 4

				//********************************************************************
				/*each face is opposite to the corresponding node number so in 2D triangle
	      0 ----- 1 2
	      1 ----- 2 0
	      2 ----- 0 1
	    */

				/*each face is opposite to the corresponding node number so in 3D tetrahedron
	      0 ----- 1 2 3
	      1 ----- 2 0 3
	      2 ----- 0 1 3
	      3 ----- 0 2 1
	    */
				//********************************************************************

				//finding boundaries and creating the "skin"
				DenseMatrix<unsigned int> lpofa; //connectivities of points defining faces
				DenseVector<unsigned int> lnofa; //number of points defining faces

				ElementWeakPtrVectorType &rE = ie->GetValue(NEIGHBOUR_ELEMENTS);

				//get matrix nodes in faces
				rElementGeometry.NodesInFaces(lpofa);
				rElementGeometry.NumberNodesInFaces(lnofa);

				//loop on neighbour elements of an element
				unsigned int iface = 0;
				for (ElementWeakPtrVectorType::iterator ne = rE.begin(); ne != rE.end(); ++ne)
				{
					unsigned int NumberNodesInFace = lnofa[iface];

					if ((ne)->Id() == ie->Id())
					{
						//if no neighbour is present => the face is free surface
						for (unsigned int j = 1; j <= NumberNodesInFace; ++j)
						{
							rElementGeometry[lpofa(j, iface)].Set(BOUNDARY);
							//std::cout<<" node ["<<j<<"]"<<rElementGeometry[lpofa(j,iface)].Id()<<std::endl;
						}

						//1.- create geometry: points array and geometry type
						Condition::NodesArrayType FaceNodes;
						Condition::GeometryType::Pointer ConditionVertices;

						FaceNodes.reserve(NumberNodesInFace);

						for (unsigned int j = 1; j <= NumberNodesInFace; ++j)
						{
							FaceNodes.push_back(rElementGeometry(lpofa(j, iface)));
						}

						if (NumberNodesInFace == 2)
						{
							if (dimension == 2)
								ConditionVertices = Kratos::make_shared<Line2D2<Node<3>>>(FaceNodes);
							else
								ConditionVertices = Kratos::make_shared<Line3D2<Node<3>>>(FaceNodes);
						}
						else if (NumberNodesInFace == 3)
						{
							if (dimension == 2)
								ConditionVertices = Kratos::make_shared<Line2D3<Node<3>>>(FaceNodes);
							else
								ConditionVertices = Kratos::make_shared<Triangle3D3<Node<3>>>(FaceNodes);
						}

						rConditionId += 1;

						//Create a composite condition
						CompositeCondition::Pointer p_cond = Kratos::make_intrusive<CompositeCondition>(rConditionId, ConditionVertices, properties);

						bool condition_found = false;
						bool point_condition = false;

						// Search for existing conditions: start
						for (ModelPart::ConditionsContainerType::iterator i_cond = rTemporaryConditions.begin(); i_cond != rTemporaryConditions.end(); ++i_cond)
						{
							Geometry<Node<3>> &rConditionGeometry = i_cond->GetGeometry();

							MesherUtilities MesherUtils;
							condition_found = MesherUtils.FindCondition(rConditionGeometry, rElementGeometry, lpofa, lnofa, iface);

							if (condition_found)
							{

								p_cond->AddChild(*(i_cond.base()));

								rPreservedConditions[i_cond->Id()] += 1;

								if (rConditionGeometry.PointsNumber() == 1)
									point_condition = true;
							}
						}
						// Search for existing conditions: end

						if (!point_condition)
						{
							// usually one MasterElement and one MasterNode for 2D and 3D simplex
							// can be more than one in other geometries -> it has to be extended to that cases

							// std::cout<<" ID "<<p_cond->Id()<<" MASTER ELEMENT "<<ie->Id()<<std::endl;
							// std::cout<<" MASTER NODE "<<rElementGeometry[lpofa(0,iface)].Id()<<" or "<<rElementGeometry[lpofa(NumberNodesInFace,iface)].Id()<<std::endl;

							ElementWeakPtrVectorType &MasterElements = p_cond->GetValue(MASTER_ELEMENTS);
							MasterElements.push_back((*(ie.base())));
							p_cond->SetValue(MASTER_ELEMENTS, MasterElements);

							NodeWeakPtrVectorType &MasterNodes = p_cond->GetValue(MASTER_NODES);
							MasterNodes.push_back(rElementGeometry(lpofa(0, iface)));
							p_cond->SetValue(MASTER_NODES, MasterNodes);
						}

						rModelPart.Conditions().push_back(Condition::Pointer(p_cond));
						// Set new conditions: end

					} //end face condition

					iface += 1;
				} //end loop neighbours
			}
			else
			{
				//set nodes to BOUNDARY for elements outside of the working space dimension
				for (unsigned int j = 0; j < rElementGeometry.size(); ++j)
				{
					rElementGeometry[j].Set(BOUNDARY);
				}
			}
		}

		return true;

		KRATOS_CATCH("")
	}

	//**************************************************************************
	//**************************************************************************

	bool FindConditionID(ModelPart &rModelPart, Geometry<Node<3>> &rConditionGeometry)
	{
		KRATOS_TRY

		//check if the condition exists and belongs to the modelpart checking node Ids
		for (unsigned int i = 0; i < rConditionGeometry.size(); ++i)
		{
			for (ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(); in != rModelPart.NodesEnd(); ++in)
			{
				if (rConditionGeometry[i].Id() == in->Id())
					return true;
			}
		}

		return false;

		KRATOS_CATCH("")
	}

	//**************************************************************************
	//**************************************************************************

	void PrintSkin(ModelPart &rModelPart)
	{

		KRATOS_TRY

		//PRINT SKIN:
		std::cout << " CONDITIONS: geometry nodes (" << rModelPart.Conditions().size() << ")" << std::endl;

		ConditionsContainerType &rCond = rModelPart.Conditions();
		for (ConditionsContainerType::iterator i_cond = rCond.begin(); i_cond != rCond.end(); ++i_cond)
		{

			Geometry<Node<3>> &rConditionGeometry = i_cond->GetGeometry();
			std::cout << "[" << i_cond->Id() << "]:" << std::endl;
			//i_cond->PrintInfo(std::cout);
			std::cout << "( ";
			for (unsigned int i = 0; i < rConditionGeometry.size(); ++i)
			{
				std::cout << rConditionGeometry[i].Id() << ", ";
			}
			std::cout << " ): ";

			i_cond->GetValue(MASTER_ELEMENTS)[0].PrintInfo(std::cout);

			std::cout << std::endl;
		}
		std::cout << std::endl;

		KRATOS_CATCH("")
	}

	//**************************************************************************
	//**************************************************************************

	bool AddOtherConditions(ModelPart &rModelPart, ModelPart::ConditionsContainerType &rTemporaryConditions, std::vector<int> &rPreservedConditions, unsigned int &rConditionId)
	{

		KRATOS_TRY

		unsigned int counter = 0;

		//add all previous conditions not found in the skin search:
		for (ModelPart::ConditionsContainerType::iterator i_cond = rTemporaryConditions.begin(); i_cond != rTemporaryConditions.end(); ++i_cond)
		{

			if (rPreservedConditions[i_cond->Id()] == 0)
			{ //I have not used the condition and any node of the condition

				if (this->CheckAcceptedCondition(rModelPart, *i_cond))
				{

					Geometry<Node<3>> &rConditionGeometry = i_cond->GetGeometry();

					Condition::NodesArrayType FaceNodes;

					FaceNodes.reserve(rConditionGeometry.size());

					for (unsigned int j = 0; j < rConditionGeometry.size(); ++j)
					{
						FaceNodes.push_back(rConditionGeometry(j));
					}

					rPreservedConditions[i_cond->Id()] += 1;

					rConditionId += 1;

					Condition::Pointer p_cond = i_cond->Clone(rConditionId, FaceNodes);
					//p_cond->Data() = i_cond->Data();

					this->AddConditionToModelPart(rModelPart, p_cond);

					counter++;
				}
			}
		}

		//std::cout<<"   recovered conditions "<<counter<<std::endl;

		//control if all previous conditions have been added:
		bool all_assigned = true;
		unsigned int lost_conditions = 0;
		for (unsigned int i = 1; i < rPreservedConditions.size(); ++i)
		{
			if (rPreservedConditions[i] == 0)
			{
				all_assigned = false;
				lost_conditions++;
			}
		}

		if (mEchoLevel >= 1)
		{

			std::cout << "   Final Conditions   : " << rModelPart.NumberOfConditions() << std::endl;

			if (all_assigned == true)
				std::cout << "   ALL_PREVIOUS_CONDITIONS_RELOCATED " << std::endl;
			else
				std::cout << "   SOME_PREVIOUS_CONDITIONS_ARE_LOST [lost_conditions:" << lost_conditions << "]" << std::endl;
		}

		return all_assigned;

		KRATOS_CATCH("")
	}

	//**************************************************************************
	//**************************************************************************

	virtual bool CheckAcceptedCondition(ModelPart &rModelPart, Condition &rCondition)
	{
		KRATOS_TRY

		Geometry<Node<3>> &rConditionGeometry = rCondition.GetGeometry();

		return FindConditionID(rModelPart, rConditionGeometry);

		KRATOS_CATCH("")
	}

	//**************************************************************************
	//**************************************************************************

	virtual void AddConditionToModelPart(ModelPart &rModelPart, Condition::Pointer pCondition)
	{
		KRATOS_TRY

		rModelPart.AddCondition(pCondition);
		//rModelPart.Conditions().push_back(pCondition);

		KRATOS_CATCH("")
	}

	///@}
	///@name Protected  Access
	///@{

	///@}
	///@name Protected Inquiry
	///@{

	///@}
	///@name Protected LifeCycle
	///@{

	///@}

private:
	///@name Static Member Variables
	///@{

	///@}
	///@name Member Variables
	///@{

	///@}
	///@name Private Operators
	///@{

	///@}
	///@name Private Operations
	///@{

	//**************************************************************************
	//**************************************************************************

	void SetComputingModelPart()
	{

		KRATOS_TRY

		std::string ComputingModelPartName;
		for (ModelPart::SubModelPartIterator i_mp = mrModelPart.SubModelPartsBegin(); i_mp != mrModelPart.SubModelPartsEnd(); ++i_mp)
		{
			if (i_mp->Is(ACTIVE))
				ComputingModelPartName = i_mp->Name();
		}

		ModelPart &rModelPart = mrModelPart.GetSubModelPart(ComputingModelPartName);

		if (mEchoLevel >= 1)
		{
			std::cout << " [" << ComputingModelPartName << " :: CONDITIONS [OLD:" << rModelPart.NumberOfConditions();
		}

		ModelPart::ConditionsContainerType KeepConditions;

		for (ModelPart::SubModelPartIterator i_mp = mrModelPart.SubModelPartsBegin(); i_mp != mrModelPart.SubModelPartsEnd(); ++i_mp)
		{

			if (i_mp->NumberOfElements() && ComputingModelPartName != i_mp->Name())
			{ //conditions of model_parts with elements only

				for (ModelPart::ConditionsContainerType::iterator i_cond = i_mp->ConditionsBegin(); i_cond != i_mp->ConditionsEnd(); ++i_cond)
				{
					// i_cond->PrintInfo(std::cout);
					// std::cout<<" -- "<<std::endl;

					KeepConditions.push_back(*(i_cond.base()));

					// KeepConditions.back().PrintInfo(std::cout);
					// std::cout<<std::endl;
				}
			}
		}

		rModelPart.Conditions().swap(KeepConditions);

		if (mEchoLevel >= 1)
		{
			std::cout << " / NEW:" << rModelPart.NumberOfConditions() << "] " << std::endl;
		}

		KRATOS_CATCH("")
	}

	//**************************************************************************
	//**************************************************************************

	void SetMainModelPartConditions()
	{

		KRATOS_TRY

		if (mEchoLevel >= 1)
		{
			std::cout << " [" << mrModelPart.Name() << " :: CONDITIONS [OLD:" << mrModelPart.NumberOfConditions();
		}

		//contact conditions are located on Mesh_0
		ModelPart::ConditionsContainerType KeepConditions;

		unsigned int condId = 1;
		if (mModelPartName == mrModelPart.Name())
		{

			for (ModelPart::SubModelPartIterator i_mp = mrModelPart.SubModelPartsBegin(); i_mp != mrModelPart.SubModelPartsEnd(); ++i_mp)
			{
				if (!(i_mp->Is(ACTIVE)) && !(i_mp->Is(CONTACT)))
				{
					//std::cout<<" ModelPartName "<<i_mp->Name()<<" conditions "<<i_mp->NumberOfConditions()<<std::endl;
					for (ModelPart::ConditionsContainerType::iterator i_cond = i_mp->ConditionsBegin(); i_cond != i_mp->ConditionsEnd(); ++i_cond)
					{
						// i_cond->PrintInfo(std::cout);
						// std::cout<<" -- "<<std::endl;

						KeepConditions.push_back(*(i_cond.base()));
						KeepConditions.back().SetId(condId);
						condId += 1;

						// KeepConditions.back().PrintInfo(std::cout);
						// std::cout<<std::endl;
					}
				}
			}
		}

		for (ModelPart::ConditionsContainerType::iterator i_cond = mrModelPart.ConditionsBegin(); i_cond != mrModelPart.ConditionsEnd(); ++i_cond)
		{
			if (i_cond->Is(CONTACT))
			{
				KeepConditions.push_back(*(i_cond.base()));
				KeepConditions.back().SetId(condId);
				condId += 1;

				//std::cout<<" -- "<<std::endl;
				//KeepConditions.back().PrintInfo(std::cout);
				//std::cout<<std::endl;
			}
		}

		mrModelPart.Conditions().swap(KeepConditions);

		if (mEchoLevel >= 1)
		{
			std::cout << " / NEW:" << mrModelPart.NumberOfConditions() << "] " << std::endl;
		}

		KRATOS_CATCH("")
	}

	void ResetNodesBoundaryFlag(ModelPart &rModelPart)
	{
		//reset the boundary flag in all nodes
		for (ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(); in != rModelPart.NodesEnd(); ++in)
		{
			in->Reset(BOUNDARY);
		}
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
	BuildModelPartBoundaryForFluidsProcess &operator=(BuildModelPartBoundaryForFluidsProcess const &rOther);

	/// Copy constructor.
	//BuildModelPartBoundaryForFluidsProcess(BuildModelPartBoundaryForFluidsProcess const& rOther);

	///@}

}; // Class BuildModelPartBoundaryForFluidsProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
								BuildModelPartBoundaryForFluidsProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
								const BuildModelPartBoundaryForFluidsProcess &rThis)
{
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_BUILD_MODEL_PART_BOUNDARY_FOR_FLUIDS_PROCESS_H_INCLUDED  defined
