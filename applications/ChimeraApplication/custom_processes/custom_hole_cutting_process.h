// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
//
// ==============================================================================
//

#if !defined(KRATOS_CUSTOM_HOLE_CUTTING_PROCESS_H_INCLUDED)
#define KRATOS_CUSTOM_HOLE_CUTTING_PROCESS_H_INCLUDED

// System includes
// Please put system includes in custom_hole_cutting_process.h

// External includes
// Please put external includes in custom_hole_cutting_process.h

// Project includes

// System includes
#include <iostream>
#include <string>
#include <algorithm>
#include "math.h"

// External includes
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/node.h"
#include "includes/model_part.h"
#include "includes/process_info.h"
// Application includes
#include "processes/process.h"
#include "processes/find_nodal_neighbours_process.h"	  // To find node neighbours using elements
#include "processes/find_conditions_neighbours_process.h" // To find node neighbours using conditions
#include "processes/node_erase_process.h"				  // To delete empty nodes
#include "utilities/normal_calculation_utils.h"			  // To calculate element's normal
#include "geometries/triangle_3d_3.h"					  // Skin face geometry template
#include "geometries/line_2d_2.h"

#include "custom_utilities/vtk_output.hpp"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

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

class CustomHoleCuttingProcess
{
  public:
	// Needed structures for the ExtractSurfaceMesh operation
	struct KeyComparor
	{
		bool operator()(const vector<std::size_t> &lhs, const vector<std::size_t> &rhs) const
		{
			if (lhs.size() != rhs.size())
				return false;

			for (std::size_t i = 0; i < lhs.size(); i++)
			{
				if (lhs[i] != rhs[i])
					return false;
			}

			return true;
		}
	};

	struct KeyHasher
	{
		std::size_t operator()(const vector<int> &k) const
		{
			std::size_t seed = 0.0;
			std::hash<int> hasher;

			for (std::size_t i = 0; i < k.size(); i++)
			{
				seed ^= hasher(k[i]) + 0x9e3779b9 + (seed<<6) + (seed>>2);
			}

			return seed;
		}

	};

/*



*/


	///@name Type Definitions
	///@{

	///@}
	///@name Pointer Definitions
	/// Pointer definition of CustomHoleCuttingProcess
	KRATOS_CLASS_POINTER_DEFINITION(CustomHoleCuttingProcess);

	///@}
	///@name Life Cycle
	///@{

	CustomHoleCuttingProcess()
	{
	}

	/// Destructor.
	virtual ~CustomHoleCuttingProcess()
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

	/// For CHIMERA boundary condition purposes: Extracts a  mesh with a certain threshold value
	void ExtractMeshBetweenLimits(ModelPart &rModelPart, ModelPart &rExtractedModelPart, double lLimit, double uLimit)
	{
		KRATOS_TRY;
		KRATOS_INFO("\n::[Mesh Extraction]::")<< std::endl;
		// Initializing mesh nodes
		// Extracting mesh elements which are only above the threshold value
		KRATOS_INFO("Extracting elements between")<< lLimit << " and " << uLimit << std::endl;
		Element::Pointer pElem;
		std::vector<std::size_t> vector_of_node_ids;
		for (ModelPart::ElementsContainerType::iterator it = rModelPart.ElementsBegin(); it != rModelPart.ElementsEnd(); ++it)
		{
			double elementDistance = 0.0;
			int numPointsInside = 0;
			std::size_t j = 0;
			for (j = 0; j < it->GetGeometry().PointsNumber(); j++)
			{
				elementDistance = it->GetGeometry()[j].FastGetSolutionStepValue(DISTANCE);
				if (elementDistance > lLimit && elementDistance < uLimit)
				{
					numPointsInside++;
				}
			}
			if (numPointsInside > 0)
			{
				pElem = Element::Pointer(new Element(*it));
				rExtractedModelPart.Elements().push_back(pElem);

				//Adding node all the node Ids of the elements satisfying the condition
				for (j = 0; j < pElem->GetGeometry().PointsNumber(); j++)
					vector_of_node_ids.push_back(pElem->GetGeometry()[j].Id());
			}
		}

		//sorting and making unique list of node ids
		std::set<std::size_t> s(vector_of_node_ids.begin(), vector_of_node_ids.end());
		vector_of_node_ids.assign(s.begin(), s.end());

		// Add unique nodes in the ModelPart
		for (auto it = vector_of_node_ids.begin(); it != vector_of_node_ids.end(); it++)
		{
			Node<3>::Pointer pnode = rModelPart.Nodes()(*it);
			rExtractedModelPart.AddNode(pnode);
		}
		KRATOS_CATCH("");
	}

	void CreateHoleAfterDistance(ModelPart &rModelPart, ModelPart &rExtractedModelPart, ModelPart &rExtractedBoundaryModelPart, double distance)
	{
		KRATOS_TRY;
		KRATOS_INFO("\n::[Creating Hole]::")<< std::endl;
		std::vector<std::size_t> vector_of_node_ids;
		//For signed distance
		distance *= -1;

		for (ModelPart::ElementsContainerType::iterator it = rModelPart.ElementsBegin(); it != rModelPart.ElementsEnd(); ++it)
		{
			double elementDistance = 0.0;
			std::size_t numPointsOutside = 0;
			std::size_t j = 0;
			Geometry<Node<3>> &geom = it->GetGeometry();

			for (j = 0; j < geom.size(); j++)
			{
				elementDistance = it->GetGeometry()[j].FastGetSolutionStepValue(DISTANCE);
				if (elementDistance < distance)
				{
					numPointsOutside++;
				}
			}

			if (numPointsOutside == geom.size())
			{
				it->Set(ACTIVE, false);
				Element::Pointer pElem = *(it.base());
				std::size_t numNodesPerElem = pElem->GetGeometry().PointsNumber();
				rExtractedModelPart.Elements().push_back(pElem);
				//Adding node all the node Ids of the elements satisfying the condition
				for (j = 0; j <numNodesPerElem; j++)
				{
					pElem->GetGeometry()[j].GetDof(VELOCITY_X).GetSolutionStepValue(0) = 0.0;
					pElem->GetGeometry()[j].GetDof(VELOCITY_Y).GetSolutionStepValue(0) = 0.0;
					if(numNodesPerElem-1 > 2)
						pElem->GetGeometry()[j].GetDof(VELOCITY_Z).GetSolutionStepValue(0) = 0.0;
					pElem->GetGeometry()[j].GetDof(PRESSURE).GetSolutionStepValue(0) = 0.0;
					pElem->GetGeometry()[j].GetDof(VELOCITY_X).GetSolutionStepValue(1) = 0.0;
					pElem->GetGeometry()[j].GetDof(VELOCITY_Y).GetSolutionStepValue(1) = 0.0;
					if(numNodesPerElem-1 > 2)
						pElem->GetGeometry()[j].GetDof(VELOCITY_Z).GetSolutionStepValue(1) = 0.0;
					pElem->GetGeometry()[j].GetDof(PRESSURE).GetSolutionStepValue(1) = 0.0;
					vector_of_node_ids.push_back(pElem->GetGeometry()[j].Id());
				}
			}
		}

		//sorting and making unique list of node ids
		std::set<std::size_t> s(vector_of_node_ids.begin(), vector_of_node_ids.end());
		vector_of_node_ids.assign(s.begin(), s.end());

		// Add unique nodes in the ModelPart
		for (auto it = vector_of_node_ids.begin(); it != vector_of_node_ids.end(); it++)
		{
			Node<3>::Pointer pnode = rModelPart.Nodes()(*it);
			rExtractedModelPart.AddNode(pnode);
		}

		std::size_t n_nodes = rModelPart.ElementsBegin()->GetGeometry().size();

		if (n_nodes == 3)
		{

			ExtractBoundaryMesh(rExtractedModelPart, rExtractedBoundaryModelPart);
		}

		else if (n_nodes == 4)
		{
			ExtractSurfaceMesh(rExtractedModelPart, rExtractedBoundaryModelPart);
		}

		else
			KRATOS_INFO("Hole cutting process is only supported for tetrahedral and triangular elements")<< std::endl;

		KRATOS_CATCH("");
	}

	void RemoveOutOfDomainPatchAndReturnModifiedPatch(ModelPart &rModelPart,ModelPart &rInsideBoundary, ModelPart &rExtractedModelPart, ModelPart &rExtractedBoundaryModelPart,int MainDomainOrNot)
	{
		KRATOS_TRY;
		KRATOS_INFO("\n:: Removing Out Of Domain Patch with Inside boundary Given ::") << std::endl;
		std::vector<std::size_t> vector_of_node_ids;

		int count = 0;

		for (ModelPart::ElementsContainerType::iterator it = rModelPart.ElementsBegin(); it != rModelPart.ElementsEnd(); ++it)
		{
			double elementDistance = 0.0;
			std::size_t numPointsOutside = 0;
			std::size_t j = 0;
			Geometry<Node<3>> &geom = it->GetGeometry();

			for (j = 0; j < geom.size(); j++)
			{
				elementDistance = it->GetGeometry()[j].FastGetSolutionStepValue(DISTANCE);

				elementDistance = elementDistance*MainDomainOrNot;
				if (elementDistance < 0)
				{
					numPointsOutside++;
				}
			}

			/* Any node goes out of the domain means the element need to be INACTIVE , otherwise the modified patch boundary
			 wont find any nodes on background */
			if (numPointsOutside > 0 )
			{
				it->Set(ACTIVE, false);
				Element::Pointer pElem = *(it.base());
				std::size_t numNodesPerElem = pElem->GetGeometry().PointsNumber();
				for (j = 0; j <numNodesPerElem; j++)
				{
					pElem->GetGeometry()[j].GetDof(VELOCITY_X).GetSolutionStepValue(0) = 0.0;
					pElem->GetGeometry()[j].GetDof(VELOCITY_Y).GetSolutionStepValue(0) = 0.0;
					if(numNodesPerElem-1 > 2)
						pElem->GetGeometry()[j].GetDof(VELOCITY_Z).GetSolutionStepValue(0) = 0.0;
					pElem->GetGeometry()[j].GetDof(PRESSURE).GetSolutionStepValue(0) = 0.0;
					pElem->GetGeometry()[j].GetDof(VELOCITY_X).GetSolutionStepValue(1) = 0.0;
					pElem->GetGeometry()[j].GetDof(VELOCITY_Y).GetSolutionStepValue(1) = 0.0;
					if(numNodesPerElem-1 > 2)
						pElem->GetGeometry()[j].GetDof(VELOCITY_Z).GetSolutionStepValue(1) = 0.0;
					pElem->GetGeometry()[j].GetDof(PRESSURE).GetSolutionStepValue(1) = 0.0;
				}
			}
			 else
			{
				count++;
				Element::Pointer pElem = *(it.base());
				std::size_t numNodesPerElem = pElem->GetGeometry().PointsNumber(); // Size()
				rExtractedModelPart.Elements().push_back(pElem); //AddElement()
				for (j = 0; j <numNodesPerElem; j++)
					vector_of_node_ids.push_back(pElem->GetGeometry()[j].Id());
			}
		}

		rExtractedModelPart.Nodes() = rModelPart.Nodes();

		KRATOS_INFO("Number of elements added to the modified patch")<<count<<std::endl;
		//sorting and making unique list of node ids
		std::set<std::size_t> s(vector_of_node_ids.begin(), vector_of_node_ids.end());
		vector_of_node_ids.assign(s.begin(), s.end());

		// Add unique nodes in the ModelPart
		for (auto it = vector_of_node_ids.begin(); it != vector_of_node_ids.end(); it++)
		{
			Node<3>::Pointer pnode = rModelPart.Nodes()(*it);
			rExtractedModelPart.AddNode(pnode);
		}

		std::size_t n_nodes = rModelPart.ElementsBegin()->GetGeometry().size();

		if (n_nodes == 3)
		{			
			KRATOS_INFO("::[Boundary line extraction of modified patch]::") << std::endl;
			ExtractOutsideBoundaryMesh(rInsideBoundary,rExtractedModelPart, rExtractedBoundaryModelPart);
		}

		else if (n_nodes == 4)
		{	
			KRATOS_INFO("::[3D surface extraction of the modified patch]::") << std::endl;
			ExtractOutsideSurfaceMesh(rInsideBoundary,rExtractedModelPart, rExtractedBoundaryModelPart);
		}

		else
			KRATOS_INFO("Hole cutting process is only supported for tetrahedral and triangular elements")<< std::endl;

		KRATOS_CATCH("");
	}

	/// Extracts a surface mesh from a provided volume mesh
	void ExtractSurfaceMesh(ModelPart &rExtractedVolumeModelPart, ModelPart &rExtractedSurfaceModelPart)
	{
		KRATOS_TRY;
		KRATOS_INFO ("::[Surface Mesh Extraction]::") << std::endl;

		// Some type-definitions
		typedef std::unordered_map<vector<std::size_t>, std::size_t, KeyHasher, KeyComparor> hashmap;
		typedef std::unordered_map<vector<std::size_t>, vector<std::size_t>, KeyHasher, KeyComparor> hashmap_vec;

		// Create map to ask for number of faces for the given set of node ids representing on face in the model part
		hashmap n_faces_map;

		// Fill map that counts number of faces for given set of nodes
		for (ModelPart::ElementIterator itElem = rExtractedVolumeModelPart.ElementsBegin(); itElem != rExtractedVolumeModelPart.ElementsEnd(); itElem++)
		{
			Element::GeometryType::GeometriesArrayType faces = itElem->GetGeometry().Faces();

			for (std::size_t face = 0; face < faces.size(); face++)
			{
				// Create vector that stores all node is of current face
				vector<std::size_t> ids(faces[face].size());

				// Store node ids
				for (std::size_t i = 0; i < faces[face].size(); i++)
					ids[i] = faces[face][i].Id();

				//*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
				std::sort(ids.begin(), ids.end());

				// Fill the map
				n_faces_map[ids] += 1;
			}
		}

		// Create a map to get nodes of skin face in original order for given set of node ids representing that face
		// The given set of node ids may have a different node order
		hashmap_vec ordered_skin_face_nodes_map;

		// Fill map that gives original node order for set of nodes
		for (ModelPart::ElementIterator itElem = rExtractedVolumeModelPart.ElementsBegin(); itElem != rExtractedVolumeModelPart.ElementsEnd(); itElem++)
		{
			Element::GeometryType::GeometriesArrayType faces = itElem->GetGeometry().Faces();

			for (std::size_t face = 0; face < faces.size(); face++)
			{
				// Create vector that stores all node is of current face
				vector<std::size_t> ids(faces[face].size());
				vector<std::size_t> unsorted_ids(faces[face].size());

				// Store node ids
				for (std::size_t i = 0; i < faces[face].size(); i++)
				{
					ids[i] = faces[face][i].Id();
					unsorted_ids[i] = faces[face][i].Id();
				}

				//*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
				std::sort(ids.begin(), ids.end());

				if (n_faces_map[ids] == 1)
					ordered_skin_face_nodes_map[ids] = unsorted_ids;
			}
		}
		// First assign to skin model part all nodes from original model_part, unnecessary nodes will be removed later
		std::size_t id_condition = 1;
		//rExtractedSurfaceModelPart.Nodes() = rExtractedVolumeModelPart.Nodes();

		// Add skin faces as triangles to skin-model-part (loop over all node sets)
		KRATOS_INFO("  Extracting surface mesh and computing normals") << std::endl;
		std::vector<std::size_t> vector_of_node_ids;
		for (typename hashmap::const_iterator it = n_faces_map.begin(); it != n_faces_map.end(); it++)
		{
			// If given node set represents face that is not overlapping with a face of another element, add it as skin element
			if (it->second == 1)
			{
				// If skin face is a triangle store triangle in with its original orientation in new skin model part
				if (it->first.size() == 3)
				{
					// Getting original order is important to properly reproduce skin face including its normal orientation
					vector<std::size_t> original_nodes_order = ordered_skin_face_nodes_map[it->first];
					Node<3>::Pointer pnode1 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[0]);
					Node<3>::Pointer pnode2 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[1]);
					Node<3>::Pointer pnode3 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[2]);

					//Storing the node ids list
					vector_of_node_ids.push_back(original_nodes_order[0]);
					vector_of_node_ids.push_back(original_nodes_order[1]);
					vector_of_node_ids.push_back(original_nodes_order[2]);
					Properties::Pointer properties = rExtractedSurfaceModelPart.rProperties()(0);
					Condition const &rReferenceTriangleCondition = KratosComponents<Condition>::Get("Condition3D");

					// Skin faces are added as conditions
					Triangle3D3<Node<3>> triangle1(pnode1, pnode2, pnode3);
					Condition::Pointer p_condition1 = rReferenceTriangleCondition.Create(id_condition++, triangle1, properties);
					rExtractedSurfaceModelPart.Conditions().push_back(p_condition1);
				}
				// If skin face is a quadrilateral then divide in two triangles and store them with their original orientation in new skin model part
				if (it->first.size() == 4)
				{
					// Getting original order is important to properly reproduce skin including its normal orientation
					vector<std::size_t> original_nodes_order = ordered_skin_face_nodes_map[it->first];

					Node<3>::Pointer pnode1 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[0]);
					Node<3>::Pointer pnode2 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[1]);
					Node<3>::Pointer pnode3 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[2]);
					Node<3>::Pointer pnode4 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[3]);
					//Storing the node ids list
					vector_of_node_ids.push_back(original_nodes_order[0]);
					vector_of_node_ids.push_back(original_nodes_order[1]);
					vector_of_node_ids.push_back(original_nodes_order[2]);
					vector_of_node_ids.push_back(original_nodes_order[3]);
					Properties::Pointer properties = rExtractedSurfaceModelPart.rProperties()(0);
					Condition const &rReferenceTriangleCondition = KratosComponents<Condition>::Get("Condition3D");

					// Add triangle one as condition
					Triangle3D3<Node<3>> triangle1(pnode1, pnode2, pnode3);
					Condition::Pointer p_condition1 = rReferenceTriangleCondition.Create(id_condition++, triangle1, properties);
					rExtractedSurfaceModelPart.Conditions().push_back(p_condition1);

					// Add triangle two as condition
					Triangle3D3<Node<3>> triangle2(pnode1, pnode3, pnode4);
					Condition::Pointer p_condition2 = rReferenceTriangleCondition.Create(id_condition++, triangle2, properties);
					rExtractedSurfaceModelPart.Conditions().push_back(p_condition2);
				}
			}
		}

		//sorting and making unique list of node ids

		std::set<std::size_t> s(vector_of_node_ids.begin(), vector_of_node_ids.end());
		vector_of_node_ids.assign(s.begin(), s.end());

		for (auto it = vector_of_node_ids.begin(); it != vector_of_node_ids.end(); it++)
		{
			//Adding the nodes to the rExtractedSurfaceModelPart
			Node<3>::Pointer pnode = rExtractedVolumeModelPart.Nodes()(*it);
			rExtractedSurfaceModelPart.AddNode(pnode);
		}

		KRATOS_INFO ("Successful extraction of the Surface ") << rExtractedSurfaceModelPart.GetMesh() << std::endl;

		KRATOS_CATCH("");
	}

	//give the outside surface of a model part given inside boundary
	void ExtractOutsideSurfaceMesh(ModelPart &rInsideBoundaryModelPart,ModelPart &rExtractedVolumeModelPart, ModelPart &rExtractedSurfaceModelPart)
	{
		KRATOS_TRY;

		KRATOS_INFO (":: [Surface Mesh Extraction]::") << std::endl;

			Parameters parameters= Parameters(R"({
						"result_file_configuration" : {
							"gidpost_flags"       : {
								"GiDPostMode"           : "GiD_PostAscii",
								"WriteDeformedMeshFlag" : "WriteDeformed",
								"WriteConditionsFlag"   : "WriteConditions",
								"MultiFileFlag"         : "SingleFile"
							},
							"file_label"          : "time",
							"output_control_type" : "time",
							"output_frequency"    : 1.0,
							"body_output"         : true,
							"node_output"         : false,
							"skin_output"         : false,
							"plane_output"        : [],
							"nodal_results"       : ["VELOCITY","PRESSURE","DISTANCE"],
							"gauss_point_results" : []
						},
						"point_data_configuration"  : []})" );

		VtkOutput VtkOutput_InsideBoundary = VtkOutput(rInsideBoundaryModelPart,"nnn",parameters);
		VtkOutput_InsideBoundary.PrintOutput();

		// Some type-definitions
		typedef std::unordered_map<vector<std::size_t>, std::size_t, KeyHasher, KeyComparor> hashmap;
		typedef std::unordered_map<vector<std::size_t>, vector<std::size_t>, KeyHasher, KeyComparor> hashmap_vec;

		// Create map to ask for number of faces for the given set of node ids representing on face in the model part
		hashmap n_faces_map;

		// Fill map that counts number of faces for given set of nodes
		for (ModelPart::ElementIterator itElem = rExtractedVolumeModelPart.ElementsBegin(); itElem != rExtractedVolumeModelPart.ElementsEnd(); itElem++)
		{
			Element::GeometryType::GeometriesArrayType faces = itElem->GetGeometry().Faces();

			for (std::size_t face = 0; face < faces.size(); face++)
			{
				// Create vector that stores all node is of current face
				vector<std::size_t> ids(faces[face].size());

				// Store node ids
				for (std::size_t i = 0; i < faces[face].size(); i++)
					ids[i] = faces[face][i].Id();

				//*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
				std::sort(ids.begin(), ids.end());

				// Fill the map
				n_faces_map[ids] += 1;
			}
		}

		KRATOS_INFO(":: comparing with inside boundary ")<<rInsideBoundaryModelPart<<std::endl;

		for (ModelPart::ConditionIterator itCond = rInsideBoundaryModelPart.ConditionsBegin(); itCond != rInsideBoundaryModelPart.ConditionsEnd(); ++itCond)
		{
			vector<std::size_t> ids(itCond->GetGeometry().size());
			for (std::size_t i = 0; i < itCond->GetGeometry().size();i++)
				ids[i] = itCond->GetGeometry()[i].Id();

			std::sort(ids.begin(), ids.end());
			if(n_faces_map[ids]==1)
				n_faces_map[ids] += 1;

		}

		// Create a map to get nodes of skin face in original order for given set of node ids representing that face
		// The given set of node ids may have a different node order
		hashmap_vec ordered_skin_face_nodes_map;

		// Fill map that gives original node order for set of nodes
		for (ModelPart::ElementIterator itElem = rExtractedVolumeModelPart.ElementsBegin(); itElem != rExtractedVolumeModelPart.ElementsEnd(); itElem++)
		{
			Element::GeometryType::GeometriesArrayType faces = itElem->GetGeometry().Faces();

			for (std::size_t face = 0; face < faces.size(); face++)
			{
				// Create vector that stores all node is of current face
				vector<std::size_t> ids(faces[face].size());
				vector<std::size_t> unsorted_ids(faces[face].size());

				// Store node ids
				for (std::size_t i = 0; i < faces[face].size(); i++)
				{
					ids[i] = faces[face][i].Id();
					unsorted_ids[i] = faces[face][i].Id();
				}

				//*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
				std::sort(ids.begin(), ids.end());
				if (n_faces_map[ids] == 1)
					ordered_skin_face_nodes_map[ids] = unsorted_ids;

			}
		}
		// First assign to skin model part all nodes from original model_part, unnecessary nodes will be removed later
		std::size_t id_condition = 1;
		//rExtractedSurfaceModelPart.Nodes() = rExtractedVolumeModelPart.Nodes();

		// Add skin faces as triangles to skin-model-part (loop over all node sets)
		std::vector<std::size_t> vector_of_node_ids;
		for (typename hashmap::const_iterator it = n_faces_map.begin(); it != n_faces_map.end(); it++)
		{
			// If given node set represents face that is not overlapping with a face of another element, add it as skin element
			if (it->second == 1)
			{
				// If skin face is a triangle store triangle in with its original orientation in new skin model part
				if (it->first.size() == 3)
				{
					// Getting original order is important to properly reproduce skin face including its normal orientation
					vector<std::size_t> original_nodes_order = ordered_skin_face_nodes_map[it->first];
					Node<3>::Pointer pnode1 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[0]);
					Node<3>::Pointer pnode2 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[1]);
					Node<3>::Pointer pnode3 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[2]);

					//Storing the node ids list
					vector_of_node_ids.push_back(original_nodes_order[0]);
					vector_of_node_ids.push_back(original_nodes_order[1]);
					vector_of_node_ids.push_back(original_nodes_order[2]);
					Properties::Pointer properties = rExtractedSurfaceModelPart.rProperties()(0);
					Condition const &rReferenceTriangleCondition = KratosComponents<Condition>::Get("Condition3D");

					// Skin faces are added as conditions
					Triangle3D3<Node<3>> triangle1(pnode1, pnode2, pnode3);
					Condition::Pointer p_condition1 = rReferenceTriangleCondition.Create(id_condition++, triangle1, properties);
					rExtractedSurfaceModelPart.Conditions().push_back(p_condition1);
				}
				// If skin face is a quadrilateral then divide in two triangles and store them with their original orientation in new skin model part
				if (it->first.size() == 4)
				{
					// Getting original order is important to properly reproduce skin including its normal orientation
					vector<std::size_t> original_nodes_order = ordered_skin_face_nodes_map[it->first];

					Node<3>::Pointer pnode1 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[0]);
					Node<3>::Pointer pnode2 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[1]);
					Node<3>::Pointer pnode3 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[2]);
					Node<3>::Pointer pnode4 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[3]);
					//Storing the node ids list
					vector_of_node_ids.push_back(original_nodes_order[0]);
					vector_of_node_ids.push_back(original_nodes_order[1]);
					vector_of_node_ids.push_back(original_nodes_order[2]);
					vector_of_node_ids.push_back(original_nodes_order[3]);
					Properties::Pointer properties = rExtractedSurfaceModelPart.rProperties()(0);
					Condition const &rReferenceTriangleCondition = KratosComponents<Condition>::Get("Condition3D");

					// Add triangle one as condition
					Triangle3D3<Node<3>> triangle1(pnode1, pnode2, pnode3);
					Condition::Pointer p_condition1 = rReferenceTriangleCondition.Create(id_condition++, triangle1, properties);
					rExtractedSurfaceModelPart.Conditions().push_back(p_condition1);

					// Add triangle two as condition
					Triangle3D3<Node<3>> triangle2(pnode1, pnode3, pnode4);
					Condition::Pointer p_condition2 = rReferenceTriangleCondition.Create(id_condition++, triangle2, properties);
					rExtractedSurfaceModelPart.Conditions().push_back(p_condition2);
				}
			}
		}
		//sorting and making unique list of node ids

		std::set<std::size_t> s(vector_of_node_ids.begin(), vector_of_node_ids.end());
		vector_of_node_ids.assign(s.begin(), s.end());

		for (auto it = vector_of_node_ids.begin(); it != vector_of_node_ids.end(); it++)

		{
			//Adding the nodes to the rExtractedSurfaceModelPart
			Node<3>::Pointer pnode = rExtractedVolumeModelPart.Nodes()(*it);
			rExtractedSurfaceModelPart.AddNode(pnode);
		}

		KRATOS_INFO ("Successful extraction of the Surface" )<< rExtractedSurfaceModelPart.GetMesh() << std::endl;

		VtkOutput VtkOutput_Extracted_Surface = VtkOutput(rExtractedSurfaceModelPart,"ExtracetedModelPart",parameters);
		VtkOutput_Extracted_Surface.PrintOutput();

		KRATOS_CATCH("");
	}

	void ExtractBoundaryMesh(ModelPart &rSurfaceModelPart, ModelPart &rExtractedBoundaryModelPart)
	{
		KRATOS_TRY;

		KRATOS_INFO ("::[Boundary Mesh Extraction]::")<< std::endl;

		// Some type-definitions
		typedef std::unordered_map<vector<std::size_t>, std::size_t, KeyHasher, KeyComparor> hashmap;
		typedef std::unordered_map<vector<std::size_t>, vector<std::size_t>, KeyHasher, KeyComparor> hashmap_vec;

		// Create map to ask for number of edges for the given set of node ids representing on edge in the model part
		hashmap n_edges_map;

		// Fill map that counts number of edges for given set of nodes
		for (ModelPart::ElementIterator itElem = rSurfaceModelPart.ElementsBegin(); itElem != rSurfaceModelPart.ElementsEnd(); itElem++)
		{
			Element::GeometryType::GeometriesArrayType edges = itElem->GetGeometry().Edges();

			for (std::size_t edge = 0; edge < edges.size(); edge++)
			{
				// Create vector that stores all node ids of current edge
				vector<std::size_t> ids(edges[edge].size());

				// Store node ids
				for (std::size_t i = 0; i < edges[edge].size(); i++)
					ids[i] = edges[edge][i].Id();

				//*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
				std::sort(ids.begin(), ids.end());

				// Fill the map
				n_edges_map[ids] += 1;
			}
		}

		// Create a map to get nodes of skin edge in original order for given set of node ids representing that edge
		// The given set of node ids may have a different node order
		hashmap_vec ordered_skin_edge_nodes_map;

		// Fill map that gives original node order for set of nodes
		for (ModelPart::ElementIterator itElem = rSurfaceModelPart.ElementsBegin(); itElem != rSurfaceModelPart.ElementsEnd(); itElem++)
		{
			Element::GeometryType::GeometriesArrayType edges = itElem->GetGeometry().Edges();

			for (std::size_t edge = 0; edge < edges.size(); edge++)
			{
				// Create vector that stores all node is of current edge
				vector<std::size_t> ids(edges[edge].size());
				vector<std::size_t> unsorted_ids(edges[edge].size());

				// Store node ids
				for (std::size_t i = 0; i < edges[edge].size(); i++)
				{
					ids[i] = edges[edge][i].Id();
					unsorted_ids[i] = edges[edge][i].Id();
				}

				//*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
				std::sort(ids.begin(), ids.end());

				if (n_edges_map[ids] == 1)
					ordered_skin_edge_nodes_map[ids] = unsorted_ids;
			}
		}
		// First assign to skin model part all nodes from original model_part, unnecessary nodes will be removed later
		std::size_t id_condition = 1;
		//rExtractedBoundaryModelPart.Nodes() = rSurfaceModelPart.Nodes();

		// Add skin edges as triangles to skin-model-part (loop over all node sets)
		KRATOS_INFO("  Extracting boundary mesh and computing normals") << std::endl;
		std::vector<std::size_t> vector_of_node_ids;
		for (typename hashmap::const_iterator it = n_edges_map.begin(); it != n_edges_map.end(); it++)
		{
			// If given node set represents edge that is not overlapping with a edge of another element, add it as skin element
			if (it->second == 1)
			{
				// If skin edge is a triangle store triangle in with its original orientation in new skin model part
				//KRATOS_INFO("size of the ordered pair : ")<<it->first.size()<<std::endl;
				if (it->first.size() == 2)
				{
					// Getting original order is important to properly reproduce skin edge including its normal orientation
					vector<std::size_t> original_nodes_order = ordered_skin_edge_nodes_map[it->first];

					//KRATOS_INFO("First Node: ")<<original_nodes_order[0]<<std::endl;
					//KRATOS_INFO("Second Node: ")<<original_nodes_order[1]<<std::endl;

					Node<3>::Pointer pnode1 = rSurfaceModelPart.Nodes()(original_nodes_order[0]);
					Node<3>::Pointer pnode2 = rSurfaceModelPart.Nodes()(original_nodes_order[1]);

					//Storing the node ids list
					vector_of_node_ids.push_back(original_nodes_order[0]);
					vector_of_node_ids.push_back(original_nodes_order[1]);

					Properties::Pointer properties = rExtractedBoundaryModelPart.rProperties()(0);
					Condition const &rReferenceLineCondition = KratosComponents<Condition>::Get("Condition2D");

					// Skin edges are added as conditions
					Line2D2<Node<3>> line1(pnode1, pnode2);
					Condition::Pointer p_condition1 = rReferenceLineCondition.Create(id_condition++, line1, properties);
					rExtractedBoundaryModelPart.Conditions().push_back(p_condition1);
				}
			}
		}

		//sorting and making unique list of node ids

		std::set<std::size_t> s(vector_of_node_ids.begin(), vector_of_node_ids.end());
		vector_of_node_ids.assign(s.begin(), s.end());

		for (auto it = vector_of_node_ids.begin(); it != vector_of_node_ids.end(); it++)

		{
			//Adding the nodes to the rExtractedSurfaceModelPart
			Node<3>::Pointer pnode = rSurfaceModelPart.Nodes()(*it);
			rExtractedBoundaryModelPart.AddNode(pnode);
		}

		KRATOS_INFO("Successful extraction of the Boundary") << rExtractedBoundaryModelPart.GetMesh() << std::endl;

		KRATOS_CATCH("");
	}


	void ExtractOutsideBoundaryMesh(ModelPart &rInsideBoundaryModelPart,ModelPart &rSurfaceModelPart, ModelPart &rExtractedBoundaryModelPart)
	{
		KRATOS_TRY;



		// Some type-definitions
		typedef std::unordered_map<vector<std::size_t>, std::size_t, KeyHasher, KeyComparor> hashmap;
		typedef std::unordered_map<vector<std::size_t>, vector<std::size_t>, KeyHasher, KeyComparor> hashmap_vec;

		// Create map to ask for number of edges for the given set of node ids representing on edge in the model part
		hashmap n_edges_map;

		// Fill map that counts number of edges for given set of nodes
		for (ModelPart::ElementIterator itElem = rSurfaceModelPart.ElementsBegin(); itElem != rSurfaceModelPart.ElementsEnd(); itElem++)
		{
			Element::GeometryType::GeometriesArrayType edges = itElem->GetGeometry().Edges();

			for (std::size_t edge = 0; edge < edges.size(); edge++)
			{
				// Create vector that stores all node ids of current edge
				vector<std::size_t> ids(edges[edge].size());

				// Store node ids
				for (std::size_t i = 0; i < edges[edge].size(); i++)
					ids[i] = edges[edge][i].Id();

				//*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
				std::sort(ids.begin(), ids.end());

				// Fill the map
				n_edges_map[ids] += 1;
			}
		}

		for (ModelPart::ConditionIterator itCond = rInsideBoundaryModelPart.ConditionsBegin(); itCond != rInsideBoundaryModelPart.ConditionsEnd(); ++itCond)
		{
			vector<std::size_t> ids(itCond->GetGeometry().size());
			for (std::size_t i = 0; i < itCond->GetGeometry().size();i++)
				ids[i] = itCond->GetGeometry()[i].Id();

			std::sort(ids.begin(), ids.end());
			n_edges_map[ids] += 1;
		}

		// Create a map to get nodes of skin edge in original order for given set of node ids representing that edge
		// The given set of node ids may have a different node order
		hashmap_vec ordered_skin_edge_nodes_map;

		// Fill map that gives original node order for set of nodes
		for (ModelPart::ElementIterator itElem = rSurfaceModelPart.ElementsBegin(); itElem != rSurfaceModelPart.ElementsEnd(); itElem++)
		{
			Element::GeometryType::GeometriesArrayType edges = itElem->GetGeometry().Edges();

			for (std::size_t edge = 0; edge < edges.size(); edge++)
			{
				// Create vector that stores all node is of current edge
				vector<std::size_t> ids(edges[edge].size());
				vector<std::size_t> unsorted_ids(edges[edge].size());

				// Store node ids
				for (std::size_t i = 0; i < edges[edge].size(); i++)
				{
					ids[i] = edges[edge][i].Id();
					unsorted_ids[i] = edges[edge][i].Id();
				}

				//*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
				std::sort(ids.begin(), ids.end());

				if (n_edges_map[ids] == 1)
					ordered_skin_edge_nodes_map[ids] = unsorted_ids;
			}
		}
		// First assign to skin model part all nodes from original model_part, unnecessary nodes will be removed later
		std::size_t id_condition = 1;
		//rExtractedBoundaryModelPart.Nodes() = rSurfaceModelPart.Nodes();

		// Add skin edges as triangles to skin-model-part (loop over all node sets)
		KRATOS_INFO("  Extracting boundary mesh and computing normals")<< std::endl;
		std::vector<std::size_t> vector_of_node_ids;
		for (typename hashmap::const_iterator it = n_edges_map.begin(); it != n_edges_map.end(); it++)
		{
			// If given node set represents edge that is not overlapping with a edge of another element, add it as skin element
			if (it->second == 1)
			{
				// If skin edge is a triangle store triangle in with its original orientation in new skin model part
				//KRATOS_INFO("size of the ordered pair : ")<<it->first.size()<<std::endl;
				if (it->first.size() == 2)
				{
					// Getting original order is important to properly reproduce skin edge including its normal orientation
					vector<std::size_t> original_nodes_order = ordered_skin_edge_nodes_map[it->first];

					//KRATOS_INFO("First Node: ")<<original_nodes_order[0]<<std::endl;
					//KRATOS_INFO("Second Node: ")<<original_nodes_order[1]<<std::endl;
					Node<3>::Pointer pnode1 = rSurfaceModelPart.Nodes()(original_nodes_order[0]);
					Node<3>::Pointer pnode2 = rSurfaceModelPart.Nodes()(original_nodes_order[1]);

					//Storing the node ids list
					vector_of_node_ids.push_back(original_nodes_order[0]);
					vector_of_node_ids.push_back(original_nodes_order[1]);

					Properties::Pointer properties = rExtractedBoundaryModelPart.rProperties()(0);
					Condition const &rReferenceLineCondition = KratosComponents<Condition>::Get("Condition2D");

					// Skin edges are added as conditions
					Line2D2<Node<3>> line1(pnode1, pnode2);
					Condition::Pointer p_condition1 = rReferenceLineCondition.Create(id_condition++, line1, properties);
					rExtractedBoundaryModelPart.Conditions().push_back(p_condition1);
				}
			}
		}

		//sorting and making unique list of node ids

		std::set<std::size_t> s(vector_of_node_ids.begin(), vector_of_node_ids.end());
		vector_of_node_ids.assign(s.begin(), s.end());

		for (auto it = vector_of_node_ids.begin(); it != vector_of_node_ids.end(); it++)

		{
			//Adding the nodes to the rExtractedSurfaceModelPart
			Node<3>::Pointer pnode = rSurfaceModelPart.Nodes()(*it);
			rExtractedBoundaryModelPart.AddNode(pnode);
		}

		KRATOS_INFO("Successful extraction of the Boundary ")<< rExtractedBoundaryModelPart.GetMesh() << std::endl;

		KRATOS_CATCH("");
	}

	virtual void Execute()
	{
	}

	virtual void Clear()
	{
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
	virtual std::string Info() const
	{
		return "CustomHoleCuttingProcess";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "CustomHoleCuttingProcess";
	}

	/// Print object's data.
	virtual void PrintData(std::ostream &rOStream) const
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

	///@}
	///@name Protected Operators
	///@{

	///@}
	///@name Protected Operations
	///@{

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
	CustomHoleCuttingProcess &operator=(CustomHoleCuttingProcess const &rOther);

	/// Copy constructor.
	//CustomHoleCuttingProcess(CustomHoleCuttingProcess const& rOther);

	///@}

}; // Class CustomHoleCuttingProcess

} // namespace Kratos.

#endif // KRATOS_CUSTOM_HOLE_CUTTING_PROCESS_H_INCLUDED  defined
