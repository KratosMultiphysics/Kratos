// ==============================================================================
//  KratosTopologyOptimizationApplication
//
//  License:         BSD License
//                   license: TopologyOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//
// ==============================================================================

#if !defined(KRATOS_TOPOLOGY_EXTRACTOR_UTILITIES_H_INCLUDED)
#define  KRATOS_TOPOLOGY_EXTRACTOR_UTILITIES_H_INCLUDED

// System includes
#include <iostream>
#include <string>
#include <algorithm>

// External includes
#include <pybind11/pybind11.h>   
#include <boost/numeric/ublas/storage.hpp>

#include "../custom_elements/small_displacement_simp_element.h"
// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/process_info.h"
#include "includes/key_hash.h"

// Application includes
#include "topology_optimization_application.h"
#include "processes/find_global_nodal_neighbours_process.h" // To find node neighbours
#include "processes/find_conditions_neighbours_process.h" // To find condition neighbours
#include "processes/node_erase_process.h" // To delete empty nodes
#include "utilities/normal_calculation_utils.h" // To calculate element's normal
#include "geometries/triangle_3d_3.h" // Skin face geometry template
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

/// Solution utility to extract optimized meshes.
/** Detail class definition.

 */

class TopologyExtractorUtilities
{
public:

	///@name Type Definitions
	///@{

	/// Pointer definition of TopologyExtractorUtilities
	KRATOS_CLASS_POINTER_DEFINITION(TopologyExtractorUtilities);

	// Needed structures for the ExtractSurfaceMesh operation
	struct KeyComparor
	{
		bool operator()(const vector<unsigned int>& lhs, const vector<unsigned int>& rhs) const
		{
			if(lhs.size() != rhs.size())
				return false;

			for(unsigned int i=0; i<lhs.size(); i++)
			{
				if(lhs[i] != rhs[i]) return false;
			}

			return true;
		}
	};

	struct KeyHasher
	{
		std::size_t operator()(const vector<int>& k) const
		{
			return HashRange(k.begin(), k.end());
		}
	};

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	TopologyExtractorUtilities()
	{
	}

	/// Destructor.
	virtual ~TopologyExtractorUtilities()
	{
	}


	///@}
	///@name Operators
	///@{


	///@}
	///@name Operations
	///@{

	// ---------------------------------------------------------------------------------------------------------------------------------------------
	// --------------------------------- EXTRACT OPTIMIZED VOLUME MESH  ----------------------------------------------------------------------------
	// ---------------------------------------------------------------------------------------------------------------------------------------------

	/// For Topology Optimization purposes: Extracts a volume mesh with a certain threshold value
	void ExtractVolumeMesh(ModelPart& rModelPart, double threshold, ModelPart& rExtractedModelPart)
	{

		KRATOS_TRY;

		std::cout<<"\n::[Volume Mesh Extraction]::"<<std::endl;

		// Initializing mesh nodes
		rExtractedModelPart.Nodes() = rModelPart.Nodes();
		const DataCommunicator& r_comm = rModelPart.GetCommunicator().GetDataCommunicator();

		// Extracting mesh elements which are only above the threshold value
		std::cout<<"  Extracting elements with X_PHYS > " << threshold <<std::endl;
		Element::Pointer pElem;
		for(ModelPart::ElementsContainerType::iterator it = rModelPart.ElementsBegin(); it != rModelPart.ElementsEnd(); ++it)
		{
			if(it->GetValue(X_PHYS) > threshold)
			{
				pElem = Element::Pointer(new SmallDisplacementSIMPElement(
						(*it).Id(),
						(*it).pGetGeometry(),
						(*it).pGetProperties() ) );
				rExtractedModelPart.Elements().push_back(pElem);
			}
		}

		// Erase free nodes
		std::cout<<"  Erasing free nodes and re-numbering kept nodes" <<std::endl;
		FindGlobalNodalNeighboursProcess nodal_finder = FindGlobalNodalNeighboursProcess(r_comm,rExtractedModelPart,10);
		nodal_finder.Execute();

		for(NodesContainerType::iterator node_i=rExtractedModelPart.NodesBegin(); node_i!=rExtractedModelPart.NodesEnd(); node_i++)
		{
			GlobalPointersVector<Element >& ng_elem = node_i->GetValue(NEIGHBOUR_ELEMENTS);
			if(ng_elem.size()==0)
				node_i->Set(TO_ERASE,true);
		}
		(NodeEraseProcess(rExtractedModelPart)).Execute();

		// Renumber nodes in extracted volume mesh (to start from 1)
		unsigned int new_id = 1;
		for(ModelPart::NodesContainerType::iterator node_i = rExtractedModelPart.NodesBegin();  node_i!=rExtractedModelPart.NodesEnd(); node_i++)
			node_i->SetId(new_id++);

		std::cout<<"  Successful extraction of the Volume "<<rExtractedModelPart.GetMesh()<<"\b"<<std::endl;

		KRATOS_CATCH("");

	}

	// ---------------------------------------------------------------------------------------------------------------------------------------------
	// --------------------------------- EXTRACT OPTIMIZED SURFACE MESH  ---------------------------------------------------------------------------
	// ---------------------------------------------------------------------------------------------------------------------------------------------

	/// For Topology Optimization purposes: Extracts a surface mesh from a provided volume mesh
	void ExtractSurfaceMesh(ModelPart& rExtractedVolumeModelPart, ModelPart& rExtractedSurfaceModelPart)
	{
		KRATOS_TRY;

		std::cout<<"::[Surface Mesh Extraction]::"<<std::endl;

		// Some type-definitions
		typedef std::unordered_map<vector<unsigned int>, unsigned int, KeyHasher, KeyComparor > hashmap;
		typedef std::unordered_map<vector<unsigned int>, vector<unsigned int>, KeyHasher, KeyComparor > hashmap_vec;

		// Some working variable
		unsigned int domain_size = 3;

		// Create map to ask for number of faces for the given set of node ids representing on face in the model part
		hashmap n_faces_map;

		// Fill map that counts number of faces for given set of nodes
		for (ModelPart::ElementIterator itElem = rExtractedVolumeModelPart.ElementsBegin(); itElem != rExtractedVolumeModelPart.ElementsEnd(); itElem++)
		{
			Element::GeometryType::GeometriesArrayType faces = itElem->GetGeometry().Faces();

			for(unsigned int face=0; face<faces.size(); face++)
			{
				// Create vector that stores all node is of current face
				vector<unsigned int> ids(faces[face].size());

				// Store node ids
				for(unsigned int i=0; i<faces[face].size(); i++)
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

			for(unsigned int face=0; face<faces.size(); face++)
			{
				// Create vector that stores all node is of current face
				vector<unsigned int> ids(faces[face].size());
				vector<unsigned int> unsorted_ids(faces[face].size());

				// Store node ids
				for(unsigned int i=0; i<faces[face].size(); i++)
				{
					ids[i] = faces[face][i].Id();
					unsorted_ids[i] = faces[face][i].Id();
				}

				//*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
				std::sort(ids.begin(), ids.end());

				if(n_faces_map[ids] == 1)
					ordered_skin_face_nodes_map[ids] = unsorted_ids;
			}
		}

		// First assign to skin model part all nodes from original model_part, unnecessary nodes will be removed later
		unsigned int face_id = 1;
		rExtractedSurfaceModelPart.Nodes() = rExtractedVolumeModelPart.Nodes();

		// Create reference conditions and elements to be assigned to a new skin-model-part
		Condition const& rReferenceTriangleCondition = KratosComponents<Condition>::Get("SurfaceCondition3D3N");
		Element const& rReferenceTriangleElement = KratosComponents<Element>::Get("SmallDisplacementSIMPElement3D3N");

		// Add skin faces as triangles to skin-model-part (loop over all node sets)
		// Add skin face both as condition and dummy element, since different useful Kratos utilities work with both elements and conditions
		std::cout<<"  Extracting surface mesh and computing normals" <<std::endl;
		for(typename hashmap::const_iterator it=n_faces_map.begin(); it!=n_faces_map.end(); it++)
		{	
			std::cout<<"  Schleife 1" <<std::endl;
			// If given node set represents face that is not overlapping with a face of another element, add it as skin element
			if(it->second == 1)
			{	std::cout<<"  Schleife 2" <<std::endl;
				// If skin face is a triangle store triangle in with its original orientation in new skin model part
				if(it->first.size()==3)
				{	std::cout<<"  Schleife 3" <<std::endl;
					// Getting original order is important to properly reproduce skin face including its normal orientation
					std::cout<<"  before constructing vector" <<std::endl;
					vector<unsigned int> original_nodes_order = ordered_skin_face_nodes_map[it->first];
					std::cout<<"  after constructing vector" <<std::endl;

					Node < 3 >::Pointer pnode1 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[0]);
					std::cout<<"  pnode1 done" <<std::endl;
					Node < 3 >::Pointer pnode2 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[1]);
					std::cout<<"  pnode2 done" <<std::endl;
					Node < 3 >::Pointer pnode3 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[2]);
					std::cout<<"  pnode3 done" <<std::endl;
					Properties::Pointer properties = rExtractedSurfaceModelPart.rProperties()(0);
					std::cout<<"  Surfacemesh done" <<std::endl;

					// Add skin face as condition
					Triangle3D3< Node<3> > triangle_c(pnode1, pnode2, pnode3);
					Condition::Pointer p_condition1 = rReferenceTriangleCondition.Create(face_id++, triangle_c, properties);
					rExtractedSurfaceModelPart.Conditions().push_back(p_condition1);

					// Add skin face as element
					Triangle3D3< Node<3> > triangle_e(pnode1, pnode2, pnode3);
					Element::Pointer p_element_1 = rReferenceTriangleElement.Create(face_id++, triangle_e, properties);
					rExtractedSurfaceModelPart.Elements().push_back(p_element_1);
				}

				// If skin face is a quadrilateral then divide in two triangles and store them with their original orientation in new skin model part
				if(it->first.size()==4)
				{	
					std::cout<<"  Schleife 4" <<std::endl;
					// Getting original order is important to properly reproduce skin including its normal orientation
					std::cout<<"  Step 1" <<std::endl;
					vector<unsigned int> original_nodes_order = ordered_skin_face_nodes_map[it->first];
					std::cout<<"  Step 2" <<std::endl;

					Node < 3 >::Pointer pnode1 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[0]);
					std::cout<<"  pnode1 done" <<std::endl;
					Node < 3 >::Pointer pnode2 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[1]);
					std::cout<<"  pnode2 done" <<std::endl;
					Node < 3 >::Pointer pnode3 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[2]);
					std::cout<<"  pnode3 done" <<std::endl;
					Node < 3 >::Pointer pnode4 = rExtractedVolumeModelPart.Nodes()(original_nodes_order[3]);
					std::cout<<"  pnode4 done" <<std::endl;
					Properties::Pointer properties = rExtractedSurfaceModelPart.rProperties()(0);
					std::cout<<"  Surfacemesh done" <<std::endl;

					// Add triangle one as condition
					std::cout<<"  Step 3" <<std::endl;
					Triangle3D3< Node<3> > triangle1_c(pnode1, pnode2, pnode3);
					std::cout<<"  Step 4" <<std::endl;
					Condition::Pointer p_condition1 = rReferenceTriangleCondition.Create(face_id++, triangle1_c, properties);
					std::cout<<"  Step 5" <<std::endl;
					rExtractedSurfaceModelPart.Conditions().push_back(p_condition1);

					// Add triangle two as condition
					Triangle3D3< Node<3> > triangle2_c(pnode1, pnode3, pnode4);
					Condition::Pointer p_condition2 = rReferenceTriangleCondition.Create(face_id++, triangle2_c, properties);
					rExtractedSurfaceModelPart.Conditions().push_back(p_condition2);

					// Add triangle one as element
					Triangle3D3< Node<3> > triangle1_e(pnode1, pnode2, pnode3);
					Element::Pointer p_elem1 = rReferenceTriangleElement.Create(face_id++, triangle1_e, properties);
					rExtractedSurfaceModelPart.Elements().push_back(p_elem1);

					// Add triangle two as element
					Triangle3D3< Node<3> > triangle2_e(pnode1, pnode3, pnode4);
					Element::Pointer p_elem2 = rReferenceTriangleElement.Create(face_id++, triangle2_e, properties);
					rExtractedSurfaceModelPart.Elements().push_back(p_elem2);
				}
			}
		}
		std::cout<<"  Bis hierher gehts" <<std::endl;
		// Remove free nodes (nodes which do not have any link to other nodes through a triangle, hence to not belong to the skin)
		(FindConditionsNeighboursProcess(rExtractedSurfaceModelPart,domain_size, 10)).Execute();
		for(ModelPart::NodesContainerType::iterator node_i =  rExtractedSurfaceModelPart.NodesBegin(); node_i !=rExtractedSurfaceModelPart.NodesEnd(); node_i++)
		{
			GlobalPointersVector<Condition >& ng_cond = node_i->GetValue(NEIGHBOUR_CONDITIONS);
			if(ng_cond.size()==0)
				node_i->Set(TO_ERASE,true);
			std::cout<<"  Schleife 5" <<std::endl;
		}
		(NodeEraseProcess(rExtractedSurfaceModelPart)).Execute();

		// Compute normalized surface normals of skin model part
		NormalCalculationUtils normal_util = NormalCalculationUtils();
		normal_util.CalculateOnSimplex(rExtractedSurfaceModelPart,domain_size);
		for ( ModelPart::ConditionIterator cond_i = rExtractedSurfaceModelPart.ConditionsBegin(); cond_i != rExtractedSurfaceModelPart.ConditionsEnd(); ++cond_i )
		{
			// Normalize normal and assign to solution step value
			noalias(cond_i->GetValue(NORMAL)) = cond_i->GetValue(NORMAL) / norm_2(cond_i->GetValue(NORMAL));
		}

		// Renumber nodes in skin model part (to start from 1)
		std::cout<<"  Erasing free nodes and re-numbering kept nodes" <<std::endl;
		unsigned int new_id = 1;
		for(ModelPart::NodesContainerType::iterator node_i = rExtractedSurfaceModelPart.NodesBegin();  node_i!=rExtractedSurfaceModelPart.NodesEnd(); node_i++)
			node_i->SetId(new_id++);

		std::cout<<"  Successful extraction of the Surface "<<rExtractedSurfaceModelPart.GetMesh()<<std::endl;

		KRATOS_CATCH("");
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
		return "TopologyExtractorUtilities";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream& rOStream) const
	{
		rOStream << "TopologyExtractorUtilities";
	}

	/// Print object's data.
	virtual void PrintData(std::ostream& rOStream) const
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
	//TopologyExtractorUtilities& operator=(TopologyExtractorUtilities const& rOther);

	/// Copy constructor.
	//TopologyExtractorUtilities(TopologyExtractorUtilities const& rOther);


	///@}

}; // Class TopologyExtractorUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif	/* KRATOS_TOPOLOGY_EXTRACTOR_UTILITIES_H_INCLUDED */
