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

// External includes

// Project includes

// System includes
#include <iostream>
#include <string>
#include <algorithm>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/process_info.h"
#include "containers/model.h"
// Application includes
#include "processes/node_erase_process.h"                 // To delete empty nodes
#include "geometries/triangle_3d_3.h"                     // Skin face geometry template
#include "geometries/line_2d_2.h"

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

class KRATOS_API(CHIMERA_APPLICATION) ChimeraHoleCuttingUtility
{
public:
    // Needed structures for the ExtractSurfaceMesh operation
    struct KeyComparator
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
                seed ^= hasher(k[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }

            return seed;
        }
    };

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of ChimeraHoleCuttingUtility
    KRATOS_CLASS_POINTER_DEFINITION(ChimeraHoleCuttingUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    ChimeraHoleCuttingUtility()
    {
    }

    /// Destructor.
    virtual ~ChimeraHoleCuttingUtility()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a hole on the given rModelPart at a Distance (-ve) from the zero distance layer.
     * @param rModelPart The modelpart where the hole is to be cut.
     * @param rHoleModelPart The modelpart containing the nodes and corresponding elements from the cut hole. (deactivated elements)
     * @param rHoleBoundaryModelPart Boundary modelpart of the rHoleModelPart
     * @param Distance is the the distance (magnitude) at which hole is to be cut from the zero distance layer.
     */
    void CreateHoleAfterDistance(ModelPart &rModelPart, ModelPart &rHoleModelPart, ModelPart &rHoleBoundaryModelPart, const double Distance)
    {
        KRATOS_TRY;
        KRATOS_INFO("\n::[Creating Hole]::") << std::endl;
        RemoveOutOfDomainElements(rModelPart, rHoleModelPart,1,Distance, true);
        ExtractBoundaryMesh(rHoleModelPart, rHoleBoundaryModelPart);
        KRATOS_CATCH("");
    }


    /**
     * @brief Removes the elements which are out of the domain.
     *          An element is removed even if one of its nodes is out of the domain (-ve or +ve) as indicated by GetInside and MainDomainOrNot
     * @param rModelPart The modelpart From where the elements are to be removed.
     * @param rModifiedModelPart The modified modelpart without the the elements which are out side.
     * @param MainDomainOrNot says which sign (-ve or +ve) is inside
     * @param OverLapDistance is the the distance (magnitude) at which hole is to be cut from the zero distance layer.
     * @param GetInside works in combination with MainDomainOrNot to get the feeling of what is inside or what is outside.
     */
    void RemoveOutOfDomainElements(ModelPart &rModelPart,
                                   ModelPart &rModifiedModelPart,
                                   const int MainDomainOrNot,
                                   const double OverLapDistance=0.0,
                                   const bool GetInside=false)
    {
        KRATOS_TRY;
        KRATOS_INFO("\n:: Removing Out Of Domain Patch with Inside boundary Given ::") << std::endl;
        std::vector<std::size_t> vector_of_node_ids;

        int count = 0;

        for (auto& i_element :rModelPart.Elements())
        {
            double nodal_distance = 0.0;
            std::size_t numPointsOutside = 0;
            std::size_t j = 0;
            Geometry<Node<3>> &geom = i_element.GetGeometry();

            for (j = 0; j < geom.size(); j++)
            {
                nodal_distance = i_element.GetGeometry()[j].FastGetSolutionStepValue(DISTANCE);

                nodal_distance = nodal_distance * MainDomainOrNot;
                if (nodal_distance < -1*OverLapDistance)
                {
                    numPointsOutside++;
                }
            }

            /* Any node goes out of the domain means the element need to be INACTIVE , otherwise the modified patch boundary
			 wont find any nodes on background */
            if (numPointsOutside > 0)
            {
                i_element.Set(ACTIVE, false);
                std::size_t num_nodes_per_elem = i_element.GetGeometry().PointsNumber();
                if(GetInside)
                    rModifiedModelPart.AddElement(rModelPart.pGetElement(i_element.Id()));
                for (j = 0; j < num_nodes_per_elem; j++)
                {
                    i_element.GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_X, 0) = 0.0;
                    i_element.GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_Y, 0) = 0.0;
                    if (num_nodes_per_elem - 1 > 2)
                        i_element.GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_Z, 0) = 0.0;
                    i_element.GetGeometry()[j].FastGetSolutionStepValue(PRESSURE, 0) = 0.0;
                    i_element.GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_X, 1) = 0.0;
                    i_element.GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_Y, 1) = 0.0;
                    if (num_nodes_per_elem - 1 > 2)
                        i_element.GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_Z, 1) = 0.0;
                    i_element.GetGeometry()[j].FastGetSolutionStepValue(PRESSURE, 1) = 0.0;
                    if(GetInside)
                        vector_of_node_ids.push_back(i_element.GetGeometry()[j].Id());
                }
            }
            else
            {
                if(!GetInside){
                    count++;
                    std::size_t num_nodes_per_elem = i_element.GetGeometry().PointsNumber(); // Size()
                    rModifiedModelPart.AddElement(rModelPart.pGetElement(i_element.Id()));//AddElement()
                    for (j = 0; j < num_nodes_per_elem; j++)
                        vector_of_node_ids.push_back(i_element.GetGeometry()[j].Id());
                }
            }
        }

        //rModifiedModelPart.Nodes() = rModelPart.Nodes();

        KRATOS_INFO("Number of elements added to the modified patch") << count << std::endl;
        //sorting and making unique list of node ids
        std::set<std::size_t> s(vector_of_node_ids.begin(), vector_of_node_ids.end());
        vector_of_node_ids.assign(s.begin(), s.end());

        // Add unique nodes in the ModelPart
        for (auto i_node_id = vector_of_node_ids.begin(); i_node_id != vector_of_node_ids.end(); i_node_id++)
        {
            Node<3>::Pointer pnode = rModelPart.Nodes()(*i_node_id);
            rModifiedModelPart.AddNode(pnode);
        }

        KRATOS_CATCH("");
    }



    /**
     * @brief Extracts the outside surface/edges of a modelpart.This uses the flag CHIMERA_INTERNAL_BOUNDARY
     *                  to check if there is an internal boundary in the given ModelPart. The flag GetInternal
     *                  specifies weather to get the internal boundary marked by CHIMERA_INTERNAL_BOUNDARY or the outside one.
     * @param rVolumeModelPart The modelpart on which the boundary is to be found.
     * @param rExtractedBoundaryModelPart The extracted surface/edge modelpart.
     * @param GetInternal A bool sepecifying which surface/edge extracted. The one marked by CHIMERA_INTERNAL_BOUNDARY or the outside one.
     */
    //
    void ExtractBoundaryMesh( ModelPart &rVolumeModelPart, ModelPart &rExtractedBoundaryModelPart, bool GetInternal = false)
    {
        KRATOS_TRY;

        KRATOS_INFO(":: [Boundary mesh extraction]::") << std::endl;

        std::size_t n_nodes = rVolumeModelPart.ElementsBegin()->GetGeometry().size();
        KRATOS_ERROR_IF(!(n_nodes!=3 || n_nodes!=4))<<"Hole cutting process is only supported for tetrahedral and triangular elements" <<Info()<< std::endl;

        // Some type-definitions
        typedef std::unordered_map<vector<std::size_t>, std::size_t, KeyHasher, KeyComparator> hashmap;
        typedef std::unordered_map<vector<std::size_t>, vector<std::size_t>, KeyHasher, KeyComparator> hashmap_vec;

        // Create map to ask for number of faces for the given set of node ids representing on face in the model part
        hashmap n_faces_map;

        // Fill map that counts number of faces for given set of nodes
        for (ModelPart::ElementIterator it_elem = rVolumeModelPart.ElementsBegin(); it_elem != rVolumeModelPart.ElementsEnd(); it_elem++)
        {
            Element::GeometryType::GeometriesArrayType faces = it_elem->GetGeometry().GenerateEdges();

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
        for (ModelPart::ElementIterator it_elem = rVolumeModelPart.ElementsBegin(); it_elem != rVolumeModelPart.ElementsEnd(); it_elem++)
        {
            Element::GeometryType::GeometriesArrayType faces = it_elem->GetGeometry().GenerateEdges();

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

        // Add skin faces as triangles to skin-model-part (loop over all node sets)
        std::vector<std::size_t> vector_of_node_ids;
        for (typename hashmap::const_iterator it = n_faces_map.begin(); it != n_faces_map.end(); it++)
        {
            // If given node set represents face that is not overlapping with a face of another element, add it as skin element
            if (it->second == 1)
            {
                // If skin edge is a triangle store triangle in with its original orientation in new skin model part
                if (it->first.size() == 2)
                {
                    // Getting original order is important to properly reproduce skin edge including its normal orientation
                    vector<std::size_t> original_nodes_order = ordered_skin_face_nodes_map[it->first];

                    Node<3>::Pointer pnode1 = rVolumeModelPart.Nodes()(original_nodes_order[0]);
                    Node<3>::Pointer pnode2 = rVolumeModelPart.Nodes()(original_nodes_order[1]);

                    //Storing the node ids list
                    vector_of_node_ids.push_back(original_nodes_order[0]);
                    vector_of_node_ids.push_back(original_nodes_order[1]);

                    Properties::Pointer properties = rExtractedBoundaryModelPart.rProperties()(0);
                    Condition const &rReferenceLineCondition = KratosComponents<Condition>::Get("LineCondition2D2N"); // Condition2D

                    // Skin edges are added as conditions
                    Line2D2<Node<3>> line1(pnode1, pnode2);
                    Condition::Pointer p_condition1 = rReferenceLineCondition.Create(id_condition++, line1, properties);
                    rExtractedBoundaryModelPart.Conditions().push_back(p_condition1);
                }
                // If skin face is a triangle store triangle in with its original orientation in new skin model part
                if (it->first.size() == 3)
                {
                    // Getting original order is important to properly reproduce skin face including its normal orientation
                    vector<std::size_t> original_nodes_order = ordered_skin_face_nodes_map[it->first];
                    Node<3>::Pointer pnode1 = rVolumeModelPart.Nodes()(original_nodes_order[0]);
                    Node<3>::Pointer pnode2 = rVolumeModelPart.Nodes()(original_nodes_order[1]);
                    Node<3>::Pointer pnode3 = rVolumeModelPart.Nodes()(original_nodes_order[2]);

                    //Storing the node ids list
                    vector_of_node_ids.push_back(original_nodes_order[0]);
                    vector_of_node_ids.push_back(original_nodes_order[1]);
                    vector_of_node_ids.push_back(original_nodes_order[2]);
                    Properties::Pointer properties = rExtractedBoundaryModelPart.rProperties()(0);
                    Condition const &rReferenceTriangleCondition = KratosComponents<Condition>::Get("SurfaceCondition3D3N"); // Condition3D

                    // Skin faces are added as conditions
                    Triangle3D3<Node<3>> triangle1(pnode1, pnode2, pnode3);
                    Condition::Pointer p_condition1 = rReferenceTriangleCondition.Create(id_condition++, triangle1, properties);
                    rExtractedBoundaryModelPart.Conditions().push_back(p_condition1);
                }
                // If skin face is a quadrilateral then divide in two triangles and store them with their original orientation in new skin model part
                if (it->first.size() == 4)
                {
                    // Getting original order is important to properly reproduce skin including its normal orientation
                    vector<std::size_t> original_nodes_order = ordered_skin_face_nodes_map[it->first];

                    Node<3>::Pointer pnode1 = rVolumeModelPart.Nodes()(original_nodes_order[0]);
                    Node<3>::Pointer pnode2 = rVolumeModelPart.Nodes()(original_nodes_order[1]);
                    Node<3>::Pointer pnode3 = rVolumeModelPart.Nodes()(original_nodes_order[2]);
                    Node<3>::Pointer pnode4 = rVolumeModelPart.Nodes()(original_nodes_order[3]);
                    //Storing the node ids list
                    vector_of_node_ids.push_back(original_nodes_order[0]);
                    vector_of_node_ids.push_back(original_nodes_order[1]);
                    vector_of_node_ids.push_back(original_nodes_order[2]);
                    vector_of_node_ids.push_back(original_nodes_order[3]);
                    Properties::Pointer properties = rExtractedBoundaryModelPart.rProperties()(0);
                    Condition const &rReferenceTriangleCondition = KratosComponents<Condition>::Get("SurfaceCondition3D3N"); // Condition3D

                    // Add triangle one as condition
                    Triangle3D3<Node<3>> triangle1(pnode1, pnode2, pnode3);
                    Condition::Pointer p_condition1 = rReferenceTriangleCondition.Create(id_condition++, triangle1, properties);
                    rExtractedBoundaryModelPart.Conditions().push_back(p_condition1);

                    // Add triangle two as condition
                    Triangle3D3<Node<3>> triangle2(pnode1, pnode3, pnode4);
                    Condition::Pointer p_condition2 = rReferenceTriangleCondition.Create(id_condition++, triangle2, properties);
                    rExtractedBoundaryModelPart.Conditions().push_back(p_condition2);
                }
            }
        }

        //sorting and making unique list of node ids
        std::set<std::size_t> sort_set(vector_of_node_ids.begin(), vector_of_node_ids.end());
        vector_of_node_ids.assign(sort_set.begin(), sort_set.end());

        for (const auto& i_node_id : vector_of_node_ids)
        {
            //Adding the nodes to the rExtractedBoundaryModelPart
            Node<3>::Pointer pnode = rVolumeModelPart.Nodes()(i_node_id);
            rExtractedBoundaryModelPart.AddNode(pnode);
        }

        for(auto& i_condition : rExtractedBoundaryModelPart.Conditions())
        {
            auto& geo = i_condition.GetGeometry();
            for(auto& node : geo )
                node.Set(TO_ERASE, false);
        }

        for(auto& i_condition : rExtractedBoundaryModelPart.Conditions())
        {
            auto& geo = i_condition.GetGeometry();
            bool is_internal = true;
            for(const auto& node : geo )
                is_internal = is_internal && node.Is(CHIMERA_INTERNAL_BOUNDARY);
            if(is_internal){
                if(!GetInternal){
                    i_condition.Set(TO_ERASE);
                    for(auto& node : geo )
                        node.Set(TO_ERASE);
                }
            } else {
                if(GetInternal){
                    i_condition.Set(TO_ERASE);
                    for(auto& node : geo )
                        node.Set(TO_ERASE);
                }
            }
        }

        rExtractedBoundaryModelPart.RemoveConditions(TO_ERASE);
        rExtractedBoundaryModelPart.RemoveNodes(TO_ERASE);

        KRATOS_INFO("Successful extraction of the boundary mesh") << std::endl;
        KRATOS_CATCH("");
    }


    /// Assignment operator.
    ChimeraHoleCuttingUtility &operator=(ChimeraHoleCuttingUtility const &rOther) = delete;

    /// Copy constructor.
    ChimeraHoleCuttingUtility(ChimeraHoleCuttingUtility const& rOther) = delete;

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
        return "ChimeraHoleCuttingUtility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "ChimeraHoleCuttingUtility";
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

    ///@}

}; // Class ChimeraHoleCuttingUtility

} // namespace Kratos.

#endif // KRATOS_CUSTOM_HOLE_CUTTING_PROCESS_H_INCLUDED  defined
