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

class ChimeraHoleCuttingUtility
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

    void CreateHoleAfterDistance(ModelPart &rModelPart, ModelPart &rExtractedModelPart, ModelPart &rExtractedBoundaryModelPart, const double Distance)
    {
        KRATOS_TRY;
        KRATOS_INFO("\n::[Creating Hole]::") << std::endl;
        Model &current_model = rModelPart.GetModel();
        ModelPart &r_dummy_internal_boundary = current_model.CreateModelPart("dummy_chimera_internal_boundary");
        RemoveOutOfDomainElements(rModelPart, rExtractedModelPart,1,Distance);
        FindOutsideBoundaryOfModelPartGivenInside(rExtractedModelPart, r_dummy_internal_boundary, rExtractedBoundaryModelPart);
        current_model.DeleteModelPart("dummy_chimera_internal_boundary");
        KRATOS_CATCH("");
    }

    void FindOutsideBoundaryOfModelPartGivenInside(ModelPart &rModelPart, ModelPart &rInsideBoundary, ModelPart &rExtractedBoundaryModelPart)
    {
        std::size_t n_nodes = rModelPart.ElementsBegin()->GetGeometry().size();
        KRATOS_ERROR_IF(!(n_nodes!=3 || n_nodes!=4))<<"Hole cutting process is only supported for tetrahedral and triangular elements" <<Info()<< std::endl;

        ExtractBoundaryMesh(rModelPart, rInsideBoundary, rExtractedBoundaryModelPart);
    }

    void RemoveOutOfDomainElements(ModelPart &rModelPart,
                                   ModelPart &rModifiedModelPart,
                                   const int MainDomainOrNot,
                                   const double Distance=0.0)
    {
        KRATOS_TRY;
        KRATOS_INFO("\n:: Removing Out Of Domain Patch with Inside boundary Given ::") << std::endl;
        std::vector<std::size_t> vector_of_node_ids;

        int count = 0;

        for (ModelPart::ElementsContainerType::iterator i_element = rModelPart.ElementsBegin(); i_element != rModelPart.ElementsEnd(); ++i_element)
        {
            double element_distance = 0.0;
            std::size_t numPointsOutside = 0;
            std::size_t j = 0;
            Geometry<Node<3>> &geom = i_element->GetGeometry();

            for (j = 0; j < geom.size(); j++)
            {
                element_distance = i_element->GetGeometry()[j].FastGetSolutionStepValue(DISTANCE);

                element_distance = element_distance * MainDomainOrNot;
                if (element_distance < -1*Distance)
                {
                    numPointsOutside++;
                }
            }

            /* Any node goes out of the domain means the element need to be INACTIVE , otherwise the modified patch boundary
			 wont find any nodes on background */
            if (numPointsOutside > 0)
            {
                i_element->Set(ACTIVE, false);
                Element::Pointer p_elem = *(i_element.base());
                std::size_t num_nodes_per_elem = p_elem->GetGeometry().PointsNumber();
                for (j = 0; j < num_nodes_per_elem; j++)
                {
                    p_elem->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_X, 0) = 0.0;
                    p_elem->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_Y, 0) = 0.0;
                    if (num_nodes_per_elem - 1 > 2)
                        p_elem->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_Z, 0) = 0.0;
                    p_elem->GetGeometry()[j].FastGetSolutionStepValue(PRESSURE, 0) = 0.0;
                    p_elem->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_X, 1) = 0.0;
                    p_elem->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_Y, 1) = 0.0;
                    if (num_nodes_per_elem - 1 > 2)
                        p_elem->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_Z, 1) = 0.0;
                    p_elem->GetGeometry()[j].FastGetSolutionStepValue(PRESSURE, 1) = 0.0;
                }
            }
            else
            {
                count++;
                Element::Pointer p_elem = *(i_element.base());
                std::size_t num_nodes_per_elem = p_elem->GetGeometry().PointsNumber(); // Size()
                rModifiedModelPart.Elements().push_back(p_elem);                   //AddElement()
                for (j = 0; j < num_nodes_per_elem; j++)
                    vector_of_node_ids.push_back(p_elem->GetGeometry()[j].Id());
            }
        }

        rModifiedModelPart.Nodes() = rModelPart.Nodes();

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

    //give the outside surface of a model part given inside boundary
    void ExtractBoundaryMesh( ModelPart &rVolumeModelPart, ModelPart &rInsideBoundaryModelPart, ModelPart &rExtractedBoundaryModelPart)
    {
        KRATOS_TRY;

        KRATOS_INFO(":: [Boundary mesh extraction]::") << std::endl;

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

        KRATOS_INFO("Comparing with inside boundary ") << rInsideBoundaryModelPart.Name() << std::endl;
        for (ModelPart::ConditionIterator it_cond = rInsideBoundaryModelPart.ConditionsBegin(); it_cond != rInsideBoundaryModelPart.ConditionsEnd(); ++it_cond)
        {
            vector<std::size_t> ids(it_cond->GetGeometry().size());
            for (std::size_t i = 0; i < it_cond->GetGeometry().size(); i++)
                ids[i] = it_cond->GetGeometry()[i].Id();

            std::sort(ids.begin(), ids.end());
            if (n_faces_map[ids] == 1)
                n_faces_map[ids] += 1;
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
                    Condition const &rReferenceLineCondition = KratosComponents<Condition>::Get("Condition2D");

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
                    Condition const &rReferenceTriangleCondition = KratosComponents<Condition>::Get("Condition3D");

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
                    Condition const &rReferenceTriangleCondition = KratosComponents<Condition>::Get("Condition3D");

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

        for (auto it = vector_of_node_ids.begin(); it != vector_of_node_ids.end(); it++)
        {
            //Adding the nodes to the rExtractedBoundaryModelPart
            Node<3>::Pointer pnode = rVolumeModelPart.Nodes()(*it);
            rExtractedBoundaryModelPart.AddNode(pnode);
        }
        KRATOS_INFO("Successful extraction of the boundary mesh") << std::endl;
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

    /// Assignment operator.
    ChimeraHoleCuttingUtility &operator=(ChimeraHoleCuttingUtility const &rOther);

    /// Copy constructor.
    ChimeraHoleCuttingUtility(ChimeraHoleCuttingUtility const& rOther);

    ///@}

}; // Class ChimeraHoleCuttingUtility

} // namespace Kratos.

#endif // KRATOS_CUSTOM_HOLE_CUTTING_PROCESS_H_INCLUDED  defined
