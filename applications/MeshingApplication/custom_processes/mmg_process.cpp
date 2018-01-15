// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____ 
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _ 
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

// System includes
#include <set>
#include <unordered_set>

// External includes
// The includes related with the MMG library
#include "mmg/libmmg.h"
#include "mmg/mmg2d/libmmg2d.h" 
#include "mmg/mmg3d/libmmg3d.h"
// #include "mmg/mmgs/libmmgs.h"

// Project includes
#include "custom_processes/mmg_process.h"
// We indlude the internal variable interpolation process
#include "custom_processes/nodal_values_interpolation_process.h"
#include "custom_processes/internal_variables_interpolation_process.h"
// Include the point locator
#include "utilities/binbased_fast_point_locator.h"
// Include the spatial containers needed for search
#include "spatial_containers/spatial_containers.h" // kd-tree 
#include "includes/io.h"
#include "includes/model_part_io.h"


// NOTE: The following contains the license of the MMG library
/* =============================================================================
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/

namespace Kratos
{
// The member variables related with the MMG library
MMG5_pMesh mmgMesh;
MMG5_pSol  mmgSol;
    
/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

template< unsigned int TDim>
MmgProcess<TDim>::MmgProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    )
    :mrThisModelPart(rThisModelPart),
     mThisParameters(ThisParameters)
{       
    Parameters DefaultParameters = Parameters(R"(
    {
        "filename"                             : "out",
        "framework"                            : "Eulerian",
        "internal_variables_parameters"        :
        {
            "allocation_size"                      : 1000, 
            "bucket_size"                          : 4, 
            "search_factor"                        : 2, 
            "interpolation_type"                   : "LST",
            "internal_variable_interpolation_list" :[]
        },
        "save_external_files"              : false,
        "save_mdpa_file"                   : false,
        "max_number_of_searchs"            : 1000,
        "echo_level"                       : 3,
        "step_data_size"                   : 0,
        "remesh_at_non_linear_iteration"   : false,
        "buffer_size"                      : 0
    })" );
    
    mThisParameters.ValidateAndAssignDefaults(DefaultParameters);
    
    mStdStringFilename = mThisParameters["filename"].GetString();
    mEchoLevel = mThisParameters["echo_level"].GetInt();
    
    mFilename = new char [mStdStringFilename.length() + 1];
    std::strcpy (mFilename, mStdStringFilename.c_str());
    
    mFramework = ConvertFramework(mThisParameters["framework"].GetString());
    
    mpRefElement.clear();
    mpRefCondition.clear();
}

/*************************************** EXECUTE ***********************************/
/***********************************************************************************/

template< unsigned int TDim>
void MmgProcess<TDim>::Execute()
{       
    KRATOS_TRY;
    
    const bool safe_to_file = mThisParameters["save_external_files"].GetBool();
    
    /* We restart the MMG mesh and solution */       
    InitMesh();
    
    /* We print the original model part */
    if (mEchoLevel > 0)
    {
        std::cout << "//---------------------------------------------------//" << std::endl;
        std::cout << "//---------------------------------------------------//" << std::endl;
        std::cout << "//---------------  BEFORE REMESHING   ---------------//" << std::endl;
        std::cout << "//---------------------------------------------------//" << std::endl;
        std::cout << "//---------------------------------------------------//" << std::endl << std::endl;
        
        KRATOS_WATCH(mrThisModelPart);
    }
    
    // We initialize the mesh and solution data
    InitializeMeshData();
    InitializeSolData();
    
    // Check if the number of given entities match with mesh size 
    CheckMeshData();
    
    // Save to file
    if (safe_to_file == true) SaveSolutionToFile(false);
    
    // We execute the remeshing
    ExecuteRemeshing();
    
    /* We print the resulting model part */
    if (mEchoLevel > 0)
    {
        std::cout << "//---------------------------------------------------//" << std::endl;
        std::cout << "//---------------------------------------------------//" << std::endl;
        std::cout << "//---------------   AFTER REMESHING   ---------------//" << std::endl;
        std::cout << "//---------------------------------------------------//" << std::endl;
        std::cout << "//---------------------------------------------------//" << std::endl << std::endl;
        
        KRATOS_WATCH(mrThisModelPart);
    }
    
    KRATOS_CATCH("");
}

/************************************* OPERATOR() **********************************/
/***********************************************************************************/

template< unsigned int TDim>
void MmgProcess<TDim>::operator()()
{
    Execute();
}
    
/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void MmgProcess<TDim>::InitializeMeshData()
{                
    // First we compute the colors
    std::unordered_map<int,int> nodes_colors, cond_colors, elem_colors;
    ComputeColors(nodes_colors, cond_colors, elem_colors);
    
    /////////* MESH FILE */////////
    // Build mesh in MMG5 format //
    
    // Iterate in the nodes
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    const SizeType num_nodes = nodes_array.end() - nodes_array.begin();
    
    // Iterate in the conditions
    ConditionsArrayType& conditions_array = mrThisModelPart.Conditions();
    const SizeType num_conditions = conditions_array.end() - conditions_array.begin();
    
    // Iterate in the elements
    ElementsArrayType& elements_array = mrThisModelPart.Elements();
    const SizeType num_elements = elements_array.end() - elements_array.begin();
    
    /* Manually set of the mesh */
    array_1d<SizeType, TDim - 1> num_array_elements;
    array_1d<SizeType, TDim - 1> num_array_conditions;
    if (TDim == 2)
    {
        num_array_conditions[0] = num_conditions;
        num_array_elements[0]   = num_elements;
    }
    else
    {
        // We initialize the values
        num_array_elements[0] = 0; // Tetrahedron
        num_array_elements[1] = 0; // Prisms
        
        num_array_conditions[0] = 0; // Triangles
        num_array_conditions[1] = 0; // Quadrilaterals
        
        /* Elements */
        #pragma omp parallel for
        for(SizeType i = 0; i < num_elements; ++i) 
        {
            auto it_elem = elements_array.begin() + i;
            
            if ((it_elem->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) // Tetrahedron
            {
                #pragma omp atomic
                num_array_elements[0] += 1;
            }
            else if ((it_elem->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Prism3D6) // Prisms
            {
                #pragma omp atomic
                num_array_elements[1] += 1;
            }
            else
                std::cout << "WARNING: YOUR GEOMETRY CONTAINS HEXAEDRON THAT CAN NOT BE REMESHED" << std::endl;
        }
        
        if (((num_array_elements[0] + num_array_elements[1]) < num_elements) && mEchoLevel > 0)
            std::cout << "Number of Elements: " << num_elements << " Number of Tetrahedron: " << num_array_elements[0] << " Number of Prisms: " << num_array_elements[1] << std::endl;
        
        /* Conditions */
        #pragma omp parallel for
        for(SizeType i = 0; i < num_conditions; ++i) 
        {
            auto it_cond = conditions_array.begin() + i;
            
            if ((it_cond->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) // Triangles
            {
                #pragma omp atomic
                num_array_conditions[0] += 1;
            }
            else if ((it_cond->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4)  // Quadrilaterals
            {
                #pragma omp atomic
                num_array_conditions[1] += 1;
            }
        }
    }
    
    SetMeshSize(num_nodes, num_array_elements, num_array_conditions);
    
    /* Nodes */
    // We copy the DOF from the first node (after we release, to avoid problem with previous conditions)
    mDofs = nodes_array.begin()->GetDofs();
    for (typename Node<3>::DofsContainerType::const_iterator it_dof = mDofs.begin(); it_dof != mDofs.end(); ++it_dof)
        it_dof->FreeDof();
    
    if (mFramework == Lagrangian) // NOTE: The code is repeated due to performance reasons
    {
        #pragma omp parallel for firstprivate(nodes_colors)
        for(SizeType i = 0; i < num_nodes; ++i) 
        {
            auto it_node = nodes_array.begin() + i;
            
            SetNodes(it_node->X0(), it_node->Y0(), it_node->Z0(), nodes_colors[it_node->Id()], i + 1);
            
            bool blocked = false;
            if (it_node->IsDefined(BLOCKED) == true)
                blocked = it_node->Is(BLOCKED);
            if (TDim == 3 && blocked == true)
                BlockNode(i + 1);
            
            // RESETING THE ID OF THE NODES (important for non consecutive meshes)
            it_node->SetId(i + 1);
        }
    }
    else
    {
        #pragma omp parallel for firstprivate(nodes_colors)
        for(SizeType i = 0; i < num_nodes; ++i) 
        {
            auto it_node = nodes_array.begin() + i;
            
            SetNodes(it_node->X(), it_node->Y(), it_node->Z(), nodes_colors[it_node->Id()], i + 1);
            
            bool blocked = false;
            if (it_node->IsDefined(BLOCKED) == true)
                blocked = it_node->Is(BLOCKED);
            if (TDim == 3 && blocked == true)
                BlockNode(i + 1);
            
            // RESETING THE ID OF THE NODES (important for non consecutive meshes)
            it_node->SetId(i + 1);
        }
    }
    
    /* Conditions */
    #pragma omp parallel for firstprivate(cond_colors)
    for(SizeType i = 0; i < num_conditions; ++i) 
    {
        auto it_cond = conditions_array.begin() + i;
        SetConditions(it_cond->GetGeometry(), cond_colors[it_cond->Id()], i + 1);
    }
    
    /* Elements */
    #pragma omp parallel for firstprivate(elem_colors)
    for(SizeType i = 0; i < num_elements; ++i) 
    {
        auto it_elem = elements_array.begin() + i;
        SetElements(it_elem->GetGeometry(), elem_colors[it_elem->Id()], i + 1);
    }
    
    /* We clone the first condition and element of each type (we will assume that each sub model part has just one kind of condition, in my opinion it is quite reccomended to create more than one sub model part if you have more than one element or condition) */
    // First we add the main model part
    bool to_check_cond = false;
    bool to_check_elem = false;
    if (num_conditions > 0)
    {
        const std::string type_name = (TDim == 2) ? "Condition2D2N" : "Condition3D";
        Condition const& r_clone_condition = KratosComponents<Condition>::Get(type_name);
        mpRefCondition[0] = r_clone_condition.Create(0, r_clone_condition.GetGeometry(), conditions_array.begin()->pGetProperties());
//         mpRefCondition[0] = conditions_array.begin()->Create(0, conditions_array.begin()->GetGeometry(), conditions_array.begin()->pGetProperties());
        to_check_cond = true;
    }
    if (num_elements > 0)
    {
        mpRefElement[0] = elements_array.begin()->Create(0, elements_array.begin()->GetGeometry(), elements_array.begin()->pGetProperties());
        to_check_elem = true;
    }
    // Now we iterate over the model parts
    for (auto & color_list : mColors)
    {
        const int key = color_list.first;
        
        if (((to_check_cond == false) && (to_check_elem == false)) == true) break;
        
        if (key != 0) // NOTE: key == 0 is the MainModelPart
        {
            bool cond_added = false;
            bool elem_added = false;
            
            for (auto sub_model_part_name : color_list.second)
            {      
                ModelPart& r_sub_model_part = mrThisModelPart.GetSubModelPart(sub_model_part_name); 
                
                if (to_check_cond == true)
                {
                    ConditionsArrayType& conditions_array_sub_model_part = r_sub_model_part.Conditions();
                    const SizeType num_conditions_sub_model_part = conditions_array_sub_model_part.end() - conditions_array_sub_model_part.begin();
                    
                    if (num_conditions_sub_model_part > 0)
                    {
                        mpRefCondition[key] = conditions_array_sub_model_part.begin()->Create(0, conditions_array_sub_model_part.begin()->GetGeometry(), conditions_array_sub_model_part.begin()->pGetProperties());
                        cond_added = true;
                    }
                }
                if (to_check_elem == true)
                {
                    ElementsArrayType& elements_array_sub_model_part = r_sub_model_part.Elements();
                    const SizeType num_elements_sub_model_part = elements_array_sub_model_part.end() - elements_array_sub_model_part.begin();
                    
                    if (num_elements_sub_model_part > 0)
                    {
                        mpRefElement[key] = elements_array_sub_model_part.begin()->Create(0, elements_array_sub_model_part.begin()->GetGeometry(), elements_array_sub_model_part.begin()->pGetProperties());
                        elem_added = true;
                    }
                }
                
                if ((cond_added && elem_added) == true) break;
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim>
void MmgProcess<TDim>::InitializeSolData()
{
    ////////* SOLUTION FILE *////////
    
    // Iterate in the nodes
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    const SizeType num_nodes = nodes_array.end() - nodes_array.begin();
    
    SetSolSizeTensor(num_nodes);

    #pragma omp parallel for 
    for(SizeType i = 0; i < num_nodes; ++i) 
    {
        auto it_node = nodes_array.begin() + i;
        
    #ifdef KRATOS_DEBUG 
        KRATOS_ERROR_IF(it_node->Has(MMG_METRIC) == false) <<  " MMG_METRIC not defined for node " << it_node->Id();
    #endif     
        
        // We get the metric
        const Vector& metric = it_node->GetValue(MMG_METRIC);
        
    #ifdef KRATOS_DEBUG 
        KRATOS_ERROR_IF((metric.size() != TDim * 3 - 3) ) << "Wrong size of vector MMG_METRIC found for node " << it_node->Id() << " size is " << metric.size() << " expected size was " << TDim * 3 - 3;
    #endif
        
        // We set the metric
        SetMetricTensor(metric, i + 1);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <unsigned int TDim>
void MmgProcess<TDim>::ExecuteRemeshing()
{
    // Getting the parameters
    const bool save_to_file = mThisParameters["save_external_files"].GetBool();
    
    // We initialize some values
    const SizeType step_data_size = mrThisModelPart.GetNodalSolutionStepDataSize();
    const SizeType buffer_size   = mrThisModelPart.NodesBegin()->GetBufferSize();
    
    mThisParameters["step_data_size"].SetInt(step_data_size);
    mThisParameters["buffer_size"].SetInt(buffer_size);
    
    if (mEchoLevel > 0)      
        std::cout << "Step data size: " << step_data_size << " Buffer size: " << buffer_size << std::endl; 
    
    ////////* MMG LIBRARY CALL *////////
    if (mEchoLevel > 0)
        std::cout << "////////* MMG LIBRARY CALL *////////" << std::endl; 
    
    MMGLibCall();
    
    const unsigned int n_nodes = mmgMesh->np;
    array_1d<unsigned int, 2> n_conditions;
    if (TDim == 2)
    {
        n_conditions[0] = mmgMesh->na;
        n_conditions[1] = 0;
    }
    else
    {
        n_conditions[0] = mmgMesh->nt;
        n_conditions[1] = mmgMesh->nquad;
    }
    array_1d<unsigned int, 2> n_elements;
    if (TDim == 2)
    {
        n_elements[0] = mmgMesh->nt;
        n_elements[1] = 0;
    }
    else
    {
        n_elements[0] = mmgMesh->ne;
        n_elements[1] = mmgMesh->nprism;
    }
    
    if (mEchoLevel > 0)
    {
        std::cout << "     Nodes created: " << n_nodes << std::endl;
        if (TDim == 2) // 2D
        {
            std::cout << "Conditions created: " << n_conditions[0] << std::endl;
            std::cout << "Elements created: " << n_elements[0] << std::endl;
        }
        else // 3D
        {
            std::cout << "Conditions created: " << n_conditions[0] + n_conditions[1] << std::endl;
            std::cout << "\tTriangles: " << n_conditions[0] << "\tQuadrilaterals: " << n_conditions[1]<< std::endl;
            std::cout << "Elements created: " << n_elements[0] + n_elements[1] << std::endl;
            std::cout << "\tTetrahedron: " << n_elements[0] << "\tPrisms: " << n_elements[1] << std::endl;
        }
    }
    
    ////////* EMPTY AND BACKUP THE MODEL PART *////////
    
    ModelPart r_old_model_part;
    
    // First we empty the model part
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    const SizeType num_nodes = nodes_array.end() - nodes_array.begin();
    
    #pragma omp parallel for 
    for(SizeType i = 0; i < num_nodes; ++i) 
    {
        auto it_node = nodes_array.begin() + i;
        it_node->Set(TO_ERASE, true);
    }
    r_old_model_part.AddNodes( mrThisModelPart.NodesBegin(), mrThisModelPart.NodesEnd() );
    mrThisModelPart.RemoveNodesFromAllLevels(TO_ERASE);  
    
    ConditionsArrayType& conditions_array = mrThisModelPart.Conditions();
    const SizeType num_conditions = conditions_array.end() - conditions_array.begin();
    
    #pragma omp parallel for 
    for(SizeType i = 0; i < num_conditions; ++i) 
    {
        auto it_cond = conditions_array.begin() + i;
        it_cond->Set(TO_ERASE, true);
    }
    r_old_model_part.AddConditions( mrThisModelPart.ConditionsBegin(), mrThisModelPart.ConditionsEnd() );
    mrThisModelPart.RemoveConditionsFromAllLevels(TO_ERASE); 
    
    ElementsArrayType& elements_array = mrThisModelPart.Elements();
    const SizeType num_elements = elements_array.end() - elements_array.begin();
    
    #pragma omp parallel for 
    for(SizeType i = 0; i < num_elements; ++i) 
    {
        auto it_elem = elements_array.begin() + i;
        it_elem->Set(TO_ERASE, true);
    }
    r_old_model_part.AddElements( mrThisModelPart.ElementsBegin(), mrThisModelPart.ElementsEnd() );
    mrThisModelPart.RemoveElementsFromAllLevels(TO_ERASE);  
    
    // Create a new model part // TODO: Use a different kind of element for each submodelpart (in order to be able of remeshing more than one kind o element or condition)
    std::unordered_map<int, std::vector<IndexType>>color_nodes, color_cond_0, color_cond_1, color_elem_0, color_elem_1;
    
    /* NODES */ // TODO: ADD OMP
    for (unsigned int i_node = 1; i_node <= n_nodes; ++i_node)
    {
        int ref, is_required;
        NodeType::Pointer p_node = CreateNode(i_node, ref, is_required);
        
        // Set the DOFs in the nodes 
        for (typename NodeType::DofsContainerType::const_iterator it_dof = mDofs.begin(); it_dof != mDofs.end(); ++it_dof)
            p_node->pAddDof(*it_dof);
        
        if (ref != 0) color_nodes[ref].push_back(i_node);// NOTE: ref == 0 is the MainModelPart
    }
    
    // Auxiliar values
    int prop_id, is_required;
    
    /* CONDITIONS */ // TODO: ADD OMP
    if (mpRefCondition.size() > 0)
    {
        unsigned int cond_id = 1;
        
        unsigned int counter_cond_0 = 0;
        const std::vector<unsigned int> condition_to_remove_0 = CheckConditions0();
        for (unsigned int i_cond = 1; i_cond <= n_conditions[0]; ++i_cond)
        {
            bool skip_creation = false;
            if (counter_cond_0 < condition_to_remove_0.size())
            {
                if (condition_to_remove_0[counter_cond_0] == i_cond)
                {
                    skip_creation = true;
                    counter_cond_0 += 1;
                }
            }
            ConditionType::Pointer p_condition = CreateCondition0(cond_id, prop_id, is_required, skip_creation);
            
            if (p_condition != nullptr)
            {
                mrThisModelPart.AddCondition(p_condition);
                if (prop_id != 0) color_cond_0[prop_id].push_back(cond_id);// NOTE: prop_id == 0 is the MainModelPart
                cond_id += 1;
            }
        }
        
        unsigned int counter_cond_1 = 0;
        const std::vector<unsigned int> condition_to_remove_1 = CheckConditions1();
        for (unsigned int i_cond = 1; i_cond <= n_conditions[1]; ++i_cond)
        {                    
            bool skip_creation = false;
            if (counter_cond_1 < condition_to_remove_1.size())
            {
                if (condition_to_remove_1[counter_cond_1] == i_cond)
                {
                    skip_creation = true;
                    counter_cond_1 += 1;
                }
            }
            ConditionType::Pointer p_condition = CreateCondition1(cond_id, prop_id, is_required, skip_creation);
            
            if (p_condition != nullptr)
            {
                mrThisModelPart.AddCondition(p_condition);
                if (prop_id != 0) color_cond_1[prop_id].push_back(cond_id);// NOTE: prop_id == 0 is the MainModelPart
                cond_id += 1;
            }
        }
    }
    
    /* ELEMENTS */ // TODO: ADD OMP
    if (mpRefElement.size() > 0)
    {
        unsigned int elem_id = 1;
        
        unsigned int counter_elem_0 = 0;
        const std::vector<unsigned int> elements_to_remove_0 = CheckElements0();
        for (unsigned int i_elem = 1; i_elem <= n_elements[0]; ++i_elem)
        {  
            bool skip_creation = false;
            if (counter_elem_0 < elements_to_remove_0.size())
            {
                if (elements_to_remove_0[counter_elem_0] == i_elem)
                {
                    skip_creation = true;
                    counter_elem_0 += 1;
                }
            }
            
            ElementType::Pointer p_element = CreateElement0(elem_id, prop_id, is_required, skip_creation);
            
            if (p_element != nullptr)
            {
                mrThisModelPart.AddElement(p_element);
                if (prop_id != 0) color_elem_0[prop_id].push_back(elem_id);// NOTE: prop_id == 0 is the MainModelPart
                elem_id += 1;
            }
        }
        
        unsigned int counter_elem_1 = 0;
        const std::vector<unsigned int> elements_to_remove_1 = CheckElements1();
        for (unsigned int i_elem = 1; i_elem <= n_elements[1]; ++i_elem)
        {
            bool skip_creation = false;  
            if (counter_elem_1 < elements_to_remove_1.size())
            {
                if (elements_to_remove_1[counter_elem_1] == i_elem)
                {
                    skip_creation = true;
                    counter_elem_1 += 1;
                }
            }
            
            ElementType::Pointer p_element = CreateElement1(elem_id, prop_id, is_required,skip_creation);
            
            if (p_element != nullptr)
            {
                mrThisModelPart.AddElement(p_element);
                if (prop_id != 0) color_elem_1[prop_id].push_back(elem_id);// NOTE: prop_id == 0 is the MainModelPart
                elem_id += 1;
            }
        }
    }
    
    // We add nodes, conditions and elements to the sub model parts
    for (auto & color_list : mColors)
    {
        const int key = color_list.first;
        
        if (key != 0) // NOTE: key == 0 is the MainModelPart
        {
            for (auto sub_model_part_name : color_list.second)
            {      
                ModelPart& r_sub_model_part = mrThisModelPart.GetSubModelPart(sub_model_part_name);
                
                if (color_nodes.find(key) != color_nodes.end()) r_sub_model_part.AddNodes(color_nodes[key]);
                if (color_cond_0.find(key) != color_cond_0.end()) r_sub_model_part.AddConditions(color_cond_0[key]);
                if (color_cond_1.find(key) != color_cond_1.end()) r_sub_model_part.AddConditions(color_cond_1[key]);
                if (color_elem_0.find(key) != color_elem_0.end()) r_sub_model_part.AddElements(color_elem_0[key]);
                if (color_elem_1.find(key) != color_elem_1.end()) r_sub_model_part.AddElements(color_elem_1[key]);
            }
        }
    }
    
    // TODO: Add OMP
    // NOTE: We add the nodes from the elements and conditions to the respective submodelparts
    const std::vector<std::string> sub_model_part_names = mrThisModelPart.GetSubModelPartNames();

    for (auto sub_model_part_name : sub_model_part_names)
    {
        ModelPart& r_sub_model_part = mrThisModelPart.GetSubModelPart(sub_model_part_name);
        
        std::unordered_set<IndexType> node_ids;
        
        ConditionsArrayType& sub_conditions_array = r_sub_model_part.Conditions();
        const SizeType sub_num_conditions = sub_conditions_array.end() - sub_conditions_array.begin();
        
        for(IndexType i = 0; i < sub_num_conditions; ++i) 
        {
            auto it_cond = sub_conditions_array.begin() + i;
            auto& cond_geom = it_cond->GetGeometry();
            
            for (SizeType i_node = 0; i_node < cond_geom.size(); ++i_node)
                node_ids.insert(cond_geom[i_node].Id());
        }
        
        ElementsArrayType& sub_elements_array = r_sub_model_part.Elements();
        const SizeType sub_num_elements = sub_elements_array.end() - sub_elements_array.begin();
        
        for(IndexType i = 0; i < sub_num_elements; ++i) 
        {
            auto it_elem = sub_elements_array.begin() + i;
            auto& elem_geom = it_elem->GetGeometry();
            
            for (SizeType i_node = 0; i_node < elem_geom.size(); ++i_node)
                node_ids.insert(elem_geom[i_node].Id());
        }
        
        std::vector<IndexType> vector_ids;
        std::copy(node_ids.begin(), node_ids.end(), std::back_inserter(vector_ids));
        r_sub_model_part.AddNodes(vector_ids);
    }
    
    /* Save to file */
    if (save_to_file == true) SaveSolutionToFile(true);

    /* Free memory */
    FreeMemory();
    
    /* After that we reorder nodes, conditions and elements: */
    ReorderAllIds();
    
    /* Unmoving the original mesh to be able to interpolate */
    if (mFramework == Lagrangian) 
    {
        NodesArrayType& old_nodes_array = r_old_model_part.Nodes();
        const int old_num_nodes = static_cast<int>(old_nodes_array.size());
        
        #pragma omp parallel for
        for(int i = 0; i < old_num_nodes; ++i)
        {
            auto it_node = old_nodes_array.begin() + i;

            noalias(it_node->Coordinates()) = it_node->GetInitialPosition().Coordinates();
        }
    }
    
    /* We interpolate all the values */
    Parameters InterpolateParameters = Parameters(R"({"echo_level": 1, "framework": "Eulerian", "max_number_of_searchs": 1000, "step_data_size": 0, "buffer_size": 0})" );
    InterpolateParameters["echo_level"].SetInt(mThisParameters["echo_level"].GetInt());
    InterpolateParameters["framework"].SetString(mThisParameters["framework"].GetString());
    InterpolateParameters["max_number_of_searchs"].SetInt(mThisParameters["max_number_of_searchs"].GetInt());
    InterpolateParameters["step_data_size"].SetInt(mThisParameters["step_data_size"].GetInt());
    InterpolateParameters["buffer_size"].SetInt(mThisParameters["buffer_size"].GetInt());
    NodalValuesInterpolationProcess<TDim> InterpolateNodalValues = NodalValuesInterpolationProcess<TDim>(r_old_model_part, mrThisModelPart, InterpolateParameters);
    InterpolateNodalValues.Execute();
    
    /* We initialize elements and conditions */
    InitializeElementsAndConditions();
    
    /* We do some operations related with the Lagrangian framework */
    if (mFramework == Lagrangian) 
    {
        // If we remesh during non linear iteration we just move to the previous displacement, to the last displacement otherwise
        const int step = mThisParameters["remesh_at_non_linear_iteration"].GetBool() ? 1 : 0;
        
        /* We move the mesh */
        nodes_array = mrThisModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size());

        #pragma omp parallel for
        for(int i = 0; i < num_nodes; ++i)
        {
            auto it_node = nodes_array.begin() + i;

            noalias(it_node->Coordinates())  = it_node->GetInitialPosition().Coordinates();
            noalias(it_node->Coordinates()) += it_node->FastGetSolutionStepValue(DISPLACEMENT, step);
        }
        
        /* We interpolate the internal variables */
        InternalVariablesInterpolationProcess InternalVariablesInterpolation = InternalVariablesInterpolationProcess(r_old_model_part, mrThisModelPart, mThisParameters["internal_variables_parameters"]);
        InternalVariablesInterpolation.Execute();
    }
}

/***********************************************************************************/
/***********************************************************************************/
    
template<unsigned int TDim>
void MmgProcess<TDim>::ReorderAllIds()
{
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    const SizeType num_nodes = nodes_array.end() - nodes_array.begin();

    for(SizeType i = 0; i < num_nodes; ++i) 
    {
        auto it_node = nodes_array.begin() + i;
        it_node->SetId(i + 1);
    }

    ConditionsArrayType& condition_array = mrThisModelPart.Conditions();
    const SizeType num_conditions = condition_array.end() - condition_array.begin();
    
    for(SizeType i = 0; i < num_conditions; ++i) 
    {
        auto it_condition = condition_array.begin() + i;
        it_condition->SetId(i + 1);
    }

    ElementsArrayType& element_array = mrThisModelPart.Elements();
    const SizeType num_elements = element_array.end() - element_array.begin();

    for(SizeType i = 0; i < num_elements; ++i) 
    {
        auto it_element = element_array.begin() + i;
        it_element->SetId(i + 1);
    }
}

/***********************************************************************************/
/***********************************************************************************/
    
template<unsigned int TDim>
void MmgProcess<TDim>::InitializeElementsAndConditions()
{
    ConditionsArrayType& condition_array = mrThisModelPart.Conditions();
    const SizeType num_conditions = condition_array.end() - condition_array.begin();
    
    for(SizeType i = 0; i < num_conditions; ++i) 
    {
        auto it_condition = condition_array.begin() + i;
        it_condition->Initialize();
    }

    ElementsArrayType& element_array = mrThisModelPart.Elements();
    const SizeType num_elements = element_array.end() - element_array.begin();

    for(SizeType i = 0; i < num_elements; ++i) 
    {
        auto it_element = element_array.begin() + i;
        it_element->Initialize();
    }
}

/***********************************************************************************/
/***********************************************************************************/
    
template<unsigned int TDim>
std::vector<unsigned int> MmgProcess<TDim>::CheckNodes()
{
    typedef std::unordered_map<std::vector<double>, unsigned int, KeyHasherRange<std::vector<double>>, KeyComparorRange<std::vector<double>> > HashMap;
    HashMap node_map;
    
    std::vector<unsigned int> nodes_to_remove_ids;
    
    std::vector<double> coords(TDim);
    
    NodesArrayType& nodes_array = mrThisModelPart.Nodes();
    const SizeType num_nodes = nodes_array.end() - nodes_array.begin();
    
    for(SizeType i = 0; i < num_nodes; ++i) 
    {
        auto it_node = nodes_array.begin() + i;
        
        const array_1d<double, 3> coordinates = it_node->Coordinates();
        
        for(unsigned int i_coord = 0; i_coord < TDim; i_coord++)
            coords[i_coord] = coordinates[i_coord];
        
        node_map[coords] += 1;
        
        if (node_map[coords] > 1)
        {
            nodes_to_remove_ids.push_back(it_node->Id());
            if (mEchoLevel > 0)
                std::cout << "The mode " << it_node->Id() <<  " is repeated"<< std::endl;
        }
    }
    
    return nodes_to_remove_ids;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
std::vector<unsigned int> MmgProcess<2>::CheckConditions0()
{
    typedef std::unordered_map<std::vector<unsigned int>, unsigned int, KeyHasherRange<std::vector<unsigned int>>, KeyComparorRange<std::vector<unsigned int>> > HashMap;
    HashMap edge_map;

    std::vector<unsigned int> ids(2);

    std::vector<unsigned int> conditions_to_remove;
    
    // Iterate in the conditions
    for(int i = 0; i < mmgMesh->na; ++i) 
    {
        int edge_0, edge_1, prop_id, is_ridge, is_required;
        
        if (MMG2D_Get_edge(mmgMesh, &edge_0, &edge_1, &prop_id, &is_ridge, &is_required) != 1 )
            exit(EXIT_FAILURE);
        
        ids[0] = edge_0;
        ids[1] = edge_1;

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids.begin(), ids.end());

        edge_map[ids] += 1;
        
        if (edge_map[ids] > 1)
            conditions_to_remove.push_back(i + 1);
    }
    
    return conditions_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
std::vector<unsigned int> MmgProcess<3>::CheckConditions0()
{
    typedef std::unordered_map<std::vector<unsigned int>, unsigned int, KeyHasherRange<std::vector<unsigned int>>, KeyComparorRange<std::vector<unsigned int>> > HashMap;
    HashMap triangle_map;

    std::vector<unsigned int> ids_triangles(3);

    std::vector<unsigned int> conditions_to_remove;
            
    for(int i = 0; i < mmgMesh->nt; ++i) 
    {
        int vertex_0, vertex_1, vertex_2, prop_id, is_required;

        if (MMG3D_Get_triangle(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &prop_id, &is_required) != 1 )
            exit(EXIT_FAILURE);

        ids_triangles[0] = vertex_0;
        ids_triangles[1] = vertex_1;
        ids_triangles[2] = vertex_2;
        
        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids_triangles.begin(), ids_triangles.end());

        triangle_map[ids_triangles] += 1;
        
        if (triangle_map[ids_triangles] > 1)
            conditions_to_remove.push_back(i + 1);
    }
    
    return conditions_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
std::vector<unsigned int> MmgProcess<2>::CheckConditions1()
{
    std::vector<unsigned int> conditions_to_remove(0);
    
    return conditions_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
std::vector<unsigned int> MmgProcess<3>::CheckConditions1()
{
    typedef std::unordered_map<std::vector<unsigned int>, unsigned int, KeyHasherRange<std::vector<unsigned int>>, KeyComparorRange<std::vector<unsigned int>> > HashMap;
    HashMap quadrilateral_map;

    std::vector<unsigned int> ids_quadrialteral(4);

    std::vector<unsigned int> conditions_to_remove;
            
    for(int i = 0; i < mmgMesh->nquad; ++i) 
    {
        int vertex_0, vertex_1, vertex_2, vertex_3, prop_id, is_required;

        if (MMG3D_Get_quadrilateral(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &prop_id, &is_required) != 1 )
            exit(EXIT_FAILURE);

        ids_quadrialteral[0] = vertex_0;
        ids_quadrialteral[1] = vertex_1;
        ids_quadrialteral[2] = vertex_2;
        ids_quadrialteral[3] = vertex_3;
        
        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids_quadrialteral.begin(), ids_quadrialteral.end());

        quadrilateral_map[ids_quadrialteral] += 1;
        
        if (quadrilateral_map[ids_quadrialteral] > 1)
            conditions_to_remove.push_back(i + 1);
    }
    
    return conditions_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
std::vector<unsigned int> MmgProcess<2>::CheckElements0()
{
    typedef std::unordered_map<std::vector<unsigned int>, unsigned int, KeyHasherRange<std::vector<unsigned int>>, KeyComparorRange<std::vector<unsigned int>> > HashMap;
    HashMap triangle_map;

    std::vector<unsigned int> ids_triangles(3);

    std::vector<unsigned int> elements_to_remove;
    
    // Iterate in the elements
    for(int i = 0; i < mmgMesh->nt; ++i) 
    {
        int vertex_0, vertex_1, vertex_2, prop_id, is_required;
        
        if (MMG2D_Get_triangle(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &prop_id, &is_required) != 1 )
            exit(EXIT_FAILURE);
        
        ids_triangles[0] = vertex_0;
        ids_triangles[1] = vertex_1;
        ids_triangles[2] = vertex_2;

        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids_triangles.begin(), ids_triangles.end());

        triangle_map[ids_triangles] += 1;
        
        if (triangle_map[ids_triangles] > 1)
            elements_to_remove.push_back(i + 1);
    }
    
    return elements_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
std::vector<unsigned int> MmgProcess<3>::CheckElements0()
{
    typedef std::unordered_map<std::vector<unsigned int>, unsigned int, KeyHasherRange<std::vector<unsigned int>>, KeyComparorRange<std::vector<unsigned int>> > HashMap;
    HashMap triangle_map;

    std::vector<unsigned int> ids_tetrahedron(4);

    std::vector<unsigned int> elements_to_remove;
            
    for(int i = 0; i < mmgMesh->ne; ++i) 
    {
        int vertex_0, vertex_1, vertex_2, vertex_3, prop_id, is_required;

        if (MMG3D_Get_tetrahedron(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &prop_id, &is_required) != 1 )
            exit(EXIT_FAILURE);

        ids_tetrahedron[0] = vertex_0;
        ids_tetrahedron[1] = vertex_1;
        ids_tetrahedron[2] = vertex_2;
        ids_tetrahedron[3] = vertex_3;
        
        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids_tetrahedron.begin(), ids_tetrahedron.end());

        triangle_map[ids_tetrahedron] += 1;
        
        if (triangle_map[ids_tetrahedron] > 1)
            elements_to_remove.push_back(i + 1);
    }
    
    return elements_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
std::vector<unsigned int> MmgProcess<2>::CheckElements1()
{
    std::vector<unsigned int> elements_to_remove(0);
    return elements_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
std::vector<unsigned int> MmgProcess<3>::CheckElements1()
{
    typedef std::unordered_map<std::vector<unsigned int>, unsigned int, KeyHasherRange<std::vector<unsigned int>>, KeyComparorRange<std::vector<unsigned int>> > HashMap;
    HashMap prism_map;

    std::vector<unsigned int> ids_prisms(6);

    std::vector<unsigned int> elements_to_remove;
            
    for(int i = 0; i < mmgMesh->nprism; ++i) 
    {
        int vertex_0, vertex_1, vertex_2, vertex_3, vertex_4, vertex_5, prop_id, is_required;

        if (MMG3D_Get_prism(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &vertex_4, &vertex_5, &prop_id, &is_required) != 1 )
            exit(EXIT_FAILURE);

        ids_prisms[0] = vertex_0;
        ids_prisms[1] = vertex_1;
        ids_prisms[2] = vertex_2;
        ids_prisms[3] = vertex_3;
        ids_prisms[4] = vertex_4;
        ids_prisms[5] = vertex_5;
        
        //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
        std::sort(ids_prisms.begin(), ids_prisms.end());

        prism_map[ids_prisms] += 1;
        
        if (prism_map[ids_prisms] > 1)
            elements_to_remove.push_back(i + 1);
    }
    
    return elements_to_remove;
}

/***********************************************************************************/
/***********************************************************************************/

// template<>  // NOTE: Not yet avalaible in the official API
// void MmgProcess<2>::BlockNode(unsigned int iNode)
// {
//     if (MMG2D_Set_requiredVertex(mmgMesh, iNode) != 1 )
//         exit(EXIT_FAILURE);
// }

/***********************************************************************************/
/***********************************************************************************/


template<>  
void MmgProcess<3>::BlockNode(unsigned int iNode)
{
    if (MMG3D_Set_requiredVertex(mmgMesh, iNode) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
NodeType::Pointer MmgProcess<2>::CreateNode(
    unsigned int iNode,
    int& Ref, 
    int& IsRequired
    )
{
    double coord_0, coord_1;
    int is_corner;
    
    if (MMG2D_Get_vertex(mmgMesh, &coord_0, &coord_1, &Ref, &is_corner, &IsRequired) != 1 )
        exit(EXIT_FAILURE);
    
    NodeType::Pointer p_node = mrThisModelPart.CreateNewNode(iNode, coord_0, coord_1, 0.0);
    
    return p_node;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
NodeType::Pointer MmgProcess<3>::CreateNode(
    unsigned int iNode,
    int& Ref, 
    int& IsRequired
    )
{
    double coord_0, coord_1, coord_2;
    int is_corner;
    
    if (MMG3D_Get_vertex(mmgMesh, &coord_0, &coord_1, &coord_2, &Ref, &is_corner, &IsRequired) != 1 )
        exit(EXIT_FAILURE);
    
    NodeType::Pointer p_node = mrThisModelPart.CreateNewNode(iNode, coord_0, coord_1, coord_2);
    
    return p_node;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
ConditionType::Pointer MmgProcess<2>::CreateCondition0(        
    const unsigned int CondId,
    int& PropId, 
    int& IsRequired, 
    bool SkipCreation
    )
{
    ConditionType::Pointer p_condition = nullptr;
    
    int edge_0, edge_1, is_ridge;
    
    if (MMG2D_Get_edge(mmgMesh, &edge_0, &edge_1, &PropId, &is_ridge, &IsRequired) != 1 )
        exit(EXIT_FAILURE);
    
    // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
    if (edge_0 == 0) SkipCreation = true;
    if (edge_1 == 0) SkipCreation = true;
    
    if (SkipCreation == false)
    {
        std::vector<NodeType::Pointer> condition_nodes (2);
        condition_nodes[0] = mrThisModelPart.pGetNode(edge_0);
        condition_nodes[1] = mrThisModelPart.pGetNode(edge_1);    
        
        p_condition = mpRefCondition[PropId]->Create(CondId, condition_nodes, mpRefCondition[PropId]->pGetProperties());
    }
    else if (mEchoLevel > 2)
        std::cout << "Condition creation avoided" << std::endl;
    
    return p_condition;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
ConditionType::Pointer MmgProcess<3>::CreateCondition0(
    const unsigned int CondId,
    int& PropId, 
    int& IsRequired,
    bool SkipCreation
    )
{
    ConditionType::Pointer p_condition = nullptr;
    
    int vertex_0, vertex_1, vertex_2;

    if (MMG3D_Get_triangle(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &PropId, &IsRequired) != 1 )
        exit(EXIT_FAILURE);
    
    // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
    if (vertex_0 == 0) SkipCreation = true;
    if (vertex_1 == 0) SkipCreation = true;
    if (vertex_2 == 0) SkipCreation = true;
    
    if (SkipCreation == false)
    {
        std::vector<NodeType::Pointer> condition_nodes (3);
        condition_nodes[0] = mrThisModelPart.pGetNode(vertex_0);
        condition_nodes[1] = mrThisModelPart.pGetNode(vertex_1);
        condition_nodes[2] = mrThisModelPart.pGetNode(vertex_2);
        
        p_condition = mpRefCondition[PropId]->Create(CondId, condition_nodes, mpRefCondition[PropId]->pGetProperties());
    }
    else if (mEchoLevel > 2)
        std::cout << "Condition creation avoided" << std::endl;
    
    return p_condition;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
ConditionType::Pointer MmgProcess<2>::CreateCondition1(
    const unsigned int CondId,
    int& PropId, 
    int& IsRequired,
    bool SkipCreation
    )
{
    return nullptr;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
ConditionType::Pointer MmgProcess<3>::CreateCondition1(
    const unsigned int CondId,
    int& PropId, 
    int& IsRequired,
    bool SkipCreation
    )
{
    ConditionType::Pointer p_condition = nullptr;
    
    int vertex_0, vertex_1, vertex_2, vertex_3;

    if (MMG3D_Get_quadrilateral(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &PropId, &IsRequired) != 1 )
        exit(EXIT_FAILURE);
    
    // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
    if (vertex_0 == 0) SkipCreation = true;
    if (vertex_1 == 0) SkipCreation = true;
    if (vertex_2 == 0) SkipCreation = true;
    if (vertex_3 == 0) SkipCreation = true;
    
    if (SkipCreation == false)
    {
        std::vector<NodeType::Pointer> condition_nodes (4);
        condition_nodes[0] = mrThisModelPart.pGetNode(vertex_0);
        condition_nodes[1] = mrThisModelPart.pGetNode(vertex_1);
        condition_nodes[2] = mrThisModelPart.pGetNode(vertex_2);
        condition_nodes[3] = mrThisModelPart.pGetNode(vertex_3);
        
        p_condition = mpRefCondition[PropId]->Create(CondId, condition_nodes, mpRefCondition[PropId]->pGetProperties());
    }
    else if (mEchoLevel > 2)
        std::cout << "Condition creation avoided" << std::endl;
    
    return p_condition;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
ElementType::Pointer MmgProcess<2>::CreateElement0(        
    const unsigned int ElemId,
    int& PropId, 
    int& IsRequired,
    bool SkipCreation
    )
{
    ElementType::Pointer p_element = nullptr;

    int vertex_0, vertex_1, vertex_2;
    
    if (MMG2D_Get_triangle(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &PropId, &IsRequired) != 1 )
        exit(EXIT_FAILURE);

    // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
    if (vertex_0 == 0) SkipCreation = true;
    if (vertex_1 == 0) SkipCreation = true;
    if (vertex_2 == 0) SkipCreation = true;
    
    if (SkipCreation == false)
    {
        std::vector<NodeType::Pointer> element_nodes (3);
        element_nodes[0] = mrThisModelPart.pGetNode(vertex_0);
        element_nodes[1] = mrThisModelPart.pGetNode(vertex_1);
        element_nodes[2] = mrThisModelPart.pGetNode(vertex_2);
        
        p_element = mpRefElement[PropId]->Create(ElemId, element_nodes, mpRefElement[PropId]->pGetProperties());
    }
    else if (mEchoLevel > 2)
        std::cout << "Element creation avoided" << std::endl;
    
    return p_element;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
ElementType::Pointer MmgProcess<3>::CreateElement0(
    const unsigned int ElemId,
    int& PropId, 
    int& IsRequired,
    bool SkipCreation
    )
{
    ElementType::Pointer p_element = nullptr;
    
    int vertex_0, vertex_1, vertex_2, vertex_3;
    
    if (MMG3D_Get_tetrahedron(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &PropId, &IsRequired) != 1 )
        exit(EXIT_FAILURE);
    
    // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
    if (vertex_0 == 0) SkipCreation = true;
    if (vertex_1 == 0) SkipCreation = true;
    if (vertex_2 == 0) SkipCreation = true;
    if (vertex_3 == 0) SkipCreation = true;
    
    if (SkipCreation == false)
    {
        std::vector<NodeType::Pointer> element_nodes (4);
        element_nodes[0] = mrThisModelPart.pGetNode(vertex_0);
        element_nodes[1] = mrThisModelPart.pGetNode(vertex_1);
        element_nodes[2] = mrThisModelPart.pGetNode(vertex_2);
        element_nodes[3] = mrThisModelPart.pGetNode(vertex_3);
        
        p_element = mpRefElement[PropId]->Create(ElemId, element_nodes, mpRefElement[PropId]->pGetProperties());
    }
    else if (mEchoLevel > 2)
        std::cout << "Element creation avoided" << std::endl;
    
    return p_element;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
ElementType::Pointer MmgProcess<2>::CreateElement1(
    const unsigned int ElemId,
    int& PropId, 
    int& IsRequired,
    bool SkipCreation
    )
{
    return nullptr;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
ElementType::Pointer MmgProcess<3>::CreateElement1(
    const unsigned int ElemId,
    int& PropId, 
    int& IsRequired,
    bool SkipCreation
    )
{
    ElementType::Pointer p_element = nullptr;
                
    int vertex_0, vertex_1, vertex_2, vertex_3, vertex_4, vertex_5;
    
    if (MMG3D_Get_prism(mmgMesh, &vertex_0, &vertex_1, &vertex_2, &vertex_3, &vertex_4, &vertex_5, &PropId, &IsRequired) != 1 )
        exit(EXIT_FAILURE);
    
    // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
    if (vertex_0 == 0) SkipCreation = true;
    if (vertex_1 == 0) SkipCreation = true;
    if (vertex_2 == 0) SkipCreation = true;
    if (vertex_3 == 0) SkipCreation = true;
    if (vertex_4 == 0) SkipCreation = true;
    if (vertex_5 == 0) SkipCreation = true;
    
    if (SkipCreation == false)
    {
        std::vector<NodeType::Pointer> element_nodes (6);
        element_nodes[0] = mrThisModelPart.pGetNode(vertex_0);
        element_nodes[1] = mrThisModelPart.pGetNode(vertex_1);
        element_nodes[2] = mrThisModelPart.pGetNode(vertex_2);
        element_nodes[3] = mrThisModelPart.pGetNode(vertex_3);
        element_nodes[4] = mrThisModelPart.pGetNode(vertex_4);
        element_nodes[5] = mrThisModelPart.pGetNode(vertex_5);
    
        p_element = mpRefElement[PropId]->Create(ElemId, element_nodes, mpRefElement[PropId]->pGetProperties());
    }
    else if (mEchoLevel > 2)
        std::cout << "Element creation avoided" << std::endl;
    
    return p_element;
}

/***********************************************************************************/
/***********************************************************************************/
    
template<unsigned int TDim>
void MmgProcess<TDim>::SaveSolutionToFile(const bool PostOutput)
{
    /* GET RESULTS */

    const unsigned int& step = mrThisModelPart.GetProcessInfo()[STEP];
    
    // Automatically save the mesh 
    OutputMesh(PostOutput, step);

    // Automatically save the solution 
    OutputSol(PostOutput, step);
    
    // Save the mesh in an .mdpa format 
    const bool save_mdpa_file = mThisParameters["save_mdpa_file"].GetBool(); 
    if(save_mdpa_file == true) OutputMdpa(); 
}

/***********************************************************************************/
/***********************************************************************************/
    
template<unsigned int TDim>
void MmgProcess<TDim>::FreeMemory()
{
    // Free the MMG structures 
    FreeAll();

    // Free filename (NOTE: Problems with more that one iteration)
//     free(mFilename);
//     mFilename = nullptr;
    
    // Free reference std::unordered_map
    mpRefElement.clear();
    mpRefCondition.clear();
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<2>::InitMesh()
{  
    mmgMesh = nullptr;
    mmgSol = nullptr;
    
    // We init the MMG mesh and sol
    MMG2D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end); 
    
    InitVerbosity();
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<3>::InitMesh()
{   
    mmgMesh = nullptr;
    mmgSol = nullptr;
    
    // We init the MMG mesh and sol
    MMG3D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end); 
    
    InitVerbosity();
}

/***********************************************************************************/
/***********************************************************************************/
    
template<unsigned int TDim>
void MmgProcess<TDim>::InitVerbosity()
{
    /* We set the MMG verbosity */
    int verbosity_mmg;
    if (mEchoLevel == 0)
        verbosity_mmg = 0;
    else if (mEchoLevel == 1)
        verbosity_mmg = 0; // NOTE: This way just the essential info from MMG will be printed, but the custom message will appear
    else if (mEchoLevel == 2)
        verbosity_mmg = 3;
    else if (mEchoLevel == 3)
        verbosity_mmg = 5;
    else
        verbosity_mmg = 10;
    
    InitVerbosityParameter(verbosity_mmg);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<2>::InitVerbosityParameter(const int& VerbosityMMG)
{  
    if ( !MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_verbose, VerbosityMMG) )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<3>::InitVerbosityParameter(const int& VerbosityMMG)
{       
    if ( !MMG3D_Set_iparameter(mmgMesh,mmgSol,MMG3D_IPARAM_verbose, VerbosityMMG) )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<2>::SetMeshSize(
    const SizeType NumNodes,
    const array_1d<SizeType, 1> NumArrayElements, 
    const array_1d<SizeType, 1> NumArrayConditions
    )
{
    //Give the size of the mesh: NumNodes vertices, num_elements triangles, num_conditions edges (2D) 
    if ( MMG2D_Set_meshSize(mmgMesh, NumNodes, NumArrayElements[0], NumArrayConditions[0]) != 1 ) 
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<3>::SetMeshSize(
    const SizeType NumNodes,
    const array_1d<SizeType, 2> NumArrayElements,  // NOTE: We do this tricky thing to take into account the prisms
    const array_1d<SizeType, 2> NumArrayConditions // NOTE: We do this tricky thing to take into account the quadrilaterals
    )
{
    //Give the size of the mesh: NumNodes Vertex, num_elements tetra and prism, NumArrayConditions triangles and quadrilaterals, 0 edges (3D) 
    if ( MMG3D_Set_meshSize(mmgMesh, NumNodes, NumArrayElements[0], NumArrayElements[1], NumArrayConditions[0], NumArrayConditions[1], 0) != 1 ) 
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<2>::SetSolSizeScalar(const int NumNodes)
{
    if ( MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,NumNodes,MMG5_Scalar) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<3>::SetSolSizeScalar(const int NumNodes)
{
    if ( MMG3D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,NumNodes,MMG5_Scalar) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<2>::SetSolSizeVector(const int NumNodes)
{
    KRATOS_ERROR << "WARNING:: Vector metric not avalaible in 2D" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<3>::SetSolSizeVector(const int NumNodes)
{
    if ( MMG3D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,NumNodes,MMG5_Vector) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<2>::SetSolSizeTensor(const int NumNodes)
{
    if ( MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,NumNodes,MMG5_Tensor) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<3>::SetSolSizeTensor(const int NumNodes)
{
    if ( MMG3D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,NumNodes,MMG5_Tensor) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<2>::CheckMeshData()
{
    if ( MMG2D_Chk_meshData(mmgMesh, mmgSol) != 1 ) 
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<3>::CheckMeshData()
{
    if ( MMG3D_Chk_meshData(mmgMesh, mmgSol) != 1 ) 
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<2>::OutputMesh(
    const bool PostOutput,
    const unsigned int Step
    )
{
    std::string mesh_name;
    if (PostOutput == true)
        mesh_name = mStdStringFilename+"_step="+std::to_string(Step)+".o.mesh";
    else
        mesh_name = mStdStringFilename+"_step="+std::to_string(Step)+".mesh";
    
    auto  mesh_file = new char [mesh_name.length() + 1];
    std::strcpy (mesh_file, mesh_name.c_str());
    
    // a)  Give the ouptut mesh name using MMG2D_Set_outputMeshName (by default, the mesh is saved in the "mesh.o.mesh" file  
    MMG2D_Set_outputMeshName(mmgMesh,mesh_file);

    // b) function calling 
    if ( MMG2D_saveMesh(mmgMesh,mesh_file) != 1) 
        std::cout << "UNABLE TO SAVE MESH" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<3>::OutputMesh(
    const bool PostOutput,
    const unsigned int Step
    )
{
    std::string mesh_name;
    if (PostOutput == true)
        mesh_name = mStdStringFilename+"_step="+std::to_string(Step)+".o.mesh";
    else
        mesh_name = mStdStringFilename+"_step="+std::to_string(Step)+".mesh";

    auto  mesh_file = new char [mesh_name.length() + 1];
    std::strcpy (mesh_file, mesh_name.c_str());
    
    // a)  Give the ouptut mesh name using MMG3D_Set_outputMeshName (by default, the mesh is saved in the "mesh.o.mesh" file 
    MMG3D_Set_outputMeshName(mmgMesh,mesh_file);

    // b) function calling 
    if ( MMG3D_saveMesh(mmgMesh,mesh_file) != 1) 
        std::cout << "UNABLE TO SAVE MESH" << std::endl;
}

/***********************************************************************************/ 
/***********************************************************************************/ 
 
template<unsigned int TDim>   
void MmgProcess<TDim>::OutputMdpa() 
{ 
    std::ofstream output_file; 
    ModelPartIO model_part_io("output", IO::WRITE); 
    model_part_io.WriteModelPart(mrThisModelPart); 
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<2>::OutputSol(
    const bool PostOutput,
    const unsigned int Step
    )
{
    std::string sol_name;
    if (PostOutput == true)
        sol_name = mStdStringFilename+"_step="+std::to_string(Step)+".o.sol";
    else
        sol_name = mStdStringFilename+"_step="+std::to_string(Step)+".sol";
    
    auto  sol_file = new char [sol_name.length() + 1];
    std::strcpy (sol_file, sol_name.c_str());
    
    // a)  Give the ouptut sol name using MMG2D_Set_outputSolName (by default, the mesh is saved in the "mesh.o.sol" file 
    MMG2D_Set_outputSolName(mmgMesh, mmgSol, sol_file);

    // b) Function calling 
    if ( MMG2D_saveSol(mmgMesh, mmgSol, sol_file) != 1) 
        std::cout << "UNABLE TO SAVE SOL" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<3>::OutputSol(
    const bool PostOutput,
    const unsigned int Step
    )
{
    std::string sol_name;
    if (PostOutput == true)
        sol_name = mStdStringFilename+"_step="+std::to_string(Step)+".o.sol";
    else
        sol_name = mStdStringFilename+"_step="+std::to_string(Step)+".sol";
    
    auto  sol_file = new char [sol_name.length() + 1];
    std::strcpy (sol_file, sol_name.c_str());
    
    // a)  Give the ouptut sol name using MMG3D_Set_outputSolName (by default, the mesh is saved in the "mesh.o.sol" file 
    MMG3D_Set_outputSolName(mmgMesh, mmgSol, sol_file);

    // b) Function calling 
    if ( MMG3D_saveSol(mmgMesh,mmgSol, sol_file) != 1) 
        std::cout << "UNABLE TO SAVE SOL" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<2>::MMGLibCall()
{
    const int ier = MMG2D_mmg2dlib(mmgMesh, mmgSol);

    if ( ier == MMG5_STRONGFAILURE ) 
        KRATOS_ERROR << "WARNING: BAD ENDING OF MMG2DLIB: UNABLE TO SAVE MESH. ier: " << ier << std::endl;
    else if ( ier == MMG5_LOWFAILURE )
        KRATOS_ERROR << "WARNING: BAD ENDING OF MMG2DLIB. ier: " << ier << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<3>::MMGLibCall()
{
    const int ier = MMG3D_mmg3dlib(mmgMesh, mmgSol);

    if ( ier == MMG5_STRONGFAILURE ) 
        KRATOS_ERROR << "WARNING: BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH. ier: " << ier << std::endl;
    else if ( ier == MMG5_LOWFAILURE )
        KRATOS_ERROR << "WARNING: BAD ENDING OF MMG3DLIB. ier: " << ier << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<2>::FreeAll()
{
    MMG2D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,MMG5_ARG_end);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<3>::FreeAll()
{
    MMG3D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,MMG5_ARG_end);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<2>::SetNodes(
    const double X,
    const double Y,
    const double Z,
    const int Color,
    const int Index
    )
{
    if ( MMG2D_Set_vertex(mmgMesh, X, Y, Color, Index) != 1 )  
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<3>::SetNodes(
    const double X,
    const double Y,
    const double Z,
    const int Color,
    const int Index
    )
{
    if ( MMG3D_Set_vertex(mmgMesh, X, Y, Z, Color, Index) != 1 )  
        exit(EXIT_FAILURE); 
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<2>::SetConditions(
    Geometry<Node<3> > & Geom,
    const int Color,
    const int Index
    )
{
    if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Point2D) // Point
        KRATOS_ERROR << "WARNING:: Nodal condition, will be meshed with the node. Condition existence after meshing not guaranteed" << std::endl;
    else if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line2D2) // Line
    {
        const int& id_1 = Geom[0].Id(); // First node id
        const int& id_2 = Geom[1].Id(); // Second node id

        if ( MMG2D_Set_edge(mmgMesh, id_1, id_2, Color, Index) != 1 ) 
        {
            exit(EXIT_FAILURE);
        }
        
        // Set fixed boundary
        bool blocked_1 = false;
        if (Geom[0].IsDefined(BLOCKED) == true)
        {
            blocked_1 = Geom[0].Is(BLOCKED);
        }
        bool blocked_2 = false;
        if (Geom[1].IsDefined(BLOCKED) == true)
        {
            blocked_2 = Geom[1].Is(BLOCKED);
        }

        if ((blocked_1 && blocked_2) == true)
        {
            if ( MMG2D_Set_requiredEdge(mmgMesh, Index) != 1 ) 
            {
                exit(EXIT_FAILURE); 
            }   
        }
    }
    else
    {
        const unsigned int size_geometry = Geom.size();
        KRATOS_ERROR << "WARNING: I DO NOT KNOW WHAT IS THIS. Size: " << size_geometry << " Type: " << Geom.GetGeometryType() << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<3>::SetConditions(
    Geometry<Node<3> > & Geom,
    const int Color,
    const int Index
    )
{
    if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Point3D) // Point
        KRATOS_ERROR << "WARNING:: Nodal condition, will be meshed with the node. Condition existence after meshing not guaranteed" << std::endl;
    else if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line3D2) // Line
    {
        KRATOS_ERROR << "Kratos_Line3D2 remeshing pending to be implemented" << std::endl;
//         const int id1 = Geom[0].Id(); // First node id
//         const int id2 = Geom[1].Id(); // Second node id
// 
//         if ( MMG3D_Set_edge(mmgMesh, id1, id2, Color, Index) != 1 ) 
//         {
//             exit(EXIT_FAILURE);
//         }
//         
//         // Set fixed boundary
//         bool blocked_1 = false;
//         if (Geom[0].IsDefined(BLOCKED) == true)
//         {
//             blocked_1 = Geom[0].Is(BLOCKED);
//         }
//         bool blocked_2 = false;
//         if (Geom[1].IsDefined(BLOCKED) == true)
//         {
//             blocked_2 = Geom[1].Is(BLOCKED);
//         }
// 
//         if ((blocked_1 && blocked_2) == true)
//         {
//             if ( MMG3D_Set_requiredEdge(mmgMesh, Index) != 1 ) 
//             {
//                 exit(EXIT_FAILURE); 
//             }   
//         }
    }
    else if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) // Triangle
    {
        const int& id_1 = Geom[0].Id(); // First node Id
        const int& id_2 = Geom[1].Id(); // Second node Id
        const int& id_3 = Geom[2].Id(); // Third node Id
    
        if ( MMG3D_Set_triangle(mmgMesh, id_1, id_2, id_3, Color, Index) != 1 )  
            exit(EXIT_FAILURE); 
        
        // Set fixed boundary
        bool blocked_1 = false;
        if (Geom[0].IsDefined(BLOCKED) == true)
            blocked_1 = Geom[0].Is(BLOCKED);
        bool blocked_2 = false;
        if (Geom[1].IsDefined(BLOCKED) == true)
            blocked_2 = Geom[1].Is(BLOCKED);
        bool blocked_3 = false;
        if (Geom[2].IsDefined(BLOCKED) == true)
            blocked_3 = Geom[2].Is(BLOCKED);
        
        if ((blocked_1 && blocked_2 && blocked_3) == true)
        {
            if ( MMG3D_Set_requiredTriangle(mmgMesh, Index) != 1 ) 
                exit(EXIT_FAILURE);  
        }
    }
    else if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4) // Quadrilaterals
    {
        const int& id_1 = Geom[0].Id(); // First node Id
        const int& id_2 = Geom[1].Id(); // Second node Id
        const int& id_3 = Geom[2].Id(); // Third node Id
        const int& id_4 = Geom[3].Id(); // Fourth node Id
        
        if ( MMG3D_Set_quadrilateral(mmgMesh, id_1, id_2, id_3, id_4, Color, Index) != 1 )  
            exit(EXIT_FAILURE); 
    }
    else
    {
        const unsigned int size_geometry = Geom.size();
        KRATOS_ERROR << "WARNING: I DO NOT KNOW WHAT IS THIS. Size: " << size_geometry << " Type: " << Geom.GetGeometryType() << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<2>::SetElements(
    Geometry<Node<3> > & Geom,
    const int Color,
    const int Index
    )
{
    const int& id_1 = Geom[0].Id(); // First node Id
    const int& id_2 = Geom[1].Id(); // Second node Id
    const int& id_3 = Geom[2].Id(); // Third node Id
    
    if ( MMG2D_Set_triangle(mmgMesh, id_1, id_2, id_3, Color, Index) != 1 ) 
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<3>::SetElements(
    Geometry<Node<3> > & Geom,
    const int Color,
    const int Index
    )
{
    const int& id_1 = Geom[0].Id(); // First node Id
    const int& id_2 = Geom[1].Id(); // Second node Id
    const int& id_3 = Geom[2].Id(); // Third node Id
    const int& id_4 = Geom[3].Id(); // Fourth node Id
    
    if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) // Tetrahedron
    {
        if ( MMG3D_Set_tetrahedron(mmgMesh, id_1, id_2, id_3, id_4, Color, Index) != 1 )  
            exit(EXIT_FAILURE); 
    }
    else if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Prism3D6) // Prisms
    {
        const int& id_5 = Geom[4].Id(); // 5th node Id
        const int& id_6 = Geom[5].Id(); // 6th node Id
        
        if ( MMG3D_Set_prism(mmgMesh, id_1, id_2, id_3, id_4, id_5, id_6, Color, Index) != 1 )  
            exit(EXIT_FAILURE); 
    }
    else if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8) // Hexaedron
    {
//         const int id_5 = Geom[4].Id(); // 5th node Id
//         const int id_6 = Geom[5].Id(); // 6th node Id
//         const int id_6 = Geom[7].Id(); // 7th node Id
//         const int id_6 = Geom[8].Id(); // 8th node Id
        
        const unsigned int size_geometry = Geom.size();
        KRATOS_ERROR << "WARNING: HEXAEDRON NON IMPLEMENTED IN THE LIBRARY " << size_geometry << std::endl;
    }
    else
    {
        const unsigned int size_geometry = Geom.size();
        KRATOS_ERROR << "WARNING: I DO NOT KNOW WHAT IS THIS. Size: " << size_geometry << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/
    
template<unsigned int TDim>
void MmgProcess<TDim>::ComputeColors(
    std::unordered_map<int,int>& NodesColors,
    std::unordered_map<int,int>& CondColors,
    std::unordered_map<int,int>& ElemColors
    )
{        
    // Initialize and create the auxiliar maps
    const std::vector<std::string> sub_model_part_names = mrThisModelPart.GetSubModelPartNames();
    std::unordered_map<int,std::set<int>> aux_nodes_colors, aux_cond_colors, aux_elem_colors;
    
    std::vector<std::string> model_part_names;
    model_part_names.push_back(mrThisModelPart.Name());
    for (const auto & sub_model_part_name : sub_model_part_names)
        model_part_names.push_back(sub_model_part_name);
    
    // Initialize Colors
    int color = 0;
    for (SizeType i_sub_model_part = 0; i_sub_model_part < model_part_names.size(); ++i_sub_model_part)
    {
        mColors[i_sub_model_part].push_back(model_part_names[i_sub_model_part]);
        
        if (color > 0)
        {            
            ModelPart& r_sub_model_part = mrThisModelPart.GetSubModelPart(model_part_names[i_sub_model_part]);
            
            // Iterate in the nodes
            NodesArrayType& nodes_array = r_sub_model_part.Nodes();
            const SizeType num_nodes = nodes_array.end() - nodes_array.begin();
            
            // Iterate in the conditions
            ConditionsArrayType& conditions_array = r_sub_model_part.Conditions();
            const SizeType num_conditions = conditions_array.end() - conditions_array.begin();
            
            // Iterate in the elements
            ElementsArrayType& elements_array = r_sub_model_part.Elements();
            const SizeType num_elements = elements_array.end() - elements_array.begin();

            /* Nodes */
            for(SizeType i = 0; i < num_nodes; ++i) 
            {
                auto it_node = nodes_array.begin() + i;
                aux_nodes_colors[it_node->Id()].insert(color);
            }
            
            /* Conditions */
            for(SizeType i = 0; i < num_conditions; ++i) 
            {
                auto it_cond = conditions_array.begin() + i;
                aux_cond_colors[it_cond->Id()].insert(color);
            }
            
            /* Elements */
            for(SizeType i = 0; i < num_elements; ++i) 
            {
                auto it_elem = elements_array.begin() + i;
                aux_elem_colors[it_elem->Id()].insert(color);
            }
        }
        
        color += 1;
    }
    
    // Now detect all the cases in which a node or a cond belongs to more than one part simultaneously 
    std::unordered_map<std::set<int>, int, KeyHasherRange<std::set<int>>, KeyComparorRange<std::set<int>> > combinations;
    
    /* Nodes */
    for(auto & aux_nodes_color : aux_nodes_colors) 
    {
        const std::set<int> value = aux_nodes_color.second;
        if (value.size() > 1) combinations[value] = -1;
    }
    
    /* Conditions */
    for(auto & aux_cond_color : aux_cond_colors) 
    {
        const std::set<int> value = aux_cond_color.second;
        if (value.size() > 1) combinations[value] = -1;
    }

    /* Elements */
    for(auto & aux_elem_color : aux_elem_colors) 
    {
        const std::set<int> value = aux_elem_color.second;
        if (value.size() > 1) combinations[value] = -1;
    }
    
    /* Combinations */
    for(auto & combination : combinations) 
    {
        const std::set<int> key = combination.first;
        for(int it : key) 
            mColors[color].push_back(mColors[it][0]);
        combinations[key] = color;
        color += 1;
    }
    
    // The final maps are created
    /* Nodes */
    for(auto & aux_nodes_color : aux_nodes_colors) 
    {
        const int key = aux_nodes_color.first;
        const std::set<int> value = aux_nodes_color.second;
        
        if (value.size() == 0)
            NodesColors[key] = 0; // Main Model Part
        else if (value.size() == 1) // Another Model Part
            NodesColors[key] = *value.begin();
        else // There is a combination
            NodesColors[key] = combinations[value];
    }
    
    /* Conditions */
    for(auto & aux_cond_color : aux_cond_colors) 
    {
        const int key = aux_cond_color.first;
        const std::set<int> value = aux_cond_color.second;
        
        if (value.size() == 0)
            CondColors[key] = 0; // Main Model Part
        else if (value.size() == 1) // Another Model Part
            CondColors[key] = *value.begin();
        else // There is a combination
            CondColors[key] = combinations[value];
    }
    
    /* Elements */
    for(auto & aux_elem_color : aux_elem_colors) 
    {
        const int key = aux_elem_color.first;
        const std::set<int> value = aux_elem_color.second;
        
        if (value.size() == 0)
            ElemColors[key] = 0; // Main Model Part
        else if (value.size() == 1) // Another Model Part
            ElemColors[key] = *value.begin();
        else // There is a combination
            ElemColors[key] = combinations[value];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<2>::SetMetricScalar(
    const double& Metric,
    const int NodeId 
    )
{
    if ( MMG2D_Set_scalarSol(mmgSol, Metric, NodeId) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<3>::SetMetricScalar(
    const double& Metric,
    const int NodeId 
    )
{
    if ( MMG3D_Set_scalarSol(mmgSol, Metric, NodeId) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<2>::SetMetricVector(
    const array_1d<double, 3>& Metric,
    const int NodeId 
    )
{
    KRATOS_ERROR << "WARNING:: Vector metric not avalaible in 2D" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<3>::SetMetricVector(
    const array_1d<double, 3>& Metric,
    const int NodeId 
    )
{
    if ( MMG3D_Set_vectorSol(mmgSol, Metric[0], Metric[1], Metric[2], NodeId) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<2>::SetMetricTensor(
    const Vector& Metric,
    const int NodeId 
    )
{
    if ( MMG2D_Set_tensorSol(mmgSol, Metric[0],  Metric[1], Metric[2], NodeId) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
void MmgProcess<3>::SetMetricTensor(
    const Vector& Metric,
    const int NodeId 
    )
{
    if ( MMG3D_Set_tensorSol(mmgSol, Metric[0], Metric[1], Metric[2], Metric[3], Metric[4], Metric[5], NodeId) != 1 )
        exit(EXIT_FAILURE);
}

/***********************************************************************************/
/***********************************************************************************/
    
template<unsigned int TDim>
FrameworkEulerLagrange MmgProcess<TDim>::ConvertFramework(const std::string& str)
{
    if(str == "Lagrangian") 
        return Lagrangian;
    else if(str == "Eulerian") 
        return Eulerian;
    else
        return Eulerian;
}

/***********************************************************************************/
/***********************************************************************************/

template class MmgProcess<2>;
template class MmgProcess<3>;

}// namespace Kratos.
