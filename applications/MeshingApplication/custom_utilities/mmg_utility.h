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

#if !defined(KRATOS_MMG_UTILITY)
#define KRATOS_MMG_UTILITY

// Project includes
#include "utilities/math_utils.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "input_output/logger.h"
#include <set>
#include <map>
// The includes related with the MMG library
#include "mmg/libmmg.h"
#include "mmg/mmg2d/libmmg2d.h" 
#include "mmg/mmg3d/libmmg3d.h"
#include "mmg/mmgs/libmmgs.h"
// Include the point locator
#include "utilities/binbased_fast_point_locator.h"

// NOTE: Inspired in the documentation from:
// https://github.com/MmgTools/mmg/blob/master/libexamples/

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
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

    enum cond_geometries_2d {Line = 0};
    
    enum elem_geometries_2d {Triangle2D = 0};
    
    enum cond_geometries_3d {Triangle3D = 0, Quadrilateral3D = 1};
    
    enum elem_geometries_3d {Tetrahedra = 0, Prism = 1};
    
///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

template< unsigned int TDim>  
class MmgUtility
{
public:

    ///@name Type Definitions
    ///@{
    
    typedef ModelPart::NodesContainerType                        NodesArrayType;
    typedef ModelPart::ElementsContainerType                  ElementsArrayType;
    typedef ModelPart::ConditionsContainerType              ConditionsArrayType;
    typedef Node <3>                                                   NodeType;
    typedef Properties                                           PropertiesType;
    typedef Element                                                 ElementType;
    typedef Condition                                             ConditionType;
    typedef std::size_t                                               IndexType;
    typedef std::size_t                                                SizeType;
    typedef Dof<double>                                                 DofType;
    typedef Mesh<NodeType, PropertiesType, ElementType, ConditionType> MeshType;
    typedef MeshType::PropertiesContainerType           PropertiesContainerType;
    typedef MeshType::NodeConstantIterator                 NodeConstantIterator;
    typedef MeshType::ConditionConstantIterator       ConditionConstantIterator;
    typedef MeshType::ElementConstantIterator           ElementConstantIterator;
    
    ///@}
    ///@name Life Cycle
    ///@{
     
    /** ------------------------------ STEP   I -------------------------- */
    // Constructor
    
    /**
     * This is the default constructor, which is used to read the input files 
     * @param Filename: The input name of the 
     */
    
    MmgUtility(
        const std::string Filename,
        const unsigned int echo_level = 3
        )
        :mStdStringFilename(Filename)
    {       
       mFilename = new char [Filename.length() + 1];
       std::strcpy (mFilename, Filename.c_str());
       
       mpRefElement.resize(TDim - 1);
       mpRefCondition.resize(TDim - 1);
       mInitRefCondition.resize(TDim - 1);
       mInitRefElement.resize(TDim - 1);
       for (unsigned int i_dim = 0; i_dim < TDim - 1; i_dim++)
       {
           mpRefElement[i_dim] = nullptr;   
           mInitRefCondition[i_dim] = false;   
           mpRefCondition[i_dim] = nullptr;   
           mInitRefElement[i_dim] = false; 
       }
       
       mmgMesh = NULL;
       mmgSol = NULL;
       
       InitMesh();
       
       int verbosityMMG; 
       if (echo_level == 0)
       {
           verbosityMMG = -10;
       }
       else if (echo_level == 1)
       {
           verbosityMMG = -5;
       }
       else if (echo_level == 2)
       {
           verbosityMMG = 0;
       }
       else if (echo_level == 3)
       {
           verbosityMMG = 5;
       }
       else
       {
           verbosityMMG = 10;
       }
       
       if ( !MMG3D_Set_iparameter(mmgMesh,mmgSol,MMG3D_IPARAM_verbose, verbosityMMG) ) // From -10 to 10
       {
           exit(EXIT_FAILURE);
       }
    }
    
    /// Destructor.
    ~MmgUtility() {}
    
    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    /** ------------------------------ STEP  II -------------------------- */
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * Instead of using an files already created we read an existing model part
     * @param rThisModelPart: The original model part, an auxiliary model part will be created an the input model part will be modified
     */
    
    void ComputeExistingModelPart(
        ModelPart& rThisModelPart,
        const Variable<double> & rVariable,
        const Variable<array_1d<double,3>> & rVariableGradient,
        const double& elementary_length = 0.1,
        const double& initial_alpha_parameter = 0.01,
        const double& value_threshold = 1.0,
        const std::string& interpolation = "Constant",
        const bool& save_to_file = false
        )
    {
        // NOTE: Step one in the constructor
        
        std::cout << "//---------------------------------------------------//" << std::endl;
        std::cout << "//---------------------------------------------------//" << std::endl;
        std::cout << "//---------------  BEFORE REMESHING   ---------------//" << std::endl;
        std::cout << "//---------------------------------------------------//" << std::endl;
        std::cout << "//---------------------------------------------------//" << std::endl;
        std::cout << std::endl;
        
        KRATOS_WATCH(rThisModelPart);
        
        // First we compute the colors
        std::map<int,int> node_colors, cond_colors, elem_colors;
        ComputeColors(rThisModelPart, node_colors, cond_colors, elem_colors);
        
        ////////* MESH FILE *////////
                
        /** 
         * 2) Build mesh in MMG5 format 
         */
        
        // Iterate in the nodes
        NodesArrayType& pNode = rThisModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        // Iterate in the conditions
        ConditionsArrayType& pConditions = rThisModelPart.Conditions();
        auto numConditions = pConditions.end() - pConditions.begin();
        
        // Iterate in the elements
        ElementsArrayType& pElements = rThisModelPart.Elements();
        auto numElements = pElements.end() - pElements.begin();
        
        /* Manually set of the mesh */
        array_1d<int, TDim - 1> numArrayElements;
        array_1d<int, TDim - 1> numArrayConditions;
        if (TDim == 2)
        {
            numArrayConditions[0] = numConditions;
            numArrayElements[0]   = numElements;
        }
        else
        {
            // We initialize the values
            numArrayElements[0] = 0; // Tetrahedron
            numArrayElements[1] = 0; // Prisms
            
            numArrayConditions[0] = 0; // Triangles
            numArrayConditions[1] = 0; // Quadrilaterals
            
            /* Elements */
            for(unsigned int i = 0; i < numElements; i++) 
            {
                auto itElem = pElements.begin() + i;
                
                const unsigned int size_geom = itElem->GetGeometry().size();
                
                if (size_geom == 4) // Tetrahedron
                {
                    numArrayElements[0] += 1;
                }
                else if (size_geom == 6) // Prisms
                {
                    numArrayElements[1] += 1;
                }
            }
            
            if ((numArrayElements[0] + numArrayElements[1]) < numElements)
            {
                std::cout << "WARNING: YOUR GEOMETRY CONTAINS HEXAEDRON THAT CAN NOT BE REMESHED" << std::endl;
                std::cout << "Number of Elements: " << numElements << " Number of Tetrahedron: " << numArrayElements[0] << " Number of Prisms: " << numArrayElements[1] << std::endl;
            }
            
            /* Conditions */
            for(unsigned int i = 0; i < numConditions; i++) 
            {
                auto itCond = pConditions.begin() + i;
                
                const unsigned int size_geom = itCond->GetGeometry().size();
                
                if (size_geom == 3) // Triangles
                {
                    numArrayConditions[0] += 1;
                }
                else if (size_geom == 4)  // Quadrilaterals
                {
                    numArrayConditions[1] += 1;
                }
            }
        }
        
        SetMeshSize(numNodes, numArrayElements, numArrayConditions);
        
        /* Nodes */
//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            // We copy the DOF from the fisrt node
            if (i == 0)
            {
                mDofs = itNode->GetDofs();
            }
            
            SetNodes(itNode->X(), itNode->Y(), itNode->Z(), node_colors[itNode->Id()], i + 1);
            
            // RESETING THE ID OF THE NODES (important for non consecutive meshes)
            itNode->SetId(i + 1);
        }
        
        /* Conditions */
//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCond = pConditions.begin() + i;
            
            // We clone the first condition of each type
            if (TDim == 2)
            {
                if (i == 0)
                {
                    mpRefCondition[0] = itCond->Create(0, itCond->GetGeometry(), itCond->pGetProperties());
                }
            }
            else
            {
                const unsigned int size_geom = itCond->GetGeometry().size();
                
                if (size_geom == 3 && mInitRefCondition[0] == false) // Triangle
                {
                    mpRefCondition[0] = itCond->Create(0, itCond->GetGeometry(), itCond->pGetProperties());
                    mInitRefCondition[0] = true;
                }
                else if (size_geom == 4 && mInitRefCondition[1] == false) // Quadrilateral
                {
                    mpRefCondition[1] = itCond->Create(0, itCond->GetGeometry(), itCond->pGetProperties());
                    mInitRefCondition[1] = true;
                }
            }
            
            SetConditions(itCond->GetGeometry(), cond_colors[itCond->Id()], i + 1);
        }
        
        /* Elements */
//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numElements; i++) 
        {
            auto itElem = pElements.begin() + i;
            
            // We clone the first element of each type
            if (TDim == 2)
            {
                if (i == 0)
                {
                    mpRefElement[0] = itElem->Create(0, itElem->GetGeometry(), itElem->pGetProperties());
                }
            }
            else
            {
                const unsigned int size_geom = itElem->GetGeometry().size();
                
                if (size_geom == 4 && mInitRefElement[0] == false) // Tetrahedra
                {
                    mpRefElement[0] = itElem->Create(0, itElem->GetGeometry(), itElem->pGetProperties());
                    mInitRefElement[0] = true;
                }
                else if (size_geom == 6 && mInitRefElement[1] == false) // Prism
                {
                    mpRefElement[1] = itElem->Create(0, itElem->GetGeometry(), itElem->pGetProperties());
                    mInitRefElement[1] = true;
                }
            }
            
            SetElements(itElem->GetGeometry(), elem_colors[itElem->Id()], i + 1);
        }
        
        ////////* SOLUTION FILE *////////
        
        SetSolSize(numNodes);

//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            const double scalar_value = itNode->FastGetSolutionStepValue(rVariable, 0);
            array_1d<double, 3> gradient_value = itNode->FastGetSolutionStepValue(rVariableGradient, 0);
            
            double element_size = itNode->FastGetSolutionStepValue(NODAL_H, 0);
            if (element_size > elementary_length)
            {
                element_size = elementary_length;
            }
            
            const double tolerance = 1.0e-12;
            const double norm = norm_2(gradient_value); // Consider isotropic when zero
            double ratio = 1.0; // NOTE: Isotropic mesh
            if (norm > tolerance)
            {
                gradient_value /= norm;
//                 if (scalar_value < value_threshold)
                if (std::abs(scalar_value) < value_threshold)
                {
                    if (interpolation.find("Constant") != std::string::npos)
                    {
                        ratio = initial_alpha_parameter;
                    }
                    else if (interpolation.find("Linear") != std::string::npos)
                    {
                        ratio = initial_alpha_parameter * (1.0 - scalar_value/value_threshold);
                    }
                    else if (interpolation.find("Exponential") != std::string::npos)
                    {
                        ratio = - std::log(scalar_value/value_threshold) * initial_alpha_parameter + tolerance;
                        if (ratio > 1.0)
                        {
                            ratio = 1.0;
                        }
                    }
                    else
                    {
                        std::cout << "No interpolation defined, considering constant" << std:: endl;
                        ratio = initial_alpha_parameter;
                    }
                }
            }
            
            ComputeTensorH(scalar_value, gradient_value, ratio, element_size, value_threshold, i + 1);
        }
        
        /** 
         * 4) Check if the number of given entities match with mesh size 
         */
        
        CheckMeshData();
        
        // Save to file
        if (save_to_file == true)
        {
            SaveSolutionToFile(false);
        }
       
        //NOTE: Don't free memmory, it will be used in the Execute
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * We consider as input an already existing files , so this reads the input files 
     */
    
    void ReadFiles()
    {
       // NOTE: Step one in the constructor
        
       /** 
        * 2) Build mesh in MMG5 format 
        */
       
       SetInputMeshName();
       
       /** 
        * 3) Build sol in MMG5 format 
        */
       
       SetInputSolName();
       
       /** 
        * 4) Check if the number of given entities match with mesh size
        */
       
       CheckMeshData();
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * The operation related with the libray is executed
     * @param rThisModelPart: The model part to modify (first we empty, later we fill)
     * @param mColors: The external dictionary needed to know how to assign the submodelparts
     * @param save_to_file: Save the solution to a *.sol and *.mesh files
     */

    void Execute(
        ModelPart& rThisModelPart,
        const bool save_to_file = false,
        const unsigned int MaxNumberOfResults = 1000
        )
    {   
        // We initialize some values
        unsigned int step_data_size = rThisModelPart.GetNodalSolutionStepDataSize();
        unsigned int buffer_size    = rThisModelPart.NodesBegin()->GetBufferSize();
                
        std::cout << "Step data size: " << step_data_size << " Buffer size: " << buffer_size << std::endl; 
        
        ////////* MMG LIBRARY CALL *////////
        
        MMGLibCall();
        const unsigned int nNodes = mmgMesh->np;
        array_1d<unsigned int, TDim - 1> nConditions;
        if (TDim == 2)
        {
            nConditions[0] = mmgMesh->na;
        }
        else
        {
            nConditions[0] = mmgMesh->nt;
            nConditions[1] = mmgMesh->nquad;
        }
        array_1d<unsigned int, TDim - 1> nElements;
        if (TDim == 2)
        {
            nElements[0] = mmgMesh->nt;
        }
        else
        {
            nElements[0] = mmgMesh->ne;
            nElements[1] = mmgMesh->nprism;
        }
        
        std::cout << "     Nodes created: " << nNodes << std::endl;
        if (TDim == 2) // 2D
        {
            std::cout << "Conditions created: " << nConditions[0] << std::endl;
            std::cout << "Elements created: " << nElements[0] << std::endl;
        }
        else // 3D
        {
            std::cout << "Conditions created: " << nConditions[0] + nConditions[1] << std::endl;
            std::cout << "\tTriangles: " << nConditions[0] << "\tQuadrilaterals: " << nConditions[1] << std::endl;
            std::cout << "Elements created: " << nElements[0] + nElements[1] << std::endl;
            std::cout << "\tTetrahedron: " << nElements[0] << "\tPrisms: " << nElements[1] << std::endl;
        }
        
        ////////* EMPTY AND BACKUP THE MODEL PART *////////
        
        ModelPart OldModelPart;
        
        // First we empty the model part
        for (NodeConstantIterator node_iterator = rThisModelPart.NodesBegin(); node_iterator != rThisModelPart.NodesEnd(); node_iterator++)
        {
            node_iterator->Set(TO_ERASE, true);
            OldModelPart.AddNode(*(node_iterator.base()));
        }
        rThisModelPart.RemoveNodesFromAllLevels(TO_ERASE);  
        
        for (ConditionConstantIterator condition_iterator = rThisModelPart.ConditionsBegin(); condition_iterator != rThisModelPart.ConditionsEnd(); condition_iterator++)
        {
            condition_iterator->Set(TO_ERASE, true);
        }
        rThisModelPart.RemoveConditionsFromAllLevels(TO_ERASE); 
        
        for (ElementConstantIterator elem_iterator = rThisModelPart.ElementsBegin(); elem_iterator != rThisModelPart.ElementsEnd(); elem_iterator++)
        {
            elem_iterator->Set(TO_ERASE, true);
            OldModelPart.AddElement(*(elem_iterator.base()));
        }
        rThisModelPart.RemoveElementsFromAllLevels(TO_ERASE);  
        
        // NOTE: Technically not necessary
//         // Create iterator
//         typedef std::map<int,std::vector<std::string>>::iterator it_type;
//         
//         // Adding the submodelparts
//         for(it_type iterator = mColors.begin(); iterator != mColors.end(); iterator++) 
//         {
//             const int key = iterator->first;
//             const std::vector<std::string> value = iterator->second;
//             
//             if (key != 0)
//             {
//                 for (unsigned int it_name = 0; it_name < value.size(); it_name++)
//                 {
//                     if ((value[it_name].find(rThisModelPart.Name()) != std::string::npos) == false) // TODO: Check this!!!
//                     {
//                         if (rThisModelPart.HasSubModelPart(value[it_name]) == false)
//                         {
//                             rThisModelPart.CreateSubModelPart(value[it_name]);
//                         }
//                     }
//                 }
//             }
//         }
        
        // Create a new model part
        /* NODES */
        unsigned int node_id = 0;
        for (int unsigned i_node = 1; i_node <= nNodes; i_node++)
        {
            node_id += 1;
            NodeType::Pointer pNode;
            
            if (TDim == 2)
            {
                pNode = rThisModelPart.CreateNewNode(node_id, mmgMesh->point[i_node].c[0], mmgMesh->point[i_node].c[1], 0.0);
            }
            else
            {
                pNode = rThisModelPart.CreateNewNode(node_id, mmgMesh->point[i_node].c[0], mmgMesh->point[i_node].c[1], mmgMesh->point[i_node].c[2]);
            }
            
            // Set the DOFs in the nodes 
            for (typename Node<3>::DofsContainerType::const_iterator i_dof = mDofs.begin(); i_dof != mDofs.end(); i_dof++)
            {
                pNode->pAddDof(*i_dof);
            }
        }
        
        /* CONDITIONS */
        unsigned int cond_id = 0;
        if (mpRefCondition[0] != nullptr)
        {
            unsigned int cond_id = 0;

            for (int unsigned i_cond = 1; i_cond <= nConditions[0]; i_cond++)
            {
                cond_id += 1;
                ConditionType::Pointer pCondition;
                unsigned int prop_id;
                
                if (TDim == 2) // Lines
                {
                    const cond_geometries_2d index_geom = Line;
                    
                    std::vector<NodeType::Pointer> ConditionNodes (2);
                    ConditionNodes[0] = rThisModelPart.pGetNode(mmgMesh->edge[i_cond].a);
                    ConditionNodes[1] = rThisModelPart.pGetNode(mmgMesh->edge[i_cond].b);
                    
                    prop_id = mmgMesh->edge[i_cond].ref;
                    
                    pCondition = mpRefCondition[index_geom]->Create(cond_id, ConditionNodes, mpRefCondition[index_geom]->pGetProperties());
                }
                else // Triangles
                {
                    const cond_geometries_3d index_geom = Triangle3D;
                    
                    std::vector<NodeType::Pointer> ConditionNodes (3);
                    ConditionNodes[0] = rThisModelPart.pGetNode(mmgMesh->tria[i_cond].v[0]);
                    ConditionNodes[1] = rThisModelPart.pGetNode(mmgMesh->tria[i_cond].v[1]);
                    ConditionNodes[2] = rThisModelPart.pGetNode(mmgMesh->tria[i_cond].v[2]);
                    
                    prop_id = mmgMesh->tria[i_cond].ref;
                    
                    pCondition = mpRefCondition[index_geom]->Create(cond_id, ConditionNodes, mpRefCondition[index_geom]->pGetProperties());
                }
                
                pCondition->Initialize();
                rThisModelPart.AddCondition(pCondition);
                                    
                if (prop_id != 0) // NOTE: prop_id == 0 is the MainModelPart
                {
                    std::vector<std::string> ColorList = mColors[prop_id];
                    for (unsigned int colors = 0; colors < ColorList.size(); colors++)
                    {
                        std::string SubModelPartName = ColorList[colors];
                        ModelPart& SubModelPart = rThisModelPart.GetSubModelPart(SubModelPartName);
                        SubModelPart.AddCondition(pCondition);
                    }
                }
            }
        }
        if (TDim == 3)
        {
            if (mpRefCondition[1] != nullptr) // Quadrilateral
            {
                for (int unsigned i_cond = 1; i_cond <= nConditions[1]; i_cond++)
                {
                    cond_id += 1;
                    ConditionType::Pointer pCondition;
                    unsigned int prop_id;
                    
                    const cond_geometries_3d index_geom = Quadrilateral3D;
                    
                    std::vector<NodeType::Pointer> ConditionNodes (4);
                    ConditionNodes[0] = rThisModelPart.pGetNode(mmgMesh->quad[i_cond].v[0]);
                    ConditionNodes[1] = rThisModelPart.pGetNode(mmgMesh->quad[i_cond].v[1]);
                    ConditionNodes[2] = rThisModelPart.pGetNode(mmgMesh->quad[i_cond].v[2]);
                    ConditionNodes[3] = rThisModelPart.pGetNode(mmgMesh->quad[i_cond].v[3]);
                    
                    prop_id = mmgMesh->quad[i_cond].ref;
                    
                    pCondition = mpRefCondition[index_geom]->Create(cond_id, ConditionNodes, mpRefCondition[index_geom]->pGetProperties());
                    
                    pCondition->Initialize();
                    rThisModelPart.AddCondition(pCondition);
                                        
                    if (prop_id != 0) // NOTE: prop_id == 0 is the MainModelPart
                    {
                        std::vector<std::string> ColorList = mColors[prop_id];
                        for (unsigned int colors = 0; colors < ColorList.size(); colors++)
                        {
                            std::string SubModelPartName = ColorList[colors];
                            ModelPart& SubModelPart = rThisModelPart.GetSubModelPart(SubModelPartName);
                            SubModelPart.AddCondition(pCondition);
                        }
                    }
                }
            }
        }
        
        /* ELEMENTS */
        unsigned int elem_id = 0;
        if (mpRefElement[0] != nullptr)
        {
            for (int unsigned i_elem = 1; i_elem <= nElements[0]; i_elem++)
            {
                elem_id += 1;
                ElementType::Pointer pElement;
                unsigned int prop_id;
                
                if (TDim == 2) // Triangle
                {
                    const elem_geometries_2d index_geom = Triangle2D;
                    
                    std::vector<NodeType::Pointer> ElementNodes (3);
                    ElementNodes[0] = rThisModelPart.pGetNode(mmgMesh->tria[i_elem].v[0]);
                    ElementNodes[1] = rThisModelPart.pGetNode(mmgMesh->tria[i_elem].v[1]);
                    ElementNodes[2] = rThisModelPart.pGetNode(mmgMesh->tria[i_elem].v[2]);
                    
                    prop_id = mmgMesh->tria[i_elem].ref;
                    
                    pElement = mpRefElement[index_geom]->Create(elem_id, ElementNodes, mpRefElement[index_geom]->pGetProperties());
                }
                else // Tetrahedra
                {
                    const elem_geometries_3d index_geom = Tetrahedra;
                    
                    std::vector<NodeType::Pointer> ElementNodes (4);
                    ElementNodes[0] = rThisModelPart.pGetNode(mmgMesh->tetra[i_elem].v[0]);
                    ElementNodes[1] = rThisModelPart.pGetNode(mmgMesh->tetra[i_elem].v[1]);
                    ElementNodes[2] = rThisModelPart.pGetNode(mmgMesh->tetra[i_elem].v[2]);
                    ElementNodes[3] = rThisModelPart.pGetNode(mmgMesh->tetra[i_elem].v[3]);
                    
                    prop_id = mmgMesh->tetra[i_elem].ref;
                    
                    pElement = mpRefElement[index_geom]->Create(elem_id, ElementNodes, mpRefElement[index_geom]->pGetProperties());
                }
                
                pElement->Initialize();
                rThisModelPart.AddElement(pElement);
                
                if (prop_id != 0) // NOTE: prop_id == 0 is the MainModelPart
                {
                    std::vector<std::string> ColorList = mColors[prop_id];
                    for (unsigned int colors = 0; colors < ColorList.size(); colors++)
                    {
                        std::string SubModelPartName = ColorList[colors];
                        ModelPart& SubModelPart = rThisModelPart.GetSubModelPart(SubModelPartName);
                        SubModelPart.AddElement(pElement);
                    }
                }
            }
        }
        if (TDim == 3)
        {
            if (mpRefElement[1] != nullptr) // Prism
            {
                for (int unsigned i_elem = 1; i_elem <= nElements[1]; i_elem++)
                {
                    elem_id += 1;
                    ElementType::Pointer pElement;
                    unsigned int prop_id;
                    
                    const elem_geometries_3d index_geom = Prism;
                    
                    std::vector<NodeType::Pointer> ElementNodes (6);
                    ElementNodes[0] = rThisModelPart.pGetNode(mmgMesh->prism[i_elem].v[0]);
                    ElementNodes[1] = rThisModelPart.pGetNode(mmgMesh->prism[i_elem].v[1]);
                    ElementNodes[2] = rThisModelPart.pGetNode(mmgMesh->prism[i_elem].v[2]);
                    ElementNodes[3] = rThisModelPart.pGetNode(mmgMesh->prism[i_elem].v[3]);
                    ElementNodes[4] = rThisModelPart.pGetNode(mmgMesh->prism[i_elem].v[4]);
                    ElementNodes[5] = rThisModelPart.pGetNode(mmgMesh->prism[i_elem].v[5]);
                    
                    prop_id = mmgMesh->prism[i_elem].ref;
                    
                    pElement = mpRefElement[index_geom]->Create(elem_id, ElementNodes, mpRefElement[index_geom]->pGetProperties());
                    
                    pElement->Initialize();
                    rThisModelPart.AddElement(pElement);
                    
                    if (prop_id != 0) // NOTE: prop_id == 0 is the MainModelPart
                    {
                        std::vector<std::string> ColorList = mColors[prop_id];
                        for (unsigned int colors = 0; colors < ColorList.size(); colors++)
                        {
                            std::string SubModelPartName = ColorList[colors];
                            ModelPart& SubModelPart = rThisModelPart.GetSubModelPart(SubModelPartName);
                            SubModelPart.AddElement(pElement);
                        }
                    }
                }
            }
        }
        
        //  Get the list of submodelparts names
        const std::vector<std::string> SubModelPartNames = rThisModelPart.GetSubModelPartNames();
       
        // Add the nodes to the differents submodelparts
        for (unsigned int i_model_part = 0; i_model_part < rThisModelPart.NumberOfSubModelParts(); i_model_part++)
        {
            ModelPart& rSubModelPart = rThisModelPart.GetSubModelPart(SubModelPartNames[i_model_part]);
           
            std::set<int> aux_set;
           
            for (ElementConstantIterator elem_iterator = rSubModelPart.ElementsBegin(); elem_iterator != rSubModelPart.ElementsEnd(); elem_iterator++)
            {
                for (unsigned int i_node = 0; i_node < elem_iterator->GetGeometry().size(); i_node++)
                {
                    aux_set.insert(elem_iterator->GetGeometry()[i_node].Id());
                }
            }
           
            for (ConditionConstantIterator condition_iterator = rSubModelPart.ConditionsBegin(); condition_iterator != rSubModelPart.ConditionsEnd(); condition_iterator++)
            {
                for (unsigned int i_node = 0; i_node < condition_iterator->GetGeometry().size(); i_node++)
                {
                    aux_set.insert(condition_iterator->GetGeometry()[i_node].Id());
                }
            }
           
            // Clean duplicated nodes
            std::vector<IndexType> NodesIds;
            for( auto it = aux_set.begin(); it != aux_set.end(); ++it ) 
            {
                NodesIds.push_back(*it);
            }
           
            rSubModelPart.AddNodes(NodesIds);
        }
        
        /* We interpolate all the values */
        
        // We create the locator
        BinBasedFastPointLocator<TDim> PointLocator = BinBasedFastPointLocator<TDim>(OldModelPart);
        PointLocator.UpdateSearchDatabase();
        
        // Iterate in the nodes
        NodesArrayType& pNode = rThisModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        /* Nodes */
//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            Vector shape_functions;
            ElementType::Pointer pElement;
            
            const bool found = PointLocator.FindPointOnMeshSimplified(itNode->Coordinates(), shape_functions, pElement, MaxNumberOfResults);
            
            if (found == false)
            {
                std::cout << "WARNING: Node "<< itNode->Id() << " not found (interpolation not posible)" << std::endl;
                std::cout << "\t X:"<< itNode->X() << "\t Y:"<< itNode->Y() << "\t Z:"<< itNode->Z() << std::endl;
            }
            else
            {
                for(unsigned int step = 0; step < buffer_size; step++)
                {
                    double* step_data = (itNode)->SolutionStepData().Data(step);
                    
                    if (TDim == 2)
                    {
                        double* node0_data = pElement->GetGeometry()[0].SolutionStepData().Data(step);
                        double* node1_data = pElement->GetGeometry()[1].SolutionStepData().Data(step);
                        double* node2_data = pElement->GetGeometry()[2].SolutionStepData().Data(step);
                        
                        for (unsigned int j = 0; j < step_data_size; j++)
                        {
                            step_data[j] = shape_functions[0] * node0_data[j]
                                         + shape_functions[1] * node1_data[j]
                                         + shape_functions[2] * node2_data[j];
                        }
                    }
                    else // NOTE: This just works with tetrahedron (you are going to have poblems with anything else)
                    {
                        double* node0_data = pElement->GetGeometry()[0].SolutionStepData().Data(step);
                        double* node1_data = pElement->GetGeometry()[1].SolutionStepData().Data(step);
                        double* node2_data = pElement->GetGeometry()[2].SolutionStepData().Data(step);
                        double* node3_data = pElement->GetGeometry()[3].SolutionStepData().Data(step);
                        
                        for (unsigned int j = 0; j < step_data_size; j++)
                        {
                            step_data[j] = shape_functions[0] * node0_data[j]
                                         + shape_functions[1] * node1_data[j]
                                         + shape_functions[2] * node2_data[j]
                                         + shape_functions[3] * node3_data[j];
                        }
                    }
                }
            }
        }
        
        /* Save to file */
        if (save_to_file == true)
        {
            SaveSolutionToFile(true);
        }
       
        /* Free memory */
        FreeMemory();
        
        /* We print the resulting model part */
        std::cout << "//---------------------------------------------------//" << std::endl;
        std::cout << "//---------------------------------------------------//" << std::endl;
        std::cout << "//---------------   AFTER REMESHING   ---------------//" << std::endl;
        std::cout << "//---------------------------------------------------//" << std::endl;
        std::cout << "//---------------------------------------------------//" << std::endl;
        std::cout << std::endl;
        
        KRATOS_WATCH(rThisModelPart);
    }
    
    /** ------------------------------ STEP III -------------------------- */
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * It saves the solution and mesh to files (for debugging pourpose g.e)
     * @param post_output: If the file to save is after or before remeshing
     */
    
    void SaveSolutionToFile(const bool post_output)
    {
        /* GET RESULTS */

        /** 
         *1) Automatically save the mesh 
         */
        
        OutputMesh(post_output);

        /** 
         * 2) Automatically save the solution 
         */
        
        OutputSol(post_output);
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * It frees the memory used during all the process
     */
    
    void FreeMemory()
    {
        /** 
         * 3) Free the MMG3D5 structures 
         */
        
        FreeAll();

        free(mFilename);
        mFilename = NULL;
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
    
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{
    
    // Storage for the dof of the node
    Node<3>::DofsContainerType  mDofs;
    
    // I/O information
    char* mFilename;
    std::string mStdStringFilename;
    
    // The member variables related with the MMG library
    MMG5_pMesh mmgMesh;
    MMG5_pSol  mmgSol;
    
    // Where the sub model parts IDs are stored
    std::map<int,std::vector<std::string>> mColors;
    
    // Reference element and condition
    std::vector<Element::Pointer>   mpRefElement;
    std::vector<Condition::Pointer> mpRefCondition;
    std::vector<bool> mInitRefElement;
    std::vector<bool> mInitRefCondition;
    
    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{
    
    /**
     * This inits the mesh
     */
    
    void InitMesh()
    {
       /** 1) Initialisation of mesh and sol structures 
        * args of InitMesh:
        * MMG5_ARG_start: we start to give the args of a variadic func
        * MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
        * &mmgMesh: pointer toward your MMG5_pMesh (that store your mesh)
        * MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
        * &mmgSol: pointer toward your MMG5_pSol (that store your metric) 
        */
       
        if (TDim == 2)
        {
            MMG2D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end); 
        }
        else
        {
            MMG3D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end); 
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This sets the size of the mesh
     * @param numNodes: Number of nodes
     * @param numElements: Number of Elements
     * @param numConditions: Number of Conditions
     */
    
    void SetMeshSize(
        const int numNodes,
        const array_1d<int, TDim - 1> numArrayElements,  // NOTE: We do this tricky thing to take into account the prisms
        const array_1d<int, TDim - 1> numArrayConditions // NOTE: We do this tricky thing to take into account the quadrilaterals
        )
    {
        if (TDim == 2)
        {
            //Give the size of the mesh: numNodes vertices, numElements triangles, numConditions edges (2D) 
            if ( MMG2D_Set_meshSize(mmgMesh, numNodes, numArrayElements[0], numArrayConditions[0]) != 1 ) 
            {
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            //Give the size of the mesh: numNodes vertex, numElements tetra and prism, numArrayConditions triangles and quadrilaterals, 0 edges (3D) 
            if ( MMG3D_Set_meshSize(mmgMesh, numNodes, numArrayElements[0], numArrayElements[1], numArrayConditions[0], numArrayConditions[1], 0) != 1 ) 
            {
                exit(EXIT_FAILURE);
            }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This sets the size of the solution
     * @param numNodes: Number of nodes
     */
    
    void SetSolSize(const int numNodes)
    {
        if (TDim == 2)
        {
            if ( MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,numNodes,MMG5_Tensor) != 1 )
            {
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            if ( MMG3D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,numNodes,MMG5_Tensor) != 1 )
            {
                exit(EXIT_FAILURE);
            }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This checks the mesh data
     */
    
    void CheckMeshData()
    {
        if (TDim == 2)
        {
            if ( MMG2D_Chk_meshData(mmgMesh, mmgSol) != 1 ) 
            {
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            if ( MMG3D_Chk_meshData(mmgMesh, mmgSol) != 1 ) 
            {
                exit(EXIT_FAILURE);
            }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This set the input mesh name
     */
    
    void SetInputMeshName()
    {
        if (TDim == 2)
        {
            if ( MMG2D_Set_inputMeshName(mmgMesh, mFilename) != 1 )
            {
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            if ( MMG3D_Set_inputMeshName(mmgMesh, mFilename) != 1 )
            {
                exit(EXIT_FAILURE);
            }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This sets the output mesh
     */
    
    void OutputMesh(const bool post_output)
    {
        std::string MeshName;
        if (post_output == true)
        {
            MeshName = mStdStringFilename+".o.mesh";
        }
        else
        {
            MeshName = mStdStringFilename+".mesh";
        }
        
        char* MeshFile = new char [MeshName.length() + 1];
        std::strcpy (MeshFile, MeshName.c_str());
        
        if (TDim == 2)
        { 
            // a)  Give the ouptut mesh name using MMG2D_Set_outputMeshName (by default, the mesh is saved in the "mesh.o.mesh" file  
            MMG2D_Set_outputMeshName(mmgMesh,MeshFile);

            // b) function calling 
            if ( MMG2D_saveMesh(mmgMesh,MeshFile) != 1 ) 
            {
                std::cout << "UNABLE TO SAVE MESH" << std::endl;
            }
        }
        else
        {
            // a)  Give the ouptut mesh name using MMG3D_Set_outputMeshName (by default, the mesh is saved in the "mesh.o.mesh" file 
            MMG3D_Set_outputMeshName(mmgMesh,MeshFile);

            // b) function calling 
            if ( MMG3D_saveMesh(mmgMesh,MeshFile) != 1 ) 
            {
                std::cout << "UNABLE TO SAVE MESH" << std::endl;
            }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This sets the input sol name
     */
    
    void SetInputSolName()
    {
        if (TDim == 2)
        {
            // a) Give the sol name (by default, the "mesh.sol" file is oppened)
            if ( MMG2D_Set_inputSolName(mmgMesh, mmgSol, mFilename) != 1 )
            {
                exit(EXIT_FAILURE);
            }

            // b) Function calling 
            if ( MMG2D_loadSol(mmgMesh, mmgSol, mFilename) != 1 )
            {
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            // a) Give the sol name (by default, the "mesh.sol" file is oppened)
            if ( MMG3D_Set_inputSolName(mmgMesh, mmgSol, mFilename) != 1 )
            {
                exit(EXIT_FAILURE);
            }

            // b) Function calling 
            if ( MMG3D_loadSol(mmgMesh, mmgSol, mFilename) != 1 )
            {
                exit(EXIT_FAILURE);
            }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This sets the output sol
     */
    
    void OutputSol(const bool post_output)
    {
        std::string SolName;
        if (post_output == true)
        {
            SolName = mStdStringFilename+".o.sol";
        }
        else
        {
            SolName = mStdStringFilename+".sol";
        }
        
        char* SolFile = new char [SolName.length() + 1];
        std::strcpy (SolFile, SolName.c_str());
        
        if (TDim == 2)
        { 
            // a)  Give the ouptut sol name using MMG2D_Set_outputSolName (by default, the mesh is saved in the "mesh.o.sol" file 
            MMG2D_Set_outputSolName(mmgMesh, mmgSol, SolFile);

            // b) Function calling 
            if ( MMG2D_saveSol(mmgMesh, mmgSol, SolFile) != 1 ) 
            {
                std::cout << "UNABLE TO SAVE SOL" << std::endl;
            }
        }
        else
        {
            // a)  Give the ouptut sol name using MMG3D_Set_outputSolName (by default, the mesh is saved in the "mesh.o.sol" file 
            MMG3D_Set_outputSolName(mmgMesh, mmgSol, SolFile);

            // b) Function calling 
            if ( MMG3D_saveSol(mmgMesh,mmgSol, SolFile) != 1 ) 
            {
                std::cout << "UNABLE TO SAVE SOL" << std::endl;
            }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This loads the solution
     */
    
    void MMGLibCall()
    {
        if (TDim == 2)
        {
            const int ier = MMG2D_mmg2dlib(mmgMesh, mmgSol);

            if ( ier == MMG5_STRONGFAILURE ) 
            {
                std::cout << "BAD ENDING OF MMG2DLIB: UNABLE TO SAVE MESH" << std::endl;
            }
            else if ( ier == MMG5_LOWFAILURE )
            {
                std::cout << "BAD ENDING OF MMG2DLIB" << std::endl;
            }
        }
        else
        {
            const int ier = MMG3D_mmg3dlib(mmgMesh, mmgSol);

            if ( ier == MMG5_STRONGFAILURE ) 
            {
                std::cout << "BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH" << std::endl;
            }
            else if ( ier == MMG5_LOWFAILURE )
            {
                std::cout << "BAD ENDING OF MMG3DLIB" << std::endl;
            }
        }

    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This frees the MMG structures
     */
    
    void FreeAll()
    {
        if (TDim == 2)
        {
            MMG2D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,MMG5_ARG_end);
        }
        else
        {
            MMG3D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,MMG5_ARG_end);
        }

    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This sets the nodes of the mesh
     * @param X: Coordinate X
     * @param Y: Coordinate Y
     * @param Z: Coordinate Z
     * @param color: Reference of the node(submodelpart)
     * @param index: The index number of the node 
     */
    
    void SetNodes(
        const double X,
        const double Y,
        const double Z,
        const int color,
        const int index
        )
    {
        // Using API
        if (TDim == 2)
        {
            if ( MMG2D_Set_vertex(mmgMesh, X, Y, color, index) != 1 )  
            {
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            if ( MMG3D_Set_vertex(mmgMesh, X, Y, Z, color, index) != 1 )  
            {
                exit(EXIT_FAILURE); 
            }
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This sets the conditions of the mesh
     * @param Geom: The geometry of the condition
     * @param color: Reference of the node(submodelpart)
     * @param index: The index number of the node 
     */
    
    void SetConditions(
        const Geometry<Node<3> > & Geom,
        const int color,
        const int index
        )
    {
        const int id1 = Geom[0].Id(); // First node id
        const int id2 = Geom[1].Id(); // Second node id
        
        // Using API
        if (TDim == 2)
        {
            if ( MMG2D_Set_edge(mmgMesh, id1, id2, color, index) != 1 ) 
            {
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            const int id3 = Geom[2].Id(); // Third node id
            
            const unsigned int size_geom = Geom.size();
            
            if (size_geom == 3) // Triangle
            {
                if ( MMG3D_Set_triangle(mmgMesh, id1, id2, id3, color, index) != 1 )  
                {
                    exit(EXIT_FAILURE); 
                }
            }
            else if (size_geom == 4) // Quadrilaterals
            {
                const int id4 = Geom[3].Id(); // Fourth node id
                
                if ( MMG3D_Set_quadrilateral(mmgMesh, id1, id2, id3, id4, color, index) != 1 )  
                {
                    exit(EXIT_FAILURE); 
                }
            }
            else
            {
                KRATOS_THROW_ERROR( std::logic_error, "WARNING: I DO NOT KNOW WHAT IS THIS. Size: ", size_geom );
            }
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This sets elements of the mesh
     * @param Geom: The geometry of the element
     * @param color: Reference of the node(submodelpart)
     * @param index: The index number of the node 
     */
    
    void SetElements(
        const Geometry<Node<3> > & Geom,
        const int color,
        const int index
        )
    {
        const int id1 = Geom[0].Id(); // First node id
        const int id2 = Geom[1].Id(); // Second node id
        const int id3 = Geom[2].Id(); // Third node id
        
        // Using API
        if (TDim == 2)
        {
            if ( MMG2D_Set_triangle(mmgMesh, id1, id2, id3, color, index) != 1 ) 
            {
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            const int id4 = Geom[3].Id(); // Fourth node id
            
            const unsigned int size_geom = Geom.size();
            
            if (size_geom == 4) // Tetrahedron
            {
                if ( MMG3D_Set_tetrahedron(mmgMesh, id1, id2, id3, id4, color, index) != 1 )  
                {
                    exit(EXIT_FAILURE); 
                }
            }
            else if (size_geom == 6) // Prisms
            {
                const int id5 = Geom[4].Id(); // 5th node id
                const int id6 = Geom[5].Id(); // 6th node id
                
                if ( MMG3D_Set_prism(mmgMesh, id1, id2, id3, id4, id5, id6, color, index) != 1 )  
                {
                    exit(EXIT_FAILURE); 
                }
            }
            else if (size_geom == 8) // Hexaedron
            {
//                 const int id5 = Geom[4].Id(); // 5th node id
//                 const int id6 = Geom[5].Id(); // 6th node id
//                 const int id6 = Geom[7].Id(); // 7th node id
//                 const int id6 = Geom[8].Id(); // 8th node id
                
                KRATOS_THROW_ERROR( std::logic_error, "WARNING: HEXAEDRON NON IMPLEMENTED IN THE LIBRARY. Size: ", size_geom );
            }
            else
            {
                KRATOS_THROW_ERROR( std::logic_error, "WARNING: I DO NOT KNOW WHAT IS THIS. Size: ", size_geom );
            }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This functions gets the "colors", parts of a model part to process
     * @param rThisModelPart: The model part to get the colors
     */
    
    void ComputeColors(
        ModelPart& rThisModelPart,
        std::map<int,int>& node_colors,
        std::map<int,int>& cond_colors,
        std::map<int,int>& elem_colors
        )
    {        
        // Initialize and create the auxiliar maps
        const std::vector<std::string> SubModelPartNames = rThisModelPart.GetSubModelPartNames();
        std::map<int,std::set<int>> aux_node_colors, aux_cond_colors, aux_elem_colors;
        
        // TODO: Add when subsubmodelparts work correctly
//         const unsigned int initsub = SubModelPartNames.size(); 
//         for (unsigned int i_sub = 0; i_sub < initsub; i_sub++)
//         {
//             ModelPart& rSubModelPart = rThisModelPart.GetSubModelPart(SubModelPartNames[i_sub]);
            
//             const std::vector<std::string> SubSubModelPartNames = rSubModelPart.GetSubModelPartNames();
//             
//             for (unsigned int i_sub_sub = 0; i_sub_sub < SubSubModelPartNames.size(); i_sub_sub++)
//             {
//                 SubModelPartNames.push_back(SubSubModelPartNames[i_sub_sub]);
//             }
//         }
        
        std::vector<std::string> ModelPartNames;
        ModelPartNames.push_back(rThisModelPart.Name());
        for (unsigned int i_sub = 0; i_sub < SubModelPartNames.size(); i_sub++)
        {
            ModelPartNames.push_back(SubModelPartNames[i_sub]);
        }
        
        // Initialize colors
        int color = 0;
        for (unsigned int i_sub = 0; i_sub < ModelPartNames.size(); i_sub++)
        {
            mColors[i_sub].push_back(ModelPartNames[i_sub]);
            
            if (color > 0)
            {
                ModelPart& rSubModelPart = rThisModelPart.GetSubModelPart(ModelPartNames[i_sub]);
                
                // Iterate in the nodes
                NodesArrayType& pNode = rSubModelPart.Nodes();
                auto numNodes = pNode.end() - pNode.begin();
                
                // Iterate in the conditions
                ConditionsArrayType& pConditions = rSubModelPart.Conditions();
                auto numConditions = pConditions.end() - pConditions.begin();
                
                // Iterate in the elements
                ElementsArrayType& pElements = rSubModelPart.Elements();
                auto numElements = pElements.end() - pElements.begin();
                
                /* Nodes */
        //         #pragma omp parallel for 
                for(unsigned int i = 0; i < numNodes; i++) 
                {
                    auto itNode = pNode.begin() + i;
                    aux_node_colors[itNode->Id()].insert(color);
                }
                
                /* Conditions */
        //         #pragma omp parallel for 
                for(unsigned int i = 0; i < numConditions; i++) 
                {
                    auto itCond = pConditions.begin() + i;
                    aux_cond_colors[itCond->Id()].insert(color);
                }
                
                /* Elements */
        //         #pragma omp parallel for 
                for(unsigned int i = 0; i < numElements; i++) 
                {
                    auto itElem = pElements.begin() + i;
                    aux_elem_colors[itElem->Id()].insert(color);
                }
            }
            
            color += 1;
        }
        
        // The iterator for the auxiliar maps is created
        typedef std::map<int,std::set<int>>::iterator it_type;
        
        // Now detect all the cases in which a node or a cond belongs to more than one part simultaneously 
        std::map<std::set<int>, int> combinations;
        
        /* Nodes */
        for(it_type iterator = aux_node_colors.begin(); iterator != aux_node_colors.end(); iterator++) 
        {
//             const int key = iterator->first;
            const std::set<int> value = iterator->second;
            
            if (value.size() > 1)
            {
                combinations[value] = -1;
            }
        }
        
        /* Conditions */
        for(it_type iterator = aux_cond_colors.begin(); iterator != aux_cond_colors.end(); iterator++) 
        {
//             const int key = iterator->first;
            const std::set<int> value = iterator->second;
            
            if (value.size() > 1)
            {
                combinations[value] = -1;
            }
        }

        /* Elements */
        for(it_type iterator = aux_elem_colors.begin(); iterator != aux_elem_colors.end(); iterator++) 
        {
//             const int key = iterator->first;
            const std::set<int> value = iterator->second;
            
            if (value.size() > 1)
            {
                combinations[value] = -1;
            }
        }
        
        /* Combinations */
        typedef std::map<std::set<int>,int>::iterator comb_type;
        for(comb_type iterator = combinations.begin(); iterator != combinations.end(); iterator++) 
        {
            const std::set<int> key = iterator->first;
//             const int value = iterator->second;
            
            for( auto it = key.begin(); it != key.end(); ++it ) 
            {
                mColors[color].push_back(mColors[*it][0]);
            }
            combinations[key] = color;
            color += 1;
            
        }
        
        // The final maps are created
        /* Nodes */
        for(it_type iterator = aux_node_colors.begin(); iterator != aux_node_colors.end(); iterator++) 
        {
            const int key = iterator->first;
            const std::set<int> value = iterator->second;
            
            if (value.size() == 0)
            {
                node_colors[key] = 0; // Main Model Part
            }
            else if (value.size() == 1) // Another Model Part
            {
                node_colors[key] = *value.begin();
            }
            else // There is a combination
            {
                node_colors[key] = combinations[value];
            }
        }
        
        /* Conditions */
        for(it_type iterator = aux_cond_colors.begin(); iterator != aux_cond_colors.end(); iterator++) 
        {
            const int key = iterator->first;
            const std::set<int> value = iterator->second;
            
            if (value.size() == 0)
            {
                cond_colors[key] = 0; // Main Model Part
            }
            else if (value.size() == 1) // Another Model Part
            {
                cond_colors[key] = *value.begin();
            }
            else // There is a combination
            {
                cond_colors[key] = combinations[value];
            }
        }
        
        /* Elements */
        for(it_type iterator = aux_elem_colors.begin(); iterator != aux_elem_colors.end(); iterator++) 
        {
            const int key = iterator->first;
            const std::set<int> value = iterator->second;
            
            if (value.size() == 0)
            {
                elem_colors[key] = 0; // Main Model Part
            }
            else if (value.size() == 1) // Another Model Part
            {
                elem_colors[key] = *value.begin();
            }
            else // There is a combination
            {
                elem_colors[key] = combinations[value];
            }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * It calculates the tensor of the scalar, necessary to get the solution before remeshing
     * @param scalar_value: The scalar used as reference to remesh
     * @param gradient_value: The gradient of the scalar to remesh
     * @param ratio: The alpha parameter used to remesh
     * @param element_size: The minimum size of the elements
     * @param value_threshold: The minimum value to consider for remesh
     */
        
    void ComputeTensorH(
        const double& scalar_value, 
        const array_1d<double, 3>& gradient_value,
        const double& ratio,
        const double& element_size,
        const double& value_threshold,
        const int node_id // NOTE: This can be a problem if the nodes are not correctly defined
    )
    {
        const double coeff0 = 1.0/(element_size * element_size);
        
        const double threshold = value_threshold; // We don't reach the minimum
//         const double threshold = element_size; // NOTE: In the case that the scalar value is related with a distance, better to consider an independent variable
        
        if (TDim == 2) // 2D: The order of the metric is m11,m12,m22
        {
            if (scalar_value > threshold) 
//             if (scalar_value > element_size) 
            {
                // Using API
                if ( MMG2D_Set_tensorSol(mmgSol, coeff0, 0.0, coeff0, node_id) != 1 )
                {
                    exit(EXIT_FAILURE);
                }
            }
            else
            {
                const double aux = ratio + (scalar_value/value_threshold)*(1.0 - ratio);
//                 const double aux = ratio + (scalar_value/element_size)*(1.0 - ratio);
                const double coeff1 = coeff0/(aux * aux);
                
                const double v0v0 = gradient_value[0]*gradient_value[0];
                const double v0v1 = gradient_value[0]*gradient_value[1];
                const double v1v1 = gradient_value[1]*gradient_value[1];
                
                // Using API
                if ( MMG2D_Set_tensorSol(mmgSol, 
                                        coeff0*(1.0 - v0v0) + coeff1*v0v0, 
                                        coeff0*(      v0v1) + coeff1*v0v1,  
                                        coeff0*(1.0 - v1v1) + coeff1*v1v1,
                                        node_id) != 1 )
                {
                    exit(EXIT_FAILURE);
                }
            }
        }
        else // 3D: The order of the metric is m11,m12,m13,m22,m23,m33
        {
            if (scalar_value > threshold)
            {
                // Using API
                if ( MMG3D_Set_tensorSol(mmgSol, coeff0, 0.0, 0.0, coeff0, 0.0, coeff0, node_id) != 1 )
                {
                    exit(EXIT_FAILURE);
                }
            }
            else
            {
                const double aux = ratio + (scalar_value/element_size)*(1.0 - ratio);
                const double coeff1 = coeff0/(aux * aux);
                
                const double v0v0 = gradient_value[0]*gradient_value[0];
                const double v0v1 = gradient_value[0]*gradient_value[1];
                const double v0v2 = gradient_value[0]*gradient_value[2];
                const double v1v1 = gradient_value[1]*gradient_value[1];
                const double v1v2 = gradient_value[1]*gradient_value[2];
                const double v2v2 = gradient_value[2]*gradient_value[2];
                
                // Using API
                if ( MMG3D_Set_tensorSol(mmgSol, 
                                        coeff0*(1.0 - v0v0) + coeff1*v0v0, 
                                        coeff0*(      v0v1) + coeff1*v0v1, 
                                        coeff0*(      v0v2) + coeff1*v0v2, 
                                        coeff0*(1.0 - v1v1) + coeff1*v1v1, 
                                        coeff0*(      v1v2) + coeff1*v1v2, 
                                        coeff0*(1.0 - v2v2) + coeff1*v2v2, 
                                        node_id) != 1 )
                {
                    exit(EXIT_FAILURE);
                }
            }
        }
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
};// class MmgUtility
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}// namespace Kratos.
#endif /* KRATOS_MMG_UTILITY defined */
