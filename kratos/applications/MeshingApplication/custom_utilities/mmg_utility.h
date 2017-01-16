// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____ 
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _ 
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                                       license: structural_mechanics_application/license.txt
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
#include "mmg/mmg2d/libmmg2d.h" // TODO: Version in 2D
#include "mmg/mmg3d/libmmg3d.h"
#include "mmg/mmgs/libmmgs.h" // TODO: what is this library doing?

// NOTE: Inspired in the documentation from:
// https://github.com/MmgTools/mmg/blob/master/libexamples/

// NOTE: The following contains the license of the MMG library //TODO: Ask Riccardo about this
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

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

// template< unsigned int TDim>  // TODO: Use a template to define the 3D or 2D case
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
    typedef MeshType::ElementConstantIterator           ElementConstantIterator;
    typedef MeshType::ConditionConstantIterator       ConditionConstantIterator;
    
    ///@}
    ///@name Life Cycle
    ///@{
     
    /** ------------------------------ STEP   I -------------------------- */
    // Constructor
    
    /**
     * This is the default constructor, which is used to read the input files 
     * @param Filename: The input name of the 
     */
    
    MmgUtility(const std::string Filename)
        :mStdStringFilename(Filename)
    {
       /** 1) Initialisation of mesh and sol structures */
       /* args of InitMesh:
        * MMG5_ARG_start: we start to give the args of a variadic func
        * MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
        * &mmgMesh: pointer toward your MMG5_pMesh (that store your mesh)
        * MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
        * &mmgSol: pointer toward your MMG5_pSol (that store your metric) 
        */
       
       mFilename = new char [Filename.length() + 1];
       std::strcpy (mFilename, Filename.c_str());
       
       mmgMesh = NULL;
       mmgSol = NULL;
       
       MMG3D_Init_mesh(
                  MMG5_ARG_start,
                  MMG5_ARG_ppMesh,
                  &mmgMesh,
                  MMG5_ARG_ppMet,
                  &mmgSol,
                  MMG5_ARG_end
                  );
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
    
    int ComputeExistingModelPart(
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
        
        // First we compute the colors
        std::map<int,int> node_colors, cond_colors, elem_colors;
        ComputeColors(rThisModelPart, node_colors, cond_colors, elem_colors);
        
        ////////* MESH FILE *////////
                
        /** 2) Build mesh in MMG5 format */
        /** Two solutions: just use the MMG3D_loadMesh function that will read a .mesh(b)
        file formatted or manually set your mesh using the MMG3D_Set* functions */
        
        // Iterate in the nodes
        NodesArrayType& pNode = rThisModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        // Iterate in the conditions
        ConditionsArrayType& pConditions = rThisModelPart.Conditions();
        auto numConditions = pConditions.end() - pConditions.begin();
        
        // Iterate in the elements
        ElementsArrayType& pElements = rThisModelPart.Elements();
        auto numElements = pElements.end() - pElements.begin();
        
        /** Manually set of the mesh */
        /** a) give the size of the mesh: numNodes vertices, numElements tetra, numConditions triangles, 0 edges (3D) */
        if ( MMG3D_Set_meshSize(mmgMesh,numNodes,numElements,numConditions,0) != 1 ) 
        {
            exit(EXIT_FAILURE);
        }
        
        /* Nodes */
//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            /** b) give the vertices: for each vertex, give the coordinates, the reference
             *and the position in mesh of the vertex */
            
            mmgMesh->point[i + 1].c[0] = itNode->X();  
            mmgMesh->point[i + 1].c[1] = itNode->Y(); 
            mmgMesh->point[i + 1].c[2] = itNode->Z(); 
            mmgMesh->point[i + 1].ref  = node_colors[itNode->Id()];
            
            /* or with the api function :*/
//             if ( MMG3D_Set_vertex(mmgMesh, itNode.X(), itNode.Y(), itNode.Z(), node_colors[itNode.Id()], i  + 1) != 1 )  
//             {
//                 exit(EXIT_FAILURE); 
//             }
        }
        
        /* Conditions */
//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCond = pConditions.begin() + i;
           /* tria*/
           mmgMesh->tria[i + 1].v[0] = itCond->GetGeometry()[0].Id();  
           mmgMesh->tria[i + 1].v[1] = itCond->GetGeometry()[1].Id();  
           mmgMesh->tria[i + 1].v[2] = itCond->GetGeometry()[2].Id();
           mmgMesh->tria[i + 1].ref  = cond_colors[itCond->Id()];
           
           /* or with the api function :*/
//            if ( MMG3D_Set_triangle(mmgMesh, itElem.GetGeometry()[0].Id() ,itElem.GetGeometry()[1].Id() ,itElem.GetGeometry()[2].Id(), cond_colors[itCond.Id()], i + 1) != 1 )  
//            {
//                exit(EXIT_FAILURE); 
//            }
        }
        
        /* Elements */
//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numElements; i++) 
        {
            auto itElem = pElements.begin() + i;
           /* tetra*/
           mmgMesh->tetra[i + 1].v[0] = itElem->GetGeometry()[0].Id();  
           mmgMesh->tetra[i + 1].v[1] = itElem->GetGeometry()[1].Id();  
           mmgMesh->tetra[i + 1].v[2] = itElem->GetGeometry()[2].Id();
           mmgMesh->tetra[i + 1].v[3] = itElem->GetGeometry()[3].Id();
           mmgMesh->tetra[i + 1].ref  = elem_colors[itElem->Id()];
           
           /* or with the api function :*/
//            if ( MMG3D_Set_tetrahedra(mmgMesh, itElem.GetGeometry()[0].Id() ,itElem.GetGeometry()[1].Id() ,itElem.GetGeometry()[2].Id() ,itElem.GetGeometry()[3].Id(), elem_colors[itCond.Id()], i + 1) != 1 )  
//            {
//                exit(EXIT_FAILURE); 
//            }
        }
        
        ////////* SOLUTION FILE *////////
        
        if ( MMG3D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,numNodes,MMG5_Tensor) != 1 )
        {
            exit(EXIT_FAILURE);
        }

//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            const double scalar_value = itNode->FastGetSolutionStepValue(rVariable);
            const array_1d<double, 3> gradient_value = itNode->FastGetSolutionStepValue(rVariableGradient);
            double element_size = itNode->GetValue(NODAL_VOLUME);
            if (element_size > elementary_length)
            {
                element_size = elementary_length;
            }
        
            double ratio = 1.0; // NOTE: Isotropic mesh
            if (scalar_value < value_threshold)
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
                    ratio = - std::log(scalar_value/value_threshold) * initial_alpha_parameter + 1.0e-12;
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
            
            ComputeTensorH(scalar_value, gradient_value, ratio, element_size, value_threshold, i + 1);
        }
        
        /** 4) If you don't use the API functions, you MUST call
        * the MMG3D_Set_handGivenMesh() function. Don't call it if you use
        * the API functions */
        MMG3D_Set_handGivenMesh(mmgMesh);
        
        /** 5) (not mandatory): check if the number of given entities match with mesh size */
        if ( MMG3D_Chk_meshData(mmgMesh,mmgSol) != 1 ) 
        {
            exit(EXIT_FAILURE);
        }
        
       // Save to file
       int ier = 0;
       if (save_to_file == true)
       {
           ier = SaveSolutionToFile(false);
       }
       
       //NOTE: Don't free memmory, it will be used in the Execute
       
       return ier;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * We consider as input an already existing files , so this reads the input files 
     */
    
    void ReadFiles()
    {
       // NOTE: Step one in the constructor
        
       /** 2) Build mesh in MMG5 format */
       if ( MMG3D_Set_inputMeshName(mmgMesh,mFilename) != 1 )
       {
           exit(EXIT_FAILURE);
       }
       
       /** 3) Build sol in MMG5 format */
       /* With MMG3D_loadSol function */
       /* a) (not mandatory): give the sol name
        * (by default, the "mesh.sol" file is oppened)
        */
       if ( MMG3D_Set_inputSolName(mmgMesh,mmgSol,mFilename) != 1 )
       {
           exit(EXIT_FAILURE);
       }
       
       /** b) function calling */
       if ( MMG3D_loadSol(mmgMesh,mmgSol,mFilename) != 1 )
       {
           exit(EXIT_FAILURE);
       }
       
       /** 4) (not mandatory): check if the number of given entities match with mesh size */if ( MMG3D_Chk_meshData(mmgMesh,mmgSol) != 1 ) 
       {
           exit(EXIT_FAILURE);
       }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * The operation related with the libray is executed
     * @param rThisModelPart: The original model part, an auxiliary model part will be created an the input model part will be modified
     * @param mColors: The external dictionary needed to know how to assign the submodelparts
     * @param save_to_file: Save the solution to a *.sol and *.mesh files
     */
    
    // FIXME: Core here
    int Execute( // NOTE: Properties are going to give problems (we need to define new properies copying the old ones)
        ModelPart& rThisModelPart, //NOTE: The model part shoul be empty
        Properties::Pointer pProperties, // TODO: Solve the problem, just one property
//         std::vector<Properties::Pointer> pPropertiesVector,
        const bool save_to_file = false
        )
    {        
        // Create iterator
        typedef std::map<int,std::vector<std::string>>::iterator it_type;
        for(it_type iterator = mColors.begin(); iterator != mColors.end(); iterator++) 
        {
            const int key = iterator->first;
            const std::vector<std::string> value = iterator->second;
            
            if (key != 0)
            {
                for (unsigned int it_name = 0; it_name < value.size(); it_name++)
                {
                    if ((value[it_name].find(rThisModelPart.Name()) != std::string::npos) == false) // TODO: Check this!!!
                    {
                        if (rThisModelPart.HasSubModelPart(value[it_name]) == false)
                        {
                            rThisModelPart.CreateSubModelPart(value[it_name]);
                        }
                    }
                }
            }
        }
        
        // Create a new model part
        /* NODES */
        unsigned int nnodes  = sizeof(mmgMesh->point);
        NodesArrayType& rNodes = rThisModelPart.Nodes();
        unsigned int node_id = rNodes.size();
        for (unsigned int i_node = 0; i_node < nnodes; i_node++)
        {
            node_id += 1;
            rThisModelPart.CreateNewNode(node_id, mmgMesh->point[i_node + 1].c[0], mmgMesh->point[i_node + 1].c[1], mmgMesh->point[i_node + 1].c[2]);
        }
       
        /* TRIANGLES */
        unsigned int ncond = sizeof(mmgMesh->tria);
        ConditionsArrayType& rConditions = rThisModelPart.Conditions();
        unsigned int cond_id = rConditions.size();
        std::string ConditionName = "Condition3D"; // NOTE: The condition should change the name to Condition3D3N for coherence
        std::vector<IndexType> ConditionNodeIds (3, 0);
        for (unsigned int i_cond = 0; i_cond < ncond; i_cond++)
        {
            cond_id += 1;
            ConditionNodeIds[0] = mmgMesh->tria[i_cond + 1].v[0];
            ConditionNodeIds[1] = mmgMesh->tria[i_cond + 1].v[1];
            ConditionNodeIds[2] = mmgMesh->tria[i_cond + 1].v[2];
            const unsigned int prop_id = mmgMesh->tria[i_cond + 1].ref;
//             rThisModelPart.AddProperties(pPropertiesVector[prop_id]);
//             PropertiesType::Pointer pProperties = rThisModelPart.pGetProperties(prop_id);
            rThisModelPart.CreateNewCondition(ConditionName ,cond_id, ConditionNodeIds, pProperties, 0);
           
            if (prop_id != 0) // NOTE: prop_id == 0 is the MainModelPart
            {
                std::vector<std::string> ColorList = mColors[prop_id];
                for (unsigned int colors = 0; colors < ColorList.size(); colors++)
                {
                    std::string SubModelPartName = ColorList[colors];
                    std::vector<IndexType> ConditionsIds (1, cond_id);
                    (rThisModelPart.GetSubModelPart(SubModelPartName)).AddConditions(ConditionsIds);
                }
            }
        }
       
        /* TETRAHEDRAS */
        unsigned int nelem = sizeof(mmgMesh->tetra);
        ElementsArrayType& rElements = rThisModelPart.Elements();
        unsigned int elem_id = rElements.size();
        std::string ElementName = "Element3D4N";
        std::vector<IndexType> ElementNodeIds (4, 0);
        for (unsigned int i_elem = 0; i_elem < nelem; i_elem++)
        {
            elem_id += 1;
            ElementNodeIds[0] = mmgMesh->tetra[i_elem + 1].v[0];
            ElementNodeIds[1] = mmgMesh->tetra[i_elem + 1].v[1];
            ElementNodeIds[2] = mmgMesh->tetra[i_elem + 1].v[2];
            ElementNodeIds[3] = mmgMesh->tetra[i_elem + 1].v[3];
            const unsigned int prop_id = mmgMesh->tetra[i_elem + 1].ref;
//             rThisModelPart.AddProperties(pPropertiesVector[prop_id]);
//             PropertiesType::Pointer pProperties = rThisModelPart.pGetProperties(prop_id);
            rThisModelPart.CreateNewElement(ElementName ,elem_id, ElementNodeIds, pProperties, 0);
           
            if (prop_id != 0) // NOTE: prop_id == 0 is the MainModelPart
            {
                std::vector<std::string> ColorList = mColors[prop_id];
                for (unsigned int colors = 0; colors < ColorList.size(); colors++)
                {
                    std::string SubModelPartName = ColorList[colors];
                    std::vector<IndexType> ElementsIds (1, elem_id);
                    (rThisModelPart.GetSubModelPart(SubModelPartName)).AddElements(ElementsIds);
                }
            }
        }
       
        //  Get the list of submodelparts names
        const std::vector<std::string> SubModelPartNames = rThisModelPart.GetSubModelPartNames();
       
        // Add the nodes to the differents submodelparts
        for (unsigned int i_model_part = 0; i_model_part < rThisModelPart.NumberOfSubModelParts(); i_model_part++)
        {
            ModelPart rSubModelPart = rThisModelPart.GetSubModelPart(SubModelPartNames[i_model_part]);
           
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
       
        // Save to file
        int ier = 0;
        if (save_to_file == true)
        {
            ier = SaveSolutionToFile(true);
        }
       
        // Free memory
        FreeMemory();
       
        return(ier);
    }
    
    /** ------------------------------ STEP III -------------------------- */
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * It saves the solution and mesh to files (for debugging pourpose g.e)
     * @param post_output: If the file to save is after or before remeshing
     */
    
    int SaveSolutionToFile(const bool post_output)
    {
        /** MMG library call */
        int ier = MMG3D_mmg3dlib(mmgMesh,mmgSol);
       
        if ( ier == MMG5_STRONGFAILURE ) 
        {
            std::cout << "BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH" << std::endl;
        }
        else if ( ier == MMG5_LOWFAILURE )
        {
            std::cout << "BAD ENDING OF MMG3DLIB" << std::endl;
        }
       
        /** get results */
        /** Two solutions: just use the MMG3D_saveMesh/MMG3D_saveSol functions
        that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
        using the MMG3D_getMesh/MMG3D_getSol functions */

        /** 1) Automatically save the mesh */
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
       
        /** a)  (not mandatory): give the ouptut mesh name using MMG3D_Set_outputMeshName
        (by default, the mesh is saved in the "mesh.o.mesh" file */
        MMG3D_Set_outputMeshName(mmgMesh,MeshFile);
        /** b) function calling */
        if ( MMG3D_saveMesh(mmgMesh,MeshFile) != 1 ) 
        {
            std::cout << "UNABLE TO SAVE MESH" << std::endl;
            return(MMG5_STRONGFAILURE);
        }

        /** 2) Automatically save the solution */
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
        /** a)  (not mandatory): give the ouptut sol name using MMG3D_Set_outputSolName
        (by default, the mesh is saved in the "mesh.o.sol" file */
        MMG3D_Set_outputSolName(mmgMesh,mmgSol,SolFile);
        /** b) function calling */
        if ( MMG3D_saveSol(mmgMesh,mmgSol,SolFile) != 1 ) 
        {
            std::cout << "UNABLE TO SAVE SOL" << std::endl;
            return(MMG5_LOWFAILURE);
        }
        
        return ier;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * It frees the memory used during all the process
     */
    
    void FreeMemory()
    {
        /** 3) Free the MMG3D5 structures */
        MMG3D_Free_all(
            MMG5_ARG_start,
            MMG5_ARG_ppMesh,
            &mmgMesh,
            MMG5_ARG_ppMet,
            &mmgSol,
            MMG5_ARG_end
            );

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
    
    // I/O information
    char* mFilename;
    std::string mStdStringFilename;
    
    // The member variables related with the MMG library
    MMG5_pMesh mmgMesh;
    MMG5_pSol  mmgSol;
    
    // Where the sub model parts IDs are stored
    std::map<int,std::vector<std::string>> mColors;
//     std::map<int,std::set<std::string>> mColors; // TODO: Maybe I need to change to set
    
    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{
    
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
        const int initial_index = 6 * (node_id - 1);
        const double coeff0 = 1.0/(element_size * element_size);
        
        if (scalar_value > value_threshold) // We don't reach the minimum
//         if (scalar_value > element_size) // NOTE: In the case that the scalar value is related with a distance, better to consider an independent variable
        {
            mmgSol->m[initial_index + 1] = coeff0;
            mmgSol->m[initial_index + 2] = 0.0;
            mmgSol->m[initial_index + 3] = coeff0;
            mmgSol->m[initial_index + 4] = 0.0;
            mmgSol->m[initial_index + 5] = 0.0;
            mmgSol->m[initial_index + 6] = coeff0;
        }
        else
        {
            const double aux = ratio + (scalar_value * element_size)*(1.0 - ratio);
            const double coeff1 = coeff0/(aux * aux);
            
            mmgSol->m[initial_index + 1] = coeff0*(1-gradient_value[0]*gradient_value[0]) + coeff1*gradient_value[0]*gradient_value[0];
            mmgSol->m[initial_index + 2] = coeff0*(gradient_value[0]*gradient_value[1]) + coeff1*gradient_value[0]*gradient_value[1];
            mmgSol->m[initial_index + 3] = coeff0*(1-gradient_value[1]*gradient_value[1]) + coeff1*gradient_value[1]*gradient_value[1];
            mmgSol->m[initial_index + 4] = coeff0*(gradient_value[0]*gradient_value[2]) + coeff1*gradient_value[0]*gradient_value[2];
            mmgSol->m[initial_index + 5] = coeff0*(gradient_value[1]*gradient_value[2]) + coeff1*gradient_value[1]*gradient_value[2];
            mmgSol->m[initial_index + 6] = coeff0*(1-gradient_value[2]*gradient_value[2]) + coeff1*gradient_value[2]*gradient_value[2];
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
