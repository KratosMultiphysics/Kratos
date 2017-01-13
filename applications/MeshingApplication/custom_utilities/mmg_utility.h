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
// #include "mmg/mmg2d/libmmg2d.h" // TODO: Version in 2D
#include "mmg/mmg3d/libmmg3d.h"
// #include "mmg/mmgs/libmmgs.h" // TODO: what is this libarry doing?

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
    
    /**
     * Instead of using an files already created we read an existing model part
     * @param rThisModelPart: The original model part, an auxiliary model part will be created an the input model part will be modified
     */
    void ComputeExistingModelPart(ModelPart& rThisModelPart)
    {
        // TODO: Add the reading of the model part (not bottleneck right now)
    }
    
    /**
     * We consider as input an already existing files , so this reads the input files 
     */
    void ReadFiles()
    {
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
    
    /**
     * The operation related with the libray is executed
     * @param rThisModelPart: The original model part, an auxiliary model part will be created an the input model part will be modified
     * @param rColors: The external dictionary needed to know how to assign the submodelparts
     * @param save_to_file: Save the solution to a *.sol and *.mesh files
     */
    bool Execute( // NOTE: Properties are going to give problems (we need to define new properies copying the old ones)
        ModelPart& rThisModelPart, //NOTE: The model part shoul be empty
        Properties::Pointer pProperties, // Solve the problem, just one property
//         std::vector<Properties::Pointer> pPropertiesVector,
        std::map<int,std::vector<std::string>> rColors,
        const bool save_to_file = false
        )
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
//            rThisModelPart.AddProperties(pPropertiesVector[prop_id]);
//            PropertiesType::Pointer pProperties = rThisModelPart.pGetProperties(prop_id);
           rThisModelPart.CreateNewCondition(ConditionName ,cond_id, ConditionNodeIds, pProperties, 0);
           
           if (prop_id != 0) // NOTE: prop_id == 0 is the MainModelPart
           {
               std::vector<std::string> ColorList = rColors[prop_id];
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
//            rThisModelPart.AddProperties(pPropertiesVector[prop_id]);
//            PropertiesType::Pointer pProperties = rThisModelPart.pGetProperties(prop_id);
           rThisModelPart.CreateNewElement(ElementName ,elem_id, ElementNodeIds, pProperties, 0);
           
           if (prop_id != 0) // NOTE: prop_id == 0 is the MainModelPart
           {
               std::vector<std::string> ColorList = rColors[prop_id];
               for (unsigned int colors = 0; colors < ColorList.size(); colors++)
               {
                   std::string SubModelPartName = ColorList[colors];
                   std::vector<IndexType> ElementsIds (1, elem_id);
                   (rThisModelPart.GetSubModelPart(SubModelPartName)).AddElements(ElementsIds);
               }
           }
       }
       
       const std::vector<std::string> SubModelPartNames = rThisModelPart.GetSubModelPartNames();
       
       // Add the nodes to the differents submodelparts //TODO: Get the list of submodelparts names
       for (unsigned int i_model_part = 0; i_model_part < rThisModelPart.NumberOfSubModelParts(); i_model_part++)
       {
           ModelPart rSubModelPart = rThisModelPart.GetSubModelPart(SubModelPartNames[i_model_part]);
           
           std::vector<IndexType> NodesIds;
           
           for (ElementConstantIterator elem_iterator = rSubModelPart.ElementsBegin(); elem_iterator != rSubModelPart.ElementsEnd(); elem_iterator++)
           {
               for (unsigned int i_node = 0; i_node < elem_iterator->GetGeometry().size(); i_node++)
               {
                   NodesIds.push_back(elem_iterator->GetGeometry()[i_node].Id());
               }
           }
           
           for (ConditionConstantIterator condition_iterator = rSubModelPart.ConditionsBegin(); condition_iterator != rSubModelPart.ConditionsEnd(); condition_iterator++)
           {
               for (unsigned int i_node = 0; i_node < condition_iterator->GetGeometry().size(); i_node++)
               {
                   NodesIds.push_back(condition_iterator->GetGeometry()[i_node].Id());
               }
           }
           
           // Clean duplicated nodes // https://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector
           std::set<int> s;
           const unsigned size = NodesIds.size();
           for( unsigned i = 0; i < size; ++i ) s.insert( NodesIds[i] );
           {
               NodesIds.assign( s.begin(), s.end() );
           }
           
           rSubModelPart.AddNodes(NodesIds);
       }
       
       // Save to file
       if (save_to_file == true)
       {
           ier = SaveSolutionToFile(ier);
       }
       
       // Free memory
       FreeMemory();
       
       return(ier);
    }
    
    /** ------------------------------ STEP III -------------------------- */
    int SaveSolutionToFile(int ier)
    {
        /** get results */
        /** Two solutions: just use the MMG3D_saveMesh/MMG3D_saveSol functions
        that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
        using the MMG3D_getMesh/MMG3D_getSol functions */

        /** 1) Automatically save the mesh */
        std::string MeshName = mStdStringFilename+".o.mesh";
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
        std::string SolName = mStdStringFilename+".o.sol";
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
    
    char* mFilename;
    std::string mStdStringFilename;
    
    // The member variables related with the MMG library
    MMG5_pMesh mmgMesh;
    MMG5_pSol  mmgSol;
    
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
};// class MmgUtility
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}// namespace Kratos.
#endif /* KRATOS_MMG_UTILITY defined */
