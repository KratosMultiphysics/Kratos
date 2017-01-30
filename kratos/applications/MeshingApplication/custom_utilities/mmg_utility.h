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

    /**
     * This enums are used to simplify the computation of the std::vector containing the conditions and elements
     */
    
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

//// This class is a remesher which uses the MMG library 
// The class uses a class for the 2D and 3D cases 

template<unsigned int TDim>  
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
     
    // Constructor
    
    /**
     * This is the default constructor, which is used to read the input files 
     * @param Filename: The input name of the output file
     * @param echo_level: The level of verbosity
     */
    
    MmgUtility(
        const std::string Filename = "output_remesh",
        const unsigned int echo_level = 3
        )
        :mStdStringFilename(Filename),
        mEchoLevel(echo_level)
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
       
       int verbosityMMG;  // TODO: Discover WTF means each value of the verbosity
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
    
    /**
     * Instead of using an files already created we read an existing model part
     * @param rThisModelPart: The original model part, an auxiliary model part will be created an the input model part will be modified
     */
    
    void RemeshModelPart(
        ModelPart& rThisModelPart,
        const bool save_to_file = false,
        const unsigned int MaxNumberOfResults = 1000
        )
    {
        if (mEchoLevel > 0)
        {
            std::cout << "//---------------------------------------------------//" << std::endl;
            std::cout << "//---------------------------------------------------//" << std::endl;
            std::cout << "//---------------  BEFORE REMESHING   ---------------//" << std::endl;
            std::cout << "//---------------------------------------------------//" << std::endl;
            std::cout << "//---------------------------------------------------//" << std::endl;
            std::cout << std::endl;
            
            KRATOS_WATCH(rThisModelPart);
        }
        
        // We initialize the mesh and solution data
        InitializeMeshData(rThisModelPart);
        InitializeSolData(rThisModelPart);
        
        // Check if the number of given entities match with mesh size 
        CheckMeshData();
        
        // Save to file
        if (save_to_file == true)
        {
            SaveSolutionToFile(false);
        }
        
        // We execute the remeshing
        ExecuteRemeshing(rThisModelPart, save_to_file, MaxNumberOfResults);
        
        /* We print the resulting model part */
        if (mEchoLevel > 0)
        {
            std::cout << "//---------------------------------------------------//" << std::endl;
            std::cout << "//---------------------------------------------------//" << std::endl;
            std::cout << "//---------------   AFTER REMESHING   ---------------//" << std::endl;
            std::cout << "//---------------------------------------------------//" << std::endl;
            std::cout << "//---------------------------------------------------//" << std::endl;
            std::cout << std::endl;
            
            KRATOS_WATCH(rThisModelPart);
        }
        
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
    unsigned int mEchoLevel;
    
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
     * This function generates the mesh MMG5 structure from a Kratos Model Part
     * @param rThisModelPart: The original model part, an auxiliary model part will be created an the input model part will be modified
     */
    
    void InitializeMeshData(ModelPart& rThisModelPart)
    {
        // First we compute the colors
        std::map<int,int> node_colors, cond_colors, elem_colors;
        ComputeColors(rThisModelPart, node_colors, cond_colors, elem_colors);
        
        /////////* MESH FILE */////////
        // Build mesh in MMG5 format //
        
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
                
                if ((itElem->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) // Tetrahedron
                {
                    numArrayElements[0] += 1;
                }
                else if ((itElem->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Prism3D6) // Prisms
                {
                    numArrayElements[1] += 1;
                }
                else
                {
                    std::cout << "WARNING: YOUR GEOMETRY CONTAINS HEXAEDRON THAT CAN NOT BE REMESHED" << std::endl;
                }
            }
            
            if (((numArrayElements[0] + numArrayElements[1]) < numElements) && mEchoLevel > 0)
            {
                std::cout << "Number of Elements: " << numElements << " Number of Tetrahedron: " << numArrayElements[0] << " Number of Prisms: " << numArrayElements[1] << std::endl;
            }
            
            /* Conditions */
            for(unsigned int i = 0; i < numConditions; i++) 
            {
                auto itCond = pConditions.begin() + i;
                
                if ((itCond->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) // Triangles
                {
                    numArrayConditions[0] += 1;
                }
                else if ((itCond->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4)  // Quadrilaterals
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
        // We clone the first condition of each type
        if (TDim == 2)
        {
            const cond_geometries_2d index_geom0 = Line;
            mpRefCondition[index_geom0] = pConditions.begin()->Create(0, pConditions.begin()->GetGeometry(), pConditions.begin()->pGetProperties());
        }
//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCond = pConditions.begin() + i;
            
            // We clone the first condition of each type
            if (TDim == 3)
            {                
                const cond_geometries_3d index_geom0 = Triangle3D;
                const cond_geometries_3d index_geom1 = Quadrilateral3D;

                if ((itCond->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3 && mInitRefCondition[index_geom0] == false) // Triangle
                {
                    mpRefCondition[index_geom0] = itCond->Create(0, itCond->GetGeometry(), itCond->pGetProperties());
                    mInitRefCondition[index_geom0] = true;
                }
                else if ((itCond->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4 && mInitRefCondition[index_geom1] == false) // Quadrilateral
                {
                    mpRefCondition[index_geom1] = itCond->Create(0, itCond->GetGeometry(), itCond->pGetProperties());
                    mInitRefCondition[index_geom1] = true;
                }
            }
            
            SetConditions(itCond->GetGeometry(), cond_colors[itCond->Id()], i + 1);
        }
        
        /* Elements */
        // We clone the first element of each type
        if (TDim == 2)
        {
            const elem_geometries_2d index_geom0 = Triangle2D;
            mpRefElement[index_geom0] = pElements.begin()->Create(0, pElements.begin()->GetGeometry(), pElements.begin()->pGetProperties());
        }
//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numElements; i++) 
        {
            auto itElem = pElements.begin() + i;
            
            // We clone the first element of each type
            if (TDim == 3)
            {                
                const elem_geometries_3d index_geom0 = Tetrahedra;
                const elem_geometries_3d index_geom1 = Prism;
                
                if ((itElem->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4 && mInitRefElement[index_geom0] == false) // Tetrahedra
                {
                    mpRefElement[index_geom0] = itElem->Create(0, itElem->GetGeometry(), itElem->pGetProperties());
                    mInitRefElement[index_geom0] = true;
                }
                else if ((itElem->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4 && mInitRefElement[index_geom1] == false) // Prism
                {
                    mpRefElement[index_geom1] = itElem->Create(0, itElem->GetGeometry(), itElem->pGetProperties());
                    mInitRefElement[index_geom1] = true;
                }
            }
            
            SetElements(itElem->GetGeometry(), elem_colors[itElem->Id()], i + 1);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * We initialize the metrics of the MMG sol using a level set approach
     * @param rThisModelPart: The original model part where we compute the metric
     */
    
    void InitializeSolData(ModelPart& rThisModelPart)
    {
        ////////* SOLUTION FILE *////////
        
        // Iterate in the nodes
        NodesArrayType& pNode = rThisModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        SetSolSize(numNodes);

//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            #ifdef KRATOS_DEBUG 
            if( itNode->Has(MMG_METRIC) == false) 
            {
                KRATOS_ERROR <<  " MMG_METRIC not defined for node " << itNode->Id();
            }
            #endif     
            
            // We get the metric
            const Vector& metric = itNode->GetValue(MMG_METRIC);
            
            #ifdef KRATOS_DEBUG 
            if(metric.size() != TDim * 3 - 3) 
            {
                KRATOS_ERROR << "Wrong size of vector MMG_METRIC found for node " << itNode->Id() << " size is " << metric.size() << " expected size was " << TDim * 3 - 3;
            }
            #endif
            
            // We set the metric
            SetMetricTensor(metric, i + 1);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * We execute the MMg library and build the new model part from the old model part
     * @param rThisModelPart: The original model part, an auxiliary model part will be created an the input model part will be modified
     * @param MaxNumberOfResults: The maxim number of results to consider in the search
     */
    
    void ExecuteRemeshing(
        ModelPart& rThisModelPart,
        const bool save_to_file = false,
        const unsigned int MaxNumberOfResults = 1000
        )
    {
        // We initialize some values
        unsigned int step_data_size = rThisModelPart.GetNodalSolutionStepDataSize();
        unsigned int buffer_size    = rThisModelPart.NodesBegin()->GetBufferSize();
        
        if (mEchoLevel > 0)
        {        
            std::cout << "Step data size: " << step_data_size << " Buffer size: " << buffer_size <<     std::endl; 
        }
        
        ////////* MMG LIBRARY CALL *////////
        if (mEchoLevel > 0)
        {
            std::cout << "////////* MMG LIBRARY CALL *////////" << std::endl; 
        }
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
        
        if (mEchoLevel > 0)
        {
            std::cout << "     Nodes created: " << nNodes << std::endl;
            if (TDim == 2) // 2D
            {
                std::cout << "Conditions created: " << nConditions[0] << std::endl;
                std::cout << "Elements created: " << nElements[0] << std::endl;
            }
            else // 3D
            {
                std::cout << "Conditions created: " << nConditions[0] + nConditions[1] << std::endl;
                std::cout << "\tTriangles: " << nConditions[0] << "\tQuadrilaterals: " << nConditions[1]<< std::endl;
                std::cout << "Elements created: " << nElements[0] + nElements[1] << std::endl;
                std::cout << "\tTetrahedron: " << nElements[0] << "\tPrisms: " << nElements[1] << std::endl;
            }
        }
        
        ////////* EMPTY AND BACKUP THE MODEL PART *////////
        
        ModelPart rOldModelPart;
        
        // First we empty the model part
        for (NodeConstantIterator node_iterator = rThisModelPart.NodesBegin(); node_iterator != rThisModelPart.NodesEnd(); node_iterator++)
        {
            node_iterator->Set(TO_ERASE, true);
            rOldModelPart.AddNode(*(node_iterator.base()));
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
            rOldModelPart.AddElement(*(elem_iterator.base()));
        }
        rThisModelPart.RemoveElementsFromAllLevels(TO_ERASE);  
        
        // Create a new model part
        /* NODES */
        for (int unsigned i_node = 1; i_node <= nNodes; i_node++)
        {
            NodeType::Pointer pNode;
            
            if (TDim == 2)
            {
                pNode = rThisModelPart.CreateNewNode(i_node, mmgMesh->point[i_node].c[0], mmgMesh->point[i_node].c[1], 0.0);
            }
            else
            {
                pNode = rThisModelPart.CreateNewNode(i_node, mmgMesh->point[i_node].c[0], mmgMesh->point[i_node].c[1], mmgMesh->point[i_node].c[2]);
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
            ConditionType::Pointer pCondition;
            int prop_id, isRequired;
            for (int unsigned i_cond = 1; i_cond <= nConditions[0]; i_cond++)
            {
                cond_id += 1;
                
                if (TDim == 2) // Lines // TODO: Move this check outside the loop
                {
                    const cond_geometries_2d index_geom = Line;
                    
                    int edge0, edge1, isRidge;
                    
                    if (MMG2D_Get_edge(mmgMesh, &edge0, &edge1, &prop_id, &isRidge, &isRequired) != 1 )
                    {
                        exit(EXIT_FAILURE);
                    }
                    
                    std::vector<NodeType::Pointer> ConditionNodes (2);
                    ConditionNodes[0] = rThisModelPart.pGetNode(edge0);
                    ConditionNodes[1] = rThisModelPart.pGetNode(edge1);    
                    
                    pCondition = mpRefCondition[index_geom]->Create(cond_id, ConditionNodes, mpRefCondition[index_geom]->pGetProperties());
                }
                else // Triangles
                {
                    const cond_geometries_3d index_geom = Triangle3D;
                    
                    int vertex0, vertex1, vertex2;

                    if (MMG3D_Get_triangle(mmgMesh, &vertex0, &vertex1, &vertex2, &prop_id, &isRequired) != 1 )
                    {
                        exit(EXIT_FAILURE);
                    }
                    
                    std::vector<NodeType::Pointer> ConditionNodes (3);
                    ConditionNodes[0] = rThisModelPart.pGetNode(vertex0);
                    ConditionNodes[1] = rThisModelPart.pGetNode(vertex1);
                    ConditionNodes[2] = rThisModelPart.pGetNode(vertex2);
                    
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
                ConditionType::Pointer pCondition;
                int prop_id, isRequired;
                for (int unsigned i_cond = 1; i_cond <= nConditions[1]; i_cond++)
                {
                    cond_id += 1;
                    
                    const cond_geometries_3d index_geom = Quadrilateral3D;
                    
                    int vertex0, vertex1, vertex2, vertex3;

                    if (MMG3D_Get_quadrilateral(mmgMesh, &vertex0, &vertex1, &vertex2, &vertex3, &prop_id, &isRequired) != 1 )
                    {
                        exit(EXIT_FAILURE);
                    }
                    
                    std::vector<NodeType::Pointer> ConditionNodes (4);
                    ConditionNodes[0] = rThisModelPart.pGetNode(vertex0);
                    ConditionNodes[1] = rThisModelPart.pGetNode(vertex1);
                    ConditionNodes[2] = rThisModelPart.pGetNode(vertex2);
                    ConditionNodes[3] = rThisModelPart.pGetNode(vertex3);
                    
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
            ElementType::Pointer pElement;
            int prop_id, isRequired;
            for (int unsigned i_elem = 1; i_elem <= nElements[0]; i_elem++)
            {
                elem_id += 1;
                
                if (TDim == 2) // Triangle
                {
                    const elem_geometries_2d index_geom = Triangle2D;
                    
                    int vertex0, vertex1, vertex2;
                    
                    if (MMG2D_Get_triangle(mmgMesh, &vertex0, &vertex1, &vertex2, &prop_id, &isRequired) != 1 )
                    {
                        exit(EXIT_FAILURE);
                    }

                    std::vector<NodeType::Pointer> ElementNodes (3);
                    ElementNodes[0] = rThisModelPart.pGetNode(vertex0);
                    ElementNodes[1] = rThisModelPart.pGetNode(vertex1);
                    ElementNodes[2] = rThisModelPart.pGetNode(vertex2);
                    
                    pElement = mpRefElement[index_geom]->Create(elem_id, ElementNodes, mpRefElement[index_geom]->pGetProperties());
                }
                else // Tetrahedra
                {
                    const elem_geometries_3d index_geom = Tetrahedra;
                    
                    int vertex0, vertex1, vertex2, vertex3;
                    
                    if (MMG3D_Get_tetrahedron(mmgMesh, &vertex0, &vertex1, &vertex2, &vertex3, &prop_id, &isRequired) != 1 )
                    {
                        exit(EXIT_FAILURE);
                    }
                    
                    std::vector<NodeType::Pointer> ElementNodes (4);
                    ElementNodes[0] = rThisModelPart.pGetNode(vertex0);
                    ElementNodes[1] = rThisModelPart.pGetNode(vertex1);
                    ElementNodes[2] = rThisModelPart.pGetNode(vertex2);
                    ElementNodes[3] = rThisModelPart.pGetNode(vertex3);
                    
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
                ElementType::Pointer pElement;
                int prop_id, isRequired;
                for (int unsigned i_elem = 1; i_elem <= nElements[1]; i_elem++)
                {
                    elem_id += 1;
                    
                    const elem_geometries_3d index_geom = Prism;
                    
                    int vertex0, vertex1, vertex2, vertex3, vertex4, vertex5;
                    
                    if (MMG3D_Get_prism(mmgMesh, &vertex0, &vertex1, &vertex2, &vertex3, &vertex4, &vertex5, &prop_id, &isRequired) != 1 )
                    {
                        exit(EXIT_FAILURE);
                    }
                    
                    std::vector<NodeType::Pointer> ElementNodes (6);
                    ElementNodes[0] = rThisModelPart.pGetNode(vertex0);
                    ElementNodes[1] = rThisModelPart.pGetNode(vertex1);
                    ElementNodes[2] = rThisModelPart.pGetNode(vertex2);
                    ElementNodes[3] = rThisModelPart.pGetNode(vertex3);
                    ElementNodes[4] = rThisModelPart.pGetNode(vertex4);
                    ElementNodes[5] = rThisModelPart.pGetNode(vertex5);
                    
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
        
        /* Save to file */
        if (save_to_file == true)
        {
            SaveSolutionToFile(true);
        }
       
        /* Free memory */
        FreeMemory();
        
        /* We interpolate all the values */
        InterpolateValues(rThisModelPart, rOldModelPart, MaxNumberOfResults, step_data_size, buffer_size);
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * It interpolates the values in the new model part using the old model part
     * @param rThisModelPart: The new model part
     * @param rOldModelPart: The old model part
     * @param MaxNumberOfResults: The maxim number of results to consider in the search
     * @param step_data_size: The size of the database
     * @param buffer_size: The size of the buffer
     */
    
    void InterpolateValues(
        ModelPart& rThisModelPart,
        ModelPart& rOldModelPart,
        const unsigned int MaxNumberOfResults,
        unsigned int step_data_size,
        unsigned int buffer_size  
        )
    {
        // We create the locator
        BinBasedFastPointLocator<TDim> PointLocator = BinBasedFastPointLocator<TDim>(rOldModelPart);
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
            
            if ((found == false) && mEchoLevel > 0)
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
    }
    
    /**
     * It saves the solution and mesh to files (for debugging pourpose g.e)
     * @param post_output: If the file to save is after or before remeshing
     */
    
    void SaveSolutionToFile(const bool post_output)
    {
        /* GET RESULTS */

        // Automatically save the mesh 
        OutputMesh(post_output);

        // Automatically save the solution 
        OutputSol(post_output);
    }
    
    /**
     * It frees the memory used during all the process
     */
    
    void FreeMemory()
    {
        // Free the MMG3D5 structures 
        FreeAll();

        free(mFilename);
        mFilename = NULL;
    }
    
    /** Initialisation of mesh and sol structures 
    * args of InitMesh:
    * MMG5_ARG_start: we start to give the args of a variadic func
    * MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
    * &mmgMesh: pointer toward your MMG5_pMesh (that store your mesh)
    * MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
    * &mmgSol: pointer toward your MMG5_pSol (that store your metric) 
    */
    
    void InitMesh();
    
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
        );
    
    /**
     * This sets the size of the solution
     * @param numNodes: Number of nodes
     */
    
    void SetSolSize(const int numNodes);
    
    /**
     * This checks the mesh data and prints if it is OK
     */
    
    void CheckMeshData();
    
    /**
     * This sets the output mesh
     */
    
    void OutputMesh(const bool post_output);
    
    /**
     * This sets the output sol
     */
    
    void OutputSol(const bool post_output);
    
    /**
     * This loads the solution
     */
    
    void MMGLibCall();
    
    /**
     * This frees the MMG structures
     */
    
    void FreeAll();
    
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
        );
    
    /**
     * This sets the conditions of the mesh
     * @param Geom: The geometry of the condition
     * @param color: Reference of the node(submodelpart)
     * @param index: The index number of the node 
     */
    
    void SetConditions(
        Geometry<Node<3> > & Geom,
        const int color,
        const int index
        );
    
    /**
     * This sets elements of the mesh
     * @param Geom: The geometry of the element
     * @param color: Reference of the node(submodelpart)
     * @param index: The index number of the node 
     */
    
    void SetElements(
        Geometry<Node<3> > & Geom,
        const int color,
        const int index
        );
    
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

    /**
     * This function is used to compute the Hessian metric tensor, note that when using the Hessian, more than one metric can be defined simultaneously, so in consecuence we need to define the elipsoid which defines the volume of maximal intersection
     * @param metric: This array contains the components of the metric tensor in the MMG defined order
     */

    void SetMetricTensor(
        const Vector& metric,
        const int node_id 
        );
    
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

///@name Explicit Specializations
///@{

    template<>  
    void MmgUtility<2>::InitMesh()
    {       
        MMG2D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end); 
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<3>::InitMesh()
    {       
        MMG3D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end); 
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<2>::SetMeshSize(
        const int numNodes,
        const array_1d<int, 1> numArrayElements, 
        const array_1d<int, 1> numArrayConditions
        )
    {
        //Give the size of the mesh: numNodes vertices, numElements triangles, numConditions edges (2D) 
        if ( MMG2D_Set_meshSize(mmgMesh, numNodes, numArrayElements[0], numArrayConditions[0]) != 1 ) 
        {
            exit(EXIT_FAILURE);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<3>::SetMeshSize(
        const int numNodes,
        const array_1d<int, 2> numArrayElements,  // NOTE: We do this tricky thing to take into account the prisms
        const array_1d<int, 2> numArrayConditions // NOTE: We do this tricky thing to take into account the quadrilaterals
        )
    {
        //Give the size of the mesh: numNodes vertex, numElements tetra and prism, numArrayConditions triangles and quadrilaterals, 0 edges (3D) 
        if ( MMG3D_Set_meshSize(mmgMesh, numNodes, numArrayElements[0], numArrayElements[1], numArrayConditions[0], numArrayConditions[1], 0) != 1 ) 
        {
            exit(EXIT_FAILURE);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    void MmgUtility<2>::SetSolSize(const int numNodes)
    {
        if ( MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,numNodes,MMG5_Tensor) != 1 )
        {
            exit(EXIT_FAILURE);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    void MmgUtility<3>::SetSolSize(const int numNodes)
    {
        if ( MMG3D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,numNodes,MMG5_Tensor) != 1 )
        {
            exit(EXIT_FAILURE);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    void MmgUtility<2>::CheckMeshData()
    {
        if ( MMG2D_Chk_meshData(mmgMesh, mmgSol) != 1 ) 
        {
            exit(EXIT_FAILURE);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    void MmgUtility<3>::CheckMeshData()
    {
        if ( MMG3D_Chk_meshData(mmgMesh, mmgSol) != 1 ) 
        {
            exit(EXIT_FAILURE);
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<2>::OutputMesh(const bool post_output)
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
        
        // a)  Give the ouptut mesh name using MMG2D_Set_outputMeshName (by default, the mesh is saved in the "mesh.o.mesh" file  
        MMG2D_Set_outputMeshName(mmgMesh,MeshFile);

        // b) function calling 
        if ( MMG2D_saveMesh(mmgMesh,MeshFile) != 1) 
        {
            std::cout << "UNABLE TO SAVE MESH" << std::endl;
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<3>::OutputMesh(const bool post_output)
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
        
        // a)  Give the ouptut mesh name using MMG3D_Set_outputMeshName (by default, the mesh is saved in the "mesh.o.mesh" file 
        MMG3D_Set_outputMeshName(mmgMesh,MeshFile);

        // b) function calling 
        if ( MMG3D_saveMesh(mmgMesh,MeshFile) != 1) 
        {
            std::cout << "UNABLE TO SAVE MESH" << std::endl;
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    void MmgUtility<2>::OutputSol(const bool post_output)
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
        
        // a)  Give the ouptut sol name using MMG2D_Set_outputSolName (by default, the mesh is saved in the "mesh.o.sol" file 
        MMG2D_Set_outputSolName(mmgMesh, mmgSol, SolFile);

        // b) Function calling 
        if ( MMG2D_saveSol(mmgMesh, mmgSol, SolFile) != 1) 
        {
            std::cout << "UNABLE TO SAVE SOL" << std::endl;
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    void MmgUtility<3>::OutputSol(const bool post_output)
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
        
        // a)  Give the ouptut sol name using MMG3D_Set_outputSolName (by default, the mesh is saved in the "mesh.o.sol" file 
        MMG3D_Set_outputSolName(mmgMesh, mmgSol, SolFile);

        // b) Function calling 
        if ( MMG3D_saveSol(mmgMesh,mmgSol, SolFile) != 1) 
        {
            std::cout << "UNABLE TO SAVE SOL" << std::endl;
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    void MmgUtility<2>::MMGLibCall()
    {
        const int ier = MMG2D_mmg2dlib(mmgMesh, mmgSol);

        if ( ier == MMG5_STRONGFAILURE ) 
        {
            KRATOS_THROW_ERROR( std::logic_error, "WARNING: BAD ENDING OF MMG2DLIB: UNABLE TO SAVE MESH. ier: ", ier );
        }
        else if ( ier == MMG5_LOWFAILURE )
        {
            KRATOS_THROW_ERROR( std::logic_error, "WARNING: BAD ENDING OF MMG2DLIB. ier: ", ier );
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    void MmgUtility<3>::MMGLibCall()
    {
        const int ier = MMG3D_mmg3dlib(mmgMesh, mmgSol);

        if ( ier == MMG5_STRONGFAILURE ) 
        {
            KRATOS_THROW_ERROR( std::logic_error, "WARNING: BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH. ier: ", ier );
        }
        else if ( ier == MMG5_LOWFAILURE )
        {
            KRATOS_THROW_ERROR( std::logic_error, "WARNING: BAD ENDING OF MMG3DLIB. ier: ", ier );
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    void MmgUtility<2>::FreeAll()
    {
        MMG2D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,MMG5_ARG_end);
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    void MmgUtility<3>::FreeAll()
    {
        MMG3D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,MMG5_ARG_end);
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    void MmgUtility<2>::SetNodes(
        const double X,
        const double Y,
        const double Z,
        const int color,
        const int index
        )
    {
        if ( MMG2D_Set_vertex(mmgMesh, X, Y, color, index) != 1 )  
        {
            exit(EXIT_FAILURE);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    void MmgUtility<3>::SetNodes(
        const double X,
        const double Y,
        const double Z,
        const int color,
        const int index
        )
    {
        if ( MMG3D_Set_vertex(mmgMesh, X, Y, Z, color, index) != 1 )  
        {
            exit(EXIT_FAILURE); 
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<2>::SetConditions(
        Geometry<Node<3> > & Geom,
        const int color,
        const int index
        )
    {
        const int id1 = Geom[0].Id(); // First node id
        const int id2 = Geom[1].Id(); // Second node id

        if ( MMG2D_Set_edge(mmgMesh, id1, id2, color, index) != 1 ) 
        {
            exit(EXIT_FAILURE);
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<3>::SetConditions(
        Geometry<Node<3> > & Geom,
        const int color,
        const int index
        )
    {
        const int id1 = Geom[0].Id(); // First node id
        const int id2 = Geom[1].Id(); // Second node id
        const int id3 = Geom[2].Id(); // Third node id
        
        if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) // Triangle
        {
            if ( MMG3D_Set_triangle(mmgMesh, id1, id2, id3, color, index) != 1 )  
            {
                exit(EXIT_FAILURE); 
            }
        }
        else if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4) // Quadrilaterals
        {
            const int id4 = Geom[3].Id(); // Fourth node id
            
            if ( MMG3D_Set_quadrilateral(mmgMesh, id1, id2, id3, id4, color, index) != 1 )  
            {
                exit(EXIT_FAILURE); 
            }
        }
        else
        {
            const unsigned int size_geom = Geom.size();
            KRATOS_THROW_ERROR( std::logic_error, "WARNING: I DO NOT KNOW WHAT IS THIS. Size: ", size_geom );
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<2>::SetElements(
        Geometry<Node<3> > & Geom,
        const int color,
        const int index
        )
    {
        const int id1 = Geom[0].Id(); // First node id
        const int id2 = Geom[1].Id(); // Second node id
        const int id3 = Geom[2].Id(); // Third node id
        
        if ( MMG2D_Set_triangle(mmgMesh, id1, id2, id3, color, index) != 1 ) 
        {
            exit(EXIT_FAILURE);
        }

    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<3>::SetElements(
        Geometry<Node<3> > & Geom,
        const int color,
        const int index
        )
    {
        const int id1 = Geom[0].Id(); // First node id
        const int id2 = Geom[1].Id(); // Second node id
        const int id3 = Geom[2].Id(); // Third node id
        const int id4 = Geom[3].Id(); // Fourth node id
        
        if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) // Tetrahedron
        {
            if ( MMG3D_Set_tetrahedron(mmgMesh, id1, id2, id3, id4, color, index) != 1 )  
            {
                exit(EXIT_FAILURE); 
            }
        }
        else if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Prism3D6) // Prisms
        {
            const int id5 = Geom[4].Id(); // 5th node id
            const int id6 = Geom[5].Id(); // 6th node id
            
            if ( MMG3D_Set_prism(mmgMesh, id1, id2, id3, id4, id5, id6, color, index) != 1 )  
            {
                exit(EXIT_FAILURE); 
            }
        }
        else if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8) // Hexaedron
        {
//                 const int id5 = Geom[4].Id(); // 5th node id
//                 const int id6 = Geom[5].Id(); // 6th node id
//                 const int id6 = Geom[7].Id(); // 7th node id
//                 const int id6 = Geom[8].Id(); // 8th node id
            
            const unsigned int size_geom = Geom.size();
            KRATOS_THROW_ERROR( std::logic_error, "WARNING: HEXAEDRON NON IMPLEMENTED IN THE LIBRARY", size_geom);
        }
        else
        {
            const unsigned int size_geom = Geom.size();
            KRATOS_THROW_ERROR( std::logic_error, "WARNING: I DO NOT KNOW WHAT IS THIS. Size: ", size_geom );
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<2>::SetMetricTensor(
        const Vector& metric,
        const int node_id 
        )
    {
        if ( MMG2D_Set_tensorSol(mmgSol, metric[0],  metric[1], metric[2], node_id) != 1 )
        {
            exit(EXIT_FAILURE);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<3>::SetMetricTensor(
        const Vector& metric,
        const int node_id 
        )
    {
        if ( MMG3D_Set_tensorSol(mmgSol, metric[0], metric[1], metric[2], metric[3], metric[4], metric[5], node_id) != 1 )
        {
            exit(EXIT_FAILURE);
        }
    }


    
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}// namespace Kratos.
#endif /* KRATOS_MMG_UTILITY defined */
