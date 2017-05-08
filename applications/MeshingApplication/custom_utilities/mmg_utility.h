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
#include "includes/kratos_parameters.h"
#include "containers/variables_list.h"
#include "utilities/openmp_utils.h"
#include "input_output/logger.h"
#include <set>
#include <map>
// The includes related with the MMG library
#include "mmg/libmmg.h"
#include "mmg/mmg2d/libmmg2d.h" 
#include "mmg/mmg3d/libmmg3d.h"
#include "mmg/mmgs/libmmgs.h"
// We indlude the internal variable interpolation process
#include "custom_processes/nodal_values_interpolation_process.h"
#include "custom_processes/internal_variables_interpolation_process.h"
// Include the point locator
#include "utilities/binbased_fast_point_locator.h"
// Include the spatial containers needed for search
#include "spatial_containers/spatial_containers.h" // kd-tree 

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

    // Containers definition
    typedef ModelPart::NodesContainerType                        NodesArrayType;
    typedef ModelPart::ElementsContainerType                  ElementsArrayType;
    typedef ModelPart::ConditionsContainerType              ConditionsArrayType;
    
    // Components definition
    typedef Node <3>                                                   NodeType;
    typedef Properties                                           PropertiesType;
    typedef Element                                                 ElementType;
    typedef Condition                                             ConditionType;
    
    // Index defintion
    typedef std::size_t                                               IndexType;
    typedef std::size_t                                                SizeType;
    
    // DoF definition
    typedef Dof<double>                                                 DofType;
    
    // Mesh definition
    typedef Mesh<NodeType, PropertiesType, ElementType, ConditionType> MeshType;
    typedef MeshType::PropertiesContainerType           PropertiesContainerType;
    typedef MeshType::NodeConstantIterator                 NodeConstantIterator;
    typedef MeshType::ConditionConstantIterator       ConditionConstantIterator;
    typedef MeshType::ElementConstantIterator           ElementConstantIterator;
    
    #if !defined(KEY_COMPAROR_VECTOR)
    #define KEY_COMPAROR_VECTOR
    template<class TClassType>
    struct KeyComparorVector
    {
        bool operator()(const vector<TClassType>& lhs, const vector<TClassType>& rhs) const
        {
            if(lhs.size() != rhs.size())
            {
                return false;
            }

            for(unsigned int i=0; i<lhs.size(); i++)
            {
                if(lhs[i] != rhs[i]) 
                {
                    return false;
                }
            }

            return true;
        }
    };
    #endif
    
    #if !defined(KEY_HASHER_VECTOR)
    #define KEY_HASHER_VECTOR
    template<class TClassType>
    struct KeyHasherVector
    {
        std::size_t operator()(const vector<TClassType>& k) const
        {
            return boost::hash_range(k.begin(), k.end());
        }
    };
    #endif
    
///@}
///@name  Enum's
///@{

    /**
     * This enums are used to simplify the computation of the std::vector containing the conditions and elements
     */
    
    enum CondGeometries2D {Line = 0};
    
    enum ElemGeometries2D {Triangle2D = 0};
    
    enum CondGeometries3D {Triangle3D = 0, Quadrilateral3D = 1};
    
    enum ElemGeometries3D {Tetrahedra = 0, Prism = 1};
    
    #if !defined(FRAMEWORK_EULER_LAGRANGE)
    #define FRAMEWORK_EULER_LAGRANGE
        enum FrameworkEulerLagrange {Eulerian = 0, Lagrangian = 1};
    #endif
    
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
    
    ///@}
    ///@name Life Cycle
    ///@{
     
    // Constructor
    
    /**
     * This is the default constructor, which is used to read the input files 
     * @param rThisModelPart: The model part
     * @param ThisParameters: The parameters
     */
    
    MmgUtility(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        )
        :mrThisModelPart(rThisModelPart),
         mThisParameters(ThisParameters)
    {       
        Parameters DefaultParameters = Parameters(R"(
            {
                "filename"                         : "out",
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
                "max_number_of_searchs"            : 1000,
                "echo_level"                       : 3,
                "step_data_size"                   : 0,
                "buffer_size"                      : 0
            })" );
        
        mThisParameters.ValidateAndAssignDefaults(DefaultParameters);
        
        mStdStringFilename = mThisParameters["filename"].GetString();
        mEchoLevel = mThisParameters["echo_level"].GetInt();
        
        mFilename = new char [mStdStringFilename.length() + 1];
        std::strcpy (mFilename, mStdStringFilename.c_str());
       
        mFramework = ConvertFramework(mThisParameters["framework"].GetString());
       
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
     */
    
    void RemeshModelPart()
    {        
        const bool SaveToFile = mThisParameters["save_external_files"].GetBool();
        
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
        if (SaveToFile == true)
        {
            SaveSolutionToFile(false);
        }
        
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
    
    // The model part to compute
    ModelPart& mrThisModelPart;                      
    
    // The parameters (can be used for general pourposes)
    Parameters mThisParameters;
    
    // Storage for the dof of the node
    Node<3>::DofsContainerType  mDofs;
    
    // I/O information
    char* mFilename;
    std::string mStdStringFilename;
    unsigned int mEchoLevel;
    
    // The framework
    FrameworkEulerLagrange mFramework;
    
    // The member variables related with the MMG library
    MMG5_pMesh mmgMesh;
    MMG5_pSol  mmgSol;
    
    // Where the sub model parts IDs are stored
    boost::unordered_map<int,std::vector<std::string>> mColors;
    
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
     */
    
    void InitializeMeshData()
    {                
        // First we compute the colors
        boost::unordered_map<int,int> NodeColors, CondColors, ElemColors;
        ComputeColors(NodeColors, CondColors, ElemColors);
        
        /////////* MESH FILE */////////
        // Build mesh in MMG5 format //
        
        // Iterate in the nodes
        NodesArrayType& pNode = mrThisModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        // Iterate in the conditions
        ConditionsArrayType& pConditions = mrThisModelPart.Conditions();
        auto numConditions = pConditions.end() - pConditions.begin();
        
        // Iterate in the elements
        ElementsArrayType& pElements = mrThisModelPart.Elements();
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
        // We copy the DOF from the fisrt node (after we release, to avoid problem with previous conditions)
        mDofs = pNode.begin()->GetDofs();
        for (typename Node<3>::DofsContainerType::const_iterator itDoF = mDofs.begin(); itDoF != mDofs.end(); itDoF++)
        {
            itDoF->FreeDof();
        }
        
//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            SetNodes(itNode->X(), itNode->Y(), itNode->Z(), NodeColors[itNode->Id()], i + 1);
            
            bool blocked = false;
            if (itNode->IsDefined(BLOCKED) == true)
            {
                blocked = itNode->Is(BLOCKED);
            }
            if (TDim == 3 && blocked == true)
            {
                BlockNode(i + 1);
            }
            
            // RESETING THE ID OF THE NODES (important for non consecutive meshes)
            itNode->SetId(i + 1);
        }
        
        /* Conditions */
        // We clone the first condition of each type
        if (TDim == 2 && numConditions > 0)
        {
            const CondGeometries2D IndexGeom0 = Line;
            mpRefCondition[IndexGeom0] = pConditions.begin()->Create(0, pConditions.begin()->GetGeometry(), pConditions.begin()->pGetProperties());
        }
        else
        {
            const CondGeometries3D IndexGeom0 = Triangle3D;
            const CondGeometries3D IndexGeom1 = Quadrilateral3D;
            
            for(unsigned int i = 0; i < numConditions; i++) 
            {
                auto itCond = pConditions.begin() + i;

                if ((itCond->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3 && mInitRefCondition[IndexGeom0] == false) // Triangle
                {
                    mpRefCondition[IndexGeom0] = itCond->Create(0, itCond->GetGeometry(), itCond->pGetProperties());
                    mInitRefCondition[IndexGeom0] = true;
                }
                else if ((itCond->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4 && mInitRefCondition[IndexGeom1] == false) // Quadrilateral
                {
                    mpRefCondition[IndexGeom1] = itCond->Create(0, itCond->GetGeometry(), itCond->pGetProperties());
                    mInitRefCondition[IndexGeom1] = true;
                }
                
                if (mInitRefCondition[IndexGeom0] == true && mInitRefCondition[IndexGeom1] == true)
                {
                    break;
                }
            }
            
        }
//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCond = pConditions.begin() + i;
            
            SetConditions(itCond->GetGeometry(), CondColors[itCond->Id()], i + 1);
        }
        
        /* Elements */
        // We clone the first element of each type
        if (TDim == 2 && numElements > 0)
        {
            const ElemGeometries2D IndexGeom0 = Triangle2D;
            mpRefElement[IndexGeom0] = pElements.begin()->Create(0, pElements.begin()->GetGeometry(), pElements.begin()->pGetProperties());
        }
        else
        {
            const ElemGeometries3D IndexGeom0 = Tetrahedra;
            const ElemGeometries3D IndexGeom1 = Prism;
            
            for(unsigned int i = 0; i < numElements; i++) 
            {
                auto itElem = pElements.begin() + i;
                
                if ((itElem->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4 && mInitRefElement[IndexGeom0] == false) // Tetrahedra
                {
                    mpRefElement[IndexGeom0] = itElem->Create(0, itElem->GetGeometry(), itElem->pGetProperties());
                    mInitRefElement[IndexGeom0] = true;
                }
                else if ((itElem->GetGeometry()).GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4 && mInitRefElement[IndexGeom1] == false) // Prism
                {
                    mpRefElement[IndexGeom1] = itElem->Create(0, itElem->GetGeometry(), itElem->pGetProperties());
                    mInitRefElement[IndexGeom1] = true;
                }
                
                if (mInitRefElement[IndexGeom0] == true && mInitRefElement[IndexGeom1] == true)
                {
                    break;
                }
            }
            
        }
//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numElements; i++) 
        {
            auto itElem = pElements.begin() + i;
            
            SetElements(itElem->GetGeometry(), ElemColors[itElem->Id()], i + 1);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * We initialize the metrics of the MMG sol using a level set approach
     */
    
    void InitializeSolData()
    {
        ////////* SOLUTION FILE *////////
        
        // Iterate in the nodes
        NodesArrayType& pNode = mrThisModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        SetSolSizeTensor(numNodes);

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
            const Vector& Metric = itNode->GetValue(MMG_METRIC);
            
            #ifdef KRATOS_DEBUG 
            if(Metric.size() != TDim * 3 - 3) 
            {
                KRATOS_ERROR << "Wrong size of vector MMG_METRIC found for node " << itNode->Id() << " size is " << Metric.size() << " expected size was " << TDim * 3 - 3;
            }
            #endif
            
            // We set the metric
            SetMetricTensor(Metric, i + 1);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * We execute the MMg library and build the new model part from the old model part
     */
    
    void ExecuteRemeshing()
    {
        // Getting the parameters
        const bool SaveToFile = mThisParameters["save_external_files"].GetBool();
        
        // We initialize some values
        const unsigned int StepDataSize = mrThisModelPart.GetNodalSolutionStepDataSize();
        const unsigned int BufferSize   = mrThisModelPart.NodesBegin()->GetBufferSize();
        
        mThisParameters["step_data_size"].SetInt(StepDataSize);
        mThisParameters["buffer_size"].SetInt(BufferSize);
        
        if (mEchoLevel > 0)
        {        
           std::cout << "Step data size: " << StepDataSize << " Buffer size: " << BufferSize << std::endl; 
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
        for (NodeConstantIterator itNode = mrThisModelPart.NodesBegin(); itNode != mrThisModelPart.NodesEnd(); itNode++)
        {
            itNode->Set(TO_ERASE, true);
            rOldModelPart.AddNode(*(itNode.base()));
        }
        mrThisModelPart.RemoveNodesFromAllLevels(TO_ERASE);  
        
        for (ConditionConstantIterator itCond = mrThisModelPart.ConditionsBegin(); itCond != mrThisModelPart.ConditionsEnd(); itCond++)
        {
            itCond->Set(TO_ERASE, true);
        }
        mrThisModelPart.RemoveConditionsFromAllLevels(TO_ERASE); 
        
        for (ElementConstantIterator itElem = mrThisModelPart.ElementsBegin(); itElem != mrThisModelPart.ElementsEnd(); itElem++)
        {
            itElem->Set(TO_ERASE, true);
            rOldModelPart.AddElement(*(itElem.base()));
        }
        mrThisModelPart.RemoveElementsFromAllLevels(TO_ERASE);  
        
        // Create a new model part
        /* NODES */
        for (int unsigned iNode = 1; iNode <= nNodes; iNode++)
        {
            int ref, isRequired;
            NodeType::Pointer pNode = CreateNode(iNode, ref, isRequired);
            
            // Set the DOFs in the nodes 
            for (typename Node<3>::DofsContainerType::const_iterator itDoF = mDofs.begin(); itDoF != mDofs.end(); itDoF++)
            {
                pNode->pAddDof(*itDoF);
            }
            
            if (ref != 0) // NOTE: ref == 0 is the MainModelPart
            {
                std::vector<std::string> ColorList = mColors[ref];
                for (unsigned int colors = 0; colors < ColorList.size(); colors++)
                {
                    std::string SubModelPartName = ColorList[colors];
                    ModelPart& SubModelPart = mrThisModelPart.GetSubModelPart(SubModelPartName);
                    SubModelPart.AddNode(pNode);
                }
            }
        }
        
        /* CONDITIONS */
        unsigned int CondId = 1;
        if (mpRefCondition[0] != nullptr)
        {
            unsigned int CounterCond0 = 0;
            const std::vector<unsigned int> ConditionToRemove0 = CheckConditions0();
            int PropId, isRequired;
            for (int unsigned iCond = 1; iCond <= nConditions[0]; iCond++)
            {
                bool SkipCreation = false;
                if (CounterCond0 < ConditionToRemove0.size())
                {
                    if (ConditionToRemove0[CounterCond0] == iCond)
                    {
                        SkipCreation = true;
                        CounterCond0 += 1;
                    }
                }
                ConditionType::Pointer pCondition = CreateCondition0(CondId, PropId, isRequired, SkipCreation);
                
                if (pCondition != nullptr)
                {
                    mrThisModelPart.AddCondition(pCondition);
                                        
                    if (PropId != 0) // NOTE: PropId == 0 is the MainModelPart
                    {
                        std::vector<std::string> ColorList = mColors[PropId];
                        for (unsigned int iColors = 0; iColors < ColorList.size(); iColors++)
                        {
                            std::string SubModelPartName = ColorList[iColors];
                            ModelPart& SubModelPart = mrThisModelPart.GetSubModelPart(SubModelPartName);
                            SubModelPart.AddCondition(pCondition);
                        }
                    }
                    
                    CondId += 1;
                }
            }
        }
        if (TDim == 3)
        {
            if (mpRefCondition[1] != nullptr) // Quadrilateral
            {
                unsigned int CounterCond1 = 0;
                const std::vector<unsigned int> ConditionToRemove1 = CheckConditions1();
                int PropId, isRequired;
                for (int unsigned iCond = 1; iCond <= nConditions[1]; iCond++)
                {                    
                    bool SkipCreation = false;
                    if (CounterCond1 < ConditionToRemove1.size())
                    {
                        if (ConditionToRemove1[CounterCond1] == iCond)
                        {
                            SkipCreation = true;
                            CounterCond1 += 1;
                        }
                    }
                    ConditionType::Pointer pCondition = CreateCondition1(CondId, PropId, isRequired, SkipCreation);
                    
                    if (pCondition != nullptr)
                    {
                        mrThisModelPart.AddCondition(pCondition);
                                            
                        if (PropId != 0) // NOTE: PropId == 0 is the MainModelPart
                        {
                            std::vector<std::string> ColorList = mColors[PropId];
                            for (unsigned int iColors = 0; iColors < ColorList.size(); iColors++)
                            {
                                std::string SubModelPartName = ColorList[iColors];
                                ModelPart& SubModelPart = mrThisModelPart.GetSubModelPart(SubModelPartName);
                                SubModelPart.AddCondition(pCondition);
                            }
                        }
                        
                        CondId += 1;
                    }
                }
            }
        }
        
        /* ELEMENTS */
        unsigned int ElemId = 1;
        if (mpRefElement[0] != nullptr)
        {
            unsigned int CounterElem0 = 0;
            const std::vector<unsigned int> ElementsToRemove0 = CheckElements0();
            int PropId, isRequired;
            for (int unsigned i_elem = 1; i_elem <= nElements[0]; i_elem++)
            {  
                bool SkipCreation = false;
                if (CounterElem0 < ElementsToRemove0.size())
                {
                    if (ElementsToRemove0[CounterElem0] == i_elem)
                    {
                        SkipCreation = true;
                        CounterElem0 += 1;
                    }
                }
                ElementType::Pointer pElement = CreateElement0(ElemId, PropId, isRequired, SkipCreation);
                
                if (pElement != nullptr)
                {
                    mrThisModelPart.AddElement(pElement);
                    
                    if (PropId != 0) // NOTE: PropId == 0 is the MainModelPart
                    {
                        std::vector<std::string> ColorList = mColors[PropId];
                        for (unsigned int iColors = 0; iColors < ColorList.size(); iColors++)
                        {
                            std::string SubModelPartName = ColorList[iColors];
                            ModelPart& SubModelPart = mrThisModelPart.GetSubModelPart(SubModelPartName);
                            SubModelPart.AddElement(pElement);
                        }
                    }
                    
                    ElemId += 1;
                }
            }
        }
        if (TDim == 3)
        {
            if (mpRefElement[1] != nullptr) // Prism
            {
                unsigned int CounterElem1 = 0;
                const std::vector<unsigned int> ElementsToRemove1 = CheckElements1();
                int PropId, isRequired;
                for (int unsigned i_elem = 1; i_elem <= nElements[1]; i_elem++)
                {
                    bool SkipCreation = false;  
                    if (CounterElem1 < ElementsToRemove1.size())
                    {
                        if (ElementsToRemove1[CounterElem1] == i_elem)
                        {
                            SkipCreation = true;
                            CounterElem1 += 1;
                        }
                    }
                    ElementType::Pointer pElement = CreateElement1(ElemId, PropId, isRequired,SkipCreation);
                    
                    if (pElement != nullptr)
                    {
                        mrThisModelPart.AddElement(pElement);
                        
                        if (PropId != 0) // NOTE: PropId == 0 is the MainModelPart
                        {
                            std::vector<std::string> ColorList = mColors[PropId];
                            for (unsigned int iColors = 0; iColors < ColorList.size(); iColors++)
                            {
                                std::string SubModelPartName = ColorList[iColors];
                                ModelPart& SubModelPart = mrThisModelPart.GetSubModelPart(SubModelPartName);
                                SubModelPart.AddElement(pElement);
                            }
                        }
                        
                        ElemId += 1;
                    }
                }
            }
        }
        
        //  Get the list of submodelparts names
        const std::vector<std::string> SubModelPartNames = mrThisModelPart.GetSubModelPartNames();
       
        // Add the nodes to the differents submodelparts
        for (unsigned int iModelPart = 0; iModelPart < mrThisModelPart.NumberOfSubModelParts(); iModelPart++)
        {
            ModelPart& rSubModelPart = mrThisModelPart.GetSubModelPart(SubModelPartNames[iModelPart]);
           
            std::set<int> AuxSet;
           
            for (ElementConstantIterator itElem = rSubModelPart.ElementsBegin(); itElem != rSubModelPart.ElementsEnd(); itElem++)
            {
                for (unsigned int iNode = 0; iNode < itElem->GetGeometry().size(); iNode++)
                {
                    AuxSet.insert(itElem->GetGeometry()[iNode].Id());
                }
            }
           
            for (ConditionConstantIterator itCond = rSubModelPart.ConditionsBegin(); itCond != rSubModelPart.ConditionsEnd(); itCond++)
            {
                for (unsigned int iNode = 0; iNode < itCond->GetGeometry().size(); iNode++)
                {
                    AuxSet.insert(itCond->GetGeometry()[iNode].Id());
                }
            }
           
            // Clean duplicated nodes
            std::vector<IndexType> NodesIds;
            for( auto it = AuxSet.begin(); it != AuxSet.end(); ++it ) 
            {
                NodesIds.push_back(*it);
            }
           
            rSubModelPart.AddNodes(NodesIds);
        }
        
        /* Save to file */
        if (SaveToFile == true)
        {
            SaveSolutionToFile(true);
        }
       
        /* Free memory */
        FreeMemory();
        
        /* After that we reorder nodes, conditions and elements: */
        ReorderAllIds();
        
        /* We interpolate all the values */
        Parameters InterpolateParameters = Parameters(R"({"echo_level": 1, "framework": "Eulerian", "max_number_of_searchs": 1000, "step_data_size": 0, "buffer_size": 0})" );
        InterpolateParameters["echo_level"].SetInt(mThisParameters["echo_level"].GetInt());
        InterpolateParameters["framework"].SetString(mThisParameters["framework"].GetString());
        InterpolateParameters["max_number_of_searchs"].SetInt(mThisParameters["max_number_of_searchs"].GetInt());
        InterpolateParameters["step_data_size"].SetInt(mThisParameters["step_data_size"].GetInt());
        InterpolateParameters["buffer_size"].SetInt(mThisParameters["buffer_size"].GetInt());
        NodalValuesInterpolationProcess<TDim> InterpolateNodalValues = NodalValuesInterpolationProcess<TDim>(rOldModelPart, mrThisModelPart, InterpolateParameters);
        InterpolateNodalValues.Execute();
        
        /* We initialize elements and conditions */
        InitializeElementsAndConditions();
        
        /* We interpolate the internal variables */
        if (mFramework == Lagrangian) 
        {
            InternalVariablesInterpolationProcess InternalVariablesInterpolation = InternalVariablesInterpolationProcess(rOldModelPart, mrThisModelPart, mThisParameters["internal_variables_parameters"]);
            InternalVariablesInterpolation.Execute();
        }
    }
    
    /**
     * This function reorder the nodes, conditions and elements to avoid problems with non-consecutive ids
     */
    
    void ReorderAllIds()
    {
        NodesArrayType& pNode = mrThisModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();

        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            itNode->SetId(i + 1);
        }

        ConditionsArrayType& pCondition = mrThisModelPart.Conditions();
        auto numConditions = pCondition.end() - pCondition.begin();
        
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCondition = pCondition.begin() + i;
            itCondition->SetId(i + 1);
        }

        ElementsArrayType& pElement = mrThisModelPart.Elements();
        auto numElements = pElement.end() - pElement.begin();

        for(unsigned int i = 0; i < numElements; i++) 
        {
            auto itElement = pElement.begin() + i;
            itElement->SetId(i + 1);
        }
    }
    
    /**
     * After we have transfer the information from the previous modelpart we initilize the elements and conditions
     */
    
    void InitializeElementsAndConditions()
    {
        ConditionsArrayType& pCondition = mrThisModelPart.Conditions();
        auto numConditions = pCondition.end() - pCondition.begin();
        
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCondition = pCondition.begin() + i;
            itCondition->Initialize();
        }

        ElementsArrayType& pElement = mrThisModelPart.Elements();
        auto numElements = pElement.end() - pElement.begin();

        for(unsigned int i = 0; i < numElements; i++) 
        {
            auto itElement = pElement.begin() + i;
            itElement->Initialize();
        }
    }
    
    /**
     * It checks if the nodes are repeated and remove the repeated ones
     */
    
    std::vector<unsigned int> CheckNodes()
    {
        typedef boost::unordered_map<vector<double>, unsigned int, KeyHasherVector<double>, KeyComparorVector<double> > HashMap;
        HashMap NodeMap;
        
        std::vector<unsigned int> NodesToRemoveIds;
        
        vector<double> Coords(TDim);
        
        NodesArrayType& pNode = mrThisModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            const array_1d<double, 3> Coordinates = itNode->Coordinates();
            
            for(unsigned int iCoord = 0; iCoord < TDim; iCoord++)
            {
                Coords[iCoord] = Coordinates[iCoord];
            }
            
            NodeMap[Coords] += 1;
            
            if (NodeMap[Coords] > 1)
            {
                NodesToRemoveIds.push_back(itNode->Id());
                if (mEchoLevel > 0)
                {
                    std::cout << "The mode " << itNode->Id() <<  " is repeated"<< std::endl;
                }
            }
        }
        
       return NodesToRemoveIds;
    }
    
    /**
     * It checks if the conditions are repeated and remove the repeated ones
     */
    
    std::vector<unsigned int> CheckConditions0();
    std::vector<unsigned int> CheckConditions1();
    
    /**
     * It checks if the elemenst are removed and remove the repeated ones
     */
    
    std::vector<unsigned int> CheckElements0();
    std::vector<unsigned int> CheckElements1();
    
    /**
     * It blocks certain nodes before remesh the model
     * @param iNode: The index of the noode
     */
    
    void BlockNode(unsigned int iNode);
    
    /**
     * It creates the new node
     * @param iNode: The index of the new noode
     * @param ref: The submodelpart id
     * @param isRequired: MMG value (I don't know that it does)
     * @return pNode: The pointer to the new node created
     */
    
    NodeType::Pointer CreateNode(
        unsigned int iNode,
        int& ref, 
        int& isRequired
        );
    
    /**
     * It creates the new condition
     * @param CondId: The id of the condition
     * @param PropId: The submodelpart id
     * @param isRequired: MMG value (I don't know that it does)
     * @return pCondition: The pointer to the new condition created
     */
    
    ConditionType::Pointer CreateCondition0(
        const unsigned int CondId,
        int& PropId, 
        int& isRequired,
        bool SkipCreation
        );
    
    ConditionType::Pointer CreateCondition1(
        const unsigned int CondId,
        int& PropId, 
        int& isRequired,
        bool SkipCreation
        );
    
    /**
     * It creates the new element
     * @param CondId: The id of the element
     * @param PropId: The submodelpart id
     * @param isRequired: MMG value (I don't know that it does)
     * @return pElement: The pointer to the new condition created
     */
    
    ElementType::Pointer CreateElement0(
        const unsigned int ElemId,
        int& PropId, 
        int& isRequired,
        bool SkipCreation
        );
    
    ElementType::Pointer CreateElement1(
        const unsigned int ElemId,
        int& PropId, 
        int& isRequired,
        bool SkipCreation
        );
    
    /**
     * It saves the solution and mesh to files (for debugging pourpose g.e)
     * @param PostOutput: If the file to save is after or before remeshing
     */
    
    void SaveSolutionToFile(const bool PostOutput)
    {
        /* GET RESULTS */

        const unsigned int step = mrThisModelPart.GetProcessInfo()[TIME_STEPS];
        
        // Automatically save the mesh 
        OutputMesh(PostOutput, step);

        // Automatically save the solution 
        OutputSol(PostOutput, step);
    }
    
    /**
     * It frees the memory used during all the process
     */
    
    void FreeMemory()
    {
        // Free the MMG structures 
        FreeAll();

        // Free filename (NOTE: Problems with more that one iteration)
//         free(mFilename);
//         mFilename = NULL;
       
        // Free reference std::vectors
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
    }
    
    /** 
     * Initialisation of mesh and sol structures 
     * args of InitMesh:
     * MMG5_ARG_start: we start to give the args of a variadic func
     * MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
     * &mmgMesh: pointer toward your MMG5_pMesh (that store your mesh)
     * MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
     * &mmgSol: pointer toward your MMG5_pSol (that store your metric) 
     */
    
    void InitMesh();
    
    /** 
     * Here the verbosity is set 
     */
    
    void InitVerbosity()
    {
       /* We set the MMG verbosity */
       int verbosityMMG;
       if (mEchoLevel == 0)
       {
           verbosityMMG = 0;
       }
       else if (mEchoLevel == 1)
       {
           verbosityMMG = 0; // NOTE: This way just the essential info from MMG will be printed, but the custom message will appear
       }
       else if (mEchoLevel == 2)
       {
           verbosityMMG = 3;
       }
       else if (mEchoLevel == 3)
       {
           verbosityMMG = 5;
       }
       else if (mEchoLevel == 4)
       {
           std::cout << "       `@@@@'  .@@@@+     :@`    @'  `@@@@@  :@@@@@+     @@@@@@@@ @,    @:  @@@@@@      `@@    @@,  @@@@@@   @@@@@   @     @       '@    `@:        "<< std::endl;
           std::cout << "       `@   @' .@   @'    @;@    @'  @,   ,  :@             @:    @,    @:  @`          `@@    @@,  @       @+   `   @     @       '@    @;@        "<< std::endl;
           std::cout << "       `@   '@ .@   +@    @ @    @'  @       :@             @:    @,    @:  @`          `@+;  ,#@,  @       @        @     @       '@    @ @        "<< std::endl;
           std::cout << "       `@   #@ .@   @@   ## @;   @'  @:      :@             @:    @,    @:  @`          `@ @  @ @,  @       @#       @     @       '@   ;@ ##       "<< std::endl;
           std::cout << "       `@``:@, .@::#@    @  `@   @'  .@@@#   :@@@@@.        @:    @@@@@@@:  @@@@@@      `@ @  @ @,  @@@@@@   @@@@`   @@@@@@@       '@   @`  @       "<< std::endl;
           std::cout << "       `@###.  .@'+@#   `@   @   @'     ,@@  :@             @:    @,    @:  @`          `@ :#+' @,  @          .#@`  @     @  ;;;. '@   @   @.      "<< std::endl;
           std::cout << "       `@      .@   @:  @@@@@@#  @'       @, :@             @:    @,    @:  @`          `@  @@  @,  @            @@  @     @  ''', '@  #@@@@@@      "<< std::endl;
           std::cout << "       `@      .@   `@  @     @  @'       @, :@             @:    @,    @:  @`          `@  ##  @,  @            @@  @     @       '@  @     @      "<< std::endl;
           std::cout << "       `@      .@    @;;@     @, @'  @,``#@  :@,,,,,        @:    @,    @:  @:,,,,      `@      @,  @,,,,,  @:``;@`  @     @       '@ .@     @'     "<< std::endl;
           std::cout << "       `#      .#    `##,     ;# #;  ,###+   ,######        #,    #,    #,  ######      `#      #.  ######  .####    #     #       ;# #'     ,#     "<< std::endl;
           std::cout << "                                                                                                                                                    "<< std::endl;
           std::cout << "                                                                                                                                                    "<< std::endl;
           std::cout << "                                                                                                                                                    "<< std::endl;
           std::cout << "                                                                                                                                                    "<< std::endl;
           std::cout << "                                                                                                                                                    "<< std::endl;
           std::cout << "                                                                      ,.                                                                            "<< std::endl;
           std::cout << "                                                                     ;;;,                                                                           "<< std::endl;
           std::cout << "                                                                     ;;;:                                                                           "<< std::endl;
           std::cout << "                                                                      ::                                                                            "<< std::endl;
           std::cout << "                                                                     ;  ;                                                                           "<< std::endl;
           std::cout << "                                                                         .                                                                          "<< std::endl;
           std::cout << "                                                                    ;                                                                               "<< std::endl;
           std::cout << "                                                                   `      :                                                                         "<< std::endl;
           std::cout << "                                                                   .       :                                                                        "<< std::endl;
           std::cout << "                                                                  :                                                                                 "<< std::endl;
           std::cout << "                                                                            .                                                                       "<< std::endl;
           std::cout << "                                                                 ;           ;                                                                      "<< std::endl;
           std::cout << "                                                                              `                                                                     "<< std::endl;
           std::cout << "                                                                :             `                                                                     "<< std::endl;
           std::cout << "                                                               .               ;                                                                    "<< std::endl;
           std::cout << "                                                              :;;               ;;`                                                                 "<< std::endl;
           std::cout << "                                                              ;;;,..,,,,,,,,,,,;;;;                                                                 "<< std::endl;
           std::cout << "                                                              :;;              .;;,                                                                 "<< std::endl;
           std::cout << "                                                             ; :`              '`;.                                                                 "<< std::endl;
           std::cout << "                                                            `    '            ;  '                                                                  "<< std::endl;
           std::cout << "                                                            :     `          .   . ;                                                                "<< std::endl;
           std::cout << "                                                           ,      ,         `       .                                                               "<< std::endl;
           std::cout << "                                                           .  :    ;        .                                                                       "<< std::endl;
           std::cout << "                                                          ;   ;            :      :  ;                                                              "<< std::endl;
           std::cout << "                                                              .     ;     ;       '   .                                                             "<< std::endl;
           std::cout << "                                                         ;           ,   ;        :   `                                                             "<< std::endl;
           std::cout << "                                                        `            ` `:              ;                                                            "<< std::endl;
           std::cout << "                                                        :    .       `;;,               .                                                           "<< std::endl;
           std::cout << "                                                      ;;`    ;       ;;;;          .    :;;                                                         "<< std::endl;
           std::cout << "                                                     :;;;    ,      ,.;;.,,`       '    ;;;`                                                        "<< std::endl;
           std::cout << "                                                     .;;:        `:         ,,     ;   `;;;                                                         "<< std::endl;
           std::cout << "             `;+                                      `.  ,`;, .,       ';     :.  :: ;  .,                                      ##:                "<< std::endl;
           std::cout << "             `@@@                                    '`    ;;;.        '''        ;;;;     ,                                    @@@;                "<< std::endl;
           std::cout << "              @@@@                                  . :    ;;;`       .''';       :;;:     `                                   +@@@                 "<< std::endl;
           std::cout << "              #@@@                                  : '    ::.        '':''        ,,:    ` ;                                  @@@@                 "<< std::endl;
           std::cout << "               @@@:                                ;  '   :          `'` ''           .   :                                    @@@.                 "<< std::endl;
           std::cout << "               @@@+                                `  ;  ,           ''  ,':           `  ;  :                                .@@@                  "<< std::endl;
           std::cout << "               #@@+                               '   , .            ''   ''           `  '   .                               .@@@                  "<< std::endl;
           std::cout << "               ,@@:                              `     `             ''   ''            . ;   `                                @@#                  "<< std::endl;
           std::cout << "               #@@                               ;   ;:`             ';   ''             :;:   :                               @@@                  "<< std::endl;
           std::cout << "              ,@@@                              ,   ;;;,             ';   ''             ;;;,                                  @@@#                 "<< std::endl;
           std::cout << " ;           '@@@@                              .   ;;;,             ';  `':             ;;;:   :                              @@@@@`          ',   "<< std::endl;
           std::cout << "+@@'`      :@@@@@@'                           `' `;` ,,,             ''  :'`             :,, ,,  `                            :@@@@@@#.     .+@@@   "<< std::endl;
           std::cout << "+@@@@@@#@@@@@@@@@@@.                         ;;;;       .            ''  ''             :      .:;;;`                        `@@@@@@@@@@@@@@@@@@@   "<< std::endl;
           std::cout << " #@@@@@@@@@@@@@@@@@@.                      ;@;;;`        `           ,' `''            ,        .;;;@@`                     .@@@@@@@@@@@@@@@@@@#    "<< std::endl;
           std::cout << "  ,@@@@@@@@; @@@@@@@@+`                  :@@@;;; .:;,.   `.:`         ''''          `:.    .,;:. ;;;@@@#                  `@@@@@@@@@;`+@@@@@@@,     "<< std::endl;
           std::cout << "    `..`    '@@@@ :@@@@+.              `@@@@@;`,      `,;,;;;         ''''          ;;;,;,`       `:@@@@@:              ,+@@@@: @@@@@               "<< std::endl;
           std::cout << "            @@@@@   +@@@@@@',.      `:@@@@@@.   '        `;;;.        `''          .;;;         :  ` +@@@@@+,`   `.,;+@@@@@@'   @@@@@`              "<< std::endl;
           std::cout << "           #@@@@@     '@@@@@@@@@@@@@@@@@@#, :    .        ,:`  ,.      '.       .,` ,:.        ;    :  '#@@@@@@@@@@@@@@@@@;     :@@@@@`             "<< std::endl;
           std::cout << "     ++:;#@@@;@@,        :+@@@@@@@@@#',          .        :       ,,    ``   .,`   ,          `     `     `:+#@@@@@@@#'.         @@`@@@@+;#@        "<< std::endl;
           std::cout << "    `@@@@@@@; @@              ```          :      '       '          ,,:;;.,`     ,   :       ,      ,                           @@ `@@@@@@@'       "<< std::endl;
           std::cout << "    ,@@@@@#` .@@                                   ,      `            ;;;:      :    '      ;       `                           '@+  '@@@@@+       "<< std::endl;
           std::cout << "     :'':    @@,                          ;        `                 `;,;;      :     `     .         ,                           @@:   `::,        "<< std::endl;
           std::cout << "            @@@                                     ;    :         .:          ;            .         `                           @@@.              "<< std::endl;
           std::cout << "           '@@#                          ;           :   '       ,,      :    ;        :   '           ,                          ,@@@              "<< std::endl;
           std::cout << "           '@@                                           `     ,,        '   ;         '  .            `                           '@@              "<< std::endl;
           std::cout << "            ;                           ;             :      ;.          ;  ;          `  `             ,                           .               "<< std::endl;
           std::cout << "                                      ,:               :;.`;`            :.;           ,;;              ,:.                                         "<< std::endl;
           std::cout << "                                     .;;;``````````````;;;```           :;;:          `;;;```````       ;;;                                         "<< std::endl;
           std::cout << "                                     .;;;              ;;;        ```...:;;;..````     ;;;              ;;;                                         "<< std::endl;
           std::cout << "                                      ::               ;,.               ;;            ;:,              ,:;                                         "<< std::endl;
           std::cout << "                                        ;                ;              . .            . ,              ; `                                         "<< std::endl;
           std::cout << "                                     :  `             '   `             :  `          :   ,                :                                        "<< std::endl;
           std::cout << "                                    `    :                ;            ,   '          .   .            '                                            "<< std::endl;
           std::cout << "                                    ,    .           '                 ,   .         :     :                :                                       "<< std::endl;
           std::cout << "                                   `      ,                '          ,              .     .          '                                             "<< std::endl;
           std::cout << "                                   ,      :         ;                 .     '       ,       ;                ;                                      "<< std::endl;
           std::cout << "                                  .        .       `        '        :      .       ,       `        ;                                              "<< std::endl;
           std::cout << "                                  .        ;       ;                 .             ,         ;       `        ;                                     "<< std::endl;
           std::cout << "                                 .          `     `          '      ;        '     ,         `      :                                               "<< std::endl;
           std::cout << "                                 .          ;     :                 `        .    ,           ;     .          ;                                    "<< std::endl;
           std::cout << "                                ,                .            '    ;              ,                ,                                                "<< std::endl;
           std::cout << "                                `            '   :                 `          '  .             '   ,            :                                   "<< std::endl;
           std::cout << "                               ,                ,              '  ;           .  :                .                                                 "<< std::endl;
           std::cout << "                              ``              '.,                               .               '.:              :`                                 "<< std::endl;
           std::cout << "                            `;;:              ;;;              .;;`           ,;;               ;;;              ;;;                                "<< std::endl;
           std::cout << "                            :;;;:::,,,,,,,,,,,;;;,:::::::::;:;:;;;;:::::::::::;;;::::;:;::::::::;;;,,::::::::::::;;;                                "<< std::endl;
           std::cout << "                             ;;.              ;;,              .;;.           :;;               ;;:              :;;                                "<< std::endl;
           std::cout << "                               ;                ;              .``            ```               . ,              ; `                                "<< std::endl;
           std::cout << "                            ,  `             '   .             ,  '           ;  '             ;   :            .   :                               "<< std::endl;
           std::cout << "                           .    '           .    .            :              `    `                `            :                                   "<< std::endl;
           std::cout << "                           `     `          ,     ;           `    '         ;    :           ;     '          :     :                              "<< std::endl;
           std::cout << "                          ;      ;         :                 '              `      ,         .       `         `      `                             "<< std::endl;
           std::cout << "                                  ,        `       '        `       '       ;      .         ,       :        ;       .                             "<< std::endl;
           std::cout << "                         :        .       '         .       :              `        ;       ;         ,                :                            "<< std::endl;
           std::cout << "                        .          ;     `          ,      :         '     ;                          .      ;                                      "<< std::endl;
           std::cout << "                        .          `     :           ;     `                         ;     '           ;    `           ;                           "<< std::endl;
           std::cout << "                       :            '   ,            `    '           '   '           .   .                 :            `                          "<< std::endl;
           std::cout << "                                     `  `             '                               .   ,             ;  :             ,                          "<< std::endl;
           std::cout << "                    `,:              ;,;               .`;             '`'             ;`;               ,,`              ,,                        "<< std::endl;
           std::cout << "                    ;;;              ;;;               ;;:             ;;;             ;;;              `;;;              ;;;                       "<< std::endl;
           std::cout << "                    ;;;..............;;;..........,,,::;;;:::::::::::::;;;:::::::::::::;;;::::,,,.......,;;;.............,;;;                       "<< std::endl;
           std::cout << "                    :;;              ,;;               ;;,             :;;             ;;;               ;;`              ;;,                       "<< std::endl;
           std::cout << "                    ,                ;  :              . ;             ' '             . '              '  ;             ,  ;                       "<< std::endl;
           std::cout << "                   `   :                `             :   ,                           ;                .   `             `                          "<< std::endl;
           std::cout << "                   ,    ;           '    :            .   .           '   '               '            `    :           ;    ;                      "<< std::endl;
           std::cout << "                  `      `                ;          :     '                         :                ;     .          :                            "<< std::endl;
           std::cout << "                  ,      ,         '       `         .      `        '     '        :      '         :       :         `      ;                     "<< std::endl;
           std::cout << "                 `        '                .        ,       :              `        `                        .        :                             "<< std::endl;
           std::cout << "                 ,         .      '         '       .        ;      '       ;      '        '       ,         :      ;         ;                    "<< std::endl;
           std::cout << "                `          .                 ,     ,                        `     .                ;          .     `                               "<< std::endl;
           std::cout << "                ,           '    '           `     ,          ;    '         ;    .          '    .            :    ,           ;                   "<< std::endl;
           std::cout << "               `             ,                :   ,            ,             .   ;                `            .   ;                                "<< std::endl;
           std::cout << "               ,             `  '              ;  ,            `  '           :               '  ;              : .              ;                  "<< std::endl;
           std::cout << "             ;;,              ;;:              ,::              ;;.           ;:;             .::               ;;;              :;:                "<< std::endl;
           std::cout << "            .;;;::::::::::::::;;;:::::::,,,,,..;;;;.............;;;.....``````;;;`````````````;;;,....,,,,,:::::;;;::::::::::::::;;;                "<< std::endl;
           std::cout << "             ;;;              ;;;              ,;;:             ;;;           ;;;             ;;;               ;;;              ;;;                "<< std::endl;
           std::cout << "             ``                `                +'              `,            `,`             :'+                `                .                 "<< std::endl;
           std::cout << "                                                @@                                            ,@@                                                   "<< std::endl;
           std::cout << "                                                @@                                            `@@                                                   "<< std::endl;
           std::cout << "                                               .@@                                             @@                                                   "<< std::endl;
           std::cout << "                                               :@@                                             @@                                                   "<< std::endl;
           std::cout << "                                               '@@                                             @@                                                   "<< std::endl;
           std::cout << "                                               #@@                                             @@                                                   "<< std::endl;
           std::cout << "                                               @@@                                             @@                                                   "<< std::endl;
           std::cout << "                                               @@@                                             @@                                                   "<< std::endl;
           std::cout << "                                               @@#                                             @@                                                   "<< std::endl;
           std::cout << "                                               @@+                                             @@`                                                  "<< std::endl;
           std::cout << "                                               @@+                                             @@,                                                  "<< std::endl;
           std::cout << "                                               @@'                                             @@'                                                  "<< std::endl;
           std::cout << "                                               @@;                                             @@#                                                  "<< std::endl;
           std::cout << "                                               @@:                                             @@@                                                  "<< std::endl;
           std::cout << "                                               @@:                                             @@@                                                  "<< std::endl;
           std::cout << "                                              `@@,                                             @@@                                                  "<< std::endl;
           std::cout << "                                              `@@,                                             @@@                                                  "<< std::endl;
           std::cout << "                                              .@@,                                             @@@                                                  "<< std::endl;
           std::cout << "                                              `@@,                                             @@@                                                  "<< std::endl;
           std::cout << "                                               @@,                                             @@@                                                  "<< std::endl;
           std::cout << "                                               @@:                                             @@@                                                  "<< std::endl;
           std::cout << "                                               @@:                                             @@@                                                  "<< std::endl;
           std::cout << "                                               @@'                                             @@#                                                  "<< std::endl;
           std::cout << "                                               @@#                                             @@'                                                  "<< std::endl;
           std::cout << "                                               @@@                                             @@:                                                  "<< std::endl;
           std::cout << "                                               @@@                                            `@@`                                                  "<< std::endl;
           std::cout << "                                               @@@                                            #@@                                                   "<< std::endl;
           std::cout << "                                               :@@@                                           @@@                                                   "<< std::endl;
           std::cout << "                                                @@@,                                         '@@@                                                   "<< std::endl;
           std::cout << "                                                 @@@                                        `@@@                                                    "<< std::endl;
           std::cout << "                                                 `@@@                                       @@@                                                     "<< std::endl;
           std::cout << "                                                  '@@@                                     @@@'                                                     "<< std::endl;
           std::cout << "                                                   #@@@                                   @@@#                                                      "<< std::endl;
           std::cout << "                                                    #@@@:                                @@@@                                                       "<< std::endl;
           std::cout << "                                                     #@@@@`                            +@@@@                                                        "<< std::endl;
           std::cout << "                                                      '@@@@@                         :@@@@#                                                         "<< std::endl;
           std::cout << "                                                       ,@@@@@#                     ,@@@@@'                                                          "<< std::endl;
           std::cout << "                                                         @@@@@@#                 ,@@@@@@.                                                           "<< std::endl;
           std::cout << "                                                          `@@@@@@#             ,@@@@@@+                                                             "<< std::endl;
           std::cout << "                                                             +@@@@            #@@@@@:                                                               "<< std::endl;
           std::cout << "                                                              @@@@            @@@@'                                                                 "<< std::endl;
           std::cout << "                                                              @@@@            +@@@                                                                  "<< std::endl;
           std::cout << "                                                              @@@@            .@@@                                                                  "<< std::endl;
           std::cout << "                                                              @@@;             @@@                                                                  "<< std::endl;
           std::cout << "                                                              @@@              @@@                                                                  "<< std::endl;
           std::cout << "                                                              @@@              @@@                                                                  "<< std::endl;
           std::cout << "                                                              @@@              ,@@.                                                                 "<< std::endl;
           std::cout << "                                                              @@+               @@,                                                                 "<< std::endl;
           std::cout << "                                                              @@.               @@:                                                                 "<< std::endl;
           std::cout << "                                                              @@                @@:                                                                 "<< std::endl;
           std::cout << "                                                              @@                @@:                                                                 "<< std::endl;
           std::cout << "                                                              @@                '@:                                                                 "<< std::endl;
           std::cout << "                                                              @@                `@,                                                                 "<< std::endl;
           std::cout << "                                                              @,                 @.                                                                 "<< std::endl;
           std::cout << "                                                              @                  @                                                                  "<< std::endl;
           std::cout << "                                                                                 .                                                                  "<< std::endl;
           std::cout << "                                                                                                                                                    "<< std::endl;
           std::cout << "                                                                                                                                                    "<< std::endl;
           verbosityMMG = 5;
       }
       else
       {
           verbosityMMG = 10;
       }
       
       InitVerbosityParameter(verbosityMMG);
    }
    
    /** 
     * Here the verbosity is set using the API
     * @param verbosityMMG: The equivalent verbosity level in the MMG API
     */
        
    void InitVerbosityParameter(int verbosityMMG);
    
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
     * This sets the size of the solution for the scalar case
     * @param numNodes: Number of nodes
     */
    
    void SetSolSizeScalar(const int numNodes);
    
    /**
     * This sets the size of the solution for the vector case
     * @param numNodes: Number of nodes
     */
    
    void SetSolSizeVector(const int numNodes);
    
    /**
     * This sets the size of the solution for the tensor case
     * @param numNodes: Number of nodes
     */
    
    void SetSolSizeTensor(const int numNodes);
    
    /**
     * This checks the mesh data and prints if it is OK
     */
    
    void CheckMeshData();
    
    /**
     * This sets the output mesh
     */
    
    void OutputMesh(
        const bool PostOutput, 
        const unsigned int step
        );
    
    /**
     * This sets the output sol
     */
    
    void OutputSol(
        const bool PostOutput, 
        const unsigned int step
        );
    
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
     * @param NodeColors: Map where the submodelparts and nodes are stored
     * @param CondColors: Map where the submodelparts and conditions are stored
     * @param ElemColors: Map where the submodelparts and elements are stored
     */
    
    void ComputeColors(
        boost::unordered_map<int,int>& NodeColors,
        boost::unordered_map<int,int>& CondColors,
        boost::unordered_map<int,int>& ElemColors
        )
    {        
        // Initialize and create the auxiliar maps
        const std::vector<std::string> SubModelPartNames = mrThisModelPart.GetSubModelPartNames();
        boost::unordered_map<int,std::set<int>> AuxNodeColors, AuxCondColors, AuxElemColors;
        
        std::vector<std::string> ModelPartNames;
        ModelPartNames.push_back(mrThisModelPart.Name());
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
                ModelPart& rSubModelPart = mrThisModelPart.GetSubModelPart(ModelPartNames[i_sub]);
                
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
                    AuxNodeColors[itNode->Id()].insert(color);
                }
                
                /* Conditions */
        //         #pragma omp parallel for 
                for(unsigned int i = 0; i < numConditions; i++) 
                {
                    auto itCond = pConditions.begin() + i;
                    AuxCondColors[itCond->Id()].insert(color);
                }
                
                /* Elements */
        //         #pragma omp parallel for 
                for(unsigned int i = 0; i < numElements; i++) 
                {
                    auto itElem = pElements.begin() + i;
                    AuxElemColors[itElem->Id()].insert(color);
                }
            }
            
            color += 1;
        }
        
        // The iterator for the auxiliar maps is created
        typedef boost::unordered_map<int,std::set<int>>::iterator itType;
        
        // Now detect all the cases in which a node or a cond belongs to more than one part simultaneously 
        boost::unordered_map<std::set<int>, int> Combinations;
        
        /* Nodes */
        for(itType iterator = AuxNodeColors.begin(); iterator != AuxNodeColors.end(); iterator++) 
        {
//             const int key = iterator->first;
            const std::set<int> Value = iterator->second;
            
            if (Value.size() > 1)
            {
                Combinations[Value] = -1;
            }
        }
        
        /* Conditions */
        for(itType iterator = AuxCondColors.begin(); iterator != AuxCondColors.end(); iterator++) 
        {
//             const int key = iterator->first;
            const std::set<int> Value = iterator->second;
            
            if (Value.size() > 1)
            {
                Combinations[Value] = -1;
            }
        }

        /* Elements */
        for(itType iterator = AuxElemColors.begin(); iterator != AuxElemColors.end(); iterator++) 
        {
//             const int key = iterator->first;
            const std::set<int> Value = iterator->second;
            
            if (Value.size() > 1)
            {
                Combinations[Value] = -1;
            }
        }
        
        /* Combinations */
        typedef boost::unordered_map<std::set<int>,int>::iterator CombType;
        for(CombType iterator = Combinations.begin(); iterator != Combinations.end(); iterator++) 
        {
            const std::set<int> key = iterator->first;
//             const int Value = iterator->second;
            
            for( auto it = key.begin(); it != key.end(); ++it ) 
            {
                mColors[color].push_back(mColors[*it][0]);
            }
            Combinations[key] = color;
            color += 1;
            
        }
        
        // The final maps are created
        /* Nodes */
        for(itType iterator = AuxNodeColors.begin(); iterator != AuxNodeColors.end(); iterator++) 
        {
            const int key = iterator->first;
            const std::set<int> Value = iterator->second;
            
            if (Value.size() == 0)
            {
                NodeColors[key] = 0; // Main Model Part
            }
            else if (Value.size() == 1) // Another Model Part
            {
                NodeColors[key] = *Value.begin();
            }
            else // There is a combination
            {
                NodeColors[key] = Combinations[Value];
            }
        }
        
        /* Conditions */
        for(itType iterator = AuxCondColors.begin(); iterator != AuxCondColors.end(); iterator++) 
        {
            const int key = iterator->first;
            const std::set<int> Value = iterator->second;
            
            if (Value.size() == 0)
            {
                CondColors[key] = 0; // Main Model Part
            }
            else if (Value.size() == 1) // Another Model Part
            {
                CondColors[key] = *Value.begin();
            }
            else // There is a combination
            {
                CondColors[key] = Combinations[Value];
            }
        }
        
        /* Elements */
        for(itType iterator = AuxElemColors.begin(); iterator != AuxElemColors.end(); iterator++) 
        {
            const int key = iterator->first;
            const std::set<int> Value = iterator->second;
            
            if (Value.size() == 0)
            {
                ElemColors[key] = 0; // Main Model Part
            }
            else if (Value.size() == 1) // Another Model Part
            {
                ElemColors[key] = *Value.begin();
            }
            else // There is a combination
            {
                ElemColors[key] = Combinations[Value];
            }
        }
    }

    /**
     * This function is used to compute the metric scalar
     * @param Metric: The inverse of the size node
     */

    void SetMetricScalar(
        const double& Metric,
        const int NodeId 
        );
    
    /**
     * This function is used to compute the metric vector (x, y, z)
     * @param Metric: This array contains the components of the metric vector
     */

    void SetMetricVector(
        const array_1d<double, 3>& Metric,
        const int NodeId 
        );
    
    /**
     * This function is used to compute the Hessian metric tensor, note that when using the Hessian, more than one metric can be defined simultaneously, so in consecuence we need to define the elipsoid which defines the volume of maximal intersection
     * @param Metric: This array contains the components of the metric tensor in the MMG defined order
     */

    void SetMetricTensor(
        const Vector& Metric,
        const int NodeId 
        );
    
    /**
     * This converts the framework string to an enum
     * @param str: The string
     * @return FrameworkEulerLagrange: The equivalent enum
     */
        
    FrameworkEulerLagrange ConvertFramework(const std::string& str)
    {
        if(str == "Lagrangian") 
        {
            return Lagrangian;
        }
        else if(str == "Eulerian") 
        {
            return Eulerian;
        }
        else
        {
            return Eulerian;
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

///@name Explicit Specializations
///@{

    template<>  
    std::vector<unsigned int> MmgUtility<2>::CheckConditions0()
    {
        typedef boost::unordered_map<vector<unsigned int>, unsigned int, KeyHasherVector<unsigned int>, KeyComparorVector<unsigned int> > HashMap;
        HashMap edge_map;

        vector<unsigned int> ids(2);

        std::vector<unsigned int> ConditionsToRemove;
        
        // Iterate in the conditions
        for(int i = 0; i < mmgMesh->na; i++) 
        {
            int edge0, edge1, PropId, isRidge, isRequired;
            
            if (MMG2D_Get_edge(mmgMesh, &edge0, &edge1, &PropId, &isRidge, &isRequired) != 1 )
            {
                exit(EXIT_FAILURE);
            }
            
            ids[0] = edge0;
            ids[1] = edge1;

            //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
            std::sort(ids.begin(), ids.end());

            edge_map[ids] += 1;
            
            if (edge_map[ids] > 1)
            {
                ConditionsToRemove.push_back(i + 1);
            }
        }
        
        return ConditionsToRemove;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    std::vector<unsigned int> MmgUtility<3>::CheckConditions0()
    {
        typedef boost::unordered_map<vector<unsigned int>, unsigned int, KeyHasherVector<unsigned int>, KeyComparorVector<unsigned int> > HashMap;
        HashMap TriangleMap;

        vector<unsigned int> IdsTriangles(3);

        std::vector<unsigned int> ConditionsToRemove;
                
        for(int i = 0; i < mmgMesh->nt; i++) 
        {
            int Vertex0, Vertex1, Vertex2, PropId, isRequired;

            if (MMG3D_Get_triangle(mmgMesh, &Vertex0, &Vertex1, &Vertex2, &PropId, &isRequired) != 1 )
            {
                exit(EXIT_FAILURE);
            }

            IdsTriangles[0] = Vertex0;
            IdsTriangles[1] = Vertex1;
            IdsTriangles[2] = Vertex2;
            
            //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
            std::sort(IdsTriangles.begin(), IdsTriangles.end());

            TriangleMap[IdsTriangles] += 1;
            
            if (TriangleMap[IdsTriangles] > 1)
            {
                ConditionsToRemove.push_back(i + 1);
            }
        }
        
        return ConditionsToRemove;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    std::vector<unsigned int> MmgUtility<3>::CheckConditions1()
    {
        typedef boost::unordered_map<vector<unsigned int>, unsigned int, KeyHasherVector<unsigned int>, KeyComparorVector<unsigned int> > HashMap;
        HashMap QuadrilateralMap;

        vector<unsigned int> IdsQuadrilateral(4);

        std::vector<unsigned int> ConditionsToRemove;
                
        for(int i = 0; i < mmgMesh->nquad; i++) 
        {
            int Vertex0, Vertex1, Vertex2, Vertex3, PropId, isRequired;

            if (MMG3D_Get_quadrilateral(mmgMesh, &Vertex0, &Vertex1, &Vertex2, &Vertex3, &PropId, &isRequired) != 1 )
            {
                exit(EXIT_FAILURE);
            }

            IdsQuadrilateral[0] = Vertex0;
            IdsQuadrilateral[1] = Vertex1;
            IdsQuadrilateral[2] = Vertex2;
            IdsQuadrilateral[3] = Vertex3;
            
            //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
            std::sort(IdsQuadrilateral.begin(), IdsQuadrilateral.end());

            QuadrilateralMap[IdsQuadrilateral] += 1;
            
            if (QuadrilateralMap[IdsQuadrilateral] > 1)
            {
                ConditionsToRemove.push_back(i + 1);
            }
        }
        
        return ConditionsToRemove;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    std::vector<unsigned int> MmgUtility<2>::CheckElements0()
    {
        typedef boost::unordered_map<vector<unsigned int>, unsigned int, KeyHasherVector<unsigned int>, KeyComparorVector<unsigned int> > HashMap;
        HashMap TriangleMap;

        vector<unsigned int> IdsTriangles(3);

        std::vector<unsigned int> ElementsToRemove;
        
        // Iterate in the elements
        for(int i = 0; i < mmgMesh->nt; i++) 
        {
            int Vertex0, Vertex1, Vertex2, PropId, isRequired;
            
            if (MMG2D_Get_triangle(mmgMesh, &Vertex0, &Vertex1, &Vertex2, &PropId, &isRequired) != 1 )
            {
                exit(EXIT_FAILURE);
            }
            
            IdsTriangles[0] = Vertex0;
            IdsTriangles[1] = Vertex1;
            IdsTriangles[2] = Vertex2;

            //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
            std::sort(IdsTriangles.begin(), IdsTriangles.end());

            TriangleMap[IdsTriangles] += 1;
            
            if (TriangleMap[IdsTriangles] > 1)
            {
                ElementsToRemove.push_back(i + 1);
            }
        }
        
        return ElementsToRemove;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    std::vector<unsigned int> MmgUtility<3>::CheckElements0()
    {
        typedef boost::unordered_map<vector<unsigned int>, unsigned int, KeyHasherVector<unsigned int>, KeyComparorVector<unsigned int> > HashMap;
        HashMap TriangleMap;

        vector<unsigned int> IdsTetrahedron(4);

        std::vector<unsigned int> ElementsToRemove;
                
        for(int i = 0; i < mmgMesh->ne; i++) 
        {
            int Vertex0, Vertex1, Vertex2, Vertex3, PropId, isRequired;

            if (MMG3D_Get_tetrahedron(mmgMesh, &Vertex0, &Vertex1, &Vertex2, &Vertex3, &PropId, &isRequired) != 1 )
            {
                exit(EXIT_FAILURE);
            }

            IdsTetrahedron[0] = Vertex0;
            IdsTetrahedron[1] = Vertex1;
            IdsTetrahedron[2] = Vertex2;
            IdsTetrahedron[3] = Vertex3;
            
            //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
            std::sort(IdsTetrahedron.begin(), IdsTetrahedron.end());

            TriangleMap[IdsTetrahedron] += 1;
            
            if (TriangleMap[IdsTetrahedron] > 1)
            {
                ElementsToRemove.push_back(i + 1);
            }
        }
        
        return ElementsToRemove;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    std::vector<unsigned int> MmgUtility<3>::CheckElements1()
    {
        typedef boost::unordered_map<vector<unsigned int>, unsigned int, KeyHasherVector<unsigned int>, KeyComparorVector<unsigned int> > HashMap;
        HashMap PrismMap;

        vector<unsigned int> IdsPrisms(6);

        std::vector<unsigned int> ElementsToRemove;
                
        for(int i = 0; i < mmgMesh->nprism; i++) 
        {
            int Vertex0, Vertex1, Vertex2, Vertex3, Vertex4, Vertex5, PropId, isRequired;

            if (MMG3D_Get_prism(mmgMesh, &Vertex0, &Vertex1, &Vertex2, &Vertex3, &Vertex4, &Vertex5, &PropId, &isRequired) != 1 )
            {
                exit(EXIT_FAILURE);
            }

            IdsPrisms[0] = Vertex0;
            IdsPrisms[1] = Vertex1;
            IdsPrisms[2] = Vertex2;
            IdsPrisms[3] = Vertex3;
            IdsPrisms[4] = Vertex4;
            IdsPrisms[5] = Vertex5;
            
            //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
            std::sort(IdsPrisms.begin(), IdsPrisms.end());

            PrismMap[IdsPrisms] += 1;
            
            if (PrismMap[IdsPrisms] > 1)
            {
                ElementsToRemove.push_back(i + 1);
            }
        }
        
        return ElementsToRemove;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
//     template<>  // NOTE: Not yet avalaible in the official API
//     void MmgUtility<2>::BlockNode(unsigned int iNode)
//     {
//         if (MMG2D_Set_requiredVertex(mmgMesh, iNode) != 1 )
//         {
//             exit(EXIT_FAILURE);
//         }
//     }

    /***********************************************************************************/
    /***********************************************************************************/
    

    template<>  
    void MmgUtility<3>::BlockNode(unsigned int iNode)
    {
        if (MMG3D_Set_requiredVertex(mmgMesh, iNode) != 1 )
        {
            exit(EXIT_FAILURE);
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    NodeType::Pointer MmgUtility<2>::CreateNode(
        unsigned int iNode,
        int& ref, 
        int& isRequired
        )
    {
        double Coord0, Coord1;
        int isCorner;
        
        if (MMG2D_Get_vertex(mmgMesh, &Coord0, &Coord1, &ref, &isCorner, &isRequired) != 1 )
        {
            exit(EXIT_FAILURE);
        }
        
        NodeType::Pointer pNode = mrThisModelPart.CreateNewNode(iNode, Coord0, Coord1, 0.0);
        
        return pNode;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    NodeType::Pointer MmgUtility<3>::CreateNode(
        unsigned int iNode,
        int& ref, 
        int& isRequired
        )
    {
        double Coord0, Coord1, Coord2;
        int isCorner;
        
        if (MMG3D_Get_vertex(mmgMesh, &Coord0, &Coord1, &Coord2, &ref, &isCorner, &isRequired) != 1 )
        {
            exit(EXIT_FAILURE);
        }
        
        NodeType::Pointer pNode = mrThisModelPart.CreateNewNode(iNode, Coord0, Coord1, Coord2);
        
        return pNode;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    ConditionType::Pointer MmgUtility<2>::CreateCondition0(        
        const unsigned int CondId,
        int& PropId, 
        int& isRequired, 
        bool SkipCreation
        )
    {
        ConditionType::Pointer pCondition = nullptr;
        
        const CondGeometries2D IndexGeom = Line;
        
        int Edge0, Edge1, isRidge;
        
        if (MMG2D_Get_edge(mmgMesh, &Edge0, &Edge1, &PropId, &isRidge, &isRequired) != 1 )
        {
            exit(EXIT_FAILURE);
        }
        
        // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
        if (Edge0 == 0) SkipCreation = true;
        if (Edge1 == 0) SkipCreation = true;
        
        if (SkipCreation == false)
        {
            std::vector<NodeType::Pointer> ConditionNodes (2);
            ConditionNodes[0] = mrThisModelPart.pGetNode(Edge0);
            ConditionNodes[1] = mrThisModelPart.pGetNode(Edge1);    
            
            pCondition = mpRefCondition[IndexGeom]->Create(CondId, ConditionNodes, mpRefCondition[IndexGeom]->pGetProperties());
        }
        else if (mEchoLevel > 0)
        {
            std::cout << "Condition creation avoided" << std::endl;
        }
        
        return pCondition;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    ConditionType::Pointer MmgUtility<3>::CreateCondition0(
        const unsigned int CondId,
        int& PropId, 
        int& isRequired,
        bool SkipCreation
        )
    {
        ConditionType::Pointer pCondition = nullptr;
        
        const CondGeometries3D IndexGeom = Triangle3D;
        
        int Vertex0, Vertex1, Vertex2;

        if (MMG3D_Get_triangle(mmgMesh, &Vertex0, &Vertex1, &Vertex2, &PropId, &isRequired) != 1 )
        {
            exit(EXIT_FAILURE);
        }
        
        // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
        if (Vertex0 == 0) SkipCreation = true;
        if (Vertex1 == 0) SkipCreation = true;
        if (Vertex2 == 0) SkipCreation = true;
        
        if (SkipCreation == false)
        {
            std::vector<NodeType::Pointer> ConditionNodes (3);
            ConditionNodes[0] = mrThisModelPart.pGetNode(Vertex0);
            ConditionNodes[1] = mrThisModelPart.pGetNode(Vertex1);
            ConditionNodes[2] = mrThisModelPart.pGetNode(Vertex2);
        
            pCondition = mpRefCondition[IndexGeom]->Create(CondId, ConditionNodes, mpRefCondition[IndexGeom]->pGetProperties());
        }
        else if (mEchoLevel > 0)
        {
            std::cout << "Condition creation avoided" << std::endl;
        }
        
        return pCondition;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    ConditionType::Pointer MmgUtility<3>::CreateCondition1(
        const unsigned int CondId,
        int& PropId, 
        int& isRequired,
        bool SkipCreation
        )
    {
        ConditionType::Pointer pCondition = nullptr;
        
        const CondGeometries3D IndexGeom = Quadrilateral3D;
        
        int Vertex0, Vertex1, Vertex2, Vertex3;

        if (MMG3D_Get_quadrilateral(mmgMesh, &Vertex0, &Vertex1, &Vertex2, &Vertex3, &PropId, &isRequired) != 1 )
        {
            exit(EXIT_FAILURE);
        }
        
        // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
        if (Vertex0 == 0) SkipCreation = true;
        if (Vertex1 == 0) SkipCreation = true;
        if (Vertex2 == 0) SkipCreation = true;
        if (Vertex3 == 0) SkipCreation = true;
        
        if (SkipCreation == false)
        {
            std::vector<NodeType::Pointer> ConditionNodes (4);
            ConditionNodes[0] = mrThisModelPart.pGetNode(Vertex0);
            ConditionNodes[1] = mrThisModelPart.pGetNode(Vertex1);
            ConditionNodes[2] = mrThisModelPart.pGetNode(Vertex2);
            ConditionNodes[3] = mrThisModelPart.pGetNode(Vertex3);
            
            pCondition = mpRefCondition[IndexGeom]->Create(CondId, ConditionNodes, mpRefCondition[IndexGeom]->pGetProperties());
        }
        else if (mEchoLevel > 0)
        {
            std::cout << "Condition creation avoided" << std::endl;
        }
        
        return pCondition;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    ElementType::Pointer MmgUtility<2>::CreateElement0(        
        const unsigned int ElemId,
        int& PropId, 
        int& isRequired,
        bool SkipCreation
        )
    {
        ElementType::Pointer pElement = nullptr;
        
        const ElemGeometries2D IndexGeom = Triangle2D;
        
        int Vertex0, Vertex1, Vertex2;
        
        if (MMG2D_Get_triangle(mmgMesh, &Vertex0, &Vertex1, &Vertex2, &PropId, &isRequired) != 1 )
        {
            exit(EXIT_FAILURE);
        }

        // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
        if (Vertex0 == 0) SkipCreation = true;
        if (Vertex1 == 0) SkipCreation = true;
        if (Vertex2 == 0) SkipCreation = true;
        
        if (SkipCreation == false)
        {
            std::vector<NodeType::Pointer> ElementNodes (3);
            ElementNodes[0] = mrThisModelPart.pGetNode(Vertex0);
            ElementNodes[1] = mrThisModelPart.pGetNode(Vertex1);
            ElementNodes[2] = mrThisModelPart.pGetNode(Vertex2);
            
            pElement = mpRefElement[IndexGeom]->Create(ElemId, ElementNodes, mpRefElement[IndexGeom]->pGetProperties());
        }
        else if (mEchoLevel > 0)
        {
            std::cout << "Element creation avoided" << std::endl;
        }
        
        return pElement;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    ElementType::Pointer MmgUtility<3>::CreateElement0(
        const unsigned int ElemId,
        int& PropId, 
        int& isRequired,
        bool SkipCreation
        )
    {
        ElementType::Pointer pElement = nullptr;
        
        const ElemGeometries3D IndexGeom = Tetrahedra;
        
        int Vertex0, Vertex1, Vertex2, Vertex3;
        
        if (MMG3D_Get_tetrahedron(mmgMesh, &Vertex0, &Vertex1, &Vertex2, &Vertex3, &PropId, &isRequired) != 1 )
        {
            exit(EXIT_FAILURE);
        }
        
        // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
        if (Vertex0 == 0) SkipCreation = true;
        if (Vertex1 == 0) SkipCreation = true;
        if (Vertex2 == 0) SkipCreation = true;
        if (Vertex3 == 0) SkipCreation = true;
        
        if (SkipCreation == false)
        {
            std::vector<NodeType::Pointer> ElementNodes (4);
            ElementNodes[0] = mrThisModelPart.pGetNode(Vertex0);
            ElementNodes[1] = mrThisModelPart.pGetNode(Vertex1);
            ElementNodes[2] = mrThisModelPart.pGetNode(Vertex2);
            ElementNodes[3] = mrThisModelPart.pGetNode(Vertex3);
            
            pElement = mpRefElement[IndexGeom]->Create(ElemId, ElementNodes, mpRefElement[IndexGeom]->pGetProperties());
        }
        else if (mEchoLevel > 0)
        {
            std::cout << "Element creation avoided" << std::endl;
        }
        
        return pElement;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    ElementType::Pointer MmgUtility<3>::CreateElement1(
        const unsigned int ElemId,
        int& PropId, 
        int& isRequired,
        bool SkipCreation
        )
    {
        ElementType::Pointer pElement = nullptr;
        
        const ElemGeometries3D IndexGeom = Prism;
                    
        int Vertex0, Vertex1, Vertex2, Vertex3, Vertex4, Vertex5;
        
        if (MMG3D_Get_prism(mmgMesh, &Vertex0, &Vertex1, &Vertex2, &Vertex3, &Vertex4, &Vertex5, &PropId, &isRequired) != 1 )
        {
            exit(EXIT_FAILURE);
        }
        
        // FIXME: This is not the correct solution to the problem, I asked in the MMG Forum
        if (Vertex0 == 0) SkipCreation = true;
        if (Vertex1 == 0) SkipCreation = true;
        if (Vertex2 == 0) SkipCreation = true;
        if (Vertex3 == 0) SkipCreation = true;
        if (Vertex4 == 0) SkipCreation = true;
        if (Vertex5 == 0) SkipCreation = true;
        
        if (SkipCreation == false)
        {
            std::vector<NodeType::Pointer> ElementNodes (6);
            ElementNodes[0] = mrThisModelPart.pGetNode(Vertex0);
            ElementNodes[1] = mrThisModelPart.pGetNode(Vertex1);
            ElementNodes[2] = mrThisModelPart.pGetNode(Vertex2);
            ElementNodes[3] = mrThisModelPart.pGetNode(Vertex3);
            ElementNodes[4] = mrThisModelPart.pGetNode(Vertex4);
            ElementNodes[5] = mrThisModelPart.pGetNode(Vertex5);
        
            pElement = mpRefElement[IndexGeom]->Create(ElemId, ElementNodes, mpRefElement[IndexGeom]->pGetProperties());
        }
        else if (mEchoLevel > 0)
        {
            std::cout << "Element creation avoided" << std::endl;
        }
        
        return pElement;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<2>::InitMesh()
    {  
        mmgMesh = NULL;
        mmgSol = NULL;
       
        // We init the MMG mesh and sol
        MMG2D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end); 
        
        InitVerbosity();
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<3>::InitMesh()
    {   
        mmgMesh = NULL;
        mmgSol = NULL;
        
        MMG3D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol, MMG5_ARG_end); 
        
        InitVerbosity();
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<2>::InitVerbosityParameter(int verbosityMMG)
    {  
       if ( !MMG2D_Set_iparameter(mmgMesh,mmgSol,MMG2D_IPARAM_verbose, verbosityMMG) )
       {
           exit(EXIT_FAILURE);
       }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<3>::InitVerbosityParameter(int verbosityMMG)
    {       
       if ( !MMG3D_Set_iparameter(mmgMesh,mmgSol,MMG3D_IPARAM_verbose, verbosityMMG) )
       {
           exit(EXIT_FAILURE);
       }
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
        //Give the size of the mesh: numNodes Vertex, numElements tetra and prism, numArrayConditions triangles and quadrilaterals, 0 edges (3D) 
        if ( MMG3D_Set_meshSize(mmgMesh, numNodes, numArrayElements[0], numArrayElements[1], numArrayConditions[0], numArrayConditions[1], 0) != 1 ) 
        {
            exit(EXIT_FAILURE);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    void MmgUtility<2>::SetSolSizeScalar(const int numNodes)
    {
        if ( MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,numNodes,MMG5_Scalar) != 1 )
        {
            exit(EXIT_FAILURE);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    void MmgUtility<3>::SetSolSizeScalar(const int numNodes)
    {
        if ( MMG3D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,numNodes,MMG5_Scalar) != 1 )
        {
            exit(EXIT_FAILURE);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    void MmgUtility<2>::SetSolSizeVector(const int numNodes)
    {
        KRATOS_ERROR << "WARNING:: Vector metric not avalaible in 2D" << std::endl;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    void MmgUtility<3>::SetSolSizeVector(const int numNodes)
    {
        if ( MMG3D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,numNodes,MMG5_Vector) != 1 )
        {
            exit(EXIT_FAILURE);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    void MmgUtility<2>::SetSolSizeTensor(const int numNodes)
    {
        if ( MMG2D_Set_solSize(mmgMesh,mmgSol,MMG5_Vertex,numNodes,MMG5_Tensor) != 1 )
        {
            exit(EXIT_FAILURE);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    void MmgUtility<3>::SetSolSizeTensor(const int numNodes)
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
    void MmgUtility<2>::OutputMesh(
        const bool PostOutput,
        const unsigned int step
        )
    {
        std::string MeshName;
        if (PostOutput == true)
        {
            MeshName = mStdStringFilename+"_step="+std::to_string(step)+".o.mesh";
        }
        else
        {
            MeshName = mStdStringFilename+"_step="+std::to_string(step)+".mesh";
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
    void MmgUtility<3>::OutputMesh(
        const bool PostOutput,
        const unsigned int step
        )
    {
        std::string MeshName;
        if (PostOutput == true)
        {
            MeshName = mStdStringFilename+"_step="+std::to_string(step)+".o.mesh";
        }
        else
        {
            MeshName = mStdStringFilename+"_step="+std::to_string(step)+".mesh";
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
    void MmgUtility<2>::OutputSol(
        const bool PostOutput,
        const unsigned int step
        )
    {
        std::string SolName;
        if (PostOutput == true)
        {
            SolName = mStdStringFilename+"_step="+std::to_string(step)+".o.sol";
        }
        else
        {
            SolName = mStdStringFilename+"_step="+std::to_string(step)+".sol";
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
    void MmgUtility<3>::OutputSol(
        const bool PostOutput,
        const unsigned int step
        )
    {
        std::string SolName;
        if (PostOutput == true)
        {
            SolName = mStdStringFilename+"_step="+std::to_string(step)+".o.sol";
        }
        else
        {
            SolName = mStdStringFilename+"_step="+std::to_string(step)+".sol";
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
            KRATOS_ERROR << "WARNING: BAD ENDING OF MMG2DLIB: UNABLE TO SAVE MESH. ier: " << ier << std::endl;
        }
        else if ( ier == MMG5_LOWFAILURE )
        {
            KRATOS_ERROR << "WARNING: BAD ENDING OF MMG2DLIB. ier: " << ier << std::endl;
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
            KRATOS_ERROR << "WARNING: BAD ENDING OF MMG3DLIB: UNABLE TO SAVE MESH. ier: " << ier << std::endl;
        }
        else if ( ier == MMG5_LOWFAILURE )
        {
            KRATOS_ERROR << "WARNING: BAD ENDING OF MMG3DLIB. ier: " << ier << std::endl;
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
        
        // Set fixed boundary
        bool blocked1 = false;
        if (Geom[0].IsDefined(BLOCKED) == true)
        {
            blocked1 = Geom[0].Is(BLOCKED);
        }
        bool blocked2 = false;
        if (Geom[1].IsDefined(BLOCKED) == true)
        {
            blocked2 = Geom[1].Is(BLOCKED);
        }

        if ((blocked1 && blocked2) == true)
        {
            if ( MMG2D_Set_requiredEdge(mmgMesh, index) != 1 ) 
            {
                exit(EXIT_FAILURE); 
            }   
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
        const int Id1 = Geom[0].Id(); // First node Id
        const int Id2 = Geom[1].Id(); // Second node Id
        const int Id3 = Geom[2].Id(); // Third node Id
        
        if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) // Triangle
        {
            if ( MMG3D_Set_triangle(mmgMesh, Id1, Id2, Id3, color, index) != 1 )  
            {
                exit(EXIT_FAILURE); 
            }
            
            // Set fixed boundary
            bool blocked1 = false;
            if (Geom[0].IsDefined(BLOCKED) == true)
            {
                blocked1 = Geom[0].Is(BLOCKED);
            }
            bool blocked2 = false;
            if (Geom[1].IsDefined(BLOCKED) == true)
            {
                blocked2 = Geom[1].Is(BLOCKED);
            }
            bool blocked3 = false;
            if (Geom[2].IsDefined(BLOCKED) == true)
            {
                blocked3 = Geom[2].Is(BLOCKED);
            }
            
            if ((blocked1 && blocked2 && blocked3) == true)
            {
                if ( MMG3D_Set_requiredTriangle(mmgMesh, index) != 1 ) 
                {
                    exit(EXIT_FAILURE); 
                }   
            }
        }
        else if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4) // Quadrilaterals
        {
            const int Id4 = Geom[3].Id(); // Fourth node Id
            
            if ( MMG3D_Set_quadrilateral(mmgMesh, Id1, Id2, Id3, Id4, color, index) != 1 )  
            {
                exit(EXIT_FAILURE); 
            }
        }
        else
        {
            const unsigned int SizeGeometry = Geom.size();
            KRATOS_ERROR << "WARNING: I DO NOT KNOW WHAT IS THIS. Size: " << SizeGeometry << std::endl;
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
        const int Id1 = Geom[0].Id(); // First node Id
        const int Id2 = Geom[1].Id(); // Second node Id
        const int Id3 = Geom[2].Id(); // Third node Id
        
        if ( MMG2D_Set_triangle(mmgMesh, Id1, Id2, Id3, color, index) != 1 ) 
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
        const int Id1 = Geom[0].Id(); // First node Id
        const int Id2 = Geom[1].Id(); // Second node Id
        const int Id3 = Geom[2].Id(); // Third node Id
        const int Id4 = Geom[3].Id(); // Fourth node Id
        
        if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) // Tetrahedron
        {
            if ( MMG3D_Set_tetrahedron(mmgMesh, Id1, Id2, Id3, Id4, color, index) != 1 )  
            {
                exit(EXIT_FAILURE); 
            }
        }
        else if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Prism3D6) // Prisms
        {
            const int Id5 = Geom[4].Id(); // 5th node Id
            const int Id6 = Geom[5].Id(); // 6th node Id
            
            if ( MMG3D_Set_prism(mmgMesh, Id1, Id2, Id3, Id4, Id5, Id6, color, index) != 1 )  
            {
                exit(EXIT_FAILURE); 
            }
        }
        else if (Geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8) // Hexaedron
        {
//                 const int Id5 = Geom[4].Id(); // 5th node Id
//                 const int Id6 = Geom[5].Id(); // 6th node Id
//                 const int Id6 = Geom[7].Id(); // 7th node Id
//                 const int Id6 = Geom[8].Id(); // 8th node Id
            
            const unsigned int SizeGeometry = Geom.size();
            KRATOS_ERROR << "WARNING: HEXAEDRON NON IMPLEMENTED IN THE LIBRARY " << SizeGeometry << std::endl;
        }
        else
        {
            const unsigned int SizeGeometry = Geom.size();
            KRATOS_ERROR << "WARNING: I DO NOT KNOW WHAT IS THIS. Size: " << SizeGeometry << std::endl;
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<2>::SetMetricScalar(
        const double& Metric,
        const int NodeId 
        )
    {
        if ( MMG2D_Set_scalarSol(mmgSol, Metric, NodeId) != 1 )
        {
            exit(EXIT_FAILURE);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<3>::SetMetricScalar(
        const double& Metric,
        const int NodeId 
        )
    {
        if ( MMG3D_Set_scalarSol(mmgSol, Metric, NodeId) != 1 )
        {
            exit(EXIT_FAILURE);
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<2>::SetMetricVector(
        const array_1d<double, 3>& Metric,
        const int NodeId 
        )
    {
        KRATOS_ERROR << "WARNING:: Vector metric not avalaible in 2D" << std::endl;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<3>::SetMetricVector(
        const array_1d<double, 3>& Metric,
        const int NodeId 
        )
    {
        if ( MMG3D_Set_vectorSol(mmgSol, Metric[0], Metric[1], Metric[2], NodeId) != 1 )
        {
            exit(EXIT_FAILURE);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<2>::SetMetricTensor(
        const Vector& Metric,
        const int NodeId 
        )
    {
        if ( MMG2D_Set_tensorSol(mmgSol, Metric[0],  Metric[1], Metric[2], NodeId) != 1 )
        {
            exit(EXIT_FAILURE);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void MmgUtility<3>::SetMetricTensor(
        const Vector& Metric,
        const int NodeId 
        )
    {
        if ( MMG3D_Set_tensorSol(mmgSol, Metric[0], Metric[1], Metric[2], Metric[3], Metric[4], Metric[5], NodeId) != 1 )
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
