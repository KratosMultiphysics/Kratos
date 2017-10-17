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

#if !defined(KRATOS_NODAL_VALUES_INTERPOLATION_PROCESS )
#define  KRATOS_NODAL_VALUES_INTERPOLATION_PROCESS

// System includes

// External includes

// Project includes
#include "meshing_application.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
// Include the point locator
#include "utilities/binbased_fast_point_locator.h"

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

/** \brief NodalValuesInterpolationProcess
 * This utilitiy has as objective to interpolate the values inside elements (and conditions?) in a model part, using as input the original model part and the new one
 * The process employs the projection.h from MeshingApplication, which works internally using a kd-tree 
 */

template<unsigned int TDim>
class NodalValuesInterpolationProcess 
    : public Process
{
public:
    ///@name Type Definitions
    ///@{
    
    // General type definitions
    typedef ModelPart::NodesContainerType                    NodesArrayType;
    typedef ModelPart::ElementsContainerType              ElementsArrayType;
    typedef ModelPart::ConditionsContainerType          ConditionsArrayType;
    typedef Node<3>                                                NodeType;
    typedef Geometry<NodeType>                                 GeometryType;

    /// Pointer definition of NodalValuesInterpolationProcess
    KRATOS_CLASS_POINTER_DEFINITION( NodalValuesInterpolationProcess );
      
    ///@}
    ///@name Life Cycle
    ///@{

    // Class Constructor
    
    /**
     * The constructor of the search utility uses the following inputs:
     * @param rOriginMainModelPart: The model part from where interpolate values
     * @param rDestinationMainModelPart: The model part where we want to interpolate the values
     * @param ThisParameters: The parameters containing all the information needed
     */
    
    NodalValuesInterpolationProcess(
        ModelPart& rOriginMainModelPart,
        ModelPart& rDestinationMainModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        )
    :mrOriginMainModelPart(rOriginMainModelPart),
     mrDestinationMainModelPart(rDestinationMainModelPart)
     {
         Parameters DefaultParameters = Parameters(R"(
         {
            "echo_level"            : 1, 
            "framework"             : "Eulerian", 
            "max_number_of_searchs" : 1000, 
            "step_data_size"        : 0, 
            "buffer_size"           : 0
         })");
         ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
         
         mEchoLevel = ThisParameters["echo_level"].GetInt();
         mFramework = ConvertFramework(ThisParameters["framework"].GetString());
         mMaxNumberOfResults = ThisParameters["max_number_of_searchs"].GetInt();
         mStepDataSize = ThisParameters["step_data_size"].GetInt();
         mBufferSize   = ThisParameters["buffer_size"].GetInt();
        
         if (mEchoLevel > 0)
         {
             std::cout << "Step data size: " << mStepDataSize << " Buffer size: " << mBufferSize << std::endl;
         }
     }
    
    ~NodalValuesInterpolationProcess() override= default;;

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
    
    /**
     * We execute the search relative to the old and new model part
     */
    
    void Execute() override
    {
        // We create the locator
        BinBasedFastPointLocator<TDim> PointLocator = BinBasedFastPointLocator<TDim>(mrOriginMainModelPart);
        PointLocator.UpdateSearchDatabase();
        
        // Iterate in the nodes
        NodesArrayType& pNode = mrDestinationMainModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        /* Nodes */
//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            Vector ShapeFunctions;
            Element::Pointer pElement;
            
            const bool IsFound = PointLocator.FindPointOnMeshSimplified(itNode->Coordinates(), ShapeFunctions, pElement, mMaxNumberOfResults, 5.0e-2);
            
            if (IsFound == false)
            {
                if (mEchoLevel > 0 || mFramework == Lagrangian) // NOTE: In the case we are in a Lagrangian framework this is serious and should print a message
                {
                   std::cout << "WARNING: Node "<< itNode->Id() << " not found (interpolation not posible)" << std::endl;
                   std::cout << "\t X:"<< itNode->X() << "\t Y:"<< itNode->Y() << "\t Z:"<< itNode->Z() << std::endl;
                   
                   if (mFramework == Lagrangian)
                   {
                       std::cout << "WARNING: YOU ARE IN A LAGRANGIAN FRAMEWORK THIS IS DANGEROUS" << std::endl;
                   }
                }
            }
            else
            {
                for(unsigned int iStep = 0; iStep < mBufferSize; iStep++)
                {
                    CalculateStepData(*(itNode.base()), pElement, ShapeFunctions, iStep);
                }
                
                // After we interpolate the DISPLACEMENT we interpolate the initial coordinates to ensure a functioning of the simulation
                if (mFramework == Lagrangian)
                {
                    if ( itNode->SolutionStepsDataHas( DISPLACEMENT ) == false ) // Fisrt we check if we have the displacement variable
                    {
                        KRATOS_ERROR << "Missing DISPLACEMENT on node " << itNode->Id() << std::endl;
                    }
                    
                    CalculateInitialCoordinates(*(itNode.base()));
                }
            }
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

    /************************************ GET INFO *************************************/
    /***********************************************************************************/
    
    std::string Info() const override
    {
        return "NodalValuesInterpolationProcess";
    }

    /************************************ PRINT INFO ***********************************/
    /***********************************************************************************/
    
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
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
    
    ModelPart& mrOriginMainModelPart;                    // The origin model part
    ModelPart& mrDestinationMainModelPart;               // The destination model part
    unsigned int mMaxNumberOfResults;                    // The maximum number of results to consider in the search
    unsigned int mStepDataSize;                          // The size of the database
    unsigned int mBufferSize;                            // The size of the buffer
    FrameworkEulerLagrange mFramework;                   // The framework
    unsigned int mEchoLevel;                             // The level of verbosity
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    /**
     * It calculates the initial coordinates interpolated to the node
     * @return itNode: The node pointer
     */
    
    void CalculateInitialCoordinates(NodeType::Pointer& pNode);
    
    /**
     * It calculates the Step data interpolated to the node
     * @return itNode: The node pointer
     * @param pElement: The element pointer
     */
    
    void CalculateStepData(
        NodeType::Pointer& pNode,
        const Element::Pointer& pElement,
        const Vector ShapeFunctions,
        const unsigned int Step
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
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class NodalValuesInterpolationProcess

///@name Explicit Specializations
///@{
    
    template<>  
    void NodalValuesInterpolationProcess<2>::CalculateInitialCoordinates(NodeType::Pointer& pNode)
    {
        // We interpolate the initial coordinates (X = X0 + DISPLACEMENT), then X0 = X - DISPLACEMENT        
        pNode->X0() = pNode->X() - pNode->FastGetSolutionStepValue(DISPLACEMENT_X);
        pNode->Y0() = pNode->Y() - pNode->FastGetSolutionStepValue(DISPLACEMENT_Y);
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void NodalValuesInterpolationProcess<3>::CalculateInitialCoordinates(NodeType::Pointer& pNode)
    {        
        // We interpolate the initial coordinates (X = X0 + DISPLACEMENT), then X0 = X - DISPLACEMENT        
        pNode->X0() = pNode->X() - pNode->FastGetSolutionStepValue(DISPLACEMENT_X);
        pNode->Y0() = pNode->Y() - pNode->FastGetSolutionStepValue(DISPLACEMENT_Y);
        pNode->Z0() = pNode->Z() - pNode->FastGetSolutionStepValue(DISPLACEMENT_Z);
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void NodalValuesInterpolationProcess<2>::CalculateStepData(
        NodeType::Pointer& pNode,
        const Element::Pointer& pElement,
        const Vector ShapeFunctions,
        const unsigned int Step
        )
    {
        double* StepData = pNode->SolutionStepData().Data(Step);
        
        double* NodeData0 = pElement->GetGeometry()[0].SolutionStepData().Data(Step);
        double* NodeData1 = pElement->GetGeometry()[1].SolutionStepData().Data(Step);
        double* NodeData2 = pElement->GetGeometry()[2].SolutionStepData().Data(Step);
        
        for (unsigned int j = 0; j < mStepDataSize; j++)
        {
            StepData[j] = ShapeFunctions[0] * NodeData0[j]
                        + ShapeFunctions[1] * NodeData1[j]
                        + ShapeFunctions[2] * NodeData2[j];
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    void NodalValuesInterpolationProcess<3>::CalculateStepData(
        NodeType::Pointer& pNode,
        const Element::Pointer& pElement,
        const Vector ShapeFunctions,
        const unsigned int Step
        )
    {
        double* StepData = pNode->SolutionStepData().Data(Step);
        
        // NOTE: This just works with tetrahedron (you are going to have problems with anything else)
        double* NodeData0 = pElement->GetGeometry()[0].SolutionStepData().Data(Step);
        double* NodeData1 = pElement->GetGeometry()[1].SolutionStepData().Data(Step);
        double* NodeData2 = pElement->GetGeometry()[2].SolutionStepData().Data(Step);
        double* NodeData3 = pElement->GetGeometry()[3].SolutionStepData().Data(Step);
        
        for (unsigned int j = 0; j < mStepDataSize; j++)
        {
            StepData[j] = ShapeFunctions[0] * NodeData0[j]
                        + ShapeFunctions[1] * NodeData1[j]
                        + ShapeFunctions[2] * NodeData2[j]
                        + ShapeFunctions[3] * NodeData3[j];
        }
    }

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/****************************** INPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

// template<class TPointType, class TPointerType>
// inline std::istream& operator >> (std::istream& rIStream,
//                                   NodalValuesInterpolationProcess& rThis);

/***************************** OUTPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

// template<class TPointType, class TPointerType>
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const NodalValuesInterpolationProcess& rThis)
// {
//     return rOStream;
// }

///@}

}  // namespace Kratos.

#endif // KRATOS_NODAL_VALUES_INTERPOLATION_PROCESS  defined 
