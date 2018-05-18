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
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
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
    
///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/** 
 * @class NodalValuesInterpolationProcess
 * @ingroup MeshingApplication
 * @brief NodalValuesInterpolationProcess
 * @details This utilitiy has as objective to interpolate the values inside elements (and conditions?) in a model part, using as input the original model part and the new one
 * The process employs the projection.h from MeshingApplication, which works internally using a kd-tree 
 * @author Vicente Mataix Ferrandiz
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
    ///@name  Enum's
    ///@{
    
    /**
     * @brief This enums allows to differentiate the working framework
     */
    enum class FrameworkEulerLagrange {EULERIAN = 0, LAGRANGIAN = 1, ALE = 2};
    
    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief The constructor of the search utility uses the following inputs:
     * @param rOriginMainModelPart The model part from where interpolate values
     * @param rDestinationMainModelPart The model part where we want to interpolate the values
     * @param ThisParameters The parameters containing all the information needed
     */
    
    NodalValuesInterpolationProcess(
        ModelPart& rOriginMainModelPart,
        ModelPart& rDestinationMainModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );
    
    /// Destructor
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
     * @brief We execute the search relative to the old and new model part
     */
    
    void Execute() override;
    
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
    
    ModelPart& mrOriginMainModelPart;                    /// The origin model part
    ModelPart& mrDestinationMainModelPart;               /// The destination model part
    unsigned int mMaxNumberOfResults;                    /// The maximum number of results to consider in the search
    unsigned int mStepDataSize;                          /// The size of the database
    unsigned int mBufferSize;                            /// The size of the buffer
    FrameworkEulerLagrange mFramework;                   /// The framework
    unsigned int mEchoLevel;                             /// The level of verbosity
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    /**
     * @brief This converts the framework string to an enum
     * @param Str The string
     * @return FrameworkEulerLagrange: The equivalent enum
     */
        
    static inline FrameworkEulerLagrange ConvertFramework(const std::string& Str)
    {
        if(Str == "Lagrangian" || Str == "LAGRANGIAN") 
            return FrameworkEulerLagrange::LAGRANGIAN;
        else if(Str == "Eulerian" || Str == "EULERIAN") 
            return FrameworkEulerLagrange::EULERIAN;
        else if(Str == "ALE") 
            return FrameworkEulerLagrange::ALE;
        else
            return FrameworkEulerLagrange::EULERIAN;
    }
    
    /**
     * @brief It calculates the Step data interpolated to the node
     * @param pNode The node pointer
     * @param pElement The element pointer
     * @param rShapeFunctions The shape functions
     * @param Step The current time step
     */
    
    void CalculateStepData(
        NodeType::Pointer pNode,
        const Element::Pointer& pElement,
        const Vector& rShapeFunctions,
        const unsigned int Step
        );
    
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
