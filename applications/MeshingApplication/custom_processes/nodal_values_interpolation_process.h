// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_NODAL_VALUES_INTERPOLATION_PROCESS )
#define  KRATOS_NODAL_VALUES_INTERPOLATION_PROCESS

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "meshing_application.h"
#include "processes/process.h"

/* Several includes */
#include "includes/model_part.h"
#include "includes/key_hash.h"
#include "includes/kratos_parameters.h"

/* Utilities */
#include "utilities/binbased_fast_point_locator.h"

/* Tree structures */
// #include "spatial_containers/bounding_volume_tree.h" // k-DOP
#include "spatial_containers/spatial_containers.h" // kd-tree

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// Defining the integers
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

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
 * @brief This utilitiy has as objective to interpolate the values inside elements (and conditions?) in a model part, using as input the original model part and the new one
 * @details The process employs the projection.h from MeshingApplication, which works internally using a kd-tree. Additionally if it can't found the node inside the reference mesh it will try to extrapolate from the skin (if the option is activated)
 * @author Vicente Mataix Ferrandiz
 */
template<SizeType TDim>
class KRATOS_API(KRATOS_MESHING_APPLICATION) NodalValuesInterpolationProcess
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

    ModelPart& mrOriginMainModelPart;      /// The origin model part
    ModelPart& mrDestinationMainModelPart; /// The destination model part
    Parameters mThisParameters;            /// Here the configuration parameters are stored
    std::unordered_set<Variable<double>, VariableHasher<Variable<double>>, VariableComparator<Variable<double>>> mListDoublesVariables;             /// List of double non-historical variables
    std::unordered_set<Variable<array_1d<double, 3>>, VariableHasher<Variable<array_1d<double, 3>>>, VariableComparator<Variable<array_1d<double, 3>>>> mListArraysVariables; /// List of array non-historical variables
    std::unordered_set<Variable<Vector>, VariableHasher<Variable<Vector>>, VariableComparator<Variable<Vector>>> mListVectorVariables;              /// List of vector non-historical variables
    std::unordered_set<Variable<Matrix>, VariableHasher<Variable<Matrix>>, VariableComparator<Variable<Matrix>>> mListMatrixVariables;              /// List of matrix non-historical variables

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
     * @brief It calculates the data (DataContainer) interpolated to the node
     * @param pNode The node pointer
     * @param pElement The element pointer
     * @param rShapeFunctions The shape functions
     */

    void CalculateData(
        NodeType::Pointer pNode,
        const Element::Pointer& pElement,
        const Vector& rShapeFunctions
        );

    /**
     * @brief This methoid creates the list of non-historical variables fro nodal interpolation
     */
    void GetListNonHistoricalVariables();

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
        const IndexType Step
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
