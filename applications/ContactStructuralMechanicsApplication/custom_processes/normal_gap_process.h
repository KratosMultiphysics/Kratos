// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:		 BSD License
//					 license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_NORMAL_GAP_PROCESS_H_INCLUDED )
#define  KRATOS_NORMAL_GAP_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "processes/simple_mortar_mapper_process.h"

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
 * @class NormalGapProcess
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This process computes the normal gap
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster = TNumNodes>
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) NormalGapProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// The type of mapper considered
    typedef SimpleMortarMapperProcess<TDim, TNumNodes, Variable<array_1d<double, 3>>, TNumNodesMaster> MapperType;

    /// General type definitions
    typedef ModelPart::NodesContainerType                    NodesArrayType;

    /// The definition of zero tolerance
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    /// Pointer definition of NormalGapProcess
    KRATOS_CLASS_POINTER_DEFINITION( NormalGapProcess );

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief The constructor of the normal gap process uses the following inputs:
     * @param rMasterModelPart The master model part to be considered
     * @param rSlaveModelPart The slave model part to be considered
     */
    NormalGapProcess(
        ModelPart& rMasterModelPart,
        ModelPart& rSlaveModelPart,
        const bool SearchOrientation = true
        ) : mrMasterModelPart(rMasterModelPart),
            mrSlaveModelPart(rSlaveModelPart),
            mSearchOrientation(SearchOrientation)
    {
    }

    virtual ~NormalGapProcess()= default;;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Execute method is used to execute the Process algorithms.
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
        return "NormalGapProcess";
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

    ModelPart& mrMasterModelPart;  /// The master model part to be considered
    ModelPart& mrSlaveModelPart;   /// The slave model part to be considered
    const bool mSearchOrientation; /// The orientation of the search (inverted or not)

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method switchs the flag of an array of nodes
     * @param rNodes The set of nodes where the flags are reset
     */
    static inline void SwitchFlagNodes(NodesArrayType& rNodes)
    {
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(rNodes.size()); ++i) {
            auto it_node = rNodes.begin() + i;
            it_node->Flip(SLAVE);
            it_node->Flip(MASTER);
        }
    }

    /**
     * @brief This method computes the normal gap
     * @param rNodes The set of nodes where the gap is computed
     */
    void ComputeNormalGap(NodesArrayType& rNodes);

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

    ///@}

}; // Class NormalGapProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/****************************** INPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
inline std::istream& operator >> (std::istream& rIStream,
                                  NormalGapProcess<TDim, TNumNodes, TNumNodesMaster>& rThis);

/***************************** OUTPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const NormalGapProcess<TDim, TNumNodes, TNumNodesMaster>& rThis)
{
    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_NORMAL_GAP_PROCESS_H_INCLUDED  defined
