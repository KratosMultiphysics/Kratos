//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "containers/model.h"

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
 * @brief This struct is used in order to identify when using the hitorical and non historical variables
 */
struct CalculateDistanceToPathSettings
{
    // Defining clearer options
    constexpr static bool SaveAsHistoricalVariable = true;
    constexpr static bool SaveAsNonHistoricalVariable = false;
};

/**
 * @class CalculateDistanceToPathProcess
 * @ingroup KratosCore 
 * @brief Computes DISTANCE from a path model part to a given model part
 * @author Vicente Mataix Ferrandiz
 * @tparam THistorical If the distance is computed as historical or non historical variable
 */
template<bool THistorical = true>
class KRATOS_API(KRATOS_CORE) CalculateDistanceToPathProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Index type definition
    typedef std::size_t IndexType;

    /// Size type definition
    typedef std::size_t SizeType;

    /// The definition of the node
    typedef Node NodeType;

    /// Pointer definition of CalculateDistanceToPathProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateDistanceToPathProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param rModel The model containing the model parts
     */
    CalculateDistanceToPathProcess(
        Model& rModel,
        Parameters ThisParameters
        );

    /// Destructor.
    ~CalculateDistanceToPathProcess() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method creates an pointer of the process
     * @details We consider as input a Model and a set of Parameters for the sake of generality
     * @warning Must be overrided in each process implementation
     * @param rModel The model to be consider
     * @param ThisParameters The configuration parameters
     */
    Process::Pointer Create(
        Model& rModel,
        Parameters ThisParameters
        ) override
    {
        return Kratos::make_shared<CalculateDistanceToPathProcess>(rModel, ThisParameters);
    }

    /**
     * @brief Execute method is used to execute the Process algorithms.
     */
    void Execute() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CalculateDistanceToPathProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CalculateDistanceToPathProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
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
    
    Model& mrModel;                             /// The model containing the ModelParts of the path and the distance to be computed
    Parameters mThisParameters;                 /// The configuration parameters
    const Variable<double>* mpDistanceVariable; /// The distance variable to be computed

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method computes the distance to the path
     * @param rModelPart The model part to compute the distance
     * @param rVectorSegments The vector of segments
     */
    void CalculateDistance(
        ModelPart& rModelPart,
        std::vector<Geometry<NodeType>::Pointer>& rVectorSegments
        );

    /**
     * @brief This method computes the distance to the path (by brute force)
     * @param rModelPart The model part to compute the distance
     * @param rVectorSegments The vector of segments
     */
    void CalculateDistanceByBruteForce(
        ModelPart& rModelPart,
        std::vector<Geometry<NodeType>::Pointer>& rVectorSegments
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

    /// Assignment operator.
    CalculateDistanceToPathProcess& operator=(CalculateDistanceToPathProcess const& rOther);

    /// Copy constructor.
    //CalculateDistanceToPathProcess(CalculateDistanceToPathProcess const& rOther);

    ///@}

}; // Class CalculateDistanceToPathProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<bool THistorical>
inline std::istream& operator >> (std::istream& rIStream,
                                  CalculateDistanceToPathProcess<THistorical>& rThis);

/// output stream function
template<bool THistorical>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CalculateDistanceToPathProcess<THistorical>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@}
} /// namespace Kratos