//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#pragma once

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"

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
 * @class FindConditionsNeighboursProcess
 * @ingroup KratosCore
 * @brief This process finds the neighboring conditions for each node and condition in a given model part.
 * @details This process iterates over all conditions in the model part and assigns neighboring conditions to each node and condition based on the connectivity information. It also provides functionality to clear the neighbor data if needed.
 * @author Riccardo Rossi
 * @todo In the future for MPI compatibility and reduce code duplication FindGlobalNodalEntityNeighboursProcess can be used as base class. Unfortunately it only computes for the nodes, not conditions, so something must be done in that class before properly computing.
 */
class KRATOS_API(KRATOS_CORE) FindConditionsNeighboursProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Defining the Index type
    using IndexType = std::size_t;

    /// Pointer definition of FindConditionsNeighboursProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindConditionsNeighboursProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor with parameters.
     * @details The better the guess for the quantities above the less memory occupied and the fastest the algorithm
     * @param rModel The model containing the target model part.
     * @param ThisParameters The setting parameters for the process.
     */
    FindConditionsNeighboursProcess(
        Model& rModel,
        Parameters ThisParameters
        );

    /**
     * @brief Default constructor.
     * @details The better the guess for the quantities above the less memory occupied and the fastest the algorithm
     * @param rModelPart The model part containing the conditions.
     * @param Dim The dimension of the problem.
     * @param AverageConditions The expected number of neighboring conditions per node.
     */
    FindConditionsNeighboursProcess(
        ModelPart& rModelPart,
        const int Dim = -1,
        const unsigned int AverageConditions = 10
        );

    /// Destructor.
    ~FindConditionsNeighboursProcess() override
    {
    }

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
     * @brief Executes the process.
     */
    void Execute() override;

    /**
     * @brief This method clears the assignation of the conditions
     */
    void Clear() override;

    /**
     * @brief Clears the neighbor data for all nodes and conditions.
     */
    void ClearNeighbours();

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
        return "FindConditionsNeighboursProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FindConditionsNeighboursProcess";
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

    ModelPart& mrModelPart;               /// Reference to the model part.
    unsigned int mAverageConditions = 10; /// Expected number of neighboring conditions per node.
    int mDim = -1;                        /// Dimension of the problem. NOTE: Should be a template argument

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief Compute the dimension for the FindConditionsNeighboursProcess.
     * @details This function computes the dimension of the process if it's not already defined. It retrieves the dimension from the first condition's geometry if it's not defined. If the dimension is not 2 or 3, it throws an error.
     */
    void ComputeDimension();

    /**
     * @brief Adds a unique weak pointer to a GlobalPointersVector if it does not already exist.
     * @details This function iterates through the GlobalPointersVector to check if the candidate weak pointer's ID already exists. If it does not, the candidate is added to the vector.
     * @tparam TDataType The data type of the elements stored it_node the GlobalPointersVector.
     * @param v Reference to the GlobalPointersVector to which the unique weak pointer is to be added.
     * @param candidate The weak pointer candidate to be added if it's unique based on its ID.
     */
    template< class TDataType >
    void AddUniqueWeakPointer(
        GlobalPointersVector<TDataType>& v,
        const typename TDataType::WeakPointer candidate
        )
    {
        auto it = v.begin();
        auto it_end = v.end();
        while ( it != it_end && (it)->Id() != (candidate.lock())->Id()) {
            it++;
        }
        if(it == it_end) {
            v.push_back(candidate);
        }

    }

    /**
    * @brief Checks for neighbouring faces around a given node that do not match a specified face ID.
    * @details This function iterates through the neighbour faces of a node (Id1) to find a face that has a node with Id2 and does not match the specified face ID. This is used to find adjacent or connected elements it_node a mesh, based on their nodes.
    * @param Id1 The ID of the node around which neighbouring faces are checked.
    * @param Id2 The ID of the second node to match it_node the neighbouring faces.
    * @param rNeighbourFace The GlobalPointersVector containing the neighbour faces to check.
    * @param Face The face ID that should not match the found neighbouring face.
    * @return Condition::WeakPointer to the first neighbouring face found that matches the criteria or an empty WeakPointer if no such face is found.
    */
    Condition::WeakPointer CheckForNeighbourFaces(
        const unsigned int Id1,
        const unsigned int Id2,
        GlobalPointersVector<Condition>& rNeighbourFace,
        const unsigned int Face
        );

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

    /// Assignment operator.
    FindConditionsNeighboursProcess& operator=(FindConditionsNeighboursProcess const& rOther);

    /// Copy constructor.
    //FindConditionsNeighboursProcess(FindConditionsNeighboursProcess const& rOther);

    ///@}
}; // Class FindConditionsNeighboursProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  FindConditionsNeighboursProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FindConditionsNeighboursProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.


