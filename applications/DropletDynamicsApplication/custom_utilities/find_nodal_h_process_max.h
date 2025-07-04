//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"

namespace Kratos
{

/**
 * @brief This struct is used in order to identify when using the hitorical and non historical variables
 */
struct FindNodalHSettings
{
    // Defining clearer options
    constexpr static bool SaveAsHistoricalVariable = true;
    constexpr static bool SaveAsNonHistoricalVariable = false;
};

/**
 * @class FindNodalHProcessMax
 * @ingroup Kratos.DropletDynamicsApplication
 * @brief Computes NODAL_H_MAX
 * @details Calculate the NODAL_H_MAX for all the nodes by means of the element sides maximum length
 * @author 
 */
template<bool THistorical = true>
class KRATOS_API(DROPLET_DYNAMICS_APPLICATION) FindNodalHProcessMax
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

    /// The definition of the node iterator
    typedef ModelPart::NodeIterator NodeIterator;

    /// Pointer definition of FindNodalHProcessMax
    KRATOS_CLASS_POINTER_DEFINITION(FindNodalHProcessMax);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit FindNodalHProcessMax(ModelPart& rModelPart)
        : mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    ~FindNodalHProcessMax() override = default;

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

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "FindNodalHProcessMax";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FindNodalHProcessMax";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;  /// The model part were to compute the NODAL_H_MAX

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method gets the current value of the NODAL_H_MAX
     * @param rNode The node iterator to be get
     * @return The current value of NODAL_H_MAX
     */
    double& GetHValue(NodeType& rNode);

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
    FindNodalHProcessMax& operator=(FindNodalHProcessMax const& rOther);

    /// Copy constructor.
    //FindNodalHProcessMax(FindNodalHProcessMax const& rOther);

    ///@}
}; // Class FindNodalHProcessMax

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<bool THistorical>
inline std::istream& operator >> (std::istream& rIStream,
                                  FindNodalHProcessMax<THistorical>& rThis);

/// output stream function
template<bool THistorical>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FindNodalHProcessMax<THistorical>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.
