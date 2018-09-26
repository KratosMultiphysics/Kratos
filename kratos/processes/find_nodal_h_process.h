//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//  Collaborator:    Vicente Mataix Ferrandiz
//                    
//

#if !defined(KRATOS_FIND_NODAL_H_PROCESS_INCLUDED )
#define  KRATOS_FIND_NODAL_H_PROCESS_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
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
 * @class FindNodalHProcess
 * @ingroup KratosCore 
 * @brief Computes NODAL_H
 * @details Calculate the NODAL_H for all the nodes by means of the element sides minimum length
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */
template<bool THistorical = true>
class KRATOS_API(KRATOS_CORE) FindNodalHProcess 
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
    typedef Node<3> NodeType;
    
    /// The definition of the node iterator
    typedef ModelPart::NodeIterator NodeIterator;
    
    /// Pointer definition of FindNodalHProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindNodalHProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FindNodalHProcess(ModelPart& rModelPart) 
        : mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    ~FindNodalHProcess() override
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
        return "FindNodalHProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FindNodalHProcess";
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
    
    ModelPart& mrModelPart;  /// The model part were to compute the NODAL_H

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method gets the current value of the NODAL_H
     * @param rNode The node iterator to be get
     * @return The current value of NODAL_H
     */
    double GetValue(NodeType& rNode);
    
    /**
     * @brief This method sets the current value of the NODAL_H to the given one
     * @param rNode The node iterator to be get
     * @param Value The current value of NODAL_H
     */
    void SetValue(
        NodeType& rNode,
        const double Value
        );
    
    /**
     * @brief This method sets the current value of the NODAL_H to the maximum
     * @param itNode The node iterator to be set
     */
    void SetInitialValue(NodeIterator itNode);

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
    FindNodalHProcess& operator=(FindNodalHProcess const& rOther);

    /// Copy constructor.
    //FindNodalHProcess(FindNodalHProcess const& rOther);


    ///@}

}; // Class FindNodalHProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<bool THistorical>
inline std::istream& operator >> (std::istream& rIStream,
                                  FindNodalHProcess<THistorical>& rThis);

/// output stream function
template<bool THistorical>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FindNodalHProcess<THistorical>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_FIND_NODAL_H_PROCESS_INCLUDED  defined 


