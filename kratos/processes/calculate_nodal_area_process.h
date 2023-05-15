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
//  Collaborators:   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_CALCULATE_NODAL_AREA_PROCESS_H_INCLUDED )
#define  KRATOS_CALCULATE_NODAL_AREA_PROCESS_H_INCLUDED

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
 * @brief This struct is used in order to identify when using the hitorical and non historical variables
 */
struct CalculateNodalAreaSettings
{
    // Defining clearer options
    constexpr static bool SaveAsHistoricalVariable = true;
    constexpr static bool SaveAsNonHistoricalVariable = false;
};
    
/** 
 * @class CalculateNodalAreaProcess
 * @ingroup KratosCore 
 * @brief Computes NODAL_AREA
 * @details Calculate the NODAL_AREA for computing the weighted area in each node
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */
template<bool THistorical = true>
class KRATOS_API(KRATOS_CORE) CalculateNodalAreaProcess
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
    
    /// Pointer definition of CalculateNodalAreaProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateNodalAreaProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param DomainSize The size of the space, if the value is not provided will compute from the model part
     */
    CalculateNodalAreaProcess(
        ModelPart& rModelPart, 
        const SizeType DomainSize = 0
        ): mrModelPart(rModelPart),
           mDomainSize(DomainSize)
    {
        // In case is not provided we will take from the model part
        if (mDomainSize == 0) {
            KRATOS_ERROR_IF_NOT(rModelPart.GetProcessInfo().Has(DOMAIN_SIZE)) << "\"DOMAIN_SIZE\" has to be specified in the ProcessInfo" << std::endl;
            mDomainSize = rModelPart.GetProcessInfo()[DOMAIN_SIZE];
        }
    }

    /// Destructor.
    ~CalculateNodalAreaProcess() override
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
        return "CalculateNodalAreaProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CalculateNodalAreaProcess";
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
    
    ModelPart& mrModelPart;  /// The model part where the nodal area is computed
    SizeType mDomainSize;    /// The dimension of the space

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method gets the current value of the NODAL_AREA
     * @param rNode The node iterator to be get
     * @return The current value of NODAL_AREA
     */
    double& GetAreaValue(NodeType& rNode);

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
    CalculateNodalAreaProcess& operator=(CalculateNodalAreaProcess const& rOther);

    /// Copy constructor.
    //CalculateNodalAreaProcess(CalculateNodalAreaProcess const& rOther);


    ///@}

}; // Class CalculateNodalAreaProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<bool THistorical = true>
inline std::istream& operator >> (std::istream& rIStream,
                                  CalculateNodalAreaProcess<THistorical>& rThis);

/// output stream function
template<bool THistorical = true>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CalculateNodalAreaProcess<THistorical>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_CALCULATE_NODAL_AREA_PROCESS_H_INCLUDED  defined 


