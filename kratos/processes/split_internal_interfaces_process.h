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
//  Collaborators:   Miguel Angel Celigueta
//                   Ruben Zorrilla
//

#if !defined(KRATOS_SPLIT_INTERNAL_INTERFACES_PROCESS_H_INCLUDED )
#define  KRATOS_SPLIT_INTERNAL_INTERFACES_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "processes/find_elements_neighbours_process.h"

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
 * @class SplitInternalInterfacesProcess
 * @ingroup KratosCore
 * @brief Computes NODAL_AREA
 * @details splits a domain across changes of property and generates a condition at the splitting positions
 * @author Riccardo Rossi
 * @author Miguel Angel Celigueta
 */
class KRATOS_API(KRATOS_CORE) SplitInternalInterfacesProcess : public Process {

public:
    ///@name Type Definitions
    ///@{

    /// Index type definition
    typedef std::size_t IndexType;

    /// Size type definition
    typedef std::size_t SizeType;

    /// The definition of the node
    typedef Node<3> NodeType;

    /// Pointer definition of SplitInternalInterfacesProcess
    KRATOS_CLASS_POINTER_DEFINITION(SplitInternalInterfacesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param DomainSize The size of the space, if the value is not provided will compute from the model part
     */
    SplitInternalInterfacesProcess(
        Model& rModel,
        Parameters rParameters)
        : Process()
        , mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString()))
    {
        KRATOS_TRY

        // Validate input parameters
        const Parameters default_parameters = GetDefaultParameters();
        rParameters.ValidateAndAssignDefaults(default_parameters);

        // Get or create the interfaces sub model part
        std::string interfaces_sub_model_part_name = rParameters["interfaces_sub_model_part_name"].GetString();
        if (mrModelPart.HasSubModelPart(interfaces_sub_model_part_name)) {
            mpInterfacesSubModelPart = &(mrModelPart.GetSubModelPart(interfaces_sub_model_part_name));
        } else {
            mpInterfacesSubModelPart = &(mrModelPart.CreateSubModelPart(interfaces_sub_model_part_name));
        }

        // Get the registering name of the conditions to be created
        mConditionName = rParameters["condition_name"].GetString();

        KRATOS_CATCH("");
    }

    /// Destructor.
    virtual ~SplitInternalInterfacesProcess() override { }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override;

    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters( R"(
        {
            "model_part_name" :"MODEL_PART_NAME",
            "interfaces_sub_model_part_name" : "internal_interfaces",
            "condition_name"   : "CONDITION_NAME"
        }  )" );
        return default_parameters;
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

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "SplitInternalInterfacesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SplitInternalInterfacesProcess";
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

    void SplitBoundary(
        const std::size_t PropertyIdBeingProcessed,
        const std::size_t InterfaceConditionsPropertyId,
        const std::unordered_map<IndexType,IndexType>& rConditionsParentMap,
        ModelPart& rModelPart);

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

    ModelPart& mrModelPart;  /// The origin model part from which the interfaces are to be created
    ModelPart* mpInterfacesSubModelPart;  /// The sub model part in which the interfaces are to be created
    std::string mConditionName;    /// The dimension of the space

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

    /// Assignment operator.
    //SplitInternalInterfacesProcess& operator=(SplitInternalInterfacesProcess const& rOther);

    /// Copy constructor.
    //SplitInternalInterfacesProcess(SplitInternalInterfacesProcess const& rOther);


    ///@}

}; // Class SplitInternalInterfacesProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, SplitInternalInterfacesProcess& rThis) {
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const SplitInternalInterfacesProcess& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_SPLIT_INTERNAL_INTERFACES_PROCESS_H_INCLUDED  defined
