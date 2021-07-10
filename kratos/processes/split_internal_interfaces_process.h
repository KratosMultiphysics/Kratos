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
    SplitInternalInterfacesProcess(Model& rModel, Parameters rParameters):Process(Flags()), mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString())) {
        KRATOS_TRY
        const Parameters default_parameters = GetDefaultParameters();
        rParameters.ValidateAndAssignDefaults(default_parameters);
        mConditionName = rParameters["condition_name"].GetString();
        KRATOS_CATCH("");
    }

    /// Destructor.
    virtual ~SplitInternalInterfacesProcess() override { }


    const Parameters GetDefaultParameters() const override {
        const Parameters default_parameters( R"(
        {
            "model_part_name" :"MODEL_PART_NAME",
            "condition_name"   : "CONDITION_NAME"
        }  )" );
        return default_parameters;
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override;

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
    std::string Info() const override {
        return "SplitInternalInterfacesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << "SplitInternalInterfacesProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
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
    void SplitBoundary(const std::size_t PropertyIdBeingProcessed, ModelPart& rModelPart);

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


