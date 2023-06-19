//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_NITSCHE_STABILIZATION_MODEL_PART_PROCESS_H_INCLUDED )
#define  KRATOS_NITSCHE_STABILIZATION_MODEL_PART_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "processes/process.h"

// Application includes
#include "iga_application_variables.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/* @class NitscheStabilizationModelPartProcess
 * @ingroup IgaApplication
 * @brief This class outputs the location of the quadrature points within the local space of the containing geometry. */
class KRATOS_API(IGA_APPLICATION) NitscheStabilizationModelPartProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NitscheStabilizationModelPartProcess
    KRATOS_CLASS_POINTER_DEFINITION(NitscheStabilizationModelPartProcess);

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    NitscheStabilizationModelPartProcess(
        ModelPart& rThisModelPart);

    /// Destructor.
    ~NitscheStabilizationModelPartProcess() = default;

    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "NitscheStabilizationModelPartProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NitscheStabilizationModelPartProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

private:
    ///@name Member Variables
    ///@{

    /// Model part and different settings
    ModelPart& mrThisModelPart; /// The model part to compute

    ///@}

}; // Class NitscheStabilizationModelPartProcess

}  // namespace Kratos.

#endif // KRATOS_NITSCHE_STABILIZATION_MODEL_PART_PROCESS_H_INCLUDED  defined
