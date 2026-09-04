//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_IGA_VTK_OUTPUT_PROCESS_H_INCLUDED)
#define KRATOS_IGA_VTK_OUTPUT_PROCESS_H_INCLUDED

// System includes
#include <memory>

// Project includes
#include "containers/model.h"
#include "controllers/output_controller.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "processes/output_process.h"

namespace Kratos
{

///@class IgaVtkOutputProcess
///@brief Writes sampled IGA surfaces and curves in the VTKHDF format.
class KRATOS_API(IGA_APPLICATION) IgaVtkOutputProcess : public OutputProcess
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(IgaVtkOutputProcess);

    IgaVtkOutputProcess(Model& rModel, Parameters ThisParameters);

    ~IgaVtkOutputProcess() override = default;

    void ExecuteInitialize() override;
    void ExecuteBeforeSolutionLoop() override;
    bool IsOutputStep() override;
    void PrintOutput() override;

    const Parameters GetDefaultParameters() const override;

    std::string Info() const override
    {
        return "IgaVtkOutputProcess";
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IgaVtkOutputProcess";
    }

    void PrintData(std::ostream& rOStream) const override
    {
    }

private:
    Model& mrModel;
    Parameters mThisParameters;
    ModelPart* mpModelPart = nullptr;
    std::unique_ptr<OutputController> mpOutputController;
};

} // namespace Kratos

#endif // KRATOS_IGA_VTK_OUTPUT_PROCESS_H_INCLUDED
