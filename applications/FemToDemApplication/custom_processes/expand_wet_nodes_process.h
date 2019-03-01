//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:             BSD License
//                               Kratos default license:
//kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#if !defined(KRATOS_EXPAND_WET_NODES_PROCESS)
#define KRATOS_EXPAND_WET_NODES_PROCESS

#include "includes/model_part.h"
#include "processes/process.h"
#include "fem_to_dem_application_variables.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include <list>

namespace Kratos {

typedef std::size_t SizeType;

class ExpandWetNodesProcess : public Process 
{

public:
    /// Pointer definition of ExpandWetNodesProcess
    KRATOS_CLASS_POINTER_DEFINITION(ExpandWetNodesProcess);

    // Constructor
    ExpandWetNodesProcess(ModelPart &r_model_part);

    // Destructor
    ~ExpandWetNodesProcess() override = default;

    void operator()() { Execute(); }

    void Execute() override;

    bool ElementHasWetNodes(
        ModelPart::ElementsContainerType::ptr_iterator itElem,
        int& rPressureId);

    void ExpandWetNodes(
        ModelPart::ElementsContainerType::ptr_iterator itElem,
        const int PressureId);

protected:
    // Member Variables
    ModelPart& mrModelPart;
};
}
#endif /* KRATOS_EXPAND_WET_NODES_PROCESS defined */