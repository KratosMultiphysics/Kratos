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

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "processes/skin_detection_process.h"
#include "processes/process.h"
#include "fem_to_dem_application_variables.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"

namespace Kratos {

typedef std::size_t SizeType;

/** 
 * @class ExpandWetNodesProcess
 * @ingroup FemToDemApplication 
 * @brief Assigns the pressure Id to the nodes of the elements that are going to be removed
 * @details When one element is going to be removed, this process assigns the pressure_id
 * to the rest of the nodes (expand the pressure load)
 * @author Alejandro Cornejo
 */
class ExpandWetNodesProcess : public Process 
{

public:

    typedef Node <3> NodeType;
    typedef Properties PropertiesType;
    typedef Element ElementType;
    typedef Condition ConditionType;
    typedef Mesh<NodeType, PropertiesType, ElementType, ConditionType> MeshType;
    typedef PointerVector<MeshType> MeshesContainerType;
    typedef MeshType::ElementIterator ElementIterator;

    /// Pointer definition of ExpandWetNodesProcess
    KRATOS_CLASS_POINTER_DEFINITION(ExpandWetNodesProcess);

    // Constructor
    ExpandWetNodesProcess(ModelPart &rModelPart);

    // Destructor
    ~ExpandWetNodesProcess() override = default;

    void operator()() { Execute(); }

    void Execute() override;

    bool ElementHasWetNodes(
        ElementIterator itElem,
        int& rPressureId,
        int& rNumberOfWetNodes);

    bool ElementHasWetNodes2(
        ElementIterator itElem,
        int& rPressureId,
        int& rNumberOfWetNodes);

    void ExpandWetNodes(
        ElementIterator itElem,
        const int PressureId);

    void ExpandWetNodesIfTheyAreSkin();

    void ExpandWetNodesWithLatestPressureId();

protected:
    // Member Variables
    ModelPart& mrModelPart;
    std::string mPressureName;
};
}
#endif /* KRATOS_EXPAND_WET_NODES_PROCESS defined */