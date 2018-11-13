//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velázquez
//

#if !defined(KRATOS_EXTEND_PRESSURE_PROCESS)
#define KRATOS_EXTEND_PRESSURE_PROCESS

#include "includes/model_part.h"
#include "processes/process.h"
#include <pybind11/pybind11.h>
#include <list>
#include "fem_to_dem_application_variables.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"

namespace Kratos {

typedef std::size_t SizeType;

template <SizeType TDim = 2>
class ExtendPressureConditionProcess : public Process 
{

public:
    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ExtendPressureConditionProcess);

    // Constructor
    ExtendPressureConditionProcess(ModelPart &r_model_part);

    // Destructor
    ~ExtendPressureConditionProcess() override = default;

    void operator()() { Execute(); }

    void Execute() override;

    void CreateAndAddPressureConditions2(
        ModelPart::ElementsContainerType::ptr_iterator itElem,
        const unsigned int LocalId,
        const int PressureId,
        int& MaximumConditionId,
        std::vector<IndexType>& ToEraseConditionsId);

    void CreateAndAddPressureConditions3(
        ModelPart::ElementsContainerType::ptr_iterator itElem,
        const int PressureId,
        int& MaximumConditionId,
        std::vector<IndexType>& ToEraseConditionsId);

    void GetMaximumConditionIdOnSubmodelPart(
        int& MaximumConditionId);

    void CalculateNumberOfElementsOnNodes();
    bool CheckIfHasConditionId(const IndexType Id);
protected:
    // Member Variables
    ModelPart &mr_model_part;
    std::vector<IndexType> mNodeIdContainer;
    std::vector<int> mNodePressureIdContainer;

};  // Class

}  // namespace Kratos
#endif /* KRATOS_EXTEND_PRESSURE_PROCESS defined */