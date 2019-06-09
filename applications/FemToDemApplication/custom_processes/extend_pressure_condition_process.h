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
//  Main authors:    Alejandro Cornejo Velazquez
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

/** 
 * @class ExtendPressureConditionProcess
 * @ingroup FemToDemApplication 
 * @brief Creates the new presure line loads after removing some elements
 * @details when several elements are removed this methods generates the line loads
 * in order to adapt to the new geometry
 * @author Alejandro Cornejo
 */
template <SizeType TDim = 2>
class ExtendPressureConditionProcess : public Process 
{

public:
    /// Pointer definition of ExtendPressureConditionProcess
    KRATOS_CLASS_POINTER_DEFINITION(ExtendPressureConditionProcess);

    // Constructor
    ExtendPressureConditionProcess(ModelPart &r_model_part);

    // Destructor
    ~ExtendPressureConditionProcess() override = default;

    void operator()() { Execute(); }

    void Execute() override;

    void RemovePreviousLineLoads();

    void CreateNewConditions();

    void GenerateLineLoads2Nodes(
        const int NonWetLocalIdNode,
        const int PressureId,
        int& rMaximumConditionId,
        ModelPart::ElementsContainerType::ptr_iterator itElem);

    void GenerateLineLoads3Nodes(
        const int PressureId,
        int& rMaximumConditionId,
        ModelPart::ElementsContainerType::ptr_iterator itElem);

    void GeneratePressureLoads3WetNodes(
        const int NonWetLocalIdNode,
        const int PressureId,
        int& rMaximumConditionId,
        ModelPart::ElementsContainerType::ptr_iterator itElem);

    void CreatePressureLoads(
        const int Id1,
        const int Id2,
        const int Id3,
        ModelPart::ElementsContainerType::ptr_iterator itElem,
        ModelPart& rSubModelPart,
        ModelPart::PropertiesType::Pointer pProperties,
        int& rMaximumConditionId);

    void CreateLineLoads(
        const int Id1,
        const int Id2,
        ModelPart::ElementsContainerType::ptr_iterator itElem,
        ModelPart& rSubModelPart,
        ModelPart::PropertiesType::Pointer pProperties,
        int& rMaximumConditionId);

    void GeneratePressureLoads4WetNodes(
        const int PressureId,
        int& rMaximumConditionId,
        ModelPart::ElementsContainerType::ptr_iterator itElem);

    void GetPressureId(
        ModelPart::ElementsContainerType::ptr_iterator itElem,
        int& rPressureId);

    void GetMaximumConditionIdOnSubmodelPart(
        int& MaximumConditionId);

    void CalculateNumberOfElementsOnNodes();
    bool CheckIfHasConditionId(const IndexType Id);

    void ResetFlagOnElements();

protected:
    // Member Variables
    ModelPart& mrModelPart;
    std::string mPressureName;

};  // Class

}  // namespace Kratos
#endif /* KRATOS_EXTEND_PRESSURE_PROCESS defined */