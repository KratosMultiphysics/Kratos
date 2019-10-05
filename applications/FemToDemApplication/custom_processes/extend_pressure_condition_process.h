//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license:
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
    ExtendPressureConditionProcess(ModelPart& r_model_part);

    // Destructor
    ~ExtendPressureConditionProcess() override = default;

    void operator()() { Execute(); }

    /**
     * @brief Regenerates the prssure load according to ne new boundary
     */
    void Execute() override;

    /**
     * @brief This erases all the line loads 
     */
    void RemovePreviousLineLoads();

    /**
     * @brief This generates the new line/surface loads
     */
    void CreateNewConditions();

    /**
     * @brief Creates line loads for an element if there are 2 wet nodes
     * @param NonWetLocalIdNode the local id of the non wet node
     * @param PressureId the pressure id
     * @param rMaximumConditionId the maximum condition id just to not repeat
     * @param itElem the element analysed
     */
    void GenerateLineLoads2Nodes(
        const int NonWetLocalIdNode,
        const int PressureId,
        int& rMaximumConditionId,
        ModelPart::ElementsContainerType::ptr_iterator itElem);

    /**
     * @brief Creates line loads for an element if there are 3 wet nodes
     * @param PressureId the pressure id
     * @param rMaximumConditionId the maximum condition id just to not repeat
     * @param itElem the element analysed
     */
    void GenerateLineLoads3Nodes(
        const int PressureId,
        int& rMaximumConditionId,
        ModelPart::ElementsContainerType::ptr_iterator itElem);

    /**
     * @brief Creates surface loads for an element if there are 3 wet nodes
     * @param PressureId the pressure id
     * @param rMaximumConditionId the maximum condition id just to not repeat
     * @param itElem the element analysed
     */
    void GeneratePressureLoads3WetNodes(
        const int NonWetLocalIdNode,
        const int PressureId,
        int& rMaximumConditionId,
        ModelPart::ElementsContainerType::ptr_iterator itElem);

    /**
     * @brief Creates surface loads 
     * @param Id1 The node 1 Id
     * @param Id2 The node 2 Id
     * @param Id1 The node 3 Id
     * @param itElem the element analysed
     * @param rMaximumConditionId the maximum condition id just to not repeat
     * @param itElem the element analysed
     * @param rSubModelPart the submodel part to add the condition
     * @param pProperties the properties of the condition
     */
    void CreatePressureLoads(
        const int Id1,
        const int Id2,
        const int Id3,
        ModelPart::ElementsContainerType::ptr_iterator itElem,
        ModelPart& rSubModelPart,
        ModelPart::PropertiesType::Pointer pProperties,
        int& rMaximumConditionId);

    /**
     * @brief Creates surface loads 
     * @param Id1 The node 1 Id
     * @param Id2 The node 2 Id
     * @param itElem the element analysed
     * @param rMaximumConditionId the maximum condition id just to not repeat
     * @param itElem the element analysed
     * @param rSubModelPart the submodel part to add the condition
     * @param pProperties the properties of the condition
     */
    void CreateLineLoads(
        const int Id1,
        const int Id2,
        ModelPart::ElementsContainerType::ptr_iterator itElem,
        ModelPart& rSubModelPart,
        ModelPart::PropertiesType::Pointer pProperties,
        int& rMaximumConditionId);

    /**
     * @brief Creates surface loads for an element if there are 4 wet nodes
     * @param PressureId the pressure id
     * @param rMaximumConditionId the maximum condition id just to not repeat
     * @param itElem the element analysed
     */
    void GeneratePressureLoads4WetNodes(
        const int PressureId,
        int& rMaximumConditionId,
        ModelPart::ElementsContainerType::ptr_iterator itElem);

    /**
     * @brief returns the maximum condition Id
     * @param PressureId the pressure id
     * @param itElem the element analysed
     */
    void GetPressureId(
        ModelPart::ElementsContainerType::ptr_iterator itElem,
        int& rPressureId);

    /**
     * @brief returns the maximum condition Id inside the submodel
     * @param rMaximumConditionId the maximum condition id just to not repeat
     */
    void GetMaximumConditionIdOnSubmodelPart(
        int& MaximumConditionId);

    /**
     * @brief calculates the number of elements that share a node
     */
    void CalculateNumberOfElementsOnNodes();

    /**
     * @brief checks if has a condition with a certain Id
     */
    bool CheckIfHasConditionId(const IndexType Id);

    /**
     * @brief Sets the Flag on all the elements to the next step
     */
    void ResetFlagOnElements();

    /**
     * @brief fills the mpPropertiesVector and mPropertiesId
     */
    void SavePreviousProperties();

protected:
    // Member Variables
    ModelPart& mrModelPart;
    std::string mPressureName;
    std::vector<ModelPart::PropertiesType::Pointer> mpPropertiesVector;

};  // Class

}  // namespace Kratos
#endif /* KRATOS_EXTEND_PRESSURE_PROCESS defined */