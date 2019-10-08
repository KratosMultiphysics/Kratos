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

#if !defined(KRATOS_REGENERATE_PFEM_PRESSURE_CONDITIONS_PROCESS)
#define KRATOS_REGENERATE_PFEM_PRESSURE_CONDITIONS_PROCESS

#include "processes/process.h"
#include "fem_to_dem_application_variables.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"

namespace Kratos {

typedef std::size_t SizeType;
    typedef Node <3> NodeType;
    typedef Properties PropertiesType;
    typedef Element ElementType;
    typedef Condition ConditionType;
    typedef Mesh<NodeType, PropertiesType, ElementType, ConditionType> MeshType;
    typedef PointerVector<MeshType> MeshesContainerType;
    typedef MeshType::ElementIterator ElementIterator;

/** 
 * @class RegeneratePfemPressureConditionsProcess
 * @ingroup FemToDemApplication 
 * @brief Regenerates the pressure conditions for the PFEM coupling
 * @details when several elements are removed this methods generates the line loads
 * in order to adapt to the new geometry
 * @author Alejandro Cornejo
 */
template <SizeType TDim = 3>
class RegeneratePfemPressureConditionsProcess : public Process 
{


public:
    /// Pointer definition of RegeneratePfemPressureConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(RegeneratePfemPressureConditionsProcess);

    // Constructor
    RegeneratePfemPressureConditionsProcess(ModelPart& r_model_part);

    // Destructor
    ~RegeneratePfemPressureConditionsProcess() override = default;

    void operator()() { Execute(); }

    /**
     * @brief Regenerates the prssure load according to ne new boundary
     */
    void Execute() override;

    /**
     * @brief returns the maximum condition Id inside the submodel
     * @param rMaximumConditionId the maximum condition id just to not repeat
     */
    void GetMaximumConditionIdOnSubmodelPart(
        int& MaximumConditionId);

    /**
     * @brief Creates the Pressure loads of the gometry
     */
    void CreatePressureLoads(
        const int Id1,
        const int Id2,
        const int Id3,
        ElementIterator itElem,
        ModelPart &rSubModelPart,
        ModelPart::PropertiesType::Pointer pProperties,
        int &rMaximumConditionId);

    /**
     * @brief Resets the Flag SMOOTHING used to stop iterations
     */
    void ResetFlagOnElements();

    /**
     * @brief Erase the old pressure conditions
     */
    void RemovePreviousPressureLoads();

    /**
     * @brief Create the new set of pfem pressure conditions
     */
    void CreateNewConditions();

    /**
     * @brief Creates a pressure condition when 3 nodes are wet
     */
    void GeneratePressureLoads3WetNodes(
        const int NonWetLocalIdNode,
        int &rMaximumConditionId,
        ElementIterator itElem);

    /**
     * @brief Creates a pressure condition when 4 nodes are wet
     */
    void GeneratePressureLoads4WetNodes(
        int &rMaximumConditionId,
        ElementIterator itElem);

    /**
     * @brief Creates line loads for an element if there are 2 wet nodes
     * @param NonWetLocalIdNode the local id of the non wet node
     * @param PressureId the pressure id
     * @param rMaximumConditionId the maximum condition id just to not repeat
     * @param itElem the element analysed
     */
    void GenerateLineLoads2Nodes(
        const int NonWetLocalIdNode,
        int& rMaximumConditionId,
        ElementIterator itElem);

    /**
     * @brief Creates line loads for an element if there are 3 wet nodes
     * @param PressureId the pressure id
     * @param rMaximumConditionId the maximum condition id just to not repeat
     * @param itElem the element analysed
     */
    void GenerateLineLoads3Nodes(
        int& rMaximumConditionId,
        ElementIterator itElem);

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
        ElementIterator itElem,
        ModelPart& rSubModelPart,
        ModelPart::PropertiesType::Pointer pProperties,
        int& rMaximumConditionId);

protected:
    // Member Variables
    ModelPart& mrModelPart;
    
}; // Class
}  // namespace Kratos
#endif /* KRATOS_REGENERATE_PFEM_PRESSURE_CONDITIONS_PROCESS defined */