// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Aron Noordam
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

#include "utilities/function_parser_utility.h"
#include "utilities/mortar_utilities.h"
namespace Kratos {

///@name Kratos Globals
///@{

///@name Kratos Classes
///@{


/**
 * @class SetMovingLoadProcess
 * @ingroup StructuralMechanicsApplication
 * @brief Process to set the moving load 
 * @details This process sorts the moving load conditions, it calculates the value and velocity of the moving load. And it places the load on the right position per
 * solution step.
 * @author Aron Noordam
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SetMovingLoadProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{
    ///
    

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of SetMovingLoadProcess
    KRATOS_CLASS_POINTER_DEFINITION(SetMovingLoadProcess);

    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    SetMovingLoadProcess(ModelPart& rModelPart,
                                  Parameters Parameters);

    ///@}
    ///@name Operations
    ///@{

    /**
     * \brief  Initializes the set moving load process. Check if load functions and a velocity function are present in the parameters.
     * Sort vector of conditions, and find the start position of the moving load, within the conditions vector.
     */
    void ExecuteInitialize() override;

    /**
     * \brief Initialize solution step. Calculate the load based on the load functions if present, else retrieve the load from the input parameters.
     * Loop over the conditions and find, on which condition the load is located. Then set the load on the condition element, if the load is located
     * within the element. If the moving load is not located on the condition element, set the load to zero.
     */
    void ExecuteInitializeSolutionStep() override;

    /**
     * \brief Finalizes solution step. Sets load velocity based on load velocity function if present, else load velocity is retrieved from the input values.
     * Then move the load based on the current position and the load velocity.
     */
    void ExecuteFinalizeSolutionStep() override;


	///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override {
        return "SetMovingLoadProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << "SetMovingLoadProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    Parameters mParameters;

    // vector of sorted conditions
    std::vector<Condition> mSortedConditions;

    // vector of booleans which indicate if the nodes in the condition are in reversed order with respect to the direction of the moving load
    std::vector<bool> mIsCondReversedVector;

    // point of origin of the moving load
    array_1d<double, 3> mOriginPoint;

    // array of directions, of the moving load. This array indicates if the load moves in positive or negative direction with respect to the global
    // axis
    array_1d<int,3> mDirection;

    // distance of the load, with respect to the first node of the first condition in the sorted conditions vector
    double mCurrentDistance;

    // bool which indicates if a load function is used to determine the value of the load
    bool mUseLoadFunction;

    // bool which indicates of velocity function is used to determine the value of the load velocity
    bool mUseVelocityFunction;

    // vector of load functions
    std::vector<BasicGenericFunctionUtility> mLoadFunctions;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * \brief Finds indices of the nodes within the conditions model part, which do not share are condition within the model part.
     * \param IndicesVector vector of all the node indices, which are found in the condition elements.
     * \return vector of non repeating node indices within the model part
     */
    static std::vector<IndexType> FindNonRepeatingIndices(const std::vector<IndexType> IndicesVector);

    /**
     * \brief Finds condition elements which are at the spatial ends of the conditions vector. This function checks which condition's points are not repeated
     * in this conditions model part. If nodes are not repeated, it means that the condition is at one of the spatial ends.
     * \return Vector of the two end conditions.
     */
    std::vector<Condition> FindEndConditions();

    /**
     * \brief Sorts conditions vector. Sorts conditions in the direction of the moving load. The sorting starts from the first condition. Repeated nodes
     * ares searched in the conditions vector. If a node is repeated, it means the condition is connected to the previous condition. 
     * \param rUnsortedConditions container of the unsorted conditions within the model part
     * \param rFirstCondition first condition in the to be sorted conditions vector
     * \return vector of sorted conditions
     */
    std::vector<Condition> SortConditions(ModelPart::ConditionsContainerType& rUnsortedConditions, Condition& rFirstCondition);

    /**
     * \brief Check if conditions points are to be flipped in the direction of the moving load. The points are sorted based on the x-coordinates; if the x-coordinates are equal,
     * sorting is done based on the y-coordinates; if y-coordinates are equal, sorting is done based on the z-coordinates
     * \param rCondition reference of the current condition
     * \param Direction direction of the moving load
     * \return bool which indicates if the condition's points are to be flipped
     */
    static bool IsConditionReversed(const Condition& rCondition, const array_1d<int, 3> Direction);

    /**
     * \brief Finds the first condition within the model part, based on the direction of the moving load.
     * \param  FirstPoint, one of the end points of the condition model part
     * \param  SecondPoint, one of the end points of the condition model part
     * \param  Direction, vector of the direction of the moving load, positive values indicate positive direction with respect to the global axis
     * \param  rEndConditions vector of the two end conditions within the model part
     * \return the first condition within the model part with respect to the direction of the moving load
     */
    static Condition& GetFirstCondition(const Point FirstPoint, const Point SecondPoint, const array_1d<int, 3> Direction, std::vector<Condition>& rEndConditions);

    /**
     * \brief Finds the first condition within the model part, based on the direction of the moving load. In this function, a single direction is checked, either
     * x,y or z.
     * \param FirstCoord one of the end x,y or z coordinate of the condition model part
     * \param SecondCoord one of the end x,y or z coordinate of the condition model part
     * \param Direction direction of the moving load, positive or negative with respect to the x,y or z axis
     * \param EndConditions vector of the two end conditions within the model part
     * \return the first condition within the model part with respect to the direction of the moving load
     */
    static Condition& GetFirstConditionFromCoord(const double FirstCoord, const double SecondCoord, const int Direction, std::vector<Condition>& EndConditions);

    /**
     * \brief Checks if the points within a condition are to be swapped, based on the direction of the moving load
     * \param FirstCoord  x,y or z coordinate of the first point within the condition element
     * \param SecondCoord x,y or z coordinate of the second point within the condition element
     * \param Direction  direction of the moving load, positive or negative with respect to the x,y or z axis
     * \return bool which indicates if the condition's points are to be swapped
     */
    static bool IsSwapPoints(const double FirstCoord, const double SecondCoord, const int Direction);

    /**
     * \brief Initializes the global distance of the load based on the origin point within the sorted condition vector, starting from the first
     * coordinate of the first condition in the sorted conditions vector
     */
    void InitializeDistanceLoadInSortedVector();

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

    ///@}

}; // Class SetMovingLoadProcess

///@}

}  // namespace Kratos.
