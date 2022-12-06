//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aron Noordam
//

#if !defined(KRATOS_SET_MOVING_LOAD_PROCESS_H_INCLUDED )
#define  KRATOS_SET_MOVING_LOAD_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

namespace Kratos {

///@name Kratos Globals
///@{

///@name Kratos Classes
///@{

/// Process to create the animated Eigenvectors
/** This process distributes a load on surface load conditions belonging to a modelpart.
 *  The load is distributed according to the surface area.
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SetMovingLoadProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

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

    void ExecuteInitialize() override;

    void ExecuteInitializeSolutionStep() override;

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

    std::vector<Condition> mSortedConditions;
    std::vector<bool> mIsCondReversedVector;
    array_1d<double,3> mLoad;
    array_1d<double, 3> mOriginPoint;
    array_1d<int,3> mDirection;
    double mLoadVelocity;
    double mCurrentDistance;

    ///@}
    ///@name Private Operations
    ///@{
    static std::vector<int> FindNonRepeatingIndices(std::vector<int> arr);

    std::vector<Condition> FindEndConditions();

    std::vector<Condition> SortConditions(ModelPart::ConditionsContainerType& unsorted_conditions, Condition& first_condition);

    static bool SortConditionPoints(Condition& rCondition, vector<int> direction);

    static Condition& GetFirstCondition(Point first_point, Point second_point, vector<int> direction, std::vector<Condition>& end_conditions);

    static Condition& GetFirstConditionFromCoord(double first_coord, double second_coord, int direction, std::vector<Condition>& end_conditions);

    static bool SwapPoints(double first_coord, double second_coord, int direction);

    void InitializeDistanceLoadInSortedVector();
    ///@}

}; // Class SetMovingLoadProcess

///@}

}  // namespace Kratos.

#endif // KRATOS_SET_MOVING_LOAD_PROCESS_H_INCLUDED  defined
