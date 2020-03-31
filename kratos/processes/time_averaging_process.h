//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//                   Riccardo Tosi
//

#if !defined(KRATOS_TIME_AVERAGING_PROCESS_H_INCLUDED)
#define KRATOS_TIME_AVERAGING_PROCESS_H_INCLUDED

// System includes
#include <functional>
#include <string>
#include <vector>

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"

namespace Kratos
{
///@addtogroup Kratos Core
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) TimeAveragingProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;
    using ElementType = ModelPart::ElementType;

    /// Pointer definition of TimeAveragingProcess
    KRATOS_CLASS_POINTER_DEFINITION(TimeAveragingProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    TimeAveragingProcess(Model& rModel, Parameters rParameters);

    /// Destructor.
    ~TimeAveragingProcess() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void ExecuteInitialize() override;

    void ExecuteFinalizeSolutionStep() override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    Model& mrModel;
    Parameters mrParameters;

    std::string mModelPartName;
    std::string mIntegrationControlVariableName;
    enum TimeAveragingContainers
    {
        NodalHistorical,
        NodalNonHistorical,
        ElementalNonHistorical
    };
    TimeAveragingContainers mTimeAveragingContainer;
    enum TimeAveragingMethods
    {
        Average,
        RootMeanSquare
    };
    TimeAveragingMethods mTimeAveragingMethod;

    std::vector<std::string> mVariableNamesList;
    std::vector<std::string> mAveragedVariableNamesList;

    int mEchoLevel;

    double mCurrentTime;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    bool IsIntegrationStep() const;

    template <typename TDataType>
    void CalculateTimeIntegratedNodalHistoricalQuantity(ModelPart::NodesContainerType& rNodes,
                                                        const Variable<TDataType>& rVariable,
                                                        const Variable<TDataType>& rAveragedVariable,
                                                        const double DeltaTime) const;

    template <typename TDataType, typename TContainerType>
    void CalculateTimeIntegratedNonHistoricalQuantity(
        TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        const Variable<TDataType>& rAveragedVariable,
        const double DeltaTime) const;

    template <typename TDataType>
    void AverageMethod(const TDataType& rTemporalVariable,
                       TDataType& rAveragedVariable,
                       const double DeltaTime) const;

    template <typename TDataType>
    void RootMeanSquareMethod(const TDataType& rTemporalVariable,
                              TDataType& rAveragedVariable,
                              const double DeltaTime) const;

    template <typename TDataType>
    std::function<void(const TDataType&, TDataType&, const double)> GetTimeAveragingMethod() const;

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TimeAveragingProcess& operator=(TimeAveragingProcess const& rOther);

    /// Copy constructor.
    TimeAveragingProcess(TimeAveragingProcess const& rOther);

    ///@}

}; // Class TimeAveragingProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const TimeAveragingProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_TIME_AVERAGING_PROCESS_H_INCLUDED defined