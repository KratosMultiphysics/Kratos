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
//

#if !defined(KRATOS_RANS_LINE_OUTPUT_PROCESS_H_INCLUDED)
#define KRATOS_RANS_LINE_OUTPUT_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"

namespace Kratos
{
///@addtogroup RANSApplication
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

class KRATOS_API(RANS_APPLICATION) RansLineOutputProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = ModelPart::NodeType;

    /// Pointer definition of RansLineOutputProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansLineOutputProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansLineOutputProcess(Model& rModel, Parameters rParameters);

    /// Destructor.
    ~RansLineOutputProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void ExecuteInitialize() override;

    void Execute() override;

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
    std::vector<std::string> mVariableNames;
    array_1d<double, 3> mStartPoint;
    array_1d<double, 3> mEndPoint;

    double mOutputStepInterval;
    double mCurrentStepCount = 0.0;
    double mPreviousStepValue;
    std::string mOutputFileName;
    std::string mOutputStepControlVariableName;

    bool mIsHistoricalValue;

    int mNumberOfSamplingPoints;
    std::vector<double> mSamplingPoints;
    std::vector<int> mSamplingPointElementIds;
    std::vector<Vector> mSamplingPointElementShapeFunctions;

    std::vector<int> mSamplePointLocalIndexList;
    std::vector<std::vector<int>> mSamplePointLocalIndexListMaster;
    std::vector<int> mFoundGlobalPoints;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    bool IsOutputStep();

    void WriteOutputFile();

    void WriteOutputFileHeader(std::ofstream& rOutputFileStream) const;

    double InterpolateVariable(const Variable<double>& rVariable, const int SamplingIndex) const;

    array_1d<double, 3> InterpolateVariable(const Variable<array_1d<double, 3>>& rVariable,
                                            const int SamplingIndex) const;

    std::string GetOutputFileName() const;

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
    RansLineOutputProcess& operator=(RansLineOutputProcess const& rOther);

    /// Copy constructor.
    RansLineOutputProcess(RansLineOutputProcess const& rOther);

    ///@}

}; // Class RansLineOutputProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const RansLineOutputProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_LINE_OUTPUT_PROCESS_H_INCLUDED defined
