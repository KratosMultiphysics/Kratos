//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
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

///@name Kratos Classes
///@{

/**
 * @brief Line output process
 *
 * This process outputs double/array_1d<double, 3> variables in historical/non-historical nodal
 * data value containers to files. Values will be sampled along a lint from starting point to end point
 * with given number of sampling points. All the sampled values will be written into one file per time step.
 *
 * Output frequency can be controlled by output control parameters. But, if a string is used as output control
 * variable, then it is assumed every step is output step.
 *
 * Output is given in Comma Seperated Values (CSV) format.
 *
 */
class KRATOS_API(RANS_APPLICATION) RansLineOutputProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansLineOutputProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansLineOutputProcess);

    using NodeType = ModelPart::NodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansLineOutputProcess(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~RansLineOutputProcess() override = default;

    /// Assignment operator.
    RansLineOutputProcess& operator=(RansLineOutputProcess const& rOther) = delete;

    /// Copy constructor.
    RansLineOutputProcess(RansLineOutputProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void ExecuteInitialize() override;

    void ExecuteFinalizeSolutionStep() override;

    const Parameters GetDefaultParameters() const override;

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

private:
    ///@name Member Variables
    ///@{

    Model& mrModel;

    std::string mModelPartName;
    std::vector<std::string> mVariableNames;
    bool mWriteHeader;
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
    ///@name Private Operations
    ///@{

    bool IsOutputStep();

    void WriteOutputFile();

    void WriteOutputFileHeader(
        std::ofstream& rOutputFileStream) const;

    template<class TDataType>
    TDataType InterpolateVariable(
        const Variable<TDataType>& rVariable,
        const int SamplingIndex) const;

    std::string GetOutputFileName() const;

    ///@}

}; // Class RansLineOutputProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansLineOutputProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_LINE_OUTPUT_PROCESS_H_INCLUDED defined
