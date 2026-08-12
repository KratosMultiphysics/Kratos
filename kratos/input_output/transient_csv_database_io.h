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

#pragma once

// System includes
#include <string>
#include <tuple>
#include <variant>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/data_communicator.h"
#include "input_output/transient_database_io.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) TransientCSVDatabaseIO: public TransientDatabaseIO
{
public:
    ///@name Type definitions
    ///@{

    using BaseType = TransientDatabaseIO;

    using IndexType = std::size_t;

    using ValueType = std::variant<std::monostate, bool, int, double, std::string>;

    KRATOS_CLASS_POINTER_DEFINITION(TransientCSVDatabaseIO);

    ///@}
    ///@name Life cycle
    ///@{

    TransientCSVDatabaseIO(
        const std::string& rFilename,
        const bool WriteKratosVersion,
        const bool WriteTimeStamp,
        const IndexType EchoLevel,
        Parameters FormatParameters,
        const DataCommunicator& rDataCommunicator);

    ~TransientCSVDatabaseIO() override = default;

    ///@}
    ///@name Public operations
    ///@{

    void Initialize() override;

    void Finalize() override;

    void SetHeaderInformation(const std::string& rHeaderInformation);

    ///@}
    ///@name Input / output
    ///@{

    void Read(
        bool& rValue,
        const int Step,
        const std::string& rKey) override;

    void Read(
        int& rValue,
        const int Step,
        const std::string& rKey) override;

    void Read(
        double& rValue,
        const int Step,
        const std::string& rKey) override;

    void Read(
        std::string& rValue,
        const int Step,
        const std::string& rKey) override;

    void Write(
        const bool Value,
        const int Step,
        const std::string& rKey) override;

    void Write(
        const int Value,
        const int Step,
        const std::string& rKey) override;

    void Write(
        const double Value,
        const int Step,
        const std::string& rKey) override;

    void Write(
        const std::string& rValue,
        const int Step,
        const std::string& rKey) override;

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

private:
    ///@name Private class definitions
    ///@{

    enum FileAccessMode
    {
        NOT_INITIALIZED = 0,
        READ_ONLY = 1,
        WRITE_ONLY = 2
    };

    struct FormatSettings
    {
        IndexType mIntLength;
        IndexType mFloatPrecision;
        IndexType mStringLength;
        std::vector<std::string> mBooleanValues;
    };

    class Column
    {
    public:
        ///@name Life cycle
        ///@{

        Column(
            const std::string& rHeaderName,
            const ValueType& rValue,
            std::shared_ptr<FormatSettings> pFormatSettings);

        ///@}
        ///@name Public operations
        ///@{

        std::string GetHeader() const;

        std::string GetFormattedHeader() const;

        std::string GetFormattedValue(const ValueType& rValue) const;

        ///@}

    private:
        ///@name Private member variables
        ///@{

        std::string mHeader;

        IndexType mColumnWidth;

        std::string mFormatString;

        ValueType mDummyValue;

        std::shared_ptr<FormatSettings> mpFormatSettings;

        ///@}
    };

    ///@}
    ///@name Private common member variables
    ///@{

    IndexType mEchoLevel;

    std::string mFileName;

    const DataCommunicator& mrDataCommunicator;

    FileAccessMode mFileAccessMode;

    ///@}
    ///@name Private member variables for writing
    ///@{

    std::vector<std::pair<Column, ValueType>> mWritingData;

    std::shared_ptr<FormatSettings> mpFormatSettings;

    std::string mHeaderInformation;

    bool mWriteKratosVersion;

    bool mWriteTimeStamp;

    int mCurrentStep = -1;

    int mLastWrittenStep = -1;

    ///@}
    ///@name Private member variables for reading
    ///@{

    std::vector<std::pair<std::string, std::variant<std::vector<bool>, std::vector<int>, std::vector<double>, std::vector<std::string>>>> mReadData;

    ///@}
    ///@name Private operations for both reading and writing
    ///@{

    bool IsMasterRank() const;

    ///@}
    ///@name Private operations for reading
    ///@{

    template<class TDataType>
    void GenericRead(
        TDataType& rValue,
        const int Step,
        const std::string& rKey);

    void ReadDataFromCSVFile();

    ///@}
    ///@name Private operations for writing
    ///@{

    template<class TDataType>
    void GenericWrite(
        const TDataType& rValue,
        const int Step,
        const std::string& rKey);


    void WriteHeaders(std::ofstream& rOutputFile) const;

    void WriteData(std::ofstream& rOutputFile);

    ///@}
};

///@}

} // namespace Kratos