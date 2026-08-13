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
#include "includes/data_communicator.h"
#include "input_output/tabular_database_io.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) CSVDatabaseIO: public TabularDatabaseIO
{
public:
    ///@name Type definitions
    ///@{

    using BaseType = TabularDatabaseIO;

    using IndexType = std::size_t;

    using ValueType = std::variant<std::monostate, bool, int, double, std::string>;

    KRATOS_CLASS_POINTER_DEFINITION(CSVDatabaseIO);

    ///@}
    ///@name Life cycle
    ///@{

    /**
     * @brief Construct a new CSVDatabaseIO object for reading csv data file
     * @details This constructor is used to construct a CSV database IO for
     *          reading only. All the method related to writing will throw errors.
     *
     *          The first column should have the column header which is equal to the @p rRowIdName
     *          having integer values in the column. This first column is always used to identify the row
     *          when reading data.
     *
     * @param rInputFileName        Input CSV filename for reading.
     * @param rRowIdName            Column header for row identifier.
     * @param rDataCommunicator     Data communicator.
     * @param rBooleanTrueValue     The string to be checked if the boolean column values are true. If the boolean column has some other value then it will be treated as false.
     * @param EchoLevel             Echo level.
     */
    explicit CSVDatabaseIO(
        const std::string& rInputFileName,
        const DataCommunicator& rDataCommunicator,
        const std::string& rRowIdName = "STEP",
        const std::string& rBooleanTrueValue = "1",
        const IndexType EchoLevel = 0);

    /**
     * @brief Construct a new CSVDatabaseIO object for writing data to a CSV file
     * @details This constructor is used to construct a CSV database IO for
     *          writing only. All the methods related to reading will throw errors.
     *
     *          The first column will be with the column header name @p rRowIdName
     *          having integer values in the column.
     *
     *          @p rTitle , Kratos version (if @p WriteKratosVersion is true), time stamp (if @p WriteTimeStamp is true)
     *          @p rHeaderInformation will be written to the csv file at the construction of the object.
     *
     * @param rOutputFileName           Output csv file name.
     * @param rDataCommunicator         Data communicator.
     * @param rTitle                    Title of the CSV file.
     * @param rHeaderInformation        Header information.
     * @param rRowIdName                Column header for row identifier.
     * @param WriteKratosVersion        If true, then the kratos version will be written in the header.
     * @param WriteTimeStamp            If true, then the construction time at the header and destruction time at the footer will be written.
     * @param IntLength                 Length of the ints (maximum of this or length of the column header will be used for formatting the csv data.)
     * @param FloatPrecision            Precision of the floats written to the csv (number of decimal places.)
     * @param StringLength              Maximum length of the strings to be expected in data (maximum of this or length of the column header will be used for formatting the csv data.)
     * @param rBooleanFalseValue        The string to be written if the boolean data is true.
     * @param rBooleanTrueValue         The string to be written if the boolean data is false.
     * @param EchoLevel                 Echo level.
     */
    explicit CSVDatabaseIO(
        const std::string& rOutputFileName,
        const DataCommunicator& rDataCommunicator,
        const std::string& rTitle,
        const std::string& rHeaderInformation,
        const std::string& rRowIdName = "STEP",
        const bool WriteKratosVersion = true,
        const bool WriteTimeStamp = true,
        const IndexType IntLength = 7,
        const IndexType FloatPrecision = 9,
        const IndexType StringLength = 10,
        const std::string& rBooleanFalseValue = "0",
        const std::string& rBooleanTrueValue = "1",
        const IndexType EchoLevel = 0);

    ~CSVDatabaseIO() override;

    ///@}
    ///@name Input / output
    ///@{

    std::vector<int> GetRowIds() const override;

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

        ValueType GetDummyValue() const;

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

    const std::string mFileName;

    const bool mIsReadOnly;

    const DataCommunicator& mrDataCommunicator;

    const IndexType mEchoLevel;

    ///@}
    ///@name Private member variables for writing
    ///@{

    std::vector<std::pair<Column, ValueType>> mWritingData;

    std::shared_ptr<FormatSettings> mpFormatSettings;

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