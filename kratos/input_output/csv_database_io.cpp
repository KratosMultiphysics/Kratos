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

// System includes
#include <sstream>
#include <fstream>
#include <format>
#include <chrono>
#include <string_view>
#include <ranges>

// External includes

// Project includes
#include "includes/kernel.h"

// Include base h
#include "csv_database_io.h"

namespace Kratos
{

namespace CSVDatabaseIOUtils
{

std::string Trim(const std::string& rInput) {
    const IndexType left  = std::distance(rInput.begin(), std::find_if(rInput.begin(), rInput.end(), [](const auto C) { return C != ' '; }));
    const IndexType right = std::distance(rInput.rbegin(), std::find_if(rInput.rbegin(), rInput.rend(), [](const auto C) { return C != ' '; }));
    const IndexType length = rInput.size() - left - right;
    return rInput.substr(left, length);
}

std::string GetType(const CSVDatabaseIO::ValueType& rValue)
{
    KRATOS_TRY

    switch (rValue.index()) {
        case 0 : return "uninitialized";
        case 1 : return "boolean";
        case 2 : return "integer";
        case 3 : return "float";
        case 4 : return "string";
        default: return "Unsupported type";
    }

    KRATOS_CATCH("");
}

} // namespace CSVDatabaseIOUtils

CSVDatabaseIO::Column::Column(
    const std::string& rHeaderName,
    const ValueType& rValue,
    const FormatSettings& rFormatSettings)
    : mHeader(rHeaderName),
      mDummyValue(rValue),
      mrFormatSettings(rFormatSettings)
{
    KRATOS_TRY

    switch (mDummyValue.index()) {
        case 1: {
            mColumnWidth = std::max(rHeaderName.size(), std::max(mrFormatSettings.mBooleanValues[0].size(), mrFormatSettings.mBooleanValues[1].size()));
            mFormatString = "{:>" + std::to_string(mColumnWidth) + "}";
            break;
        }
        case 2: {
            mColumnWidth = std::max(rHeaderName.size(), mrFormatSettings.mIntLength);
            mFormatString = "{: " + std::to_string(mColumnWidth) + "d}";
            break;
        }
        case 3: {
            mColumnWidth = std::max(rHeaderName.size(), mrFormatSettings.mFloatPrecision + 7);
            mFormatString = "{: " + std::to_string(mColumnWidth) + "." + std::to_string(mrFormatSettings.mFloatPrecision) + "e}";
            break;
        }
        case 4: {
            mColumnWidth = std::max(rHeaderName.size(), mrFormatSettings.mStringLength);
            mFormatString = "{:>" + std::to_string(mColumnWidth) + "}";
            break;
        }
        default: {
            KRATOS_ERROR << "The columns always needs a proper value type [ header name = \"" << rHeaderName << "\" ]." ;
        }
    }

    KRATOS_CATCH("");
}

CSVDatabaseIO::ValueType CSVDatabaseIO::Column::GetDummyValue() const
{
    return mDummyValue;
}

std::string CSVDatabaseIO::Column::GetHeader() const
{
    return mHeader;
}

std::string CSVDatabaseIO::Column::GetFormattedHeader() const
{
    return std::vformat("{:>" + std::to_string(mColumnWidth) + "}", std::make_format_args(mHeader));;
}

std::string CSVDatabaseIO::Column::GetFormattedValue(const ValueType& rValue) const
{
    if (mDummyValue.index() == rValue.index()) {
        switch (rValue.index()) {
            case 1 : return std::vformat(mFormatString, std::make_format_args(mrFormatSettings.mBooleanValues[std::get<bool>(rValue)]));
            case 2 : return std::vformat(mFormatString, std::make_format_args(std::get<int>(rValue)));
            case 3 : return std::vformat(mFormatString, std::make_format_args(std::get<double>(rValue)));
            case 4 : return std::vformat(mFormatString, std::make_format_args(std::get<std::string>(rValue)));
            default: return std::vformat("{:>" + std::to_string(mColumnWidth) + "}", std::make_format_args("n/a"));
        }
    } else {
        if (std::holds_alternative<std::monostate>(rValue)) {
            return std::vformat("{:>" + std::to_string(mColumnWidth) + "}", std::make_format_args("n/a"));;
        } else {
            KRATOS_ERROR << "Type mismatch for header \"" << mHeader << "\"."
                         << "\n\t Column type = " << CSVDatabaseIOUtils::GetType(mDummyValue)
                         << "\n\t Given value type = " << CSVDatabaseIOUtils::GetType(rValue);
        }
    }
    return "";
}

CSVDatabaseIO::CSVDatabaseIO(
    const std::string& rInputFileName,
    const DataCommunicator& rDataCommunicator,
    const std::string& rRowIdName,
    const std::string& rBooleanTrueValue,
    const IndexType EchoLevel)
    : mFileName(rInputFileName),
      mIsReadOnly(true),
      mrDataCommunicator(rDataCommunicator),
      mEchoLevel(EchoLevel),
      mRowIdName(rRowIdName),
      mFormatSettings(9, 9, 9, std::vector<std::string>{"0", rBooleanTrueValue}),
      mTitle(""),
      mHeader(""),
      mWriteKratosVersion(false),
      mWriteTimeStamp(false)
{
    KRATOS_TRY

    ReadCSVFile(mFileName);

    KRATOS_CATCH("")
}

CSVDatabaseIO::CSVDatabaseIO(
    const std::string& rOutputFileName,
    const DataCommunicator& rDataCommunicator,
    const std::string& rTitle,
    const std::string& rHeaderInformation,
    const std::string& rRowIdName,
    const bool WriteKratosVersion,
    const bool WriteTimeStamp,
    const IndexType IntLength,
    const IndexType FloatPrecision,
    const IndexType StringLength,
    const std::string& rBooleanFalseValue,
    const std::string& rBooleanTrueValue,
    const IndexType EchoLevel)
    : mFileName(rOutputFileName),
      mIsReadOnly(false),
      mrDataCommunicator(rDataCommunicator),
      mEchoLevel(EchoLevel),
      mRowIdName(rRowIdName),
      mFormatSettings(IntLength, FloatPrecision, StringLength, std::vector<std::string>{rBooleanFalseValue, rBooleanTrueValue}),
      mTitle(rTitle),
      mHeader(rHeaderInformation),
      mWriteKratosVersion(WriteKratosVersion),
      mWriteTimeStamp(WriteTimeStamp)
{
    KRATOS_TRY

    WriteTitleBlock(mFileName);

    KRATOS_CATCH("");
}

CSVDatabaseIO::~CSVDatabaseIO()
{
    if (!mIsReadOnly) {

        std::ofstream output_file(mFileName, std::ios::out | std::ios::app | std::ios::binary);
        WriteData(output_file);

        if (mWriteTimeStamp) {
            output_file << "# End of File - " << std::format("{:%Y-%m-%d %H:%M:%S}", std::chrono::system_clock::now()) << std::endl;
        } else {
            output_file << "# End of File" << std::endl;
        }

        output_file.close();
    }
}

std::vector<int> CSVDatabaseIO::GetRowIds() const
{
    return std::get<std::vector<int>>(mReadData[0].second);
}

void CSVDatabaseIO::Read(
    bool& rValue,
    const int Step,
    const std::string& rKey)
{
    GenericRead(rValue, Step, rKey);
}

void CSVDatabaseIO::Read(
    int& rValue,
    const int Step,
    const std::string& rKey)
{
    GenericRead(rValue, Step, rKey);
}

void CSVDatabaseIO::Read(
    double& rValue,
    const int Step,
    const std::string& rKey)
{
    GenericRead(rValue, Step, rKey);
}

void CSVDatabaseIO::Read(
    std::string& rValue,
    const int Step,
    const std::string& rKey)
{
    GenericRead(rValue, Step, rKey);
}

void CSVDatabaseIO::Write(
    const bool Value,
    const int Step,
    const std::string& rKey)
{
    GenericWrite(Value, Step, rKey);
}

void CSVDatabaseIO::Write(
    const int Value,
    const int Step,
    const std::string& rKey)
{
    GenericWrite(Value, Step, rKey);
}

void CSVDatabaseIO::Write(
    const double Value,
    const int Step,
    const std::string& rKey)
{
    GenericWrite(Value, Step, rKey);
}

void CSVDatabaseIO::Write(
    const std::string& rValue,
    const int Step,
    const std::string& rKey)
{
    GenericWrite(rValue, Step, rKey);
}

void CSVDatabaseIO::ReadCSVFile(const std::string& rFileName)
{
    KRATOS_TRY

    if (!IsMasterRank()) {
        return;
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) <<"Reading CSV data from \"" << rFileName << "\"...\n";

    std::ifstream input_file(rFileName);

    std::string line;
    bool found_column_information_block{false};
    bool found_header_line{false};

    while (std::getline(input_file, line)) {
        if (line[0] == '#') {
            if (line.find("<Column information>") != std::string::npos) {
                found_column_information_block = true;
                KRATOS_INFO_IF("CSVDatabaseIO", mEchoLevel > 0) << "Found the column information block." << std::endl;
                continue;
            }

            if (line.find("End of column information") != std::string::npos) {
                KRATOS_ERROR_IF(mReadData.empty())
                    << "No column information found.\n" << *this << std::endl;
                KRATOS_ERROR_IF(mReadData.front().first != mRowIdName)
                    << "The first column should be \"" << mRowIdName << "\" column with data type \"integer\" [ first column name = " << mReadData.front().first << " ].\n" << *this << std::endl;
                KRATOS_ERROR_IF(!std::holds_alternative<std::vector<int>>(mReadData.front().second))
                    << "The first column should be \"" << mRowIdName << "\" column with data type \"integer\" [ first column name = " << mReadData.front().first << " ].\n" << *this << std::endl;
                found_column_information_block = false;
                continue;
            }

            if (found_column_information_block) {
                // getting the header type information
                const IndexType delimiter = line.find(':');
                const std::string header_name = CSVDatabaseIOUtils::Trim(line.substr(1, delimiter - 1));
                const std::string header_type = CSVDatabaseIOUtils::Trim(line.substr(delimiter + 1, line.size() - 1 - delimiter));

                if (header_type == "boolean") {
                    mReadData.push_back(std::make_pair(header_name, std::vector<bool>{}));
                } else if (header_type == "integer") {
                    mReadData.push_back(std::make_pair(header_name, std::vector<int>{}));
                } else if (header_type == "float") {
                    mReadData.push_back(std::make_pair(header_name, std::vector<double>{}));
                } else if (header_type == "string") {
                    mReadData.push_back(std::make_pair(header_name, std::vector<std::string>{}));
                } else {
                    KRATOS_ERROR << "Unsupported header type = \"" << header_type << "\" for header name = \""
                                << header_name << "\" found. Only supports following types: \n\tboolean\n\tinteger\n\tfloat\n\tstring\n" << *this << std::endl;
                }

                KRATOS_INFO_IF("CSVDatabaseIO", mEchoLevel > 1) << "Found the column header \"" << header_name << "\" with data type \"" << header_type << "\"." << std::endl;
            }
        } else {
            if (CSVDatabaseIOUtils::Trim(line.substr(0, line.find(','))) == mRowIdName) {
                KRATOS_INFO_IF("CSVDatabaseIO", mEchoLevel > 0) << "Found the headers." << std::endl;
                if (mReadData.empty()) {
                    // no column information block was found. Assume all data are floats.
                    for (auto part : line | std::views::split(',')) {
                        const auto& header = CSVDatabaseIOUtils::Trim(std::string{std::string_view(part.begin(), part.end())});
                        mReadData.push_back(std::make_pair(header, std::vector<double>{}));
                    }
                    // replace the first one with a in vector because, first one always should be the step
                    if (!mReadData.empty()) {
                        mReadData[0].second = std::vector<int>{};
                    }

                    KRATOS_ERROR_IF(mReadData.empty())
                        << "No column information found.\n" << *this << std::endl;
                    KRATOS_ERROR_IF(mReadData.front().first != mRowIdName)
                        << "The first column should be \"" << mRowIdName << "\" column with data type \"integer\" [ first column name = " << mReadData.front().first << " ].\n" << *this << std::endl;
                    KRATOS_ERROR_IF(!std::holds_alternative<std::vector<int>>(mReadData.front().second))
                        << "The first column should be \"" << mRowIdName << "\" column with data type \"integer\" [ first column name = " << mReadData.front().first << " ].\n" << *this << std::endl;

                        KRATOS_INFO_IF("CSVDatabaseIO", mEchoLevel > 0) << "No column information block found. Assuming all the columns are of float type.\n";
                }
                found_header_line = true;
                continue;
            }

            if (found_header_line) {
                // split the line by commas (",")
                IndexType column_index = 0;
                for (auto part : line | std::views::split(',')) {
                    if (column_index < mReadData.size()) {
                        const std::string& value = CSVDatabaseIOUtils::Trim(std::string{std::string_view(part.begin(), part.end())});
                        switch (mReadData[column_index].second.index()) {
                            case 0: {
                                if (value != "n/a")
                                    std::get<std::vector<bool>>(mReadData[column_index].second).push_back(value == mFormatSettings.mBooleanValues[1]);
                                else
                                    std::get<std::vector<bool>>(mReadData[column_index].second).push_back(false);
                                break;
                            }
                            case 1: {
                                if (value != "n/a")
                                    std::get<std::vector<int>>(mReadData[column_index].second).push_back(std::stoi(value));
                                else
                                    std::get<std::vector<int>>(mReadData[column_index].second).push_back(0);
                                break;
                            }
                            case 2: {
                                if (value != "n/a")
                                    std::get<std::vector<double>>(mReadData[column_index].second).push_back(std::stod(value));
                                else
                                    std::get<std::vector<double>>(mReadData[column_index].second).push_back(0);
                                break;
                            }
                            case 3: {
                                std::get<std::vector<std::string>>(mReadData[column_index].second).push_back(value);
                                break;
                            }
                        }
                        KRATOS_INFO_IF("CSVDatabaseIO", mEchoLevel > 2)
                            << "Read data \"" << value << "\" for column = \"" << mReadData[column_index].first << "\"." << std::endl;
                    } else {
                        KRATOS_ERROR << "Number of columns mismatch [ line = \"" << line << " ]\n" << *this << std::endl;
                    }
                    ++column_index;
                }

                KRATOS_INFO_IF("CSVDatabaseIO", mEchoLevel > 1) << "Read data for " << mRowIdName << " = " << std::get<std::vector<int>>(mReadData[0].second).back() << "." << std::endl;
            }
        }
    }

    input_file.close();

    KRATOS_CATCH("")
}

void CSVDatabaseIO::WriteTitleBlock(const std::string& rFileName)
{
    KRATOS_TRY

    if (!IsMasterRank()) {
        return;
    }

    mWritingData.clear();

    // Adding the step column
    mWritingData.push_back(std::make_pair(Column(mRowIdName, 0, mFormatSettings), 0));

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) <<"Writing CSV header information to \"" << rFileName << "\"...\n";

    std::ofstream output_file(rFileName, std::ios::out | std::ios::trunc | std::ios::binary);

    output_file << "# ===========================================================" << std::endl;
    output_file << "# ";
    const int left_gap = (59 - mTitle.size()) / 2 - 1;
    if (left_gap > 0) {
        int i;
        for (i = 0; i < left_gap; ++i) {
            output_file << "=";
        }
        output_file << " ";
        for (; i - left_gap < static_cast<int>(mTitle.size()); ++i) {
            output_file << mTitle[i - left_gap];
        }
        output_file << " ";
        for (; i < 59 - 2 ; ++i) {
            output_file << "=";
        }
        output_file << std::endl;
    } else {
        output_file << mTitle << std::endl;
    }

    output_file << "# ===========================================================" << std::endl;

    output_file << "# ------------------- Kratos information --------------------" << std::endl;
    if (mWriteKratosVersion) {
        output_file << "# Kratos version: " << Kernel::Version() << std::endl;
    } else {
        output_file << "# Kratos version: not_given" << std::endl;
    }
    if (mWriteTimeStamp) {
        output_file << "# Timestamp     : " << std::format("{:%Y-%m-%d %H:%M:%S}", std::chrono::system_clock::now()) << std::endl;
    } else {
        output_file << "# Timestamp     : not_specified" << std::endl;
    }
    output_file << "# --------------- End of Kratos information -----------------" << std::endl;

    output_file << mHeader;

    output_file.close();

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << "Header block is written to file \"" << rFileName << "\".\n";

    KRATOS_CATCH("");
}

std::string CSVDatabaseIO::Info() const
{
    return "CSVDatabaseIO";
}

void CSVDatabaseIO::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void CSVDatabaseIO::PrintData(std::ostream& rOStream) const
{
    rOStream << "\t Filename: " << mFileName << std::endl;
    rOStream << "\t File access mode: " << (mIsReadOnly ? "read_only" : "write_only") << "\n";
    rOStream << "\t Row id name: " << mRowIdName << std::endl;

    if (mIsReadOnly) {
        rOStream << "\t Columns: " << std::endl;
        for (const auto& r_pair : mReadData) {
            switch (r_pair.second.index()) {
                case 0: {
                    rOStream << "\t\t " << r_pair.first << ": boolean\n";
                    break;
                }
                case 1: {
                    rOStream << "\t\t " << r_pair.first << ": integer\n";
                    break;
                }
                case 2: {
                    rOStream << "\t\t " << r_pair.first << ": float\n";
                    break;
                }
                case 3: {
                    rOStream << "\t\t " << r_pair.first << ": string\n";
                    break;
                }
            }
        }
    } else {
        rOStream << "\t Title: " << mTitle << std::endl;
        rOStream << "\t Write kratos version: " << (mWriteKratosVersion ? "yes" : "no") << std::endl;
        rOStream << "\t Write time stamp: " << (mWriteTimeStamp ? "yes" : "no") << std::endl;
        rOStream << "\t Current step: " << mCurrentStep << std::endl;
        rOStream << "\t Last written step: " << mLastWrittenStep << std::endl;
        rOStream << "\t Format settings: " << std::endl;
        rOStream << "\t\t Integer length: " << mFormatSettings.mIntLength << std::endl;
        rOStream << "\t\t Float precision: " << mFormatSettings.mFloatPrecision << std::endl;
        rOStream << "\t\t String length: " << mFormatSettings.mStringLength << std::endl;
        rOStream << "\t\t Boolean values: [ \"" << mFormatSettings.mBooleanValues[0] << "\", \"" << mFormatSettings.mBooleanValues[1] << "\" ]" << std::endl;
        rOStream << "\t Columns: " << std::endl;
        for (const auto& r_pair : mWritingData) {
            rOStream << "\t\t " << r_pair.first.GetHeader() << ": " << CSVDatabaseIOUtils::GetType(r_pair.first.GetDummyValue()) << std::endl;
        }
    }
}

void CSVDatabaseIO::WriteData(std::ofstream& rOutputFile)
{
    KRATOS_TRY

    rOutputFile << mWritingData.front().first.GetFormattedValue(mWritingData.front().second);
    for (IndexType i = 1; i < mWritingData.size(); ++i) {
        rOutputFile << ", " << mWritingData[i].first.GetFormattedValue(mWritingData[i].second);
        mWritingData[i].second = ValueType();
    }
    rOutputFile << std::endl;

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 1) << "Data for " << mWritingData[0].first.GetHeader() << " = " << std::get<int>(mWritingData.front().second) << " is written to file \"" << mFileName << "\".\n";

    KRATOS_CATCH("");
}

bool CSVDatabaseIO::IsMasterRank() const
{
    return mrDataCommunicator.Rank() == 0;
}

template<class TDataType>
void CSVDatabaseIO::GenericRead(
    TDataType& rValue,
    const int RowId,
    const std::string& rKey)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mIsReadOnly)
        << "The file \"" << mFileName << "\" is opened for write only access. Hence cannot read.\n" << *this;

    for (const auto& r_pair : mReadData) {
        if (r_pair.first == rKey) {
            KRATOS_ERROR_IF(!std::holds_alternative<std::vector<TDataType>>(r_pair.second))
                << "Type mismatch for key = \"" << rKey << "\".\n" << *this << std::endl;

            const auto& row_ids = std::get<std::vector<int>>(mReadData[0].second);
            const IndexType index = std::distance(row_ids.begin(), std::find(row_ids.begin(), row_ids.end(), RowId));
            KRATOS_ERROR_IF(index == row_ids.size())
                << "The row id = " << RowId << " not found in the list of row ids.\n" << *this;

            rValue = std::get<std::vector<TDataType>>(r_pair.second)[index];
            break;
        }
    }

    KRATOS_CATCH("")
}

void CSVDatabaseIO::WriteHeaders(std::ofstream& rOutputFile) const
{
    rOutputFile << "# ------------------ <Column information> -------------------" << std::endl;
    for (IndexType i = 0; i < mWritingData.size(); ++i) {
        rOutputFile << "#         " << mWritingData[i].first.GetHeader() << ": " << CSVDatabaseIOUtils::GetType(mWritingData[i].second) << std::endl;
    }
    rOutputFile << "# --------------- End of column information -----------------" << std::endl;

    rOutputFile << "# Headers:" << std::endl;

    rOutputFile << mWritingData.front().first.GetFormattedHeader();
    for (IndexType i = 1; i < mWritingData.size(); ++i) {
        rOutputFile << ", " << mWritingData[i].first.GetFormattedHeader();
    }
    rOutputFile << std::endl;

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << "Header names are written to file \"" << mFileName << "\".\n";
}

template<class TDataType>
void CSVDatabaseIO::GenericWrite(
    const TDataType& rValue,
    const int Step,
    const std::string& rKey)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mIsReadOnly)
        << "The file \"" << mFileName << "\" is opened for read only access. Hence cannot write.\n" << *this;

    if (!IsMasterRank()) {
        return;
    }

    if (mCurrentStep == -1) {
        // this is to record the first step
        mCurrentStep = Step;
    }

    if (mCurrentStep == Step) {
        // This can happen in two situations.
        // 1. At the first most step where we collect the headers
        // 2. At a sub sequent step where we need to cross check with the existing headers

        auto p_itr = std::find_if(mWritingData.begin(), mWritingData.end(), [&rKey](const auto rColumnValuePair) { return rKey == rColumnValuePair.first.GetHeader(); });
        if (mLastWrittenStep == -1) {
            // No headers have been written yet. So this is the first step.
            if (mWritingData.end() == p_itr) {
                mWritingData.push_back(std::make_pair(Column(rKey, rValue, mFormatSettings), rValue));
            } else {
                // rKey should be unique
                KRATOS_ERROR << "The keys should be unique."
                             << "\n\t Input    (key, value) = (" << rKey << ", " << rValue << " )"
                             << "\n" << *this;
            }
        } else {
            // Headers are already written. So now we have to check if the rKey is there on the list of headers written
            if (p_itr == mWritingData.end()) {
                std::stringstream msg;
                msg << "The header was not found for key = \"" << rKey << "\". Following headers are available:";
                for (const auto& r_column_value_pair : mWritingData) {
                    msg << "\n\t\"" << r_column_value_pair.first.GetHeader() << "\"";
                }
                KRATOS_ERROR << msg.str() << "\n" << *this;
            }

            KRATOS_ERROR_IF_NOT(std::holds_alternative<TDataType>(p_itr->first.GetDummyValue()))
                << "Type mismatch for the key = \"" << rKey << "\"."
                << "\n\t header (key) type = " << CSVDatabaseIOUtils::GetType(p_itr->first.GetDummyValue())
                << "\n\t given data = " << rValue
                << "\n\t given data type = " << CSVDatabaseIOUtils::GetType(ValueType(rValue)) << std::endl << *this;

            p_itr->second = rValue;
        }
    } else {
        // now information with new Step is given.
        std::ofstream output_file(mFileName, std::ios::out | std::ios::app | std::ios::binary);

        // so we first need to check if the headers are written.
        if (mLastWrittenStep == -1) {
            // writing the headers
            WriteHeaders(output_file);
        }

        WriteData(output_file);
        output_file.close();

        mLastWrittenStep = Step;
        mCurrentStep = Step;
        mWritingData[0].second = Step;
        Write(rValue, Step, rKey);
    }

    KRATOS_CATCH("");
}

// template instantiations
template void CSVDatabaseIO::GenericRead<bool>(bool&, int, const std::string&);
template void CSVDatabaseIO::GenericRead<int>(int&, int, const std::string&);
template void CSVDatabaseIO::GenericRead<double>(double&, int, const std::string&);
template void CSVDatabaseIO::GenericRead<std::string>(std::string&, int, const std::string&);

template void CSVDatabaseIO::GenericWrite<bool>(const bool&, int, const std::string&);
template void CSVDatabaseIO::GenericWrite<int>(const int&, int, const std::string&);
template void CSVDatabaseIO::GenericWrite<double>(const double&, int, const std::string&);
template void CSVDatabaseIO::GenericWrite<std::string>(const std::string&, int, const std::string&);

} // namespace Kratos