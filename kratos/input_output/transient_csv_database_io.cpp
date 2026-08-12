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
#include "transient_csv_database_io.h"

namespace Kratos
{

namespace TransientCSVDatabaseIOUtils
{

std::string Trim(const std::string& rInput) {
    const IndexType left  = std::distance(rInput.begin(), std::find_if(rInput.begin(), rInput.end(), [](const auto C) { return C != ' '; }));
    const IndexType right = std::distance(rInput.rbegin(), std::find_if(rInput.rbegin(), rInput.rend(), [](const auto C) { return C != ' '; }));
    const IndexType length = rInput.size() - left - right;
    return rInput.substr(left, length);
}

std::string GetType(const TransientCSVDatabaseIO::ValueType& rValue)
{
    KRATOS_TRY

    switch (rValue.index()) {
        case 1 : return "boolean";
        case 2 : return "integer";
        case 3 : return "float";
        case 4 : return "string";
        default: return "Unsupported type";
    }

    KRATOS_CATCH("");
}

} // namespace TransientCSVDatabaseIOUtils

TransientCSVDatabaseIO::Column::Column(
    const std::string& rHeaderName,
    const ValueType& rValue,
    std::shared_ptr<FormatSettings> pFormatSettings)
    : mHeader(rHeaderName),
      mDummyValue(rValue),
      mpFormatSettings(pFormatSettings)
{
    KRATOS_TRY

    switch (mDummyValue.index()) {
        case 1: {
            mColumnWidth = std::max(rHeaderName.size(), std::max(pFormatSettings->mBooleanValues[0].size(), pFormatSettings->mBooleanValues[1].size()));
            mFormatString = "{:>" + std::to_string(mColumnWidth) + "}";
            break;
        }
        case 2: {
            mColumnWidth = std::max(rHeaderName.size(), pFormatSettings->mIntLength);
            mFormatString = "{: " + std::to_string(mColumnWidth) + "d}";
            break;
        }
        case 3: {
            mColumnWidth = std::max(rHeaderName.size(), pFormatSettings->mFloatPrecision + 7);
            mFormatString = "{: " + std::to_string(mColumnWidth) + "." + std::to_string(pFormatSettings->mFloatPrecision) + "e}";
            break;
        }
        case 4: {
            mColumnWidth = std::max(rHeaderName.size(), pFormatSettings->mStringLength);
            mFormatString = "{:>" + std::to_string(mColumnWidth) + "}";
            break;
        }
        default: {
            KRATOS_ERROR << "The columns always needs a proper value type [ header name = \"" << rHeaderName << "\" ]." ;
        }
    }

    KRATOS_CATCH("");
}

std::string TransientCSVDatabaseIO::Column::GetHeader() const
{
    return mHeader;
}

std::string TransientCSVDatabaseIO::Column::GetFormattedHeader() const
{
    return std::vformat("{:>" + std::to_string(mColumnWidth) + "}", std::make_format_args(mHeader));;
}

std::string TransientCSVDatabaseIO::Column::GetFormattedValue(const ValueType& rValue) const
{
    if (mDummyValue.index() == rValue.index()) {
        switch (rValue.index()) {
            case 1 : return std::vformat(mFormatString, std::make_format_args(mpFormatSettings->mBooleanValues[std::get<bool>(rValue)]));
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
                         << "\n\t Column type = " << TransientCSVDatabaseIOUtils::GetType(mDummyValue)
                         << "\n\t Given value type = " << TransientCSVDatabaseIOUtils::GetType(rValue);
        }
    }
    return "";
}

TransientCSVDatabaseIO::TransientCSVDatabaseIO(
    const std::string& rFilename,
    const bool WriteKratosVersion,
    const bool WriteTimeStamp,
    const IndexType EchoLevel,
    Parameters FormatParameters,
    const DataCommunicator& rDataCommunicator)
    : mrDataCommunicator(rDataCommunicator),
      mFileAccessMode(NOT_INITIALIZED)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
    {
        "int_length"     : 7,
        "float_precision": 9,
        "bool_values"    : ["0", "1"],
        "string_length"  : 10
    })" );

    FormatParameters.ValidateAndAssignDefaults(default_parameters);

    mFileName = rFilename;
    mWriteKratosVersion = WriteKratosVersion;
    mWriteTimeStamp = WriteTimeStamp;
    mEchoLevel = EchoLevel;

    mpFormatSettings = std::make_shared<FormatSettings>(
                            FormatParameters["int_length"].GetInt(),
                            FormatParameters["float_precision"].GetInt(),
                            FormatParameters["string_length"].GetInt(),
                            FormatParameters["bool_values"].GetStringArray());

    // Adding the step column
    mWritingData.push_back(std::make_pair(Column("STEP", 0, mpFormatSettings), 0));

    KRATOS_CATCH("");
}

void TransientCSVDatabaseIO::SetHeaderInformation(const std::string& rHeaderInformation)
{
    mHeaderInformation = rHeaderInformation;
}

void TransientCSVDatabaseIO::Initialize()
{
}

void TransientCSVDatabaseIO::Finalize()
{
    if (mFileAccessMode == WRITE_ONLY) {
        std::ofstream output_file(mFileName, std::ios::app);
        WriteData(output_file);
        if (mWriteTimeStamp) {
            output_file << "# End of File - " << std::format("{:%Y-%m-%d %H:%M:%S}", std::chrono::system_clock::now()) << "\n";
        } else {
            output_file << "# End of File\n";
        }
        output_file.close();
    }
}

void TransientCSVDatabaseIO::Read(
    bool& rValue,
    const int Step,
    const std::string& rKey)
{
    GenericRead(rValue, Step, rKey);
}

void TransientCSVDatabaseIO::Read(
    int& rValue,
    const int Step,
    const std::string& rKey)
{
    GenericRead(rValue, Step, rKey);
}

void TransientCSVDatabaseIO::Read(
    double& rValue,
    const int Step,
    const std::string& rKey)
{
    GenericRead(rValue, Step, rKey);
}

void TransientCSVDatabaseIO::Read(
    std::string& rValue,
    const int Step,
    const std::string& rKey)
{
    GenericRead(rValue, Step, rKey);
}

void TransientCSVDatabaseIO::Write(
    const bool Value,
    const int Step,
    const std::string& rKey)
{
    GenericWrite(Value, Step, rKey);
}

void TransientCSVDatabaseIO::Write(
    const int Value,
    const int Step,
    const std::string& rKey)
{
    GenericWrite(Value, Step, rKey);
}

void TransientCSVDatabaseIO::Write(
    const double Value,
    const int Step,
    const std::string& rKey)
{
    GenericWrite(Value, Step, rKey);
}

void TransientCSVDatabaseIO::Write(
    const std::string& rValue,
    const int Step,
    const std::string& rKey)
{
    GenericWrite(rValue, Step, rKey);
}

std::string TransientCSVDatabaseIO::Info() const
{
    return "TransientCSVDatabaseIO";
}

void TransientCSVDatabaseIO::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void TransientCSVDatabaseIO::PrintData(std::ostream& rOStream) const
{
    rOStream << "\t Filename: " << mFileName << std::endl;
    rOStream << "\t File access mode: ";
    switch (mFileAccessMode) {
        case NOT_INITIALIZED: {
            rOStream << "not initialized" << std::endl;
            break;
        }
        case READ_ONLY: {
            rOStream << "read only" << std::endl;
            break;
        }
        case WRITE_ONLY: {
            rOStream << "write only" << std::endl;
            break;
        }
    }
    switch (mFileAccessMode) {
        case READ_ONLY: {
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
            break;
        }
        case WRITE_ONLY: {
            rOStream << "\t Write kratos version: " << (mWriteKratosVersion ? "Yes" : "No") << std::endl;
            rOStream << "\t Write time stamp: " << (mWriteTimeStamp ? "Yes" : "No") << std::endl;
            rOStream << "\t Current step: " << mCurrentStep << std::endl;
            rOStream << "\t Last written step: " << mLastWrittenStep << std::endl;
            rOStream << "\t Format settings: " << std::endl;
            rOStream << "\t\t Integer length: " << mpFormatSettings->mIntLength << std::endl;
            rOStream << "\t\t Float precision: " << mpFormatSettings->mFloatPrecision << std::endl;
            rOStream << "\t\t String length: " << mpFormatSettings->mStringLength << std::endl;
            rOStream << "\t\t Boolean values: [ \"" << mpFormatSettings->mBooleanValues[0] << "\", \"" << mpFormatSettings->mBooleanValues[1] << "\" ]" << std::endl;
            rOStream << "\t Columns: " << std::endl;
            for (const auto& r_pair : mWritingData) {
                rOStream << "\t\t " << r_pair.first.GetHeader() << ": " << TransientCSVDatabaseIOUtils::GetType(r_pair.second) << std::endl;
            }
            break;
        }
        default: break;
    }

}

void TransientCSVDatabaseIO::WriteHeaders(std::ofstream& rOutputFile) const
{
    rOutputFile << "# ------------------- Kratos information --------------------\n";
    if (mWriteKratosVersion) {
        rOutputFile << "# Kratos version: " << Kernel::Version() << "\n";
    } else {
        rOutputFile << "# Kratos version: not_given\n";
    }
    if (mWriteTimeStamp) {
        rOutputFile << "# Timestamp     : " << std::format("{:%Y-%m-%d %H:%M:%S}", std::chrono::system_clock::now()) << "\n";
    } else {
        rOutputFile << "# Timestamp     : not_specified\n";
    }
    rOutputFile << "# --------------- End of Kratos information -----------------\n";
    rOutputFile << "# ------------------ <Column information> -------------------\n";
    for (IndexType i = 0; i < mWritingData.size(); ++i) {
        rOutputFile << "#         " << mWritingData[i].first.GetHeader() << ": " << TransientCSVDatabaseIOUtils::GetType(mWritingData[i].second) << "\n";
    }
    rOutputFile << "# --------------- End of column information -----------------\n";
    for (auto part : mHeaderInformation | std::views::split('\n')) {
        rOutputFile << "# " << std::string_view(part.begin(), part.end()) << "\n";
    }

    rOutputFile << mWritingData.front().first.GetFormattedHeader();
    for (IndexType i = 1; i < mWritingData.size(); ++i) {
        rOutputFile << ", " << mWritingData[i].first.GetFormattedHeader();
    }
    rOutputFile << "\n";

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << "Header block is written.\n";
}

void TransientCSVDatabaseIO::WriteData(std::ofstream& rOutputFile)
{
    rOutputFile << mWritingData.front().first.GetFormattedValue(mWritingData.front().second);
    for (IndexType i = 1; i < mWritingData.size(); ++i) {
        rOutputFile << ", " << mWritingData[i].first.GetFormattedValue(mWritingData[i].second);
        // KRATOS_WATCH(mWritingData[i].first)
        mWritingData[i].second = ValueType();
    }
    rOutputFile << "\n";

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 1) << "Data for step " << std::get<int>(mWritingData.front().second) << " is written.\n";
}

bool TransientCSVDatabaseIO::IsMasterRank() const
{
    return mrDataCommunicator.Rank() == 0;
}

void TransientCSVDatabaseIO::ReadDataFromCSVFile()
{
    KRATOS_TRY

    if (!IsMasterRank()) {
        return;
    }

    std::ifstream input_file(mFileName);

    std::string line;
    bool found_column_information_block{false};
    bool found_header_line{false};

    while (std::getline(input_file, line)) {
        if (line[0] == '#') {
            if (line.find("<Column information>") != std::string::npos) {
                found_column_information_block = true;
                KRATOS_INFO_IF("TransientCSVDatabaseIO", mEchoLevel > 0) << "Found the column information block." << std::endl;
                continue;
            }

            if (line.find("End of column information") != std::string::npos) {
                KRATOS_ERROR_IF(mReadData.empty())
                    << "No column information found.\n" << *this << std::endl;
                KRATOS_ERROR_IF(mReadData.front().first != "STEP")
                    << "The first column always should be \"STEP\" column with data type \"integer\"\n" << *this << std::endl;
                KRATOS_ERROR_IF(!std::holds_alternative<std::vector<int>>(mReadData.front().second))
                    << "The first column always should be \"STEP\" column with data type \"integer\"\n" << *this << std::endl;
                found_column_information_block = false;
                continue;
            }

            if (found_column_information_block) {
                // getting the header type information
                const IndexType delimiter = line.find(':');
                const std::string header_name = TransientCSVDatabaseIOUtils::Trim(line.substr(1, delimiter - 1));
                const std::string header_type = TransientCSVDatabaseIOUtils::Trim(line.substr(delimiter + 1, line.size() - 1 - delimiter));

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

                KRATOS_INFO_IF("TransientCSVDatabaseIO", mEchoLevel > 1) << "Found the column header \"" << header_name << "\" with data type \"" << header_type << "\"." << std::endl;
            }
        } else {
            if (TransientCSVDatabaseIOUtils::Trim(line.substr(0, line.find(','))) == "STEP") {
                KRATOS_INFO_IF("TransientCSVDatabaseIO", mEchoLevel > 0) << "Found the headers." << std::endl;
                if (mReadData.empty()) {
                    // no column information block was found. Assume all data are floats.
                    for (auto part : line | std::views::split(',')) {
                        const auto& header = TransientCSVDatabaseIOUtils::Trim(std::string{std::string_view(part.begin(), part.end())});
                        mReadData.push_back(std::make_pair(header, std::vector<double>{}));
                    }
                    // replace the first one with a in vector because, first one always should be the step
                    if (!mReadData.empty()) {
                        mReadData[0].second = std::vector<int>{};
                    }

                    KRATOS_ERROR_IF(mReadData.empty())
                        << "No column information found.\n" << *this << std::endl;
                    KRATOS_ERROR_IF(mReadData.front().first != "STEP")
                        << "The first column always should be \"STEP\" column with data type \"integer\"\n" << *this << std::endl;
                    KRATOS_ERROR_IF(!std::holds_alternative<std::vector<int>>(mReadData.front().second))
                        << "The first column always should be \"STEP\" column with data type \"integer\"\n" << *this << std::endl;

                        KRATOS_INFO_IF("TransientCSVDatabaseIO", mEchoLevel > 0) << "No column information block found. Assuming all the columns are of float type.\n";
                }
                found_header_line = true;
                continue;
            }

            if (found_header_line) {
                // split the line by commas (",")
                IndexType column_index = 0;
                for (auto part : line | std::views::split(',')) {
                    if (column_index < mReadData.size()) {
                        const std::string& value = TransientCSVDatabaseIOUtils::Trim(std::string{std::string_view(part.begin(), part.end())});
                        switch (mReadData[column_index].second.index()) {
                            case 0: {
                                if (value != "n/a")
                                    std::get<std::vector<bool>>(mReadData[column_index].second).push_back(value == mpFormatSettings->mBooleanValues[1]);
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
                        KRATOS_INFO_IF("TransientCSVDatabaseIO", mEchoLevel > 2)
                            << "Read data \"" << value << "\" for column = \"" << mReadData[column_index].first << "\"." << std::endl;
                    } else {
                        KRATOS_ERROR << "Number of columns mismatch [ line = \"" << line << " ]\n" << *this << std::endl;
                    }
                    ++column_index;
                }

                KRATOS_INFO_IF("TransientCSVDatabaseIO", mEchoLevel > 1) << "Read data for STEP = " << std::get<std::vector<int>>(mReadData[0].second).back() << "." << std::endl;
            }
        }
    }

    input_file.close();

    KRATOS_CATCH("")
}

template<class TDataType>
void TransientCSVDatabaseIO::GenericRead(
    TDataType& rValue,
    const int Step,
    const std::string& rKey)
{
    KRATOS_TRY

    switch (mFileAccessMode) {
        case NOT_INITIALIZED: {
                KRATOS_CRITICAL_SECTION
                mFileAccessMode = READ_ONLY;
                break;
        }
        case READ_ONLY: {
            break;
        }
        default: {
            KRATOS_ERROR << "The file can either be opened for reading or writing. Not for both.\n" << *this;
        }
    }

    if (mReadData.empty()) {
        {
            // the GenericRead may be called from parallel threads. But the CSV read needs to be done in serial.
            KRATOS_CRITICAL_SECTION
            ReadDataFromCSVFile();
        }
    }

    for (const auto& r_pair : mReadData) {
        if (r_pair.first == rKey) {
            KRATOS_ERROR_IF(!std::holds_alternative<std::vector<TDataType>>(r_pair.second))
                << "Type mismatch for key = \"" << rKey << "\".\n" << *this << std::endl;
            const auto& steps = std::get<std::vector<int>>(mReadData[0].second);
            const IndexType index = std::distance(steps.begin(), std::find(steps.begin(), steps.end(), Step));
            rValue = std::get<std::vector<TDataType>>(r_pair.second)[index];
            break;
        }
    }

    KRATOS_CATCH("")
}

template<class TDataType>
void TransientCSVDatabaseIO::GenericWrite(
    const TDataType& rValue,
    const int Step,
    const std::string& rKey)
{
    KRATOS_TRY

    switch (mFileAccessMode) {
        case NOT_INITIALIZED:
        case WRITE_ONLY: {
            mFileAccessMode = WRITE_ONLY;
            break;
        }
        default: {
            KRATOS_ERROR << "The file can either be opened for reading or writing. Not for both.\n" << *this;
        }
    }

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
                mWritingData.push_back(std::make_pair(Column(rKey, rValue, mpFormatSettings), rValue));
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
            p_itr->second = rValue;
        }
    } else {
        // now information with new Step is given.
        // so we first need to check if the headers are written.
        if (mLastWrittenStep == -1) {
            std::ofstream output_file(mFileName);
            // no step data has written, so no headers are also written.
            // writing the headers
            WriteHeaders(output_file);
            // now write the data
            WriteData(output_file);
            output_file.close();
        } else {
            std::ofstream output_file(mFileName, std::ios::app);
            // now write the data
            WriteData(output_file);
            output_file.close();
        }

        mLastWrittenStep = Step;
        mCurrentStep = Step;
        mWritingData[0].second = Step;
        Write(rValue, Step, rKey);
    }

    KRATOS_CATCH("");
}

// template instantiations
template void TransientCSVDatabaseIO::GenericRead<bool>(bool&, int, const std::string&);
template void TransientCSVDatabaseIO::GenericRead<int>(int&, int, const std::string&);
template void TransientCSVDatabaseIO::GenericRead<double>(double&, int, const std::string&);
template void TransientCSVDatabaseIO::GenericRead<std::string>(std::string&, int, const std::string&);

template void TransientCSVDatabaseIO::GenericWrite<bool>(const bool&, int, const std::string&);
template void TransientCSVDatabaseIO::GenericWrite<int>(const int&, int, const std::string&);
template void TransientCSVDatabaseIO::GenericWrite<double>(const double&, int, const std::string&);
template void TransientCSVDatabaseIO::GenericWrite<std::string>(const std::string&, int, const std::string&);

} // namespace Kratos