//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes
#include <cctype>

// External includes

// Project includes
#include "includes/kratos_filesystem.h"

// Include base h
#include "file_name_data_collector.h"

namespace Kratos
{

FileNameDataCollector::FileNameDataCollector(
    const ModelPart& rModelPart,
    const std::string& rFileNamePattern,
    const std::unordered_map<std::string, std::string>& rFlagFormatMap)
    : mrModelPart(rModelPart)
{
    KRATOS_TRY

    const auto& overall_patterns = GetPatternFlagStrings(rFileNamePattern);

    for (const auto& r_pattern : overall_patterns) {
        if (r_pattern.front() == '<' && r_pattern.back() == '>') {
            KRATOS_ERROR_IF_NOT(r_pattern == "<model_part_name>" ||
                                r_pattern == "<model_part_full_name>" ||
                                r_pattern == "<rank>" ||
                                r_pattern == "<time>" ||
                                r_pattern == "<step>")
                << "Unsupported flag name found. [ r_pattern = " << r_pattern
                << " ].\n Supported flags are:"
                << "\n\t\"<model_part_name>\""
                << "\n\t\"<model_part_full_name>\""
                << "\n\t\"<rank>\""
                << "\n\t\"<time>\""
                << "\n\t\"<step>\"\n";
        }
    }

    auto file_name_pattern = rFileNamePattern;
    FindAndReplace(file_name_pattern, "<model_part_name>", rModelPart.Name());
    FindAndReplace(file_name_pattern, "<model_part_full_name>", rModelPart.FullName());
    mFilePath = std::filesystem::path(file_name_pattern).parent_path().string();

    KRATOS_ERROR_IF(mFilePath.find("<rank>") != std::string::npos)
        << "Flag \"<rank>\" is not allowed to be used inside the file path. "
           "Please use it in the file name only. [ file_path = "
        << mFilePath << " ].\n";

    KRATOS_ERROR_IF(mFilePath.find("<step>") != std::string::npos)
        << "Flag \"<step>\" is not allowed to be used inside the file path. "
           "Please use it in the file name only. [ file_path = "
        << mFilePath << " ].\n";

    KRATOS_ERROR_IF(mFilePath.find("<time>") != std::string::npos)
        << "Flag \"<step>\" is not allowed to be used inside the file path. "
           "Please use it in the file name only. [ file_path = "
        << mFilePath << " ].\n";

    const auto& process_flag_formats = [&](const std::string& rFlagName) -> std::string {
        const auto p_flag_format = rFlagFormatMap.find(rFlagName);
        if (p_flag_format != rFlagFormatMap.end()) {
            const auto& r_flag_format = p_flag_format->second;
            if (r_flag_format != "" && r_flag_format.front() != '%') {
                return "%" + r_flag_format;
            } else {
                return r_flag_format;
            }
        } else {
            return "";
        }
    };

    const auto& file_name_flags = GetPatternFlagStrings(std::filesystem::path(file_name_pattern).filename().string());
    for (const auto& r_flag : file_name_flags) {
        mFileNamePatternFlags.emplace_back(PatternFlag(r_flag, process_flag_formats(r_flag)));
    }

    // check patterns sections
    for (std::size_t i = 1; i < mFileNamePatternFlags.size(); ++i) {
        const auto& r_pattern_flag = mFileNamePatternFlags[i];
        const auto& r_previous_pattern_flag = mFileNamePatternFlags[i - 1];

        KRATOS_ERROR_IF(r_pattern_flag.IsFlag() && r_previous_pattern_flag.IsFlag())
            << "Having two flags adjacent to each other is not allowed. Please "
               "separate \""
            << r_previous_pattern_flag.GetPatternFlagString() << "\" and \""
            << r_pattern_flag.GetPatternFlagString() << "\" flags by a seperation character/string which is not a starting with a digit. [ PatternFileName = "
            << rFileNamePattern << " ]\n";

        // check patterns starting with digits
        if (!r_pattern_flag.IsFlag()) {
            const auto& pattern = r_pattern_flag.GetPatternFlagString();
            KRATOS_ERROR_IF(std::isdigit(pattern.front()))
                << "Found a string pattern starting with a digit in between "
                   "two flags. This is not allowed. [ found string pattern = "
                << pattern << ", file name pattern : " << rFileNamePattern << " ].\n";
        }
    }

    KRATOS_CATCH("");
}

std::string FileNameDataCollector::GetFileName() const
{
    KRATOS_TRY

    std::string file_name = "";

    for (const auto& r_flag : mFileNamePatternFlags) {
        file_name += r_flag.GetValueString(mrModelPart);
    }

    return (std::filesystem::path(mFilePath) / file_name).string();

    KRATOS_CATCH("");
}

std::string FileNameDataCollector::GetPath() const
{
    return mFilePath;
}

bool FileNameDataCollector::RetrieveFileNameData(
    FileNameData& rFileNameData,
    const std::string& rFileNameWithoutPath) const
{
    KRATOS_TRY

    rFileNameData.Clear();

    std::size_t current_position = 0;
    for (const auto& r_flag : mFileNamePatternFlags) {
        // checks whether rFileNameWithoutPath matterns the constructor prescribed  file name pattern,
        // if so fille it with relevant data.
        if (!r_flag.UpdateFileNameData(rFileNameData, current_position, rFileNameWithoutPath)) {
            // filename data update unsuccessfull because given rFileNameWithoutPath doesn not match
            // constructor prescribed file name pattern.
            return false;
        }
    }

    rFileNameData.SetFileName((std::filesystem::path(mFilePath) / rFileNameWithoutPath).string());

    return true;

    KRATOS_CATCH("");
}

std::vector<FileNameDataCollector::FileNameData> FileNameDataCollector::GetFileNameDataList() const
{
    KRATOS_TRY

    std::vector<FileNameData> result;

    for (const auto& current_file : FilesystemExtensions::ListDirectory(mFilePath)) {
        if (std::filesystem::is_regular_file(current_file)) {
            const auto& current_file_name_without_path = current_file.filename().string();
            FileNameData file_name_data;
            if (RetrieveFileNameData(file_name_data, current_file_name_without_path)) {
                result.push_back(file_name_data);
            }
        }
    }

    return result;

    KRATOS_CATCH("");
}

std::vector<std::string> FileNameDataCollector::GetSortedFileNamesList(
    const std::vector<std::string>& rFlagsSortingOrder) const
{
    KRATOS_TRY

    auto file_name_data_list = GetFileNameDataList();
    SortListOfFileNameData(file_name_data_list, rFlagsSortingOrder);

    std::vector<std::string> result;
    result.reserve(file_name_data_list.size());
    for (const auto& r_item : file_name_data_list) {
        result.push_back(r_item.GetFileName());
    }

    return result;

    KRATOS_CATCH("");
}

void FileNameDataCollector::SortListOfFileNameData(
    std::vector<FileNameData>& rFileNameDataList,
    const std::vector<std::string>& rFlagsSortingOrder)
{
    KRATOS_TRY

    std::vector<std::function<double(const FileNameData&)>> comparator_value_getter_list;

    for (const auto& r_flag_name : rFlagsSortingOrder) {
        if (r_flag_name == "<rank>") {
            comparator_value_getter_list.push_back([](const FileNameData& rFileNameData) -> double{
                return rFileNameData.GetRank();
            });
        } else if (r_flag_name == "<time>") {
            comparator_value_getter_list.push_back([](const FileNameData& rFileNameData) -> double {
                return rFileNameData.GetTime();
            });
        } else if (r_flag_name == "<step>") {
            comparator_value_getter_list.push_back([](const FileNameData& rFileNameData) -> double {
                return rFileNameData.GetStep();
            });
        } else {
            KRATOS_ERROR << "Unsupported flag name found. [ r_flag_name = " << r_flag_name
                         << " ].\n Supported flags are:"
                         << "\n\t\"<rank>\""
                         << "\n\t\"<time>\""
                         << "\n\t\"<step>\"\n";
        }
    }

    std::sort(rFileNameDataList.begin(), rFileNameDataList.end(),
              [&](const FileNameData& rA, const FileNameData& rB) {
                  for (const auto& r_comparator_value_getter : comparator_value_getter_list) {
                      const auto value_a = r_comparator_value_getter(rA);
                      const auto value_b = r_comparator_value_getter(rB);
                      if (value_a != value_b) {
                          return value_a < value_b;
                      }
                  }
                  return true;
              });

    KRATOS_CATCH("");
}

std::pair<bool, std::string> FileNameDataCollector::ExtractFileNamePattern(
    const std::string& rFileName,
    const std::vector<std::string>& rFlagsList)
{
    KRATOS_TRY

    int temp_int;
    double temp_double;

    std::size_t pos = 0;
    std::size_t current_flag_index = 0;
    std::string file_name_pattern = "";
    while (pos < rFileName.length()) {
        auto current_pos = pos;
        if (current_flag_index < rFlagsList.size()) {
            const auto& current_flag = rFlagsList[current_flag_index];
            if (current_flag == "<time>") {
                if (PatternFlag::RetrieveFloatingPointValue(temp_double, current_pos, rFileName)) {
                    file_name_pattern += current_flag;
                    current_flag_index++;
                    pos = current_pos;
                } else {
                    file_name_pattern += rFileName[pos];
                    pos++;
                }
            } else if (current_flag == "<step>" || current_flag == "<rank>") {
                if (PatternFlag::RetrieveIntegerValue(temp_int, current_pos, rFileName)) {
                    file_name_pattern += current_flag;
                    current_flag_index++;
                    pos = current_pos;
                } else {
                    file_name_pattern += rFileName[pos];
                    pos++;
                }
            } else if (current_flag == "<skip_float>") {
                if (PatternFlag::RetrieveFloatingPointValue(temp_double, current_pos, rFileName)) {
                    file_name_pattern += rFileName.substr(pos, current_pos - pos);
                    current_flag_index++;
                    pos = current_pos;
                } else {
                    file_name_pattern += rFileName[pos];
                    pos++;
                }
            } else if (current_flag == "<skip_int>") {
                if (PatternFlag::RetrieveIntegerValue(temp_int, current_pos, rFileName)) {
                    file_name_pattern += rFileName.substr(pos, current_pos - pos);
                    current_flag_index++;
                    pos = current_pos;
                } else {
                    file_name_pattern += rFileName[pos];
                    pos++;
                }
            } else {
                KRATOS_ERROR << "Unsupported flag name. [ flag_name = " << current_flag
                             << " ].\n"
                             << "Supported flags:"
                             << "\n\t<time>"
                             << "\n\t<rank>"
                             << "\n\t<step"
                             << "\n\t<skip_float>"
                             << "\n\t<skip_int";
            }
        } else {
            file_name_pattern += rFileName[pos++];
        }
    }

    return std::make_pair(current_flag_index == rFlagsList.size(), file_name_pattern);

    KRATOS_CATCH("");
}

void FileNameDataCollector::FindAndReplace(
    std::string& rInputString,
    const std::string& rSearchString,
    const std::string& rReplaceString)
{
    // Get the first occurrence
    size_t pos = rInputString.find(rSearchString);
    // Repeat till end is reached
    while (pos != std::string::npos) {
        // Replace this occurrence of Sub String
        rInputString.replace(pos, rSearchString.size(), rReplaceString);
        // Get the next occurrence from the current position
        pos = rInputString.find(rSearchString, pos + rReplaceString.size());
    }
}

std::vector<std::string> FileNameDataCollector::GetPatternFlagStrings(
    const std::string& rPattern)
{
    std::vector<std::string> result;

    const auto& add_section = [&](const std::string& rSection){
        if (rSection != "") {
            result.push_back(rSection);
        }
    };

    std::string current_section = "";
    for(const auto& c : rPattern) {
        if (c == '<') {
            add_section(current_section);
            current_section = "";
        }

        current_section += c;

        if (c == '>') {
            add_section(current_section);
            current_section = "";
        }
    }

    add_section(current_section);

    return result;
}

FileNameDataCollector::PatternFlag::PatternFlag(
    const std::string& rPatternFlag,
    const std::string& rPatternValueFormat)
    : mPatternFlag(rPatternFlag),
      mPatternValueFormat(rPatternValueFormat)
{
    if (rPatternFlag == "<time>") {
        this->mUpdateData = &PatternFlag::UpdateTime;
        this->mGetValueString = &PatternFlag::GetTimeString;
    } else if (rPatternFlag == "<step>") {
        this->mUpdateData = &PatternFlag::UpdateStep;
        this->mGetValueString = &PatternFlag::GetStepString;
    } else if (rPatternFlag == "<rank>") {
        this->mUpdateData = &PatternFlag::UpdateRank;
        this->mGetValueString = &PatternFlag::GetRankString;
    } else {
        this->mUpdateData = &PatternFlag::UpdateCustomString;
        this->mGetValueString = &PatternFlag::GetCustomString;
    }
}

std::string FileNameDataCollector::PatternFlag::GetRankString(
    const ModelPart& rModelPart) const
{
    return GetFormattedValueString<int>(rModelPart.GetCommunicator().MyPID());
}

std::string FileNameDataCollector::PatternFlag::GetStepString(
    const ModelPart& rModelPart) const
{
    return GetFormattedValueString<int>(rModelPart.GetProcessInfo()[STEP]);
}

std::string FileNameDataCollector::PatternFlag::GetTimeString(
    const ModelPart& rModelPart) const
{
    return GetFormattedValueString<double>(rModelPart.GetProcessInfo()[TIME]);
}

std::string FileNameDataCollector::PatternFlag::GetCustomString(
    const ModelPart& rModelPart) const
{
    return GetPatternFlagString();
}

bool FileNameDataCollector::PatternFlag::UpdateRank(
    FileNameData& rFileNameData,
    std::size_t& rCurrentPosition,
    const std::string& rData) const
{
    return RetrieveIntegerValue(rFileNameData.mRank, rCurrentPosition, rData);
}

bool FileNameDataCollector::PatternFlag::UpdateStep(
    FileNameData& rFileNameData,
    std::size_t& rCurrentPosition,
    const std::string& rData) const
{
    return RetrieveIntegerValue(rFileNameData.mStep, rCurrentPosition, rData);
}

bool FileNameDataCollector::PatternFlag::UpdateTime(
    FileNameData& rFileNameData,
    std::size_t& rCurrentPosition,
    const std::string& rData) const
{
    return RetrieveFloatingPointValue(rFileNameData.mTime, rCurrentPosition, rData);
}

bool FileNameDataCollector::PatternFlag::UpdateCustomString(
    FileNameData& rFileNameData,
    std::size_t& rCurrentPosition,
    const std::string& rData) const
{
    if (rData.substr(rCurrentPosition, mPatternFlag.size()) == mPatternFlag) {
        rCurrentPosition += mPatternFlag.size();
        return true;
    } else {
        return false;
    }
}

bool FileNameDataCollector::PatternFlag::RetrieveIntegerValue(
    int& rValue,
    std::size_t& rCurrentPosition,
    const std::string& rData)
{
    bool found_digit = false;
    std::string s_value = "";
    for (; rCurrentPosition < rData.length(); ++rCurrentPosition) {
        const auto& c = rData[rCurrentPosition];
        if (std::isdigit(c)) {
            s_value += c;
            found_digit = true;
        } else if (found_digit || c != ' ') {
            break;
        }
    }

    if (s_value != "") {
        rValue = std::stoi(s_value);
        return true;
    } else {
        return false;
    }
}

bool FileNameDataCollector::PatternFlag::RetrieveFloatingPointValue(
    double& rValue,
    std::size_t& rCurrentPosition,
    const std::string& rData)
{
    bool found_digit = false;
    bool found_point = false;
    bool found_e = false;

    std::string s_value = "";

    for (; rCurrentPosition < rData.length(); ++rCurrentPosition) {
        const auto c = rData[rCurrentPosition];
        if (std::isdigit(c)) {
            found_digit = true;
            s_value += c;
        } else if (c == '.' && !found_point && found_digit && !found_e) {
            found_point = true;
            s_value += c;
        } else if ((c == 'e' || c == 'E') && !found_e && found_digit && (rCurrentPosition + 2 < rData.length())) {
            const auto n_c = rData[rCurrentPosition + 1];
            const auto nn_c = rData[rCurrentPosition + 2];
            if ((n_c == '-' || n_c == '+') && (std::isdigit(nn_c))) {
                found_e = true;
                s_value += c;
                s_value += n_c;
                s_value += nn_c;
                rCurrentPosition += 2;
            } else {
                break;
            }
        } else {
            break;
        }
    }

    if (s_value != "") {
        rValue = std::stod(s_value);
        return true;
    } else {
        return false;
    }
}

} // namespace Kratos
