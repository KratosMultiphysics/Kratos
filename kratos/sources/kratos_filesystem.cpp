//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes
#include <algorithm>

// External includes
#include "ghc/filesystem.hpp" // TODO after moving to C++17 this can be removed since the functions can be used directly

// Project includes
#include "includes/kratos_filesystem.h"

namespace Kratos {
namespace filesystem {

bool exists(const std::string& rPath)
{
    return ghc::filesystem::exists(rPath);
}


bool is_regular_file(const std::string& rPath)
{
    return ghc::filesystem::is_regular_file(rPath);
}


bool is_directory(const std::string& rPath)
{
    return ghc::filesystem::is_directory(rPath);
}


bool create_directory(const std::string& rPath)
{
    return ghc::filesystem::create_directory(rPath);
}


bool create_directories(const std::string& rPath)
{
    return ghc::filesystem::create_directories(rPath);
}


bool remove(const std::string& rPath)
{
    return ghc::filesystem::remove(rPath);
}


std::uintmax_t remove_all(const std::string& rPath)
{
    return ghc::filesystem::remove_all(rPath);
}


void rename(const std::string& rPathFrom, const std::string& rPathTo)
{
    return ghc::filesystem::rename(rPathFrom, rPathTo);
}

PatternSection::PatternSection(
    const std::string& rPatternSection,
    const std::string& rPatternValueFormat)
    : mPatternSectionString(rPatternSection),
      mPatternValueFormat(rPatternValueFormat)
{
    if (rPatternSection == "<time>") {
        this->mUpdateData = &PatternSection::UpdateTimeStep;
        this->mGetValueString = &PatternSection::GetTimeStepString;
    } else if (rPatternSection == "<step>") {
        this->mUpdateData = &PatternSection::UpdateStep;
        this->mGetValueString = &PatternSection::GetStepString;
    } else if (rPatternSection == "<rank>") {
        this->mUpdateData = &PatternSection::UpdateRank;
        this->mGetValueString = &PatternSection::GetRankString;
    } else {
        this->mUpdateData = &PatternSection::UpdateString;
        this->mGetValueString = &PatternSection::GetString;
    }
}

std::string PatternSection::GetRankString(
    const ModelPart& rModelPart) const
{
    return std::to_string(rModelPart.GetCommunicator().MyPID());
}

std::string PatternSection::GetStepString(
    const ModelPart& rModelPart) const
{
    return std::to_string(rModelPart.GetProcessInfo()[STEP]);
}

std::string PatternSection::GetTimeStepString(
    const ModelPart& rModelPart) const
{
    const auto& current_time = rModelPart.GetProcessInfo()[TIME];
    if (mPatternValueFormat == "") {
        return std::to_string(current_time);
    } else {
        int length = std::snprintf(nullptr, 0, mPatternValueFormat.c_str(), current_time);
        assert(length >= 0);

        char* buf = new char[length + 1];
        std::snprintf(buf, length + 1, mPatternValueFormat.c_str(), current_time);

        std::string current_time_str(buf);
        delete[] buf;

        return current_time_str;
    }
}

std::string PatternSection::GetString(
    const ModelPart& rModelPart) const
{
    return GetPatternSectionString();
}

bool PatternSection::UpdateRank(
    FileNameData& rFileNameData,
    std::size_t& rCurrentPosition,
    const std::string& rData) const
{
    return RetrieveIntegerValue(rFileNameData.Rank, rCurrentPosition, rData);
}

bool PatternSection::UpdateStep(
    FileNameData& rFileNameData,
    std::size_t& rCurrentPosition,
    const std::string& rData) const
{
    return RetrieveIntegerValue(rFileNameData.Step, rCurrentPosition, rData);
}

bool PatternSection::UpdateTimeStep(
    FileNameData& rFileNameData,
    std::size_t& rCurrentPosition,
    const std::string& rData) const
{
    return RetrieveFloatingPointValue(rFileNameData.TimeStep, rCurrentPosition, rData);
}

bool PatternSection::UpdateString(
    FileNameData& rFileNameData,
    std::size_t& rCurrentPosition,
    const std::string& rData) const
{
    if (rData.substr(rCurrentPosition, mPatternSectionString.size()) == mPatternSectionString) {
        rCurrentPosition += mPatternSectionString.size();
        return true;
    } else {
        return false;
    }
}

bool PatternSection::RetrieveIntegerValue(
    int& rValue,
    std::size_t& rCurrentPosition,
    const std::string& rData)
{
    std::string s_value = "";
    for (; rCurrentPosition < rData.size(); ++rCurrentPosition) {
        const auto& c = rData[rCurrentPosition];
        if (std::isdigit(c)) {
            s_value += c;
        } else {
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

bool PatternSection::RetrieveFloatingPointValue(
    double& rValue,
    std::size_t& rCurrentPosition,
    const std::string& rData)
{
    bool found_digit = false;
    bool found_point = false;
    bool found_e = false;

    std::string s_value = "";

    for (; rCurrentPosition < rData.size(); ++rCurrentPosition) {
        const auto c = rData[rCurrentPosition];
        if (isdigit(c)) {
            found_digit = true;
            s_value += c;
        } else if (c == '.' && !found_point && found_digit && !found_e) {
            found_point = true;
            s_value += c;
        } else if ((c == 'e' || c == 'E') && !found_e && found_digit && (rCurrentPosition + 2 < rData.size())) {
            const auto n_c = rData[rCurrentPosition + 1];
            const auto nn_c = rData[rCurrentPosition + 2];
            if ((n_c == '-' || n_c == '+') && (isdigit(nn_c))) {
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

FileNameInformationCollector::FileNameInformationCollector(
    const ModelPart& rModelPart,
    const std::string& rFileNamePattern,
    const std::string& rTimeStepFormat)
    : mrModelPart(rModelPart)
{
    KRATOS_TRY

    const auto& overall_pattern_sections = GetPatternSections(rFileNamePattern);

    for (const auto& r_pattern_section : overall_pattern_sections) {
        if (r_pattern_section.front() == '<' && r_pattern_section.back() == '>') {
            KRATOS_ERROR_IF_NOT(r_pattern_section == "<model_part_name>" ||
                                r_pattern_section == "<model_part_full_name>" ||
                                r_pattern_section == "<rank>" ||
                                r_pattern_section == "<time>" || r_pattern_section == "<step>")
                << "Unsupported flag name found. [ r_pattern_section = " << r_pattern_section
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

    const auto& path_file_name_pattern = ghc::filesystem::path(file_name_pattern);
    mPatternPath = path_file_name_pattern.parent_path();

    KRATOS_ERROR_IF(mPatternPath.find("<rank>") != std::string::npos)
        << "Flag \"<rank>\" is not allowed to be used inside the file path. "
           "Please use it in the file name only. [ file_path = "
        << mPatternPath << " ].\n";

    KRATOS_ERROR_IF(mPatternPath.find("<step>") != std::string::npos)
        << "Flag \"<step>\" is not allowed to be used inside the file path. "
           "Please use it in the file name only. [ file_path = "
        << mPatternPath << " ].\n";

    KRATOS_ERROR_IF(mPatternPath.find("<time>") != std::string::npos)
        << "Flag \"<step>\" is not allowed to be used inside the file path. "
           "Please use it in the file name only. [ file_path = "
        << mPatternPath << " ].\n";

    const auto& file_name_pattern_sections = GetPatternSections(path_file_name_pattern.filename());

    for (const auto& r_pattern_section : file_name_pattern_sections) {
        if (r_pattern_section == "<time>") {
            mPatternFileNameSections.emplace_back(PatternSection(
                r_pattern_section, (rTimeStepFormat != "")
                                       ? (rTimeStepFormat[0] == '%')
                                            ? rTimeStepFormat
                                            : "%" + rTimeStepFormat
                                       : ""));
        } else {
            mPatternFileNameSections.emplace_back(PatternSection(r_pattern_section));
        }
    }

    // check patterns sections
    for (std::size_t i = 1; i < mPatternFileNameSections.size(); ++i) {
        KRATOS_ERROR_IF(mPatternFileNameSections[i].IsFlag() &&
                        mPatternFileNameSections[i - 1].IsFlag())
            << "Having two flags adjacent to each other is not allowed. Please "
               "separate \""
            << mPatternFileNameSections[i - 1].GetPatternSectionString() << "\" and \""
            << mPatternFileNameSections[i].GetPatternSectionString() << "\" flags by a seperation character which is not a digit. [ PatternFileName = "
            << rFileNamePattern << " ]\n";
    }

    KRATOS_CATCH("");
}

std::string FileNameInformationCollector::GetFileName() const
{
    std::string file_name = "";
    for (const auto& r_current_pattern : mPatternFileNameSections) {
        file_name += r_current_pattern.GetValueString(mrModelPart);
    }

    return FilesystemExtensions::JoinPaths({mPatternPath, file_name});
}

bool FileNameInformationCollector::RetrieveFileNameInformation(
    FileNameData& rFileNameData,
    const std::string& rFileNameWithoutPath) const
{
    rFileNameData.Rank = -1;
    rFileNameData.Step = -1;
    rFileNameData.TimeStep = -1.0;

    std::size_t current_position = 0;
    for (auto& r_current_pattern : mPatternFileNameSections) {
        if (!r_current_pattern.UpdateFileNameData(
                rFileNameData, current_position, rFileNameWithoutPath)) {
            return false;
        }
    }

    return true;
}

std::vector<FileNameData> FileNameInformationCollector::GetFileNameDataList() const
{
    KRATOS_TRY

    std::vector<FileNameData> result;

    for (const auto& current_file : ghc::filesystem::directory_iterator(mPatternPath)) {
        const std::string& current_file_name_without_path =
            ghc::filesystem::path(current_file).filename();
        FileNameData file_name_data;
        if (RetrieveFileNameInformation(file_name_data, current_file_name_without_path)) {
            file_name_data.FileName = current_file.path();
            result.push_back(file_name_data);
        }
    }

    return result;

    KRATOS_CATCH("");
}

std::vector<std::string> FileNameInformationCollector::GetSortedFileNamesList(
    const std::vector<std::string>& rFlagsSortingOrder) const
{
    KRATOS_TRY

    std::vector<FileNameData> file_name_data_list = this->GetFileNameDataList();
    SortListOfFileNameData(file_name_data_list, rFlagsSortingOrder);

    std::vector<std::string> result;
    for (const auto& r_item : file_name_data_list) {
        result.push_back(r_item.FileName);
    }

    return result;

    KRATOS_CATCH("");
}

void FileNameInformationCollector::SortListOfFileNameData(
    std::vector<FileNameData>& rFileNameDataList,
    const std::vector<std::string>& rFlagsSortingOrder)
{
    KRATOS_TRY

    std::vector<std::function<double(const FileNameData&)>> comparator_value_getter_list;

    for (const auto& r_flag_name : rFlagsSortingOrder) {
        if (r_flag_name == "<rank>") {
            comparator_value_getter_list.push_back([](const FileNameData& rFileNameData) -> double{
                return rFileNameData.Rank;
            });
        } else if (r_flag_name == "<time>") {
            comparator_value_getter_list.push_back([](const FileNameData& rFileNameData) -> double {
                return rFileNameData.TimeStep;
            });
        } else if (r_flag_name == "<step>") {
            comparator_value_getter_list.push_back([](const FileNameData& rFileNameData) -> double {
                return rFileNameData.Step;
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

void FileNameInformationCollector::FindAndReplace(
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

std::vector<std::string> FileNameInformationCollector::GetPatternSections(const std::string& rPattern)
{
    std::vector<std::string> result;

    std::string current_section = "";
    for(const auto& c : rPattern) {
        if (c == '<') {
            result.push_back(current_section);
            current_section = "";
        }

        current_section += c;

        if (c == '>') {
            result.push_back(current_section);
            current_section = "";
        }
    }

    if (current_section != "") {
        result.push_back(current_section);
    }

    return result;
}

} // namespace filesystem


namespace FilesystemExtensions {

std::string CurrentWorkingDirectory()
{
    return ghc::filesystem::current_path().string();
}

std::string JoinPaths(const std::vector<std::string>& rPaths)
{
    auto paths(rPaths); // create local copy

    // first remove empty paths
    paths.erase(std::remove_if(paths.begin(), paths.end(),
                         [](const std::string& s)
                         { return s.empty(); }), paths.end());

    const std::size_t num_paths = paths.size();

    if (num_paths == 0) { return ""; }

    std::string full_path = paths[0];
    if (num_paths > 1) {
        for(std::size_t i=1; i<num_paths; ++i) {
            full_path += "/" + paths[i]; // using portable separator "/"
        }
    }

    return full_path;
}

} // namespace FilesystemExtensions
} // namespace Kratos
