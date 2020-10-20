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

FileNameInformationCollector::FileNameInformationCollector(
    const ModelPart& rModelPart,
    const std::string& rFileNamePattern,
    const std::string& rTimeStepFormat)
    : mrModelPart(rModelPart)
{
    KRATOS_TRY

    mTimeStepFormat = (rTimeStepFormat != "")
                          ? (rTimeStepFormat[0] == '%')
                                ? rTimeStepFormat
                                : "%" + rTimeStepFormat
                          : "";

    // check for flags in pattern
    std::string flag_name = "";
    bool found_flag = false;
    for (const auto& c : rFileNamePattern) {
        if (c == '<') {
            found_flag = true;
            flag_name = "<";
        } else if (c == '>') {
            found_flag = false;
            flag_name += ">";
            KRATOS_ERROR_IF_NOT(flag_name == "<model_part_name>" ||
                                flag_name == "<model_part_full_name>" ||
                                flag_name == "<rank>" ||
                                flag_name == "<time>" ||
                                flag_name == "<step>")
                << "Unsupported flag name found. [ flag_name = " << flag_name
                << " ].\n Supported flags are:"
                << "\n\t\"<model_part_name>\""
                << "\n\t\"<model_part_full_name>\""
                << "\n\t\"<rank>\""
                << "\n\t\"<time>\""
                << "\n\t\"<step>\"\n";
        } else {
            if (found_flag) {
                flag_name += c;
            }
        }
    }

    auto file_name_pattern = rFileNamePattern;
    FindAndReplace(file_name_pattern, "<model_part_name>", rModelPart.Name());
    FindAndReplace(file_name_pattern, "<model_part_full_name>", rModelPart.FullName());
    mPattern = file_name_pattern;

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

    mPatternFileName = path_file_name_pattern.filename();

    KRATOS_CATCH("");
}

std::string FileNameInformationCollector::GetFileName() const
{
    const auto& r_process_info = mrModelPart.GetProcessInfo();
    auto file_name = mPattern;
    FindAndReplace(file_name, "<step>", std::to_string(r_process_info[STEP]));
    FindAndReplace(file_name, "<rank>", std::to_string(mrModelPart.GetCommunicator().MyPID()));

    if (mTimeStepFormat == "") {
        FindAndReplace(file_name, "<time>", std::to_string(r_process_info[TIME]));
    } else {
        int length = std::snprintf(nullptr, 0, mTimeStepFormat.c_str(),  r_process_info[TIME]);
        assert( length >= 0 );

        char* buf = new char[length + 1];
        std::snprintf(buf, length + 1, mTimeStepFormat.c_str(), r_process_info[TIME]);

        std::string current_time( buf );
        delete[] buf;

        FindAndReplace(file_name, "<time>", current_time);
    }
    return file_name;
}

bool FileNameInformationCollector::RetrieveFileNameInformation(
    FileNameData& rFileNameData,
    const std::string& rFileNameWithoutPath) const
{
    rFileNameData.Rank = -1;
    rFileNameData.Step = -1;
    rFileNameData.TimeStep = -1.0;

    std::size_t file_name_pos = 0;
    std::size_t file_name_pattern_pos = 0;
    bool found_flag = false;
    std::string flag_name = "";
    std::string matching_string = "";
    bool is_file_name_matching = true;

    const auto& update_for_next_match = [&]() {
        // match string until now
        const auto& found_pos = rFileNameWithoutPath.find(matching_string, file_name_pos);
        if (found_pos == std::string::npos) {
            is_file_name_matching = false;
            return;
        } else {
            const auto& data = rFileNameWithoutPath.substr(
                file_name_pos, found_pos - file_name_pos);

            is_file_name_matching = is_file_name_matching &&
                ((flag_name == "" && data == "") || (flag_name != "" && data != ""));

            if (data != "") {
                switch (GetDataType(data)) {
                    case DataType::IntegerNumber:
                        break;
                    case DataType::FloatingPointNumber:
                        break;
                    case DataType::UnknownType:
                        is_file_name_matching = false;
                        return;
                };

                if (flag_name == "<step>") {
                    rFileNameData.Step = std::stoi(data);
                } else if (flag_name == "<rank>") {
                    rFileNameData.Rank = std::stoi(data);
                } else if (flag_name == "<time>") {
                    rFileNameData.TimeStep = std::stod(data);
                } else {
                    is_file_name_matching = false;
                    return;
                }
            }

            file_name_pos = found_pos + matching_string.size();
        }
    };

    while (file_name_pattern_pos < mPatternFileName.size()) {
        const auto& c_file_name_pattern = mPatternFileName[file_name_pattern_pos];
        if (c_file_name_pattern == '<') {
            update_for_next_match();
            found_flag = true;
            matching_string = "";
            flag_name = "<";
        } else if (c_file_name_pattern == '>') {
            found_flag = false;
            flag_name += ">";
        } else {
            if (!found_flag) {
                matching_string += c_file_name_pattern;
            } else {
                flag_name += c_file_name_pattern;
            }
        }

        if (!is_file_name_matching) {
            return false;
        }

        file_name_pattern_pos++;
    }

    // check the last mathing pattern
    update_for_next_match();
    return is_file_name_matching && (file_name_pos == rFileNameWithoutPath.size());
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

    for (const auto& r_flag_name : rFlagsSortingOrder) {
        KRATOS_ERROR_IF_NOT(r_flag_name == "<rank>" ||
                            r_flag_name == "<time>" ||
                            r_flag_name == "<step>")
            << "Unsupported flag name found. [ r_flag_name = " << r_flag_name
            << " ].\n Supported flags are:"
            << "\n\t\"<rank>\""
            << "\n\t\"<time>\""
            << "\n\t\"<step>\"\n";
    }

    std::sort(rFileNameDataList.begin(), rFileNameDataList.end(),
              [&](const FileNameData& rA, const FileNameData& rB) {
                  for (const auto& r_flag : rFlagsSortingOrder) {
                      if (r_flag == "<rank>") {
                          if (rA.Rank != rB.Rank) {
                              return rA.Rank < rB.Rank;
                          }
                      } else if (r_flag == "<time>") {
                          if (rA.TimeStep != rB.TimeStep) {
                              return rA.TimeStep < rB.TimeStep;
                          }
                      } else if (r_flag == "<step>") {
                          if (rA.Step != rB.Step) {
                              return rA.Step < rB.Step;
                          }
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

FileNameInformationCollector::DataType FileNameInformationCollector::GetDataType(
    const std::string& rData)
{
    std::size_t number_of_digits = 0;
    std::size_t number_of_points = 0;
    std::size_t number_of_scientific_notations = 0;
    std::size_t number_of_dashes = 0;
    std::size_t number_of_e = 0;

    for (std::size_t i = 0; i < rData.size(); ++i) {
        const auto& c = rData[i];
        if (std::isdigit(c)) {
            number_of_digits++;
        } else if (c == '.') {
            number_of_points++;
        } else if (c == '-') {
            number_of_dashes++;
            if (i > 0) {
                if (rData[i - 1] == 'e' || rData[i - 1] == 'E') {
                    number_of_scientific_notations++;
                }
            }
        } else if (c == 'e' || c == 'E') {
            number_of_e++;
        } else {
            return DataType::UnknownType;
        }
    }

    if (number_of_digits == rData.size()) {
        return DataType::IntegerNumber;
    } else if ((number_of_dashes == 1 && number_of_e == 1 && number_of_scientific_notations == 1) ||
               (number_of_dashes == 0 && number_of_e == 0 && number_of_points == 1)) {
        return DataType::FloatingPointNumber;
    } else {
        return DataType::UnknownType;
    }
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
