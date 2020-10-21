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

#if !defined(KRATOS_FILESYSTEM)
#define KRATOS_FILESYSTEM

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos {
// wrapper functions for std::filesystem (part of C++17)
// the function signatures are identical, hence after moving to C++17 Kratos::filesystem can be replaced with std::filesystem
// please check the documentation of std::filesystem for the function documentation

// Note: the filesystem functinos have a filesystem::path as input, but currently std::string is used as filesystem::path is not available
// this should not be a problem for upgrading to std::filesystem, since filesystem::path has a constructor accepting a string
namespace filesystem {

bool KRATOS_API(KRATOS_CORE) exists(const std::string& rPath);

bool KRATOS_API(KRATOS_CORE) is_regular_file(const std::string& rPath);

bool KRATOS_API(KRATOS_CORE) is_directory(const std::string& rPath);

bool KRATOS_API(KRATOS_CORE) create_directory(const std::string& rPath);

bool KRATOS_API(KRATOS_CORE) create_directories(const std::string& rPath);

bool KRATOS_API(KRATOS_CORE) remove(const std::string& rPath);

std::uintmax_t KRATOS_API(KRATOS_CORE) remove_all(const std::string& rPath);

void KRATOS_API(KRATOS_CORE) rename(const std::string& rPathFrom, const std::string& rPathTo);

// holder for file name data retrieved from the pattern and filename
struct FileNameData {
    FileNameData()
    {
    }

    FileNameData(const int Rank, const int Step, const double TimeStep)
    {
        this->Rank = Rank;
        this->Step = Step;
        this->TimeStep = TimeStep;
    }

    std::string FileName{""};
    int Rank{-1};
    int Step{-1};
    double TimeStep{-1.0};
};

class KRATOS_API(KRATOS_CORE) PatternSection
{
public:

    ///@name Life Cycle
    ///@{

    // Constructor
    PatternSection(
        const std::string& rPatternSection,
        const std::string& rPatternValueFormat = "");

    ///@}
    ///@ name Public operations

    bool UpdateFileNameData(
        FileNameData& rFileNameData,
        std::size_t& rCurrentPosition,
        const std::string& rData) const
    {
        return (this->*mUpdateData)(rFileNameData, rCurrentPosition, rData);
    }

    std::string GetValueString(
        const ModelPart& rModelPart) const
    {
        return (this->*mGetValueString)(rModelPart);
    }

    bool IsFlag() const
    {
        return (mPatternSectionString.front() == '<' && mPatternSectionString.back() == '>');
    }

    std::string GetPatternSectionString() const
    {
        return mPatternSectionString;
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    const std::string mPatternSectionString;
    const std::string mPatternValueFormat;
    bool (PatternSection::*mUpdateData)(FileNameData& rFileNameData,  std::size_t& rCurrentPosition, const std::string& rData) const;
    std::string (PatternSection::*mGetValueString)(const ModelPart& rModelPart) const;

    ///@}
    ///@name Private operations
    ///@{

    bool UpdateRank(FileNameData& rFileNameData, std::size_t& rCurrentPosition, const std::string& rData) const;
    bool UpdateStep(FileNameData& rFileNameData, std::size_t& rCurrentPosition, const std::string& rData) const;
    bool UpdateTimeStep(FileNameData& rFileNameData, std::size_t& rCurrentPosition, const std::string& rData) const;
    bool UpdateString(FileNameData& rFileNameData, std::size_t& rCurrentPosition, const std::string& rData) const;

    std::string GetRankString(const ModelPart& rModelPart) const;
    std::string GetStepString(const ModelPart& rModelPart) const;
    std::string GetTimeStepString(const ModelPart& rModelPart) const;
    std::string GetString(const ModelPart& rModelPart) const;

    static bool RetrieveIntegerValue(
        int& rValue,
        std::size_t& rCurrentPosition,
        const std::string& rData);

    static bool RetrieveFloatingPointValue(
        double& rValue,
        std::size_t& rCurrentPosition,
        const std::string& rData);

    ///@}
};

class KRATOS_API(KRATOS_CORE) FileNameInformationCollector
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Mesh
    KRATOS_CLASS_POINTER_DEFINITION(FileNameInformationCollector);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new File Name Information Collector object
     *
     * Construct an object which will handle given filename pattern. This rFileNamePattern
     * can incldue folders, subfolders. It must include the file name. Followings are
     * the list of flags supported.
     *
     * Accepted flags are:
     *      "<rank>"                : ranks of the model part
     *      "<step>"                : STEP value in model parts's process info
     *      "<time>"                : TIME value in model part's process info
     *      "<model_part_name>"     : name of the model part
     *      "<model_part_full_name>": full name of the model part.
     *
     * Flags "<rank>", "<step>" and "<time>" is only allowed be in the file name only. (Those flags are not
     * allowed to be present in file path)
     *
     * Example rFileNamePattern:
     *      1. "test_cases/<model_part_name>/<model_part_full_name>-<time>.h5"
     *      2. "test/test_1_<rank>_<time>.h5"
     *
     * Example rTimeStepFormat:
     *      1. "0.4f"
     *      2. "%0.4f"
     *      3. "0.3e"
     *
     * @param rModelPart
     * @param rFileNamePattern
     * @param rTimeStepFormat
     */

    FileNameInformationCollector(
        const ModelPart& rModelPart,
        const std::string& rFileNamePattern,
        const std::string& rTimeStepFormat = "");

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Get the File Name string
     *
     * This method returns the file name string according to file name pattern
     * and time step format provided in the constructor.
     *
     * @return std::string
     */
    std::string GetFileName() const;

    /**
     * @brief Gets sorted list of files
     *
     * This method returns list of file names matching the pattern given in the constructor
     * The list of file names are sorted in the flags order given in rFlagsSortingOrder.
     *
     * Accepted flags are:
     *      "<rank>" : ranks of the model part
     *      "<step>" : STEP value in model parts's process info
     *      "<time>" : TIME value in model part's process info
     *
     * @param rFlagsSortingOrder            Sorting order list of strings
     * @return std::vector<std::string>     List of sorted file names
     */
    std::vector<std::string> GetSortedFileNamesList(
        const std::vector<std::string>& rFlagsSortingOrder) const;

    /**
     * @brief Get filename data.
     *
     * This will return FileNameData object containing rank, step, time
     * if they are available in the given file name pattern. rFileNameWithoutPath
     * should be only filename without the path.
     *
     * If the corrensponding rFileNameWithoutPath does not math the pattern given in the
     * constructor then this returns false; otherwise true.
     *
     * @param rFileNameData             Data retrieved from rFileNameWithoutPath
     * @param rFileNameWithoutPath      File name without file path
     * @return true                     rFileNameWithoutPath matches the given pattern
     * @return false                    rFileNameWithoutPath does not match the given pattern
     */
    bool RetrieveFileNameInformation(
        FileNameData& rFileNameData,
        const std::string& rFileNameWithoutPath) const;

    /**
     * @brief Retrieves list of file name data objects
     *
     * This iterates through the constructor given file name pattern path files, and
     * returns file name data vector containing file name data for matching file names.
     *
     * @return std::vector<FileNameData>
     */
    std::vector<FileNameData> GetFileNameDataList() const;

    /**
     * @brief Sorts given list of file name data vector
     *
     * This method sorts given list of file name data vector
     * according to sorting flags order given in rFlagsSortingOrder.
     *
     * Accepted flags are:
     *      "<rank>" : ranks of the model part
     *      "<step>" : STEP value in model parts's process info
     *      "<time>" : TIME value in model part's process info
     *
     * @param rFileNameDataList
     * @param rFlagsSortingOrder
     */
    static void SortListOfFileNameData(
        std::vector<FileNameData>& rFileNameDataList,
        const std::vector<std::string>& rFlagsSortingOrder);

    ///@}

private:

    ///@name Private member variables
    ///@{

    const ModelPart& mrModelPart;
    std::string mPatternPath;
    std::vector<PatternSection> mPatternFileNameSections;

    ///@}
    ///@name Private member operations
    ///@{

    static std::vector<std::string> GetPatternSections(const std::string& rPattern);

    static void FindAndReplace(
        std::string& rInputString,
        const std::string& rSearchString,
        const std::string& rReplaceString);

    ///@}
};

} // namespace filesystem


namespace FilesystemExtensions {
// helper functions related to filesystem

std::string KRATOS_API(KRATOS_CORE) CurrentWorkingDirectory();

std::string KRATOS_API(KRATOS_CORE) JoinPaths(const std::vector<std::string>& rPaths);

} // namespace FilesystemExtensions
} // namespace Kratos

#endif // KRATOS_FILESYSTEM defined