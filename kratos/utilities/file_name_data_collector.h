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

#if !defined(KRATOS_FILE_NAME_DATA_COLLECTOR)
#define KRATOS_FILE_NAME_DATA_COLLECTOR

// System includes
#include <string>
#include <unordered_map>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) FileNameDataCollector
{
    // private forward declarations
    class PatternFlag;

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FileNameDataCollector
    KRATOS_CLASS_POINTER_DEFINITION(FileNameDataCollector);

    ///@}
    ///@name Public classes
    ///@{

    class FileNameData
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of FileNameData
        KRATOS_CLASS_POINTER_DEFINITION(FileNameData);

        ///@}
        ///@name Life Cycle
        ///@{

        // Default constructor
        FileNameData()
        {}

        FileNameData(
            const std::string& rFileName,
            const int Rank,
            const int Step,
            const double Time)
            : mFileName(rFileName),
              mRank(Rank),
              mStep(Step),
              mTime(Time)
        {}

        ///@}
        ///@name Public member access
        ///@{

        void SetFileName(const std::string& rFileName) { mFileName = rFileName; }

        std::string GetFileName() const { return mFileName; }

        void SetRank(const int Rank) { mRank = Rank; }

        int GetRank() const { return mRank; }

        void SetStep(const int Step) { mStep = Step; }

        int GetStep() const { return mStep; }

        void SetTime(const double Time) { mTime = Time; }

        double GetTime() const { return mTime; }

        ///@}
        ///@name Public member operations
        ///@{

        void Clear()
        {
            mFileName = "";
            mRank = -1;
            mStep = -1;
            mTime = -1.0;
        }

        ///}
        ///@name Public operators
        ///@{

        bool operator==(const FileNameData& rRHS)
        {
            return (mFileName == rRHS.mFileName && mRank == rRHS.mRank &&
                    mStep == rRHS.mStep && mTime == rRHS.mTime);
        }

        ///@}

    private:
        ///@name Private members variables
        ///@{

        std::string mFileName = "";
        int mRank = -1;
        int mStep = -1;
        double mTime = -1.0;

        friend class PatternFlag;

        ///@}
    };

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
     * rFileNamePattern should always have a seperation character/string in between two flags. This seperation
     * character/string can not start with a digit.
     *
     * rFlagFormatMap holds format for each flag. If no format is specified, then
     * value for that flag will be written without any formatting.
     * Example map:
     *          {
     *              "<time>": "0.4f",
     *              "<step>": "%4d"
     *          }
     *
     * Example formats:
     *      1. "0.4f"
     *      2. "%0.4f"
     *      3. "0.3e"
     *
     * @param rModelPart            Model part on which this file name data collection is based on
     * @param rFileNamePattern      File name pattern (can have file name with or without file path)
     * @param rFlagFormatMap        Map of formats (key : flag, value: format [both as strings])
     */

    FileNameDataCollector(
        const ModelPart& rModelPart,
        const std::string& rFileNamePattern,
        const std::unordered_map<std::string, std::string>& rFlagFormatMap);

    ///@}
    ///@name Public member operations
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
     * @brief Get the Path
     *
     * This method returns path obtained from the constructor specified
     * pattern
     *
     * @return std::string              Path from the pattern
     */
    std::string GetPath() const;

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
     * If retrieval is successfully, FileNameData::mFileName will be filled with
     * full file name including the path (specified by the constructor pattern). It will
     * not check for existence of the file in the same path.
     *
     * @param rFileNameData             Data retrieved from rFileNameWithoutPath
     * @param rFileNameWithoutPath      File name without file path
     * @return true                     rFileNameWithoutPath matches the given pattern
     * @return false                    rFileNameWithoutPath does not match the given pattern
     */
    bool RetrieveFileNameData(
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

    ///@}
    ///@name Public static operations
    ///@{

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

    /**
     * @brief Extracts file name pattern from a given file name
     *
     * This method returns the file name pattern matching given order of flags from left to right.
     * If all the given flags are found, this returns true, otherwise false.
     *
     * Example patterns:
     *      Check test_file_name_data_collector.py::test_ExtractFileNamePattern
     *
     * @param rFileName                         Filename with or without file path
     * @param rFlagsList                        List of flags exists in the file name
     * @return std::pair<bool, std::string>     bool: true if all flags are found, false otherwise; string: found pattern
     */
    static std::pair<bool, std::string> ExtractFileNamePattern(
        const std::string& rFileName,
        const std::vector<std::string>& rFlagsList);

    ///@}

private:
    ///@name Private classes
    ///@{

    /**
     * @brief This class holds pattern flags
     *
     * This class is used to hold different sections of the pattern and the flags.
     *
     */
    class KRATOS_API(KRATOS_CORE) PatternFlag
    {
    public:
        ///@name Life Cycle
        ///@{

        // Constructor
        /**
         * @brief Construct a new Pattern Flag object
         *
         * This construct a PatterFlag by checking the passed rPatternFlag.
         * So string checking for rPatternFlag is done only here.
         *
         * @param rPatternFlag          String pattern
         * @param rPatternValueFormat   Pattern value format used in GetValueString
         */
        PatternFlag(
            const std::string& rPatternFlag,
            const std::string& rPatternValueFormat = "");

        ///@}
        ///@ name Public member operations
        ///@{

        /**
         * @brief Updates FileNameData container from rData string
         *
         * This method updates file name data from rData for specified
         * rPatternFlag given in the constructor.
         *
         * @param rFileNameData         FileName data container which is filled from rData
         * @param rCurrentPosition      Current position from where this pattern flag should try to match
         * @param rData                 File name without the path
         * @return true                 If rPatternFlag matches given rData starting at rCurrentPosition
         * @return false                If rPatternFlag does not match given rData starting at rCurrentPosition
         */
        bool UpdateFileNameData(
            FileNameData& rFileNameData,
            std::size_t& rCurrentPosition,
            const std::string& rData) const
        {
            return (this->*mUpdateData)(rFileNameData, rCurrentPosition, rData);
        }

        /**
         * @brief Constructs file name from the pattern
         *
         * @param rModelPart            Model part to get information to construct file name
         * @return std::string          Constructed file name (full path will be returned)
         */
        std::string GetValueString(
            const ModelPart& rModelPart) const
        {
            return (this->*mGetValueString)(rModelPart);
        }


        /**
         * @brief Checks whether the PatternFlag is a Flag or not
         *
         * @return true                 If PatternFlag is a flag
         * @return false                If Pattern flag is not a flag
         */
        bool IsFlag() const
        {
            return (mPatternFlag.front() == '<' && mPatternFlag.back() == '>');
        }

        /**
         * @brief Get the Pattern Flag String
         *
         * @return std::string
         */
        std::string GetPatternFlagString() const
        {
            return mPatternFlag;
        }

        ///@}
        ///@name Public static operations
        ///{

        /**
         * @brief Retrieves integer value from rData
         *
         * This method tries to retrieve an integer value from
         * rData starting from rCurrenPosition of rData string.
         *
         * Following string are recognized as integers
         *      1. "123"
         *      2. "   123"
         *
         * @param rValue                        Output integer value
         * @param rCurrentPosition              Starting position of rData string to look for integer
         * @param rData                         String to look for integer
         * @return true                         If integer retrieval is successfull
         * @return false                        If integer retrieval is not successfull
         */
        static bool RetrieveIntegerValue(
            int& rValue,
            std::size_t& rCurrentPosition,
            const std::string& rData);

        /**
         * @brief Retrieves floating point value from rData
         *
         * This method tries to retrieve a floating point value from
         * rData starting from rCurrenPosition of rData string.
         *
         * Followings patterns currently accepted as floating point values:
         *      1. "10"
         *      2. "12.2"
         *      3. "10e-5"
         *      4. "13e+2"
         *      5. "1.4e-3"
         *      6. "12.4E+10"
         *      7. "1.4E-3"
         *
         * @param rValue                        Output floating value
         * @param rCurrentPosition              Starting position of rData string to look for floating value
         * @param rData                         String to look for floating value
         * @return true                         If floating value retrieval is successfull
         * @return false                        If floating value retrieval is not successfull
         */
        static bool RetrieveFloatingPointValue(
            double& rValue,
            std::size_t& rCurrentPosition,
            const std::string& rData);

        ///@}

    private:
        ///@name Private member variables
        ///@{

        const std::string mPatternFlag;
        const std::string mPatternValueFormat;
        bool (PatternFlag::*mUpdateData)(FileNameData& rFileNameData, std::size_t& rCurrentPosition, const std::string& rData) const;
        std::string (PatternFlag::*mGetValueString)(const ModelPart& rModelPart) const;

        ///@}
        ///@name Private member operations
        ///@{

        bool UpdateRank(FileNameData& rFileNameData, std::size_t& rCurrentPosition, const std::string& rData) const;
        bool UpdateStep(FileNameData& rFileNameData, std::size_t& rCurrentPosition, const std::string& rData) const;
        bool UpdateTime(FileNameData& rFileNameData, std::size_t& rCurrentPosition, const std::string& rData) const;
        bool UpdateCustomString(FileNameData& rFileNameData, std::size_t& rCurrentPosition, const std::string& rData) const;

        std::string GetRankString(const ModelPart& rModelPart) const;
        std::string GetStepString(const ModelPart& rModelPart) const;
        std::string GetTimeString(const ModelPart& rModelPart) const;
        std::string GetCustomString(const ModelPart& rModelPart) const;

        /**
         * @brief Get formatted value string
         *
         * @tparam TDataType                    Data type
         * @param rData                         Value to be formatted
         * @return std::string                  Formatting string
         */
        template <class TDataType>
        std::string GetFormattedValueString(
            const TDataType& rData) const
        {
            if (mPatternValueFormat == "") {
                return std::to_string(rData);
            } else {
                int length = std::snprintf(nullptr, 0, mPatternValueFormat.c_str(), rData);
                assert(length >= 0);

                char* buf = new char[length + 1];
                std::snprintf(buf, length + 1, mPatternValueFormat.c_str(), rData);

                std::string data_str(buf);
                delete[] buf;

                return data_str;
            }
        }

        ///@}
    };

    ///@}
    ///@name Private member variables
    ///@{

    const ModelPart& mrModelPart;
    std::string mFilePath;
    std::vector<PatternFlag> mFileNamePatternFlags;

    ///@}
    ///@name Private static operations
    ///@{

    /**
     * @brief Get the pattern sections and flags
     *
     * @param rPattern                      File name pattern
     * @return std::vector<std::string>     Sections of file name pattern
     */
    static std::vector<std::string> GetPatternFlagStrings(
        const std::string& rPattern);

    /**
     * @brief Simple recursive find and replace method for strings
     *
     * @param rInputString                  String to find and replace
     * @param rSearchString                 Search string which is replaced
     * @param rReplaceString                Replacement string
     */
    static void FindAndReplace(
        std::string& rInputString,
        const std::string& rSearchString,
        const std::string& rReplaceString);

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_FILE_NAME_DATA_COLLECTOR
