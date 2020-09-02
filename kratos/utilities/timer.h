//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//
//

#if !defined(KRATOS_TIMER_H_INCLUDED )
#define  KRATOS_TIMER_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <unordered_map>
#include <chrono>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class Timer
 * @ingroup KratosCore
 * @brief This utility can be used to compute the time employed on computations
 * @author Pooyan Dadvand
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) Timer
{
    /**
    * @class TimerData
    * @ingroup KratosCore
    * @brief This is an internal class used to manage the timer data
    * @author Pooyan Dadvand
    * @author Riccardo Rossi
    */
    class TimerData
    {
        int mRepeatNumber;
        double mStartTime;
        double mTotalElapsedTime;
        double mMaximumTime;
        double mMinimumTime;
    public:
        TimerData() : mRepeatNumber(int()), mStartTime(double()), mTotalElapsedTime(double()), mMaximumTime(double()), mMinimumTime(double()) {}
        double GetStartTime()
        {
            return mStartTime;
        }
        void SetStartTime(double StartTime)
        {
            mStartTime = StartTime;
        }
        double GetTotalElapsedTime()
        {
            return mTotalElapsedTime;
        }
        void SetTotalElapsedTime(double TotalElapsedTime)
        {
            mTotalElapsedTime = TotalElapsedTime;
        }
        void Update(double StopTime)
        {
            double elapsed = StopTime - mStartTime;
            if(mRepeatNumber == 0)
                mMinimumTime = elapsed;
            mTotalElapsedTime += elapsed;
            if(mMaximumTime < elapsed)
                mMaximumTime = elapsed;

            if((mMinimumTime > elapsed))
                mMinimumTime = elapsed;

            mRepeatNumber++;

        }
        /// Print object's data.
        void PrintData(std::ostream& rOStream, double GlobalElapsedTime = -1.00) const;
    };

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Timer
    KRATOS_CLASS_POINTER_DEFINITION(Timer);

    /// The type of float used to store the time
    typedef double TimeType;

    /// The timer data container type (map)
    typedef std::map<std::string, TimerData> ContainerType;

    /// This is used to know the internal name used for the time table
    typedef std::unordered_map<std::string, std::string> InternalNameDatabaseType;

    /// The number of 0 in the internal name
    static constexpr std::size_t NumberOfZeros = 3;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Timer();

    /// Destructor.
    virtual ~Timer()
    {

    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method starts the timer meassures
     * @param rIntervalName The internal name that will store the timing data
     */
    static void Start(std::string const& rIntervalName);

    /**
     * @brief This method stops the timer meassures
     * @param rIntervalName The internal name that will store the timing data
     */
    static void Stop(std::string const& rIntervalName);

    /**
     * @brief This method returns the resulting time
     */
    static inline double GetTime()
    {
        const auto current_time = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<std::chrono::duration<double>>(current_time.time_since_epoch()).count();
    }

    /**
     * @brief This method returns the resulting time
     */
    static inline double ElapsedSeconds(const std::chrono::steady_clock::time_point StartTime)
    {
        const auto current_time = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<std::chrono::duration<double>>(current_time - StartTime).count();
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This method sets the output file *.time that will store the timing
     * @param rOutputFileName The name of the output file
     */
    static int SetOuputFile(std::string const& rOutputFileName);

    /**
     * @brief This method closes the output file
     */
    static int CloseOuputFile();

    /**
     * @brief This method gets the variable which stores if the information is printed on screen
     * @return True if the information is printed on screen, false otherwise
     */
    static bool GetPrintOnScreen();

    /**
     * @brief This method sets the variable which stores if the information is printed on screen
     * @param PrintOnScreen True if the information is printed on screen, false otherwise
     */
    static void SetPrintOnScreen(bool const PrintOnScreen);

    /**
     * @brief This method gets the variable which stores if the information is printed on each interval
     * @return True if the information is printed on each interval, false otherwise
     */
    static bool GetPrintIntervalInformation();

    /**
     * @brief This method sets the variable which stores if the information is printed on each interval
     * @param PrintIntervalInformation True if the information is printed on each interval, false otherwise
     */
    static void SetPrintIntervalInformation(bool const PrintIntervalInformation);

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief This method prints the internal information in a given stream
     * @param rOStream The strem considered
     * @param rIntervalName The internal name that will store the timing data
     * @param StartTime The starting time
     * @param StopTime The stoping time
     */
    static void PrintIntervalInformation(
        std::ostream& rOStream,
        std::string const& rIntervalName,
        const double StartTime,
        const double StopTime
        );

    /**
     * @brief This method prints the internal information
     * @param rIntervalName The internal name that will store the timing data
     * @param StartTime The starting time
     * @param StopTime The stoping time
     */
    static void PrintIntervalInformation(
        std::string const& rIntervalName,
        const double StartTime,
        const double StopTime
        );

    /**
     * @brief This method prints the timing information
     */
    static void PrintTimingInformation();

    /**
     * @brief This method prints the timing information in a giving stream
     * @param rOStream The strem considered
     */
    static void PrintTimingInformation(std::ostream& rOStream);

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Timer";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        PrintTimingInformation(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:

    ///@name Static Member Variables
    ///@{

    static InternalNameDatabaseType msInternalNameDatabase;        /// The names used on the time tables

    static ContainerType msTimeTable;                              /// The time tables

    static std::ofstream msOutputFile;                             /// The file to be written

    static std::size_t msCounter;                                  /// Counter of the instances

    static bool msPrintOnScreen;                                   /// If the information is printed on screen

    static bool msPrintIntervalInformation;                        /// If the information of the interval is printed

    static const std::chrono::steady_clock::time_point mStartTime; /// The starting time

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief This method returns the internal name used on the table
     * @param rName The base name
     * @return The internal name
     */
    static std::string GetInternalName(const std::string& rName);

    static std::vector<std::string>& GetLabelsStackInstance()
    {
      static std::vector<std::string> instance;
      return instance;
    }

    static std::string CreateFullLabel(){
        const auto& r_labels_stack = GetLabelsStackInstance();
        std::string result;
        for(const auto& r_label : r_labels_stack){
        result += "/" + r_label;
        }
        return result;
    }

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    Timer& operator=(Timer const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    /*       Timer(Timer const& rOther); */

    ///@}

}; // Class Timer

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{


/// input stream function
//   inline std::istream& operator >> (std::istream& rIStream,
// 				    Timer& rThis){}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Timer& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_TIMER_H_INCLUDED  defined


