//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:       BSD License
//                Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//


#if !defined(KRATOS_LOGGER_TABLE_H_INCLUDED )
#define  KRATOS_LOGGER_TABLE_H_INCLUDED

// System includes
#include <map> // NOTE: MUST BE ORDERED

// External includes

// Project includes
#include "includes/table_stream.h" // TODO: Think about move it to input_output folder instead
#include "input_output/logger_message.h"

namespace Kratos
{
    ///@addtogroup KratosCore
    ///@{

    ///@name Kratos Classes
    ///@{

    /// LoggerTable class holdes message and the properties of the message in a table form.
    /** LoggerTable holds the origin of the message, severity, level and 
        the category of it.
        Most of the methods are defined in header to be inlined in order to
        increase the performance. 
    */
    class KRATOS_API(KRATOS_CORE) LoggerTable : public LoggerMessage
    {
    public:
        ///@name Type Definitions
        ///@{

        typedef LoggerMessage BaseType;
        
        typedef BaseType::Severity Severity;
        
        typedef BaseType::Category Category;
        
        typedef std::map<std::string, bool> MapStringsType;
        
        typedef TableStream TableStreamType;
        
        using TimePointType = std::chrono::steady_clock::time_point;

        ///@}
        ///@name Enums
        ///@{

        ///@}
        ///@name Life Cycle
        ///@{

        LoggerTable(std::string const& TheLabel) : BaseType(TheLabel){
            mTableIsStarted = false;
        }

        LoggerTable(LoggerTable const& Other) : BaseType(Other){
            mTableIsStarted = false;
        }

        /// Destructor.
        virtual ~LoggerTable() {}

        ///@}
        ///@name Operators
        ///@{

        LoggerTable& operator=(LoggerTable const& Other){
            mLabels = Other.mLabels;
            mMessage = Other.mMessage;
            mTableIsStarted = Other.mTableIsStarted;

            return *this;
        }

        ///@}
        ///@name Operations
        ///@{

        /**
         * It clears the table // TODO: Maybe a std::vector containing several tables can be craetes for each category and location. Ask Pooyan
         */
        void ClearTable()
        {
            mLabels.clear();
        }
        
        /**
         * It starts the table
         */
        void StartTable()
        {
            if(!mTableIsStarted){
                mMessage.PrintHeader();
                for (auto& label : mLabels){
                    label.second = false; // We reset the "print"
                }
                mTableIsStarted = true;
            }
        }
        
        /**
         * It finishs the table
         */
        void EndTable()
        {
            mMessage.PrintFooter();
            mTableIsStarted = false;
        }
        
        ///@}
        ///@name Access
        ///@{

        /**
         * It sets the label. For that pouporse a set is created in order to check that a label is just printed once and it is printed in a correct order TODO: See how to do this
         * @param TheLabel The label that is set
         */
        void SetLabel(std::string const& TheLabel) {
            bool existing_label = false;
            MapStringsType::iterator set = mLabels.find(TheLabel);
            if(set != mLabels.end())
                existing_label = true;
            if (!existing_label){
                mLabels.insert({TheLabel, false});
                size_t n = 0;
                for (char charac : TheLabel) n++;
                mMessage.AddColumn(TheLabel, n);
            }
            BaseType::SetLabel(TheLabel);
        }
        
        /**
         * Adds to the current string a new message. It check if the current label has been streamed already
         * @param TheMessage The message to be streamed
         */
        void SetMessage(std::string const& TheMessage) {
            if (!mLabels[GetLabel()]){
                mMessage << TheMessage;
                mLabels[GetLabel()] = true; // NOTE: The bool says if already streamed
            }
        }

        /**
         * It gets you the current stream
         */
        std::string const& GetMessage() const {
            return mMessage.GetCurrentStream();
        }
    
        ///@}
        ///@name Inquiry
        ///@{

        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a string.
        virtual std::string Info() const;

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const;

        /// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const;

        /// string stream function
        template<class StreamValueType>
        LoggerTable& operator << (StreamValueType const& rValue)
        {
            if (!mLabels[GetLabel()]){
                mMessage << rValue;
                mLabels[GetLabel()] = true; // NOTE: The bool says if already streamed
            }

            return *this;
        }

        /// Manipulator stream function
        LoggerTable& operator << (std::ostream& (*pf)(std::ostream&));

        /// char stream function
        LoggerTable& operator << (const char * rString);
        
        /// Location stream function
        LoggerTable& operator << (CodeLocation const& TheLocation);

        /// Severity stream function
        LoggerTable& operator << (Severity const& TheSeverity);

        /// Category stream function
        LoggerTable& operator << (Category const& TheCategory);

        ///@}

    private:
        ///@name Life Cycle
        ///@{

        ///@}
        ///@name Member Variables
        ///@{

        MapStringsType mLabels;     // The labels used (used to create the header of the table)
        TableStreamType mMessage;   // The table used to print data
        bool mTableIsStarted;       // A "flag" used to check that the table has beed started

        ///@}
    }; // Class LoggerTable

    ///@}

    ///@name Input and output
    ///@{

    /// output stream function
    std::ostream& operator << (std::ostream& rOStream,
        const LoggerTable& rThis);

    ///@}
    ///@name macros
    ///@{


    ///@}

    ///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_LOGGER_TABLE_H_INCLUDED  defined 
