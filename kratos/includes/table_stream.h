//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//  Inspired in bprinter, work of Dat Chu: https://github.com/dattanchu/bprinter
//  Removing all the dependencies of boost::karma and adding bold fonts and additional functionalities needed
//

#ifndef KRATOS_TABLE_STREAM_H_INCLUDED
#define KRATOS_TABLE_STREAM_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

namespace Kratos 
{
///@addtogroup Kratos Core
///@{

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
    
class endl{};

/** \class TableStream
  This is a fancy table to stream data in a fancy way
  */
class TableStream
{
public:
    
    ///@name Type Definitions
    ///@{
    
    ///@}
    ///@name Life Cycle
    ///@{
    
    /// Default constructor.
    TableStream(
        std::ostream * Output, 
        const std::string Separator = "|", 
        const bool UseBoldFont = true
        ) : mOutStream(Output),
            mSeparator(Separator),
            mBoldFont(UseBoldFont)
    {
        // Initialize values
        mIndexRow    = 0;
        mIndexColumn = 0;
        mTableWidth  = 0;
        mFlushLeft = false;
    }
    
    /// Destructor.
    virtual ~TableStream() = default;

    ///@}
    ///@name Operators
    ///@{

    /**
     * This is the operator << for any kind of input
     * @param Input: The input considered
     * @return The updated table stream
     */
    template<typename TClass> 
    TableStream& operator<<(TClass Input)
    {
        if (typeid(TClass) == typeid(endl))
        {
            while (mIndexColumn != 0)
            {
                *this << "";
            }
        }
        else
        {
            if (mIndexColumn == 0)
            {
                *mOutStream << "|";
            }
            
            if(mFlushLeft)
            {
                *mOutStream << std::left;
            }
            else
            {
                *mOutStream << std::right; 
            }

            // Leave 3 extra space: One for negative sign, one for zero, one for decimal
            *mOutStream << std::setw(mColumnWidths.at(mIndexColumn)) << Input;

            if (mIndexColumn == GetNumColumns()-1)
            {
                *mOutStream << "|\n";
                mIndexRow = mIndexRow + 1;
                mIndexColumn = 0;
            } 
            else 
            {
                *mOutStream << mSeparator;
                mIndexColumn = mIndexColumn + 1;
            }
        }
        
        return *this;
    }
    
    /**
     * This is the operator << just for floats
     * @param Input: The float considered
     * @return The updated table stream
     */
    TableStream& operator<<(float Input)
    {
        OutputDecimalNumber<float>(Input);
        return *this;
    }
    
    /**
     * This is the operator << just for doubles
     * @param Input: The double considered
     * @return The updated table stream
     */
    TableStream& operator<<(double Input)
    {
        OutputDecimalNumber<double>(Input);
        return *this;
    }
    
    ///@}
    ///@name Operations
    ///@{
    
    /**
     * It returns the number of columns 
     * @return The size of mColumnHeaders (the column headers)
     */
    unsigned int GetNumColumns() const
    {
        return mColumnHeaders.size();
    }

    /**
     * It returns the table width
     * @return mTableWidth: The table width
     */
    unsigned int GetTableWidth() const
    {
        return mTableWidth;
    }
    
    /**
     * Set the separator used for the table
     */
    void SetSeparator(const std::string& Separator)
    {
        mSeparator = Separator;
    }
    
    /**
     * Set if the bold fonts are used for the table
     */
    void SetBold(const bool& UseBoldFont)
    {
        mBoldFont = UseBoldFont;
    }
    
    /**
     * Set the flush orientation to the left
     */
    void SetFlushLeft()
    {
        mFlushLeft = true;
    }
    
    /**
     * Set the flush orientation to the right
     */
    void SetFlushRight()
    {
        mFlushLeft = false;
    }

    /** Add a column to our table
     * @param HeaderName: Name to be print for the header
     * @param ColumnWidth: The width of the column must be at least 4 spaces
     */
    void AddColumn(
        const std::string& HeaderName, 
        const int ColumnWidth
        )
    {
        if (ColumnWidth < 4)
        {
            KRATOS_ERROR << "Column size has to be >= 4" << std::endl;
        }

        mColumnHeaders.push_back(HeaderName);
        mColumnWidths.push_back(ColumnWidth);
        mTableWidth += ColumnWidth + mSeparator.size(); // for the separator  
    }
    
    /**
     * This function prints the header of the stream
     */
    void PrintHeader()
    {
        PrintHorizontalLine();
        
        if (mBoldFont == true)
        {
        #if !defined(_WIN32)
            *mOutStream << "\e[1m";
        #endif
        }

        *mOutStream << "|";
            
        for (unsigned int i = 0; i < GetNumColumns(); ++i)
        {
            if(mFlushLeft)
            {
                *mOutStream << std::left;
            }
            else
            {
                *mOutStream << std::right; 
            }

            *mOutStream << std::setw(mColumnWidths.at(i)) << mColumnHeaders.at(i).substr(0, mColumnWidths.at(i));
            
            if (i != GetNumColumns()-1)
            {
                *mOutStream << mSeparator;
            }
        }

        *mOutStream << "|";

        if (mBoldFont == true)
        {
        #if !defined(_WIN32)
            *mOutStream << "\e[0m";
        #endif
        }

        *mOutStream << "\n";
            
        PrintHorizontalLine();
    }
    
    /**
     * This function prints the footer of the stream
     */
    void PrintFooter()
    {
        PrintHorizontalLine();
    }
  
    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{
    
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


    ///@}
    ///@name Member Variables
    ///@{
    
    // Stream related variables
    std::ostream* mOutStream;                // The stream considered 
    std::vector<std::string> mColumnHeaders; // This vector contains the header to print in each column
    std::vector<int> mColumnWidths;          // This vector containts the spaces of each column
    std::string mSeparator;                  // The separator considered in each column

    // The indexes currently used
    unsigned int mIndexRow;                  // Index of current row
    unsigned int mIndexColumn;               // Index of current column

    // Other variables related with the table output
    unsigned int mTableWidth;                // The table width
    bool mFlushLeft;                         // If the flush is aligned to the left or the right
    bool mBoldFont;                          // If the bold fonts are considered to use in the 
    
    ///@}
    ///@name Private Operators
    ///@{



    ///@}
    ///@name Private Operations
    ///@{
    
    /**
     * This functions prints an horizontal line
     */
    void PrintHorizontalLine()
    {
        *mOutStream << "+"; // the left bar

        for (unsigned int i = 0; i< mTableWidth-1; ++i)
        {
            *mOutStream << "-";
        }

        *mOutStream << "+"; // the right bar
        *mOutStream << "\n";
    }
    
    /**
     * This functions prints into the stream a double or float value with scientific notation (respeting the spaces asigned)
     * @param Input: The double or float to print
     */
    template<typename TClass> 
    void OutputDecimalNumber(TClass Input)
    {
        // If we cannot handle this number, indicate so
        if (Input < 10*(mColumnWidths.at(mIndexColumn)-1) || Input > 10*mColumnWidths.at(mIndexColumn))
        {
            std::stringstream string_out;
            string_out 
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(3)
            << std::uppercase
            << std::setw(mColumnWidths.at(mIndexColumn))
            << Input;

            std::string string_to_print = string_out.str();

            *mOutStream << string_to_print;
        } 
        else 
        {
            *mOutStream 
            << std::setiosflags(std::ios::scientific)
            << std::setprecision(3)
            << std::uppercase
            << std::setw(mColumnWidths.at(mIndexColumn))
            << Input;
        }

        if (mIndexColumn == GetNumColumns()-1)
        {
            *mOutStream << "|\n";
            mIndexRow = mIndexRow + 1;
            mIndexColumn = 0;
        } 
        else 
        {
            *mOutStream << mSeparator;
            mIndexColumn = mIndexColumn + 1;
        }
    }
    
    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{
}; // Class TableStream

} // namespace Kratos.
#endif
