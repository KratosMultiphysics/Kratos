//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                    
//

#if !defined(KRATOS_TABLE_STREAM_UTILITY)
#define KRATOS_TABLE_STREAM_UTILITY

// System includes

// External includes

// Project includes
#include "includes/serializer.h"
#include "includes/table_stream.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    
    typedef TableStream TableStreamType;
    
///@}
///@name  Enum's
///@{
    
///@}
///@name  Functions
///@{
    
class TableStreamUtility
{
public:
    ///@name Type Definitions
    ///@{
    
    /// Counted pointer of TableStreamUtility
    KRATOS_CLASS_POINTER_DEFINITION( TableStreamUtility );
    
    ///@}
    ///@name Life Cycle
    ///@{

    // Class Constructor
    /**
     * The default constructor
     */
    
    TableStreamUtility(const bool UseBoldFont = true):
        mTable(&std::cout, "|", UseBoldFont)
    {
    }
    
    virtual ~TableStreamUtility()= default;
    
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
    ///@name Operations
    ///@{
    
    /**
     * It prints the header of the table
     */
    
    void PrintHeader()
    {
        mTable.PrintHeader();
    }
    
    /**
     * It prints the footer of the table 
     */
    
    void PrintFooter()
    {
        mTable.PrintFooter();
    }
    
    /**
     * It adds a value to the ouput of the current row 
     */
    
    template<class TClass>
    void AddToRow(TClass ThisClass)
    {
        mTable << ThisClass;
    }
    
    /**
     * It adds a column to the table
     * @param ThisName The name of the variable
     * @param ThisSpaces The number of spaces to consider
     */
        
    void AddColumn(std::string ThisName, unsigned int ThisSpaces)
    {
        mTable.AddColumn(ThisName, ThisSpaces);
    }
    
    /**
     * It returns the table of BPrinter
     * @return mTable: The table stream table
     */
        
    TableStreamType& GetTable()
    {
        return mTable;
    }
    
    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "TableStreamUtility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info() << std::endl;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << Info() << std::endl;
    }

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
    
    TableStreamType mTable;
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
        rSerializer.save("Table", mTable);
    }

    void load(Serializer& rSerializer)
    {
        rSerializer.load("Table", mTable);
    }

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
};// class TableStreamUtility

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.
#endif /* KRATOS_TABLE_STREAM_UTILITY defined */
