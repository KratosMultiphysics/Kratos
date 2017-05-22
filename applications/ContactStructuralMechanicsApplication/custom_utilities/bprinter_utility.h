// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_BPRINTER_UTILITY)
#define KRATOS_BPRINTER_UTILITY

// System includes

// External includes

// Project includes
#include "custom_external_libraries/bprinter/table_printer.h"

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
    
class BprinterUtility
{
public:
    ///@name Type Definitions
    ///@{
    
    ///@}
    ///@name Life Cycle
    ///@{

    // Class Constructor
    /**
     * The default constructor
     */
    
    BprinterUtility()
    {        
    }
    
    virtual ~BprinterUtility(){};
    
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
    
    /*
     * It prints the header of the table
     */
    
    void PrintHeader()
    {
        mTable.PrintHeader();
    }
    
    /*
     * It prints the footer of the table 
     */
    
    void PrintFooter()
    {
        mTable.PrintFooter();
    }
    
    /*
     * It adds a value to the ouput of the current row 
     */
    
    template<class TClass>
    void AddToRow(TClass ThisClass)
    {
        mTable << ThisClass;
    }
    
    /*
     * It adds a column to the table
     */
        
    void AddColumn(std::string ThisName, unsigned int ThisSpaces)
    {
        std::string ThisNameBold = "\e[1m";
        ThisNameBold.append(ThisName);
        ThisNameBold.append("\e[0m");
        mTable.AddColumn(ThisNameBold, ThisSpaces);
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
    
    bprinter::TablePrinter mTable{&std::cout, "|"};
    
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

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
};// class BprinterUtility

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.
#endif /* KRATOS_BPRINTER_UTILITY defined */
