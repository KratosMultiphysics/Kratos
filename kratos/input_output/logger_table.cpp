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

// System includes
#include <sstream>

// External includes 

// Project includes
#include "input_output/logger_table.h"

namespace Kratos
{
std::string LoggerTable::Info() const {
    return "LoggerTable";
}

/***********************************************************************************/
/***********************************************************************************/

void LoggerTable::PrintInfo(std::ostream& rOStream) const {
    rOStream << Info();
}

/***********************************************************************************/
/***********************************************************************************/

void LoggerTable::PrintData(std::ostream& rOStream) const {
    rOStream << GetMessage();
}

/***********************************************************************************/
/***********************************************************************************/

LoggerTable& LoggerTable::operator << (const char * rString) {
    mMessage << rString;

    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

LoggerTable& LoggerTable::operator << (std::ostream& (*pf)(std::ostream&)) {
    std::stringstream buffer;
    pf(buffer);
    mMessage << buffer.str();

    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

LoggerTable& LoggerTable::operator << (CodeLocation const& TheLocation) {
    if (typeid(TheLocation) != typeid(GetLocation())){
        SetLocation(TheLocation);
        EndTable();
        ClearTable();
    }

    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

LoggerTable& LoggerTable::operator << (Severity const& TheSeverity) {
    if (typeid(TheSeverity) != typeid(GetSeverity())){
        SetSeverity(TheSeverity);
        EndTable();
        ClearTable();
    }

    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

LoggerTable& LoggerTable::operator << (Category const& TheCategory) {
    if (typeid(TheCategory) != typeid(GetCategory())){
        SetCategory(TheCategory);
        EndTable();
        ClearTable();
    }

    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

std::ostream& operator << (std::ostream& rOStream,
    const LoggerTable& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos.


