//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

#ifndef LOGGER_INCLUDE_HPP
#define LOGGER_INCLUDE_HPP

//// STL includes
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iomanip>

namespace queso {

///@name QuESo Classes
///@{

///
/**
 * @class  Logger
 * @author Manuel Messmer
*/
class Logger {
public:
    ///@}
    ///@name Life cycle
    ///@{

    /// Default Constructor
    Logger() {}

    /// Constructor
    Logger(std::string rWhere ) {
        AppendMessage( rWhere );
        AppendMessage( " -- " );

    }

    /// Destructor
    ~Logger() {
        std::cout << mMessage;
    }

    ///@}
    ///@name Operations
    ///@{

    Logger& operator << (std::ostream& (*pf)(std::ostream&)) {
       	std::stringstream buffer;
		pf(buffer);

        AppendMessage(buffer.str());

        return *this;
    }

    template<class StreamValueType>
    Logger& operator << (StreamValueType const& rValue)
    {
        std::stringstream buffer;
        buffer << std::setprecision(10) << rValue;

        AppendMessage(buffer.str());

        return *this;
    }

    void AppendMessage(std::string const& rMessage) {
		mMessage.append(rMessage);
	}

    Logger& operator << (const char * rString) {
        AppendMessage(rString);
        return *this;
    }
    ///@}

private:
    ///@}
    ///@name Private Members
    ///@{
    std::string mMessage;

    ///@}
};

#define QuESo_INFO Logger()
#define QuESo_INFO_IF(Conditional) if(Conditional) Logger()

///@} // End QuESo Classes
} // End namespace queso

#endif // LOGGER_INCLUDE_HPP