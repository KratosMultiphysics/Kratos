//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

#ifndef EXCEPTION_INCLUDE_HPP
#define EXCEPTION_INCLUDE_HPP

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
 * @class  Exception
 * @author Manuel Messmer
*/
class Exception : public std::exception {
public:
    ///@}
    ///@name Life cycle
    ///@{

    /// Constructor
    Exception(std::string const& rFileName, std::string const& rFunctionName, std::size_t LineNumber) {
        AppendMessage( "\nWhere:" );
        AppendMessage( "\tFilename: " );
        AppendMessage( rFileName );
        AppendMessage( ", Line: " );
        AppendMessage( std::to_string(LineNumber) );
        AppendMessage( "\n\tFunction: " );
        AppendMessage( rFunctionName );
        AppendMessage( "\nWhat:\t" );
    }

    ///@}
    ///@name Operations
    ///@{

    template<class StreamValueType>
    Exception& operator << (StreamValueType const& rValue)
    {
        std::stringstream buffer;
        buffer << rValue;

        AppendMessage(buffer.str());

        return *this;
    }

    Exception& operator << (std::ostream& (*pf)(std::ostream&)) {
       	std::stringstream buffer;
		pf(buffer);

        AppendMessage(buffer.str());

        return *this;
    }

    void AppendMessage(std::string const& rMessage) {
		mMessage.append(rMessage);
	}

    Exception& operator << (const char * rString) {
        AppendMessage(rString);
        return *this;
    }

    const char* what() const noexcept override {
        return mMessage.c_str();
    }

private:
    ///@}
    ///@name Private Members
    ///@{
    std::string mMessage;
    ///@}
};

#if defined(__PRETTY_FUNCTION__)
#define QuESo_CURRENT_FUNCTION __PRETTY_FUNCTION__
#elif defined(__GNUC__)
#define QuESo_CURRENT_FUNCTION __PRETTY_FUNCTION__
#elif defined(__FUNCTION__)
#define QuESo_CURRENT_FUNCTION __FUNCTION__
#elif defined(__func__)
#define QuESo_CURRENT_FUNCTION __func__
#else
#define QuESo_CURRENT_FUNCTION "unknown function"
#endif

#define QuESo_ERROR throw Exception(__FILE__, QuESo_CURRENT_FUNCTION, __LINE__)
#define QuESo_ERROR_IF(Conditional) if(Conditional) throw Exception(__FILE__, QuESo_CURRENT_FUNCTION, __LINE__)

///@} // End QuESo Classes
} // End namespace queso

#endif // EXCEPTION_INCLUDE_HPP